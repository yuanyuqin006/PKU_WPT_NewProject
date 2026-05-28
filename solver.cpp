// solver.cpp —— 升级版：基于矩阵主动解耦的 4 线圈阵列 2D 扫描
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <limits>

#ifndef __cpp_lib_math_special_functions
namespace std {
    inline double comp_ellint_1(double k) {
        constexpr double pi_v = 3.14159265358979323846;
        const double ak = (k < 0 ? -k : k);
        if (ak >= 1.0) return std::numeric_limits<double>::infinity();
        if (ak == 0.0) return pi_v * 0.5;
        double a = 1.0;
        double b = sqrt(1.0 - ak * ak);
        for (int i = 0; i < 100; ++i) {
            const double an = 0.5 * (a + b);
            const double bn = sqrt(a * b);
            if ((a > b ? a - b : b - a) <= 1e-16 * an) { a = an; break; }
            a = an; b = bn;
        }
        return pi_v / (2.0 * a);
    }
    inline double comp_ellint_2(double k) {
        constexpr double pi_v = 3.14159265358979323846;
        const double ak = (k < 0 ? -k : k);
        if (ak == 0.0) return pi_v * 0.5;
        if (ak >= 1.0) return 1.0;
        double a = 1.0;
        double b = sqrt(1.0 - ak * ak);
        double c = ak;
        double pow2 = 1.0;
        double s = c * c;
        for (int i = 0; i < 100; ++i) {
            const double an = 0.5 * (a + b);
            const double bn = sqrt(a * b);
            c    = 0.5 * (a - b);
            pow2 *= 2.0;
            s    += pow2 * c * c;
            if ((a > b ? a - b : b - a) <= 1e-16 * an) { a = an; break; }
            a = an; b = bn;
        }
        const double K = pi_v / (2.0 * a);
        return K * (1.0 - 0.5 * s);
    }
}
#endif

constexpr double PI  = 3.14159265358979323846;
constexpr double MU0 = 4.0 * PI * 1e-7;

// ---- 物理参数 ----
constexpr double A_INNER = 0.045;
constexpr double A_OUTER = 0.055;
constexpr int    N_A     = 5;

constexpr double B_INNER = 0.045;
constexpr double B_OUTER = 0.055;
constexpr int    N_B     = 5;

constexpr int    N_PHI   = 360;
constexpr double H_AXIAL = 5.0e-3; // Qi 气隙 5mm

// 高斯-约旦法 4x4 矩阵求逆助手
bool invert_4x4(const double src[4][4], double dst[4][4]) {
    double mat[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mat[i][j] = src[i][j];
            mat[i][j+4] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < 4; ++i) {
        double pivot = mat[i][i];
        if (std::abs(pivot) < 1e-15) return false;
        for (int j = 0; j < 8; ++j) mat[i][j] /= pivot;
        for (int k = 0; k < 4; ++k) {
            if (k != i) {
                double factor = mat[k][i];
                for (int j = 0; j < 8; ++j) mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) dst[i][j] = mat[i][j+4];
    }
    return true;
}

// 核心积分求解器
double calc_M_ring(double a, double b, double h, double d) {
    const double dphi = 2.0 * PI / static_cast<double>(N_PHI);
    double sum = 0.0;
    for (int i = 0; i < N_PHI; ++i) {
        const double phi  = (i + 0.5) * dphi;
        const double cphi = std::cos(phi);
        const double rho  = std::sqrt(b * b + d * d + 2.0 * b * d * cphi);
        if (rho < 1e-9) continue;

        const double denom = (a + rho) * (a + rho) + h * h;
        const double k2    = 4.0 * a * rho / denom;
        const double k     = std::sqrt(k2);

        const double K = std::comp_ellint_1(k);
        const double E = std::comp_ellint_2(k);
        const double bracket    = (1.0 - 0.5 * k2) * K - E;
        const double Aphi_per_I = (MU0 / (PI * k)) * std::sqrt(a / rho) * bracket;

        sum += Aphi_per_I * b * (b + d * cphi) / rho;
    }
    return sum * dphi;
}

double calc_M_total(double a_in, double a_out, int Na, double b_in, double b_out, int Nb, double h, double d) {
    double M = 0.0;
    for (int i = 0; i < Na; ++i) {
        const double ai = a_in + (a_out - a_in) * (i + 0.5) / static_cast<double>(Na);
        for (int j = 0; j < Nb; ++j) {
            const double bj = b_in + (b_out - b_in) * (j + 0.5) / static_cast<double>(Nb);
            M += calc_M_ring(ai, bj, h, d);
        }
    }
    return M;
}

int main() {
    // ---- 1. 坐标与工艺叠层设定 ----
    constexpr double R_TX_MM = 1.0e3 * 0.5 * (A_INNER + A_OUTER); // 50 mm
    constexpr double TX_POS_MM[4][2] = {
        {0.0,  R_TX_MM},   // Tx1 上 (Layer 0)
        {0.0, -R_TX_MM},   // Tx2 下 (Layer 0)
        {-R_TX_MM, 0.0},   // Tx3 左 (Layer 1)
        { R_TX_MM, 0.0},   // Tx4 右 (Layer 1)
    };

    // ---- 2. 预计算发射阵列自感与互感全矩阵 L_tx ----
    double L_tx[4][4] = {0.0};
    // 计算单线圈自感 (利用等效线径高斯微偏置 0.8mm 替代奇异点)
    double L_self = calc_M_total(A_INNER, A_OUTER, N_A, A_INNER, A_OUTER, N_A, 0.8e-3, 0.0);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i == j) {
                L_tx[i][j] = L_self + 1.0e-7;
            } else {
                double dx = (TX_POS_MM[i][0] - TX_POS_MM[j][0]) * 1.0e-3;
                double dy = (TX_POS_MM[i][1] - TX_POS_MM[j][1]) * 1.0e-3;
                double d_tx = std::sqrt(dx * dx + dy * dy);
                // 工艺错层：如果跨层(一正交一垂直)，引入 1.2mm PCB板厚间距
                double h_layer = ((i < 2 && j >= 2) || (i >= 2 && j < 2)) ? 1.2e-3 : 0.0;
                L_tx[i][j] = calc_M_total(A_INNER, A_OUTER, N_A, A_INNER, A_OUTER, N_A, h_layer, d_tx);
            }
        }
    }

    // 矩阵求逆
    double inv_L_tx[4][4] = {0.0};
    if (!invert_4x4(L_tx, inv_L_tx)) {
        std::cerr << "Fatal Error: L_tx matrix is singular!" << std::endl;
        return 1;
    }

    // ---- 3. 2D 网格扫描 ----
    constexpr double GRID_MIN  = -50.0;
    constexpr double GRID_MAX  =  50.0;
    constexpr double GRID_STEP =   0.5;
    const int N = static_cast<int>(std::round((GRID_MAX - GRID_MIN) / GRID_STEP)) + 1;
    const long long TOTAL = static_cast<long long>(N) * N;

    std::ofstream f("array_data.json");
    f << std::scientific << std::setprecision(8);
    f << "{\"grid_size\": " << GRID_STEP << ", \"data\": [";

    bool first = true;
    for (int iy = 0; iy < N; ++iy) {
        const double y_mm = GRID_MIN + GRID_STEP * iy;
        const double y_m  = y_mm * 1.0e-3;

        for (int ix = 0; ix < N; ++ix) {
            const double x_mm = GRID_MIN + GRID_STEP * ix;
            const double x_m  = x_mm * 1.0e-3;

            // 基础互感向量 M_txrx
                        double M[4];
                        for (int k = 0; k < 4; ++k) {
                            const double dx = x_m - TX_POS_MM[k][0] * 1.0e-3;
                            const double dy = y_m - TX_POS_MM[k][1] * 1.0e-3;
                            
                            // 为 Tx3 (k=2) 和 Tx4 (k=3) 引入 1.2mm 的工艺错层对空高度修正
                            const double h_curr = (k >= 2) ? (H_AXIAL + 1.2e-3) : H_AXIAL;
                            
                            M[k] = calc_M_total(A_INNER, A_OUTER, N_A, B_INNER, B_OUTER, N_B, h_curr, std::sqrt(dx*dx + dy*dy));
                        }

            // ============================================================
                        // A. 原始算法：最大比合并（物理修正：回归真实的等效耦合互感模长）
                        // ============================================================
                        const double s2 = M[0]*M[0] + M[1]*M[1] + M[2]*M[2] + M[3]*M[3];
                        double w[4];
                        for (int k = 0; k < 4; ++k) w[k] = (s2 > 0.0) ? (M[k]*M[k]) / s2 : 0.25;
                        
                        // 【物理修正】MRC驱动下，归一化总电流后的真实等效互感就是 M 向量的几何模长
                        const double eff_total = std::sqrt(s2);

                        // ============================================================
                        // B. 新算法：主动矩阵零迫解耦权重（物理修正：考虑电流方向/相位线性叠加）
                        // ============================================================
                        double w_dec[4] = {0.0};
                        for (int i = 0; i < 4; ++i) {
                            for (int j = 0; j < 4; ++j) {
                                w_dec[i] += inv_L_tx[i][j] * M[j];
                            }
                        }
                        // 将解耦驱动流转换为归一化功率分布占比（保持传给前端的百分比全为正数，供 UI 进度条使用）
                        double s2_dec = w_dec[0]*w_dec[0] + w_dec[1]*w_dec[1] + w_dec[2]*w_dec[2] + w_dec[3]*w_dec[3];
                        double w_dn[4];
                        for (int k = 0; k < 4; ++k) w_dn[k] = (s2_dec > 0.0) ? (w_dec[k]*w_dec[k]) / s2_dec : 0.25;

                        // 【物理修正】核心核心：计算真实物理场的线性叠加（携带 w_dec 的正负号！）
                        double linear_field_sum = w_dec[0]*M[0] + w_dec[1]*M[1] + w_dec[2]*M[2] + w_dec[3]*M[3];
                        
                        // 真实物理耦合效率 = 矢量叠加后的总磁链 / 发射总电流有效值
                        const double eff_decoupled = (s2_dec > 0.0) ? (std::abs(linear_field_sum) / std::sqrt(s2_dec)) : 0.0;
            // ---- 4. 输出到 JSON ----
            if (!first) f << ",";
            first = false;
            f << "[" << x_mm << "," << y_mm
              << "," << M[0] << "," << M[1] << "," << M[2] << "," << M[3]     // 2~5
              << "," << w[0] << "," << w[1] << "," << w[2] << "," << w[3]     // 6~9
              << "," << eff_total                                             // 10
              << "," << w_dn[0] << "," << w_dn[1] << "," << w_dn[2] << "," << w_dn[3] // 11~14 (新)
              << "," << eff_decoupled << "]";                                 // 15 (新)
        }
    }
    f << "]}\n";
    f.close();
    std::cout << "Successfully generated decoupled matrix scan!" << std::endl;
    return 0;
}
