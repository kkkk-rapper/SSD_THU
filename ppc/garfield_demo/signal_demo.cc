// =============================================================================
// signal_demo.cc
//
// 最小可运行的 Garfield++ 半导体信号仿真示例。
// 加载 SSD 导出的电场 + P+ 权势，在 (r=5 mm, z=5 mm) 放一个 122 keV 全沉积事件，
// 用 AvalancheMC 漂移 e/h 对，输出 P+ 上的诱导电荷 Q(t)。
//
// 与 SSD 的对照点：
//   /root/ssd_projects/ppc/output/signals/wf_axis_z5.csv  (r=0, z=5mm)
//   注意位置不完全一致 —— 本程序用 r=5mm 是为了和 SSD 的 mid_r6 (r=6mm) 对比
//
// !! 重要警告 !!
//   Garfield++ 没有原生 MediumGermanium。本程序用 MediumSilicon 占位，
//   并手动覆盖低场迁移率为 Ge @ 78K 的近似值。结果只能做"工作流验证"，
//   不能直接做 quantitative HPGe 信号研究。
//
// 编译运行：
//   source /root/gpp/install/share/Garfield/setupGarfield.sh
//   make
//   ./signal_demo
// =============================================================================

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"

using namespace Garfield;

int main() {

    // -------------------------------------------------------------------------
    // 1. Medium：用 MediumSilicon 占位，灌入 Ge @ 78K 的迁移率
    //    Garfield++ 迁移率单位是 cm^2 / (V * ns)
    //    Ge 78K:  mu_e ≈ 36000 cm^2/(V*s) = 3.6e-5 cm^2/(V*ns)
    //             mu_h ≈ 42000 cm^2/(V*s) = 4.2e-5 cm^2/(V*ns)
    //    Saturation velocity ≈ 1e7 cm/s = 1e-2 cm/ns
    // -------------------------------------------------------------------------
    MediumSilicon ge;
    ge.SetTemperature(78.0);
    ge.SetLowFieldMobility(3.6e-5, 4.2e-5);          // Ge 78K 近似
    ge.SetSaturationVelocity(1.0e-2, 1.0e-2);        // cm/ns
    ge.SetHighFieldMobilityModelCanali();            // 高场饱和模型
    std::cout << "[info] Medium configured as Ge-like Si placeholder.\n";

    // -------------------------------------------------------------------------
    // 2. 加载 SSD 导出的电场 + 静电势
    //    网格参数必须和导出文件一致 (181 x 1 x 141, cm 单位)
    // -------------------------------------------------------------------------
    // 注：Garfield++ 的 SetMesh 拒绝 ymin==ymax。轴对称下用 ny=1 + ymax=2π
    // 表示"单 phi 切片代表整圈"，Garfield 会自动启用 theta 周期性。
    constexpr double TWO_PI = 6.283185307179586;
    ComponentGrid efield;
    efield.SetCylindricalCoordinates();
    efield.SetMesh(181, 1, 141,
                   0.0, 1.8,        // r ∈ [0, 1.8] cm
                   0.0, TWO_PI,     // phi 完整一圈，ny=1 表示轴对称
                   -0.2, 1.2);      // z ∈ [-0.2, 1.2] cm
    if (!efield.LoadElectricField("../output/garfield/efield_1000V.txt",
                                  "xyz", true, false)) {
        std::cerr << "[fatal] LoadElectricField failed\n";
        return 1;
    }
    efield.SetMedium(&ge);
    std::cout << "[info] Loaded electric field map (181x141 nodes).\n";

    // -------------------------------------------------------------------------
    // 3. 加载 P+ 权势 (Shockley-Ramo)
    // -------------------------------------------------------------------------
    ComponentGrid wfield;
    wfield.SetCylindricalCoordinates();
    wfield.SetMesh(181, 1, 141,
                   0.0, 1.8, 0.0, TWO_PI, -0.2, 1.2);
    if (!wfield.LoadWeightingField("../output/garfield/wpot_Pplus.txt",
                                   "xyz", true)) {
        std::cerr << "[fatal] LoadWeightingField failed\n";
        return 1;
    }
    wfield.SetMedium(&ge);
    std::cout << "[info] Loaded P+ weighting field.\n";

    // -------------------------------------------------------------------------
    // 4. Sensor + 时间窗口
    // -------------------------------------------------------------------------
    Sensor sensor;
    sensor.AddComponent(&efield);
    sensor.AddElectrode(&wfield, "Pplus");

    const double t_start  = 0.0;     // ns
    const double t_step   = 1.0;     // ns
    const int    n_bins   = 1000;    // -> 总长 1000 ns = 1 μs
    sensor.SetTimeWindow(t_start, t_step, n_bins);
    std::cout << "[info] Sensor time window: 0..1000 ns, dt = 1 ns.\n";

    // -------------------------------------------------------------------------
    // 5. AvalancheMC：漂移 e/h 对
    //    位置：r=5mm, z=5mm  (cartesian: x=0.5cm, y=0, z=0.5cm)
    //    电荷对数：122 keV / 2.96 eV/pair ≈ 41216
    // -------------------------------------------------------------------------
    AvalancheMC mc;
    mc.SetSensor(&sensor);
    mc.EnableSignalCalculation(true);
    mc.SetTimeSteps(0.05);           // 50 ps step
    mc.EnableRKFSteps(true);         // Runge-Kutta-Fehlberg 自适应步长

    const double x0 = 0.5;           // cm
    const double y0 = 0.0;
    const double z0 = 0.5;           // cm
    const double t0 = 0.0;
    const int    n_pairs = 41216;

    // 让一个虚拟电子/空穴的信号代表 n_pairs 对
    mc.SetElectronSignalScalingFactor(static_cast<double>(n_pairs));
    mc.SetHoleSignalScalingFactor(static_cast<double>(n_pairs));

    std::cout << "[info] Drifting e/h pair at (x=" << x0
              << ", y=" << y0 << ", z=" << z0 << ") cm, scale = "
              << n_pairs << " pairs.\n";

    bool ok_e = mc.DriftElectron(x0, y0, z0, t0);
    if (!ok_e) std::cerr << "[warn] DriftElectron returned false\n";
    bool ok_h = mc.DriftHole(x0, y0, z0, t0);
    if (!ok_h) std::cerr << "[warn] DriftHole returned false\n";

    // -------------------------------------------------------------------------
    // 6. 提取波形并写 CSV
    //    Sensor::GetSignal 返回 induced current（每 bin），单位 fC/ns 或类似
    //    累积求和得到 Q(t)，并把 Garfield 内部电荷归一化到电子数
    // -------------------------------------------------------------------------
    std::ofstream out("signal_garfield.csv");
    out << "time_ns,current,Q_cum\n";
    double q_cum = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        const double t = t_start + (i + 0.5) * t_step;
        const double s = sensor.GetSignal("Pplus", i);   // induced current
        q_cum += s * t_step;
        out << std::fixed << std::setprecision(3) << t << ","
            << std::scientific << std::setprecision(6) << s << ","
            << q_cum << "\n";
    }
    out.close();
    std::cout << "[info] Wrote signal_garfield.csv\n";
    std::cout << "[result] Q_final ≈ " << q_cum
              << "  (expected ~" << n_pairs << " e-)\n";
    return 0;
}
