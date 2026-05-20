# =============================================================================
# 任务 06：把 SSD 算出的电场 + 电势 + 权势导出为 Garfield++ 的
#         ComponentGrid 格式（柱坐标、cm、V/cm、V）
#
# 用途：在 Garfield++ 里跑 Heed/AvalancheMicroscopic，跟 SSD 做信号对照
#
# Garfield++ ComponentGrid 在柱坐标 + "xyz" 格式下每行：
#     r  phi  z  Er  Et  Ez  V
# 我们是轴对称 -> ny = 1, phi = 0；Et 始终为 0
#
# 单位约定：
#   - 坐标 cm
#   - 电场 V/cm
#   - 电势 V
#
# 输出：output/garfield/
#   efield_1000V.txt        E-field + 静电势
#   wpot_Pplus.txt          P+ (contact 1) 的权势
#   load_ssd_field.C        ROOT/Garfield++ 加载示例 (C++)
# =============================================================================

using SolidStateDetectors
using SolidStateDetectors: interpolated_scalarfield, interpolated_vectorfield,
                           CylindricalPoint, CartesianPoint
using Unitful
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG       = joinpath(PROJECT_ROOT, "config", "ppc.yaml")
const OUT_DIR      = joinpath(PROJECT_ROOT, "output", "garfield")
isdir(OUT_DIR) || mkdir(OUT_DIR)

const V_OP = 1000.0u"V"

# ----- 1. 跑 SSD，得到电势、电场、P+ 权势 -----
@info "加载探测器 + 设置 V_op = $V_OP"
sim = Simulation{Float32}(CONFIG)
sim.detector = SolidStateDetectors.SolidStateDetector(sim.detector;
    contact_id = 2, contact_potential = V_OP)

@info "求解电势..."
calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1, 0.05, 0.01])
@info "求解电场..."
calculate_electric_field!(sim)
@info "求解 P+ 权势..."
calculate_weighting_potential!(sim, 1, refinement_limits = [0.2, 0.1, 0.05, 0.01])

# 用 SSD 自带的插值器在任意 (r, z) 评估
pot_itp = interpolated_scalarfield(sim.electric_potential)
efld_itp = interpolated_vectorfield(sim.electric_field.data, sim.electric_field.grid)
wpot_itp = interpolated_scalarfield(sim.weighting_potentials[1])

# ----- 2. 定义导出用的均匀网格 -----
# 内部单位是 米；ComponentGrid 需要 cm。Garfield++ 网格是均匀的。
# r: 0 .. 18 mm，步长 0.1 mm  -> 181 点
# z: -2 .. 12 mm，步长 0.1 mm -> 141 点
const DR_mm = 0.1
const DZ_mm = 0.1
const RMIN_mm, RMAX_mm = 0.0, 18.0
const ZMIN_mm, ZMAX_mm = -2.0, 12.0
r_mm = collect(RMIN_mm:DR_mm:RMAX_mm)
z_mm = collect(ZMIN_mm:DZ_mm:ZMAX_mm)
Nr = length(r_mm)
Nz = length(z_mm)
@info "导出网格大小" Nr=Nr Nz=Nz total=(Nr*Nz)

# 网格点坐标（米，给插值器用）
r_m = r_mm .* 1e-3
z_m = z_mm .* 1e-3

# 转 cm，给 Garfield++ 用
r_cm = r_mm .* 0.1
z_cm = z_mm .* 0.1

# ----- 3. 评估场 + 电势在每个网格点 -----
# CylindricalPoint{T}(r, phi, z) in 米
@info "评估场..."
Er_Vpm  = zeros(Float64, Nr, Nz)
Ez_Vpm  = zeros(Float64, Nr, Nz)
V_V     = zeros(Float64, Nr, Nz)
Wpot    = zeros(Float64, Nr, Nz)

for i in 1:Nr, k in 1:Nz
    pt = CylindricalPoint{Float32}(Float32(r_m[i]), 0f0, Float32(z_m[k]))
    Er_Vpm[i, k] = efld_itp(r_m[i], 0.0, z_m[k])[1]
    Ez_Vpm[i, k] = efld_itp(r_m[i], 0.0, z_m[k])[3]
    V_V[i, k]    = pot_itp(r_m[i], 0.0, z_m[k])
    Wpot[i, k]   = wpot_itp(r_m[i], 0.0, z_m[k])
end

# V/m -> V/cm
Er_Vpcm = Er_Vpm ./ 100.0
Ez_Vpcm = Ez_Vpm ./ 100.0

# ----- 4. 写 Garfield++ ComponentGrid 文件 -----
# 用 "xyz" 格式（柱坐标下 x=r, y=phi, z=z），每行：
#   r_cm  phi(=0)  z_cm  Er  Et(=0)  Ez  V

function write_field_file(fname::String, header_lines::Vector{String},
                         field_columns::Function)
    open(fname, "w") do f
        for line in header_lines
            println(f, "# " * line)
        end
        # 头部参数（Garfield++ 用它们 sanity check，可选）
        @printf(f, "# XMIN = %.6f, XMAX = %.6f, NX = %d\n",
                r_cm[1], r_cm[end], Nr)
        @printf(f, "# YMIN = %.6f, YMAX = %.6f, NY = %d\n", 0.0, 0.0, 1)
        @printf(f, "# ZMIN = %.6f, ZMAX = %.6f, NZ = %d\n",
                z_cm[1], z_cm[end], Nz)
        # 数据行：(i, j=1, k) 嵌套，注意 Garfield++ 期望 i 慢、k 快
        # 实际查看源码 LoadElectricField 的循环 -> for k for j for i
        # 但 SaveElectricField 写的是 for i for j for k -> 写时按 i,j,k 顺序
        # 保险起见我们按 i, j, k 三层循环（j 只有 1）
        for i in 1:Nr
            for k in 1:Nz
                line = field_columns(i, k)
                println(f, line)
            end
        end
    end
end

@info "写 efield_1000V.txt..."
efield_file = joinpath(OUT_DIR, "efield_1000V.txt")
write_field_file(efield_file,
    ["SSD-computed PPC HPGe field map, V_op = 1000 V on N+ contact (id=2)",
     "Cylindrical coordinates: r, phi=0, z   (axisymmetric)",
     "Units: cm for coords, V/cm for fields, V for potential",
     "Format string for Garfield++ LoadElectricField: \"xyz\"",
     "withPotential = true, withFlag = false",
     "Columns: r[cm]  phi[rad]  z[cm]  Er[V/cm]  Et[V/cm]  Ez[V/cm]  V[V]"],
    (i, k) -> @sprintf("%.5f  %.5f  %.5f  %+.6e  %+.6e  %+.6e  %+.6e",
        r_cm[i], 0.0, z_cm[k],
        Er_Vpcm[i, k], 0.0, Ez_Vpcm[i, k], V_V[i, k]))

@info "写 wpot_Pplus.txt..."
wpot_file = joinpath(OUT_DIR, "wpot_Pplus.txt")
# 权势的 "电场" 是其负梯度。这里我们写权势本身 + 零电场。
# Garfield++ 加载权势时只需要电势数据；它会自己求负梯度得到加权场。
# 但 ComponentGrid 的 LoadWeightingField 在 withPotential=true 时仍要求场列。
# 简单做法：手动数值求负梯度填进去。这里先填 0，让 Garfield++ 内部从权势求导。
# 实际上 LoadWeightingField 必须要场（Garfield++ 不内插得到），所以我们要自己算 grad(wpot)。
function central_diff_r(A, dr)
    g = zeros(size(A))
    for i in 2:size(A,1)-1, k in axes(A, 2)
        g[i, k] = (A[i+1, k] - A[i-1, k]) / (2 * dr)
    end
    g[1, :]   .= (A[2, :]   .- A[1, :])     ./ dr
    g[end, :] .= (A[end, :] .- A[end-1, :]) ./ dr
    return g
end
function central_diff_z(A, dz)
    g = zeros(size(A))
    for i in axes(A, 1), k in 2:size(A,2)-1
        g[i, k] = (A[i, k+1] - A[i, k-1]) / (2 * dz)
    end
    g[:, 1]   .= (A[:, 2]   .- A[:, 1])     ./ dz
    g[:, end] .= (A[:, end] .- A[:, end-1]) ./ dz
    return g
end
# wpot 是无量纲（0~1），梯度单位 1/cm，对应 "权场" V/cm 时按惯例直接用（电极电压归一化）
Wr_pcm = -central_diff_r(Wpot, DR_mm * 0.1)  # 注意 dr 用 cm
Wz_pcm = -central_diff_z(Wpot, DZ_mm * 0.1)

write_field_file(wpot_file,
    ["SSD-computed P+ contact (id=1) weighting potential & weighting field",
     "Cylindrical coordinates: r, phi=0, z   (axisymmetric)",
     "Units: cm for coords, 1/cm for weighting field, dimensionless for potential",
     "Format string for Garfield++ LoadWeightingField: \"xyz\"",
     "withPotential = true",
     "Columns: r[cm]  phi[rad]  z[cm]  Wr[1/cm]  Wt[1/cm]  Wz[1/cm]  Wpot[dimensionless]"],
    (i, k) -> @sprintf("%.5f  %.5f  %.5f  %+.6e  %+.6e  %+.6e  %+.6e",
        r_cm[i], 0.0, z_cm[k],
        Wr_pcm[i, k], 0.0, Wz_pcm[i, k], Wpot[i, k]))

# ----- 5. 写一个 ROOT 宏作为 Garfield++ 加载示例 -----
cpp_example = """
// =============================================================================
// load_ssd_field.C — ROOT 宏，演示如何在 Garfield++ 里加载 SSD 导出的场图
//
// 编译运行：
//     cd /root/ssd_projects/ppc/output/garfield
//     root -l load_ssd_field.C
//
// 注意 Garfield++ 必须先 source /root/gpp/install/share/Garfield/setupGarfield.sh
// 或类似设置好 LD_LIBRARY_PATH 和 GARFIELD_INSTALL
// =============================================================================

#include <iostream>
#include "Garfield/ComponentGrid.hh"
#include "Garfield/MediumSilicon.hh"   // 没有现成的 MediumGermanium；用 Si 占位或自定义
#include "Garfield/Sensor.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

void load_ssd_field() {
    // 1. 电场 + 静电势
    ComponentGrid efield;
    efield.SetCylindricalCoordinates();
    // 网格尺寸要跟导出文件的头部一致
    efield.SetMesh(
        $(Nr), 1, $(Nz),                            // nr, nphi, nz
        $(r_cm[1]),  $(r_cm[end]),                  // r 范围 (cm)
        0.0,         0.0,                            // phi (单平面 axisym)
        $(z_cm[1]),  $(z_cm[end])                   // z 范围 (cm)
    );
    if (!efield.LoadElectricField("efield_1000V.txt", "xyz", true, false))
        std::cerr << "Failed to load electric field\\n";

    // 2. P+ 权势（用于 Shockley-Ramo 信号仿真）
    ComponentGrid wfield;
    wfield.SetCylindricalCoordinates();
    wfield.SetMesh($(Nr), 1, $(Nz),
                   $(r_cm[1]),  $(r_cm[end]), 0.0, 0.0,
                   $(z_cm[1]),  $(z_cm[end]));
    if (!wfield.LoadWeightingField("wpot_Pplus.txt", "xyz", true))
        std::cerr << "Failed to load weighting field\\n";

    // 3. 装入 Sensor（你后面用 AvalancheMicroscopic 等需要 sensor）
    Sensor sensor;
    sensor.AddComponent(&efield);
    sensor.AddElectrode(&wfield, "Pplus");

    // 4. 可视化检查（可选）
    ViewField view;
    view.SetComponent(&efield);
    view.SetPlane(0, -1, 0, 0, 0, 0);   // r-z 平面
    view.SetArea($(r_cm[1]), $(z_cm[1]),
                 $(r_cm[end]), $(z_cm[end]));
    view.PlotContour("v");
}
"""

cpp_file = joinpath(OUT_DIR, "load_ssd_field.C")
write(cpp_file, cpp_example)
@info "写 ROOT 宏示例" cpp_file

println("\n完成。导出文件:")
println("  $efield_file")
println("  $wpot_file")
println("  $cpp_file")
println("\n下一步:")
println("  cd $OUT_DIR && root -l load_ssd_field.C")
