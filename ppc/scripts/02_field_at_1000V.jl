# =============================================================================
# 任务 02：在工作电压 V_op = 1000 V 下计算电势和电场
#
# 流程：
#   1. 加载 YAML（里面 contact 2 默认是 2000 V，仅用于任务 01）
#   2. 用 SolidStateDetector(...; contact_id=2, contact_potential=1000V) 覆盖电位
#   3. 求解 Poisson 方程 -> 电势
#   4. 由电势数值微分 -> 电场矢量 (Er, Ephi, Ez)
#   5. 把柱坐标 (phi=1 切片) 的 2D 电势/电场矩阵写成 CSV
#
# 输出：output/potential_V_rz.csv, Efield_*_rz.csv, field_summary_1000V.txt
# =============================================================================

using SolidStateDetectors
using Unitful
using DelimitedFiles

# ----- 路径定义 -----
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG       = joinpath(PROJECT_ROOT, "config", "ppc.yaml")
const OUTPUT_DIR   = joinpath(PROJECT_ROOT, "output")
isdir(OUTPUT_DIR) || mkdir(OUTPUT_DIR)

# ----- 工作电压：施加在 N+ 包裹电极 (contact id=2) -----
const V_OP = 1000.0u"V"

# ----- 1. 从 YAML 加载探测器 -----
@info "加载探测器配置" CONFIG
sim = Simulation{Float32}(CONFIG)

# ----- 2. 把 contact 2 的电位覆盖为 V_OP（不修改 YAML 文件本身） -----
# SolidStateDetector 的关键字构造器，只改 contact 2 的 potential，其他不变
sim.detector = SolidStateDetectors.SolidStateDetector(sim.detector;
    contact_id = 2, contact_potential = V_OP)
@info "覆盖后的电极电位" map(c -> (c.id, c.potential), sim.detector.contacts)

# ----- 3. 求解 Poisson 方程，得到电势 -----
@info "求解 Poisson 方程..."
calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1, 0.05, 0.01])

# ----- 4. 检查 V_op = 1000 V 时是否完全耗尽 -----
# 任务 01 算出 V_dep ≈ 120 V，所以 1000 V 应当远超耗尽，结果应为 true
depleted = is_depleted(sim.point_types)
@info "完全耗尽检查（V_op = $V_OP）" depleted
if !depleted
    @warn "在 $V_OP 下未完全耗尽 -> 电场中将存在未耗尽（无场）区域"
end

# ----- 5. 由电势计算电场矢量 -----
@info "计算电场..."
calculate_electric_field!(sim)

# ----- 6. 提取轴对称切片（phi index = 1） -----
# 注意：SSD 内部长度单位是 米；这里把网格轴转换为 mm 输出更直观
pot  = sim.electric_potential
efld = sim.electric_field

r_m = collect(pot.grid.axes[1].ticks)   # r 轴（米）
z_m = collect(pot.grid.axes[3].ticks)   # z 轴（米）
r_mm = Float64.(r_m) .* 1000.0           # 转 mm
z_mm = Float64.(z_m) .* 1000.0           # 转 mm

# 电势数据：pot.data[r, phi, z]，标量，单位 V
pot_2d = Float64.(pot.data[:, 1, :])

# 电场数据：efld.data[r, phi, z] 每个元素是 SVector{3} = (Er, Ephi, Ez)
# 单位 V/m
Nr, Nz = length(r_mm), length(z_mm)
Er = zeros(Float64, Nr, Nz)             # 径向分量
Ez = zeros(Float64, Nr, Nz)             # 轴向分量
Em = zeros(Float64, Nr, Nz)             # 矢量模长 |E|
for i in 1:Nr, k in 1:Nz
    v = efld.data[i, 1, k]
    Er[i, k] = v[1]
    Ez[i, k] = v[3]
    Em[i, k] = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
end

# ----- 7. 写 CSV：矩阵的 行 = r 轴方向，列 = z 轴方向 -----
writedlm(joinpath(OUTPUT_DIR, "r_axis_mm.csv"),               r_mm,   ',')
writedlm(joinpath(OUTPUT_DIR, "z_axis_mm.csv"),               z_mm,   ',')
writedlm(joinpath(OUTPUT_DIR, "potential_V_rz.csv"),          pot_2d, ',')
writedlm(joinpath(OUTPUT_DIR, "Efield_Er_Vpm_rz.csv"),        Er,     ',')
writedlm(joinpath(OUTPUT_DIR, "Efield_Ez_Vpm_rz.csv"),        Ez,     ',')
writedlm(joinpath(OUTPUT_DIR, "Efield_Emag_Vpm_rz.csv"),      Em,     ',')

# ----- 8. 限制在晶体内的统计量（剔除边界外、电极内部的网格点） -----
# 晶体范围：r <= 15 mm，0 <= z <= 10 mm
inside = falses(Nr, Nz)
for i in 1:Nr, k in 1:Nz
    if r_mm[i] <= 15.0 && 0.0 <= z_mm[k] <= 10.0
        inside[i, k] = true
    end
end
Emag_in = Em[inside]
pot_in  = pot_2d[inside]

# ----- 9. 输出文本摘要 -----
open(joinpath(OUTPUT_DIR, "field_summary_1000V.txt"), "w") do f
    println(f, "PPC 电场/电势计算结果 @ V_op = $V_OP (施加在 contact id=2, N+ 包裹)")
    println(f, "=====================================================================")
    println(f, "V_op 下是否完全耗尽: $depleted")
    println(f, "")
    println(f, "网格 r (mm): N=$Nr, $(round(first(r_mm), digits=3)) .. $(round(last(r_mm), digits=3))")
    println(f, "网格 z (mm): N=$Nz, $(round(first(z_mm), digits=3)) .. $(round(last(z_mm), digits=3))")
    println(f, "")
    println(f, "晶体内统计 (r<=15, 0<=z<=10):")
    println(f, "  电势  min/max  : $(minimum(pot_in)) / $(maximum(pot_in)) V")
    println(f, "  |E|   min/max  : $(minimum(Emag_in)) / $(maximum(Emag_in)) V/m")
    println(f, "  |E|   median   : $(sort(Emag_in)[length(Emag_in)÷2]) V/m")
    println(f, "")
    println(f, "已写入文件（CSV 矩阵：行 = r 轴, 列 = z 轴）：")
    println(f, "  r_axis_mm.csv, z_axis_mm.csv")
    println(f, "  potential_V_rz.csv")
    println(f, "  Efield_Er_Vpm_rz.csv, Efield_Ez_Vpm_rz.csv, Efield_Emag_Vpm_rz.csv")
end
println("完成。输出目录: $OUTPUT_DIR")
