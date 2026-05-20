# =============================================================================
# 任务 01：计算 PPC 探测器的耗尽电压 V_dep
#
# 原理：在 YAML 里给一个足够高的初始偏压（这里是 2000 V），求解 Poisson 方程，
#       确认探测器在该偏压下已经完全耗尽，然后用 SSD 内置的
#       estimate_depletion_voltage() 通过二分法搜索 V_dep。
#
# 输出：output/depletion_voltage.txt
# =============================================================================

using SolidStateDetectors
using Unitful

# ----- 路径定义（脚本可独立运行，路径自动相对于 scripts/） -----
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG       = joinpath(PROJECT_ROOT, "config", "ppc.yaml")
const OUTPUT_DIR   = joinpath(PROJECT_ROOT, "output")
isdir(OUTPUT_DIR) || mkdir(OUTPUT_DIR)

# ----- 1. 从 YAML 加载探测器，使用 Float32 精度 -----
@info "加载探测器配置" CONFIG
sim = Simulation{Float32}(CONFIG)

# ----- 2. 求解 Poisson 方程，得到电势分布 -----
# refinement_limits：网格自适应加密阈值；数值越小越精细，但收敛慢
# 这里 4 级加密足够；2D 柱坐标问题计算量很小
@info "在初始过偏压下求解 Poisson 方程 (contact 2 = $(sim.detector.contacts[2].potential))"
calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1, 0.05, 0.01])

# ----- 3. 检查是否完全耗尽 -----
# is_depleted 检查 point_types 数组中是否所有体内点都是 "bulk + depleted" 状态
depleted = is_depleted(sim.point_types)
@info "完全耗尽检查（初始偏压下）" depleted
@assert depleted "在初始偏压下探测器未完全耗尽，请提高 YAML 里 contact 2 的电位"

# ----- 4. 用二分法估计耗尽电压 -----
# 算法：基于权势 + 杂质势的线性叠加，扫描 (Umin, Umax) 找出 V_dep
# 默认 Umin=0, Umax=initial_bias，tolerance=0.1 V
V_dep = estimate_depletion_voltage(sim; verbose = true)
@info "估计的耗尽电压" V_dep

# ----- 5. 把结果写入文本文件 -----
open(joinpath(OUTPUT_DIR, "depletion_voltage.txt"), "w") do f
    println(f, "PPC 耗尽电压计算结果 (SolidStateDetectors)")
    println(f, "===========================================")
    println(f, "晶体        : 直径 30 mm × 高 10 mm, p 型 HPGe")
    println(f, "杂质浓度    : -2.8e6/mm^3 (z=0, N+ 面)  →  -3.0e6/mm^3 (z=10, P+ 面)")
    println(f, "施压电极    : id=2 (N+ 包裹)，初始过偏压 = $(sim.detector.contacts[2].potential) V")
    println(f, "")
    println(f, "V_dep = $V_dep")
end
println("V_dep = $V_dep")
println("结果已写入: ", joinpath(OUTPUT_DIR, "depletion_voltage.txt"))
