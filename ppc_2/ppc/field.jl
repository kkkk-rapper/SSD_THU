# field.jl — 1000V 下计算电场，并画「电势 / 电场」图 (只含晶体部分)
# 用法: julia -t auto --project=<SSD根目录> ppc/field.jl
#
# 对 config.yaml(有划痕) 与 config_noscratch.yaml(无划痕) 各做:
#   · calculate_electric_potential! (CUDA 加速) + calculate_electric_field!
#   · 电势图、电场强度图 —— r-z 截面 @ φ=90° (穿过划痕最深处)，
#     裁剪到 r∈[0,15]mm、z∈[0,10]mm —— 即「只含晶体部分」
# 结果写 output/<tag>/{potential,field}.png

ENV["GKSwstype"] = "100"                       # GR 离屏渲染

using SolidStateDetectors
using Unitful
using LinearAlgebra
using Plots
using Plots.PlotMeasures              # 提供 mm 等边距单位 (Plots 默认未导出)
using CUDA
gr()

SCRIPT_DIR = @__DIR__
ROOT_DIR   = dirname(SCRIPT_DIR)
T = Float32
DEVICE = CUDA.functional() ? CuArray : Array

# 先跑便宜的无划痕基线 (探路)，再跑有划痕
CONFIGS = [("config_noscratch.yaml", "noscratch"), ("config.yaml", "scratch")]

# 晶体范围: r 0→15mm, z 0→10mm —— 绘图裁剪到「只含晶体」(坐标轴用 mm)
CROP = (xlims = (0u"mm", 15u"mm"), ylims = (0u"mm", 10u"mm"))

# 求解参数: max_tick_distance 用元组 (Δr, Δφ, Δz) —— 细化 r,z，φ 上限 30°
# (★ 切勿用标量 max_tick_distance：那会按弧长强制 φ 加密、令 3D 划痕模型爆炸)
SOLVER_KW = (
    convergence_limit = 1e-5,
    refinement_limits = [0.2, 0.1, 0.05],
    max_tick_distance = (0.4u"mm", 30u"°", 0.4u"mm"),
    device_array_type = DEVICE,
    verbose           = true,
)

println("="^60)
println("PPC HPGe — 1000V 电场计算 + 电势/电场图")
println("  线程 = ", Threads.nthreads(), "   设备 = ", DEVICE === CuArray ? "GPU" : "CPU")
println("="^60)

for (cfg, tag) in CONFIGS
    println("\n", "█"^58)
    println("█ ", cfg, "   (", tag, ")")
    println("█"^58); flush(stdout)
    outdir = joinpath(SCRIPT_DIR, "output", tag)
    mkpath(outdir)

    sim = Simulation{T}(joinpath(SCRIPT_DIR, cfg))

    println("[1] 求解电势 (CUDA)..."); flush(stdout)
    calculate_electric_potential!(sim; SOLVER_KW...)

    println("[2] 求解电场..."); flush(stdout)
    calculate_electric_field!(sim)

    Emax = maximum(norm.(sim.electric_field))                 # V/m
    println("    完全耗尽? ", is_depleted(sim.point_types))
    println("    最大场强 |E|max ≈ ", round(Emax / 1e5, digits = 2), " kV/cm")
    flush(stdout)

    # --- 电势图: r-z 截面 @ φ=90°，坐标轴 mm，裁剪到只含晶体 ---
    plot(sim.electric_potential, φ = 90u"°";
         title = "Electric Potential — $tag  (bias 1000 V)",
         xunit = u"mm", yunit = u"mm",
         size = (880, 620), right_margin = 12mm, CROP...)
    savefig(joinpath(outdir, "potential.png"))
    println("    → output/$tag/potential.png")

    # --- 电场强度图 ---
    # clims 截断到 6e5 V/m (=6 kV/cm)：接触电极边缘是近奇异尖峰、会吃掉整个
    # 色标，截断后才能看清晶体体内的场分布。
    plot(sim.electric_field, φ = 90u"°";
         title = "Electric Field |E| — $tag  (1000 V, color clipped @ 6e5 V/m)",
         clims = (0.0, 6e5),
         xunit = u"mm", yunit = u"mm",
         size = (880, 620), right_margin = 12mm, CROP...)
    savefig(joinpath(outdir, "field.png"))
    println("    → output/$tag/field.png")
    flush(stdout)

    GC.gc()
    CUDA.functional() && CUDA.reclaim()           # 释放显存给下个配置
end

println("\n", "="^60)
println("完成 — 4 张图: output/{scratch,noscratch}/{potential,field}.png")
println("="^60)
