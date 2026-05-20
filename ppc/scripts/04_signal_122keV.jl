# =============================================================================
# 任务 04：在 V_op = 1000 V 下，仿真 122 keV 全沉积事例在不同位置的信号
#
# 物理设定：
#   - 122 keV 来自 57Co 主峰，在 Ge 中以光电效应为主，假定单点全沉积
#   - 122 keV → 122000 eV / 2.96 (eV/对) ≈ 41200 个电子-空穴对
#   - 信号是 P+ 点电极 (contact id=1) 上的诱导电荷（Shockley-Ramo）
#
# 流程：
#   1. 加载几何，覆盖 N+ 电位为 1000 V
#   2. 解电势 + 电场
#   3. 解 P+ 电极的权势 (weighting potential)
#   4. 在晶体内布置一系列 (r, z) 测试点
#   5. 每个点放一个 122 keV 单点事件，调 simulate!() 漂移电荷 + 计算信号
#   6. 保存波形 CSV + 画图
#
# 输出：output/signals/*.csv，output/plots/signals_*.png
# =============================================================================

using SolidStateDetectors
using SolidStateDetectors: CartesianPoint
using Unitful
using DelimitedFiles
using Plots
gr()

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG       = joinpath(PROJECT_ROOT, "config", "ppc.yaml")
const OUTPUT_DIR   = joinpath(PROJECT_ROOT, "output")
const SIGNAL_DIR   = joinpath(OUTPUT_DIR, "signals")
const PLOTS_DIR    = joinpath(OUTPUT_DIR, "plots")
isdir(SIGNAL_DIR) || mkdir(SIGNAL_DIR)
isdir(PLOTS_DIR)  || mkdir(PLOTS_DIR)

const V_OP   = 1000.0u"V"
const E_GAMMA = 122.0u"keV"
const DT      = 1.0u"ns"           # 仿真时间步长
const READOUT_CID = 1               # P+ 点电极

# ----- 1. 加载并配置探测器 -----
@info "加载探测器配置 + 覆盖工作电压"
sim = Simulation{Float32}(CONFIG)
sim.detector = SolidStateDetectors.SolidStateDetector(sim.detector;
    contact_id = 2, contact_potential = V_OP)

# ----- 2. 电势 / 电场（信号仿真依赖电场来定义漂移速度方向） -----
@info "求解电势..."
calculate_electric_potential!(sim, refinement_limits = [0.2, 0.1, 0.05, 0.01])
@info "求解电场..."
calculate_electric_field!(sim)

# ----- 3. 权势：Shockley-Ramo 定理用，告诉电荷在每个位置贡献多少诱导信号 -----
@info "求解 P+ 电极的权势（读出电极 contact id=$READOUT_CID）..."
calculate_weighting_potential!(sim, READOUT_CID,
    refinement_limits = [0.2, 0.1, 0.05, 0.01])

# 同时把 N+ 权势也算了，方便后续看 N+ 信号（对 PPC 不是必须的）
@info "求解 N+ 权势（contact id=2）..."
calculate_weighting_potential!(sim, 2,
    refinement_limits = [0.2, 0.1, 0.05, 0.01])

# ----- 4. 测试位置：(r, z) in mm -----
# 分三组覆盖晶体：
#   (a) 中轴上不同深度  -> 看深度对脉冲形状的影响
#   (b) 中深度不同半径  -> 看径向位置的影响
#   (c) 四个角落        -> 极端漂移路径
positions_mm = [
    # (label, r, z)
    ("axis_z1",  0.0,  1.0),
    ("axis_z3",  0.0,  3.0),
    ("axis_z5",  0.0,  5.0),
    ("axis_z7",  0.0,  7.0),
    ("axis_z9",  0.0,  9.0),
    ("mid_r3",   3.0,  5.0),
    ("mid_r6",   6.0,  5.0),
    ("mid_r9",   9.0,  5.0),
    ("mid_r12", 12.0,  5.0),
    ("mid_r14", 14.0,  5.0),
    ("corner_inner_bottom", 1.0, 1.0),
    ("corner_outer_bottom",14.0, 1.0),
    ("corner_inner_top",    1.0, 9.0),
    ("corner_outer_top",   14.0, 9.0),
]

# ----- 5. 逐个仿真 -----
@info "开始仿真 $(length(positions_mm)) 个事件 @ $E_GAMMA"
results = []  # 每个元素: (label, r_mm, z_mm, time_us, signal_e)

for (label, r_mm, z_mm) in positions_mm
    # CartesianPoint 内部单位是 米；这里直接构造 m 数值
    # 中轴对称：x = r, y = 0, z = z
    pt = CartesianPoint{Float32}(Float32(r_mm/1000), 0.0f0, Float32(z_mm/1000))
    evt = Event(pt, E_GAMMA)

    # simulate! 先漂移电荷再算波形
    simulate!(evt, sim; Δt = DT, max_nsteps = 2000, verbose = false)

    wf = evt.waveforms[READOUT_CID]
    t_us  = Float64.(ustrip.(uconvert.(u"µs", collect(wf.time))))
    sig_e = Float64.(ustrip.(wf.signal))   # 默认 signal_unit = e_au (电荷数)

    push!(results, (label, r_mm, z_mm, t_us, sig_e))
    println("  $label  (r=$(r_mm) mm, z=$(z_mm) mm): drift t = $(round(last(t_us), digits=3)) µs, " *
            "Q_final = $(round(last(sig_e), digits=1)) e-")
end

# ----- 6. 保存 CSV：每个事件一行 + 单独保存波形数据 -----
# 一份总表
open(joinpath(SIGNAL_DIR, "events_summary.csv"), "w") do f
    println(f, "label,r_mm,z_mm,drift_time_us,Q_final_e")
    for (label, r, z, t, s) in results
        println(f, "$label,$r,$z,$(last(t)),$(last(s))")
    end
end

# 每个事件单独的 (时间, 信号) CSV
for (label, r, z, t, s) in results
    open(joinpath(SIGNAL_DIR, "wf_$label.csv"), "w") do f
        println(f, "time_us,signal_e")
        for i in eachindex(t)
            println(f, "$(t[i]),$(s[i])")
        end
    end
end

# ----- 7. 画图：三组分别画 -----
function plot_group(group_results, title_text, fname)
    p = plot(xlabel = "Time (µs)", ylabel = "Induced charge on P+ (e⁻)",
             title = title_text, size = (800, 500), dpi = 150, legend = :bottomright)
    cgrad_colors = palette(:viridis, length(group_results))
    for (i, (label, r, z, t, s)) in enumerate(group_results)
        plot!(p, t, s, lw = 2, color = cgrad_colors[i],
              label = "$label (r=$r, z=$z)")
    end
    savefig(p, joinpath(PLOTS_DIR, fname))
    @info "保存波形图" file=fname
end

axis_group   = filter(r -> startswith(r[1], "axis"),   results)
mid_group    = filter(r -> startswith(r[1], "mid"),    results)
corner_group = filter(r -> startswith(r[1], "corner"), results)

plot_group(axis_group,   "122 keV signals along central axis (r=0)", "signals_axis.png")
plot_group(mid_group,    "122 keV signals at z=5, varying r",        "signals_radial.png")
plot_group(corner_group, "122 keV signals at corner positions",       "signals_corner.png")

# 全部画在一张图上
p_all = plot(xlabel = "Time (µs)", ylabel = "Induced charge on P+ (e⁻)",
    title = "122 keV signals — all positions", size = (900, 550), dpi = 150,
    legend = false)
cgrad_colors = palette(:turbo, length(results))
for (i, (label, r, z, t, s)) in enumerate(results)
    plot!(p_all, t, s, lw = 1.5, color = cgrad_colors[i], alpha = 0.85)
end
savefig(p_all, joinpath(PLOTS_DIR, "signals_all.png"))
@info "保存波形图" file="signals_all.png"

println("\n完成。波形 CSV: $SIGNAL_DIR")
println("波形图: $PLOTS_DIR/signals_*.png")
