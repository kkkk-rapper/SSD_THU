# =============================================================================
# 比较 Garfield++ 输出和 SSD 输出 (相同位置 r=5mm, z=5mm)
#
# 用法：
#   1. 先确保 signal_garfield.csv 已经由 ./signal_demo 生成
#   2. /root/.juliaup/bin/julia --project=../ compare_with_ssd.jl
# =============================================================================

using DelimitedFiles
using Plots
gr()

const HERE      = @__DIR__
const SSD_DIR   = normpath(joinpath(HERE, "..", "output", "signals"))
const GARF_FILE = joinpath(HERE, "signal_garfield.csv")

# SSD: r=6mm 中深度 z=5mm 是 mid_r6
ssd_data = readdlm(joinpath(SSD_DIR, "wf_mid_r6.csv"), ',', skipstart=1)
t_ssd  = Float64.(ssd_data[:, 1])   # μs
Q_ssd  = Float64.(ssd_data[:, 2])   # e-

garf_data = readdlm(GARF_FILE, ',', skipstart=1)
t_g_ns = Float64.(garf_data[:, 1])  # ns
I_g    = Float64.(garf_data[:, 2])  # current
Q_g    = Float64.(garf_data[:, 3])  # cumulative charge

# 转单位
t_g_us = t_g_ns ./ 1000.0

p = plot(
    xlabel = "Time (μs)",
    ylabel = "Induced charge on P+ (e⁻)",
    title  = "122 keV @ (r≈5-6 mm, z=5 mm): SSD vs Garfield++",
    size = (900, 550), dpi = 150, legend = :bottomright,
)
plot!(p, t_ssd, Q_ssd,
    lw = 2, color = :blue,
    label = "SSD (r=6mm, mid_r6)")
plot!(p, t_g_us, Q_g,
    lw = 2, color = :red, linestyle = :dash,
    label = "Garfield++ (r=5mm)")

savefig(p, joinpath(HERE, "compare_ssd_garfield.png"))
println("Saved: ", joinpath(HERE, "compare_ssd_garfield.png"))
