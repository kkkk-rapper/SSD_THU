# compare_spectra.jl — 有/无划痕能谱对比 (线性坐标)
# 读 output/{noscratch,scratch}/spectrum.txt (signal.jl 生成),
# 用线性 y 轴画两条能谱叠加对比 → output/spectrum_compare_linear.png
# 用法: julia --project=<SSD根目录> ppc/compare_spectra.jl

ENV["GKSwstype"] = "100"
using Plots
gr()

SCRIPT_DIR = @__DIR__

# 读 spectrum.txt: 列 = E_keV(bin中心)  counts_measured  counts_deposited
function read_spec(tag)
    E = Float64[]; cm = Float64[]
    for line in eachline(joinpath(SCRIPT_DIR, "output", tag, "spectrum.txt"))
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        p = split(s)
        push!(E, parse(Float64, p[1])); push!(cm, parse(Float64, p[2]))
    end
    E, cm
end

E, cm_ns = read_spec("noscratch")     # 无划痕
_, cm_s  = read_spec("scratch")       # 有划痕  (bin 与无划痕一致)

# 连续谱区 (E<110 keV) 的峰值 → 决定 zoom 面板的 y 上限
cont_max = maximum(cm_ns[i] for i in eachindex(E) if E[i] < 110)

# 面板 1: 线性, 全谱 (光电峰主导)
p1 = plot(E, cm_ns, seriestype=:steppost, lw=1.5, color=:navy, label="no scratch",
          title="Linear spectrum — full range",
          xlabel="energy (keV)", ylabel="counts / 0.5 keV", legend=:topleft)
plot!(p1, E, cm_s, seriestype=:steppost, lw=1.3, color=:orangered, ls=:dash,
      label="scratch")

# 面板 2: 线性, y 放大到连续谱 (看康普顿区; 光电峰被截顶)
p2 = plot(E, cm_ns, seriestype=:steppost, lw=1.5, color=:navy, label="no scratch",
          title="Linear spectrum — zoom to Compton region",
          xlabel="energy (keV)", ylabel="counts / 0.5 keV",
          ylims=(0, cont_max*1.25), legend=:topright)
plot!(p2, E, cm_s, seriestype=:steppost, lw=1.3, color=:orangered, ls=:dash,
      label="scratch")

plt = plot(p1, p2, layout=(1,2), size=(1340,560))
savefig(plt, joinpath(SCRIPT_DIR, "output", "spectrum_compare_linear.png"))
println("→ output/spectrum_compare_linear.png")

# 简单数字对比
pk_ns = sum(cm_ns[i] for i in eachindex(E) if 120 ≤ E[i] ≤ 124)
pk_s  = sum(cm_s[i]  for i in eachindex(E) if 120 ≤ E[i] ≤ 124)
println("光电峰(120-124keV)计数:  无划痕 ", round(Int,pk_ns), "   有划痕 ", round(Int,pk_s))
