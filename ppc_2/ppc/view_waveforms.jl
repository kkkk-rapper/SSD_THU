# view_waveforms.jl — 按测得能量段调看已保存的信号波形
# 用法: julia --project=<SSD根目录> ppc/view_waveforms.jl <tag> <Emin_keV> <Emax_keV> [最多条数]
#   例: julia --project=. ppc/view_waveforms.jl noscratch 120 124
#       julia --project=. ppc/view_waveforms.jl noscratch 35 42 50
#
# 读 output/<tag>/ 下 signal.jl 保存的 waveforms.bin + events.txt,
# 画出测得能量落在 [Emin,Emax] 的所有波形 → output/<tag>/view_<Emin>-<Emax>keV.png

ENV["GKSwstype"] = "100"
using Plots
gr()

length(ARGS) >= 3 || error("用法: view_waveforms.jl <tag> <Emin> <Emax> [最多条数]")
TAG  = ARGS[1]
EMIN = parse(Float64, ARGS[2])
EMAX = parse(Float64, ARGS[3])
NMAX = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 200

SCRIPT_DIR = @__DIR__
OUTDIR = joinpath(SCRIPT_DIR, "output", TAG)

# ---- 读 events.txt: 第 6 列 = 测得能量 emeas_keV ----
emeas = Float64[]
for line in eachline(joinpath(OUTDIR, "events.txt"))
    startswith(strip(line), "#") && continue
    isempty(strip(line)) && continue
    push!(emeas, parse(Float64, split(line)[6]))
end

# ---- 读 waveforms.bin: Int64 n, Int64 L, 然后 L×n Float32 ----
N, L, W = open(joinpath(OUTDIR, "waveforms.bin")) do io
    n = read(io, Int64); l = read(io, Int64)
    m = Matrix{Float32}(undef, l, n); read!(io, m)
    n, l, m
end
@assert length(emeas) == N "events.txt 行数与 waveforms.bin 不一致"

# ---- 选出能量段内的事件 ----
sel = findall(e -> EMIN ≤ e ≤ EMAX, emeas)
println("能量段 ", EMIN, "–", EMAX, " keV: 命中 ", length(sel), " 个事件 (共 ", N, ")")
isempty(sel) && error("该能量段没有事件")
length(sel) > NMAX && (sel = sel[1:NMAX])

# ---- 画图 ----
plt = plot(size=(960,620), legend=false,
           title="$TAG — waveforms with measured energy $(EMIN)–$(EMAX) keV  ($(length(sel)) shown)",
           xlabel="time (ns)", ylabel="signal (keV)")
for i in sel
    plot!(plt, 1:L, W[:,i], lw=0.8, alpha=0.4, color=:steelblue)
end
outpng = joinpath(OUTDIR, "view_$(EMIN)-$(EMAX)keV.png")
savefig(plt, outpng)
println("→ ", outpng)
