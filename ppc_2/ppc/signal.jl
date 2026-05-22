# signal.jl — Geant4 hits → SSD 信号波形 + 能谱
# 用法: julia -t auto --project=<SSD根目录> ppc/signal.jl [配置] [事件数]
#   默认: config_noscratch.yaml, 100 事件
#
# 流程: 解电场+权重势(CUDA) → 读 Geant4 hits.csv → 多线程逐事件漂移电荷
#       → p+ 读出电极波形 → 测得能量(波形平台) → 能谱。
# 不含电荷陷阱(SSD 默认理想收集)。
#
# 输出到 output/<tag>/ :
#   events.txt               每事件: id, n_hits, r,z 质心(mm), 沉积能量, 测得能量
#   waveforms.bin            全部波形 (二进制: Int64 n, Int64 L, 然后 L×n Float32)
#   spectrum.png / .txt      测得能量能谱
#   waveforms_by_energy.png  各能量段代表波形
# 后期用 view_waveforms.jl 按能量段调看波形。

ENV["GKSwstype"] = "100"

using SolidStateDetectors
using Unitful
using Plots
using CUDA
using Statistics
gr()

# ---- 参数 ----
CONFIG   = length(ARGS) >= 1 ? ARGS[1] : "config_noscratch.yaml"
N_EVENTS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100
SCRIPT_DIR = @__DIR__
ROOT_DIR   = dirname(SCRIPT_DIR)
T = Float32
DEVICE     = CUDA.functional() ? CuArray : Array
HITS_CSV   = "/root/geant4/sim_project/output/hits.csv"
READOUT_ID = 2                                       # p+ 点电极 = 读出电极
WF_CAP     = 1500                                    # 波形保存最大长度 (ns)
TAG = CONFIG == "config.yaml"           ? "scratch"   :
      CONFIG == "config_noscratch.yaml" ? "noscratch" :
      first(splitext(basename(CONFIG)))
OUTDIR = joinpath(SCRIPT_DIR, "output", TAG)
mkpath(OUTDIR)

println("="^60)
println("信号波形 + 能谱模拟")
println("  配置=", CONFIG, "  (tag=", TAG, ")   事件数=", N_EVENTS)
println("  线程=", Threads.nthreads(), "   设备=", DEVICE === CuArray ? "GPU" : "CPU")
println("="^60)

# ===================================================================
# 1. 建探测器 + 解电场 / 权重势 (CUDA)
# ===================================================================
println("\n[1] 求解电场 + 权重势 (CUDA)..."); flush(stdout)
sim = Simulation{T}(joinpath(SCRIPT_DIR, CONFIG))
simulate!(sim;
    convergence_limit = 1e-5,
    refinement_limits = [0.2, 0.1, 0.05],
    max_tick_distance = (0.4u"mm", 30u"°", 0.4u"mm"),
    device_array_type = DEVICE,
    verbose           = true)
println("    完全耗尽? ", is_depleted(sim.point_types)); flush(stdout)

# ===================================================================
# 2. 读 Geant4 hits.csv 前 N_EVENTS 个事件
# ===================================================================
println("\n[2] 读 hits.csv 前 ", N_EVENTS, " 个事件..."); flush(stdout)
function read_events(path, n)
    evlist = Int[]
    hits = Dict{Int, Vector{NTuple{4,Float64}}}()
    open(path) do io
        readline(io)
        for line in eachline(io)
            p  = split(line, ',')
            ev = parse(Int, p[1])
            if !haskey(hits, ev)
                length(evlist) >= n && break
                push!(evlist, ev); hits[ev] = NTuple{4,Float64}[]
            end
            push!(hits[ev], (parse(Float64,p[2]), parse(Float64,p[3]),
                             parse(Float64,p[4]), parse(Float64,p[5])))
        end
    end
    evlist, hits
end
evlist, hits = read_events(HITS_CSV, N_EVENTS)
println("    读到 ", length(evlist), " 个事件, 共 ",
        sum(length(hits[e]) for e in evlist), " 个能量沉积点")

# ===================================================================
# 3. 把 hit 钳进有效半导体体区，再构造 SSD hits 表 (DHE 格式)
# ===================================================================
# Geant4 晶体实心 r≤15；但 SSD 里 n+ 侧电极占 r∈[14.9,15.1]，电荷云不能落在
# 电极内或晶体外。落在这些区域的 hit 沿径向钳到 r≤14.8、z 钳到 [0.15,9.85]
# (保留其能量)；钳后仍无效的 (如划痕腔内) 丢弃。Geant4 (x,y,z)mm → 单位 m。
const MM = 1.0f-3
function clamp_hit(x, y, z)                          # 单位 mm
    r = hypot(x, y)
    if r > 14.8
        x *= 14.8 / r; y *= 14.8 / r
    end
    return x, y, clamp(z, 0.15, 9.85)
end
# 钳位 + 过滤无效 hit：就地修改 hits / evlist，返回 (钳位数, 丢弃数)
function prepare_events!(hits, evlist, sim)
    sc_geom  = sim.detector.semiconductor.geometry
    ct_geoms = [c.geometry for c in sim.detector.contacts]
    in_active(p) = (p in sc_geom) && !any(p in g for g in ct_geoms)
    nc = 0; nd = 0
    for ev in evlist
        fixed = NTuple{4,Float64}[]
        for h in hits[ev]
            x, y, z = clamp_hit(h[1], h[2], h[3])
            (x, y, z) != (h[1], h[2], h[3]) && (nc += 1)
            if in_active(CartesianPoint{Float32}(x*MM, y*MM, z*MM))
                push!(fixed, (x, y, z, h[4]))
            else
                nd += 1
            end
        end
        hits[ev] = fixed
    end
    filter!(ev -> !isempty(hits[ev]), evlist)
    nc, nd
end
n_clamped, n_dropped = prepare_events!(hits, evlist, sim)
N = length(evlist)
println("    钳位 ", n_clamped, " 个 hit、丢弃 ", n_dropped, " 个 → 有效事件 ", N)

mcevents = NamedTuple[
    (evtno = ev, detno = 1,
     thit  = fill(0.0u"ns", length(hits[ev])),
     edep  = [h[4] for h in hits[ev]] .* u"keV",
     pos   = [CartesianPoint{T}(h[1]*MM, h[2]*MM, h[3]*MM) for h in hits[ev]])
    for ev in evlist
]

# ===================================================================
# 4. 多线程模拟波形 (轮转分块 → Threads.@threads)
# ===================================================================
println("\n[3] 多线程漂移电荷 + 生成波形 (", Threads.nthreads(), " 线程)..."); flush(stdout)
nchunks = min(Threads.nthreads(), N)
chunks  = [collect(c:nchunks:N) for c in 1:nchunks]
results = Vector{Any}(undef, nchunks)
t_sig = @elapsed begin
    Threads.@threads for c in 1:nchunks
        results[c] = simulate_waveforms(mcevents[chunks[c]], sim;
            Δt = 1u"ns", max_nsteps = 2000, signal_unit = u"keV", verbose = false)
    end
end
println("    波形模拟耗时 ", round(t_sig, digits=1), " s  (",
        round(t_sig/N*1000, digits=1), " ms/事件)"); flush(stdout)

# 把 p+ 读出波形按事件顺序排好
all_wfs = Vector{Any}(undef, N)
for c in 1:nchunks
    rc  = results[c]
    idx = findall(==(READOUT_ID), rc.chnid)
    for (k, ei) in enumerate(chunks[c])
        all_wfs[ei] = rc.waveform[idx[k]]
    end
end

# ===================================================================
# 5. 每事件: 沉积能量 / 测得能量 / 质心
# ===================================================================
edep  = [sum(h[4] for h in hits[ev]) for ev in evlist]                       # keV
emeas = [abs(Float64(ustrip(all_wfs[i].signal[end]))) for i in 1:N]          # keV (波形平台)
rc_mm = zeros(N); zc_mm = zeros(N)
for (i, ev) in enumerate(evlist)
    hs = hits[ev]; w = sum(h[4] for h in hs)
    rc_mm[i] = sum(h[4]*hypot(h[1],h[2]) for h in hs)/w
    zc_mm[i] = sum(h[4]*h[3] for h in hs)/w
end

# ===================================================================
# 6. 保存 events.txt
# ===================================================================
open(joinpath(OUTDIR, "events.txt"), "w") do io
    println(io, "# Geant4→SSD 信号模拟 — ", TAG, " — 每事件一行")
    println(io, "# event_id  n_hits  r_mm  z_mm  edep_keV  emeas_keV")
    for i in 1:N
        @views println(io, evlist[i], "  ", length(hits[evlist[i]]), "  ",
            round(rc_mm[i],digits=3), "  ", round(zc_mm[i],digits=3), "  ",
            round(edep[i],digits=4), "  ", round(emeas[i],digits=4))
    end
end
println("[4] → events.txt  (", N, " 行)")

# ===================================================================
# 7. 保存全部波形 waveforms.bin  (平台值补齐到统一长度 L)
# ===================================================================
Ls = [length(all_wfs[i].signal) for i in 1:N]
L  = min(maximum(Ls), WF_CAP)
W  = Matrix{Float32}(undef, L, N)
for i in 1:N
    s  = Float32.(ustrip.(all_wfs[i].signal))
    li = length(s)
    if li >= L
        @views W[:, i] .= s[1:L]
    else
        @views W[1:li, i] .= s
        @views W[li+1:L, i] .= s[li]                  # 补平台值
    end
end
open(joinpath(OUTDIR, "waveforms.bin"), "w") do io
    write(io, Int64(N), Int64(L)); write(io, W)
end
println("    → waveforms.bin  (", N, " 条波形 × ", L, " ns, Float32 二进制)")

# ===================================================================
# 8. 能谱 (测得能量直方图)
# ===================================================================
bw   = 0.5
bins = 0:bw:130
ctr  = collect(bins)[1:end-1] .+ bw/2
hm   = [count(e -> b ≤ e < b+bw, emeas) for b in bins[1:end-1]]              # 测得
hd   = [count(e -> b ≤ e < b+bw, edep)  for b in bins[1:end-1]]              # 沉积(参考)
open(joinpath(OUTDIR, "spectrum.txt"), "w") do io
    println(io, "# 能谱 — ", TAG, " — ", N, " 事件, 测得能量(p+读出电极波形平台)")
    println(io, "# E_keV(bin中心)  counts_measured  counts_deposited")
    for k in eachindex(ctr)
        println(io, round(ctr[k],digits=3), "  ", hm[k], "  ", hd[k])
    end
end
plt_s = plot(size=(960,600), title="Energy spectrum — $TAG ($N events)",
             xlabel="energy (keV)", ylabel="counts / $(bw) keV", yscale=:log10,
             legend=:topleft)
plot!(plt_s, ctr, max.(hd,0.1), seriestype=:steppost, lw=1.2, color=:gray,
      label="Geant4 deposited")
plot!(plt_s, ctr, max.(hm,0.1), seriestype=:steppost, lw=1.4, color=:navy,
      label="SSD measured (p+ signal)")
vline!(plt_s, [122.0], color=:red, ls=:dash, lw=1, label="122 keV")
savefig(plt_s, joinpath(OUTDIR, "spectrum.png"))
println("    → spectrum.png / spectrum.txt")

# ===================================================================
# 9. 各能量段代表波形
# ===================================================================
bands = [(119.0,124.0,:red,"photopeak ~122"),
         (105.0,116.0,:orange,"escape ~112"),
         (60.0,90.0,:green,"Compton ~75"),
         (33.0,43.0,:purple,"Compton edge ~39"),
         (8.0,22.0,:steelblue,"low ~15")]
plt_w = plot(size=(960,620), title="Representative p+ waveforms by energy — $TAG",
             xlabel="time (ns)", ylabel="signal (keV)", legend=:bottomright)
for (lo,hi,col,lab) in bands
    sel = findall(i -> lo ≤ emeas[i] < hi, 1:N)
    isempty(sel) && continue
    for (j,i) in enumerate(sel[1:min(end,25)])
        plot!(plt_w, 1:L, W[:,i], lw=0.8, alpha=0.45, color=col,
              label = j==1 ? lab : "")
    end
end
savefig(plt_w, joinpath(OUTDIR, "waveforms_by_energy.png"))
println("    → waveforms_by_energy.png")

# ===================================================================
# 10. 汇总
# ===================================================================
println("\n[5] 汇总:")
println("    测得能量范围 ", round(minimum(emeas),digits=1), " – ",
        round(maximum(emeas),digits=1), " keV")
println("    122keV 光电峰 (120–124): ", count(e->120≤e≤124,emeas), " / ", N,
        "  (", round(100*count(e->120≤e≤124,emeas)/N,digits=1), "%)")
println("    波形模拟 ", round(t_sig,digits=1), " s")
println("\n", "="^60)
println("完成: output/", TAG, "/  (events.txt, waveforms.bin, spectrum.png, waveforms_by_energy.png)")
println("="^60)
