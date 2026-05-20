# =============================================================================
# 任务 05：把任务 04 的诱导电荷 Q(t) 转换成更真实的波形
#
# 三种波形对比：
#   (1) Q_induced(t)  — Shockley-Ramo 直接给的诱导电荷
#                       物理上等价于"理想 CSP 无衰减输出"
#                       形状：从 0 阶跃上升到 Q_total = 41356 e- (122 keV)
#
#   (2) I(t) = dQ/dt  — 真正的"电流脉冲"
#                       宽度 = 漂移时间 (~60 ns ~ 446 ns)，0.5 ms 尺度看就是一根针
#
#   (3) V_CSP(t)      — 现实中电荷灵敏前放的输出
#                       模型：dV/dt + V/τ = I(t)，τ = R_f·C_f ≈ 50 μs
#                       离散化递归：V[n] = α·V[n-1] + ΔQ[n], α = exp(-dt/τ)
#                       形状：陡升 + 指数衰减
#
# 输入：output/signals/wf_*.csv (任务 04 写出)
# 输出：output/plots/waveforms_*.png
# =============================================================================

using DelimitedFiles
using Plots
gr()

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const SIGNAL_DIR   = joinpath(PROJECT_ROOT, "output", "signals")
const PLOTS_DIR    = joinpath(PROJECT_ROOT, "output", "plots")

# ----- 物理参数 -----
const TAU_US   = 50.0          # CSP 衰减时间常数 (μs)，典型 HPGe 前放值
const T_MAX_US = 500.0         # 显示窗口 = 0.5 ms
const DT_US    = 0.01          # 重采样步长 10 ns，保证早期上升细节 + 长尾衰减

# ----- 读 summary 拿到事件列表 -----
summary = readdlm(joinpath(SIGNAL_DIR, "events_summary.csv"), ',', skipstart = 1)
labels = string.(summary[:, 1])
rs_mm  = Float64.(summary[:, 2])
zs_mm  = Float64.(summary[:, 3])
N_evt  = length(labels)
@info "处理 $N_evt 个事件"

# ----- 重采样后的均匀时间轴 -----
t_uniform_us = collect(0.0:DT_US:T_MAX_US)
N_t  = length(t_uniform_us)
α    = exp(-DT_US / TAU_US)

# ----- 读入每个 wf_*.csv，重采样到均匀网格，并计算 I 和 V_CSP -----
Q_arr = Matrix{Float64}(undef, N_t, N_evt)   # 诱导电荷 (e-)
I_arr = Matrix{Float64}(undef, N_t, N_evt)   # 电流 (e-/μs)
V_arr = Matrix{Float64}(undef, N_t, N_evt)   # CSP 输出 (任意单位，相对 Q_total)

for j in 1:N_evt
    data = readdlm(joinpath(SIGNAL_DIR, "wf_$(labels[j]).csv"), ',', skipstart = 1)
    t_raw = Float64.(data[:, 1])
    Q_raw = Float64.(data[:, 2])
    Q_total = Q_raw[end]

    # 1) Q 重采样到均匀网格；漂移结束后保持平台 Q_total
    for n in 1:N_t
        t = t_uniform_us[n]
        if t <= t_raw[end]
            i = searchsortedlast(t_raw, t)
            if i < 1
                Q_arr[n, j] = 0.0
            elseif i >= length(t_raw)
                Q_arr[n, j] = Q_raw[end]
            else
                Q_arr[n, j] = Q_raw[i] +
                    (Q_raw[i+1] - Q_raw[i]) * (t - t_raw[i]) / (t_raw[i+1] - t_raw[i])
            end
        else
            Q_arr[n, j] = Q_total
        end
    end

    # 2) 电流 I = dQ/dt
    I_arr[1, j] = 0.0
    for n in 2:N_t
        I_arr[n, j] = (Q_arr[n, j] - Q_arr[n-1, j]) / DT_US
    end

    # 3) CSP：一阶 RC 递归 V[n] = α V[n-1] + ΔQ[n]
    V_arr[1, j] = 0.0
    for n in 2:N_t
        ΔQ = Q_arr[n, j] - Q_arr[n-1, j]
        V_arr[n, j] = α * V_arr[n-1, j] + ΔQ
    end
end

# ----- 颜色编码：按总漂移时间排序（越长越红） -----
drift_times = [findlast(I_arr[:, j] .> maximum(I_arr[:, j]) * 0.01) for j in 1:N_evt]
drift_times = [isnothing(t) ? 0 : t for t in drift_times] .* DT_US
sortidx = sortperm(drift_times)
colors = palette(:turbo, N_evt)

# =============================================================================
# 图 1：CSP 输出，全 0.5 ms 视野（最贴近示波器观感）
# =============================================================================
p_csp = plot(
    xlabel = "Time (μs)",
    ylabel = "CSP output  (∝ e⁻, normalized)",
    title  = "Realistic CSP waveform, τ = $(TAU_US) μs   [122 keV @ V_op = 1000 V]",
    size = (1000, 550), dpi = 150, legend = false,
)
for k in 1:N_evt
    j = sortidx[k]
    plot!(p_csp, t_uniform_us, V_arr[:, j], lw = 1.3, color = colors[k], alpha = 0.85)
end
savefig(p_csp, joinpath(PLOTS_DIR, "waveforms_CSP_500us.png"))
@info "保存" file = "waveforms_CSP_500us.png"

# 同一份数据放对数 y 轴，看清衰减
p_csp_log = plot(
    xlabel = "Time (μs)",
    ylabel = "CSP output (log)",
    title  = "CSP waveform (log scale), τ = $(TAU_US) μs",
    size = (1000, 550), dpi = 150, legend = false, yscale = :log10,
    ylims = (1.0, 1e5),
)
for k in 1:N_evt
    j = sortidx[k]
    plot!(p_csp_log, t_uniform_us, max.(V_arr[:, j], 1.0), lw = 1.3, color = colors[k], alpha = 0.85)
end
savefig(p_csp_log, joinpath(PLOTS_DIR, "waveforms_CSP_500us_log.png"))
@info "保存" file = "waveforms_CSP_500us_log.png"

# =============================================================================
# 图 2：电流脉冲 I(t)，放大到 0..2 μs（在 500 μs 上是看不见的针）
# =============================================================================
zoom_us = 2.0
zoom_n  = round(Int, zoom_us / DT_US)
p_I = plot(
    xlabel = "Time (μs)",
    ylabel = "Induced current  I(t)  (e⁻/μs)",
    title  = "Current pulse  dQ/dt  (zoom to 0..$(zoom_us) μs)",
    size = (1000, 550), dpi = 150, legend = false,
)
for k in 1:N_evt
    j = sortidx[k]
    plot!(p_I, t_uniform_us[1:zoom_n], I_arr[1:zoom_n, j],
          lw = 1.3, color = colors[k], alpha = 0.85,
          label = "$(labels[j]) (r=$(rs_mm[j]), z=$(zs_mm[j]))")
end
savefig(p_I, joinpath(PLOTS_DIR, "waveforms_current_2us.png"))
@info "保存" file = "waveforms_current_2us.png"

# =============================================================================
# 图 3：电流脉冲在 0..500 μs 上，证明确实只是一根针
# =============================================================================
p_I_full = plot(
    xlabel = "Time (μs)",
    ylabel = "Induced current  I(t)  (e⁻/μs)",
    title  = "Current pulse on 0..500 μs scale  (essentially a delta near t=0)",
    size = (1000, 550), dpi = 150, legend = false,
)
for k in 1:N_evt
    j = sortidx[k]
    plot!(p_I_full, t_uniform_us, I_arr[:, j], lw = 1.0, color = colors[k], alpha = 0.7)
end
savefig(p_I_full, joinpath(PLOTS_DIR, "waveforms_current_500us.png"))
@info "保存" file = "waveforms_current_500us.png"

# =============================================================================
# 图 4：CSP 输出放大到上升段 (0..2 μs)，看见各事件的"形状指纹"
# =============================================================================
p_csp_zoom = plot(
    xlabel = "Time (μs)",
    ylabel = "CSP output  (∝ e⁻)",
    title  = "CSP waveform — rising edge zoom (0..$(zoom_us) μs)",
    size = (1000, 550), dpi = 150, legend = :bottomright,
)
for k in 1:N_evt
    j = sortidx[k]
    plot!(p_csp_zoom, t_uniform_us[1:zoom_n], V_arr[1:zoom_n, j],
          lw = 1.5, color = colors[k], alpha = 0.85,
          label = "$(labels[j])")
end
savefig(p_csp_zoom, joinpath(PLOTS_DIR, "waveforms_CSP_2us.png"))
@info "保存" file = "waveforms_CSP_2us.png"

println("\n完成。图：")
println("  waveforms_CSP_500us.png       — 现实 CSP 波形，0.5 ms 全景")
println("  waveforms_CSP_500us_log.png   — 同上，对数 y 轴看清衰减")
println("  waveforms_CSP_2us.png         — CSP 波形上升段放大")
println("  waveforms_current_2us.png     — 电流脉冲 I(t)，2 μs 视野")
println("  waveforms_current_500us.png   — 电流脉冲在 0.5 ms 上 (只是一根针)")
