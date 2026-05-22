using SolidStateDetectors
using Unitful
using Printf
using CUDA

# ===================================================================
# PPC 高纯锗探测器 — 电场 / 耗尽电压 / 信号模拟  (有/无划痕 A/B 对照)
# ===================================================================
# 用法:
#   julia -t auto --project=<SSD根目录> ppc/simulate.jl [阶段] [配置...]
#     阶段  : depletion → 只算 电场 + 耗尽电压
#             full (默认) → 再算信号波形
#     配置  : 省略 → 默认依次跑 config.yaml 与 config_noscratch.yaml
#             也可显式给一个或多个 yaml
#   例:
#     julia -t auto --project=. ppc/simulate.jl depletion
#     julia -t auto --project=. ppc/simulate.jl full config.yaml
#
# 性能要点:
#   · 多个配置在「同一个 julia 进程」里依次跑 —— JIT 编译只付一次
#   · 用 `-t auto` 多线程启动 —— SSD 的网格 setup 等可并行部分会用上多核
#   · depletion 阶段只解「电势 + 偏压电极权重势」，跳过用不到的另一个
#     电极权重势和电场 —— 比 simulate!(全解) 少 ~1/3 求解量
#   · 结果按配置写入独立子目录 output/<tag>/
# ===================================================================

# ---- 命令行参数 ----
STAGE   = length(ARGS) >= 1 ? ARGS[1] : "full"
CONFIGS = length(ARGS) >= 2 ? ARGS[2:end] : ["config.yaml", "config_noscratch.yaml"]

SCRIPT_DIR = @__DIR__
ROOT_DIR   = dirname(SCRIPT_DIR)
T = Float32

# 网格控制: r、z 最大网格间距。主要作用是给「轴对称基线」一个够细的网格
# (基线 φ 方向只有 ~8 点、很便宜，可放心加细)。0.3mm 对 3D 划痕模型几乎无
# 影响 —— 划痕几何本身播下的网格已比这更细。
const MAX_TICK = 0.3u"mm"

# 求解设备 (一次性确定)
const DEVICE = CUDA.functional() ? CuArray : Array

# 公共求解参数 (收敛限 1e-5 ≈ 0.01V，远细于耗尽电压 0.1V 容差，足够且更快)
const SOLVER_KW = (
    convergence_limit = 1e-5,
    refinement_limits = [0.2, 0.1, 0.05],
    max_tick_distance = MAX_TICK,
    device_array_type = DEVICE,
    verbose           = true,
)

# ===================================================================
# 单个配置的完整流程
# ===================================================================
function run_config(config_file::AbstractString, stage::AbstractString)
    depletion_only = (stage == "depletion")
    config_path = isabspath(config_file) ? config_file : joinpath(SCRIPT_DIR, config_file)
    tag = config_file == "config.yaml"           ? "scratch"   :
          config_file == "config_noscratch.yaml" ? "noscratch" :
          first(splitext(basename(config_file)))
    output_dir = joinpath(SCRIPT_DIR, "output", tag)
    mkpath(output_dir)

    println("\n", "█"^60)
    println("█ 配置: ", config_file, "   (tag = ", tag, ")")
    println("█"^60)

    # --- 1. 加载 ---
    sim = Simulation{T}(config_path)
    det = sim.detector
    println("[1] 探测器: ", det.name, "  /  ", det.semiconductor.material.name,
            "  /  ", det.semiconductor.temperature)
    for c in det.contacts
        println("    电极 id$(c.id)  ", rpad(c.name, 12), c.potential, " V")
    end
    bias_id = det.contacts[argmax([abs(c.potential) for c in det.contacts])].id   # 偏压电极

    # --- 2. 求解场 ---
    if depletion_only
        # 耗尽电压只需要「电势 + 偏压电极权重势」—— 跳过另一电极权重势和电场
        println("[2] 求解 电势 + 偏压电极(id$bias_id)权重势  (设备: ",
                DEVICE === CuArray ? "GPU $(CUDA.name(CUDA.device()))" : "CPU", ")")
        calculate_electric_potential!(sim; SOLVER_KW...)
        calculate_weighting_potential!(sim, bias_id; SOLVER_KW...)
    else
        # 信号模拟需要 电场 + 所有电极权重势
        println("[2] 求解 电势 + 全部权重势 + 电场  (设备: ",
                DEVICE === CuArray ? "GPU $(CUDA.name(CUDA.device()))" : "CPU", ")")
        simulate!(sim; SOLVER_KW...)
    end
    depleted = is_depleted(sim.point_types)
    println("    完全耗尽? ", depleted)

    # --- 3. 耗尽电压 ---
    bias = maximum(c.potential for c in det.contacts)
    deplV = depleted ? estimate_depletion_voltage(sim; contact_id = bias_id, verbose = true) : nothing
    if depleted
        println("[3] 耗尽电压 ≈ ", round(typeof(1.0u"V"), deplV))
    else
        println("[3] ⚠ 在 ", bias, " V 下未完全耗尽，无法估算耗尽电压")
    end
    open(joinpath(output_dir, "depletion_voltage.txt"), "w") do io
        println(io, "# PPC HPGe — 耗尽电压估算")
        println(io, "配置     : ", config_file, "   (tag = ", tag, ")")
        println(io, "探测器   : ", det.name)
        println(io, "工作偏压 : ", bias, " V")
        println(io, "完全耗尽 : ", depleted)
        println(io, "耗尽电压 : ", isnothing(deplV) ? "无法估算 (未完全耗尽)" :
                string(round(typeof(1.0u"V"), deplV)))
    end

    depletion_only && return (tag = tag, depleted = depleted, deplV = deplV)

    # --- 4. 信号波形 (full 阶段) ---
    println("[4] 信号波形模拟")
    event_specs = [
        ("center", 0.0, 5.0), ("mid-radius", 7.0, 5.0), ("near-edge", 12.0, 5.0),
        ("near-top", 5.0, 9.0), ("near-bottom", 5.0, 1.0), ("corner", 3.0, 1.0),
    ]
    EVENT_ENERGY = 661.7u"keV"                                  # Cs-137 光电峰
    println("    粗略 (Δt=1ns):")
    for (label, r_mm, z_mm) in event_specs
        cp = CylindricalPoint{T}(r_mm * 1e-3, 0.0, z_mm * 1e-3)
        ev = Event([CartesianPoint(cp)], [EVENT_ENERGY])
        simulate!(ev, sim, Δt = 1e-9, max_nsteps = 10000)
        peaks = join(["c$(ci)=$(round(maximum(abs.(ustrip.(wf.signal))), digits=4))"
                      for (ci, wf) in enumerate(ev.waveforms)], ", ")
        println("      ", rpad(label, 14), peaks)
    end
    println("    高精度中心点 (Δt=0.5ns):")
    ev = Event([CartesianPoint(CylindricalPoint{T}(0.0, 0.0, 5e-3))], [EVENT_ENERGY])
    simulate!(ev, sim, Δt = 5e-10, max_nsteps = 20000)
    for (ci, wf) in enumerate(ev.waveforms)
        path = joinpath(output_dir, "waveform_contact$(ci).txt")
        open(path, "w") do io
            println(io, "# PPC HPGe ($tag) — Contact $ci")
            println(io, "# t(ns)  signal")
            for (t, s) in zip(ustrip.(wf.time .* 1e9), ustrip.(wf.signal))
                @printf(io, "%.4f  %.8f\n", t, s)
            end
        end
        println("      → ", relpath(path, ROOT_DIR))
    end
    return (tag = tag, depleted = depleted, deplV = deplV)
end

# ===================================================================
# 主流程 — 同一进程内依次跑所有配置 (JIT 编译只付一次)
# ===================================================================
println("="^60)
println("PPC HPGe 模拟")
println("  阶段  : ", STAGE)
println("  配置  : ", join(CONFIGS, ",  "))
println("  线程  : ", Threads.nthreads(), "    设备: ",
        DEVICE === CuArray ? "GPU" : "CPU")
println("="^60)

results = NamedTuple[]
for cfg in CONFIGS
    push!(results, run_config(cfg, STAGE))
    GC.gc()
    CUDA.functional() && CUDA.reclaim()       # 释放 GPU 显存，给下个配置腾地方
end

println("\n", "="^60)
println("汇总")
for r in results
    println("  ", rpad(r.tag, 12),
            r.depleted ? "耗尽电压 ≈ $(round(typeof(1.0u"V"), r.deplV))" : "未完全耗尽")
end
println("="^60)
