using SolidStateDetectors
using Unitful
using Printf

# ===================================================================
# PPC 高纯锗探测器 — 电场 & 信号波形模拟
# 用法: julia --project=<SSD根目录> simulate.jl
# ===================================================================

# 路径解析
SCRIPT_DIR  = @__DIR__                          # ppc/
ROOT_DIR    = dirname(SCRIPT_DIR)               # SSD 根目录
CONFIG_PATH = joinpath(SCRIPT_DIR, "config.yaml")
OUTPUT_DIR  = joinpath(SCRIPT_DIR, "output")
mkpath(OUTPUT_DIR)

# 激活根目录的 Julia 环境
import Pkg; Pkg.activate(ROOT_DIR)

T = Float32

println("="^60)
println("PPC HPGe 探测器模拟")
println("  配置: ", relpath(CONFIG_PATH, ROOT_DIR))
println("  输出: ", relpath(OUTPUT_DIR, ROOT_DIR))
println("="^60)

# ===================================================================
# 1. 加载探测器配置
# ===================================================================
sim = Simulation{T}(CONFIG_PATH)
det  = sim.detector
println("\n[1/4] 探测器: ", det.name)
println("  材料: ", det.semiconductor.material)
println("  温度: ", det.semiconductor.temperature)
println("  偏压: ", det.contacts[1].potential, " V (contact1), ",
                         det.contacts[2].potential, " V (contact2)")

# ===================================================================
# 2. 计算电势 + 电场 + 权重势
# ===================================================================
println("\n[2/4] 求解电场 & 权重势 (SOR + 网格精化)...")
simulate!(sim,
    convergence_limit = 1e-6,
    refinement_limits = [0.2, 0.1, 0.05],
    verbose = true
)
depleted = is_depleted(sim.point_types)
println("  完全耗尽? ", depleted)
if !depleted
    println("  ⚠ 探测器未完全耗尽，建议提高偏压或降低杂质浓度")
end

# ===================================================================
# 3. 估算耗尽电压
# ===================================================================
# 注意: estimate_depletion_voltage 内部 @assert is_depleted(...)，
#       探测器未完全耗尽时会抛异常。这里先判断，未耗尽则跳过，
#       否则脚本会在此中断、跑不到第 4 步的信号模拟。
println("\n[3/4] 耗尽电压估算...")
if depleted
    deplV = estimate_depletion_voltage(sim, verbose = true)
    println("  耗尽电压 ≈ ", round(typeof(1.0u"V"), deplV))
else
    println("  跳过: 探测器未完全耗尽，estimate_depletion_voltage 仅适用于完全耗尽情形。")
    println("       请提高偏压使其完全耗尽后重跑 (或调用时传 check_for_depletion = false)。")
end

# ===================================================================
# 4. 多位置信号模拟
# ===================================================================
println("\n[4/4] 信号波形模拟")

# 相互作用位置 — 格式: (标签, r, z)，r/z 单位为 mm
# (下面构造 CylindricalPoint 时用 *1e-3 转成 SSD 内部单位「米」)
event_specs = [
    ("center",         0.0,  5.0),
    ("mid-radius",     7.0,  5.0),
    ("near-edge",     12.0,  5.0),
    ("near-top",       5.0,  9.0),
    ("near-bottom",    5.0,  1.0),
    ("corner",         3.0,  1.0),
]

# 每个事件的能量沉积 — Cs-137 光电峰 661.7 keV
# (不指定能量时 SSD 默认按 1 eV 处理，信号幅度会变成 1eV/2.95eV≈0.339，无物理意义)
EVENT_ENERGY = 661.7u"keV"

# 粗略模拟 (快速)
println("  粗略模拟 (Δt=1ns):")
for (label, r_mm, z_mm) in event_specs
    cyl_pt = CylindricalPoint{T}(r_mm * 1e-3, 0.0, z_mm * 1e-3)
    evt = Event([CartesianPoint(cyl_pt)], [EVENT_ENERGY])
    simulate!(evt, sim, Δt = 1e-9, max_nsteps = 10000)
    peaks = join(["c$(ci)=$(round(maximum(abs.(ustrip.(wf.signal))), digits=4))"
                  for (ci, wf) in enumerate(evt.waveforms)], ", ")
    println("    $(rpad(label, 14))  $(peaks)")
end

# 高精度中心点波形
println("\n  高精度中心点波形 (Δt=0.5ns):")
cyl_pt = CylindricalPoint{T}(0.0, 0.0, 5e-3)
evt = Event([CartesianPoint(cyl_pt)], [EVENT_ENERGY])
simulate!(evt, sim, Δt = 5e-10, max_nsteps = 20000)

for (ci, wf) in enumerate(evt.waveforms)
    t_ns = ustrip.(wf.time .* 1e9)
    sig  = ustrip.(wf.signal)
    path = joinpath(OUTPUT_DIR, "waveform_contact$(ci).txt")
    open(path, "w") do io
        println(io, "# PPC HPGe Detector — Contact $ci")
        println(io, "# t(ns)  signal")
        for (t, s) in zip(t_ns, sig)
            @printf(io, "%.4f  %.8f\n", t, s)
        end
    end
    println("    → $(relpath(path, ROOT_DIR)) ($(length(t_ns)) 步)")
end

println("\n" * "="^60)
println("完成: $(relpath(OUTPUT_DIR, ROOT_DIR))/")
println("="^60)
