# 配置对照验证 — 构建 config.yaml(有划痕) 与 config_noscratch.yaml(基线) 两份
# 探测器，打印关键参数并核对：除「划痕」外两者应完全一致。
# 用法: julia --project=<SSD根目录> ppc/verify_configs.jl

using SolidStateDetectors

SCRIPT_DIR = @__DIR__
mm = 1.0f-3
pt = CartesianPoint{Float32}(0.0f0, 1.5f0 * mm, 9.95f0 * mm)   # 划痕沟槽腔内一点

function show_cfg(cfg)
    println("─"^60)
    println("配置: ", cfg)
    sim = Simulation{Float32}(joinpath(SCRIPT_DIR, cfg))
    det = sim.detector
    println("  name        : ", det.name)
    println("  material    : ", det.semiconductor.material.name)
    println("  temperature : ", det.semiconductor.temperature, " K")
    for c in det.contacts
        println("  contact id$(c.id) : ", rpad(c.name, 12), " ", c.potential, " V")
    end
    println("  偏压(最高电位): ", maximum(c.potential for c in det.contacts), " V")
    println("  passives    : ", det.passives === nothing ? 0 : length(det.passives), " 个 (钝化层)")
    insc = pt in det.semiconductor.geometry
    inp  = pt in det.contacts[2].geometry
    println("  点(0,1.5,9.95)mm: ∈半导体=", insc, "  ∈p+电极=", inp,
            insc ? "   → 实心晶体 (无划痕)" : "   → 沟槽腔 (有划痕)")
    return nothing
end

println("\n", "="^60)
println("配置对照验证")
println("="^60)
show_cfg("config.yaml")
show_cfg("config_noscratch.yaml")
println("─"^60)
println("预期: 两者 偏压/材料/温度/电极/钝化层 全部相同；")
println("      仅「点(0,1.5,9.95)」一项不同 —— 这正是划痕。")
