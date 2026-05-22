# 几何渲染脚本 — 构建探测器、自检、并保存结构图 (静态 PNG)
# 用法: julia --project=<SSD根目录> ppc/render_geometry.jl
#
# 说明：
#  · SSD 自带的 3D 几何绘图按「图元整面」描边，不会按 CSG 裁切结果裁面，
#    所以 intersection 裁掉的部分在 3D 图里仍会被画出来 → 易误判。
#  · 因此本脚本除 3D 总览外，另用「2D 截面分类图」：在截面上逐点调用
#    `point ∈ geometry` 判断属于 晶体/电极/真空，再上色 —— 这才如实反映
#    布尔运算后的真实几何 (能看出划痕只剩 z≤10 的下半部)。

ENV["GKSwstype"] = "100"          # GR 离屏渲染 (无显示环境也能出图)

using SolidStateDetectors
using Plots
gr()

SCRIPT_DIR  = @__DIR__
ROOT_DIR    = dirname(SCRIPT_DIR)
CONFIG_PATH = joinpath(SCRIPT_DIR, "config.yaml")
OUTPUT_DIR  = joinpath(SCRIPT_DIR, "output")
mkpath(OUTPUT_DIR)

println("="^60)
println("构建探测器几何  配置: ", relpath(CONFIG_PATH, ROOT_DIR))
println("="^60)

sim = Simulation{Float32}(CONFIG_PATH)
det = sim.detector

println("\n探测器: ", det.name, "   电极数: ", length(det.contacts))
for c in det.contacts
    println("    · contact $(c.id): $(c.name)  →  $(c.potential) V")
end

sc = det.semiconductor.geometry          # 半导体 (晶体 − 划痕)
c1 = det.contacts[1].geometry            # n+ 电极
c2 = det.contacts[2].geometry            # p+ 电极 (φ8 圆盘 ∪ 划痕)
mm = 1.0f-3

# ===================================================================
# 几何自检 — 用「点是否属于」确认划痕已裁到 z≤10
# ===================================================================
println("\n几何自检 (✓ = 符合预期):")
checks = [   # (说明, 点[mm], 期望∈半导体, 期望∈p+电极)
    ("沟槽腔内   (0, 1.5,  9.95)", (0.0f0, 1.5f0,  9.95f0), false, true ),
    ("沟槽口上方 (0, 1.5, 10.05)", (0.0f0, 1.5f0, 10.05f0), false, false),  # ★裁切关键
    ("沟槽正下方 (0, 1.5,  9.50)", (0.0f0, 1.5f0,  9.50f0), true,  false),
    ("远处晶体内 (5, 0,    9.95)", (5.0f0, 0.0f0,  9.95f0), true,  false),
]
for (label, (x, y, z), exp_sc, exp_ct) in checks
    p = CartesianPoint{Float32}(x * mm, y * mm, z * mm)
    in_sc, in_ct = (p in sc), (p in c2)
    ok = (in_sc == exp_sc) && (in_ct == exp_ct)
    println("  ", ok ? "✓" : "✗", " ", label,
            "   半导体=", in_sc, "  p+电极=", in_ct,
            ok ? "" : "   ← 期望 半导体=$exp_sc p+=$exp_ct")
end

# ===================================================================
# 逐点分类
# ===================================================================
function classify(x, y, z)               # 坐标单位 m
    p = CartesianPoint{Float32}(x, y, z)
    p in c2 && return 3                   # p+ 电极 (含划痕)
    p in c1 && return 2                   # n+ 电极
    p in sc && return 1                   # Ge 晶体
    return 0                              # 真空
end

# 离散配色: 0 真空(白) / 1 Ge(浅蓝) / 2 n+(深蓝) / 3 p+(橙红)
catcolors = cgrad([:white, :deepskyblue, :royalblue, :orangered], categorical = true)
catnames  = ["vacuum", "Ge crystal", "n+ contact", "p+ contact (scratch)"]

# ===================================================================
# 2D 截面分类图 (heatmap)
# ===================================================================
function cross_section(fixed_axis, fixed_val, urange, vrange,
                       ulabel, vlabel, ttl, fname)
    us = range(urange[1], urange[2], step = urange[3])    # 单位 mm
    vs = range(vrange[1], vrange[2], step = vrange[3])
    M  = Array{Int}(undef, length(vs), length(us))
    for (j, u) in enumerate(us), (i, v) in enumerate(vs)
        x, y, z = fixed_axis === :x ? (fixed_val, u * mm, v * mm) :
                                      (u * mm, fixed_val, v * mm)
        M[i, j] = classify(Float32(x), Float32(y), Float32(z))
    end
    plt = heatmap(us, vs, M; color = catcolors, clims = (-0.5, 3.5),
                  colorbar = false, aspect_ratio = :equal, framestyle = :box,
                  size = (920, 640), title = ttl, xlabel = ulabel, ylabel = vlabel)
    # 晶体顶面参考线 z=10
    hline!(plt, [10.0]; color = :black, ls = :dash, lw = 1.2, label = "")
    # 手动图例 (画在图外的占位散点)
    for c in 1:3
        scatter!(plt, [NaN], [NaN]; label = catnames[c + 1],
                 color = [:deepskyblue, :royalblue, :orangered][c],
                 markerstrokewidth = 0, markersize = 6, legend = :outertopright)
    end
    savefig(plt, joinpath(OUTPUT_DIR, fname))
    println("→ ", fname)
end

# ===================================================================
# 出图
# ===================================================================
axis_kw = (xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)", legend = false)

println("\n[3D 总览]")
plot(det; size = (1000, 800), camera = (45, 25),
     title = "PPC HPGe — crystal + contacts", axis_kw...)
savefig(joinpath(OUTPUT_DIR, "geometry_3d_iso.png")); println("→ geometry_3d_iso.png")

plot(det; size = (900, 850), camera = (0, 90),
     title = "Top view", axis_kw...)
savefig(joinpath(OUTPUT_DIR, "geometry_3d_top.png")); println("→ geometry_3d_top.png")

println("\n[2D 截面分类图 — 如实反映布尔运算后的真实几何]")
# x-z 截面 @ y=1.5mm：划痕沿长度方向的深度剖面
cross_section(:y, 1.5f0 * mm, (-2.2, 2.2, 0.01), (9.55, 10.25, 0.005),
              "x (mm)", "z (mm)",
              "Cross-section x–z @ y=1.5mm  (scratch length profile)",
              "geometry_cross_xz.png")
# y-z 截面 @ x=0：划痕最深处的横截面 —— 应为 z≤10 的半圆，顶面以上是真空
cross_section(:x, 0.0f0, (0.9, 2.1, 0.005), (9.55, 10.25, 0.005),
              "y (mm)", "z (mm)",
              "Cross-section y–z @ x=0  (scratch cross-section — lower half only)",
              "geometry_cross_yz.png")

println("\n" * "="^60)
println("完成: ", relpath(OUTPUT_DIR, ROOT_DIR), "/")
println("="^60)
