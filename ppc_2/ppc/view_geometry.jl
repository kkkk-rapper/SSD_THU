# 可视化脚本 — 在 Julia REPL 中交互查看探测器结构
# 用法 (在 SSD 根目录下):
#   julia --project=. -i ppc/view_geometry.jl
# 脚本用 @__DIR__ 定位 config.yaml，因此从哪个目录启动都能找到配置文件。
#
# 进入 REPL 后可交互旋转、缩放

using SolidStateDetectors
using Plots; plotlyjs()   # 交互式 3D (浏览器中打开)

# 路径解析 — 不依赖当前工作目录
CONFIG_PATH = joinpath(@__DIR__, "config.yaml")
OUTPUT_DIR  = joinpath(@__DIR__, "output")

sim = Simulation{Float32}(CONFIG_PATH)

# 1. 探测器 3D 全貌
# 注意: SSD 的几何绘图使用内部单位「米」，所以坐标轴是 m，不是 mm。
plot(sim.detector,
    title = "PPC HPGe — Scratch on p+ Surface",
    xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
    size = (1000, 800),
    legend = false,
)
# 在浏览器中打开 → 鼠标拖拽旋转、滚轮缩放

# 2. 保存静态图 (取消注释即可)
# savefig(joinpath(OUTPUT_DIR, "detector_3d.png"))
