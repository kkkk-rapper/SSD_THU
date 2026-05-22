# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

基于 [SolidStateDetectors.jl](https://github.com/JuliaPhysics/SolidStateDetectors.jl)
(简称 SSD) 的 **PPC (点电极) 高纯锗 γ 探测器** Julia 仿真项目。
计算锗晶体的电势 / 电场、耗尽电压和感应信号波形 —— 并对两种几何做并排对照:
**带表面划痕的探测器** vs. 除划痕外完全相同的**无划痕基线**。

仓库根目录是一个 Julia 工程 (`Project.toml` / `Manifest.toml`),所有代码位于 `ppc/`。

## 常用命令

所有命令均从仓库根目录运行,并激活工程环境。`-t auto` 启用多线程
(SSD 的网格 setup 与波形模拟会用到)。

```bash
# 完整流程 (电场 + 耗尽电压 + 波形),对两份配置都跑
julia -t auto --project=. ppc/simulate.jl

# 只算耗尽电压 (跳过电场和波形,求解量少约 1/3)
julia -t auto --project=. ppc/simulate.jl depletion

# 只跑单个配置
julia -t auto --project=. ppc/simulate.jl full config.yaml

# 电势 / 电场截面图 (PNG),对两份配置都跑
julia -t auto --project=. ppc/field.jl

# 用真实 Geant4 hits 算波形 —— 参数: [配置] [事件数],默认 config_noscratch.yaml 100
julia -t auto --project=. ppc/signal.jl config.yaml 200

# 几何: 输出静态 PNG + 对 CSG 结果做「点是否属于」自检
julia --project=. ppc/render_geometry.jl

# 几何: 浏览器中交互式 3D (用 -i 让 REPL 保持打开)
julia --project=. -i ppc/view_geometry.jl

# 验证两份配置「仅划痕一项不同」
julia --project=. ppc/verify_configs.jl

# 从 Geant4 hits 画能谱 (Python,不是 Julia)
cd ppc/input && python3 plot_spectrum.py
```

本项目没有测试套件、linter 或构建步骤 —— 这些脚本本身就是工作流。

## 代码架构

**双配置 A/B 对照设计。** `ppc/config.yaml` (有划痕) 与
`ppc/config_noscratch.yaml` (基线) 应当**除划痕几何外逐字节完全一致**。
修改任何物理 / 网格参数 (偏压、杂质分布、漂移模型、网格、钝化层……) 时,
**必须同时改两份文件**。`verify_configs.jl` 就是用来检查两者是否产生了偏差。

**划痕几何。** 划痕是一条带锥度的沟槽,用 `Polycone` 回转体建模,旋转到全局
x 轴后与晶体本体求交,只保留落在晶内 (z ≤ 10 mm) 的下半部。它通过 YAML 锚点
`&scratch` 定义一次,并用 `*scratch` **复用**为 p+ 电极。这样可保证「挖掉的
沟槽」与「p+ 覆盖的沟槽内壁」永远不会失配 —— 不要复制 polycone,保留锚点/别名。

**输出布局。** 脚本根据配置文件名推导出 `tag` (`config.yaml` → `scratch`,
`config_noscratch.yaml` → `noscratch`),结果写入 `ppc/output/<tag>/`。
每个脚本在**同一个 Julia 进程**里依次跑两份配置,JIT 编译只付一次。

**设备选择。** 每个脚本都自动检测 CUDA: `DEVICE = CUDA.functional() ?
CuArray : Array`。有 GPU 就用 GPU,否则用 CPU —— 无需改代码。

**探测器模型。** 78 K 下的高纯锗;电极 id 1 = n+ (底面 + 侧壁,偏压 1000 V),
电极 id 2 = p+ (顶面 φ8 mm 圆盘 + 划痕内壁,接地 = 读出电极)。顶面裸露区覆盖
一层浮空的 HPGe 钝化层。杂质浓度为 z 方向线性模型;电荷漂移用各向异性的
`ADLChargeDriftModel`。

## 关键注意事项

- **3D 模型的 `max_tick_distance` 必须是元组 `(Δr, Δφ, Δz)`。** 标量值会按弧长
  强制加密 φ 方向,令 3D 划痕网格爆炸。`simulate.jl` 故意用了标量 —— 但那只是
  为了给便宜的轴对称基线一个够细的网格;`field.jl` 和 `signal.jl` 用元组形式。
- **单位。** YAML 配置里的长度单位是 **mm**;SSD 内部 / 绘图单位是 **米**。
  脚本会显式转换 (`mm = 1.0f-3`、`xunit = u"mm"`)。
- **SSD 的 3D 几何图会画出每个图元的整面**,不按 CSG 裁切结果裁面 —— 所以
  `intersection` / `difference` 的结果在 3D 视图里看起来是错的。要判断真实
  几何,请看 `render_geometry.jl` 输出的 2D 点分类截面图。
- **`signal.jl` 读取的是硬编码的 `HITS_CSV` 路径** (`/root/geant4/sim_project/
  output/hits.csv`),不是仓库内的 `ppc/input/hits.csv`。运行 `signal.jl` 前请
  先核对 / 更新这个路径。
- **不同脚本的能量来源不同。** `simulate.jl` 注入 661.7 keV 的合成点事件
  (Cs-137);`signal.jl` 用真实 Geant4 hits (122 keV 源)。两者能标不可比较。
