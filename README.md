# PPC 高纯锗 γ 探测器仿真 —— 表面划痕影响研究

基于 [SolidStateDetectors.jl](https://github.com/JuliaPhysics/SolidStateDetectors.jl)（简称 SSD）的**点电极（PPC）高纯锗（HPGe）γ 探测器**仿真。核心目标：用 **A/B 对照**研究 **p+ 电极表面一道划痕**对探测器电场、电荷收集和能谱的影响。

## 物理背景与目标

真实 HPGe 探测器表面的机械划痕会扰动局部电场，可能影响电荷收集。本项目对**两台除划痕外完全相同**的探测器做并排仿真：

- **有划痕**（`config.yaml`）—— p+ 点电极表面有一道带锥度的沟槽
- **无划痕基线**（`config_noscratch.yaml`）—— 其余参数逐字段一致

通过对比两者的电场、信号波形与能谱，定量考察划痕的影响。

## 探测器模型

| 项 | 值 |
|---|---|
| 晶体 | 高纯锗圆柱，φ30 mm × 高 10 mm，78 K |
| 杂质 | p 型，沿 z 线性分布（约 2.8–3.0 ×10⁶ mm⁻³） |
| n+ 电极 | 底面 + 侧壁，偏压 1000 V |
| p+ 点电极 | 顶面中心 φ8 mm 圆盘，接地，**读出电极** |
| 钝化层 | 顶面裸露区的浮空 HPGe 钝化膜 |
| 电荷漂移 | 各向异性 `ADLChargeDriftModel` |
| 划痕 | 顶面一道沟槽：沿 x、偏心 1.5 mm、长 3 mm、最深 0.16 mm；用 `Polycone` 回转体建模 |

耗尽电压 ≈ **124 V**（1000 V 工作偏压下完全耗尽）。

## 仿真流程

```
Geant4（122 keV γ 源，独立项目）          本仓库：SolidStateDetectors.jl
──────────────────────────────          ──────────────────────────────
γ 在锗晶体中沉积能量    ── hits.csv ──►   解电场 / 权重势（CUDA）
(event, x, y, z, edep)                          │
                                         逐事件漂移电荷（多线程）
                                                │
                                         p+ 读出电极感应波形 ──► 能谱
```

Geant4 部分为独立项目（122 keV 单能 γ 点源，输出每步能量沉积 `hits.csv`）；本仓库负责其后的电场求解、电荷漂移与信号 / 能谱。

## 目录结构

```
ppc_2/
├── Project.toml / Manifest.toml    Julia 工程环境
├── CLAUDE.md                       项目说明（给 AI 助手）
└── ppc/
    ├── config.yaml                 有划痕配置
    ├── config_noscratch.yaml       无划痕基线配置
    ├── simulate.jl                 电场 + 耗尽电压
    ├── field.jl                    电场 / 电势截面图
    ├── signal.jl                   Geant4 hits → 波形 → 能谱  ★主信号链
    ├── render_geometry.jl          几何可视化（3D + 2D 截面分类图）
    ├── view_geometry.jl            交互式 3D 几何
    ├── verify_configs.jl           核对两配置仅划痕不同
    ├── compare_spectra.jl          能谱对比（线性坐标）
    ├── view_waveforms.jl           按能量段调看已存波形
    ├── input/                      Geant4 hits.csv 等输入
    └── output/{scratch,noscratch}/ 各配置的结果
```

## 环境依赖

- Julia 1.12+，SolidStateDetectors.jl 0.11+（精确版本见 `ppc_2/Manifest.toml`）
- 可选：NVIDIA GPU + CUDA.jl —— 脚本自动检测，有则用 GPU 解电场
- Geant4 部分单独运行，产出 `hits.csv`

## 使用方法

所有命令在 `ppc_2/` 目录下运行，`-t auto` 启用多线程：

```bash
cd ppc_2

# 耗尽电压（快，跳过波形）
julia -t auto --project=. ppc/simulate.jl depletion

# 1000V 电场 / 电势截面图（有 / 无划痕各一套）
julia -t auto --project=. ppc/field.jl

# Geant4 hits → 信号波形 → 能谱    参数：[配置] [事件数]
julia -t auto --project=. ppc/signal.jl config_noscratch.yaml 30000
julia -t auto --project=. ppc/signal.jl config.yaml 30000

# 有 / 无划痕能谱对比（线性坐标）
julia --project=. ppc/compare_spectra.jl

# 按测得能量段调看已保存的波形
julia --project=. ppc/view_waveforms.jl noscratch 120 124

# 几何可视化 / 配置校验
julia --project=. ppc/render_geometry.jl
julia --project=. ppc/verify_configs.jl
```

结果写入 `ppc/output/<tag>/`（tag 由配置名推出：`config.yaml`→`scratch`，`config_noscratch.yaml`→`noscratch`）。

## 主要结果

- **耗尽电压 ≈ 124 V**；1000 V 下完全耗尽。
- **122 keV γ 能谱**复现良好：光电峰（约 87%）、康普顿连续谱与边缘（约 39 keV）、锗 X 射线逃逸峰（约 112 keV）。
- **划痕对能谱的影响**：在当前模型（无电荷陷阱、完全耗尽）下，电荷 100% 收集 → 测得能量 ≡ 沉积能量 → **有 / 无划痕的能谱基本重合**。划痕改变的是**脉冲波形的形状与时序**，而非收集到的总电荷。若要让划痕在**能谱**上产生可见差异，需引入电荷陷阱等电荷损失机制。

## 注意事项

- **3D 模型的 `max_tick_distance` 必须用元组 `(Δr, Δφ, Δz)`** —— 标量值会按弧长强制加密 φ 方向、令网格爆炸。
- **单位**：YAML 配置用 mm；SSD 内部 / 绘图用 m，脚本会显式转换。
- **两份配置必须除划痕外逐字段一致** —— 改任何物理 / 网格参数都要同步改两份，`verify_configs.jl` 用于校验。
- **波形二进制 `waveforms.bin`（每个约 190 MB）不入库**（超 GitHub 文件大小限制），由 `signal.jl` 重算。

详细变更见 [CHANGELOG.md](CHANGELOG.md)。
