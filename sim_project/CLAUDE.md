# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目简介

这是一个 Geant4 仿真程序：单能伽马点源照射一块高纯锗（HPGe）晶体。程序记录晶
体内每一个有能量沉积的步（相互作用位置 + 沉积能量），写入 `hits.csv`。该文件
是交给 `~/ssd_projects` 中 SolidStateDetectors.jl 信号仿真（`ppc`、`ppc_2`）的
输入——后者把这些 hit 转换成探测器脉冲，并最终得到能谱。

## 编译

默认环境中没有 Geant4，每个用于编译或运行的 shell 都必须先 source 它的环境
脚本：

```bash
source /root/geant4/g4_install/bin/geant4.sh
cd /root/geant4/sim_project/build
cmake -DCMAKE_PREFIX_PATH=/root/geant4/g4_install ..   # 只需执行一次 / 修改 CMakeLists 后再执行
make -j$(nproc)
```

## 运行

`hits.csv` 写入**当前工作目录**，因此批处理运行从 `output/` 目录启动，把结果
收集到那里。每次运行都会覆盖 `hits.csv`。

```bash
source /root/geant4/g4_install/bin/geant4.sh
cd /root/geant4/sim_project/output
../build/sim_project ../build/run.mac      # 批处理模式 -> output/hits.csv

cd /root/geant4/sim_project/build
./sim_project                              # 不带宏参数 -> 交互式 Qt 可视化界面（运行 init_vis.mac）
```

后处理生成能谱图：`cd output && python3 plot_spectrum.py`
（读取 `hits.csv`，输出 `spectrum.png`）。

## 架构

标准的 Geant4 应用布局。`sim_project.cc` 通过三次 `SetUserInitialization` 把
对象注入 run manager；其余都是 Geant4 的必需类或可选用户类。所有类都位于
`namespace sim` 中。

- **`DetectorConstruction`** —— 几何：一个 20 cm 的空气立方体世界，以及一块锗
  圆柱体（`G4Tubs`，r = 15 mm，h = 10 mm）。通过 `GetCrystalVolume()` 暴露晶体
  逻辑体。
- **`PhysicsList`** —— 只注册 `G4EmStandardPhysics_option4`（不含衰变物理，因为
  源是普通伽马）。默认产生阈（production cut）0.1 mm。
- **`PrimaryGeneratorAction`** —— 一个纯粹的 `G4GeneralParticleSource`；粒子、
  能量、位置和方向完全由宏驱动（`run.mac` 中的 `/gps/...`）。
- **`ActionInitialization::Build()`** —— 构造每次运行所需的用户 action，并把
  `RunAction` 指针传给 `SteppingAction`，使两者共用同一个输出流。
- **`SteppingAction`** —— 打分（scoring）。对晶体体内每一个 `edep > 0` 的步，把
  步中点位置和能量转发给 `RunAction::RecordHit`。晶体体指针在首次使用时惰性
  解析。
- **`RunAction`** —— 持有 `hits.csv` 的 `ofstream`；在运行开始时打开并写表头，
  运行结束时关闭并打印 hit 计数。

### 坐标系 —— 关键的跨项目约定

晶体被放置成在**世界坐标 z [0, 10] mm** 范围内（`G4Tubs` 以自身原点为中心，
故放置在 `z = +height/2`）。这是有意为之：Geant4 世界坐标 `(x, y, z)`（单位
mm）可以**直接**映射到 SSD 的 `CartesianPoint(x, y, z) * 1e-3`（单位 m）——
无需任何变换。晶体半径、高度以及这个 z 坐标系必须与 `~/ssd_projects/ppc` 和
`ppc_2` 中的 `config.yaml` 保持一致。在此处修改几何而不同步更新那些 SSD 配置，
会悄无声息地破坏数据交接。

### `hits.csv` 格式

每个有能量沉积的步占一行：`event,x_mm,y_mm,z_mm,edep_keV`。`event` id 对应
一个初级伽马；共享同一 `event` 的所有行属于同一个瞬时探测器脉冲（按 `event`
分组以构建一个 SSD `Event`）。

## 约定与注意事项

- run manager **有意设为串行**（`G4RunManagerType::Serial`），使一个
  `RunAction` 独占一个输出流，`hits.csv` 不会被交错写入。改为多线程模式需要
  新的输出策略（按线程分文件、加互斥锁，或用 `G4AnalysisManager`）。
- `init_vis.mac`、`vis.mac`、`run.mac` 由 CMake 复制进 `build/`
  （`configure_file ... COPYONLY`）。修改顶层 `.mac` 后，需重新运行 `make`
  （或 `cmake`）来刷新 `build/` 中的副本；可执行文件运行的是 build 目录下的
  副本。
- 新的几何（低温恒温器、铝端窗、n+ 死层）应作为额外的体加入
  `DetectorConstruction.cc`——它们不会改变晶体坐标系。
- 若要模拟放射性同位素源而非普通伽马，需在 `PhysicsList.cc` 中注册
  `G4DecayPhysics` + `G4RadioactiveDecayPhysics`。
- `src/*.cc` 和 `include/*.hh` 由 CMake 的 `file(GLOB ...)` 收集，因此新增源
  文件需要重新运行一次 `cmake` 才会被编译。
