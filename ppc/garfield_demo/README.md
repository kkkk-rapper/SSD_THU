# Garfield++ 信号仿真最小示例

加载 SSD 导出的电场和 P+ 权势，跑 122 keV 单点漂移仿真，输出 P+ 上的诱导电荷波形，并和 SSD 的同位置结果对照。

## 前置条件

1. SSD 已经导出场图：`../output/garfield/efield_1000V.txt` + `wpot_Pplus.txt`
   (任务 06 跑完后会自动生成)
2. Garfield++ 已经编译安装到 `/root/gpp/install/`
3. ROOT 在 PATH 里 (`root-config --version` 不报错)

## 编译运行

```bash
cd /root/ssd_projects/ppc/garfield_demo

# 设置 Garfield++ 环境变量
source /root/gpp/install/share/Garfield/setupGarfield.sh

# 编译
make

# 运行
./signal_demo
# 输出: signal_garfield.csv

# 跟 SSD 输出对照画图
/root/.juliaup/bin/julia --project=../ compare_with_ssd.jl
# 输出: compare_ssd_garfield.png
```

## 已知差异/警告

| 项 | SSD | Garfield++ Demo |
|---|---|---|
| 半导体材料 | Ge (ADL 模型) | **MediumSilicon 占位**，迁移率手动设成 Ge 78K 近似 |
| 漂移模型 | ADL 各向异性 | 各向同性 (Canali 模型) |
| 电荷云 | 单点 / N-body | 单点 |
| 时间步长 | 1 ns | 50 ps (自适应 RKF) |
| 信号计算 | Shockley-Ramo + 内置场插值 | Shockley-Ramo + ComponentGrid 三线性插值 |
| 单事件总电荷 | 41356 e⁻ (固定) | 41216 e⁻ (固定 -> 取决于设定) |

**预期对照结果**：两条 Q(t) 曲线应该在**总电荷**和**漂移时间数量级**上吻合 (~100 ns 量级)，但形状会有差异，因为：
- Si 迁移率 + Ge 数据混搭，速度不完全对
- 各向异性 vs 各向同性
- 数值插值方案不同

## 文件清单

- `signal_demo.cc` — C++ 主程序
- `Makefile`       — 构建脚本
- `compare_with_ssd.jl` — Julia 对比/画图脚本
- `README.md`      — 本文件
