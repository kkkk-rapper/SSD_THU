# =============================================================================
# 任务 03：从任务 02 的 CSV 数据生成 2D 电势/电场热图
#
# 输入  ：output/ 下任务 02 写出的 CSV 文件
# 输出  ：output/plots/ 下的 PNG 图
# 后端  ：GR（无头环境也能跑，不需要 X11）
# =============================================================================

using DelimitedFiles
using Plots
gr()   # 指定 GR 后端，纯像素输出

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR   = joinpath(PROJECT_ROOT, "output")
const PLOTS_DIR    = joinpath(OUTPUT_DIR, "plots")
isdir(PLOTS_DIR) || mkdir(PLOTS_DIR)

# ----- 1. 读 CSV 数据 -----
# 单列轴文件（CSV 一列） + 矩阵文件（CSV，行=r 方向, 列=z 方向）
r_mm = vec(readdlm(joinpath(OUTPUT_DIR, "r_axis_mm.csv"), ','))
z_mm = vec(readdlm(joinpath(OUTPUT_DIR, "z_axis_mm.csv"), ','))
pot  = readdlm(joinpath(OUTPUT_DIR, "potential_V_rz.csv"),     ',')
Er   = readdlm(joinpath(OUTPUT_DIR, "Efield_Er_Vpm_rz.csv"),   ',')
Ez   = readdlm(joinpath(OUTPUT_DIR, "Efield_Ez_Vpm_rz.csv"),   ',')
Em   = readdlm(joinpath(OUTPUT_DIR, "Efield_Emag_Vpm_rz.csv"), ',')

@info "读入数据" size_pot=size(pot) size_E=size(Em) Nr=length(r_mm) Nz=length(z_mm)

# ----- 2. 晶体/电极外形线（用于在图上勾出轮廓） -----
# 晶体：矩形 r ∈ [0,15], z ∈ [0,10]
crystal_r = [0.0, 15.0, 15.0, 0.0, 0.0]
crystal_z = [0.0,  0.0, 10.0, 10.0, 0.0]
# P+ 点电极：z=10 顶面上 r ∈ [0,4]
ppos_r = [0.0, 4.0]
ppos_z = [10.0, 10.0]
# N+ 包裹：底面 + 侧面（两段线）
npos_r_bot  = [0.0, 15.0];  npos_z_bot  = [0.0, 0.0]
npos_r_side = [15.0, 15.0]; npos_z_side = [0.0, 10.0]

# ----- 3. 通用热图绘制函数 -----
# heatmap 用 (x=z_mm, y=r_mm, z=矩阵) 配置：竖向 z，水平 r 不太直观
# PPC 文献里通常 x=r, y=z，所以我们做转置：x=r, y=z, color=value
# 数据矩阵 mat[r_index, z_index] -> 需要 transpose 给 heatmap (heatmap 期望 mat[y, x])
function plot_heatmap(mat, title_text, cbar_label;
                      fname::String, colormap = :viridis,
                      clims = nothing, logscale = false)

    data = logscale ? log10.(max.(mat, 1e-3)) : mat
    label = logscale ? "log10($cbar_label)" : cbar_label

    kwargs = Dict{Symbol,Any}(
        :xlabel => "r (mm)",
        :ylabel => "z (mm)",
        :title  => title_text,
        :colorbar_title => "  " * label,
        :aspect_ratio => :equal,
        :c => colormap,
        :size => (700, 550),
        :dpi => 150,
    )
    if !isnothing(clims)
        kwargs[:clims] = clims
    end

    # heatmap(x_vec, y_vec, Matrix) 期望 Matrix 的 size = (length(y_vec), length(x_vec))
    # 我们的数据 size = (Nr, Nz)，要画 x=r, y=z -> 需要 permutedims
    p = heatmap(r_mm, z_mm, permutedims(data); kwargs...)

    # 叠加轮廓
    plot!(p, crystal_r, crystal_z, lw=2, color=:white, label="")
    plot!(p, ppos_r,    ppos_z,    lw=4, color=:red,   label="P+")
    plot!(p, npos_r_bot,  npos_z_bot,  lw=4, color=:cyan, label="N+ bottom")
    plot!(p, npos_r_side, npos_z_side, lw=4, color=:cyan, label="")

    savefig(p, joinpath(PLOTS_DIR, fname))
    @info "保存图像" file=fname
    return p
end

# ----- 4. 生成 4 张图 -----
# (1) 电势
plot_heatmap(pot, "Electric Potential @ V_op = 1000 V", "V";
             fname = "potential_V.png", colormap = :viridis)

# (2) 电场模 |E|（kV/m，方便读数）
plot_heatmap(Em ./ 1000, "Electric Field |E| @ V_op = 1000 V", "kV/m";
             fname = "Efield_magnitude_kVpm.png", colormap = :inferno)

# (3) 电场模 |E|（对数刻度，能看清点电极尖端的高场区）
plot_heatmap(Em, "Electric Field |E| (log scale)", "V/m";
             fname = "Efield_magnitude_log.png", colormap = :inferno, logscale = true)

# (4) E_z 分量（带符号，发散色板）
# 手动限色标 ±200 kV/m：覆盖体内典型值 ~150 kV/m，
# P+ 尖端 ~1800 kV/m 区域会饱和成深色，但晶体内部的方向/大小细节看得清。
Ez_kVpm = Ez ./ 1000
plot_heatmap(Ez_kVpm, "E_z @ V_op = 1000 V (saturated at ±200 kV/m)", "kV/m";
             fname = "Efield_Ez_kVpm.png", colormap = :RdBu,
             clims = (-200, 200))

# (5) E_r 分量
Er_kVpm = Er ./ 1000
plot_heatmap(Er_kVpm, "E_r @ V_op = 1000 V (saturated at ±200 kV/m)", "kV/m";
             fname = "Efield_Er_kVpm.png", colormap = :RdBu,
             clims = (-200, 200))

# ----- 5. 一维切片：沿 P+ 中轴线 (r=0) 的电势/电场 -----
i0 = argmin(abs.(r_mm .- 0.0))
mask = (z_mm .>= 0) .& (z_mm .<= 10)
z_in = z_mm[mask]
pot_axis = pot[i0, :][mask]
Em_axis  = Em[i0,  :][mask] ./ 1000

p_axis = plot(z_in, pot_axis;
    xlabel="z (mm)  [r=0, P+ at z=10]",
    ylabel="Potential (V)",
    label="Potential",
    lw=2, color=:blue,
    title="Axial profile @ r = 0",
    size=(700, 450), dpi=150,
)
p_axis_E = twinx(p_axis)
plot!(p_axis_E, z_in, Em_axis;
    ylabel="|E| (kV/m)",
    label="|E|",
    lw=2, color=:red, linestyle=:dash,
)
savefig(p_axis, joinpath(PLOTS_DIR, "axial_profile_r0.png"))
@info "保存图像" file="axial_profile_r0.png"

println("完成。图像目录: $PLOTS_DIR")
