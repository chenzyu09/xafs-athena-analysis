---
name: "xafs-athena-analysis"
description: "分析 Athena 工程文件(.prj)里的 XAFS 数据（XANES 价态、LCF、FT-EXAFS），画论文风格图，并生成 Origin 内置 Python 作图脚本。用户发来 athena.prj 或要求做 XAS/XAFS 分析与作图时使用。"
---

# XAFS（Athena .prj）分析与作图

用户常把 XAFS 误称为“核磁”。.prj 是 gzip 压缩的 Demeter/Athena 工程文件，不是 NMR。先说明这一点再分析。

## 1. 环境与读取

- 安装 xraylarch：`pip install --break-system-packages xraylarch` 可能因 bibtexparser 编译失败。绕过：
  `pip install --break-system-packages --no-deps xraylarch` 然后
  `pip install --break-system-packages lmfit asteval uncertainties pyshortcuts xraydb silx h5py peakutils pyyaml toml termcolor dill charset-normalizer packaging matplotlib`
- 解压：`zcat athena.prj > athena.txt`（放在单独目录，Python 用 `-I` 运行）。
- 读取：`from larch.io import read_athena; p = read_athena('athena.txt', do_preedge=False, do_bkg=False, do_fft=False)`；
  遍历 `dir(p)` 中带 `athena_params` 的对象，`g.label` 为谱名，`g.athena_params.bkg.e0/eshift` 为 Athena 设定。
- 预处理：按能量排序，删除间隔 <0.05 eV 的重复能量点（会在导数谱上产生假尖峰）。

## 2. 质量与标定检查（必须做，并在回复中明确列出问题）

1. 用 `grep -o "'label','[^']*'"`、`'source'`、`'bkg_e0'`、`'bkg_eshift'` 列出所有谱、来源文件、E0、能量平移。
2. 金属箔一阶导最大值对照标准值：Fe 7112、Co 7709、Ni 8333、Cu 8979 eV。偏差 >1 eV 要指出；同一批次样品的标定可能受影响。
3. 任何非零 `bkg_eshift` 都要核实是否合理（例如 Co²⁺ 标样被平移到 7709 eV 就是错的）。
4. 检查边跳（edge_step），过小表示信号弱 → 限制 k 范围（如 3–8 Å⁻¹）。
5. 检查 source 文件名与 label 是否一致（如样品 source 叫 foil 要提醒核实）。
6. 样品若与金属箔几乎重合（相关系数 >0.998），提醒可能测到的是基底（如长在泡沫铜/铜片上透射测试），结论要据此限定。
7. 标样与样品来自不同线站/批次时，跨批次能量误差约 ±1 eV，小的边位移不要过度解读。

## 3. 分析

- 归一化：`pre_edge(g, e0=athena_e0, pre1=-150, pre2=-30, norm1=150, norm2=min(400, Emax-e0-10), nnorm=2)`，作图用 `g.flat`。
- 价态：一阶导峰位、半高边位置、白线强度；与 MO / M₂O₃ / 箔做边位置线性插值。
- LCF：在 E0−20 ~ E0+40 eV 用 `scipy.optimize.nnls` 对标样拟合，报告权重、总和与 R 因子。
- EXAFS：`autobk(g, rbkg=1.0, kweight=2或3)`；`xftf(g, kmin=3, kmax=..., dk=1, kweight=..., window='hanning')`。报告各壳层峰位（未相位校正）和样品/箔峰强比。
- 结论按元素写：价态、M–O / M–M 壳层、有无金属峰、可信度；再给整体结构图像，并列出写论文前需核实的问题。机理推断要标明是推断。

## 4. 论文风格图（matplotlib）

参考常见样式：左 FT-EXAFS（k³，堆叠，样品在最下，每条曲线归一到最大值=1 后加偏移，曲线名标在右侧上方；粉/黄阴影带标 M–O、M–M 壳层），右 XANES 叠加（图例左上无边框），y 轴无刻度数字（a.u.）。
- 若用户给参考图，核对其标签是否与曲线形状匹配（金属箔：边前台阶、上升平缓；MO：有边前峰、白线较强），发现标反要指出。
- 归一化堆叠会掩盖峰强差异，必须在图注写明，并可另给未归一化版本。
- 先用低 dpi 预览检查遮挡，再输出 600 dpi PNG + PDF/SVG。

## 5. Origin 版本（originpro 脚本）

生成单个 .py：处理后的数据以 JSON 字符串嵌入脚本（避免路径/依赖问题），只用标准库 + `originpro`。

运行方式告诉用户：
- Python 控制台（提示符 `>>>`）：`exec(open(r"C:\路径\XAFS_origin_plot.py", encoding="utf-8").read())`
- LabTalk 命令窗口：`run -pyf "C:\路径\XAFS_origin_plot.py"`（在 Python 控制台输入会报 unicodeescape 错）
- 重跑前 File → New，避免重名。

Origin 写法与坑（在 Origin 中实测过的标【实测】，其余为根据报错修正后的写法、尚待更多版本验证）：
- 数据：`wks = op.new_sheet('w', lname=...)`；先设 `wks.cols = N`，再 `wks.from_list(col, data, lname=, units=, comments=)`；`wks.cols_axis('XYXY...')`。
- 【实测】作图：`gp = op.new_graph(template='line'); gl = gp[0]; p = gl.add_plot(wks, coly=i, colx=j, type='l')`；`p.color = '#RRGGBB'`；线宽 `p.set_cmd('-w 1500')`（500 = 1 pt）。
- 所有格式命令用 `gl.lt_exec(...)` 并包 try/except，失败只打印“[格式设置跳过]”。
- 文字标注：【实测】`label -n name -a x y 文本` 会把 `-a x y` 当文字显示；改为 `label -n name 文本`（**不要加 `-a x y`**，会把参数当文字显示），再用 `name.x=右端-name.dx/2; name.y=底边+name.dy/2;` 定位（.x/.y 为中心、轴单位）；`name.font=font(Arial); name.fsize=36; name.background=0;`
- 阴影带：【实测】面积图（type=204）只改 p.color 时填充为黑色。改为用一条粗竖线：X=[b,b]、Y=[ylo,yhi]，线宽 = 带宽(Å) × 图层宽(pt)/6，先画以置于底层。
- 边框与刻度：`layer.x.showAxes=3; layer.y.showAxes=3;` `layer.x.thickness=3; layer.y.thickness=3;`（x2/y2 同设）；刻度 `layer.x.ticks=10`（主+次朝外），无刻度为 0；上/右轴 `layer.x2.ticks=0; layer.y2.ticks=0;`。
- 【实测】设了 `layer.x2.showLabels=0` 等后下轴数字消失，因此不要设；用 `layer.x.showLabels=1; layer.y.showLabels=0;`。
- 字体：`layer.x.label.font=font(Arial); layer.x.label.pt=36; xb.font=...; xb.fsize=36; yl...`；图例 `legend.font/fsize/background`，删除图例 `label -r legend;`。
- 页面：`page.width=10*page.resx; page.height=8*page.resy;`，`layer.unit=1; layer.left/top/width/height=%`。
- 大字号防遮挡：在 Python 中按 字宽≈0.62×字号、行高≈1.25×字号 估算文字尺寸，自动计算每条曲线偏移、标注高度（标注底边高于本曲线在标注 x 范围内的最大值，且下一条曲线让开）以及 XANES 的 y 上限（图例区域内曲线最大值低于图例底边）。大字号时刻度间隔放宽（FT 每 2 Å，XANES 每 20 eV，起点取整）。
- 导出：`gp.save_fig(path, width=3000)`；`op.save(path + '.opju')`。
- 本地无法运行 Origin：用假的 originpro 模块做空跑检查逻辑，并用 matplotlib 按相同页面/字号模拟排版检查遮挡；如实告诉用户未在 Origin 中实测。

配色：用户给参考配色时从图片图例取色（PIL getpixel）。很浅的颜色（如 #FEE090、#E0F3F8）不要画曲线，留给阴影带。

## 6. 交付

- 发送图（PNG/PDF）、可选 Excel 数据（前三行 Long Name / Units / Comments，数据从第 4 行起）和 Origin 脚本。
- 若连接了用户电脑且给了路径，直接把脚本写到该文件夹。
- 回复简洁：结论、需核实问题、下一步（EXAFS 拟合、补标样等）。