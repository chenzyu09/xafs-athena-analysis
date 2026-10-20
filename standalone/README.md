# xafs_plot —— Athena .prj 一键分析与作图（不需要 Claude）

把 Athena 工程文件（`.prj`）里的 XAFS 数据一键处理成：

- **论文风格图**（每个元素一张）：左边是 FT-EXAFS（k³ 加权，曲线堆叠，阴影带标注 M–O / M–M 壳层），右边是 XANES 叠加图。输出 PNG（600 dpi）和 PDF。
- **Origin 作图脚本** `XAFS_origin_plot.py`：数据已嵌入，在 Origin 2021 及以上版本中运行即可生成同样的图。默认格式为 Arial 36 号、边框 3 pt，下轴刻度朝外，其余三边无刻度。
- **Excel 数据** `XAFS_data.xlsx`：可以直接导入 Origin 或其他软件。
- **检查报告** `report.txt`：
  - 金属箔能量标定是否偏移；
  - Athena 中可疑的能量平移、重复能量点、弱信号；
  - 样品是否与金属箔几乎重合（可能测到的是基底）；
  - 吸收边位置、白线强度、按标样插值的价态、线性组合拟合（LCF）、FT 峰位。

程序会自动识别元素，支持常见元素的 K 边和 W/Ir/Pt/Au 等的 L3 边。它还会自动区分样品、氧化物标样和金属箔。

## 1. 安装（只需一次）

1. 安装 Python，推荐 3.11 或 3.12：<https://www.python.org/downloads/>
   Windows 安装时请勾选 **Add python.exe to PATH**。
2. 下载本文件夹：`xafs_plot.py`、`origin_template.py`、`requirements.txt`，两个 .py 文件要放在同一个文件夹里。
3. 在这个文件夹里打开命令行（在文件夹地址栏输入 `cmd` 回车），运行：
   ```
   pip install -r requirements.txt
   ```

## 2. 使用

```
python xafs_plot.py "D:\数据\athena.prj" --sample-name "FeCoNiRu-Cu"
```

- `--sample-name`：可选，把所有样品统一显示为这个名字；不加时使用 Athena 里的谱名。
- 结果在 .prj 旁边的 `athena_xafs_output` 文件夹里。

### 修改设置后重画

第一次运行会在 .prj 旁边生成 `athena_config.json`，用记事本就能打开修改：

| 想改什么 | 改哪里 |
|---|---|
| 某条谱的显示名 | `spectra` 里的 `name` |
| 是样品还是标样 | `role`：`sample` / `reference` / `foil` |
| 不画某条谱 | `include` 改为 `false` |
| 曲线颜色 | `color`（如 `"#D73027"`） |
| FT 图从上到下的顺序 | 调整 `spectra` 里各条的先后顺序 |
| k 范围 / 加权 | `elements` → 元素 → `kmin`、`kmax`（`"auto"` 或数字）、`kweight` |
| XANES 横轴范围 | `xanes_range`，如 `[7700, 7750]` |
| 阴影带 | `bands`：`"auto"`，或 `[[1.5, "#FEE090"], [2.7, "#C8DFED"]]`，`[]` 表示不画 |
| 字体与字号 | `style` |

改完后再运行一次同样的命令即可，程序会自动读取这个配置文件。

### 在 Origin 里出图

打开 Origin，在 **Python 控制台**（提示符 `>>>`）里运行下面这行，括号里换成你的实际路径：
```python
exec(open(r"D:\数据\athena_xafs_output\XAFS_origin_plot.py", encoding="utf-8").read())
```
也可以在 **命令窗口**（Window → Command Window）输入 `run -pyf "完整路径"`。
重新运行前建议先 File → New 新建项目。如果某项格式在你的 Origin 版本中不支持，输出窗口会提示“[格式设置跳过]”，但图仍会画出来。

## 3. 注意

- 程序只做固定流程的处理和检查，**不会替你下结论**。`report.txt` 里 ⚠ 开头的条目请重点核实。
- 价态插值和 LCF 的结果取决于标样。标样和样品如果不是同一批测的，能量可能有约 ±1 eV 的误差。
- 配位数和键长需要在 Artemis 等软件中做 EXAFS 拟合。
- FT 图中每条曲线都归一化到最大值等于 1，方便比较峰位，峰强请看 Excel 里的原始 |χ(R)| 或 report 里的强度比。

## 安装问题

如果 `pip install xraylarch` 在安装 `bibtexparser` 时失败，可以改为：
```
pip install --no-deps xraylarch
pip install numpy scipy matplotlib openpyxl lmfit asteval uncertainties pyshortcuts xraydb silx h5py peakutils pyyaml toml termcolor dill charset-normalizer packaging
```
