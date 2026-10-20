# xafs-athena-analysis

**Athena 工程文件（`.prj`）XAFS 数据的一键分析与论文作图工具**

One-click analysis and publication-quality plotting for XAFS data stored in Athena (Demeter) project files: data-quality checks, XANES edge/valence analysis, linear combination fitting (LCF), FT-EXAFS, publication-style figures, and a ready-to-run Origin script.

---

## 它能做什么

输入一个 Athena 的 `.prj` 文件，自动完成：

| 输出 | 内容 |
|---|---|
| 📊 **论文风格图**（PNG 600 dpi / PDF） | 每个元素一张：左边是 FT-EXAFS（k³ 加权，曲线堆叠，阴影带标注 M–O / M–M 壳层），右边是 XANES 叠加 |
| 📈 **Origin 作图脚本** | 数据已嵌入，在 Origin 2021+ 中运行即可生成同样的图（默认 Arial 36 号、边框 3 pt，各项格式可调） |
| 📄 **Excel 数据** | 前 3 行为 Long Name / Units / Comments，可直接导入 Origin |
| 🔍 **检查与分析报告** | 数据质量检查、吸收边位置、白线强度、价态插值、LCF、FT 峰位、峰强比 |

### 自动识别

- **元素与吸收边**：常见元素的 K 边（Ti–Ce），以及 W、Ir、Pt、Au 等的 L3 边。
- **样品 / 氧化物标样 / 金属箔**：根据谱名和化学式自动分类。
- **价态**：根据标样化学式计算（如 Fe₂O₃ → +3，Co(OH)₂ → +2）。

### 自动检查的常见问题

- 金属箔能量标定偏移（与标准边位置比较）。
- Athena 中设置了可疑的能量平移（平移过大的谱默认不画）。
- 重复能量点（会在导数谱上产生假尖峰，程序会自动删除）。
- 样品信号弱（会自动缩短 EXAFS 的 k 范围）。
- 样品与金属箔几乎重合（可能测到的是基底，比如长在泡沫铜上的样品）。
- 线性组合拟合（LCF）效果差（说明现有标样可能代表不了样品）。

---

## 两种使用方式

### 方式一：独立 Python 程序（推荐，不需要任何账号）

```
pip install -r standalone/requirements.txt
python standalone/xafs_plot.py athena.prj --sample-name "你的样品名"
```

第一次运行会生成一个配置文件。改里面的显示名、颜色、k 范围、阴影带等设置，再运行一次即可重画。

👉 详细说明见 **[standalone/README.md](standalone/README.md)**

### 方式二：Claude Skill（可选）

如果你使用 [Claude](https://claude.ai)，可以把 `xafs-athena-analysis/` 文件夹作为 Skill 安装。相比独立程序，Claude 除了出图，还能：

- 结合谱图解释检查出的问题，判断可能的原因；
- 按元素写出结构结论（价态、配位环境、有无金属相等）；
- 按你给的参考图或要求调整作图风格。

安装方法：
1. 把 `xafs-athena-analysis` 文件夹压缩成 ZIP（ZIP 内顶层应是该文件夹）。
2. Claude → **Customize → Skills** → **+** → **Create skill** → **Upload a skill**，上传 ZIP。
3. 发送 `.prj` 文件，说“分析一下”或“画成 Origin 的图”即可。

---

## 仓库结构

```
├── standalone/                独立 Python 程序
│   ├── xafs_plot.py           主程序
│   ├── origin_template.py     Origin 脚本模板
│   ├── requirements.txt
│   └── README.md              使用说明
├── xafs-athena-analysis/      Claude Skill
│   └── SKILL.md
├── LICENSE
└── README.md
```

## 说明与局限

- 结果为定性 / 半定量分析。配位数与键长需在 Artemis 等软件中做 EXAFS 拟合。
- 价态插值与 LCF 依赖标样。标样与样品不是同批测试时，能量误差约 ±1 eV。
- FT 图中每条曲线都归一化到最大值等于 1，方便比较峰位；比较峰强请看报告或 Excel 中的原始 |χ(R)|。
- Origin 脚本中部分格式命令（LabTalk）在不同 Origin 版本中的支持程度可能不同；脚本做了容错，某项格式失败时只打印提示，不会中断作图。
- 数据处理基于 [xraylarch](https://github.com/xraypy/xraylarch)。
- 本仓库不包含任何实验数据。

## License

[MIT](LICENSE)
