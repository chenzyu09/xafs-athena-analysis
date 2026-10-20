# xafs-athena-analysis

一个 Claude Skill：分析 Athena 工程文件（`.prj`）中的 XAFS 数据，并输出论文风格的图和 Origin 作图脚本。

A Claude Skill for analyzing XAFS data stored in Athena (Demeter) project files: XANES oxidation-state analysis, linear combination fitting (LCF), FT-EXAFS, publication-style figures, and an Origin (originpro) plotting script.

## 功能

- 读取 `.prj`（gzip 压缩的 Athena 工程文件），基于 [xraylarch](https://github.com/xraypy/xraylarch)
- 数据质量检查：金属箔能量标定、可疑的能量平移、重复能量点、弱信号、样品与金属箔过于相似（可能测到基底）等
- 归一化、边位置/白线分析、LCF、AUTOBK 扣背景与傅里叶变换
- 论文风格图（FT 堆叠 + 阴影带标注壳层，XANES 叠加），PNG/PDF/SVG
- 生成可在 Origin 2021+ 内置 Python 中运行的作图脚本（数据嵌入，单文件）

## 安装

1. 下载本仓库中的 `xafs-athena-analysis` 文件夹，压缩成 ZIP（ZIP 内顶层应是该文件夹）。
2. Claude → **Customize → Skills** → **+** → **Create skill** → **Upload a skill**，上传 ZIP。
3. 之后发送 `athena.prj` 并说明需求（如“分析一下”“画成 Origin 的图”）即可。

## 说明与局限

- 分析结果为定性/半定量，配位数与键长需在 Artemis 等软件中做 EXAFS 拟合。
- Origin 脚本中部分格式命令（LabTalk）在不同 Origin 版本中的支持程度可能不同；脚本对格式命令做了容错，失败时只打印提示、不中断作图。
- 本仓库不包含任何实验数据。

## License

MIT
