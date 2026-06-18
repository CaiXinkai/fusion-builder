# Fusion Transcript Builder (hg38)

这是一个用于构建基因融合转录本并翻译蛋白序列的桌面小工具（Python + PyQt5），支持单条和批量（CSV）输入。断点格式为 `chr:pos`（例如 `chr1:1234567`），参考基因组使用 hg38。

主要功能
- 输入 leftGene, rightGene, leftBreakpoint, rightBreakpoint，得到融合转录本核苷酸序列（在连接处用 `/` 标注）和对应蛋白序列（在氨基酸边界插入 `/`；若有 frameshift 会在 notes 中说明）。
- 支持批量 CSV 输入（见 sample_input.csv）。
- 支持导出结果为 TSV。
- GUI 基于 PyQt5，可在 Windows 下打包为可执行文件。

依赖
- Python 3.8+
- see requirements.txt

安装
1. 建议先创建虚拟环境：
   - python -m venv venv
   - source venv/bin/activate（Linux/macOS） 或 venv\Scripts\activate（Windows）

2. 安装依赖：
   - pip install -r requirements.txt

准备参考数据（你需要提供）
- hg38 参考 FASTA（建议文件名： hg38.fa 或 hg38.fa.gz）
- 对应的 GTF（例如 GENCODE/Ensembl 的 hg38 GTF）

运行
- 运行 GUI：
  python app.py
- 在 GUI 中选择参考 FASTA 和 GTF，输入单条或批量 CSV（CSV 示例见 sample_input.csv），点击 Run，最后可导出 TSV。

CSV 格式（批量）
CSV 必须包含列：leftgene,rightgene,leftbp,rightbp（列名不区分大小写）。
示例：
leftgene,rightgene,leftbp,rightbp
TMPRSS2,ERG,chr21:42880000,chr21:39750000

打包为 Windows 可执行（可选）
- 使用 PyInstaller：
  - Windows (命令行):
    build_windows.bat
  - Linux/macOS:
    ./build_unix.sh
- 打包说明在仓库内的 build_windows.bat / build_unix.sh 文件。

注意与局限
- 当前默认每个基因选择“最长转录本”作为代表转录本；如果你需要指定 transcript_id，请告诉我，我可以修改成优先使用用户输入或在 GUI 中选择。
- 断点为基因组坐标 chr:pos；若断点实际上是转录本坐标，需要调整解析逻辑。
- 如果断点位于内含子或外显子边界，程序会尽力决定左右序列并在 notes 中给出提示。
- 本程序不包含对 NMD、启动子/终止密码子复杂生物学规则的深度处理；如需更严格的 ORF 判定、同义/非同义注释，可后续扩展。

License: MIT
