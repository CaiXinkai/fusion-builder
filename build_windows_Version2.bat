@echo off
REM Windows 打包脚本（在含 Python 环境的 cmd 中运行）
REM 请先激活你的虚拟环境并确保已安装 pyinstaller: pip install pyinstaller
REM 这个脚本将生成一个单文件可执行文件（不包含大基因组文件——建议在运行时由用户选择本地 FASTA/GTF）

pyinstaller --onefile --windowed --name FusionBuilder app.py
echo Done. 输出在 dist\FusionBuilder.exe
pause