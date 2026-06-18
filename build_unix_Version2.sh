#!/bin/bash
# Unix/macOS 打包脚本（假如你要在 macOS/Linux 上打包）
# 需要安装 pyinstaller： pip install pyinstaller
pyinstaller --onefile --windowed --name FusionBuilder app.py
echo "Done. check ./dist/FusionBuilder"