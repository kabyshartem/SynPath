#!/bin/bash
pip install -r requirements.txt
mkdir -p static
cd static
curl -L "https://github.com/epam/ketcher/releases/download/v2.12.0/ketcher-standalone-2.12.0.zip" -o ketcher.zip
unzip -o ketcher.zip
rm ketcher.zip
cd ..
echo "Build complete"
