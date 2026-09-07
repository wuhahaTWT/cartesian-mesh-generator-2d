#!/bin/sh
set -eu
PROJECT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
APP_PATH="$PROJECT_DIR/desktop/dist/mac-arm64/CartMesh2D.app"
if [ -d "$APP_PATH" ]; then
  exec /usr/bin/open "$APP_PATH"
fi
printf '%s\n' '还没有当前版本的 App。请在终端执行：' \
  "cd \"$PROJECT_DIR/desktop\"" 'npm ci' 'sh scripts/build-macos.sh' \
  '完成后再双击本文件。按回车关闭。'
read -r answer
