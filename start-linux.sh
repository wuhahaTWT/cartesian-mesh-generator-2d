#!/bin/sh
set -eu
PROJECT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
APP_PATH="$PROJECT_DIR/desktop/dist/linux-unpacked/cartmesh2d-desktop"
if [ -x "$APP_PATH" ]; then
  exec "$APP_PATH" "$@"
fi
printf '%s\n' 'Build the Linux app first:' \
  'npm ci --prefix desktop' \
  'npm --prefix desktop run build:native' \
  'npm --prefix desktop run pack:linux'
exit 1
