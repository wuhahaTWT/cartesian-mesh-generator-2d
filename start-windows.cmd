@echo off
setlocal
set "APP_PATH=%~dp0desktop\dist\win-unpacked\CartMesh2D.exe"
if exist "%APP_PATH%" (
  start "CartMesh2D" "%APP_PATH%" %*
  exit /b 0
)
echo Build the Windows app first:
echo npm ci --prefix desktop
echo npm --prefix desktop run build:native
echo npm --prefix desktop run pack:win
pause
exit /b 1
