$ErrorActionPreference = 'Stop'
$desktopDir = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

Push-Location $desktopDir
try {
    npm run build:native
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    npm test
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    npm run pack:win
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} finally {
    Pop-Location
}
