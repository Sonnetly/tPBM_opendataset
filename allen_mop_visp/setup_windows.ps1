$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
py -3.12 -m venv .venv
if ($LASTEXITCODE -ne 0) { throw 'Python 3.12 venv creation failed' }
& "$PSScriptRoot\.venv\Scripts\python.exe" -m pip install -r requirements.txt
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
& "$PSScriptRoot\.venv\Scripts\python.exe" -m pip check
if ($LASTEXITCODE -ne 0) { throw 'Dependency check failed' }
Write-Host 'Environment prepared. See README_ko.md for analysis commands.'
