param(
    [Parameter(Position = 0)]
    [string]$GramRoot = ""
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$EnvFile = Join-Path $ScriptDir "gram.env.ps1"

if (Test-Path $EnvFile) {
    . $EnvFile
}

if ($GramRoot) {
    $env:GRAM_ROOT = (Resolve-Path $GramRoot).Path
}

$JuliaBin = Get-Command julia -ErrorAction SilentlyContinue
if (-not $JuliaBin) {
    throw "Could not find julia on PATH."
}

$SmokeTest = Join-Path $ScriptDir "julia_smoke_test.jl"
& $JuliaBin.Source $SmokeTest
if ($LASTEXITCODE -ne 0) {
    throw "Smoke test failed with exit code $LASTEXITCODE"
}
