param(
    [Parameter(Position = 0)]
    [string]$GramRoot = "",
    [switch]$Clean
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

function Resolve-GramRoot {
    param(
        [string]$InputRoot,
        [string]$ScriptDirPath
    )

    if ($InputRoot) {
        return (Resolve-Path $InputRoot).Path
    }
    if ($env:GRAM_ROOT) {
        return (Resolve-Path $env:GRAM_ROOT).Path
    }

    $candidates = @(
        (Join-Path $ScriptDirPath "..\.."),
        (Join-Path $ScriptDirPath "..\..\GRAM Suite 2.0"),
        (Join-Path $ScriptDirPath "..\..\GRAM"),
        (Join-Path $ScriptDirPath "..\GRAM Suite 2.0"),
        (Join-Path $ScriptDirPath "..\GRAM")
    )

    foreach ($c in $candidates) {
        $full = [System.IO.Path]::GetFullPath($c)
        if ((Test-Path (Join-Path $full "Build")) -and (Test-Path (Join-Path $full "Julia"))) {
            return $full
        }
    }

    throw "Could not find GRAM root. Pass it as argument or set GRAM_ROOT."
}

function Get-MakeBin {
    $candidate = Get-Command mingw32-make -ErrorAction SilentlyContinue
    if ($candidate) { return $candidate.Source }

    $candidate = Get-Command make -ErrorAction SilentlyContinue
    if ($candidate) { return $candidate.Source }

    throw "Could not find mingw32-make or make on PATH. Install MSYS2 MinGW toolchain."
}

$GramRoot = Resolve-GramRoot -InputRoot $GramRoot -ScriptDirPath $ScriptDir
$BuildDir = Join-Path $GramRoot "Build"
if (-not (Test-Path $BuildDir)) {
    throw "Invalid GRAM root: $GramRoot"
}

$MakeBin = Get-MakeBin
Write-Host "GRAM_ROOT: $GramRoot"
Write-Host "MAKE_BIN:  $MakeBin"

$SetupScript = Join-Path $BuildDir "setup_cspice.sh"
$BashBin = Get-Command bash -ErrorAction SilentlyContinue
if ($BashBin) {
    & $BashBin.Source $SetupScript
    if ($LASTEXITCODE -ne 0) {
        throw "setup_cspice.sh failed with exit code $LASTEXITCODE"
    }
} else {
    $MingwArchive = Join-Path $GramRoot "common\cspice\lib\cspice_mingw64.a"
    if (-not (Test-Path $MingwArchive)) {
        throw "bash not found and bundled MinGW CSPICE archive missing: $MingwArchive"
    }
    Write-Host "Using bundled Windows CSPICE archive:"
    Write-Host "  $MingwArchive"
}

if ($Clean) {
    & $MakeBin -C $BuildDir clean
    if ($LASTEXITCODE -ne 0) {
        throw "Clean failed with exit code $LASTEXITCODE"
    }
}

$Jobs = 4
if ($env:NUMBER_OF_PROCESSORS) {
    $ParsedJobs = 0
    if ([int]::TryParse($env:NUMBER_OF_PROCESSORS, [ref]$ParsedJobs) -and $ParsedJobs -gt 0) {
        $Jobs = $ParsedJobs
    }
}

& $MakeBin -C $BuildDir shared "-j$Jobs"
if ($LASTEXITCODE -ne 0) {
    throw "Build failed with exit code $LASTEXITCODE"
}

$GramLib = Join-Path $GramRoot "Build\lib\libGRAM.dll"
if (-not (Test-Path $GramLib)) {
    throw "Build finished but shared library not found: $GramLib"
}

$EnvFile = Join-Path $ScriptDir "gram.env.ps1"
@"
`$env:GRAM_ROOT = "$GramRoot"
`$env:GRAM_LIB = "$GramLib"
"@ | Set-Content -Encoding UTF8 $EnvFile

Write-Host "Build complete."
Write-Host "GRAM_LIB: $GramLib"
Write-Host "Wrote:    $EnvFile"
