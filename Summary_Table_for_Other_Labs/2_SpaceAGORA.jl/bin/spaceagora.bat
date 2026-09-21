@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "REPO_ROOT=%%~fI"

if "%JULIA_PROJECT%"=="" (
  set "PROJECT_DIR=%REPO_ROOT%"
) else (
  set "PROJECT_DIR=%JULIA_PROJECT%"
)

julia --project="%PROJECT_DIR%" "%REPO_ROOT%\src\cli\main.jl" %*
