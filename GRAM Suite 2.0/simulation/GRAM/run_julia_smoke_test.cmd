@echo off
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_julia_smoke_test.ps1" %*
