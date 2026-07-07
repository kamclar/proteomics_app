@echo off
setlocal enabledelayedexpansion

set "APP_DIR=%~dp0"
cd /d "%APP_DIR%"

echo Checking R installation...
echo.

if exist "%APP_DIR%R-runtime\bin\Rscript.exe" (
  echo Found bundled R runtime:
  echo %APP_DIR%R-runtime\bin\Rscript.exe
  echo.
  echo No R installation is needed.
  pause
  exit /b 0
)

where winget >nul 2>nul
if errorlevel 1 (
  echo Local R-runtime was not found and winget is not available.
  echo Install R for Windows from:
  echo https://cloud.r-project.org/bin/windows/base/
  echo.
  echo Then run install_packages.bat and run_app.bat.
  pause
  exit /b 1
)

echo Local R-runtime was not found.
echo Installing R for Windows into:
echo %APP_DIR%R-runtime
echo.

winget install --id RProject.R -e --source winget --location "%APP_DIR%R-runtime" --accept-package-agreements --accept-source-agreements
if errorlevel 1 (
  echo.
  echo R installation failed.
  echo You can install R manually from:
  echo https://cloud.r-project.org/bin/windows/base/
  pause
  exit /b 1
)

if not exist "%APP_DIR%R-runtime\bin\Rscript.exe" (
  echo.
  echo R installation finished, but Rscript.exe was not found in R-runtime.
  echo Check whether winget installed R to a system location instead.
  pause
  exit /b 1
)

echo.
echo Local R-runtime installation finished.
echo Next steps:
echo 1. Run install_packages.bat
echo 2. Run run_app.bat
echo.
pause
exit /b 0
