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

set "R_SCRIPT="
for /f "delims=" %%R in ('where Rscript.exe 2^>nul') do (
  if not defined R_SCRIPT set "R_SCRIPT=%%R"
)

if not defined R_SCRIPT (
  for /f "delims=" %%D in ('dir /b /ad "%ProgramFiles%\R\R-*" 2^>nul') do (
    if exist "%ProgramFiles%\R\%%D\bin\Rscript.exe" set "R_SCRIPT=%ProgramFiles%\R\%%D\bin\Rscript.exe"
  )
)

if defined R_SCRIPT (
  echo Found installed R:
  echo %R_SCRIPT%
  echo.
  echo No R installation is needed.
  pause
  exit /b 0
)

where winget >nul 2>nul
if errorlevel 1 (
  echo R was not found and winget is not available.
  echo Install R for Windows from:
  echo https://cloud.r-project.org/bin/windows/base/
  echo.
  echo Then run install_packages.bat and run_app.bat.
  pause
  exit /b 1
)

echo R was not found.
echo Installing R for Windows with winget...
echo.

winget install --id RProject.R -e --source winget --accept-package-agreements --accept-source-agreements
if errorlevel 1 (
  echo.
  echo R installation failed.
  echo You can install R manually from:
  echo https://cloud.r-project.org/bin/windows/base/
  pause
  exit /b 1
)

echo.
echo R installation finished.
echo Next steps:
echo 1. Run install_packages.bat
echo 2. Run run_app.bat
echo.
pause
exit /b 0
