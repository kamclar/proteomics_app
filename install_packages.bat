@echo off
setlocal enabledelayedexpansion

set "APP_DIR=%~dp0"
cd /d "%APP_DIR%"


echo Hledam R...
echo.

set "R_SCRIPT="

if exist "%APP_DIR%R-runtime\bin\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R-runtime\bin\Rscript.exe"
) else if exist "%APP_DIR%R-runtime\bin\x64\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R-runtime\bin\x64\Rscript.exe"
) else (
  for /f "delims=" %%R in ('where Rscript.exe 2^>nul') do (
    if not defined R_SCRIPT set "R_SCRIPT=%%R"
  )
)

if not defined R_SCRIPT (
  for /f "delims=" %%D in ('dir /b /ad "%ProgramFiles%\R\R-*" 2^>nul') do (
    if exist "%ProgramFiles%\R\%%D\bin\Rscript.exe" set "R_SCRIPT=%ProgramFiles%\R\%%D\bin\Rscript.exe"
  )
)

if not defined R_SCRIPT (
  echo !!! CHYBA: R nenalezeno !!!
  echo Spust install_r_runtime.bat nebo nainstaluj R for Windows.
  pause
  exit /b 1
)

echo Nalezeno R:
echo %R_SCRIPT%
echo.

"%R_SCRIPT%" install_packages.R

if errorlevel 1 (
  echo.
  echo !!! CHYBA pri instalaci !!!
) else (
  echo.
  echo ===== INSTALACE USPESNA =====
)

echo.
echo Stiskni Enter pro zavreni...
pause >nul
