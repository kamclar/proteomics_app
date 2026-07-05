@echo off
setlocal enabledelayedexpansion

cd /d "%~dp0"


echo Hledam R...
echo.

if exist "R-runtime\bin\Rscript.exe" (
  echo Nalezen: R-runtime\bin\Rscript.exe
  echo.
  
  "R-runtime\bin\Rscript.exe" install_impute.R
  
  if errorlevel 1 (
    echo.
    echo !!! CHYBA pri instalaci !!!
  ) else (
    echo.
    echo ===== INSTALACE USPESNA =====
  )
) else (
  echo !!! CHYBA: R nenalezeno !!!
  echo Ocekavano: R-runtime\bin\Rscript.exe
)

echo.
echo Stiskni Enter pro zavreni...
pause >nul