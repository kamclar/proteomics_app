@echo off
setlocal

set "APP_DIR=%~dp0"
cd /d "%APP_DIR%"

set "R_SCRIPT="

if exist "%APP_DIR%R-runtime\bin\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R-runtime\bin\Rscript.exe"
) else if exist "%APP_DIR%R-runtime\bin\x64\Rscript.exe" (
  set "R_SCRIPT=%APP_DIR%R-runtime\bin\x64\Rscript.exe"
) else if exist "%ProgramFiles%\R\R-4.5.2\bin\Rscript.exe" (
  set "R_SCRIPT=%ProgramFiles%\R\R-4.5.2\bin\Rscript.exe"
) else (
  for /f "delims=" %%R in ('where Rscript.exe 2^>nul') do (
    if not defined R_SCRIPT set "R_SCRIPT=%%R"
  )
)

if not defined R_SCRIPT (
  echo Could not find Rscript.exe.
  echo Put a portable R runtime in "%APP_DIR%R-runtime" or install R and add it to PATH.
  pause
  exit /b 1
)

if not exist "%APP_DIR%packages" (
  echo The local packages folder was not found.
  echo Run prepare_standalone_packages.R before distributing the app.
  pause
  exit /b 1
)

set "PROTEOMICS_APP_R_LIB=%APPDATA%\ProteomicsApp\R\win-library"
set "R_LIBS_USER=%APP_DIR%packages;%PROTEOMICS_APP_R_LIB%"

echo Starting ProteomicsApp...
echo Using Rscript: %R_SCRIPT%

"%R_SCRIPT%" -e "shiny::runApp('.', launch.browser=TRUE, port=3838)"

if errorlevel 1 (
  echo.
  echo ProteomicsApp stopped with an error.
  pause
)

endlocal
