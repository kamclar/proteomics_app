@echo off
setlocal

cd /d "%~dp0"

where git >nul 2>nul
if errorlevel 1 (
  echo Git was not found. Install Git for Windows first.
  pause
  exit /b 1
)

git rev-parse --is-inside-work-tree >nul 2>nul
if errorlevel 1 (
  echo This folder is not a Git repository.
  pause
  exit /b 1
)

set "HAS_README_CHANGE="
for /f "delims=" %%S in ('git status --short README.md') do set "HAS_README_CHANGE=1"

if not defined HAS_README_CHANGE (
  echo README.md has no local changes to commit.
  pause
  exit /b 0
)

echo Changes in README.md:
git status --short README.md
echo.

set "COMMIT_MSG="
set /p COMMIT_MSG=Commit message [Update README]: 
if "%COMMIT_MSG%"=="" set "COMMIT_MSG=Update README"

git add README.md
if errorlevel 1 goto error

git commit -m "%COMMIT_MSG%"
if errorlevel 1 goto error

git push origin HEAD:main
if errorlevel 1 goto error

echo.
echo README.md was pushed to GitHub.
pause
exit /b 0

:error
echo.
echo Push failed. Check the message above.
pause
exit /b 1
