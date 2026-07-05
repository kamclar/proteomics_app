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

set "HAS_CHANGES="
for /f "delims=" %%S in ('git status --short') do set "HAS_CHANGES=1"

if not defined HAS_CHANGES (
  echo No local changes to commit.
  pause
  exit /b 0
)

echo Local changes:
git status --short
echo.
echo Ignored files are not committed. Check .gitignore before continuing.
echo.

set "CONFIRM="
set /p CONFIRM=Commit and push these changes to GitHub main? [y/N]: 
if /i not "%CONFIRM%"=="y" (
  echo Cancelled.
  pause
  exit /b 0
)

set "COMMIT_MSG="
set /p COMMIT_MSG=Commit message [Update app]: 
if "%COMMIT_MSG%"=="" set "COMMIT_MSG=Update app"

git add -A
if errorlevel 1 goto error

git commit -m "%COMMIT_MSG%"
if errorlevel 1 goto error

git push origin HEAD:main
if errorlevel 1 goto error

echo.
echo Changes were pushed to GitHub.
pause
exit /b 0

:error
echo.
echo Commit or push failed. Check the message above.
pause
exit /b 1
