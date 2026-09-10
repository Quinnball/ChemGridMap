@echo off
cd /d "%~dp0"
if not exist .app-venv\Scripts\python.exe (
  echo First launch: installing ChemGridMap into an isolated local environment.
  py -3 -m venv .app-venv
  if errorlevel 1 goto failed
)
if not exist .app-venv\chemgridmap-ready (
  .app-venv\Scripts\python.exe -m pip install ".[umap]"
  if errorlevel 1 goto failed
  echo ready>.app-venv\chemgridmap-ready
)
.app-venv\Scripts\python.exe -m chemgridmap.app
if errorlevel 1 goto failed
exit /b 0
:failed
echo Could not start ChemGridMap. Install Python 3.10 or newer, then try again.
pause
exit /b 1
