#!/bin/sh
set -eu
cd "$(dirname "$0")"
if [ ! -x .app-venv/bin/python ]; then
  echo "First launch: installing ChemGridMap into an isolated local environment."
  python3 -m venv .app-venv
fi
if [ ! -f .app-venv/chemgridmap-ready ]; then
  .app-venv/bin/python -m pip install '.[umap]'
  touch .app-venv/chemgridmap-ready
fi
exec .app-venv/bin/python -m chemgridmap.app
