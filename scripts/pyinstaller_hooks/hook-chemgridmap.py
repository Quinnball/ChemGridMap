# Numba resolves kernel source paths when initializing the UMAP disk cache.
# Collect these modules as source, not code objects with archive-relative paths.
# https://pyinstaller.org/en/stable/hooks.html#hook-global-variables
module_collection_mode = {"umap": "py", "pynndescent": "py"}
