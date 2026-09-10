"""Build a self-contained folder/app on the current OS; no cross-compilation."""
from pathlib import Path
import hashlib
import platform
import shutil
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    # Numba needs Python sources for cached UMAP kernels in a frozen distribution.
    args = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "--onedir",
            "--name", "ChemGridMap", "--paths", str(root / "src"),
            "--additional-hooks-dir", str(root / "scripts" / "pyinstaller_hooks"),
            "--collect-data", "chemgridmap", "--collect-all", "chembl_structure_pipeline",
            "--collect-data", "rdkit", "--collect-submodules", "rdkit.Chem",
            "--collect-all", "umap", "--collect-all", "pynndescent",
            "--hidden-import", "matplotlib.backends.backend_pdf",
            "--recursive-copy-metadata", "chemgridmap", "--recursive-copy-metadata", "umap-learn",
            "--add-data", str(root / "LICENSE") + (";" if sys.platform == "win32" else ":") + ".",
            "--exclude-module", "tkinter", "--exclude-module", "torch", "--exclude-module", "tensorflow",
            "--distpath", str(root / "dist" / "desktop"),
            "--workpath", str(root / "build" / "desktop"),
            "--specpath", str(root / "build")]
    if sys.platform == "darwin":
        args += ["--windowed", "--osx-bundle-identifier", "org.chemgridmap.local"]
    elif sys.platform == "win32":
        args += ["--noconsole"]
    args.append(str(root / "scripts" / "desktop_entry.py"))
    subprocess.run(args, cwd=root, check=True)
    from chemgridmap import __version__
    base = root / "dist" / f"ChemGridMap-{__version__}-{platform.system()}-{platform.machine()}"
    desktop = root / "dist" / "desktop"
    if sys.platform == "darwin":
        archive = Path(str(base) + ".zip")
        subprocess.run(["ditto", "-c", "-k", "--sequesterRsrc", "--keepParent",
                        str(desktop / "ChemGridMap.app"), str(archive)], check=True)
    else:
        archive = Path(shutil.make_archive(str(base), "zip" if sys.platform == "win32" else "gztar",
                                          root_dir=desktop, base_dir="ChemGridMap"))
    Path(str(archive) + ".sha256").write_text(
        hashlib.sha256(archive.read_bytes()).hexdigest() + "  " + archive.name + "\n", encoding="ascii")
    print(f"Archive: {archive}")
    print(f"Built for {platform.system()} {platform.machine()}. Test on the target platform before distributing.")


if __name__ == "__main__":
    main()
