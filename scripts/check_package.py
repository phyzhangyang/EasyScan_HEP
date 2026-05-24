#!/usr/bin/env python3
"""Build and smoke-test the EasyScan_HEP package before release."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], *, cwd: Path = ROOT, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="easyscan-package-") as tmp:
        tmp_path = Path(tmp)
        dist_dir = tmp_path / "dist"
        target_dir = tmp_path / "install"
        egg_info = ROOT / "src" / "easyscan_hep.egg-info"

        if egg_info.exists():
            shutil.rmtree(egg_info)

        run([sys.executable, "-m", "build", "--outdir", str(dist_dir)])
        artifacts = sorted(glob.glob(str(dist_dir / "*")))
        run([sys.executable, "-m", "twine", "check", *artifacts])

        wheels = sorted(dist_dir.glob("*.whl"))
        if not wheels:
            raise RuntimeError("No wheel was built.")

        run([sys.executable, "-m", "pip", "install", "--target", str(target_dir), str(wheels[-1])])

        executable = target_dir / "bin" / "easyscan"
        if not executable.exists():
            executable = target_dir / "Scripts" / "easyscan.exe"
        if not executable.exists():
            raise RuntimeError("Installed easyscan console script was not found.")

        installed_example = target_dir / "templates" / "example_random.ini"
        if not installed_example.exists():
            raise RuntimeError("Installed example template was not found in the wheel.")

        env = os.environ.copy()
        env["PYTHONPATH"] = str(target_dir)
        run([str(executable), "--check", "templates/example_random.ini"], cwd=target_dir, env=env)

    print("Package check passed.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
