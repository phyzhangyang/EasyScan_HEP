#!/usr/bin/env python3
"""Prepare the fixed-version SSM dark-matter/EWPT EasyScan example."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from pathlib import Path


EXAMPLE_DIR = Path(__file__).resolve().parent
ROOT = EXAMPLE_DIR.parents[1]
PROGRAMS_DIR = EXAMPLE_DIR / "programs"
PATCH_DIR = EXAMPLE_DIR / "patches"

MICROMEGAS_VERSION = "7.1"
MICROMEGAS_URL = "https://zenodo.org/records/20267206/files/micromegas_7.1.tgz?download=1"
MICROMEGAS_DIR = PROGRAMS_DIR / "micromegas_7.1"

PHASETRACER_REF = "2.2.1"
PHASETRACER_COMMIT = "a68b2fd2801e748d143058d828683f8f76b67ce3"
PHASETRACER_URL_TEMPLATE = "https://github.com/PhaseTracer/PhaseTracer/archive/{commit}.tar.gz"
PHASETRACER_DIR = PROGRAMS_DIR / "PhaseTracer"

PHASETRACER_DEPENDENCY_HELP = """
PhaseTracer CMake configuration failed.

Please install the normal PhaseTracer build dependencies so CMake can find
them: Boost, Eigen3, NLopt and ALGLIB. On macOS with Homebrew, the commonly
available packages are:

  brew install boost eigen nlopt pkgconf

ALGLIB must also be installed in a location visible to CMake; use the
PhaseTracer/ALGLIB installation instructions appropriate for your system.
After installing the missing dependency, rerun this bootstrap script.
""".strip()


def run(command: list[str], *, cwd: Path | None = None, shell: bool = False) -> subprocess.CompletedProcess:
    display = command if not shell else command[0]
    print("+", display if isinstance(display, str) else " ".join(display), flush=True)
    return subprocess.run(command if not shell else command[0], cwd=cwd, shell=shell, check=True, text=True)


def run_optional(command: list[str], *, cwd: Path) -> int:
    return subprocess.run(command, cwd=cwd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True).returncode


def require_tools(tools: list[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise SystemExit("Missing required command(s): " + ", ".join(missing))


def safe_extract(archive: Path, dest: Path) -> None:
    dest = dest.resolve()
    with tarfile.open(archive) as handle:
        for member in handle.getmembers():
            target = (dest / member.name).resolve()
            if not str(target).startswith(str(dest)):
                raise RuntimeError(f"Refusing unsafe archive member: {member.name}")
        handle.extractall(dest)


def download_file(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"{dest} already exists; skipping download.", flush=True)
        return
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response, dest.open("wb") as output:
        shutil.copyfileobj(response, output)


def extract_single_directory(archive: Path, dest: Path) -> None:
    dest_parent = dest.parent
    dest_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".extract-", dir=dest_parent) as tmp:
        tmp_dir = Path(tmp)
        safe_extract(archive, tmp_dir)
        candidates = [path for path in tmp_dir.iterdir() if path.is_dir()]
        if len(candidates) != 1:
            raise SystemExit(f"Could not locate a single source directory in {archive}.")
        if dest.exists():
            shutil.rmtree(dest)
        candidates[0].rename(dest)


def apply_patch_file(patch_file: Path, cwd: Path) -> None:
    dry = ["patch", "-l", "--forward", "--batch", "--dry-run", "-p1", "-i", str(patch_file)]
    apply = ["patch", "-l", "--forward", "--batch", "-p1", "-i", str(patch_file)]
    reverse = ["patch", "-l", "--reverse", "--batch", "--dry-run", "-p1", "-i", str(patch_file)]
    if run_optional(dry, cwd=cwd) == 0:
        run(apply, cwd=cwd)
    elif run_optional(reverse, cwd=cwd) == 0:
        print(f"{patch_file.name} is already applied.", flush=True)
    else:
        raise SystemExit(f"Patch can not be applied cleanly: {patch_file}")


def prepare_micromegas(args: argparse.Namespace) -> None:
    if args.force and MICROMEGAS_DIR.exists():
        shutil.rmtree(MICROMEGAS_DIR)
    PROGRAMS_DIR.mkdir(parents=True, exist_ok=True)
    if not MICROMEGAS_DIR.exists():
        archive = Path(args.micromegas_tar).expanduser().resolve() if args.micromegas_tar else EXAMPLE_DIR / "micromegas_7.1.tgz"
        if args.micromegas_tar:
            if not archive.is_file():
                raise SystemExit(f"micrOMEGAs archive was not found: {archive}")
        else:
            download_file(args.micromegas_url, archive)
        safe_extract(archive, PROGRAMS_DIR)
        if not MICROMEGAS_DIR.is_dir():
            candidates = sorted(PROGRAMS_DIR.glob("micromegas*"))
            if len(candidates) == 1:
                candidates[0].rename(MICROMEGAS_DIR)
            else:
                raise SystemExit("Could not locate extracted micromegas_7.1 directory.")
        if not args.keep_archives and not args.micromegas_tar and archive.exists():
            archive.unlink()
    apply_patch_file(PATCH_DIR / "micromegas_7.1_easyscan.patch", MICROMEGAS_DIR)
    write_micromegas_input()
    # The top-level CalcHEP/micrOMEGAs build has archive-generation steps that
    # can race under parallel make on macOS, so keep this stage serial.
    run(["make"], cwd=MICROMEGAS_DIR)
    run(["make", f"-j{args.jobs}", "main=main.c"], cwd=MICROMEGAS_DIR / "SingletDM")


def prepare_phasetracer(args: argparse.Namespace) -> None:
    if args.force and PHASETRACER_DIR.exists():
        shutil.rmtree(PHASETRACER_DIR)
    PROGRAMS_DIR.mkdir(parents=True, exist_ok=True)
    if not PHASETRACER_DIR.exists():
        archive = (
            Path(args.phasetracer_tar).expanduser().resolve()
            if args.phasetracer_tar
            else EXAMPLE_DIR / f"PhaseTracer-{args.phasetracer_ref}-{args.phasetracer_commit[:12]}.tar.gz"
        )
        if args.phasetracer_tar:
            if not archive.is_file():
                raise SystemExit(f"PhaseTracer archive was not found: {archive}")
        else:
            download_file(args.phasetracer_url, archive)
        extract_single_directory(archive, PHASETRACER_DIR)
        if not args.keep_archives and not args.phasetracer_tar and archive.exists():
            archive.unlink()
    write_phasetracer_input()
    build_dir = PHASETRACER_DIR / "build-easyscan"
    if args.clean_build and build_dir.exists():
        shutil.rmtree(build_dir)
    cmake_cmd = [
        "cmake",
        "-S",
        str(PHASETRACER_DIR),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DBUILD_WITH_FS=OFF",
        "-DBUILD_WITH_FS_ScalarSingletZ2DM=OFF",
        "-DBUILD_WITH_FS_ScalarSingletZ2DMEWSBoutputlamHEFTHiggs=OFF",
        "-DBUILD_WITH_BSMPT=OFF",
        "-DBUILD_WITH_BP=OFF",
    ]
    try:
        run(cmake_cmd)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(PHASETRACER_DEPENDENCY_HELP) from exc
    run(["cmake", "--build", str(build_dir), "--target", "run_xSM_MSbar", "--parallel", str(args.jobs)])
    if args.clean_build and build_dir.exists():
        shutil.rmtree(build_dir)


def write_micromegas_input() -> None:
    path = MICROMEGAS_DIR / "SingletDM" / "data1.par"
    path.write_text("Q      100\nMh     125\nlaS    1\nlaSH   0.30\nMdm1   90\n", encoding="utf-8")


def write_phasetracer_input() -> None:
    path = PHASETRACER_DIR / "easyscan_phasetracer.in"
    path.write_text(
        "\n".join(
            [
                "ms 90",
                "lambda_s 1",
                "lambda_hs 0.60",
                "Q 173",
                "xi 1",
                "daisy_flag 2",
                "use_1L_EWSB_in_0L_mass 0",
                "use_Goldstone_resum 1",
                "",
            ]
        ),
        encoding="utf-8",
    )


def smoke_test() -> None:
    print("Running micrOMEGAs smoke test.", flush=True)
    mm_output = MICROMEGAS_DIR / "SingletDM" / "easyscan_micromegas.out"
    run(["./main data1.par > easyscan_micromegas.out"], cwd=MICROMEGAS_DIR / "SingletDM", shell=True)
    text = mm_output.read_text(encoding="utf-8", errors="replace")
    for label in ("ES_Omega_h2", "ES_sigmaSIp_pb"):
        if label not in text:
            raise SystemExit(f"micrOMEGAs smoke test did not produce {label}.")

    print("Running PhaseTracer smoke test.", flush=True)
    run(["./bin/run_xSM_MSbar $(awk '{print $2}' easyscan_phasetracer.in)"], cwd=PHASETRACER_DIR, shell=True)
    pt_output = PHASETRACER_DIR / "output.txt"
    fields = pt_output.read_text(encoding="utf-8").strip().split()
    if len(fields) < 14:
        raise SystemExit("PhaseTracer smoke test output.txt has too few columns.")
    if abs(float(fields[2]) - 0.60) > 1e-12:
        raise SystemExit("PhaseTracer smoke test did not use lambda_hs = 2 * lambdaHS.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--micromegas-url", default=MICROMEGAS_URL)
    parser.add_argument("--micromegas-tar", default=None, help="Use a local micromegas_7.1 archive instead of downloading.")
    parser.add_argument("--phasetracer-url", default=None, help="Download a fixed PhaseTracer source archive from this URL.")
    parser.add_argument("--phasetracer-tar", default=None, help="Use a local PhaseTracer source archive instead of downloading.")
    parser.add_argument("--phasetracer-ref", default=PHASETRACER_REF)
    parser.add_argument("--phasetracer-commit", default=PHASETRACER_COMMIT)
    parser.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--force", action="store_true", help="Remove existing downloaded program directories first.")
    parser.add_argument("--keep-archives", action="store_true")
    parser.add_argument("--clean-build", action="store_true", default=True, help="Delete the PhaseTracer build directory after compiling.")
    parser.add_argument("--no-clean-build", action="store_false", dest="clean_build")
    parser.add_argument("--no-smoke", action="store_true")
    args = parser.parse_args()
    if args.phasetracer_url is None:
        args.phasetracer_url = PHASETRACER_URL_TEMPLATE.format(commit=args.phasetracer_commit)
    return args


def main() -> int:
    args = parse_args()
    require_tools(["make", "cmake", "patch"])
    prepare_micromegas(args)
    prepare_phasetracer(args)
    if not args.no_smoke:
        smoke_test()
    print("\nSSM_DM_EWPT example is ready.", flush=True)
    print("Next commands:", flush=True)
    print("  python3 bin/easyscan.py --check templates/scan_SSM_DM_EWPT.ini --json", flush=True)
    print("  python3 bin/easyscan.py templates/scan_SSM_DM_EWPT.ini", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
