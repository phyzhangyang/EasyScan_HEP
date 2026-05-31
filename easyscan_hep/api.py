"""Agent-friendly public API for EasyScan_HEP."""

from __future__ import annotations

import configparser
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from .config_io import resolve_path
from .results import read_results


VALID_OVERWRITE_ACTIONS = {"replace", "backup", "stop"}
AGENT_SKILL_REPOSITORY = "https://github.com/PhenoAgent/easyscan-skill.git"


def source_root() -> Path:
    return Path(__file__).resolve().parents[1]


def configure_import_paths() -> Path:
    root = source_root()
    src_root = root / "src"
    if src_root.is_dir() and str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))
    if root.is_dir() and str(root) not in sys.path:
        sys.path.insert(0, str(root))
    os.environ.setdefault("EASYSCAN_ROOT", str(root))
    return root


def check_config(path: str | Path, base_dir: str | Path | None = None) -> dict[str, Any]:
    configure_import_paths()
    from config_checker import check_config_file

    config_path = Path(path).expanduser()
    if base_dir is not None and not config_path.is_absolute():
        config_path = Path(base_dir).expanduser().resolve() / config_path
    return check_config_file(config_path, base_dir=base_dir)


def result_dir_from_config(config_path: str | Path, base_dir: str | Path | None = None) -> Path:
    config_path = Path(config_path).expanduser().resolve()
    launch_dir = Path(base_dir or Path.cwd()).expanduser().resolve()
    config = configparser.ConfigParser(interpolation=None)
    config.read(config_path, encoding="utf-8")
    if not config.has_option("scan", "Result folder name"):
        return launch_dir / "result"
    return resolve_path(config.get("scan", "Result folder name").strip(), launch_dir)


def run_config(
    config_path: str | Path,
    *,
    cwd: str | Path | None = None,
    overwrite: str | None = "stop",
    log_path: str | Path | None = None,
    timeout: float | None = None,
) -> dict[str, Any]:
    """Run a config through the CLI and return a structured report."""
    configure_import_paths()
    launch_dir = Path(cwd or Path.cwd()).expanduser().resolve()
    config_path = Path(config_path).expanduser()
    if not config_path.is_absolute():
        config_path = launch_dir / config_path
    config_path = config_path.resolve()
    if overwrite is not None and overwrite not in VALID_OVERWRITE_ACTIONS:
        raise ValueError('overwrite must be one of "replace", "backup", "stop", or None')
    if log_path is None:
        log_path = config_path.with_suffix(".log")
    log_path = Path(log_path).expanduser()
    if not log_path.is_absolute():
        log_path = launch_dir / log_path
    log_path = log_path.resolve()
    result_dir = result_dir_from_config(config_path, launch_dir)

    command = [sys.executable, "-m", "easyscan_hep.cli", str(config_path)]
    env = os.environ.copy()
    pythonpath = [str(source_root()), str(source_root() / "src")]
    if env.get("PYTHONPATH"):
        pythonpath.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(pythonpath)
    if overwrite:
        env["EASYSCAN_RESULT_EXISTS_ACTION"] = overwrite

    completed = subprocess.run(
        command,
        cwd=launch_dir,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=timeout,
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(completed.stdout, encoding="utf-8")
    manifest = {
        "ok": completed.returncode == 0,
        "return_code": completed.returncode,
        "command": command,
        "cwd": str(launch_dir),
        "config_path": str(config_path),
        "log_path": str(log_path),
        "result_dir": str(result_dir),
        "overwrite": overwrite,
    }
    manifest_path = None
    if completed.returncode == 0 and result_dir.exists():
        manifest_path = result_dir / "run_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    results = read_results(result_dir) if completed.returncode == 0 and result_dir.exists() else None
    return {
        **manifest,
        "manifest_path": str(manifest_path) if manifest_path else None,
        "stdout": completed.stdout,
        "results": results,
    }


def install_agent_skill(agent: str = "codex", target: str | Path | None = None) -> dict[str, str]:
    """Install the EasyScan_HEP skill from its standalone repository."""
    agent = agent.lower()
    if target is None:
        if agent != "codex":
            raise ValueError('Only "codex" has a built-in default target. Pass target=... for other agents.')
        target = Path.home() / ".codex" / "skills" / "easyscan-hep"
    target = Path(target).expanduser().resolve()
    if target.exists():
        raise FileExistsError(f"Target already exists: {target}. Update it with git pull or choose a new --target.")
    target.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        ["git", "clone", AGENT_SKILL_REPOSITORY, str(target)],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout.strip() or "Failed to clone EasyScan_HEP agent skill repository.")
    return {"source": AGENT_SKILL_REPOSITORY, "target": str(target), "agent": agent, "action": "cloned"}
