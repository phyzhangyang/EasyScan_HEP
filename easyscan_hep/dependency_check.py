"""Runtime dependency checks for EasyScan_HEP command-line entry points."""

from __future__ import annotations

import configparser
import importlib
from importlib import metadata
import os
import re
import shlex
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


MIN_PYTHON = (3, 9)


@dataclass(frozen=True)
class RuntimeRequirement:
    pip_name: str
    import_name: str
    min_version: str | None = None
    purpose: str = ""

    @property
    def pip_spec(self) -> str:
        if self.min_version:
            return f"{self.pip_name}>={self.min_version}"
        return self.pip_name


@dataclass(frozen=True)
class RequirementIssue:
    requirement: RuntimeRequirement
    message: str


CORE_REQUIREMENTS = (
    RuntimeRequirement("numpy", "numpy", "1.20", "core scan arrays and file bounds"),
    RuntimeRequirement("scipy", "scipy", "1.7", "BESTFIT and contour interpolation"),
    RuntimeRequirement("matplotlib", "matplotlib", "3.5", "plot generation"),
    RuntimeRequirement("pandas", "pandas", "1.3", "scan-result and plot data reading"),
)

UI_REQUIREMENTS = (
    RuntimeRequirement("fastapi", "fastapi", "0.95", "local Web UI backend"),
    RuntimeRequirement("uvicorn", "uvicorn", "0.20", "local Web UI server"),
    RuntimeRequirement("jinja2", "jinja2", "3.0", "local Web UI templates"),
    RuntimeRequirement("python-multipart", "multipart", "0.0.6", "form parsing support"),
)

SAMPLER_REQUIREMENTS = {
    "EMCEE": (RuntimeRequirement("emcee", "emcee", "3.1", "EMCEE ensemble sampler"),),
    "DYNESTY": (RuntimeRequirement("dynesty", "dynesty", "2.1", "DYNESTY nested sampler"),),
    "MULTINEST": (RuntimeRequirement("pymultinest", "pymultinest", "2.12", "MultiNest sampler"),),
}


def _version_key(value: str) -> tuple[int, ...]:
    parts = re.findall(r"\d+", value)
    return tuple(int(part) for part in parts[:3])


def _version_lt(current: str, required: str) -> bool:
    current_key = _version_key(current)
    required_key = _version_key(required)
    length = max(len(current_key), len(required_key))
    current_key = current_key + (0,) * (length - len(current_key))
    required_key = required_key + (0,) * (length - len(required_key))
    return current_key < required_key


def _distribution_version(requirement: RuntimeRequirement, module: object) -> str | None:
    try:
        return metadata.version(requirement.pip_name)
    except metadata.PackageNotFoundError:
        return getattr(module, "__version__", None)


def check_requirement(requirement: RuntimeRequirement) -> RequirementIssue | None:
    if requirement.import_name == "matplotlib":
        os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "easyscan-matplotlib"))
    try:
        module = importlib.import_module(requirement.import_name)
    except Exception as exc:
        first_line = str(exc).strip().splitlines()[0] if str(exc).strip() else type(exc).__name__
        return RequirementIssue(requirement, f"can not import {requirement.import_name}: {first_line}")

    if requirement.min_version:
        version = _distribution_version(requirement, module)
        if version is None:
            return RequirementIssue(requirement, f"can not determine installed version for {requirement.pip_name}")
        if _version_lt(version, requirement.min_version):
            return RequirementIssue(requirement, f"installed version is {version}, but >= {requirement.min_version} is required")
    return None


def scan_method_from_config(config_path: str | Path | None) -> str | None:
    if not config_path:
        return None
    path = Path(config_path).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    if not path.is_file():
        return None

    config = configparser.ConfigParser(interpolation=None)
    config.optionxform = str
    try:
        config.read(path, encoding="utf-8")
    except configparser.Error:
        return None
    if not config.has_option("scan", "Scan method"):
        return None
    method = config.get("scan", "Scan method").strip()
    return method.split("/", 1)[0].upper() if method else None


def requirements_for(mode: str, config_path: str | Path | None = None) -> list[RuntimeRequirement]:
    mode = mode.lower()
    requirements: list[RuntimeRequirement] = []
    if mode in {"scan", "run", "ui"}:
        requirements.extend(CORE_REQUIREMENTS)
    if mode == "ui":
        requirements.extend(UI_REQUIREMENTS)

    method = scan_method_from_config(config_path)
    if method:
        requirements.extend(SAMPLER_REQUIREMENTS.get(method, ()))

    unique: dict[str, RuntimeRequirement] = {}
    for requirement in requirements:
        unique[requirement.pip_name] = requirement
    return list(unique.values())


def install_command(requirements: Iterable[RuntimeRequirement]) -> str:
    specs = [shlex.quote(requirement.pip_spec) for requirement in requirements]
    return f"{shlex.quote(sys.executable)} -m pip install " + " ".join(specs)


def check_runtime_requirements(mode: str, config_path: str | Path | None = None) -> list[str]:
    messages: list[str] = []
    if sys.version_info < MIN_PYTHON:
        messages.append(
            "Python %d.%d or newer is required; current Python is %s at %s."
            % (MIN_PYTHON[0], MIN_PYTHON[1], sys.version.split()[0], sys.executable)
        )
        messages.append("Install Python 3.9 or newer, then reinstall EasyScan_HEP from GitHub with:")
        messages.append("  python3 -m pip install -U git+https://github.com/phyzhangyang/EasyScan_HEP.git")

    issues = [issue for issue in (check_requirement(requirement) for requirement in requirements_for(mode, config_path)) if issue]
    if issues:
        messages.append("Missing or unusable Python packages:")
        for issue in issues:
            purpose = f" ({issue.requirement.purpose})" if issue.requirement.purpose else ""
            messages.append(f"  - {issue.requirement.pip_spec}{purpose}: {issue.message}")
        messages.append("Install or update them with:")
        messages.append(f"  {install_command(issue.requirement for issue in issues)}")
        messages.append("Use the same Python shown in the command above to avoid pip/python mismatches.")

    return messages


def ensure_runtime_requirements(mode: str, config_path: str | Path | None = None) -> None:
    messages = check_runtime_requirements(mode, config_path)
    if not messages:
        return
    print("EasyScan_HEP environment check failed.", file=sys.stderr)
    print("", file=sys.stderr)
    print("\n".join(messages), file=sys.stderr)
    sys.exit(1)
