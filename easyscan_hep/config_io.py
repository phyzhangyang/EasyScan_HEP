"""Shared configuration parsing helpers for EasyScan_HEP."""

from __future__ import annotations

import csv
from pathlib import Path


def split_row(row: str) -> list[str]:
    """Split one comma-separated EasyScan row while respecting quotes."""
    try:
        return [item.strip() for item in next(csv.reader([row], skipinitialspace=True))]
    except csv.Error:
        return [item.strip() for item in row.split(",")]


def clean_items(items: list[str]) -> list[str]:
    return [item.strip() for item in items if item.strip()]


def rows_for(config, section: str, option: str) -> list[str]:
    if not config.has_option(section, option):
        return []
    return [line.strip() for line in config.get(section, option).splitlines() if line.strip()]


def resolve_path(value: str, base_dir: str | Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (Path(base_dir).expanduser().resolve() / path).resolve()


def parse_scan_method(value: str, base_dir: str | Path) -> tuple[str, str, Path | None]:
    method = value.strip()
    batch_file = ""
    if "/" in method and method.upper().startswith("ONEPOINT"):
        batch_file = method.split("/", 1)[1]
        method = "ONEPOINTBATCH"
    else:
        method = method.upper()
    batch_path = resolve_path(batch_file, base_dir) if batch_file else None
    return method, batch_file, batch_path
