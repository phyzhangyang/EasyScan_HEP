"""Read EasyScan_HEP result directories into machine-friendly summaries."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from .config_io import split_row


KNOWN_TABLES = (
    "ScanResult.txt",
    "All_ScanResult.txt",
    "EMCEEChain.txt",
    "DynestySamples.txt",
    "Previous_ScanResult.txt",
)


def _to_number(value: str) -> int | float | str | None:
    stripped = value.strip()
    if stripped == "":
        return None
    try:
        number = float(stripped)
    except ValueError:
        return stripped
    if math.isfinite(number) and number.is_integer():
        return int(number)
    return number


def _split_data_row(line: str, expected_columns: int) -> list[str]:
    parts = split_row(line)
    if len(parts) == expected_columns:
        return parts
    whitespace_parts = line.split()
    if len(whitespace_parts) == expected_columns:
        return whitespace_parts
    return parts


def read_table(path: str | Path, max_preview: int = 5) -> dict[str, Any]:
    table_path = Path(path).expanduser().resolve()
    if not table_path.is_file():
        return {"path": str(table_path), "exists": False, "rows": 0, "columns": [], "preview": []}

    lines = [line.strip() for line in table_path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if not lines:
        return {"path": str(table_path), "exists": True, "rows": 0, "columns": [], "preview": []}

    columns = split_row(lines[0])
    preview = []
    row_count = 0
    best_row = None
    best_metric = None
    best_mode = "min"
    if "Chi2" in columns:
        best_metric = "Chi2"
    elif "-2lnlike" in columns:
        best_metric = "-2lnlike"
    elif "log_probability" in columns:
        best_metric = "log_probability"
        best_mode = "max"
    elif "probability" in columns:
        best_metric = "probability"
        best_mode = "max"

    for line in lines[1:]:
        values = _split_data_row(line, len(columns))
        if len(values) != len(columns):
            continue
        row = {name: _to_number(value) for name, value in zip(columns, values)}
        row_count += 1
        if len(preview) < max_preview:
            preview.append(row)
        if best_metric:
            metric_value = row.get(best_metric)
            if isinstance(metric_value, (int, float)) and math.isfinite(metric_value):
                if best_row is None:
                    best_row = row
                else:
                    current = best_row.get(best_metric)
                    if not isinstance(current, (int, float)):
                        best_row = row
                    elif best_mode == "min" and metric_value < current:
                        best_row = row
                    elif best_mode == "max" and metric_value > current:
                        best_row = row

    return {
        "path": str(table_path),
        "exists": True,
        "rows": row_count,
        "columns": columns,
        "preview": preview,
        "best_metric": best_metric,
        "best_row": best_row,
    }


def read_results(result_dir: str | Path, max_preview: int = 5) -> dict[str, Any]:
    directory = Path(result_dir).expanduser().resolve()
    files: list[dict[str, Any]] = []
    plots: list[dict[str, Any]] = []
    tables: dict[str, dict[str, Any]] = {}

    if directory.exists():
        for path in sorted(directory.rglob("*")):
            if not path.is_file():
                continue
            item = {"name": path.name, "path": str(path), "size": path.stat().st_size}
            if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}:
                plots.append(item)
            else:
                files.append(item)

    for table_name in KNOWN_TABLES:
        table_path = directory / table_name
        if table_path.exists():
            tables[table_name] = read_table(table_path, max_preview=max_preview)

    multinest_table = directory / "MultiNestData" / ".txt"
    if multinest_table.exists():
        tables["MultiNestData/.txt"] = read_table(multinest_table, max_preview=max_preview)

    primary = tables.get("ScanResult.txt")
    return {
        "result_dir": str(directory),
        "exists": directory.exists(),
        "files": files,
        "plots": plots,
        "tables": tables,
        "row_count": primary.get("rows", 0) if primary else 0,
        "columns": primary.get("columns", []) if primary else [],
        "best_row": primary.get("best_row") if primary else None,
    }


def format_results_summary(summary: dict[str, Any]) -> str:
    lines = [f"Result directory: {summary['result_dir']}"]
    if not summary.get("exists"):
        lines.append("Result directory does not exist.")
        return "\n".join(lines)
    lines.append(f"Rows in ScanResult.txt: {summary.get('row_count', 0)}")
    if summary.get("columns"):
        lines.append("Columns: " + ", ".join(summary["columns"]))
    if summary.get("best_row"):
        lines.append("Best row:")
        for key, value in summary["best_row"].items():
            lines.append(f"  {key}: {value}")
    if summary.get("plots"):
        lines.append(f"Plots: {len(summary['plots'])}")
    if summary.get("files"):
        lines.append(f"Files: {len(summary['files'])}")
    return "\n".join(lines)
