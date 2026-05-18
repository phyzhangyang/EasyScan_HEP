from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import threading
import time
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field


REPO_ROOT = Path(__file__).resolve().parents[1]
UI_ROOT = REPO_ROOT / "ui"
RUNS_ROOT = REPO_ROOT / "ui_runs"
PYTHON_EXE = REPO_ROOT / ".venv" / "bin" / "python"
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*m")

app = FastAPI(title="EasyScan_HEP Local UI")
app.mount("/static", StaticFiles(directory=UI_ROOT / "static"), name="static")
templates = Jinja2Templates(directory=UI_ROOT / "templates")


class InputParameter(BaseModel):
    name: str = "x"
    prior: str = "flat"
    value: str = ""
    minimum: str = "0"
    maximum: str = "3.14"
    bins: str = "10"
    interval: str = "10"
    initial: str = "1.5"


class FileEntry(BaseModel):
    file_id: str = "1"
    path: str = ""


class VariableEntry(BaseModel):
    name: str = ""
    file_id: str = "1"
    method: str = "Position"
    arguments: str = ""


class ProgramEntry(BaseModel):
    program_name: str = "TestFunction"
    execute_command: str = "./TestFunction.py"
    command_path: str = "utils/"
    input_files: list[FileEntry] = Field(default_factory=lambda: [FileEntry(path="utils/TestFunction_input.dat")])
    input_variables: list[VariableEntry] = Field(
        default_factory=lambda: [
            VariableEntry(name="x", arguments="1, 1"),
            VariableEntry(name="y", arguments="1, 2"),
        ]
    )
    output_files: list[FileEntry] = Field(default_factory=lambda: [FileEntry(path="utils/TestFunction_output.dat")])
    output_variables: list[VariableEntry] = Field(default_factory=lambda: [VariableEntry(name="f", arguments="1, 2")])
    command_executor: str = ""
    time_limit: str = ""
    clean_output: str = ""
    bounds: str = ""


class GaussianEntry(BaseModel):
    variable: str = "f"
    mean: str = "1"
    deviation: str = "0.2"
    kind: str = "symm"
    label: str = ""


class FreeFormChi2Entry(BaseModel):
    variable: str = ""
    label: str = ""


class PlotEntry(BaseModel):
    kind: str = "Color"
    x: str = "x"
    y: str = "y"
    value: str = "f"
    figure_name: str = "test_run_random"


class EasyScanConfig(BaseModel):
    result_folder: str = "ui_runs/example/result"
    scan_method: str = "RANDOM"
    batch_file: str = "utils/OnePointBatch.in"
    number_of_points: str = "100"
    random_seed: str = ""
    print_interval: str = "1"
    parallel_threads: str = "1"
    parallel_folder: str = "utils"
    input_parameters: list[InputParameter] = Field(
        default_factory=lambda: [
            InputParameter(name="x", minimum="0", maximum="3.14"),
            InputParameter(name="y", minimum="-3.14", maximum="3.14"),
        ]
    )
    programs: list[ProgramEntry] = Field(default_factory=lambda: [ProgramEntry()])
    gaussian_constraints: list[GaussianEntry] = Field(default_factory=lambda: [GaussianEntry()])
    freeform_chi2: list[FreeFormChi2Entry] = Field(default_factory=list)
    plots: list[PlotEntry] = Field(
        default_factory=lambda: [
            PlotEntry(kind="Color", x="x", y="y", value="f", figure_name="test_run_random"),
            PlotEntry(kind="Color", x="x", y="y", value="Chi2", figure_name="Color_x_y_Chi2"),
        ]
    )


class RunRecord:
    def __init__(self, run_id: str, run_dir: Path, config_path: Path):
        self.run_id = run_id
        self.run_dir = run_dir
        self.config_path = config_path
        self.log_path = run_dir / "run.log"
        self.status = "queued"
        self.return_code: int | None = None
        self.process: subprocess.Popen[str] | None = None
        self.logs: list[str] = []
        self.started_at = time.time()
        self.ended_at: float | None = None
        self.lock = threading.Lock()

    def append_log(self, line: str) -> None:
        with self.lock:
            self.logs.append(line)
            if len(self.logs) > 2000:
                self.logs = self.logs[-2000:]


runs: dict[str, RunRecord] = {}


def as_payload(model: BaseModel) -> dict[str, Any]:
    if hasattr(model, "model_dump"):
        return model.model_dump()
    return model.dict()


def path_for_config(path_value: str) -> str:
    if not path_value:
        return ""
    path = Path(path_value).expanduser()
    if path.is_absolute():
        try:
            return path.resolve().relative_to(REPO_ROOT).as_posix()
        except ValueError:
            return str(path)
    return path_value


def clean_items(items: list[str]) -> list[str]:
    return [item.strip() for item in items if item.strip()]


def clean_log_line(line: str) -> str:
    return ANSI_ESCAPE_RE.sub("", line)


def split_row(row: str) -> list[str]:
    return clean_items(row.split(","))


def variable_line(item: VariableEntry) -> str:
    parts = [item.name, item.file_id, item.method]
    parts.extend(clean_items(item.arguments.split(",")))
    return ",       ".join(parts)


def input_parameter_line(item: InputParameter, scan_method: str) -> str:
    method = scan_method.upper()
    if item.prior.upper() == "FIXED" or method in {"ONEPOINT", "ONEPOINTBATCH"}:
        value = item.value or item.initial or item.minimum
        return f"{item.name},     Fixed,  {value}"
    if method == "GRID":
        return f"{item.name},     {item.prior},   {item.minimum},       {item.maximum},  {item.bins or '10'}"
    if method == "MCMC":
        return (
            f"{item.name},     {item.prior},   {item.minimum},       {item.maximum},  "
            f"{item.interval or '10'},       {item.initial or item.minimum}"
        )
    return f"{item.name},     {item.prior},   {item.minimum},       {item.maximum}"


def multiline_option(label: str, rows: list[str], indent: int = 19) -> str:
    if not rows:
        return ""
    spaces = " " * indent
    lines = [f"{label}:  {rows[0]}"]
    lines.extend(f"{spaces}{row}" for row in rows[1:])
    return "\n".join(lines)


def serialize_config(config: EasyScanConfig) -> str:
    method = config.scan_method.upper()
    scan_method = method
    if method == "ONEPOINTBATCH":
        scan_method = f"ONEPOINT/{path_for_config(config.batch_file)}"

    lines: list[str] = ["[scan]"]
    lines.append(f"Result folder name:  {path_for_config(config.result_folder)}")
    lines.append(f"Scan method:       {scan_method}")
    lines.append(multiline_option("Input parameters", [input_parameter_line(item, method) for item in config.input_parameters]))
    if method not in {"ONEPOINT", "ONEPOINTBATCH", "GRID"} and config.number_of_points:
        lines.append(f"Number of points:  {config.number_of_points}")
    elif method == "MCMC" and config.number_of_points:
        lines.append(f"Number of points:  {config.number_of_points}")
    if config.random_seed:
        lines.append(f"Random seed:       {config.random_seed}")
    if config.print_interval:
        lines.append(f"Interval of print: {config.print_interval}")
    if config.parallel_threads:
        lines.append(f"Parallel threads:  {config.parallel_threads}")
    if config.parallel_folder:
        lines.append(f"Parallel folder:   {path_for_config(config.parallel_folder)}")

    for index, program in enumerate(config.programs, start=1):
        lines.extend(["", f"[program{index}]"])
        if program.program_name:
            lines.append(f"Program name:    {program.program_name}")
        lines.append(f"Execute command: {program.execute_command}")
        lines.append(f"Command path:    {path_for_config(program.command_path)}")
        if program.input_files:
            rows = [f"{item.file_id}, {path_for_config(item.path)}" for item in program.input_files if item.file_id and item.path]
            lines.append(multiline_option("Input file", rows, indent=18))
        if program.input_variables:
            rows = [variable_line(item) for item in program.input_variables if item.name]
            lines.append(multiline_option("Input variable", rows, indent=17))
        if program.output_files:
            rows = [f"{item.file_id}, {path_for_config(item.path)}" for item in program.output_files if item.file_id and item.path]
            lines.append(multiline_option("Output file", rows, indent=18))
        if program.output_variables:
            rows = [variable_line(item) for item in program.output_variables if item.name]
            lines.append(multiline_option("Output variable", rows, indent=17))
        if program.command_executor:
            lines.append(f"Command executor: {program.command_executor}")
        if program.time_limit:
            lines.append(f"Time limit in minute: {program.time_limit}")
        if program.clean_output:
            lines.append(f"Clean output file: {program.clean_output}")
        if program.bounds.strip():
            lines.append("Bound:           " + program.bounds.strip().replace("\n", "\n                 "))

    constraint_lines: list[str] = []
    gaussian_rows = []
    for item in config.gaussian_constraints:
        if not item.variable:
            continue
        row = [item.variable, item.mean, item.deviation]
        if item.kind and item.kind != "symm":
            row.append(item.kind)
        if item.label:
            if len(row) == 3:
                row.append(item.kind or "symm")
            row.append(item.label)
        gaussian_rows.append(",      ".join(row))
    if gaussian_rows:
        constraint_lines.append(multiline_option("Gaussian", gaussian_rows, indent=12))
    freeform_rows = []
    for item in config.freeform_chi2:
        if item.variable:
            freeform_rows.append(",      ".join(clean_items([item.variable, item.label])))
    if freeform_rows:
        constraint_lines.append(multiline_option("FreeFormChi2", freeform_rows, indent=15))
    if constraint_lines:
        lines.extend(["", "[constraint]", *constraint_lines])

    plot_groups: dict[str, list[str]] = {"Histogram": [], "Scatter": [], "Color": [], "Contour": []}
    for item in config.plots:
        kind = item.kind.capitalize()
        if kind not in plot_groups:
            continue
        if kind == "Histogram":
            row = ",       ".join(clean_items([item.x, item.figure_name]))
        elif kind == "Scatter":
            row = ",       ".join(clean_items([item.x, item.y, item.figure_name]))
        else:
            row = ",       ".join(clean_items([item.x, item.y, item.value, item.figure_name]))
        if row:
            plot_groups[kind].append(row)
    plot_lines = [multiline_option(kind, rows, indent=len(kind) + 3) for kind, rows in plot_groups.items() if rows]
    if plot_lines:
        lines.extend(["", "[plot]", *plot_lines])

    return "\n".join(line for line in lines if line != "") + "\n"


def parse_easy_scan_template(path: Path) -> EasyScanConfig:
    sections: dict[str, dict[str, list[str]]] = {}
    current_section = ""
    current_key = ""
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        stripped = raw_line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current_section = stripped.strip("[]").lower()
            sections.setdefault(current_section, {})
            current_key = ""
            continue
        if not current_section:
            continue
        if not raw_line.startswith((" ", "\t")) and ":" in raw_line:
            key, value = raw_line.split(":", 1)
            current_key = key.strip()
            sections[current_section].setdefault(current_key, [])
            if value.strip():
                sections[current_section][current_key].append(value.strip())
        elif current_key and stripped:
            sections[current_section][current_key].append(stripped)

    config = EasyScanConfig()
    scan = sections.get("scan", {})
    config.result_folder = (scan.get("Result folder name") or [config.result_folder])[0]
    method = (scan.get("Scan method") or [config.scan_method])[0]
    if "/" in method and method.upper().startswith("ONEPOINT"):
        config.scan_method, config.batch_file = "ONEPOINTBATCH", method.split("/", 1)[1]
    else:
        config.scan_method = method.upper()
    config.number_of_points = (scan.get("Number of points") or [config.number_of_points])[0]
    config.random_seed = (scan.get("Random seed") or [config.random_seed])[0]
    config.print_interval = (scan.get("Interval of print") or [config.print_interval])[0]
    config.parallel_threads = (scan.get("Parallel threads") or [config.parallel_threads])[0]
    config.parallel_folder = (scan.get("Parallel folder") or [config.parallel_folder])[0]
    input_parameters = []
    for row in scan.get("Input parameters", []):
        parts = split_row(row)
        if len(parts) < 3:
            continue
        parameter = InputParameter(name=parts[0], prior=parts[1])
        if parts[1].upper() == "FIXED":
            parameter.value = parts[2]
            parameter.initial = parts[2]
        else:
            parameter.minimum = parts[2]
            parameter.maximum = parts[3] if len(parts) > 3 else ""
            if len(parts) > 4 and config.scan_method == "GRID":
                parameter.bins = parts[4]
            elif len(parts) > 4:
                parameter.interval = parts[4]
            if len(parts) > 5:
                parameter.initial = parts[5]
        input_parameters.append(parameter)
    if input_parameters:
        config.input_parameters = input_parameters

    programs = []
    for section_name in sorted(name for name in sections if name.startswith("program")):
        source = sections[section_name]
        program = ProgramEntry()
        program.program_name = (source.get("Program name") or [program.program_name])[0]
        program.execute_command = (source.get("Execute command") or [program.execute_command])[0]
        program.command_path = (source.get("Command path") or [program.command_path])[0]
        program.command_executor = (source.get("Command executor") or [program.command_executor])[0]
        program.time_limit = (source.get("Time limit in minute") or [program.time_limit])[0]
        program.clean_output = (source.get("Clean output file") or [program.clean_output])[0]
        program.bounds = "\n".join(source.get("Bound", []))
        program.input_files = [FileEntry(file_id=parts[0], path=parts[1]) for parts in map(split_row, source.get("Input file", [])) if len(parts) >= 2]
        program.output_files = [FileEntry(file_id=parts[0], path=parts[1]) for parts in map(split_row, source.get("Output file", [])) if len(parts) >= 2]
        program.input_variables = [
            VariableEntry(name=parts[0], file_id=parts[1], method=parts[2], arguments=", ".join(parts[3:]))
            for parts in map(split_row, source.get("Input variable", []))
            if len(parts) >= 3
        ]
        program.output_variables = [
            VariableEntry(name=parts[0], file_id=parts[1], method=parts[2], arguments=", ".join(parts[3:]))
            for parts in map(split_row, source.get("Output variable", []))
            if len(parts) >= 3
        ]
        programs.append(program)
    if programs:
        config.programs = programs

    constraint = sections.get("constraint", {})
    gaussian_constraints = []
    for row in constraint.get("Gaussian", []):
        parts = split_row(row)
        if len(parts) >= 3:
            gaussian_constraints.append(
                GaussianEntry(
                    variable=parts[0],
                    mean=parts[1],
                    deviation=parts[2],
                    kind=parts[3] if len(parts) > 3 else "symm",
                    label=parts[4] if len(parts) > 4 else "",
                )
            )
    if gaussian_constraints:
        config.gaussian_constraints = gaussian_constraints

    plot = sections.get("plot", {})
    plots = []
    for kind in ("Histogram", "Scatter", "Color", "Contour"):
        for row in plot.get(kind, []):
            parts = split_row(row)
            if kind == "Histogram" and parts:
                plots.append(PlotEntry(kind=kind, x=parts[0], figure_name=parts[1] if len(parts) > 1 else ""))
            elif kind == "Scatter" and len(parts) >= 2:
                plots.append(PlotEntry(kind=kind, x=parts[0], y=parts[1], figure_name=parts[2] if len(parts) > 2 else ""))
            elif len(parts) >= 3:
                plots.append(PlotEntry(kind=kind, x=parts[0], y=parts[1], value=parts[2], figure_name=parts[3] if len(parts) > 3 else ""))
    if plots:
        config.plots = plots
    return config


def default_config() -> EasyScanConfig:
    template_path = REPO_ROOT / "templates" / "example_random.ini"
    if template_path.exists():
        try:
            return parse_easy_scan_template(template_path)
        except Exception:
            return EasyScanConfig()
    return EasyScanConfig()



def allowed_roots() -> list[Path]:
    roots = [REPO_ROOT, Path.home()]
    unique: list[Path] = []
    for root in roots:
        resolved = root.resolve()
        if resolved not in unique:
            unique.append(resolved)
    return unique


def ensure_allowed(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    for root in allowed_roots():
        if resolved == root or root in resolved.parents:
            return resolved
    raise HTTPException(status_code=403, detail="Path is outside the allowed roots.")


def result_files(record: RunRecord) -> dict[str, Any]:
    result_dir = record.run_dir / "result"
    files = []
    plots = []
    if result_dir.exists():
        for path in sorted(result_dir.rglob("*")):
            if not path.is_file():
                continue
            item = {"name": path.name, "path": str(path), "url": f"/api/files?path={path}"}
            if path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg"}:
                plots.append(item)
            else:
                files.append(item)
    return {"result_dir": str(result_dir), "files": files, "plots": plots}


def run_worker(record: RunRecord) -> None:
    record.status = "running"
    command = [str(PYTHON_EXE), "bin/easyscan.py", str(record.config_path)]
    with record.log_path.open("w", encoding="utf-8") as log_file:
        log_file.write("$ " + " ".join(command) + "\n")
        log_file.flush()
        try:
            record.process = subprocess.Popen(
                command,
                cwd=REPO_ROOT,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                preexec_fn=os.setsid,
            )
            if record.process.stdin:
                record.process.stdin.write("r\n")
                record.process.stdin.flush()
            assert record.process.stdout is not None
            for line in record.process.stdout:
                line = clean_log_line(line)
                log_file.write(line)
                log_file.flush()
                record.append_log(line)
            record.return_code = record.process.wait()
            if record.status != "stopped":
                record.status = "completed" if record.return_code == 0 else "failed"
        except Exception as exc:  # pragma: no cover - surfaced through API
            record.status = "failed"
            record.append_log(f"Backend error: {exc}\n")
            log_file.write(f"Backend error: {exc}\n")
        finally:
            record.ended_at = time.time()


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> HTMLResponse:
    payload = as_payload(default_config())
    static_version = int(max((UI_ROOT / "static" / "app.js").stat().st_mtime, (UI_ROOT / "static" / "styles.css").stat().st_mtime))
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "default_config": json.dumps(payload),
            "static_version": static_version,
        },
    )


@app.post("/api/config")
def preview_config(config: EasyScanConfig) -> dict[str, str]:
    return {"ini": serialize_config(config)}


@app.post("/api/runs")
def start_run(config: EasyScanConfig) -> dict[str, Any]:
    RUNS_ROOT.mkdir(exist_ok=True)
    run_id = time.strftime("%Y%m%d_%H%M%S")
    suffix = 1
    while (RUNS_ROOT / run_id).exists():
        suffix += 1
        run_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{suffix}"
    run_dir = RUNS_ROOT / run_id
    run_dir.mkdir(parents=True)
    result_dir = run_dir / "result"
    if result_dir.exists():
        shutil.rmtree(result_dir)
    config.result_folder = str(result_dir.relative_to(REPO_ROOT))
    config_path = run_dir / "config.ini"
    config_path.write_text(serialize_config(config), encoding="utf-8")
    record = RunRecord(run_id, run_dir, config_path)
    runs[run_id] = record
    thread = threading.Thread(target=run_worker, args=(record,), daemon=True)
    thread.start()
    return {"run_id": run_id, "status": record.status, "config_path": str(config_path)}


@app.get("/api/runs/{run_id}")
def get_run(run_id: str) -> dict[str, Any]:
    record = runs.get(run_id)
    if record is None:
        run_dir = RUNS_ROOT / run_id
        config_path = run_dir / "config.ini"
        if not config_path.exists():
            raise HTTPException(status_code=404, detail="Run not found.")
        record = RunRecord(run_id, run_dir, config_path)
        record.status = "completed"
        record.log_path = run_dir / "run.log"
    with record.lock:
        log_text = "".join(record.logs)
    if not log_text and record.log_path.exists():
        log_text = record.log_path.read_text(encoding="utf-8", errors="replace")
    return {
        "run_id": run_id,
        "status": record.status,
        "return_code": record.return_code,
        "config_path": str(record.config_path),
        "log_path": str(record.log_path),
        "logs": log_text,
        "results": result_files(record),
    }


@app.post("/api/runs/{run_id}/stop")
def stop_run(run_id: str) -> dict[str, str]:
    record = runs.get(run_id)
    if record is None or record.process is None:
        raise HTTPException(status_code=404, detail="Active run not found.")
    if record.process.poll() is None:
        record.status = "stopped"
        os.killpg(os.getpgid(record.process.pid), signal.SIGTERM)
    return {"status": record.status}


@app.get("/api/browse")
def browse(
    path: str | None = Query(default=None),
    select: str = Query(default="file"),
) -> dict[str, Any]:
    roots = allowed_roots()
    current = ensure_allowed(Path(path)) if path else roots[0]
    if current.is_file():
        current = current.parent
    entries = []
    try:
        children = sorted(current.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
    except OSError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    for child in children:
        if child.name.startswith("."):
            continue
        entries.append(
            {
                "name": child.name,
                "path": str(child),
                "is_dir": child.is_dir(),
                "selectable": child.is_dir() if select == "dir" else child.is_file(),
            }
        )
    parent = current.parent if current != current.anchor else current
    return {
        "current": str(current),
        "parent": str(parent),
        "roots": [{"name": str(root), "path": str(root)} for root in roots],
        "entries": entries,
    }


@app.get("/api/files")
def get_file(path: str) -> FileResponse:
    resolved = ensure_allowed(Path(path))
    if not resolved.is_file():
        raise HTTPException(status_code=404, detail="File not found.")
    return FileResponse(resolved)
