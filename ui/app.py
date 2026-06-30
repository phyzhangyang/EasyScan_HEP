from __future__ import annotations

import json
import os
import re
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import FileResponse, HTMLResponse, Response
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field


REPO_ROOT = Path(os.environ.get("EASYSCAN_ROOT", Path(__file__).resolve().parents[1])).expanduser().resolve()
RUN_CWD = Path(os.environ.get("EASYSCAN_UI_CWD", Path.cwd())).expanduser().resolve()
UI_ROOT = Path(__file__).resolve().parent
SRC_ROOT = REPO_ROOT / "src"
PYTHON_EXE = Path(sys.executable)
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*m")

if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from easyscan_hep.config_io import split_row as split_config_row
from easyscan_hep.results import read_results
from config_checker import check_config_text, format_check_report

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
    config_file_path: str = ""
    result_folder: str = "result"
    scan_method: str = "RANDOM"
    batch_file: str = "utils/OnePointBatch.in"
    number_of_points: str = "100"
    mcmc_walkers: str = ""
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


class ImportConfigRequest(BaseModel):
    content: str


class LLMConfigRequest(BaseModel):
    provider: str = "openai"
    model: str = ""
    api_key: str = ""
    base_url: str = ""
    prompt: str
    current_config: EasyScanConfig


class RunRecord:
    def __init__(self, run_id: str, run_dir: Path, config_path: Path, log_path: Path, result_dir: Path):
        self.run_id = run_id
        self.run_dir = run_dir
        self.config_path = config_path
        self.log_path = log_path
        self.result_dir = result_dir
        self.status = "queued"
        self.return_code: Optional[int] = None
        self.process: Optional[subprocess.Popen[str]] = None
        self.logs: list[str] = []
        self.started_at = time.time()
        self.ended_at: Optional[float] = None
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
            return path.resolve().relative_to(RUN_CWD).as_posix()
        except ValueError:
            pass
        try:
            return path.resolve().relative_to(REPO_ROOT).as_posix()
        except ValueError:
            return str(path)
    return path_value


def resolve_config_path(path_value: str) -> Path:
    path = Path(path_value).expanduser()
    if path.is_absolute():
        return path.resolve()
    return (RUN_CWD / path).resolve()


def config_base_dir(config: EasyScanConfig) -> Path:
    if config.config_file_path.strip():
        path = Path(config.config_file_path).expanduser()
        if not path.is_absolute():
            path = RUN_CWD / path
        return path.resolve().parent
    return RUN_CWD


def run_artifact_paths(run_id: str) -> tuple[Path, Path]:
    return RUN_CWD / f"easyscan_{run_id}.ini", RUN_CWD / f"easyscan_{run_id}.log"


def config_target_paths(config: EasyScanConfig, run_id: str) -> tuple[Path, Path]:
    if config.config_file_path.strip():
        config_path = Path(config.config_file_path).expanduser()
        if not config_path.is_absolute():
            config_path = RUN_CWD / config_path
        config_path = config_path.resolve()
        return config_path, config_path.with_suffix(".log")
    return run_artifact_paths(run_id)


def clean_items(items: list[str]) -> list[str]:
    return [item.strip() for item in items if item.strip()]


def clean_log_line(line: str) -> str:
    return ANSI_ESCAPE_RE.sub("", line)


def split_row(row: str) -> list[str]:
    return clean_items(split_config_row(row))


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
    if method in {"MCMC", "EMCEE"}:
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
    if method == "EMCEE" and config.mcmc_walkers:
        lines.append(f"MCMC walkers:      {config.mcmc_walkers}")
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


def parse_easy_scan_text(text: str) -> EasyScanConfig:
    sections: dict[str, dict[str, list[str]]] = {}
    current_section = ""
    current_key = ""
    for raw_line in text.splitlines():
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
    config.mcmc_walkers = (scan.get("MCMC walkers") or [config.mcmc_walkers])[0]
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
    config.gaussian_constraints = gaussian_constraints
    freeform_chi2 = []
    for row in constraint.get("FreeFormChi2", []):
        parts = split_row(row)
        if parts:
            freeform_chi2.append(FreeFormChi2Entry(variable=parts[0], label=parts[1] if len(parts) > 1 else ""))
    config.freeform_chi2 = freeform_chi2

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
    config.plots = plots
    return config


def parse_easy_scan_template(path: Path) -> EasyScanConfig:
    return parse_easy_scan_text(path.read_text(encoding="utf-8"))


def load_config_file(path: Path) -> EasyScanConfig:
    config = parse_easy_scan_template(path)
    config.config_file_path = str(path)
    return config


LLM_PROVIDERS: dict[str, dict[str, str]] = {
    "openai": {
        "base_url": "https://api.openai.com/v1/chat/completions",
        "model": "gpt-4.1-mini",
        "requires_key": "1",
    },
    "deepseek": {
        "base_url": "https://api.deepseek.com/chat/completions",
        "model": "deepseek-chat",
        "requires_key": "1",
    },
    "qwen": {
        "base_url": "https://dashscope.aliyuncs.com/compatible-mode/v1/chat/completions",
        "model": "qwen-plus",
        "requires_key": "1",
    },
    "moonshot": {
        "base_url": "https://api.moonshot.cn/v1/chat/completions",
        "model": "moonshot-v1-8k",
        "requires_key": "1",
    },
    "zhipu": {
        "base_url": "https://open.bigmodel.cn/api/paas/v4/chat/completions",
        "model": "glm-4-flash",
        "requires_key": "1",
    },
    "siliconflow": {
        "base_url": "https://api.siliconflow.cn/v1/chat/completions",
        "model": "Qwen/Qwen2.5-72B-Instruct",
        "requires_key": "1",
    },
    "openrouter": {
        "base_url": "https://openrouter.ai/api/v1/chat/completions",
        "model": "openai/gpt-4o-mini",
        "requires_key": "1",
    },
    "ollama": {
        "base_url": "http://127.0.0.1:11434/v1/chat/completions",
        "model": "llama3.1",
        "requires_key": "0",
    },
    "custom": {
        "base_url": "",
        "model": "",
        "requires_key": "0",
    },
}


def llm_provider_defaults(provider: str) -> dict[str, str]:
    return LLM_PROVIDERS.get(provider.lower().strip(), LLM_PROVIDERS["custom"])


def normalize_chat_url(provider: str, base_url: str) -> str:
    value = base_url.strip()
    if not value:
        value = llm_provider_defaults(provider).get("base_url", "")
    if not value:
        raise HTTPException(status_code=400, detail="LLM API base URL is required.")
    lowered = value.rstrip("/").lower()
    if lowered.endswith("/chat/completions"):
        return value.rstrip("/")
    if lowered.endswith("/v1") or lowered.endswith("/compatible-mode/v1"):
        return value.rstrip("/") + "/chat/completions"
    if provider.lower().strip() == "deepseek":
        return value.rstrip("/") + "/chat/completions"
    return value.rstrip("/") + "/v1/chat/completions"


def llm_system_prompt() -> str:
    return """You convert natural-language EasyScan_HEP requests into a complete EasyScan_HEP INI file.
Return only INI text. Do not use Markdown fences or explanations.
Start from the current INI and change only what the user requests.
Keep existing local paths, program blocks, variable names, constraints, and plots unless the user asks to change them.

EasyScan_HEP syntax summary:
- Required [scan] options: Result folder name, Scan method, Input parameters.
- Supported Scan method values: RANDOM, GRID, BESTFIT, MCMC, EMCEE, DYNESTY, MULTINEST, ONEPOINT, ONEPOINT/path/to/batch.in.
- Input parameters rows:
  - RANDOM, BESTFIT, DYNESTY, MULTINEST: name, prior, min, max.
  - GRID: name, prior, min, max, bins. To get N grid points on an axis, use bins = N - 1.
  - MCMC, EMCEE: name, prior, min, max, interval, initial.
  - Fixed or one-point mode: name, Fixed, value.
  - prior must be flat, log, or Fixed.
- EMCEE may use optional "MCMC walkers"; it must be at least twice the number of sampled input parameters.
- Likelihood-based methods BESTFIT, MCMC, EMCEE, DYNESTY and MULTINEST need [constraint] with Gaussian or FreeFormChi2.
- [programN] blocks describe external programs with Execute command, Command path, Input/Output file and variable rows.
- [plot] rows: Histogram x, figure; Scatter x, y, figure; Color/Contour x, y, value, figure.
"""


def build_llm_messages(request: LLMConfigRequest) -> list[dict[str, str]]:
    current_ini = serialize_config(request.current_config)
    user_prompt = f"""Current EasyScan_HEP INI:

{current_ini}

User request:
{request.prompt.strip()}

Write the full updated EasyScan_HEP INI now."""
    return [
        {"role": "system", "content": llm_system_prompt()},
        {"role": "user", "content": user_prompt},
    ]


def call_openai_compatible_llm(request: LLMConfigRequest) -> str:
    provider = request.provider.lower().strip()
    defaults = llm_provider_defaults(provider)
    model = request.model.strip() or defaults.get("model", "")
    if not model:
        raise HTTPException(status_code=400, detail="LLM model is required.")
    if defaults.get("requires_key") == "1" and not request.api_key.strip():
        raise HTTPException(status_code=400, detail="LLM API key is required for this provider.")

    body = {
        "model": model,
        "messages": build_llm_messages(request),
        "temperature": 0.1,
        "stream": False,
    }
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "EasyScan_HEP-UI",
    }
    if request.api_key.strip():
        headers["Authorization"] = f"Bearer {request.api_key.strip()}"

    url = normalize_chat_url(provider, request.base_url)
    http_request = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        headers=headers,
        method="POST",
    )
    try:
        with urllib.request.urlopen(http_request, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        error_body = exc.read(3000).decode("utf-8", errors="replace").strip()
        message = error_body or str(exc.reason)
        raise HTTPException(status_code=502, detail=f"LLM API error {exc.code}: {message}") from exc
    except (urllib.error.URLError, TimeoutError, socket.timeout) as exc:
        reason = getattr(exc, "reason", exc)
        raise HTTPException(status_code=502, detail=f"Can not reach LLM API: {reason}") from exc
    except json.JSONDecodeError as exc:
        raise HTTPException(status_code=502, detail="LLM API returned invalid JSON.") from exc

    try:
        content = payload["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError) as exc:
        raise HTTPException(status_code=502, detail="LLM API response did not contain a message.") from exc
    if not isinstance(content, str) or not content.strip():
        raise HTTPException(status_code=502, detail="LLM API returned an empty message.")
    return content


def extract_ini_from_llm_text(text: str) -> str:
    fenced_blocks = re.findall(r"```(?:ini|cfg|text)?\s*(.*?)```", text, flags=re.IGNORECASE | re.DOTALL)
    candidates = fenced_blocks if fenced_blocks else [text]
    allowed_options = {
        "Result folder name",
        "Scan method",
        "Input parameters",
        "Number of points",
        "MCMC walkers",
        "Random seed",
        "Interval of print",
        "Parallel threads",
        "Parallel folder",
        "Program name",
        "Execute command",
        "Command path",
        "Input file",
        "Input variable",
        "Output file",
        "Output variable",
        "Command executor",
        "Time limit in minute",
        "Clean output file",
        "Bound",
        "Gaussian",
        "FreeFormChi2",
        "Histogram",
        "Scatter",
        "Color",
        "Contour",
    }
    for candidate in candidates:
        lines: list[str] = []
        in_ini = False
        for raw_line in candidate.splitlines():
            stripped = raw_line.strip()
            if stripped.lower().startswith("[scan]"):
                in_ini = True
            if not in_ini:
                continue
            if stripped.startswith("[") and stripped.endswith("]"):
                section_name = stripped.strip("[]").lower()
                if section_name != "scan" and not section_name.startswith("program") and section_name not in {"constraint", "plot"}:
                    break
            if (
                stripped
                and not stripped.startswith("#")
                and not (stripped.startswith("[") and stripped.endswith("]"))
                and not raw_line.startswith((" ", "\t"))
            ):
                if ":" not in raw_line:
                    break
                key = raw_line.split(":", 1)[0].strip()
                if key not in allowed_options:
                    break
            lines.append(raw_line.rstrip())
        ini_text = "\n".join(lines).strip()
        if ini_text:
            return ini_text + "\n"
    raise HTTPException(status_code=422, detail="The LLM response did not contain an EasyScan INI file.")


def applescript_string(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def choose_open_path() -> Path:
    if sys.platform == "darwin":
        default_location = applescript_string(str(RUN_CWD))
        script = [
            f'set defaultFolder to POSIX file "{default_location}"',
            'set chosenFile to choose file with prompt "Load EasyScan INI:" default location defaultFolder',
            "POSIX path of chosenFile",
        ]
        command = ["osascript"]
        for line in script:
            command.extend(["-e", line])
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode != 0:
            message = (completed.stderr or completed.stdout or "").strip()
            if "User canceled" in message or completed.returncode == 1:
                raise HTTPException(status_code=400, detail="Load canceled.")
            raise HTTPException(status_code=500, detail=message or "Failed to open file dialog.")
        selected = completed.stdout.strip()
        if not selected:
            raise HTTPException(status_code=400, detail="Load canceled.")
        return Path(selected).expanduser()

    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as exc:
        raise HTTPException(status_code=500, detail="Native open dialog is not available.") from exc

    root = tk.Tk()
    root.withdraw()
    try:
        selected = filedialog.askopenfilename(
            title="Load EasyScan INI",
            initialdir=str(RUN_CWD),
            filetypes=[("INI files", "*.ini"), ("All files", "*.*")],
        )
    finally:
        root.destroy()
    if not selected:
        raise HTTPException(status_code=400, detail="Load canceled.")
    return Path(selected).expanduser()


def choose_save_path(default_name: str = "easyscan_config.ini") -> Path:
    if sys.platform == "darwin":
        default_location = applescript_string(str(RUN_CWD))
        default_file = applescript_string(default_name)
        script = [
            f'set defaultFolder to POSIX file "{default_location}"',
            f'set chosenFile to choose file name with prompt "Save EasyScan INI as:" default name "{default_file}" default location defaultFolder',
            "POSIX path of chosenFile",
        ]
        command = ["osascript"]
        for line in script:
            command.extend(["-e", line])
        completed = subprocess.run(command, capture_output=True, text=True)
        if completed.returncode != 0:
            message = (completed.stderr or completed.stdout or "").strip()
            if "User canceled" in message or completed.returncode == 1:
                raise HTTPException(status_code=400, detail="Save canceled.")
            raise HTTPException(status_code=500, detail=message or "Failed to open save dialog.")
        selected = completed.stdout.strip()
        if not selected:
            raise HTTPException(status_code=400, detail="Save canceled.")
        return Path(selected).expanduser()

    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as exc:
        raise HTTPException(status_code=500, detail="Native save dialog is not available.") from exc

    root = tk.Tk()
    root.withdraw()
    try:
        selected = filedialog.asksaveasfilename(
            title="Save EasyScan INI as",
            defaultextension=".ini",
            initialfile=default_name,
            filetypes=[("INI files", "*.ini"), ("All files", "*.*")],
        )
    finally:
        root.destroy()
    if not selected:
        raise HTTPException(status_code=400, detail="Save canceled.")
    return Path(selected).expanduser()


def default_config() -> EasyScanConfig:
    for template_path in (RUN_CWD / "templates" / "example_random.ini", REPO_ROOT / "templates" / "example_random.ini"):
        if template_path.exists():
            try:
                return parse_easy_scan_template(template_path)
            except Exception:
                return EasyScanConfig()
    return EasyScanConfig()


def initial_config_from_request(request: Request) -> tuple[EasyScanConfig | None, dict[str, str]]:
    path_value = request.query_params.get("config") or os.environ.get("EASYSCAN_UI_INITIAL_CONFIG", "")
    if not path_value:
        return None, {}
    path = resolve_config_path(path_value)
    try:
        path = ensure_allowed(path)
        if not path.is_file():
            return None, {"message": f"Initial INI file not found: {path}", "kind": "failed"}
        return load_config_file(path), {"message": f"Loaded {path}.", "kind": "completed"}
    except Exception as exc:
        return None, {"message": f"Failed to load initial INI: {exc}", "kind": "failed"}



def allowed_roots() -> list[Path]:
    roots = [RUN_CWD, REPO_ROOT, Path.home()]
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
    result_dir = record.result_dir
    if record.status != "completed" or record.return_code not in (0, None):
        return {
            "result_dir": str(result_dir),
            "exists": result_dir.exists(),
            "files": [],
            "plots": [],
            "tables": {},
            "row_count": 0,
            "columns": [],
            "best_row": None,
            "message": "Results are shown only after a successful run.",
        }
    summary = read_results(result_dir)
    for group in ("files", "plots"):
        for item in summary[group]:
            item["url"] = f"/api/files?path={item['path']}"
    return summary


def run_record_payload(record: RunRecord) -> dict[str, Any]:
    return {
        "run_id": record.run_id,
        "status": record.status,
        "return_code": record.return_code,
        "started_at": record.started_at,
        "ended_at": record.ended_at,
        "config_path": str(record.config_path),
        "log_path": str(record.log_path),
        "result_dir": str(record.result_dir),
    }


def persisted_run_payloads() -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    seen = set(runs)
    for config_path in sorted(RUN_CWD.glob("easyscan_*.ini")):
        run_id = config_path.stem.removeprefix("easyscan_")
        if run_id in seen:
            continue
        log_path = config_path.with_suffix(".log")
        try:
            result_dir = resolve_config_path(parse_easy_scan_template(config_path).result_folder)
        except Exception:
            result_dir = RUN_CWD / "result"
        status = "completed" if log_path.exists() else "saved"
        items.append(
            {
                "run_id": run_id,
                "status": status,
                "return_code": None,
                "started_at": config_path.stat().st_mtime,
                "ended_at": log_path.stat().st_mtime if log_path.exists() else None,
                "config_path": str(config_path),
                "log_path": str(log_path),
                "result_dir": str(result_dir),
            }
        )
    return items


def run_worker(record: RunRecord) -> None:
    record.status = "running"
    command = [str(PYTHON_EXE), "-m", "easyscan_hep.cli", str(record.config_path)]
    env = os.environ.copy()
    pythonpath = [str(REPO_ROOT), str(SRC_ROOT)]
    if env.get("PYTHONPATH"):
        pythonpath.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(pythonpath)
    env["EASYSCAN_RESULT_EXISTS_ACTION"] = "replace"
    with record.log_path.open("w", encoding="utf-8") as log_file:
        log_file.write("$ " + " ".join(command) + "\n")
        log_file.flush()
        try:
            record.process = subprocess.Popen(
                command,
                cwd=RUN_CWD,
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                preexec_fn=os.setsid,
            )
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
    initial_config, initial_status = initial_config_from_request(request)
    payload = as_payload(initial_config or default_config())
    static_version = int(max((UI_ROOT / "static" / "app.js").stat().st_mtime, (UI_ROOT / "static" / "styles.css").stat().st_mtime))
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "default_config": json.dumps(payload),
            "initial_config_status": json.dumps(initial_status),
            "static_version": static_version,
        },
    )


@app.post("/api/config")
def preview_config(config: EasyScanConfig) -> dict[str, str]:
    return {"ini": serialize_config(config)}


@app.post("/api/config/import")
def import_config(request: ImportConfigRequest) -> dict[str, Any]:
    try:
        return as_payload(parse_easy_scan_text(request.content))
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to parse INI file: {exc}") from exc


@app.post("/api/config/import/open")
def import_config_from_file() -> dict[str, Any]:
    path = choose_open_path()
    if not path.is_file():
        raise HTTPException(status_code=404, detail="INI file not found.")
    try:
        config = load_config_file(path)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=f"Failed to parse INI file: {exc}") from exc
    return as_payload(config)


@app.post("/api/config/export")
def export_config(config: EasyScanConfig) -> Response:
    return Response(
        content=serialize_config(config),
        media_type="text/plain; charset=utf-8",
        headers={"Content-Disposition": 'attachment; filename="easyscan_config.ini"'},
    )


@app.post("/api/config/export/save")
def export_config_to_file(config: EasyScanConfig) -> dict[str, str]:
    default_name = Path(config.config_file_path).name if config.config_file_path.strip() else "easyscan_config.ini"
    path = choose_save_path(default_name)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(serialize_config(config), encoding="utf-8")
    return {"path": str(path)}


@app.post("/api/config/check")
def check_config(config: EasyScanConfig) -> dict[str, Any]:
    report = check_config_text(serialize_config(config), base_dir=config_base_dir(config))
    report["text"] = format_check_report(report)
    return report


@app.post("/api/config/ai/generate")
def generate_config_with_ai(request: LLMConfigRequest) -> dict[str, Any]:
    if not request.prompt.strip():
        raise HTTPException(status_code=400, detail="Natural-language request is required.")

    llm_text = call_openai_compatible_llm(request)
    ini_text = extract_ini_from_llm_text(llm_text)
    try:
        generated_config = parse_easy_scan_text(ini_text)
    except Exception as exc:
        raise HTTPException(
            status_code=422,
            detail={
                "message": f"LLM output could not be parsed as an EasyScan INI file: {exc}",
                "ini": ini_text,
            },
        ) from exc

    generated_config.config_file_path = request.current_config.config_file_path
    report = check_config_text(ini_text, base_dir=config_base_dir(request.current_config))
    check_text = format_check_report(report)
    if not report.get("ok"):
        raise HTTPException(
            status_code=422,
            detail={
                "message": "AI generated an INI file, but Check Config found errors.",
                "check_text": check_text,
                "ini": ini_text,
            },
        )

    return {
        "config": as_payload(generated_config),
        "ini": ini_text,
        "check": report,
        "check_text": check_text,
    }


@app.get("/api/runs")
def list_runs() -> dict[str, Any]:
    current = [run_record_payload(record) for record in runs.values()]
    history = current + persisted_run_payloads()
    history.sort(key=lambda item: item.get("started_at") or 0, reverse=True)
    return {"runs": history}


@app.post("/api/runs")
def start_run(config: EasyScanConfig) -> dict[str, Any]:
    run_id = time.strftime("%Y%m%d_%H%M%S")
    suffix = 1
    config_path, log_path = config_target_paths(config, run_id)
    while not config.config_file_path.strip() and (config_path.exists() or log_path.exists()):
        suffix += 1
        run_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{suffix}"
        config_path, log_path = config_target_paths(config, run_id)
    run_dir = RUN_CWD
    result_dir = resolve_config_path(config.result_folder)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(serialize_config(config), encoding="utf-8")
    record = RunRecord(run_id, run_dir, config_path, log_path, result_dir)
    runs[run_id] = record
    thread = threading.Thread(target=run_worker, args=(record,), daemon=True)
    thread.start()
    return {"run_id": run_id, "status": record.status, "config_path": str(config_path), "log_path": str(log_path)}


@app.get("/api/runs/{run_id}")
def get_run(run_id: str) -> dict[str, Any]:
    record = runs.get(run_id)
    if record is None:
        config_path, log_path = run_artifact_paths(run_id)
        if not config_path.exists():
            raise HTTPException(status_code=404, detail="Run not found.")
        try:
            result_dir = resolve_config_path(parse_easy_scan_template(config_path).result_folder)
        except Exception:
            result_dir = RUN_CWD / "result"
        record = RunRecord(run_id, RUN_CWD, config_path, log_path, result_dir)
        record.status = "completed"
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
    path: Optional[str] = Query(default=None),
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
