import configparser
import csv
import io
import os
import re
import shutil
from pathlib import Path


SCAN_METHODS = {"ONEPOINT", "ONEPOINTBATCH", "RANDOM", "GRID", "MCMC", "EMCEE", "MULTINEST", "DYNESTY", "POSTPROCESS", "PLOT", "READ"}
NO_NUMBER_OF_POINTS = {"ONEPOINT", "ONEPOINTBATCH", "GRID", "POSTPROCESS", "PLOT", "READ"}
NO_LIKE = {"ONEPOINT", "ONEPOINTBATCH", "RANDOM", "GRID", "POSTPROCESS", "PLOT", "READ"}
PARALLEL_METHODS = {"RANDOM", "GRID", "MCMC", "ONEPOINTBATCH", "MULTINEST", "POSTPROCESS"}
FORBIDDEN_NAMES = {
    "abs",
    "float",
    "int",
    "acos",
    "asin",
    "atan",
    "atan2",
    "ceil",
    "cos",
    "cosh",
    "degrees",
    "e",
    "exp",
    "fabs",
    "floor",
    "fmod",
    "frexp",
    "hypot",
    "ldexp",
    "log",
    "log10",
    "modf",
    "pi",
    "pow",
    "radians",
    "sin",
    "sinh",
    "sqrt",
    "tan",
    "tanh",
    "ONEPOINT",
    "ONEPOINTBATCH",
    "RANDOM",
    "GRID",
    "MCMC",
    "EMCEE",
    "MULTINEST",
    "DYNESTY",
    "POSTPROCESS",
    "PLOT",
    "READ",
    "dwell",
    "probability",
    "-2lnlike",
}


def split_row(row):
    try:
        return [item.strip() for item in next(csv.reader([row], skipinitialspace=True))]
    except csv.Error:
        return [item.strip() for item in row.split(",")]


def rows_for(config, section, option):
    if not config.has_option(section, option):
        return []
    return [line.strip() for line in config.get(section, option).splitlines() if line.strip()]


def parse_scan_method(value, base_dir):
    method = value.strip()
    batch_file = ""
    if "/" in method and method.upper().startswith("ONEPOINT"):
        batch_file = method.split("/", 1)[1]
        method = "ONEPOINTBATCH"
    else:
        method = method.upper()
    if batch_file:
        batch_path = resolve_path(batch_file, base_dir)
    else:
        batch_path = None
    return method, batch_file, batch_path


def resolve_path(value, base_dir):
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def add_duplicate_errors(kind, names, errors):
    seen = {}
    for name in names:
        seen[name] = seen.get(name, 0) + 1
    for name, count in sorted(seen.items()):
        if count > 1:
            errors.append(f'{kind} "{name}" is duplicated {count} times.')


def check_positive_int(config, section, option, errors, warnings, required=False, used=True):
    if not config.has_option(section, option):
        if required:
            errors.append(f'Missing "{option}" in [{section}].')
        return None
    value = config.get(section, option).strip()
    if not used:
        warnings.append(f'"{option}" is present but not used for this scan method.')
    try:
        number = int(value)
    except ValueError:
        errors.append(f'"{option}" must be an integer.')
        return None
    if number < 1:
        errors.append(f'"{option}" must be larger than 0.')
    return number


def check_float(value, context, errors):
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        errors.append(f"{context} must be a number.")
        return False


def check_config_text(text, base_dir=None):
    base_dir = Path(base_dir or os.getcwd()).expanduser().resolve()
    errors = []
    warnings = []
    info = []
    config = configparser.ConfigParser(interpolation=None)
    try:
        config.read_file(io.StringIO(text))
    except configparser.Error as exc:
        return {"ok": False, "errors": [f"Can not parse config file: {exc}"], "warnings": [], "info": []}

    if not config.has_section("scan"):
        return {"ok": False, "errors": ['Missing required section "[scan]".'], "warnings": warnings, "info": info}

    method = ""
    if not config.has_option("scan", "Scan method"):
        errors.append('Missing "Scan method" in [scan].')
    else:
        method, batch_file, batch_path = parse_scan_method(config.get("scan", "Scan method"), base_dir)
        if method not in SCAN_METHODS:
            errors.append(f'"{method}" is not a supported scan method.')
        if batch_file and batch_path and not batch_path.exists():
            errors.append(f'OnePointBatch file does not exist: {batch_file}')
        if method == "ONEPOINTBATCH" and not batch_file:
            warnings.append('ONEPOINTBATCH requires a file in "Scan method", e.g. ONEPOINT/path/to/file.in.')
        if method in {"ONEPOINT", "GRID"}:
            warnings.append(f"{method} is currently marked as not usable in the UI.")

    if not config.has_option("scan", "Result folder name"):
        errors.append('Missing "Result folder name" in [scan].')
    else:
        result_folder = config.get("scan", "Result folder name").strip()
        if " " in result_folder:
            errors.append('"Result folder name" must not contain spaces.')

    check_positive_int(
        config,
        "scan",
        "Number of points",
        errors,
        warnings,
        required=method not in NO_NUMBER_OF_POINTS,
        used=method not in NO_NUMBER_OF_POINTS,
    )
    check_positive_int(config, "scan", "Interval of print", errors, warnings, required=False, used=True)
    threads = check_positive_int(
        config,
        "scan",
        "Parallel threads",
        errors,
        warnings,
        required=False,
        used=True,
    )
    if config.has_option("scan", "Random seed") and method not in {"RANDOM", "MCMC", "EMCEE", "MULTINEST"}:
        warnings.append('"Random seed" is present but not used for this scan method.')
    if threads and threads > 1:
        if method not in PARALLEL_METHODS:
            warnings.append(f"Parallel execution is not used by {method}.")
        if not config.has_option("scan", "Parallel folder"):
            errors.append('Missing "Parallel folder" because "Parallel threads" is larger than 1.')
        else:
            parallel_folder = resolve_path(config.get("scan", "Parallel folder"), base_dir)
            if not parallel_folder.is_dir():
                errors.append(f"Parallel folder does not exist: {parallel_folder}")

    input_names = []
    fixed_names = []
    if not config.has_option("scan", "Input parameters"):
        errors.append('Missing "Input parameters" in [scan].')
    else:
        for row_index, row in enumerate(rows_for(config, "scan", "Input parameters"), start=1):
            parts = split_row(row)
            if len(parts) < 3:
                errors.append(f"Input parameters row {row_index} needs at least ID, Prior, and Value/Min.")
                continue
            name, prior = parts[0], parts[1]
            input_names.append(name)
            if name in FORBIDDEN_NAMES:
                errors.append(f'Input parameter "{name}" uses a forbidden name.')
            if prior.upper() == "FIXED" or method in {"ONEPOINT", "ONEPOINTBATCH"}:
                fixed_names.append(name)
                check_float(parts[2], f'Input parameter "{name}" value', errors)
                if len(parts) > 3:
                    warnings.append(f'Input parameter "{name}" has extra fields ignored for Fixed/one-point mode.')
                continue
            if prior.lower() not in {"flat", "log"}:
                errors.append(f'Input parameter "{name}" prior must be flat, log, or Fixed.')
            if len(parts) < 4:
                errors.append(f'Input parameter "{name}" needs Minimum and Maximum.')
                continue
            check_float(parts[2], f'Input parameter "{name}" minimum', errors)
            check_float(parts[3], f'Input parameter "{name}" maximum', errors)
            if method == "GRID":
                if len(parts) >= 5:
                    try:
                        if int(parts[4]) < 1:
                            errors.append(f'Input parameter "{name}" bins must be larger than 0.')
                    except ValueError:
                        errors.append(f'Input parameter "{name}" bins must be an integer.')
                else:
                    warnings.append(f'Input parameter "{name}" has no GRID bins; EasyScan will use default 10.')
            elif method in {"MCMC", "EMCEE"}:
                if len(parts) < 5:
                    warnings.append(f'Input parameter "{name}" has no {method} interval; EasyScan will use default 10.')
                elif not check_float(parts[4], f'Input parameter "{name}" interval', errors):
                    pass
                if len(parts) < 6:
                    warnings.append(f'Input parameter "{name}" has no {method} initial value; EasyScan will use midpoint.')
                else:
                    check_float(parts[5], f'Input parameter "{name}" initial value', errors)
            elif method in {"RANDOM", "MULTINEST", "DYNESTY"} and len(parts) > 4:
                warnings.append(f'Input parameter "{name}" has extra fields ignored by {method}.')
    add_duplicate_errors("Input parameter", input_names, errors)

    program_sections = sorted(
        [section for section in config.sections() if section.lower().startswith("program")],
        key=lambda item: int(item[7:]) if item[7:].isdigit() else 10**9,
    )
    if method != "PLOT" and not program_sections:
        errors.append('Missing required "[programN]" section.')
    output_names = []
    for section in program_sections:
        if not section[7:].isdigit():
            errors.append(f'Section "[{section}]" should be named like [program1].')
        if not config.has_option(section, "Execute command"):
            errors.append(f'Missing "Execute command" in [{section}].')
        else:
            command = config.get(section, "Execute command").strip().split()
            executable = command[0] if command else ""
            command_path = config.get(section, "Command path", fallback="").strip()
            resolved_command_path = resolve_path(command_path or ".", base_dir)
            if command_path and not resolved_command_path.is_dir():
                errors.append(f"Command path does not exist in [{section}]: {command_path}")
            if executable and not executable.startswith(("/", ".")) and shutil.which(executable) is None:
                warnings.append(f'Executable "{executable}" in [{section}] is not on PATH.')
            elif executable.startswith("./"):
                executable_path = resolved_command_path / executable[2:]
                if not executable_path.exists():
                    errors.append(f'Execute command file does not exist in [{section}]: {executable_path}')
        for option in ("Input file", "Output file"):
            for row in rows_for(config, section, option):
                parts = split_row(row)
                if len(parts) < 2:
                    errors.append(f'"{option}" row in [{section}] needs file ID and path.')
                    continue
                path = resolve_path(parts[1], base_dir)
                if option == "Input file" and not path.exists():
                    errors.append(f"Input file does not exist in [{section}]: {parts[1]}")
                if option == "Output file" and not path.parent.exists():
                    warnings.append(f"Output file parent does not exist in [{section}]: {path.parent}")
        for row in rows_for(config, section, "Output variable"):
            parts = split_row(row)
            if parts:
                output_names.append(parts[0])
    add_duplicate_errors("Output variable", output_names, errors)

    if method not in NO_LIKE:
        if not config.has_section("constraint"):
            errors.append('Missing required section "[constraint]" for likelihood-based scan.')
        elif not (config.has_option("constraint", "Gaussian") or config.has_option("constraint", "FreeFormChi2")):
            errors.append('Section [constraint] needs "Gaussian" or "FreeFormChi2".')

    known_variables = set(input_names) | set(fixed_names) | set(output_names) | {"Chi2", "probability", "-2lnlike", "dwell"}
    if config.has_section("plot"):
        for kind, min_count in (("Histogram", 1), ("Scatter", 2), ("Color", 3), ("Contour", 3)):
            for row in rows_for(config, "plot", kind):
                parts = split_row(row)
                if len(parts) < min_count:
                    errors.append(f'Plot "{kind}" row has too few fields.')
                    continue
                for variable in parts[:min_count]:
                    if variable and variable not in known_variables:
                        warnings.append(f'Plot variable "{variable}" is not defined by input/output/constraint variables.')

    if not errors:
        info.append("Config check passed.")
    return {"ok": not errors, "errors": errors, "warnings": warnings, "info": info}


def check_config_file(path):
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        return {"ok": False, "errors": [f"Config file does not exist: {path}"], "warnings": [], "info": []}
    return check_config_text(path.read_text(encoding="utf-8"), base_dir=Path.cwd())


def format_check_report(report):
    lines = []
    if report.get("errors"):
        lines.append("Errors:")
        lines.extend(f"  - {item}" for item in report["errors"])
    if report.get("warnings"):
        lines.append("Warnings:")
        lines.extend(f"  - {item}" for item in report["warnings"])
    if report.get("info"):
        lines.append("Info:")
        lines.extend(f"  - {item}" for item in report["info"])
    if not lines:
        lines.append("Config check passed.")
    return "\n".join(lines)
