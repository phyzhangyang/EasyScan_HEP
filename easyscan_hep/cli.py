"""Command-line entry point for EasyScan_HEP."""

from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
import urllib.parse
import urllib.request
import webbrowser
from pathlib import Path


def source_root() -> Path:
    root = Path(__file__).resolve().parents[1]
    return root


def configure_import_paths() -> Path:
    root = source_root()
    src_root = root / "src"
    if src_root.is_dir() and str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))
    if root.is_dir() and str(root) not in sys.path:
        sys.path.insert(0, str(root))
    os.environ.setdefault("EASYSCAN_ROOT", str(root))
    return root


def ui_is_running(url: str) -> bool:
    try:
        urllib.request.urlopen(url, timeout=1).close()
        return True
    except Exception:
        return False


def open_ui_when_ready(url: str) -> None:
    for _ in range(30):
        if ui_is_running(url):
            webbrowser.open(url)
            return
        time.sleep(1)
    print("Open this URL in your browser: %s" % url)


def resolve_ui_config_path(path_value: str | None) -> Path | None:
    if not path_value:
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    path = path.resolve()
    if not path.is_file():
        print("ERROR: UI config file not found: %s" % path)
        sys.exit(1)
    return path


def run_ui(config_path: str | None = None) -> None:
    configure_import_paths()
    from easyscan_hep.dependency_check import ensure_runtime_requirements

    ensure_runtime_requirements("ui", config_path)
    initial_config = resolve_ui_config_path(config_path)
    url = "http://127.0.0.1:8000/"
    if initial_config is not None:
        os.environ["EASYSCAN_UI_INITIAL_CONFIG"] = str(initial_config)
        url = "%s?config=%s" % (url, urllib.parse.quote(str(initial_config), safe=""))
    os.environ["EASYSCAN_UI_CWD"] = os.getcwd()

    if ui_is_running(url):
        webbrowser.open(url)
        print("EasyScan_HEP Local UI is already running.", flush=True)
        print("Opened %s" % url, flush=True)
        return

    try:
        import uvicorn
    except ImportError:
        print("FastAPI UI dependencies are not installed.")
        print("Install them with:")
        print("  %s -m pip install fastapi uvicorn jinja2 python-multipart" % sys.executable)
        sys.exit(1)

    threading.Thread(target=open_ui_when_ready, args=(url,), daemon=True).start()
    print("Starting EasyScan_HEP Local UI at %s" % url, flush=True)
    print("Keep this terminal open while using the UI.", flush=True)
    print("Press Control-C here to stop the server.", flush=True)
    print("", flush=True)
    uvicorn.run("ui.app:app", host="127.0.0.1", port=8000)


def run_check(argv: list[str]) -> None:
    configure_import_paths()
    from config_checker import format_check_report
    from easyscan_hep.api import check_config

    parser = argparse.ArgumentParser(prog="easyscan --check", description="Check an EasyScan_HEP config without running it.")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--json", action="store_true", dest="json_output")
    parser.add_argument("--base-dir", default=None, help="Base directory for resolving relative paths.")
    parser.add_argument("config", nargs="?")
    args = parser.parse_args(argv[1:])
    if not args.config:
        parser.print_usage()
        sys.exit(1)
    report = check_config(args.config, base_dir=args.base_dir)
    if args.json_output:
        print_json(report)
    else:
        print(format_check_report(report))
    sys.exit(0 if report["ok"] else 1)


def print_json(payload) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False))


def print_main_usage() -> None:
    print("Usage:")
    print("  easyscan CONFIG.ini [--overwrite replace|backup|stop]")
    print("  easyscan --ui [config.ini]")
    print("  easyscan --check CONFIG.ini [--json]")
    print("  easyscan --run CONFIG.ini [--overwrite replace|backup|stop] [--json]")
    print("  easyscan --results RESULT_DIR [--json]")


def run_results(argv: list[str]) -> None:
    from easyscan_hep.results import format_results_summary, read_results

    parser = argparse.ArgumentParser(prog="easyscan --results", description="Read an EasyScan_HEP result directory.")
    parser.add_argument("--results", dest="result_dir")
    parser.add_argument("--json", action="store_true", dest="json_output")
    parser.add_argument("--max-preview", type=int, default=5)
    args = parser.parse_args(argv[1:])
    if not args.result_dir:
        parser.print_usage()
        sys.exit(1)
    summary = read_results(args.result_dir, max_preview=args.max_preview)
    if args.json_output:
        print_json(summary)
    else:
        print(format_results_summary(summary))
    sys.exit(0 if summary["exists"] else 1)


def run_agent_scan(argv: list[str]) -> None:
    from easyscan_hep.dependency_check import ensure_runtime_requirements
    from easyscan_hep.api import run_config
    from easyscan_hep.results import format_results_summary

    parser = argparse.ArgumentParser(prog="easyscan --run", description="Run EasyScan_HEP with machine-readable options.")
    parser.add_argument("--run", dest="run_config", nargs="?")
    parser.add_argument("--json", action="store_true", dest="json_output")
    parser.add_argument("--overwrite", choices=["replace", "backup", "stop"], default="stop")
    parser.add_argument("--cwd", default=None, help="Launch directory. Defaults to the current directory.")
    parser.add_argument("--log", default=None, help="Path to write captured run output.")
    parser.add_argument("config", nargs="?")
    args = parser.parse_args(argv[1:])
    config = args.run_config or args.config
    if not config:
        parser.print_usage()
        sys.exit(1)
    ensure_runtime_requirements("run", config)
    report = run_config(config, cwd=args.cwd, overwrite=args.overwrite, log_path=args.log)
    if args.json_output:
        report = {key: value for key, value in report.items() if key != "stdout"}
        print_json(report)
    else:
        if report["stdout"]:
            print(report["stdout"], end="" if report["stdout"].endswith("\n") else "\n")
        print("Run %s." % ("completed" if report["ok"] else "failed"))
        print("Log file: %s" % report["log_path"])
        if report.get("results"):
            print(format_results_summary(report["results"]))
    sys.exit(0 if report["ok"] else report["return_code"])


def run_install_agent_skill(argv: list[str]) -> None:
    from easyscan_hep.api import AGENT_SKILL_REPOSITORY, install_agent_skill

    parser = argparse.ArgumentParser(prog="easyscan --install-agent-skill", description="Install the EasyScan_HEP agent skill from its standalone repository.")
    parser.add_argument("--install-agent-skill", dest="agent", nargs="?", const="codex", default="codex")
    parser.add_argument("--target", default=None, help="Target skill directory. Defaults to ~/.codex/skills/easyscan-hep for Codex.")
    parser.add_argument("--json", action="store_true", dest="json_output")
    args = parser.parse_args(argv[1:])
    try:
        report = install_agent_skill(args.agent, target=args.target)
        report["ok"] = True
    except Exception as exc:
        report = {"ok": False, "source": AGENT_SKILL_REPOSITORY, "agent": args.agent, "target": args.target or "", "error": str(exc)}
        if args.json_output:
            print_json(report)
        else:
            print("Failed to install EasyScan_HEP agent skill.")
            print("Source: %s" % report["source"])
            print("Error: %s" % report["error"])
        sys.exit(1)
    if args.json_output:
        print_json(report)
    else:
        print("Installed EasyScan_HEP agent skill for %s." % report["agent"])
        print("Source: %s" % report["source"])
        print("Target: %s" % report["target"])
    sys.exit(0)


def apply_overwrite_option(argv: list[str]) -> list[str]:
    cleaned = [argv[0]]
    index = 1
    while index < len(argv):
        arg = argv[index]
        if arg == "--overwrite":
            if index + 1 >= len(argv):
                print('Usage: easyscan CONFIG.ini --overwrite replace|backup|stop')
                sys.exit(1)
            action = argv[index + 1]
            if action not in {"replace", "backup", "stop"}:
                print('ERROR: --overwrite must be "replace", "backup", or "stop".')
                sys.exit(1)
            os.environ["EASYSCAN_RESULT_EXISTS_ACTION"] = action
            index += 2
            continue
        if arg.startswith("--overwrite="):
            action = arg.split("=", 1)[1]
            if action not in {"replace", "backup", "stop"}:
                print('ERROR: --overwrite must be "replace", "backup", or "stop".')
                sys.exit(1)
            os.environ["EASYSCAN_RESULT_EXISTS_ACTION"] = action
            index += 1
            continue
        cleaned.append(arg)
        index += 1
    return cleaned


def run_scan() -> None:
    configure_import_paths()
    from easyscan_hep.dependency_check import ensure_runtime_requirements

    ensure_runtime_requirements("scan", sys.argv[1] if len(sys.argv) > 1 else None)

    # Internal modules are imported only after UI/check handling because
    # initialize.py parses command-line arguments during import.
    import initialize  # noqa: F401
    import auxfun as af
    import statfun
    import scanner
    from constraint import CONSTRAINT
    from ploter import PLOTER
    from readin_config import ReadIn
    from scan_controller import CONTROLLER

    es = CONTROLLER()
    programs = {}
    constraint = CONSTRAINT()
    ploter = PLOTER()

    prog_id = ReadIn(sys.argv[1], es, programs, constraint, ploter)

    af.WriteResultInf(es.InPar, es.FixedPar, es.OutPar, constraint.Chi2, es.getFolderName(), es.getScanMethod())

    def lnlike(cube, ndim, nparams, i_process):
        physical_point = True

        for ii in prog_id:
            for name in list(programs[ii].invar):
                es.AllPar[name] = af.NaN
            for name in list(programs[ii].outvar):
                es.AllPar[name] = af.NaN
            for name in list(programs[ii].boundvar):
                es.AllPar[name] = af.NaN

        af.Info("->\n       Calculating in a new scan\n       <-")
        for i, name in enumerate(es.InPar):
            es.AllPar[name] = cube[i]
            af.Info(f"Scanned parameter: {name} = {cube[i]}")

        for ii in prog_id:
            physical_point = programs[ii].WriteInputFile(es.AllPar, i_process)
            if not physical_point:
                af.Info(f"Scanned parameters are outside domains for some math functions in input variables of {ii}. Go to next scan.")
                break
            af.Info(f"Running {ii}")
            physical_point = programs[ii].RunProgram(i_process)
            if not physical_point:
                break
            physical_point = programs[ii].ReadOutputFile(es.AllPar, es.getFolderName(), i_process)
            if not physical_point:
                af.Info(f"Missing output files or output variables could not be get from output files in {ii}. Go to next scan.")
                break
            physical_point = programs[ii].ReadBound(es.AllPar)
            if not physical_point:
                af.Info(f'Scanned parameters are outside domains for some math functions or excluded in "Bound" of {ii}. Go to next scan.')
                break

        for i, name in enumerate(es.FixedPar):
            cube[i + ndim] = es.AllPar[name]
        for i, name in enumerate(es.OutPar):
            cube[i + ndim + len(es.FixedPar)] = es.AllPar[name]

        if not physical_point:
            return af.log_zero

        loglike = -0.5 * constraint.getChisq(es.AllPar)

        for i, name in enumerate(constraint.Chi2):
            cube[i + ndim + len(es.FixedPar) + len(es.OutPar)] = constraint.Chi2[name]

        return loglike

    def prior(cube, ndim, nparams):
        for i, name in enumerate(es.InPar):
            cube[i] = statfun.prior(cube[i], es.InputPar[name])

    if es.getScanMethod() == af._onepoint:
        scanner.onepointrun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            outputfolder=es.getFolderName(),
        )
    elif es.getScanMethod() == af._onepointbatch:
        scanner.onepointbatchrun(
            LnLike=lnlike,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            scanfile=es.getScanFile(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
        )
    elif es.getScanMethod() == af._random:
        scanner.randomrun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            n_live_points=es.getPointNum(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
        )
    elif es.getScanMethod() == af._grid:
        scanner.gridrun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            bin_num=es.GridBin,
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
        )
    elif es.getScanMethod() == af._bestfit:
        scanner.bestfitrun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            maxiter=es.getPointNum(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
        )
    elif es.getScanMethod() == af._mcmc:
        scanner.mcmcrun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            n_live_points=es.getPointNum(),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            StepSize=es.getStepSize(),
            AccepRate=es.getAccepRate(),
            FlagTuneR=es.getFlagTuneR(),
            InitVal=es.getInitialValue(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
        )
    elif es.getScanMethod() == af._emcee:
        scanner.emceerun(
            LnLike=lnlike,
            Prior=prior,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            n_live_points=es.getPointNum(),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            StepSize=es.getStepSize(),
            InitVal=es.getInitialValue(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
            mcmc_walkers=es.getMCMCWalkers(),
        )
    elif es.getScanMethod() == af._multinest:
        scanner.multinestrun(
            LnLike=lnlike,
            Prior=prior,
            n_dims=len(es.InPar),
            n_params=len(es.AllPar) + len(constraint.Chi2),
            seed=es.getRandomSeed(),
            outputfiles_basename=es.MNOutputFile,
            n_live_points=es.getPointNum(),
            verbose=True,
            resume=af.resume,
            importance_nested_sampling=True,
            num_processes=es.getParallelThreads(),
        )
    elif es.getScanMethod() == af._dynesty:
        scanner.dynestyrun(
            LnLike=lnlike,
            Prior=prior,
            n_dims=len(es.InPar),
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            n_live_points=es.getPointNum(),
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
        )
    elif es.getScanMethod() == af._postprocess:
        scanner.postprocessrun(
            LnLike=lnlike,
            n_params=len(es.AllPar) + len(constraint.Chi2),
            inpar=es.InPar,
            fixedpar=es.FixedPar,
            outpar=es.OutPar,
            n_print=es.getPrintNum(),
            outputfolder=es.getFolderName(),
            num_processes=es.getParallelThreads(),
        )

    if es.getScanMethod() != af._plot:
        for ii in programs:
            if not es.getParallelMode():
                programs[ii].Recover("")
            else:
                for jj in range(es.getParallelThreads()):
                    programs[ii].Recover("p%s_" % str(jj))

    ploter.setPlotPar(es.getFolderName(), es._ScanMethod)
    ploter.getPlot(es._ScanMethod)


def main() -> None:
    configure_import_paths()
    if len(sys.argv) == 1:
        print_main_usage()
        sys.exit(1)
    if sys.argv[1] in ["-h", "--help"]:
        print_main_usage()
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] in ["-ui", "--ui"]:
        if len(sys.argv) > 3:
            print("Usage: easyscan --ui [config.ini]")
            sys.exit(1)
        run_ui(sys.argv[2] if len(sys.argv) == 3 else None)
        return
    if "--install-agent-skill" in sys.argv:
        run_install_agent_skill(sys.argv)
        return
    if "--results" in sys.argv:
        run_results(sys.argv)
        return
    if "--check" in sys.argv:
        run_check(sys.argv)
        return
    if "--run" in sys.argv or "--json" in sys.argv:
        run_agent_scan(sys.argv)
        return
    cleaned_argv = apply_overwrite_option(sys.argv)
    if len(cleaned_argv) == 1:
        print_main_usage()
        sys.exit(1)
    if cleaned_argv[1].startswith("-"):
        print_main_usage()
        sys.exit(1)
    sys.argv = cleaned_argv
    run_scan()


if __name__ == "__main__":
    main()
