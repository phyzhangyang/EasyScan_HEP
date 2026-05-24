"""Command-line entry point for EasyScan_HEP."""

from __future__ import annotations

import os
import sys
import threading
import time
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


def run_ui() -> None:
    configure_import_paths()
    url = "http://127.0.0.1:8000/"
    os.environ["EASYSCAN_UI_CWD"] = os.getcwd()

    if ui_is_running(url):
        webbrowser.open(url)
        print("EasyScan_HEP Local UI is already running.")
        print("Opened %s" % url)
        return

    try:
        import uvicorn
    except ImportError:
        print("FastAPI UI dependencies are not installed.")
        print("Install them with:")
        print("  %s -m pip install fastapi uvicorn jinja2 python-multipart" % sys.executable)
        sys.exit(1)

    threading.Thread(target=open_ui_when_ready, args=(url,), daemon=True).start()
    print("Starting EasyScan_HEP Local UI at %s" % url)
    print("Keep this terminal open while using the UI.")
    print("Press Control-C here to stop the server.")
    print("")
    uvicorn.run("ui.app:app", host="127.0.0.1", port=8000)


def run_check(argv: list[str]) -> None:
    configure_import_paths()
    from config_checker import check_config_file, format_check_report

    config_args = [arg for arg in argv[1:] if arg != "--check"]
    if len(config_args) != 1:
        print("Usage: easyscan --check CONFIG.ini")
        sys.exit(1)
    report = check_config_file(config_args[0])
    print(format_check_report(report))
    sys.exit(0 if report["ok"] else 1)


def run_scan() -> None:
    configure_import_paths()

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
    if len(sys.argv) > 1 and sys.argv[1] in ["-ui", "--ui"]:
        run_ui()
        return
    if "--check" in sys.argv:
        run_check(sys.argv)
        return
    run_scan()


if __name__ == "__main__":
    main()
