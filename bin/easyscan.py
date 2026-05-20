#!/usr/bin/env python3
#########################################################################
r"""
       ____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
       |___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
       |___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

    A tool for easily connecting programs to scan physics models.
        
    Author: Yang Zhang and Liangliang Shang
    Web: https://github.com/phyzhangyang/EasyScan_HEP
                                                                     """
##########################################################################

# External modules.
import os,sys,time,threading,urllib.request,webbrowser

EASYSCAN_ROOT = os.path.split(os.path.split(os.path.realpath(__file__))[0])[0]
sys.path.append(os.path.join(EASYSCAN_ROOT, "src"))


def ui_is_running(url):
    try:
        urllib.request.urlopen(url, timeout=1).close()
        return True
    except Exception:
        return False


def open_ui_when_ready(url):
    for _ in range(30):
        if ui_is_running(url):
            webbrowser.open(url)
            return
        time.sleep(1)
    print("Open this URL in your browser: %s" % url)


def run_ui():
    url = "http://127.0.0.1:8000/"
    os.environ["EASYSCAN_UI_CWD"] = os.getcwd()
    os.chdir(EASYSCAN_ROOT)
    if EASYSCAN_ROOT not in sys.path:
        sys.path.insert(0, EASYSCAN_ROOT)

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


if len(sys.argv) > 1 and sys.argv[1] in ["-ui", "--ui"]:
    run_ui()
    sys.exit(0)

if "--check" in sys.argv:
    from config_checker import check_config_file, format_check_report

    config_args = [arg for arg in sys.argv[1:] if arg != "--check"]
    if len(config_args) != 1:
        print("Usage: ./bin/easyscan.py --check CONFIG.ini")
        sys.exit(1)
    report = check_config_file(config_args[0])
    print(format_check_report(report))
    sys.exit(0 if report["ok"] else 1)

# Internal modules.
import initialize
import auxfun as af
import statfun
import scanner
from readin_config   import ReadIn
from scan_controller import CONTROLLER
from constraint      import CONSTRAINT
from ploter          import PLOTER

# Create basic objects
ES         = CONTROLLER()
Programs   = {}
Constraint = CONSTRAINT()
Ploter     = PLOTER()

# Set objects
ProgID     = ReadIn(sys.argv[1], ES, Programs, Constraint, Ploter)

# Write names of parameters into result file
af.WriteResultInf(ES.InPar, ES.FixedPar, ES.OutPar, Constraint.Chi2, ES.getFolderName(), ES.getScanMethod())

# Natural logarithm of likelihood function
def LnLike(cube, ndim, nparams, i_process):
    PhysicalPoint = True
    
    for ii in ProgID:
        for name in list(Programs[ii].invar):
            ES.AllPar[name] = af.NaN
        for name in list(Programs[ii].outvar):
            ES.AllPar[name] = af.NaN
        for name in list(Programs[ii].boundvar):
            ES.AllPar[name] = af.NaN
    
    af.Info(f'->\n       Calculating in a new scan\n       <-')
    # Pass input value from cube to AllPar
    for i,name in enumerate(ES.InPar):
        ES.AllPar[name]=cube[i]
        af.Info(f'Scanned parameter: {name} = {cube[i]}')
    
    # Run programs
    for ii in ProgID:
        # wirte input file
        PhysicalPoint = Programs[ii].WriteInputFile(ES.AllPar,i_process)
        if not PhysicalPoint:
            af.Info(f'Scanned parameters are outside domains for some math functions in input variables of {ii}. Go to next scan.')
            break
        # command running
        af.Info(f'Running {ii}')
        PhysicalPoint = Programs[ii].RunProgram(i_process)
        if not PhysicalPoint:
            break
        # read output file
        PhysicalPoint = Programs[ii].ReadOutputFile(ES.AllPar, ES.getFolderName(),i_process)
        if not PhysicalPoint:
            af.Info(f'Missing output files or output variables could not be get from output files in {ii}. Go to next scan.')
            break
        # Apply bound
        PhysicalPoint = Programs[ii].ReadBound(ES.AllPar)
        # If the point is unphysical, return log(0)
        if not PhysicalPoint:
            af.Info(f'Scanned parameters are outside domains for some math functions or excluded in "Bound" of {ii}. Go to next scan.')
            break
    
    # Pass fixed variables to cube
    for i,name in enumerate(ES.FixedPar):
        cube[i+ndim]   = ES.AllPar[name]    
    # Pass output variables to cube
    for i,name in enumerate(ES.OutPar):
        cube[i+ndim+len(ES.FixedPar)]   = ES.AllPar[name]
    
    if not PhysicalPoint:
        return af.log_zero

    loglike = - 0.5*Constraint.getChisq(ES.AllPar)

    # Pass constraint likelihood to cube
    for i,name in enumerate(Constraint.Chi2) :
        cube[i+ndim+len(ES.FixedPar)+len(ES.OutPar)]   = Constraint.Chi2[name]

    return loglike 

# Prior function
def Prior(cube, ndim, nparams):
    for i,name in enumerate(ES.InPar): 
        cube[i] = statfun.prior(cube[i],ES.InputPar[name])

# Load corresponding scan method
if ES.getScanMethod() == af._onepoint:
    scanner.onepointrun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
        outpar        = ES.OutPar,
        outputfolder  = ES.getFolderName())

elif ES.getScanMethod() == af._onepointbatch:
    scanner.onepointbatchrun(
        LnLike        = LnLike,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
        outpar        = ES.OutPar,
        scanfile      = ES.getScanFile(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName(),
        num_processes = ES.getParallelThreads())

elif ES.getScanMethod() == af._random:
    scanner.randomrun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
	      outpar        = ES.OutPar,
        n_live_points = ES.getPointNum(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName(),
        num_processes = ES.getParallelThreads())

elif ES.getScanMethod() == af._grid:
    scanner.gridrun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
	      outpar        = ES.OutPar,
        bin_num       = ES.GridBin,
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName(),
        num_processes = ES.getParallelThreads())

elif ES.getScanMethod() == af._bestfit:
    scanner.bestfitrun(
        LnLike       = LnLike,
        Prior        = Prior,
        n_params     = len(ES.AllPar)+len(Constraint.Chi2),
        inpar        = ES.InPar,
        fixedpar     = ES.FixedPar,
        outpar       = ES.OutPar,
        maxiter      = ES.getPointNum(),
        n_print      = ES.getPrintNum(),
        outputfolder = ES.getFolderName())

elif ES.getScanMethod() == af._mcmc:
    scanner.mcmcrun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        n_live_points = ES.getPointNum(),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
        outpar        = ES.OutPar,
        StepSize      = ES.getStepSize(),
        AccepRate     = ES.getAccepRate(),
        FlagTuneR     = ES.getFlagTuneR(),
        InitVal       = ES.getInitialValue(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName(),
        num_processes = ES.getParallelThreads())

elif ES.getScanMethod() == af._emcee:
    scanner.emceerun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        n_live_points = ES.getPointNum(),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
        outpar        = ES.OutPar,
        StepSize      = ES.getStepSize(),
        InitVal       = ES.getInitialValue(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName(),
        num_processes = ES.getParallelThreads())

elif ES.getScanMethod() == af._multinest:
    scanner.multinestrun(
        LnLike               = LnLike,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(Constraint.Chi2),
        seed                 = ES.getRandomSeed(),
        outputfiles_basename = ES.MNOutputFile,
        n_live_points        = ES.getPointNum(),
        verbose              = True,
        resume               = af.resume,
        importance_nested_sampling = True,
        num_processes        = ES.getParallelThreads())


elif ES.getScanMethod() == af._dynesty:
    scanner.dynestyrun(
        LnLike        = LnLike,
        Prior         = Prior,
        n_dims        = len(ES.InPar),
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
        outpar        = ES.OutPar,
        n_live_points = ES.getPointNum(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName())


elif ES.getScanMethod() == af._postprocess:
    scanner.postprocessrun(
            LnLike       = LnLike,
            n_params     = len(ES.AllPar)+len(Constraint.Chi2),
            inpar        = ES.InPar,
            fixedpar     = ES.FixedPar,
            outpar       = ES.OutPar,
            n_print      = ES.getPrintNum(),
            outputfolder = ES.getFolderName(),
            num_processes= ES.getParallelThreads())

## recover the modified input file(s) for external programs
if ES.getScanMethod() != af._plot:
    for ii in Programs:
        if not ES.getParallelMode():
            Programs[ii].Recover("")
        else:
            for jj in range(ES.getParallelThreads()):
                Programs[ii].Recover("p%s_"%str(jj))

""" Plot """
Ploter.setPlotPar(ES.getFolderName(), ES._ScanMethod)
Ploter.getPlot(ES._ScanMethod)
