#!/usr/bin/env python3
#########################################################################
"""
       ____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
       |___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
       |___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

    A tool for easily connecting programs to scan physics models.
        
    Author: Yang Zhang and Liangliang Shang
    Web: https://github.com/phyzhangyang/EasyScan_HEP
                                                                     """
##########################################################################

# External modules.
import os,sys
sys.path.append(os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "src"))
# Internal modules.
import initialize
import auxfun as af
import statfun
import scanner
from readin_config   import ReadIn
from scan_controller import SCANINPUT
from constraint      import CONSTRAINT
from ploter          import PLOTER

# Create basic objects
ES         = SCANINPUT()
Programs   = {}
Constraint = CONSTRAINT()
Ploter     = PLOTER()

# Set objects
ProgID     = ReadIn(sys.argv[1], ES, Programs, Constraint, Ploter)

# Write names of parameters into result file
af.WriteResultInf(ES.InPar, ES.FixedPar, ES.OutPar, Constraint.Chi2, ES.getFolderName(), ES.getScanMethod())

# Natural logarithm of likelihood function
def LnLike(cube, ndim, nparams):
    # Pass input value from cube to AllPar
    for i,name in enumerate(ES.InPar):
        ES.AllPar[name]=cube[i]
        
    # Run programs
    for ii in ProgID:
        Programs[ii].WriteInputFile(ES.AllPar)
        Programs[ii].RunProgram()
        Phy = Programs[ii].ReadOutputFile(ES.AllPar, ES.getFolderName())
        # Apply bound
        if Phy: Phy = Programs[ii].ReadBound(ES.AllPar)
        # If the point is unphysical, return log(0)
        if not Phy : return af.log_zero

    # Pass fixed variables to cube
    for i,name in enumerate(ES.FixedPar) :
        cube[i+ndim]   = ES.AllPar[name]    
    # Pass output variables to cube
    for i,name in enumerate(ES.OutPar) :
        cube[i+ndim+len(ES.FixedPar)]   = ES.AllPar[name]

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
if ES.getScanMethod() == af._random:
    scanner.randomrun(
        LnLike = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
	      outpar        = ES.OutPar,
        n_live_points = ES.getPointNum(),
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName())

elif ES.getScanMethod() == af._grid:
    scanner.gridrun(
        LnLike = LnLike,
        Prior         = Prior,
        n_params      = len(ES.AllPar)+len(Constraint.Chi2),
        inpar         = ES.InPar,
        fixedpar      = ES.FixedPar,
	      outpar        = ES.OutPar,
        bin_num       = ES.GridBin,
        n_print       = ES.getPrintNum(),
        outputfolder  = ES.getFolderName())

elif ES.getScanMethod() == af._mcmc:
    scanner.mcmcrun(
        LnLike = LnLike,
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
        outputfolder  = ES.getFolderName())

elif ES.getScanMethod() == af._multinest:
    import pymultinest
    # See https://johannesbuchner.github.io/PyMultiNest/_modules/pymultinest/run.html
    # for more settings
    pymultinest.run(
        LogLikelihood        = LnLike,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(Constraint.Chi2),
        seed                 = ES.getRandomSeed(),
        outputfiles_basename = ES.MNOutputFile,
        n_live_points        = ES.getPointNum(),
        verbose                    = True,
        resume                     = True, #TODO
        importance_nested_sampling = True)

elif ES.getScanMethod() == af._postprocess:
    scanner.postprocessrun(
            LnLike       = LnLike,
            n_params     = len(ES.AllPar)+len(Constraint.Chi2),
            inpar                = ES.InPar,
            fixedpar      = ES.FixedPar,
            outpar               = ES.OutPar,
            bin_num              = ES.GridBin,
            n_print              = ES.getPrintNum(),
            outputfolder = ES.getFolderName())

## recover the modified input file(s) for external programs
if ES.getScanMethod() != af._plot: 
    for ii in Programs: Programs[ii].Recover()

""" Plot """
Ploter.setPlotPar(ES.getFolderName(), ES._ScanMethod)
Ploter.getPlot(ES._ScanMethod)



