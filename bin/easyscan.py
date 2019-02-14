#!/usr/bin/env python
#########################################################################
"""
       ____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
       |___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
       |___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

    An Easy-to-use tool providing a comfortable way connecting programs 
    to Scan the parameter space for high energy physics(HEP) models.
        
    Author: Liangliang Shang and Yang Zhang
    Web: http://easyscanhep.hepforge.org
                                                                     """
##########################################################################

## External modules.
import os,sys,math
sys.path.append(os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "src"))
## Internal modules.
import init       as sf
import mainfun    as mf
import statfun    as stat
import readin

# define basic class object
ES       = mf.EasyScanInput()
Programs = {}
CS       = mf.constraint()
Ploter   = mf.plot()
ProgID   = readin.ReadIn(sys.argv[1],ES,Programs,CS,Ploter)

## new 20180416 liang
if ES.getScanMethod() == 'RANDOM':
    ResultFile = 'RandomData.txt'
elif ES.getScanMethod() == 'MCMC':
    ResultFile = 'MCMCData.txt'
elif ES.getScanMethod() == 'MULTINEST':
    ResultFile = 'MultiNestData/.txt'
elif ES.getScanMethod() == 'GRID':
    ResultFile = 'GridData.txt'
elif ES.getScanMethod() == 'READ':
    ResultFile = 'ReadData.txt'
if ES.getScanMethod() != 'PLOT':
    sf.WriteResultInf(ES.InPar,ES.OutPar,CS.Chi2,ES.getFileName(),ES.getScanMethod(), ResultFile)

# logarithm of likelihood function
def LogLikelihood(cube, ndim, nparams):
    # pass the input value from cube to InPar
    for i,name in enumerate(ES.InPar) :
        ES.InPar [name]=cube[i]
        ES.AllPar[name]=cube[i]
    # Run each programs
    # ES.AllPar is a dictionary involving variables and their values of scanning parameters and output variables. 
    for ii in ProgID:
        if Programs[ii].getRunFlag(ES.AllPar):
            Programs[ii].WriteInputFile(ES.AllPar)
            # new 20180416 liang
            Programs[ii].RemoveOutputFile()
            Programs[ii].RunProgram()
            Phy = Programs[ii].ReadOutputFile(ES.AllPar,ES.getFileName())
            if Phy:
                Phy = Programs[ii].ReadBound(ES.AllPar)
        else:
            Phy = Programs[ii].SetOutput(ES.AllPar)
        # if the point is unphysical, return log(0)
        if not Phy : return sf.log_zero

    # pass all output variables with values to "cube" for using in each scan method in "scanmanner.py"
    for i,name in enumerate(ES.OutPar) :
        cube[i+ndim]   = ES.AllPar[name]

    loglike = - CS.getChisq(ES.AllPar)/2.0

    ## new 20180428 liang
    for i,name in enumerate(CS.Chi2) :
        cube[i+ndim+len(ES.OutPar)]   = CS.Chi2[name]

    return loglike 


def Prior(cube, ndim, nparams):
    for i,name in enumerate(ES.InPar):
        if ES.InputPar[name][1].lower() == 'flat':
            min = float(ES.InputPar[name][2])
            max = float(ES.InputPar[name][3])
            cube[i] = cube[i] * (max - min) + min
        elif ES.InputPar[name][1].lower() == 'log':
            min = math.log10(float(ES.InputPar[name][2]))
            max = math.log10(float(ES.InputPar[name][3]))
            cube[i] = 10.0**(cube[i]*(max - min) + min )
        else:
            sf.ErrorStop( 'Not ready. Only "flat" and "log" prior can be used.' )

## Load corresponding scan method
if ES.getScanMethod() == 'RANDOM':
    from scanmanner import randomrun
    randomrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(CS.Chi2),
#        inpar                = ES.InPar,
#	outpar               = ES.OutPar,
        n_live_points        = ES.getPointNum(),
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile )

elif ES.getScanMethod() == 'MCMC':
    from scanmanner import mcmcrun
    mcmcrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(CS.Chi2),
        n_live_points        = ES.getPointNum(),
        inpar                = ES.InPar,
        outpar               = ES.OutPar,
        StepSize             = ES.getStepSize(),
        AccepRate            = ES.getAccepRate(),
        FlagTuneR            = ES.getFlagTuneR(),
        InitVal              = ES.getInitialValue(),
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile)

elif ES.getScanMethod() == 'MULTINEST':
    import pymultinest
    pymultinest.run(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(CS.Chi2),
        seed                 = ES.getRandomSeed(),
        outputfiles_basename = ES.MNOutputFile,
        n_live_points        = ES.getPointNum(),
        n_clustering_params        = 2,
        wrapped_params             = None,
        multimodal                 = True,
        const_efficiency_mode      = False,
        evidence_tolerance         = 1.0,
        sampling_efficiency        = 2.0,
        n_iter_before_update       = 1,
        null_log_evidence          = -1e+90,
        max_modes                  = 5,
        verbose                    = True,
        resume                     = False, #!!!!!!!!!
        context                    = 0,
        importance_nested_sampling = True)

elif ES.getScanMethod() == 'GRID':
    from scanmanner import gridrun
    gridrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar)+len(CS.Chi2),
        inpar                = ES.InPar,
#	outpar               = ES.OutPar,
        bin_num              = ES.GridBin,
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile )

elif ES.getScanMethod() == 'READ':
    from scanmanner import readrun
    readrun(
            LogLikelihood        = LogLikelihood,
            Prior                = Prior,
            n_dims               = len(ES.InPar),
            n_params             = len(ES.AllPar)+len(CS.Chi2),
            inpar                = ES.InPar,
            bin_num              = ES.GridBin,
            n_print              = ES.getPrintNum(),
            outputfiles_basename = ES.getFileName(),
            outputfiles_filename = ResultFile )

## recover the modified input file(s) for external programs
if ES.getScanMethod() != 'PLOT':
    for ii in Programs: Programs[ii].Recover()
#    sf.WriteResultInf(ES.InPar,ES.OutPar,CS.Chi2,ES.getFileName(),ResultFile,ES.getScanMethod())


""" Plot """
Ploter.setPlotPar(ES.getFileName(), ES._ScanMethod)
Ploter.getPlot(ES._ScanMethod)



