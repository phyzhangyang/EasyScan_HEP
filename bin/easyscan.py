#!/usr/bin/env python
#########################################################################
"""
       ____ ____ ____ _   _ ____ ____ ____ _  _ _  _ ____ ___
       |___ |__| [__   \_/  [__  |    |__| |\ | |__| |___ |__]
       |___ |  | ___]   |   ___] |___ |  | | \| |  | |___ |

    An Easy-to-use tool providing a comfortable way connecting programs 
    to Scan the parameter space for high energy physics(HEP) models.
        
    Author: Junjie Cao, Liangliang Shang, Jin Min Yang and Yang Zhang
    Web: http://easyscanhep.hepforge.org
                                                                     """
##########################################################################

## External modules.
import os,sys,math
sys.path.append(os.path.join(os.path.split(os.path.split(os.path.realpath(__file__))[0])[0], "src"))
## Internal modules.
import init       as ip
import mainfun    as mf
import statfun    as stat
import readin

# define basic class object
ES       = mf.EasyScanInput()
Programs ={}
CS       = mf.constraint()
Ploter   = mf.plot()
ProgID   = readin.ReadIn(sys.argv[1],ES,Programs,CS,Ploter)


if ES.getScanMethod() != 'PLOT':
    ip.WriteResultInf(ES.InPar,ES.OutPar,ES.getFileName(),ES.getScanMethod())

# logarithm of likelihood function
def LogLikelihood(cube, ndim, nparams):
    # pass the input value from cube to InPar
    for i,name in enumerate(ES.InPar) :
        ES.InPar [name]=cube[i]
        ES.AllPar[name]=cube[i]
    # Run each programs
    for ii in ProgID:
        if Programs[ii].getRunFlag(ES.AllPar):
            Programs[ii].WriteInputFile(ES.AllPar)
            Programs[ii].RunProgram()
            Phy = Programs[ii].ReadOutputFile(ES.AllPar,ES.getFileName())
        else:
            Phy = Programs[ii].SetOutput(ES.AllPar)
        # if the point is unphysical, return log(0)
        if not Phy : return ip.log_zero
    # pass the output value to AllPar
    for i,name in enumerate(ES.OutPar) :
        cube[i+ndim]   = ES.AllPar[name]

    return - CS.getChisq(ES.AllPar)/2.0


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
            ip.ErrorStop( 'Not ready. Only "flat" and "log" prior can be used.' )

## Load corresponding scan method
if ES.getScanMethod() == 'RANDOM':
    from scanmanner import randomrun
    ResultFile = 'RandomData.txt'
    randomrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar),
        n_live_points        = ES.getPointNum(),
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile )

elif ES.getScanMethod() == 'MCMC':
    from scanmanner import mcmcrun
    ResultFile = 'MCMCData.txt'
    mcmcrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar),
        n_live_points        = ES.getPointNum(),
        inpar                = ES.InPar,
        outpar               = ES.OutPar,
        StepSize             = ES.getStepSize(),
        AccepRate            = ES.getAccepRate(),
        FalgTune             = ES.getFalgTuneR(),
        InitVal              = ES.getInitialValue(),
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile)

elif ES.getScanMethod() == 'MULTINEST':
    import pymultinest
    ResultFile = 'MultiNestData/ev.dat'
    pymultinest.run(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar),
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
    ResultFile = 'GridData.txt'
    gridrun(
        LogLikelihood        = LogLikelihood,
        Prior                = Prior,
        n_dims               = len(ES.InPar),
        n_params             = len(ES.AllPar),
        inpar                = ES.InPar,
        bin_num              = ES.GridBin,
        n_print              = ES.getPrintNum(),
        outputfiles_basename = ES.getFileName(),
        outputfiles_filename = ResultFile )

elif ES.getScanMethod() == 'READ':
    from scanmanner import readrun
    ResultFile = 'ReadData.txt'
    readrun(
            LogLikelihood        = LogLikelihood,
            Prior                = Prior,
            n_dims               = len(ES.InPar),
            n_params             = len(ES.AllPar),
            inpar                = ES.InPar,
            bin_num              = ES.GridBin,
            n_print              = ES.getPrintNum(),
            outputfiles_basename = ES.getFileName(),
            outputfiles_filename = ResultFile )

## recover the modified input file(s) for external programs
if ES.getScanMethod() != 'PLOT':
    for ii in Programs: Programs[ii].Recover()
#    ip.WriteResultInf(ES.InPar,ES.OutPar,ES.getFileName(),ResultFile,ES.getScanMethod())


""" Plot """
Ploter.setPlotPar(ES.getFileName())
Ploter.getPlot()



