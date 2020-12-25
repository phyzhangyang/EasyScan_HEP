####################################################################
#    Funtions used for differeny scan methods.                     #
####################################################################

## External modules.
import os,sys,shutil
from random import random, gauss
from math import exp
# Internal modules
import auxfun as af
import ploter

_random = "RANDOM"
_grid = "GRID"
_mcmc = "MCMC"
_multinest = "MULTINEST"
_postprocess = "POSTPROCESS"
_plot = "PLOT"

_all = [_random, _grid, _mcmc, _multinest, _postprocess, _plot]
_no_random = [_grid, _postprocess, _plot]
_no_like   = [_random, _grid, _postprocess, _plot]
_post = [_postprocess, _plot]

def saveCube(cube, data_file, file_path, num, save_file):
  result = []
  for ii in cube:
    try:
      float(ii)
      result.append(str(ii))
    except:
      if os.path.exists(ii):
        if save_file: 
          shutil.copy(ii, os.path.join(file_path, os.path.basename(ii)+"."+num))
        result.append(str(af.NaN))
      else:
        result.append(str(ii))
                
  data_file.write(','.join(result)+'\n')
  data_file.flush()

def printPoint(Numrun, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept):
    if Numrun == 0:
        af.Info('------------ Start Point ------------')
    else:
        af.Info('------------ Num: %i ------------'%Numrun)
    for i,name in enumerate(inpar):
        af.Info('Input  - %s = %s '%(name,cube[i]))
    for i,name in enumerate(fixedpar):
        af.Info('Input  - %s = %s '%(name,cube[i+n_dims]))
    for i,name in enumerate(outpar):
        outVar = cube[i+n_dims+len(fixedpar)]
        try:
          float(outVar)
          af.Info('Output - %s = %s '%(name,outVar))
        except:
          if '/' not in outVar: 
            af.Info('Output - %s = %s '%(name,outVar))

    af.Info('LnLike   = '+str(lnlike))
    if Numrun == 0:
      af.Info('Initial lnlike = '+str(-2*lnlike))
    af.Info('Accepted Num    = '+str(Naccept))
    af.Info('Total    Num    = '+str(Numrun))

def printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar):
    af.Info('........... MCMC info .............')
    af.Info('Test     Chi^2 = '+str(Chisq))
    af.Info('Current  Chi^2 = '+str(CurChisq))
    af.Info('Mimimum  Chi^2 = '+str(MinChisq))
    af.Info('Accepted Ratio = '+str(AccRat))
    if FlagTuneR :
        af.Info('StepZize factor= '+str(exp(kcovar)))


def postprocessrun(LnLike, n_params, inpar, fixedpar, outpar, bin_num,n_print,outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a') 
    file_path = os.path.join(outputfolder,"SavedFile")
    
    n_dims = len(inpar)
        
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    # Read data using ploter
    Data = ploter.PLOTER()
    Data.setPlotPar(outputfolder, _plot, postprocess=True)
    data =Data._data
    if not Data.checkPar([i for i in inpar],n_dims, section_name="scan"):
      af.ErrorStop('Can not postprocess it.')
      
    ntotal = data.shape[0]
    # Initialise cube
    cube = [af.NaN] * n_params

    af.Info('Begin read scan ...')
    Naccept = 0
    for Nrun in range(ntotal):
        for i,name in enumerate(inpar):
            cube[i] = data[name][Nrun]
            
        lnlike = LnLike(cube, n_dims, n_params)
        if lnlike > af.log_zero:
            Naccept += 1
            saveCube(cube,data_file,file_path,str(Naccept),True)
        
        if (Nrun+1)%n_print == 0:
            printPoint(Nrun+1, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept)
            

def gridrun(LnLike, Prior, n_params, inpar, fixedpar, outpar, bin_num, n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a') 
    file_path = os.path.join(outputfolder,"SavedFile")
    
    n_dims = len(inpar)
    
    ntotal = 1
    cube = []
    interval = {}
    for i,name in enumerate(inpar):
        interval[name] = 1.0 / bin_num[name]
        bin_num[name] += 1
        ntotal     *= bin_num[name]
    for i in range(n_params): 
        cube.append(0)

    af.Info('Begin grid scan ...')
    
    Naccept = 0

    for Nrun in range(ntotal):
        iner = 1
        for i,name in enumerate(inpar):
            cube[i] = ( int(Nrun/iner) )%bin_num[name] * interval[name]
            iner   *= bin_num[name]

        for i,name in enumerate(outpar):
            cube[i+n_dims] = af.NaN
        
        Prior(cube, n_dims, n_params)
        lnlike = LnLike(cube, n_dims, n_params)

        if lnlike > af.log_zero:
            Naccept += 1
            saveCube(cube,data_file,file_path,str(Naccept),True)

        if (Nrun+1)%n_print == 0: 
          printPoint(Nrun+1, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept)

def randomrun(LnLike, Prior, n_params, inpar, fixedpar, outpar, n_live_points, n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a')
    file_path = os.path.join(outputfolder,"SavedFile")
   
    n_dims = len(inpar)
    
    cube = []
    # Initialise the cube
    for i in range(n_params): 
        cube.append(0.0)

    af.Info('Begin random scan ...')
    Naccept = 0
    for Nrun in range(n_live_points) :
        for j in range(n_dims):
          cube[j] = random()
        
        Prior(cube, n_dims, n_params)
        lnlike = LnLike(cube, n_dims, n_params)
        
        if lnlike > af.log_zero:
            Naccept += 1
            saveCube(cube, data_file, file_path, str(Naccept), True)
        
        if (Nrun+1)%n_print == 0: 
            printPoint(Nrun+1, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept)

def mcmcrun(LnLike, Prior, n_params, n_live_points, inpar, fixedpar, outpar, StepSize, AccepRate, FlagTuneR, InitVal, n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a')
    all_data_file = open(os.path.join(outputfolder, af.ResultFile_MCMC),'a')
    file_path = os.path.join(outputfolder,"SavedFile")
        
    n_dims = len(inpar)
        
    # Initialise the cube
    cube = []
    for i in range(n_params):
        cube.append(0)

    covar = [] # the sigma of gauss distribution, normalized to 1
    par  = []  # test par, normalized to 1
    CurPar=[]  # current par, normalized to 1 
    for i,name in enumerate(inpar):
        covar.append(StepSize[name])
        cube[i] = InitVal[name]
        par.append(cube[i])
        CurPar.append( cube[i] )
    n_init = 0
    while True:
        Prior(cube, n_dims, n_params) # normalized to cube to real value
        lnlike = LnLike(cube, n_dims, n_params)
        AllOutMCMC = cube.copy()
        AllOutMCMC.append(1)
        #"True" for saving files of initial physical point
        saveCube(AllOutMCMC, all_data_file, file_path, '0', False)
        if lnlike > af.log_zero / 2.0 : break
        if n_init == 0 : 
            af.WarningNoWait('The initial point is unphysical, it will find the physical initial points randmly.')
        n_init = n_init +1
        if n_init>100:
            af.ErrorStop('Can not find physical initial points with 100 tries.')
        for i in range(n_dims):
            cube[i] = random()
            CurPar[i] = cube[i]

    CurObs=[]
    CurChisq = - 2.0 * lnlike
    for i in range(n_params): CurObs.append( cube[i] )
    CurObs.append(0) # mult
    printPoint(0, cube, n_dims, inpar, fixedpar, outpar, lnlike, 0)

    # Initialize the MCMC parameters
    MinChisq = CurChisq
    Chisq = CurChisq
    Nrun= 0
    Naccept = 0
    Nout=0
    mult = 1
    kcovar = 0 
    while Naccept < n_live_points:

        #Nrun += 1
        RangeFlag = True
        for j in range(n_dims):
            par[j] = gauss(CurPar[j],exp(kcovar)*covar[j]) # normalized to 1 
            #rd = random()
            #par[j] = CurPar[j] + covar[j] * (0.5-rd)*2
        if max(par)>1 or min(par)<0 :
            RangeFlag = False
            Nout = Nout +1
            if Nout%100 == 0: 
                af.WarningNoWait("Too many points out of range!")
        else:
            Nrun += 1
            Nout=0
            for i in range(n_dims): cube[i] = par[i]
            Prior(cube, n_dims, n_params)
            lnlike = LnLike(cube, n_dims, n_params)
            AllOutMCMC = cube.copy()
            AllOutMCMC.append(1)
            saveCube(AllOutMCMC, all_data_file, file_path, '0', False)
            Chisq = - 2.0 * lnlike

        Flag_accept = RangeFlag and (Chisq < CurChisq + 20) 
        if Flag_accept: 
            if CurChisq > Chisq: 
                Flag_accept = True
            else:
                Flag_accept = random() < exp(CurChisq-Chisq) 
        if Flag_accept :
            CurObs[-1]=mult
            #"Naccept+1" due to file of Chisq have covered file of CurChisq
            saveCube(CurObs, data_file, file_path, str(Naccept+1), True)
            CurChisq = Chisq
            for i in range(n_params): CurObs[i]   = cube[i]
            for i in range(n_dims):   CurPar[i]   = par[i]
            
            if Chisq < MinChisq : MinChisq = Chisq
            Naccept += 1
            mult = 1
        else:
            if RangeFlag:
                mult +=1

        AccRat = float(Naccept)/float(Nrun)

        if FlagTuneR and Nrun < 1000: 
            kcovar = kcovar + 1.0/(float(Nrun)**0.7)*(AccRat - AccepRate)
        else: 
            kcovar = 1

        if Nrun%n_print == 0: 
            if RangeFlag:
                printPoint(Nrun, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept)
                printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar)

    # save the last point
    CurObs[-1]=mult
    saveCube(CurObs, data_file, file_path, str(Naccept), True)
