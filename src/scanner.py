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


class SCANNER:
  def __init__(self, name, no_random=False, no_like=False, only_plot=False, use_AccepRate=False, postprocess=False, multinest=False):
    self.name = name.upper()
    self.no_random = no_random
    self.no_like = no_like
    self.use_AccepRate = use_AccepRate
    self.only_plot = only_plot
    self.postprocess = postprocess
    self.multinest = multinest
    if self.only_plot:
      self.postprocess = True
    if self.postprocess:
      self.no_random = True
      self.no_like = True


scanner_random = SCANNER("RANDOM", no_like=True)
scanner_grid = SCANNER("GRID", no_random=True, no_like=True)
scanner_mcmc = SCANNER("MCMC", use_AccepRate=True)
scanner_multinest = SCANNER("MULTINEST", multinest=True)
scanner_postprocess = SCANNER("POSTPROCESS", postprocess=True)
scanner_plot = SCANNER("Plot", only_plot=True)

scanners = [scanner_random, scanner_grid, scanner_mcmc, scanner_multinest, scanner_postprocess, scanner_plot]

names = list(map(lambda x:x.name, scanners))
names_no_random = list(map(lambda x:x.name, filter(lambda y:y.no_random, scanners)))
names_no_like   = list(map(lambda x:x.name, filter(lambda y:y.no_like,   scanners)))
names_only_plot = list(map(lambda x:x.name, filter(lambda y:y.only_plot, scanners)))
names_post = list(map(lambda x:x.name, filter(lambda y:y.postprocess, scanners)))
names_multinest = list(map(lambda x:x.name, filter(lambda y:y.multinest, scanners)))
names_use_AccepRate = list(map(lambda x:x.name, filter(lambda y:y.use_AccepRate, scanners)))


def saveCube(cube, data_file, file_path, num, save_file):
  for ii in cube:
    try:
      float(ii)
      data_file.write(str(ii)+'\t')
    except:
      if os.path.exists(ii):
        # TODO Do not save path. Make sure this!
        if save_file: 
          shutil.copy(ii, os.path.join(file_path, os.path.basename(ii)+"."+num))
      else:
        data_file.write(str(ii)+'\t')
                
    data_file.write('\n')
    data_file.flush()

def printPoint(Numrun, cube, n_dims, inpar, outpar, loglike, Naccept):
    if Numrun == 0:
        af.Info('------------ Start Point ------------')
    else:
        af.Info('------------ Num: %i ------------'%Numrun)
    for i,name in enumerate(inpar):
        af.Info('Input  - %s = %s '%(name,cube[i]))
    af.Info('.................................')
    for i,name in enumerate(outpar):
        outVar = cube[i+n_dims]
        try:
            float(outVar)
            af.Info('Output - %s = %s '%(name,outVar))
        except:
            if '/' not in outVar: af.Info('Output - %s = %s '%(name,outVar))
            
    af.Info('.................................')
    af.Info('loglikelihood   = '+str(loglike))
    if Numrun == 0:
      af.Info('Initial loglike = '+str(-2*loglike))
    af.Info('Accepted Num    = '+str(Naccept))
    af.Info('Total    Num    = '+str(Numrun))

def printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar):
    af.Info('.......... MCMC info ............')
    af.Info('Test     Chi^2 = '+str(Chisq))
    af.Info('Current  Chi^2 = '+str(CurChisq))
    af.Info('Mimimum  Chi^2 = '+str(MinChisq))
    af.Info('Accepted Ratio = '+str(AccRat))
    if FlagTuneR :
        af.Info('StepZize factor= '+str(exp(kcovar)))


def readrun(LogLikelihood,Prior,n_dims,n_params,inpar,outpar,bin_num,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    f_path = os.path.join(outputfiles_basename,"SavedFile")
    if not os.path.exists(f_path):
        os.makedirs(f_path)

    Ploter = ploter.PLOTER()
    Ploter.setPlotPar(outputfiles_basename, 'READ')
    
    for i,name in enumerate(inpar):
        try:
            Ploter._data[name]
        except:
            af.ErrorStop('Input parameter "%s" could not found in ScanInf.txt.'%name)

    for i,name in enumerate(inpar):
        ntotal = len(Ploter._data[name])
        break

    cube = []
    cubePre = []
    # Initialise the cube
    for i in range(n_params): 
        cube.append(0.0)
        cubePre.append(0.0)

    af.Info('Begin read scan ...')
    
    Naccept = 0
    
    for Nrun in range(ntotal):
        iner = Nrun
        for i,name in enumerate(inpar):
            cube[i] = Ploter._data[name][iner]
        
        loglike = LogLikelihood(cube, n_dims, n_params)
        if loglike > af.log_zero:
            Naccept += 1
            saveCube(cube,f_out,f_path,str(Naccept),True)
        #saveCube(cube,f_out2,f_path,str(Naccept),False)
        
        if (Nrun+1)%n_print == 0: printPoint(Nrun+1,cube,n_dims,inpar,outpar,loglike,Naccept)

        #new 20180519 liang
        #cubeProtect = list(cube)
        #if cubePre[n_dims:n_params] == cube[n_dims:n_params]:
        #    for i in range(n_dims, n_params):
        #        cube[i]=af.NaN
        #cube = list(cubeProtect)
        #cubePre = list(cube)

def gridrun(LogLikelihood, Prior, n_dims,n_params,inpar,outpar,bin_num,n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, "ScanResult.txt"),'w') 
    file_path = os.path.join(outputfolder,"SavedFile")
    
    ntotal = 1 # TODO this is where I stopped
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
        loglike = LogLikelihood(cube, n_dims, n_params)

        if loglike > af.log_zero:
            Naccept += 1
            saveCube(cube,data_file,file_path,str(Naccept),True)

        if (Nrun+1)%n_print == 0: printPoint(Nrun+1,cube,n_dims,inpar,outpar,loglike,Naccept)

def randomrun(LogLikelihood, Prior, n_dims, n_params, inpar, outpar, n_live_points, n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, "ScanResult.txt"),'w') # TODO w -> a
    file_path = os.path.join(outputfolder,"SavedFile")
   
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
        loglike = LogLikelihood(cube, n_dims, n_params)
        
        if loglike > af.log_zero:
            Naccept += 1
            saveCube(cube, data_file, file_path, str(Naccept), True)
        
        if (Nrun+1)%n_print == 0: 
            printPoint(Nrun+1, cube, n_dims, inpar, outpar, loglike, Naccept)

def mcmcrun(LogLikelihood,Prior,n_dims,n_params,n_live_points,inpar,outpar,StepSize,AccepRate,FlagTuneR,InitVal,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    f_path = os.path.join(outputfiles_basename,"SavedFile")
    
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
        loglike = LogLikelihood(cube, n_dims, n_params)
        AllOutMCMC = cube.copy()
        AllOutMCMC.append(1)
        #"True" for saving files of initial physical point
        saveCube(AllOutMCMC,f_out2,f_path,'0',True)
        if loglike > af.log_zero / 2.0 : break
        if n_init == 0 : 
            af.WarningNoWait('The initial point is unphysical, it will find the physical initial points randmly.')
        n_init = n_init +1
        if n_init>100:
            af.ErrorStop('Can not find physical initial points with 100 tries.')
        for i in range(n_dims):
            cube[i] = random()
            CurPar[i] = cube[i]

    CurObs=[]
    CurChisq = - 2.0 * loglike
    for i in range(n_params): CurObs.append( cube[i] )
    CurObs.append(0) # mult
    printPoint(0,cube,n_dims,inpar,outpar,loglike,0)

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
            loglike = LogLikelihood(cube, n_dims, n_params)
            AllOutMCMC = cube.copy()
            AllOutMCMC.append(1)
            saveCube(AllOutMCMC,f_out2,f_path,'0',False)
            Chisq = - 2.0 * loglike

        Flag_accept = RangeFlag and (Chisq < CurChisq + 20) 
        if Flag_accept: 
            if CurChisq > Chisq: 
                Flag_accept = True
            else:
                Flag_accept = random() < exp(CurChisq-Chisq) 
        if Flag_accept :
            CurObs[-1]=mult
            #"Naccept+1" due to file of Chisq have covered file of CurChisq
            saveCube(CurObs,f_out,f_path,str(Naccept+1),True)
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
                printPoint(Nrun,cube,n_dims,inpar,outpar,loglike,Naccept)
                printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar)

    # save the last point
    CurObs[-1]=mult
    saveCube(CurObs,f_out,f_path,str(Naccept),True)
