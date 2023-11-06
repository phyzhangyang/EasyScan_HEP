####################################################################
#    Funtions used for differeny scan methods.                     #
####################################################################

## External modules.
import os,sys,shutil
from random import random, gauss
from numpy import zeros, sum
from math import exp
# Internal modules
import auxfun as af
import ploter
import time
import multiprocessing
multiprocessing.set_start_method('fork')
    
lock = multiprocessing.Lock()

def run_processes(processes):
    # Start all subprocesses
    for p in processes:
        p.start()
    
    # Wait for all subprocesses to finish
    for p in processes:
        p.join()
    
    af.Info('All processes finished')

def getFilelength(datafile):
  with open(datafile, 'r') as f:
    num_lines = sum(1 for line in f)
  return num_lines
  
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
  with lock:
      data_file.write(','.join(result)+'\n')
      data_file.flush()

def printPoint(Numrun, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept, i_process=''):
  with lock:
    if Numrun == 1:
        af.Info('------------ Start point -%s-------'%i_process)
    else:
        af.Info('------------ Num: %i -%s-----------'%(Numrun,i_process))
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
    if Numrun == 1:
      af.Info('Initial lnlike = '+str(-2*lnlike))
    af.Info('Accepted Num    = '+str(Naccept))
    af.Info('Total    Num    = '+str(Numrun))

def printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar):
  with lock:
    af.Info('........... MCMC info .............')
    af.Info('Test     Chi^2 = '+str(Chisq))
    af.Info('Current  Chi^2 = '+str(CurChisq))
    af.Info('Mimimum  Chi^2 = '+str(MinChisq))
    af.Info('Accepted Ratio = '+str(AccRat))
    if FlagTuneR :
        af.Info('StepZize factor= '+str(exp(kcovar)))

def postprocessrun(LnLike, n_params, inpar, fixedpar, outpar, n_print, outputfolder, num_processes):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a') 
    file_path = os.path.join(outputfolder,"SavedFile")
    
    n_dims = len(inpar)
        
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    # Read data using ploter
    Data = ploter.PLOTER()
    Data.setPlotPar(outputfolder, af._plot, postprocess=True)
    data =Data._data
    if not Data.checkPar([i for i in inpar],n_dims, section_name="scan"):
      af.ErrorStop('Can not postprocess it.')
      
    ntotal = data.shape[0]
    # Initialise cube
    cube = [af.NaN] * n_params

    af.Info('Begin read scan ...')
    
    def per_run(i_process, i_start, i_end):
        N_Accept = 0
        for Nrun in range(i_start, i_end) :
            for i,name in enumerate(inpar):
                cube[i] = data[name][Nrun]
            
            lnlike = LnLike(cube, n_dims, n_params,i_process)
            if lnlike > af.log_zero:
                N_Accept += 1
                saveCube(cube,data_file,file_path,i_process+str(N_Accept),True)
        
            if (Nrun+1)%n_print == 0:
                printPoint(Nrun+1-i_start, cube, n_dims, inpar, fixedpar, outpar, lnlike, N_Accept, i_process)

    if num_processes == 1:
        per_run("")
        return
        
    # Create subprocesses
    processes = []
    i_end = 0
    for ii in range(num_processes):
        i_start = i_end
        i_end += af.divide_jobs(ntotal, num_processes, ii)
        af.Info('p%i process: %i >>>> %i '%(ii, i_start, i_end))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(ii),i_start,i_end))
        processes.append(p)
    
    run_processes(processes)
    

def onepointrun(LnLike, Prior, n_params, inpar, fixedpar, outpar, outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a')
    file_path = os.path.join(outputfolder,"SavedFile")
        
    n_dims = len(inpar) # see ReadIn() in src/readin_config.py
        
    # Initialise cube
    cube = [af.NaN] * n_params

    af.Info('Begin onepoint mode ...')
    Naccept = 0  

    Prior(cube, n_dims, n_params) # normalized to cube to real value
    lnlike = LnLike(cube, n_dims, n_params)

    if lnlike > af.log_zero:   
        Naccept = Naccept +1
    if Naccept == 0: 
        af.WarningNoWait('The initial point is unphysical.')

    printPoint(1, cube, n_dims, inpar, fixedpar, outpar, lnlike, Naccept)
    saveCube(cube, data_file, file_path, str(Naccept), True)

def onepointbatchrun(LnLike, n_params, inpar, fixedpar, outpar, scanfile, n_print, outputfolder):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a') 
    file_path = os.path.join(outputfolder, "SavedFile")
   
    n_dims = len(inpar) # see ReadIn() in src/readin_config.py 
        
    if not os.path.exists(file_path):
        os.makedirs(file_path)

    # Read data using ploter
    Data = ploter.PLOTER()
    Data.setPlotPar(outputfolder, af._onepointbatch, ScanFile=scanfile, Plot=False)
    data =Data._data
    if not Data.checkPar([i for i in inpar], n_dims, section_name="scan", severe=True):
        af.ErrorStop('Can not postprocess it.')
    
    ntotal = data.shape[0]
    # Initialise cube
    cube = [af.NaN] * n_params

    af.Info('Begin one point batch scan ...')
    
    def per_run(i_process, i_start, i_end):
        N_Accept = 0
        for Nrun in range(i_start, i_end) :
            for i,name in enumerate(inpar):
                cube[i] = data[name][Nrun]
            
            lnlike = LnLike(cube, n_dims, n_params,i_process)
            if lnlike > af.log_zero:
                N_Accept += 1
                saveCube(cube,data_file,file_path,i_process+str(N_Accept),True)
        
            if (Nrun+1)%n_print == 0:
                printPoint(Nrun+1-i_start, cube, n_dims, inpar, fixedpar, outpar, lnlike, N_Accept, i_process)

    if num_processes == 1:
        per_run("")
        return
        
    # Create subprocesses
    processes = []
    i_end = 0
    for ii in range(num_processes):
        i_start = i_end
        i_end += af.divide_jobs(ntotal, num_processes, ii)
        af.Info('p%i process: %i >>>> %i '%(ii, i_start, i_end))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(ii),i_start,i_end))
        processes.append(p)
    
    run_processes(processes)
        
def gridrun(LnLike, Prior, n_params, inpar, fixedpar, outpar, bin_num, n_print, outputfolder,num_processes):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a') 
    file_path = os.path.join(outputfolder,"SavedFile")
    
    n_dims = len(inpar)
    
    ntotal = 1
    interval = {}
    for i,name in enumerate(inpar):
        try:
            interval[name] = 1.0 / bin_num[name]
        except ZeroDivisionError:
            af.ErrorStop('The number of intervals in grid scanning could not be zero.')
        bin_num[name] += 1
        ntotal     *= bin_num[name]
    # Initialise cube
    cube = [af.NaN] * n_params

    af.Info('Begin grid scan ...')
    
    Naccept = zeros(num_processes)
    if af.resume:
        try:
            for i_process in range(num_processes):
                Naccept[i_process] = int(open(os.path.join(outputfolder, "Nrun.txt"),'r').read().strip()) + 1
        except: 
            af.ErrorStop('Can not use resume mode because of no p[i]_Nrun.txt.')

        if sum(Naccept) >= ntotal:
            af.ErrorStop('There are already %s living samples in the data file.'%Naccept)
          
    def per_run(i_process, i_start, i_end):
        i_start_ = i_start
        for Nrun in range(i_start, i_end) :
            iner = 1
            for i,name in enumerate(inpar):
                cube[i] = ( int(Nrun/iner) )%bin_num[name] * interval[name]
                iner   *= bin_num[name]
            
            for i,name in enumerate(outpar):
                cube[i+n_dims] = af.NaN
            
            
            Prior(cube, n_dims, n_params)
            lnlike = LnLike(cube, n_dims, n_params, i_process)

            if lnlike > af.log_zero:
                i_start += 1
                saveCube(cube,data_file,file_path,i_process+str(i_start),True)
            
            # for resume
            open(os.path.join(outputfolder, i_process+"Nrun.txt"),'w').write(str(Nrun))
        
            if (Nrun+1-i_start_)%n_print == 0:
                printPoint(Nrun+1-i_start_, cube, n_dims, inpar, fixedpar, outpar, lnlike, i_start, i_process)
            
    if num_processes == 1:
        per_run("",int(Naccept[0]),ntotal)
        return
        
    # Create subprocesses
    processes = []
    i_end = 0
    for ii in range(num_processes):
        i_start = int(Naccept[ii] + i_end)
        i_end += af.divide_jobs(ntotal, num_processes, i)
        af.Info('p%i process: %i >>>> %i '%(ii, i_start, i_end))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(ii),i_start,i_end))
        processes.append(p)
    
    run_processes(processes)

        
def randomrun(LnLike, Prior, n_params, inpar, fixedpar, outpar, n_live_points, n_print, outputfolder, num_processes):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a')
    file_path = os.path.join(outputfolder,"SavedFile")
   
    n_dims = len(inpar)
    
    # Initialise cube
    cube = [af.NaN] * n_params

    # For resume 
    # If it is a new scan, Naccept == 1, otherwise Naccept >= 2
    Naccept = getFilelength(os.path.join(outputfolder, af.ResultFile))
    if Naccept >= n_live_points:
      af.ErrorStop('There are already %s living samples in the data file.'%Naccept)
    elif Naccept == 0 :
      af.ErrorStop('The data file is empty. Please start a new scan instead of using the resume mode.')
    Naccept -= 1 
    
    af.Info('Begin random scan ...')
    
    def per_run(i_process, i_accept, i_tot):
        for Nrun in range(i_accept, i_tot) :
            for j in range(n_dims):
                cube[j] = random()
            
            Prior(cube, n_dims, n_params)
            lnlike = LnLike(cube, n_dims, n_params, i_process)
        
            if lnlike > af.log_zero:
                i_accept += 1
                saveCube(cube, data_file, file_path, i_process+str(i_accept), True)
        
            if (Nrun+1)%n_print == 0: 
                printPoint(Nrun+1, cube, n_dims, inpar, fixedpar, outpar, lnlike, i_accept, i_process)
    
    if num_processes == 1:
        per_run("",Naccept,n_live_points)
        return
        
    # Create subprocesses
    processes = []
    for i in range(num_processes):
        i_accept = af.divide_jobs(Naccept, num_processes, i)
        i_tot = af.divide_jobs(n_live_points, num_processes, i)
        af.Info('p%i process: %i >>>> %i '%(i, i_accept, i_tot))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(i),i_accept,i_tot))
        processes.append(p)
    
    run_processes(processes)

def mcmcrun(LnLike, Prior, n_params, n_live_points, inpar, fixedpar, outpar, StepSize, AccepRate, FlagTuneR, InitVal, n_print, outputfolder, num_processes):
    data_file = open(os.path.join(outputfolder, af.ResultFile),'a')
    all_data_file = open(os.path.join(outputfolder, af.ResultFile_MCMC),'a')
    file_path = os.path.join(outputfolder,"SavedFile")
        
    n_dims = len(inpar)
        
    # Initialise cube
    cube = [af.NaN] * n_params

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
        if lnlike > af.log_zero * 0.1 : break
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
    CurObs.append(0) # dwell
    printPoint(0, cube, n_dims, inpar, fixedpar, outpar, lnlike, 0, '')

    if af.resume:
      if lnlike < af.log_zero * 10. :
            af.ErrorStop('Cannot use resume mode because the last point in the previous scan is unphysical.')
      Nrun= getFilelength(os.path.join(outputfolder, af.ResultFile_MCMC))
      # If Nrun <= 2, it will already give error in scan_controller
      with open(os.path.join(outputfolder, af.ResultFile), 'r+') as f:
        # remove the last point from previous scan
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        f.writelines(lines[:-1])
      Naccept = getFilelength(os.path.join(outputfolder, af.ResultFile))
      if Naccept >= n_live_points:
        af.ErrorStop('There are already %s living samples in the data file.'%Naccept)
    else:
      Nrun= 1
      Naccept = 1      
      
    
    def per_run(i_process, i_accept, i_run, i_tot, CurChisq):
    
      # Initialize the MCMC parameters
      MinChisq = CurChisq
      Chisq = CurChisq
      Nout=1
      dwell = 1
      kcovar = 0
      while i_accept < i_tot:

        RangeFlag = True
        for j in range(n_dims):
            par[j] = gauss(CurPar[j],exp(kcovar)*covar[j]) # normalized to 1 
            #rd = random()
            #par[j] = CurPar[j] + covar[j] * (0.5-rd)*2
        if max(par)>1 or min(par)<0 :
            RangeFlag = False
            Nout = Nout +1
            if Nout%100 == 0: 
                af.WarningNoWait("Too many points out of range! Check the range or start point.")
        else:
            i_run += 1
            Nout=0
            for i in range(n_dims): cube[i] = par[i]
            Prior(cube, n_dims, n_params)
            lnlike = LnLike(cube, n_dims, n_params,i_process)
            AllOutMCMC = cube.copy()
            AllOutMCMC.append(1)
            saveCube(AllOutMCMC, all_data_file, file_path, i_process+'0', False)
            Chisq = - 2.0 * lnlike

        Flag_accept = RangeFlag and (Chisq < CurChisq + 20) 
        if Flag_accept: 
            if CurChisq > Chisq: 
                Flag_accept = True
            else:
                Flag_accept = random() < exp(CurChisq-Chisq) 
        if Flag_accept :
            CurObs[-1]=dwell
            #"Naccept+1" due to file of Chisq have covered file of CurChisq
            saveCube(CurObs, data_file, file_path, i_process+str(i_accept+1), True)
            CurChisq = Chisq
            for i in range(n_params): CurObs[i]   = cube[i]
            for i in range(n_dims):   CurPar[i]   = par[i]
            
            if Chisq < MinChisq : MinChisq = Chisq
            i_accept += 1
            dwell = 1
        else:
            if RangeFlag:
                dwell +=1

        AccRat = float(i_accept)/float(i_run)

        if FlagTuneR and i_run < 1000 and i_run > 10 :
            kcovar = kcovar + 1.0/(float(i_run)**0.7)*(AccRat - AccepRate)
        else: 
            kcovar = 1

        if i_run%n_print == 0:
            if RangeFlag:
                printPoint(i_run, cube, n_dims, inpar, fixedpar, outpar, lnlike, i_accept-1, i_process)
                printPoint4MCMC(Chisq,CurChisq,MinChisq,AccRat,FlagTuneR,kcovar)

      # save the last point
      CurObs[-1]=dwell
      saveCube(CurObs, data_file, file_path, i_process+str(i_accept), True)

    if num_processes == 1:
        per_run("", Naccept, Nrun, n_live_points, CurChisq)
        return
        
    # Create subprocesses
    processes = []
    for i in range(num_processes):
        i_tot = af.divide_jobs(n_live_points, num_processes, i)
        i_accept = max(af.divide_jobs(Naccept, num_processes, i),1)
        i_run = max(af.divide_jobs(Nrun, num_processes, i),1)
        af.Info('p%i process has %i live points .............'%(i, i_tot))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(i), i_accept, i_run, i_tot, CurChisq))
        processes.append(p)

    run_processes(processes)


def multinestrun(LnLike, Prior, n_dims, n_params, seed, outputfiles_basename, n_live_points, verbose, resume, importance_nested_sampling, num_processes):
    import pymultinest
    # See https://johannesbuchner.github.io/PyMultiNest/_modules/pymultinest/run.html
    # for more settings
    import functools
    
    def get_i_folder(i_process):
        return outputfiles_basename[:-14] + i_process + "MultiNestData/"
    
    # TODO: check whether we can set random seed here
    def per_run(i_process, i_live_points):
        i_outputfiles_basename = get_i_folder(i_process)
        i_LnLike = functools.partial(LnLike, i_process=i_process)
        pymultinest.run(
            LogLikelihood        = i_LnLike,
            Prior                = Prior,
            n_dims               = n_dims,
            n_params             = n_params,
            seed                 = seed,
            outputfiles_basename = i_outputfiles_basename,
            n_live_points        = i_live_points,
            verbose                    = verbose,
            resume                     = resume,
            importance_nested_sampling = importance_nested_sampling)
        
    if num_processes == 1:
        per_run("",n_live_points)
        return
    
    if resume:
        for ii in range(num_processes):
            if not os.path.exists(get_i_folder("p%s_"%str(ii))):
                af.ErrorStop('Can not use resume mode because of no '+ "p%s_"%str(ii) + "MultiNestData/")
    else:
        for ii in range(num_processes):
            os.mkdir(get_i_folder("p%s_"%str(ii)))
    
    # Create subprocesses
    processes = []

    for i in range(num_processes):
        i_live_points = af.divide_jobs(n_live_points, num_processes, i)
        af.Info('p%i process has %i live points .............'%(i, i_live_points))
        p = multiprocessing.Process(target = per_run, args=("p%s_"%str(i), i_live_points))
        processes.append(p)

    run_processes(processes)
    
    # Combine results
    
    with open(get_i_folder('')+'.txt', 'w') as merged_file:
        af.Info('Merge result.')
        for ii in range(num_processes):
            with open(get_i_folder("p%s_"%str(ii))+'.txt', 'r') as file:
                filetext = file.read()
                merged_file.write(filetext)
                if not filetext.endswith('\n'):
                    merged_file.write('\n')
    

