####################################################################
#    Funtions used for differeny scan methods.                     #
####################################################################

## External modules.
import os,sys,shutil
from random import random,gauss
from math import exp
## Internal modules.
import init as sf
from ploter import PLOTER

def saveCube(cube,f,path,num,save):
    for ii in cube:
        try:
            float(ii)
            f.write(str(ii)+'\t')
        except:
            if os.path.exists(ii):
                if save: shutil.copy(ii, os.path.join(path, os.path.basename(ii)+"."+num))
            else:
                f.write(str(ii)+'\t')
                
    f.write('\n')
    f.flush()


def readrun(LogLikelihood,Prior,n_dims,n_params,inpar,outpar,bin_num,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    f_path = os.path.join(outputfiles_basename,"SavedFile")

    Ploter = PLOTER()
    Ploter.setPlotPar(outputfiles_basename)
    
    for i,name in enumerate(inpar):
        try:
            Ploter._data[name]
        except:
            sf.ErrorStop('Input parameter "%s" could not found in ScanInf.txt.'%name)

    for i,name in enumerate(inpar):
        ntotal = len(Ploter._data[name])
        break

    cube = []
    for i in range(n_params): cube.append(0)

    sf.Info('Begin read scan ...')
    
    Naccept = 0
    
    for Nrun in range(ntotal):
        iner = Nrun
        for i,name in enumerate(inpar):
            cube[i] = Ploter._data[name][iner]
        
        loglike = LogLikelihood(cube, n_dims, n_params)
        if loglike > sf.log_zero:
            Naccept += 1
            saveCube(cube,f_out,f_path,str(Naccept),True)
        saveCube(cube,f_out2,f_path,str(Naccept),False)
        
        if Nrun%n_print == 0:
            sf.Info('------------ Num: %i ------------'%(Nrun+1))
            for i,name in enumerate(inpar):
                sf.Info('Input  - %s = %s '%(name,cube[i]))
            sf.Info('.................................')
            for i,name in enumerate(outpar):
                sf.Info('Output - %s = %s '%(name,cube[i+n_dims]))
            sf.Info('.................................')
            sf.Info('     loglike   = '+str(loglike))
            sf.Info('Accepted Num   = '+str(Naccept))
            sf.Info('Total    Num   = '+str(Nrun+1))


def gridrun(LogLikelihood,Prior,n_dims,n_params,inpar,outpar,bin_num,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #new 20180420 liang
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    f_path = os.path.join(outputfiles_basename,"SavedFile")
    
    ntotal = 1
    cube = []
    interval = {}
    for i,name in enumerate(inpar):
        interval[name] = 1.0 / bin_num[name]
        bin_num[name] += 1
        ntotal     *= bin_num[name]
    for i in range(n_params): cube.append(0)

    sf.Info('Begin grid scan ...')
    
    Naccept = 0

    for Nrun in range(ntotal):
        iner = 1
        for i,name in enumerate(inpar):
            cube[i] = ( int(Nrun/iner) )%bin_num[name] * interval[name]
            iner   *= bin_num[name]
        
        Prior(cube, n_dims, n_params)
        loglike = LogLikelihood(cube, n_dims, n_params)
        if loglike > sf.log_zero:
            Naccept += 1
            saveCube(cube,f_out,f_path,str(Naccept),True)
        saveCube(cube,f_out2,f_path,str(Naccept),False)

        if Nrun%n_print == 0:
            sf.Info('------------ Num: %i ------------'%(Nrun+1))
            for i,name in enumerate(inpar):
                sf.Info('Input  - %s = %s '%(name,cube[i]))
            sf.Info('.................................')
            for i,name in enumerate(outpar):
                sf.Info('Output - %s = %s '%(name,cube[i+n_dims]))
            sf.Info('.................................')
            sf.Info('     loglike   = '+str(loglike))
            sf.Info('Accepted Num   = '+str(Naccept))
            sf.Info('Total    Num   = '+str(Nrun+1))


def randomrun(LogLikelihood,Prior,n_dims,n_params,inpar,outpar,n_live_points,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #new 20180420 liang
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    f_path = os.path.join(outputfiles_basename,"SavedFile")
    
    cube = []
    for i in range(n_params): cube.append(0.0)
    # Initialise the cube

    sf.Info('Begin random scan ...')
    Naccept = 0
    for Nrun in range(n_live_points) :
        for j in range(n_dims):
            cube[j] = random()
        
        Prior(cube, n_dims, n_params)
        loglike = LogLikelihood(cube, n_dims, n_params)
        
        if loglike > sf.log_zero:
            Naccept += 1
            saveCube(cube,f_out,f_path,str(Naccept),True)
        saveCube(cube,f_out2,f_path,str(Naccept),False)
        
        if Nrun%n_print == 0:
            sf.Info('------------ Num: %i ------------'%(Nrun+1))
            for i,name in enumerate(inpar):
                sf.Info('Input  - %s = %s '%(name,cube[i]))
            sf.Info('.................................')
            for i,name in enumerate(outpar):
                sf.Info('Output - %s = %s '%(name,cube[i+n_dims]))
            sf.Info('.................................')
            sf.Info('     loglike   = '+str(loglike))
            sf.Info('Accepted Num   = '+str(Naccept))
            sf.Info('Total    Num   = '+str(Nrun+1))




def mcmcrun(LogLikelihood,Prior,n_dims,n_params,n_live_points,inpar,outpar,StepSize,AccepRate,FalgTune,InitVal,n_print,outputfiles_basename,outputfiles_filename):
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
        saveCube(cube,f_out2,f_path,'0',False)
        if loglike > sf.log_zero / 2.0 : break
        if n_init == 0 : 
            sf.WarningNoWait('The initial point is unphysical, it will find the physical initial points randmly.')
        n_init = n_init +1
        if n_init>100:
            sf.ErrorStop('Can not find physical initial points with 100 tries.')
        for i in range(n_dims):
            cube[i] = random()
            CurPar[i] = cube[i]

    CurObs=[]
    CurChisq = - 2.0 * loglike
    for i in range(n_params): CurObs.append( cube[i] )
    CurObs.append(0) # mult
    sf.Info('------------ Start Point ------------')
    for i,name in enumerate(inpar):
       sf.Info('Input  - %s = %s '%(name,cube[i]))
    for i,name in enumerate(outpar):
       sf.Info('Output - %s = %s '%(name,cube[i+n_dims]))
    sf.Info('.................................')
    sf.Info('Current  Chi^2 = '+str(CurChisq))

    # Initialize the MCMC parameters
    MinChisq = CurChisq
    Chisq = CurChisq
    Nrun= 0
    Naccept = 0
    Nout=0
    mult = 1
    kcovar = 0 
    while Naccept < n_live_points:

        Nrun += 1
        RangeFlag = True
        for j in range(n_dims):
            rd = random()
            par[j] = gauss(CurPar[j],exp(kcovar)*covar[j]) # normalized to 1 
            #par[j] = CurPar[j] + covar[j] * (0.5-rd)*2
        if max(par)>1 or min(par)<0 :
            RangeFlag = False
            Nout = Nout +1
            if Nout%100 == 0: 
                sf.WarningNoWait("Too many points out of range!")
        else:
            Nout=0
            for i in range(n_dims): cube[i] = par[i]
            Prior(cube, n_dims, n_params)
            loglike = LogLikelihood(cube, n_dims, n_params)
            saveCube(cube,f_out2,f_path,'0',False)
            Chisq = - 2.0 * loglike

        Flag_accept = RangeFlag and (Chisq < CurChisq + 20) 
        if Flag_accept: 
            if CurChisq > Chisq: 
                Flag_accept = True
            else:
                Flag_accept = random() < exp(CurChisq-Chisq) 
        if Flag_accept :
            CurObs[-1]=mult
            saveCube(CurObs,f_out,f_path,str(Naccept),True)
            CurChisq = Chisq
            for i in range(n_params): CurObs[i]   = cube[i]
            for i in range(n_dims):   CurPar[i]   = par[i]
            
            if Chisq < MinChisq : MinChisq = Chisq
            Naccept += 1
            mult = 1

        else:
            mult +=1


        AccRat = float(Naccept)/float(Nrun)

        if FalgTune and Nrun < 1000: 
            kcovar = kcovar + 1.0/(float(Nrun)**0.7)*(AccRat - AccepRate)
        else: 
            kcovar = 1

        if Nrun%n_print == 0:
            sf.Info('------------ Num: %i ------------'%Nrun)
            for i,name in enumerate(inpar):
                sf.Info('Input  - %s = %s '%(name,cube[i]))
            if (Chisq < - 2.0 * sf.log_zero) and RangeFlag:
                sf.Info('.................................')
                for i,name in enumerate(outpar):
                    sf.Info('Output - %s = %s '%(name,cube[i+n_dims]))
                sf.Info('.................................')
                sf.Info('Test     Chi^2 = '+str(Chisq))
            sf.Info('Current  Chi^2 = '+str(CurChisq))
            sf.Info('Mimimum  Chi^2 = '+str(MinChisq))
            sf.Info('Accepted Num  = '+str(Naccept))
            sf.Info('Total    Num   = '+str(Nrun))
            sf.Info('Accepted Ratio = '+str(AccRat))
            if FalgTune :
                sf.Info('StepZize factor= '+str(exp(kcovar)))

    # save the last point
    CurObs[-1]=mult
    saveCube(CurObs,f_out,f_path,str(Naccept),True)






