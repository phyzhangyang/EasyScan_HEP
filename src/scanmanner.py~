####################################################################
#    Funtions used for differeny scan methods.
####################################################################

## External modules.
import os,sys,shutil
from random import random,gauss
from math import exp
## Internal modules.
import init as sf

def readrun(LogLikelihood,Prior,n_dims,n_params,inpar,bin_num,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #new 20180420 liang
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')

    import mainfun as mf
    Ploter = mf.plot()
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
            #new 20180420 liang
            cube, acceptFiles=sf.checkFileInList(cube)

            Naccept += 1
            f_out.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
            f_out.flush()

            #new 20180420 liang
            for File in acceptFiles:
                path = os.path.join(outputfiles_basename,"SavedFile")
                SavePath = os.path.join(path, os.path.basename(File)+"."+str(Naccept))
                shutil.copy(File, SavePath)
        
        if (Nrun+1)%n_print == 0:
            print '------------ Num: %i ------------'%(Nrun+1)
            print 'Input    par   = '+str(cube[0:n_dims])
            print 'Output   par   = '+str(cube[n_dims:n_params])
            print '     loglike   = '+str(loglike)
            print 'Accepted Num   = '+str(Naccept)
            print 'Total    Num   = '+str(Nrun+1)

        f_out2.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
        f_out2.flush()

def gridrun(LogLikelihood,Prior,n_dims,n_params,inpar,bin_num,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #new 20180420 liang
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    
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

    for Nrun in xrange(ntotal):
        iner = 1
        for i,name in enumerate(inpar):
            cube[i] = ( int(Nrun/iner) )%bin_num[name] * interval[name]
            iner   *= bin_num[name]
        
        Prior(cube, n_dims, n_params)
        loglike = LogLikelihood(cube, n_dims, n_params)
        if loglike > sf.log_zero:

            #new 20180420 liang
            cube, acceptFiles=sf.checkFileInList(cube)

            Naccept += 1
            f_out.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
            f_out.flush()

            #new 20180420 liang
            for File in acceptFiles:
                path = os.path.join(outputfiles_basename,"SavedFile")
                SavePath = os.path.join(path, os.path.basename(File)+"."+str(Naccept))
                shutil.copy(File, SavePath)

        if (Nrun+1)%n_print == 0:
            print '------------ Num: %i ------------'%(Nrun+1)
            print 'Input    par   = '+str(cube[0:n_dims])
            print 'Output   par   = '+str(cube[n_dims:n_params])
            print '     loglike   = '+str(loglike)
            print 'Accepted Num   = '+str(Naccept)
            print 'Total    Num   = '+str(Nrun+1)

        f_out2.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
        f_out2.flush()


def randomrun(LogLikelihood,Prior,n_dims,n_params,n_live_points,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    #new 20180420 liang
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')
    
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

            #new 20180420 liang
            cube, acceptFiles=sf.checkFileInList(cube)

            Naccept += 1
            f_out.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
            f_out.flush()

            #new 20180420 liang
            for File in acceptFiles:
                path = os.path.join(outputfiles_basename,"SavedFile")
                SavePath = os.path.join(path, os.path.basename(File)+"."+str(Naccept))
                shutil.copy(File, SavePath)
 

        if (Nrun+1)%n_print == 0:
            print '------------ Num: %i ------------'%(Nrun+1)
            print 'Input    par   = '+str(cube[0:n_dims])
            print 'Output   par   = '+str(cube[n_dims:n_params])
            print '     loglike   = '+str(loglike)
            print 'Accepted Num   = '+str(Naccept)
            print 'Total    Num   = '+str(Nrun+1)

        f_out2.write('\t'.join([str(x) for x in cube])+'\t'+str(loglike)+'\n')
        f_out2.flush()


def mcmcrun(LogLikelihood,Prior,n_dims,n_params,n_live_points,inpar,outpar,StepSize,AccepRate,FalgTune,InitVal,n_print,outputfiles_basename,outputfiles_filename):
    f_out = open(os.path.join(outputfiles_basename,outputfiles_filename),'w')
    f_out2 = open(os.path.join(outputfiles_basename,'All_'+outputfiles_filename),'w')

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
    print '------------ Start Point ------------'
    for i,name in enumerate(inpar):
       print 'Input  - %s = %s '%(name,cube[i])
    for i,name in enumerate(outpar):
       print 'Output - %s = %s '%(name,cube[i+n_dims])
    print '.................................'
    print 'Current  Chi^2 = '+str(CurChisq)

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
            Chisq = - 2.0 * loglike

        Flag_accept = RangeFlag and (Chisq < CurChisq + 20) 
        if Flag_accept: 
            if CurChisq > Chisq: 
                Flag_accept = True
            else:
                print random()
                raw_input()
                Flag_accept = random() < exp(CurChisq-Chisq) 

        if Flag_accept :
            #new 20180420 liang
            CurObs, acceptFiles=sf.checkFileInList(CurObs)

            f_out.write('\t'.join([str(x) for x in CurObs])+'\t'+str(-2*CurChisq)+'\t'+str(mult)+'\n')
            f_out.flush()
            CurChisq = Chisq
            for i in range(n_params): CurObs[i]   = cube[i]
            for i in range(n_dims):   CurPar[i]   = par[i]
            
            if Chisq < MinChisq : MinChisq = Chisq
            Naccept += 1
            mult = 1

            #new 20180420 liang
            for File in acceptFiles:
                path = os.path.join(outputfiles_basename,"SavedFile")
                SavePath = os.path.join(path, os.path.basename(File)+"."+str(Naccept))
                shutil.copy(File, SavePath)

        else:
            mult +=1

        AccRat = float(Naccept)/float(Nrun)

        if FalgTune and Nrun < 1000: kcovar = kcovar + 1.0/(float(Nrun)**0.7)*(AccRat - AccepRate)
        else: kcovar =1

        if Nrun%n_print == 0:
            print '------------ Num: %i ------------'%Nrun
            for i,name in enumerate(inpar):
                print 'Input  - %s = %s '%(name,cube[i])
            if (Chisq < - 2.0 * sf.log_zero) and RangeFlag:
                print '.................................'
                for i,name in enumerate(outpar):
                    print 'Output - %s = %s '%(name,cube[i+n_dims])
                print '.................................'
                print 'Test     Chi^2 = '+str(Chisq)
            print 'Current  Chi^2 = '+str(CurChisq)
            print 'Mimimum  Chi^2 = '+str(MinChisq)
            print 'Accepted Num  = '+str(Naccept)
            print 'Total    Num   = '+str(Nrun)
            print 'Accepted Ratio = '+str(AccRat)
            if FalgTune :
                print 'StepZize factor= '+str(exp(kcovar))

        if RangeFlag and (Chisq < - 2.0 * sf.log_zero) :
            f_out2.write('\t'.join([str(x) for x in CurObs])+'\t'+str(-2*CurChisq)+'\t'+str(mult)+'\n')
            f_out2.flush()







