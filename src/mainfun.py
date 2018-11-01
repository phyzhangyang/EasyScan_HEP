####################################################################
#    Read in the input configure file
####################################################################

import os,sys
import re,shutil
import subprocess
import linecache
import random
import numpy
import time
import math

import init as sf

class EasyScanInput:
    def __init__(self):
        self._FileName   = 'test'
        self._PointNum   = 100
        self._ScanMethod = 'random'
        self._RandomSeed = -1
        self._PrintNum   = 10
        self._AccepRate  = 0.25
        self._FalgTuneR  = False

        self.InputPar = {}
        
        self._Prog    = {}
        self.AllPar   = {}
        self.InPar    = {}
        self.OutPar   = {}
        
        self.GridBin = {}
        self.MCMCsw = {}  # Step width
        self.MCMCiv = {}  # Initial value

        self.MNOutputFile = 'test/MultiNestData/'
    
        self._Count   = 0
    
    def setFileName(self, name):
        # turn the result file path into absolute path
        if name.startswith('/home') or name.startswith('~'):
            self._FileName = name
        else:
            self._FileName = os.path.join(sf.CurrentPath, name)
        if self._ScanMethod not in ['READ','PLOT']:
            # deal with the satuation that the result file path already exists.
            if os.path.exists(self._FileName):
                sf.Info(("* The Result file [%s] already exists." % name ))
                while True:
                    c = input("Choose: (r)replace, (b)backup, (s)stop\n")
                    if c == "r":
                         os.system(r"rm -r %s" %self._FileName)
                         break
                    elif c == "b":
                        if not (os.path.exists(sf.CurrentPath+"/Backup")):
                            os.mkdir(sf.CurrentPath+"/Backup")
                        BackupTime = time.strftime("_%Y_%m_%d_%H_%M_%S", time.localtime())
                        BackupPath = os.path.join(sf.CurrentPath, 'Backup/'+name+BackupTime)
                        os.system(r"mv %s %s" %(self._FileName, BackupPath))
                        break
                    elif c == "s":
                        exit(1)
                    else:
                        sf.Info("Wrong input! Please type in one of ('r','b','s')")
            # generate the result file path
            os.mkdir(self._FileName)
            os.mkdir(os.path.join(self._FileName,'SavedFile'))
            if self._ScanMethod == 'MULTINEST':
                self.MNOutputFile = os.path.join(self._FileName, "MultiNestData/")
                os.mkdir(self.MNOutputFile)
            sf.Info('...............................................')
            sf.Info('Result file name  = %s'%self._FileName)
        else:
            if not os.path.exists(self._FileName):
                sf.ErrorStop("The result file %s doesn't exist."%self._FileName)
            if not os.path.exists(os.path.join(self._FileName,'ScanInf.txt')):
                sf.ErrorStop("The info filr %s/ScanInf.txt doesn't exist. Please Check the input."%self._FileName)
                    
            if self._ScanMethod == 'READ':
                sf.Info("* Now you are in READ mode. About files in \n* %s\n*, do you want to  would be delete." % os.path.join(self._FileName, "SavedFile") )
                while True:
                    c = input("Choose: (d)delete, (b)backup, (s)stop\n")
                    if c == "s":
                        exit(1)
                    elif c == "d":
                        os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'SavedFile'))
                        break
                    elif c == "d":
                        sf.Info("Not ready, please choose another way.")
                        exit(1)
                        if not (os.path.exists(sf.CurrentPath+"/Backup")):
                            os.mkdir(sf.CurrentPath+"/Backup")
                        BackupTime = time.strftime("_%Y_%m_%d_%H_%M_%S", time.localtime())
                        BackupPath = os.path.join(sf.CurrentPath, 'Backup/'+name+BackupTime)
                        os.system(r"mv %s %s" %(self._FileName, BackupPath))
                        break
                    else:
                        sf.Info("Wrong input! Please type in one of ('c', 's')")
            sf.Info('...............................................')
            sf.Info('Result file name  = %s'%self._FileName)

    def setPointNum(self, ntot):
        if int(ntot) < 1 :
            sf.ErrorStop('"Points number" should larger than 0')
        self._PointNum = int(ntot)
        sf.Info('Points number     = %s'%self._PointNum)

    def setScanMethod(self, method):
        if method.upper() not in ['RANDOM', 'MCMC', 'MULTINEST','GRID','READ','PLOT','ONEPOINT']:
            sf.ErrorStop('%s is not in supported methods'%method)
        self._ScanMethod = method.upper()
        sf.Info('Scan method       = %s'%self._ScanMethod)

    def setRandomSeed(self, iseed):
        self._RandomSeed = int(iseed)
        ## If iseed is provided in the input file, initialize the basic random number generator
        ## Otherwise, it will be initialized by current system time, and self._RandomSeed = -1,
        ## which means also initialized by current system time in MultiNest
        random.seed( self._RandomSeed )
        sf.Info('Random seed       = %s'%self._RandomSeed)

    def setAccepRate(self, AccepRate):
        AccepRate = float(AccepRate)
        if AccepRate >= 1 or AccepRate <= 0:
            sf.ErrorStop('The acceptance rate must be in [0,1]. The suggest value is 0.5 for d<=2, 0.25 otherwise.')
        self._AccepRate = AccepRate
        self._FalgTuneR = True
        sf.Info('Acceptance rate   = %s'%self._AccepRate)

    def setPrintNum(self, nprint):
        if int(nprint) < 1 :
            sf.ErrorStop('Nprint should larger than 0')
        self._PrintNum = int(nprint)
        sf.Info('Interval of print = %s'%self._PrintNum)

    def setInputPar(self, inputvar):
        inputvar = sf.string2nestlist(inputvar)
        
#        if self._ScanMethod in ['PLOT','READ']:
#            if len(inputvar) == 1 :
#                if inputvar[0][0].startswith('/home') or inputvar[0][0].startswith('~'):
#                    inputfolder = inputvar[0][0]
#                else:
#                    inputfolder = os.path.join(sf.CurrentPath, inputvar[0][0])
#                if not os.path.exists(inputfolder):
#                    sf.ErrorStop('Input folder "%s" do not exist.'%inputfolder)
#            else:
#                sf.ErrorStop( 'For the scan method you choosed, [Plot,READ], only need one parameter.' )
#            sf.Info('Input parameters   = %s'%inputfolder)

        ## inputvar is list of list of input parameters define in section [scan]
        sf.Info('Input parameters   =  ')
        for ii in inputvar:
            self.AllPar[ii[0]]=sf.NaN
            self.InPar[ii[0]]=sf.NaN
            self.InputPar[ii[0]] = ii
            lenii = len(ii)
            if self._ScanMethod in ['RANDOM', 'MULTINEST','GRID','MCMC']:
                if lenii < 4 :
                    sf.ErrorStop( 'For the scan method you choosed, the items of each input parameter should include at least 4 items ( ID, Prior distribution, Minimum, Maximum ). But the paramter [%s] missed %i of them.'%(ii[0],4-lenii) )
                elif lenii > 4 and (self._ScanMethod in ['RANDOM', 'MULTINEST']):
                    sf.WarningWait( 'For the scan method you choosed, only 4 items ( ID, Prior distribution, Minimum, Maximum ) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,lenii-4) )
                    sf.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue
                elif self._ScanMethod == 'GRID':
                    if   lenii == 4 :
                        self.GridBin[ii[0]]=20
                        sf.WarningWait('As the scan method "Grid", the bins number of the parameter [%s] is not provided, which will be set to default value, %i, in this run.'%(ii[0],self.GridBin[ii[0]]) )
                    elif lenii == 5:
                        if ii[4]<=0.0:
                            self.GridBin[ii[0]]=20
                            sf.WarningWait('As the scan method "Grid", the bins number of the parameter [%s] should be positive, which will be set to default value, %i, in this run.'%(ii[0],self.GridBin[ii[0]]) )
                            continue
                        if type(ii[4]) != int and (not ii[4].is_integer()):
                            sf.WarningNoWait('As the scan method "Grid", the bins number of the parameter [%s] shoule be integer, which will be set to int(%f)=%i in this run.'%(ii[0],ii[4],int(ii[4])) )
                        self.GridBin[ii[0]]=int(ii[4])
                    else :
                        sf.WarningWait('For the scan method you choosed, only 5 items ( ID, Prior distribution, Minimum, Maximum, Bins number ) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,lenii-5) )
                    sf.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue
                elif self._ScanMethod == 'MCMC':
                    if   lenii < 6 :
                        sf.WarningWait('As the scan method "MCMC", for the parameter [%s], the step width and initial value (or only initial value) is not provided. \nBoth of them will be set to default values, "step width = (Max-Min)/30" and "initial value = (Max-Min)/2", in this run.'%ii[0] )
                        self.MCMCsw[ii[0]] = 1.0/30.0
                        self.MCMCiv[ii[0]] = 1.0/2.0
                    elif lenii == 6:
                        self.MCMCsw[ii[0]] = 1.0/float(ii[4])      
                        if ii[1].lower() == 'flat':
                            self.MCMCiv[ii[0]] = float(ii[5]-ii[2])/float(ii[3]-ii[2])
                        elif ii[1].lower() == 'log':
                            self.MCMCiv[ii[0]] = (math.log10(ii[5])-math.log10(ii[2]))/(math.log10(ii[3]) - math.log10(ii[2]))
                    else :
                        sf.WarningWait('For the scan method you choosed, only 6 items ( ID, Prior distribution, Minimum, Maximum, Step width, Initial value) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,5-lenii) )
                    sf.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue

            if self._ScanMethod == 'OnePoint':
                if lenii < 2 :
                    sf.ErrorStop('For the scan method you choosed, the items of each input parameter should include at least 2 items ( ID, Value). But the paramter [%s] missed %i of them.'%(ii[0],2-lenii))
                elif lenii > 2:
                    sf.WarningWait('For the scan method you choosed, only 2 items ( ID, Value) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],linii,2-lenii))
                continue

    # not used for now
    def setProgID(self,progID):
        self._ProgID = progID
    
    def setProgram(self,prog):
        self._Prog = prog
        ## copy input vars into allvars
        for ii in prog:
            sf.Debug('Programe ID', ii)
            sf.Debug('Corresponding output vars', prog[ii].outvar)

            ## deep copy
            ## prog[ii] is different program
            ## and prog[ii].outvar is dictionary of name and value of output variables in each program
            for jj in prog[ii].outvar:
                ## jj is just key of each item in dictionary
                self.AllPar[jj] = prog[ii].outvar[jj]
                self.OutPar[jj] = prog[ii].outvar[jj]

            ## add "math .." in "Input variable"/"Bound"/"gaussian" to self.AllPar to avoid duplication
            for jj in prog[ii].invar:
                if jj not in list(self.AllPar.keys()):
                    self.AllPar[jj] = prog[ii].invar[jj]
                    self.OutPar[jj] = prog[ii].invar[jj]
            for jj in prog[ii].boundvar:
                if jj not in list(self.AllPar.keys()):
                    self.AllPar[jj] = prog[ii].boundvar[jj]
                    self.OutPar[jj] = prog[ii].boundvar[jj]
            for jj in prog[ii].cgauvar:
                if jj not in list(self.AllPar.keys()):
                    self.AllPar[jj] = prog[ii].cgauvar[jj]
                    self.OutPar[jj] = prog[ii].cgauvar[jj]
            for jj in prog[ii].cffchi2var:
                if jj not in list(self.AllPar.keys()):
                    self.AllPar[jj] = prog[ii].cffchi2var[jj]
                    self.OutPar[jj] = prog[ii].cffchi2var[jj]

        ## new 20180428 liang thanks to yuanfang
        self.OutPar = sf.sortDic(self.OutPar)

        sf.Debug('All vars:   ',self.AllPar)
        sf.Debug('Input Pars: ',self.InPar)
        sf.Debug('Output Pars:',self.OutPar)
                
    def getFileName(self):
        return self._FileName
    def getPointNum(self):
        return self._PointNum
    def getScanMethod(self):
        return self._ScanMethod
    def getRandomSeed(self):
        return self._RandomSeed
    def getPrintNum(self):
        return self._PrintNum
    def getDebugFlag(self):
        return self._DebugFlag
    
    def getStepSize(self):
        return self.MCMCsw
    def getInitialValue(self):
        return self.MCMCiv

    def getFalgTuneR(self):
        return self._FalgTuneR
    def getAccepRate(self):
        return self._AccepRate
    

############################################################
#################   Constraint    Class   ##################
############################################################
class constraint:
    def __init__(self):
        self._Gaussian=[]
        #self._Limit=[]
        self._FreeFormChi2=[]
        self.Chi2={}
    
    def setGaussian(self,var):
        var = sf.string2nestlist(var)
        self.Chi2['Chi2'] = sf.NaN 
        sf.Info('Gaussian Constraint:')
        for ii in var:
            if len(ii) in [3]:
                jj = ii+['symm','Gaussian_%s'%ii[0]]
                self._Gaussian.append(jj)
                self.Chi2[jj[4]] = sf.NaN 
                
                sf.Info('    varID= %s\tMean= %e\tDeviation= %e\tType= %s\tName= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
            elif len(ii) in [4,5]:
                if not ii[3] in ['symm','lower','upper']:
                    sf.ErrorStop( 'For the "Gaussian" constraint on "%s", the "Type" can only be "symm", "upper" or "lower", not "%s".'%(ii[0],ii[3]) )
                ## new 20180428 liang
                if len(ii) == 4: 
                    jj = ii+['Gaussian_%s'%ii[0]]
                else:
                    jj = ii
                self._Gaussian.append(jj)
                self.Chi2[jj[4]] = sf.NaN 
                
                sf.Info('    varID= %s\tMean= %e\tDeviation= %e\tType= %s\tName= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))

            else:
                sf.ErrorStop( 'The "Gaussian" constraint on "%s" need 4 or 5 items( ID, Mean, Deviation, Type [, Name] ).'%(ii[0]) )

        ## new 20180428 liang
        self.Chi2 = sf.sortDic(self.Chi2) 

## marked 20180430 liang
#    def setLimit(self,var):
#        var = sf.string2nestlist(var)
#        sf.Info('Upper/Lower limit:')
#        for ii in var:
#            if len(ii) == 4:
#                self._Limit.append(ii)
#                sf.Info('  varID(X)= %s\tvarID(Y)= %s\tConstraintFile= %s\tType= %s'%(ii[0],ii[1],ii[2],ii[3]))
#        ## add check the ConstraintFile exist
#        ## add check the ConstraintFile has two columns and >1 lines
#        ## very useful, you can simply let a>b
#            else:
#                sf.ErrorStop( 'The "Limit" constraint on "(%s,%s)" need 4 items( X ID, Y ID, ConstraintFile, Type ).'%(ii[0],ii[1]) )

    ## new 20180430 liang
    def setFreeFormChi2(self,var):
        var = sf.string2nestlist(var)
        sf.Info('FreeFormChi2:')
        for ii in var:
            if len(ii) in [1,2]:
                if len(ii) == 1:
                    jj = ii + ['FreeFormChi2_%s'%ii[0]]
                else:
                    jj = ii
                self._FreeFormChi2.append(jj)
                sf.Info('    varID= %s\tName= %s'%(jj[0], jj[1]))

                self.Chi2[jj[1]] = sf.NaN
            else:
                sf.ErrorStop( 'The "FreeFormChi2" constraint on "%s" need 1 item or 2 items( VarID [, Name] ).'%(ii[0]) )

            self.Chi2 = sf.sortDic(self.Chi2)

    def getChisq(self,par):
        chisq = 0.0

        ## add for "math ..." in [constrain]
        sf.parseMath(par)

        for ii in self._Gaussian:
            if   ii[3] == 'symm':
                ichisq =     (ii[1] - par[ii[0]])**2/ii[2]**2
            elif ii[3] == 'upper':
                if ii[1] > par[ii[0]]:
                    ichisq = (ii[1] - par[ii[0]])**2/ii[2]**2
            elif ii[3] == 'lower':
                if ii[1] < par[ii[0]]:
                    ichisq = (ii[1] - par[ii[0]])**2/ii[2]**2

            ## new 20180428 liang
            self.Chi2[ii[4]]=ichisq

            chisq += ichisq

        ## marked 20180430 liang
        #for ii in self._Limit:
        #    sf.ErrorStop('Line limit constraint is not ready.')
        for ii in self._FreeFormChi2:
            ichisq = par[ii[0]]

            self.Chi2[ii[1]]=ichisq

            chisq += ichisq

        ## new 20180428 liang
        self.Chi2['Chi2'] = chisq

        return chisq

############################################################
#################   Plot    Class   ########################
############################################################
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata

histconf={'bins':50, 'normed':False, 'facecolor':'green', 'alpha':0.7}
scatterconf={'s':50, 'marker':'o', 'edgecolors':'None', 'alpha':0.8}
colorconf={'s':40, 'edgecolors':'None',
             'cmap':plt.get_cmap('Reds')}
figconf={'figsize':(12,7), 'dpi':80}
labelconf={'fontsize':20}
legendconf={'fontsize':20}


class plot():
    def __init__(self):
        self._Histogram=[]
        self._Scatter=[]
        self._Color=[]
        self._Contour=[]
    
        self._data = {}
        self._dataAllTry = {}
        self._path = {}

        self._FigNames = []

    def setHistogram(self, hist):
        hist=sf.string2nestlist(hist)
        jj=0
        for ii in hist:
            if len(ii)==1 :
                self._Histogram.append([ii[0],'Histogram_%s.png'%ii[0]])
            elif len(ii)==2 :
                self._Histogram.append([ii[0],'%s.png'%ii[1]])

            self._FigNames.append(self._Histogram[jj][1])
            jj+=1
  
    def setScatter(self, scatter):
        scatter = sf.string2nestlist(scatter)
        jj=0
        for ii in scatter:
            if len(ii)==2 :
                self._Scatter.append([ii[0],ii[1],'Scatter_%s_%s.png'%(ii[0],ii[1])])
            elif len(ii)==3 :
                self._Scatter.append([ii[0],ii[1],'%s.png'%ii[2]])

            self._FigNames.append(self._Scatter[jj][2])
            jj+=1

    def setColor(self, color):
        color = sf.string2nestlist(color)
        jj=0
        for ii in color:
            if len(ii)==3 :
                self._Color.append([ii[0],ii[1],ii[2],'Color_%s_%s_%s.png'%(ii[0],ii[1],ii[2])])
            elif len(ii)==4 :
                self._Color.append([ii[0],ii[1],ii[2],'%s.png'%ii[3]])

            self._FigNames.append(self._Color[jj][3])
            jj+=1

    def setContour(self, Contour):
        Contour = sf.string2nestlist(Contour)
        jj=0
        for ii in Contour:
            if len(ii)==3 :
                self._Contour.append([ii[0],ii[1],ii[2],'Contour_%s_%s_%s.png'%(ii[0],ii[1],ii[2])])
            elif len(ii)==4 :
                self._Contour.append([ii[0],ii[1],ii[2],'%s.png'%ii[3]])

            self._FigNames.append(self._Contour[jj][3])
            jj+=1


    def setPlotPar(self,path):
        ## try this
        # self._data =  np.loadtxt(path)
        
        f_data = open(os.path.join(path,'ScanInf.txt'),'r')
        path   = list(map(str,f_data.readline().split()))
        var    = {}
        while True:
            plot_line = f_data.readline()
            if not plot_line :
                break
            plot_line = list(map(str,plot_line.split()))
            var["+".join(plot_line[:-1])] = int(plot_line[-1])
    
        self._path = os.path.join(path[0],'Figures')
        if not os.path.exists(self._path):
            os.mkdir(self._path)
        else:
            __import__("shutil").rmtree(self._path) 
            os.mkdir(self._path)

        for ii in var:
            self._data[ii] = []
        
        f_data = open(os.path.join(path[0],path[1]),'r')
        while True:
            line = f_data.readline()
            if not line :
                break
            line_par = list(map(str,line.split()))
            for ii in var:
                try:
                    self._data[ii].append(float( line_par[var[ii]] ))
                except:
                    sf.Debug('Skip parameter %s'%ii)

        ##new 20180418 liang
        for ii in var:
            self._dataAllTry[ii] = []

        f_dataAllTry = open(os.path.join(path[0],'All_%s'%path[1]), 'r')
        while True:
            line = f_dataAllTry.readline()
            if not line :
                break
            line_par = list(map(str,line.split()))
            for ii in var:
                try:
                    self._dataAllTry[ii].append(float( line_par[var[ii]] ))
                except:
                    sf.Debug('Skip parameter %s'%ii)

    def checkPar(self,par,num):                
            for jj in range(num):
                if max(self._data[par[jj]]) == min(self._data[par[jj]]):
                    sf.ErrorStop("The parameter %s=%f is a cosntant with all samples, can not creat plot for it. Please correct in [plot] in your input file!"%(par[jj], min(self._data[par[jj]]) )  )
                #new 20180416 liang
                if len(self._data[par[jj]]) == 1:
                    sf.ErrorStop("One sample (e.g., see parameter %s) only, can not creat plot for it. Please correct in [plot] in your input file!"%par[jj])
                    return False 
                try:
                    list(map(float, self._data[par[jj]]))
                except ValueError:
                    sf.WarningNoWait("The parameter %s not a number, can not creat plot for it."%par[jj])
                    return False 
            return True

    ##only debug
    def get_contour_verts(self,cn):
        contours = []
        # for each contour line
        for cc in cn.collections:
            paths = []
            # for each separate section of the contour line
            for pp in cc.get_paths():
                xy = []
                # for each segment of that section
                for vv in pp.iter_segments():
                    xy.append(vv[0])
                paths.append(np.vstack(xy))
            contours.append(paths)
    
        return contours

    def getPlot(self):
        if len(self._Histogram) + len(self._Scatter) + len(self._Color) + len(self._Contour) == 0:
            sf.Info('You have close ploting the result ... ')
            return
        sf.Info('Start to plot the result ... ')
        FigNames=[]
        for ii in self._Histogram :
            sf.Info('Generate histogram plot of parameter %s'%ii[0])
            if not self.checkPar(ii,1): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            subplot.hist(self._data[ii[0]], **histconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel('Count', **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[1]))

        for ii in self._Scatter :
            sf.Info('Generate scatter plot of parameter %s and %s'%(ii[0],ii[1]))
            if not self.checkPar(ii,2): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            subplot.scatter(self._data[ii[0]],self._data[ii[1]],**scatterconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[2]))

            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            subplot.scatter(self._dataAllTry[ii[0]],self._dataAllTry[ii[1]],label='All',**scatterconf)
            subplot.scatter(self._data[ii[0]],self._data[ii[1]],label='Surviving',**scatterconf)
            plt.legend(**legendconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, 'Compare_%s'%ii[2]))

        for ii in self._Color :
            sf.Info('Generate color plot of parameter %s and %s with color %s'%(ii[0],ii[1],ii[2]))
            if not self.checkPar(ii,3): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            sc1=subplot.scatter(self._data[ii[0]],self._data[ii[1]], c= self._data[ii[2]], **colorconf)
            cb1=plt.colorbar(sc1)
            cb1.set_label(ii[2], **labelconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[3]))

        for ii in self._Contour :
            sf.Info('Generate contour plot of parameter %s and %s with contour %s'%(ii[0],ii[1],ii[2]))
            if not self.checkPar(ii,3): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)

            x = self._data[ii[0]]
            X = np.linspace(min(x),max(x),100)
            y = self._data[ii[1]]
            Y = np.linspace(min(y),max(y),100)
            # z = [ np.log10(abs(u)) for u in self._data[ii[2]] ]  # log10
            z = self._data[ii[2]] 
            Z = griddata(x,y,z,X,Y,interp='linear')

            #C = subplot.contour(X,Y,Z, 3,linewidths=2)
            #plt.clabel(C, inline=True, fontsize=8)

            ## debug
            #Cpoint = self.get_contour_verts(C)  
            #np.savetxt(os.path.join(self._path, "contour_1_1.dat"),Cpoint[0][0])  

            C = subplot.contourf(X,Y,Z,3, cmap=plt.cm.rainbow)
            
            cb1=plt.colorbar(C)
            cb1.set_label(ii[2], **labelconf)

            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[3]))


