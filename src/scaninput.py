####################################################################
#    Class SCANINPUT: contral scan                                 #
####################################################################
# External modules
import os
import random
import time
import math
# Internal modules
import auxfun as af

class SCANINPUT:
    def __init__(self):
        self._FileName   = 'test'
        self._PointNum   = 2
        self._ScanMethod = 'random'
        self._RandomSeed = -1
        self._PrintNum   = 10
        self._AccepRate  = 0.25
        self._FlagTuneR  = False

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
            self._FileName = os.path.join(af.CurrentPath, name)
        if self._ScanMethod not in ['READ','PLOT']:
            # deal with the satuation that the result file path already exists.
            if os.path.exists(self._FileName):
                af.Info(("* The Result file [%s] already exists." % name ))
                while True:
                    c = input("Choose: (r)replace, (b)backup, (s)stop\n")
                    if c == "r":
                         os.system(r"rm -r %s" %self._FileName)
                         break
                    elif c == "b":
                        if not (os.path.exists(af.CurrentPath+"/Backup")):
                            os.mkdir(af.CurrentPath+"/Backup")
                        BackupTime = time.strftime("_%Y_%m_%d_%H_%M_%S", time.localtime())
                        BackupPath = os.path.join(af.CurrentPath, 'Backup/'+name+BackupTime)
                        os.system(r"mv %s %s" %(self._FileName, BackupPath))
                        break
                    elif c == "s":
                        exit(1)
                    else:
                        af.Info("Wrong input! Please type in one of ('r','b','s')")
            # generate the result file path
            os.mkdir(self._FileName)
            os.mkdir(os.path.join(self._FileName,'SavedFile'))
            if self._ScanMethod == 'MULTINEST':
                self.MNOutputFile = os.path.join(self._FileName, "MultiNestData/")
                os.mkdir(self.MNOutputFile)
            af.Info('...............................................')
            af.Info('Result file name  = %s'%self._FileName)
        else:
            if not os.path.exists(self._FileName):
                af.ErrorStop("The result file %s doesn't exist."%self._FileName)

            if self._ScanMethod in ['PLOT'] and not os.path.exists(os.path.join(self._FileName,'ScanInf.txt')):
                af.ErrorStop("The information file %s/ScanInf.txt doesn't exist. Please give file ScanInf.txt. Note please properly set the first line in the file."%self._FileName)

            if self._ScanMethod in ['READ'] and not os.path.exists(os.path.join(self._FileName,'ScanInfINPUT.txt')):
                af.ErrorStop("The information file %s/ScanInfINPUT.txt doesn't exist. Please copy ScanInf.txt as ScanInfINPUT.txt and properly set the first line in the file."%self._FileName)
                    
            if self._ScanMethod in ['PLOT', 'READ']:
                af.Info("* Now you are in PLOT (READ) mode:\n* Files of folder Figure (Figure and SavedFile) in\n* %s\n* would be deleted and other files in\n* %s\n* not be deleted.\n* Please choose the way how to deal with files in\n* %s\n*." % (self._FileName, self._FileName, self._FileName) )
                while True:
                    c = input("Choose: (d)elete, (b)ackup, (s)top\n")
                    if c == "s":
                        exit(1)
                    elif c == "d":
                        os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'Figure'))
                        if self._ScanMethod == 'READ':
                            os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'SavedFile'))
                        break
                    elif c == "b":
                        if not (os.path.exists(af.CurrentPath+"/Backup")):
                            os.mkdir(af.CurrentPath+"/Backup")
                        BackupTime = time.strftime("_%Y_%m_%d_%H_%M_%S", time.localtime())
                        BackupPath = os.path.join(af.CurrentPath, 'Backup/'+name+BackupTime)
                        #os.system(r"mv %s %s" %(self._FileName, BackupPath))
                        os.system(r"cp -r %s %s" %(self._FileName, BackupPath))
                        os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'Figure'))
                        if self._ScanMethod == 'READ':
                            os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'SavedFile'))
                        break
                    else:
                        af.Info("Wrong input! Please type in one of ('d', 'b', 's')")
            af.Info('...............................................')
            af.Info('Result file name  = %s'%self._FileName)

    def setPointNum(self, ntot):
        if int(ntot) < 1 :
            af.ErrorStop('"Points number" should larger than 0')
        self._PointNum = int(ntot)
        af.Info('Points number     = %s'%self._PointNum)

    def setScanMethod(self, method):
        if method.upper() not in ['RANDOM', 'MCMC', 'MULTINEST','GRID','READ','PLOT','ONEPOINT']:
            af.ErrorStop('%s is not in supported methods'%method)
        self._ScanMethod = method.upper()
        af.Info('Scan method       = %s'%self._ScanMethod)

    def setRandomSeed(self, iseed):
        self._RandomSeed = int(iseed)
        ## If iseed is provided in the input file, initialize the basic random number generator
        ## Otherwise, it will be initialized by current system time, and self._RandomSeed = -1,
        ## which means also initialized by current system time in MultiNest
        random.seed( self._RandomSeed )
        af.Info('Random seed       = %s'%self._RandomSeed)

    def setAccepRate(self, AccepRate):
        AccepRate = float(AccepRate)
        if AccepRate >= 1 or AccepRate <= 0:
            af.ErrorStop('The acceptance rate must be in [0,1]. The suggest value is 0.5 for d<=2, 0.25 otherwise.')
        self._AccepRate = AccepRate
        self._FalgTuneR = True
        af.Info('Acceptance rate   = %s'%self._AccepRate)

    def setPrintNum(self, nprint):
        if int(nprint) < 1 :
            af.ErrorStop('Nprint should larger than 0')
        self._PrintNum = int(nprint)
        af.Info('Interval of print = %s'%self._PrintNum)

    def setInputPar(self, inputvar):
        inputvar = af.string2nestlist(inputvar)
        
#        if self._ScanMethod in ['PLOT','READ']:
#            if len(inputvar) == 1 :
#                if inputvar[0][0].startswith('/home') or inputvar[0][0].startswith('~'):
#                    inputfolder = inputvar[0][0]
#                else:
#                    inputfolder = os.path.join(af.CurrentPath, inputvar[0][0])
#                if not os.path.exists(inputfolder):
#                    af.ErrorStop('Input folder "%s" do not exist.'%inputfolder)
#            else:
#                af.ErrorStop( 'For the scan method you choosed, [Plot,READ], only need one parameter.' )
#            af.Info('Input parameters   = %s'%inputfolder)

        ## inputvar is list of list of input parameters define in section [scan]
        af.Info('Input parameters   =  ')
        for ii in inputvar:
            self.AllPar[ii[0]]=af.NaN
            self.InPar[ii[0]]=af.NaN
            self.InputPar[ii[0]] = ii
            lenii = len(ii)
            if self._ScanMethod in ['RANDOM', 'MULTINEST','GRID','MCMC']:
                if lenii < 4 :
                    af.ErrorStop( 'For the scan method you choosed, the items of each input parameter should include at least 4 items ( ID, Prior distribution, Minimum, Maximum ). But the paramter [%s] missed %i of them.'%(ii[0],4-lenii) )
                elif lenii > 4 and (self._ScanMethod in ['RANDOM', 'MULTINEST']):
                    af.WarningWait( 'For the scan method you choosed, only 4 items ( ID, Prior distribution, Minimum, Maximum ) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,lenii-4) )
                    af.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue
                elif self._ScanMethod == 'GRID':
                    if   lenii == 4 :
                        self.GridBin[ii[0]]=20
                        af.WarningWait('As the scan method "Grid", the bins number of the parameter [%s] is not provided, which will be set to default value, %i, in this run.'%(ii[0],self.GridBin[ii[0]]) )
                    elif lenii == 5:
                        if ii[4]<=0.0:
                            self.GridBin[ii[0]]=20
                            af.WarningWait('As the scan method "Grid", the bins number of the parameter [%s] should be positive, which will be set to default value, %i, in this run.'%(ii[0],self.GridBin[ii[0]]) )
                            continue
                        if type(ii[4]) != int and (not ii[4].is_integer()):
                            af.WarningNoWait('As the scan method "Grid", the bins number of the parameter [%s] shoule be integer, which will be set to int(%f)=%i in this run.'%(ii[0],ii[4],int(ii[4])) )
                        self.GridBin[ii[0]]=int(ii[4])
                    else :
                        af.WarningWait('For the scan method you choosed, only 5 items ( ID, Prior distribution, Minimum, Maximum, Bins number ) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,lenii-5) )
                    af.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue
                elif self._ScanMethod == 'MCMC':
                    if   lenii < 6 :
                        af.WarningWait('As the scan method "MCMC", for the parameter [%s], the step width and initial value (or only initial value) is not provided. \nBoth of them will be set to default values, "step width = (Max-Min)/30" and "initial value = (Max-Min)/2", in this run.'%ii[0] )
                        self.MCMCsw[ii[0]] = 1.0/30.0
                        self.MCMCiv[ii[0]] = 1.0/2.0
                    elif lenii == 6:
                        self.MCMCsw[ii[0]] = 1.0/float(ii[4])      
                        if ii[1].lower() == 'flat':
                            self.MCMCiv[ii[0]] = float(ii[5]-ii[2])/float(ii[3]-ii[2])
                        elif ii[1].lower() == 'log':
                            self.MCMCiv[ii[0]] = (math.log10(ii[5])-math.log10(ii[2]))/(math.log10(ii[3]) - math.log10(ii[2]))
                    else :
                        af.WarningWait('For the scan method you choosed, only 6 items ( ID, Prior distribution, Minimum, Maximum, Step width, Initial value) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,5-lenii) )
                    af.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))
                    continue

            if self._ScanMethod == 'OnePoint':
                if lenii < 2 :
                    af.ErrorStop('For the scan method you choosed, the items of each input parameter should include at least 2 items ( ID, Value). But the paramter [%s] missed %i of them.'%(ii[0],2-lenii))
                elif lenii > 2:
                    af.WarningWait('For the scan method you choosed, only 2 items ( ID, Value) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],linii,2-lenii))
                continue

    # not used for now
    def setProgID(self,progID):
        self._ProgID = progID
    
    def setProgram(self,prog):
        self._Prog = prog
        ## copy input vars into allvars
        for ii in prog:
            af.Debug('Programe ID', ii)
            af.Debug('Corresponding output vars', prog[ii].outvar)

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
        self.OutPar = af.sortDic(self.OutPar)

        af.Debug('All vars:   ',self.AllPar)
        af.Debug('Input Pars: ',self.InPar)
        af.Debug('Output Pars:',self.OutPar)
                
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

    def getFlagTuneR(self):
        return self._FlagTuneR
    def getAccepRate(self):
        return self._AccepRate
