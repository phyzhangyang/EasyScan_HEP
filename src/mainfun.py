#!/usr/bin/env python
import os,sys
import re,shutil
import subprocess
import linecache
import random
import numpy
import time
import math

import init  as sf

class program:
    def __init__(self):
        self._ProgName  =''
        self._Command   =''
        self._ComPath   =''
        self._InputFile ={}
        self._InputVar  =[]
        self._OutputFile={}
        self._OutputVar ={}

        self._BoundVar  =[]
        
        self._InFileID =[]
        self._OutFileID =[]

        self._InFilVar  ={}
        self._InRepVar  ={}
        self._InPosVar  ={}
        self._InLabVar  ={}
        self._InSLHAVar ={}
        
        self._OutFileVar  ={}
        self._OutPosVar  ={}
        self._OutLabelVar  ={}
        self._OutSLHAVar ={}
    
        self.invar = {}
        self.outvar = {}
        self.boundvar = {}
        self.cgauvar = {}
        self.cffchi2var= {}
            
        self._Count   = 0
        self._runflag = 'True'

    def setProgName(self, name):
        self._ProgName=name
        sf.Info('...............................................')
        sf.Info('...............................................')
        sf.Info('Program name    = %s'% self._ProgName)
    def setCommand(self, command):
        self._Command=sf.string2nestlist(command)
        sf.Info('Execute command = %s'% self._Command)
    def setComPath(self, cpath):
        if cpath.startswith('/home') or cpath.startswith('~'):
            self._ComPath=cpath
        else:
            self._ComPath=os.path.join(sf.CurrentPath, cpath)
        if not os.path.exists(self._ComPath):
            sf.ErrorStop('Command path "%s" do not exist.'%self._ComPath)
        sf.Info('Command path    = %s'% self._ComPath)

    def setInputFile(self, inputfile):
        inputfile=sf.string2nestlist(inputfile)
        self._InFileID = [x[0] for x in inputfile ]
        if self._InFileID != list(set(self._InFileID)):
            sf.ErrorStop('Input file in program "%s" have same File ID.'%self._ProgName)
        sf.Info('Input file      = ')
        for ii in inputfile:
            if len(ii) != 2:
                if ii[0] == '':
                    break
                sf.ErrorStop('The input file of %s need two items (File ID, File path).'%self._ProgName)
            if not (ii[1].startswith('/home') or ii[1].startswith('~')):
                ii[1]=os.path.join(sf.CurrentPath, ii[1])
            self._InputFile[ii[0]]=ii[1]
            sf.Info('  fileID= %s \tFile= %s'%(ii[0],ii[1]))

    ## check functions for whether contents in input variable and output variable in configure is matching with contents in corresponding input file and output file with methods "file", "position", "label", "slha""
    def checkVar_file(self, fileID):
            ## For 'File' method
            ## check the input vars that use 'File' method

            ## file_flag stands for existing the checked input file which is created by user-self and not output by the previous program(s).
            ii = fileID
            file_flag = True

            self._InFilVar[ii]  = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'file')]

            for jj in self._InFilVar[ii]:
                if len(jj) != 4 :
                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "File" method, 4 items (Name, FileID, "File", Method) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['PREVIOUS', 'SAVE', 'REPLACE', 'ADD']:
                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "File" method, the 4th item must be "PREVIOUS", "SAVE", "REPLACE" or "ADD". If you can use other formats, please contact with the authors.'%(jj[0], self._ProgName) )
                if jj[3].upper() == "PREVIOUS":
                    file_flag = False
                sf.Info( 'Becasue file (ID=%s) in program "%s" is obtained from previous program(s), check this input file is ignored.'%(jj[1], self._ProgName))
                sf.Info('  Name= %s \tFileID= %s \t"File"= %s \tMethod %s'%(jj[0],jj[1],jj[2],jj[3]))
           
            return file_flag 

    def checkVar_position(self, fileID, fileFlag, status):
            ## For 'Position' method
            ## check the input vars that use 'Position' method
           
            ii=fileID
            file_flag=fileFlag 

            if status.lower()=="invar":
               self._InPosVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'position')]
            elif status.lower()=="outvar":
               self._OutPosVar[ii] = [x for x in self._OutputVar if (x[1] == ii) and (x[2].lower() == 'position')]
            else:
               return

            for jj in self._InPosVar[ii]:
                if len(jj) != 5 :
                    sf.ErrorStop('For input/output variable "%s" in program "%s" with "Position" method, 5 items ( Name ID, Input file ID, Method, Line number, Column number ) need to be provived.'%(jj[0], self._ProgName) )
                sf.Info('  varID= %s \t fileID= %s \t Method= %s \t Line num= %s \t Column num= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]) )
                
                if file_flag:
                    ## inlines is a list of all lines
                    ## invars is a list of list of words for all lines
                    if status.lower()=="invar":
                        inlines = open(self._InputFile[ii]).readlines()
                    elif status.lower()=="outvar":
                        inlines = open(self._OutputFile[ii]).readlines()
                    else:
                        return
                    invar = [ss.split() for ss in inlines]

                    sf.Debug('Position len(invar)',len(invar))
                    sf.Debug('Position line num',jj[3])
                    if len(invar) < jj[3]:
                        sf.ErrorStop('For input variable "%s" in program "%s" with "Position" method, the line number is larger than the number of lines of inoput file "%s". Please check your input file.'%(jj[0],self._ProgName,self._InputFile[ii]) )
                    sf.Debug('Position len(invar[line num])',len(invar[jj[3]-1]))
                    sf.Debug('Position Column number',jj[4])

                    if len(invar[jj[3]-1]) < jj[4]:
                        sf.ErrorStop('For input variable "%s" in program "%s" with "Position" method, the column number is larger than the number of columns in line %s of input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )


    def checkVar_label(self, fileID, fileFlag):
            ## For 'Label' method
            ## save the input vars that use 'Label' method

            ii=fileID
            file_flag=fileFlag 

            self._InLabVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'label')]
            
            for jj in self._InLabVar[ii]:
                if len(jj) != 5 :
                    sf.ErrorStop('For input variable "%s" in program "%s" with "Label" method, 5 items ( Name ID,  Input file ID,  Method,  Label name,  Input variable column number ) need to be provived.'%(jj[0],self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \tMethod= %s \tLabel= %s \tColumn= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
                
                if file_flag:
                    ## inlines is a list of all lines
                    ## invars is a list of list of words for all lines
                    inlines = open(self._InputFile[ii]).readlines()
                    invar = [ss.split() for ss in inlines]

                    ## enumerate return object that could by use by for loop by parsing list where xxi is id and xx is value in list.
                    labelinum = [xxi for xxi,xx in enumerate(inlines) if jj[3] in xx]
                    if len(labelinum)==0:
                        sf.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, there is no "%s" in the input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                    if len(labelinum)!=1:
                        sf.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, there is more than one line with label "%s" in the input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                    for kk in labelinum:
                        sf.Debug('Labeled line',inlines[kk].strip('\n'))
                        if len(invar[kk]) < jj[4]:
                            sf.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, the column number "%i" defined by user is larger than the number of columns "%i" in the corresponding labeled line "%s" of input file "%s".'%(jj[0], self._ProgName, jj[4], len(invar[labelinum[0]]), invar[labelinum[0]], self._InputFile[ii]) )
                        if invar[kk][jj[4]-1] == jj[3]:
                            sf.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, the extracting item with column ID "%i" is just the label "%s" in the corresponding labeled line "%s" of input file "%s".'%(jj[0], self._ProgName, jj[4], jj[3], invar[labelinum[0]], self._InputFile[ii]) )

    def checkVar_slha(self, fileID, fileFlag):
            ## For 'SLHA' method
            ## check the input vars that use 'SLHA' method

            ii=fileID
            file_flag=fileFlag 

            self._InSLHAVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'slha')]

            for jj in self._InSLHAVar[ii]:
                if len(jj) < 6 :
                    sf.ErrorStop('For input variable "%s" in program "%s" with "SLHA" method, at least 6 items ( Name ID,  Input file ID,  Method, BLOCK/DECAY, Block name/PDG, Keys) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['BLOCK','DECAY']:
                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, the 4th item must be "BLOCK" or "DECAY". If you can use other formats, please contact with the authors.'%(jj[0],self._ProgName) )

                sf.Info('  varID= %s \t fileID= %s \t Method= %s \t Block/Decay= %s \tName= %s \tKeys= %s'%(jj[0], jj[1], jj[2], jj[3], jj[4], jj[5:]))

                if file_flag:
                    ## inlines is a list of all lines
                    ## invars is a list of list of words for all lines
                    inlines = open(self._InputFile[ii]).readlines()
                    invar = [ss.split() for ss in inlines]

                    ## list[i:], begin with i to the end involved; list[i:j], begin with i to the end j but not involved.
                    ## string.split() could get list with seperated string]                    ## jj[4] may be \"a\" or \"a b\"
                    blk = str(jj[4]).split()
                    blk_flag = False
                    #ks  = str(jj[5]).split()
                    ks  = map(str, jj[5:])
                    ks_flag  = False
                    for kk in invar:
                        if kk[0].startswith('#'): continue

                        if blk_flag:
                            ## quick break for loop if no data in block (but not effect on decay, because decay could not have no data)
                            ## quick break for loop if finished iterating over data in block or decay
                            ## if there are two same block or decay, only return data in the first one.
                            if kk[0].upper() in ['BLOCK','DECAY']:
                                break
                            ## smart select condition, if len(kk) < len(jj)-4, this line in SLHA file could not give desired info in any case.
                            if len(kk) < len(jj)-4:
                                continue
                            if jj[3].upper() == 'BLOCK' and ''.join(ks) ==  ''.join(kk[0:len(ks)]):
                                ks_flag  = True
                                sf.Debug('SLHA match data line', kk)
                            if jj[3].upper() == 'DECAY' and ''.join(ks) ==  ''.join(kk[1:len(ks)+1]):
                                ks_flag  = True
                                sf.Debug('SLHA match data line', kk)
                        if jj[3].upper() == kk[0].upper() and ''.join(blk) == ''.join(kk[1:len(blk)+1]) :
                            blk_flag = True
                            sf.Debug('SLHA match line',kk)
                            if jj[3].upper() == 'DECAY' and jj[5] == 0:
                                if len(kk) < 3 :
                                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, there are only %i column is the head line of "%s %s" in the input file.'%(jj[0],self._ProgName,len(kk),jj[3],jj[4],self._InputFile[ii]) )
                                else:
                                    sf.Debug('SLHA match data line', kk)
                                    ks_flag  = True
                                break
                        
                    if not blk_flag:
                        sf.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, can not find "%s %s" in the input file "%s".'%(jj[0], self._ProgName, jj[3], jj[4], self._InputFile[ii]) )

                    if not ks_flag:
                        sf.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, can not find required line with key "%s" in "%s %s" of the input file "%s".'%(jj[0], self._ProgName, jj[5:], jj[3], jj[4], self._InputFile[ii]) )

    def checkVar_replace(self, fileID, fileFlag):
            ## For 'Replace' method
            ## check the input vars that use 'Replace' method

            ii=fileID
            file_flag=fileFlag 

            ## if the input file is not obtained from previous program(s), open the file and check it.
            if file_flag:
                try :
                    ## Note infile and infile_bk stands for content of file
                    infile = open(self._InputFile[ii]).read()
                    ## if 'Replace' method is used and the last run of the same program stop by accident, there will exist a ESbackup file which is generated by easyscan and contains the replaced item.
                    if os.path.exists(self._InputFile[ii]+'.ESbackup'):
                        try:
                            infile_bk = open(self._InputFile[ii]+'.ESbackup').read()
                            BackupFlag = True
                        except:
                            BackupFlag = False
                    else:
                        BackupFlag = False
                ## Stop if easyscan could not read input file or input file.ESbackup.
                except:
                    sf.ErrorStop('Can not find/open the input file "%s" in program "%s".'%(self._InputFile[ii],self._ProgName))

            self._InRepVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'replace')]
           
            for jj in self._InRepVar[ii]:
                if len(jj) != 4 :
                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, 4 items ( Name ID,  Input file ID,  Method, Name of replaced parameter ) need to be provived.'%(jj[0], self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \tMethod= %s \tName= %s'%(jj[0],jj[1],jj[2],jj[3]))
               
                if file_flag:
                    ## re.findall used for returning a list of matching object.
                    match = re.findall(r"\b%s\b"%jj[3], infile)
                    if len(match)==0:
                        if not BackupFlag:
                            sf.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s" and its ESbackup file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                        else:
                            bk_match = re.findall(r"\b%s\b"%jj[3], infile_bk)
                            if len(bk_match)==0:
                                sf.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s" and its ESbackup file.'%(jj[0], self._ProgName, jj[3], self._InputFile[ii]) )
                            else:
                                ## the input file is wrong, use the ESbackup file
                                infile = infile_bk
                                match = re.findall(r"\b%s\b"%jj[3], infile)
                                ## there is no backup file for the following par
                                BackupFlag = False
                                sf.WarningNoWait( 'For input variable "%s" in program "%s" with "Replace" method, the input file "%s" does not contain "%s", but ESbackup file exists and contain it. In the following, ESbackup file will replace the content of the original input file.'%(jj[0],self._ProgName,self._InputFile[ii],jj[3]) )
                    if len(match)>1:
                        sf.WarningNoWait( 'For input variable "%s" in program "%s" with "Replace" method, find %i "%s" in coressponding input file "%s". They will all be replaced by variable "%s" in the following program.'%(jj[0],self._ProgName,len(match),jj[3],self._InputFile[ii],jj[0]) )
                    ## if the fist var do not use Backup file, the next var can not use
                    ## if the fist var use Backup, BackupFlag is already False
                    BackupFlag = False
            
            ## auto backup input file which involving replaced items, because the replaced items in origin input file would be replaced by values.
            if file_flag:
                if len(self._InRepVar[ii])>0:
                    open(self._InputFile[ii]+'.ESbackup','w').write(infile)
                    open(self._InputFile[ii],'w').write(infile)


    def setInputVar(self, inputvar):
        ## inputvar is a string of all content in input variable section of configure file
        self._InputVar=sf.string2nestlist(inputvar)
        
        for ii in self._InputVar:
            if len(ii) <4:
                if ii[0] == '': ## Program can have zero input parameters
                    return
                sf.ErrorStop( 'The input variables in program "%s" must have at least 4 items (Name ID,  Input file ID,  Method,  Note).'%self._ProgName )
            ## self._InFileID is file ID.
            if not ii[1] in self._InFileID:
                sf.ErrorStop( 'For input variable "%s" in program "%s", There is no input file with ID "%s"'%(ii[0],self._ProgName, ii[1]) )

            ## add for "math ..." in "Input variable" in [programX]
            self.invar[ii[0]] = sf.NaN

        sf.Info('Input variable = ')
        #file_flag = True
        for ii in self._InFileID:
            ## For 'File' method 
            ## check the input vars that use 'File' method
            file_flag=self.checkVar_file(ii)

            ## For 'Replace' method
            ## check the input vars that use 'Replace' method
            self.checkVar_replace(ii, file_flag)
            
            ## For 'Position' method
            ## check the input vars that use 'Position' method
            self.checkVar_position(ii, file_flag, "invar")
                         
            ## For 'Label' method
            ## check the input vars that use 'Label' method
            self.checkVar_label(ii, file_flag)
                                
            ## For 'SLHA' method
            ## check the input vars that use 'SLHA' method
            self.checkVar_slha(ii, file_flag)

    def setOutputFile(self, outputfile):
        outputfile=sf.string2nestlist(outputfile)
        self._OutFileID = [x[0] for x in outputfile ]
        sf.Debug('OutFileID',self._OutFileID)
        sf.Info('Output file     = ')
        for ii in outputfile:
            if not (ii[1].startswith('/home') or ii[1].startswith('~')):
                ii[1] = os.path.join(sf.CurrentPath, ii[1])
            self._OutputFile[ii[0]]=ii[1]
            sf.Info('  ID= %s \tFile= %s'%(ii[0],ii[1]))

    def setOutputVar(self, outputvar):
        outputvar=sf.string2nestlist(outputvar)
        for ii in outputvar:
            if not ii[1] in self._OutFileID:
                sf.ErrorStop( 'For output variable "%s" in program "%s", There is no output file with ID="%s".'%(ii[0],self._ProgName, ii[1]) )
            if not ii[2].upper() in ['FILE', 'POSITION', 'LABEL', 'SLHA']:
                sf.ErrorStop( 'For output variable "%s" in program "%s", EasyScan_HEP not supporting method="%s" now!'%(ii[0],self._ProgName, ii[2]) )
            self.outvar[ii[0]] = sf.NaN
        sf.Info('Output variable = ')
        for ii in self._OutFileID:
            self._OutputVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() != 'file')]

            ## For 'File' method
            self._OutFileVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'file')]
            if len(self._OutFileVar[ii])>1:
                sf.ErrorStop( 'In program "%s", there is no need to use more than one vars to stand the output file "%s" where you use %s.'%(self._ProgName, self._OutputFile[ii],' '.join(self._OutFileVar[ii])) )
            for jj in self._OutFileVar[ii]:
                if len(jj) != 4 :
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "File" method, 4 items ( Name, FileID, "File", Method ) need to be provived.'%(jj[0],self._ProgName) )
                if jj[3].upper() != 'SAVE' :
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "File" method, 4 items ( Name, FileID, "File", Method ) need to be provived. Note Method=save only supporting now.'%(jj[0],self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \t"File"= %s     Method=%s'%(jj[0],jj[1],jj[2],jj[3]))

            ## For 'Position' method
            self._OutPosVar[ii]  = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'position')]
            for jj in self._OutPosVar[ii]:
                if len(jj) != 5 :
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "Position" method, 5 items ( Name ID,  Input file ID,  Method, Line number,  Column number ) need to be provived.'%(jj[0],self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \tMethod= %s \tLine num= %s \tColumn num= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
        
            ## For 'label' method
            self._OutLabelVar[ii]  = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'label')]
            for jj in self._OutLabelVar[ii]:
                if len(jj) != 5 :
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "Label" method, 5 items ( Name ID,  Input file ID,  Method, Label name,  Input variable column number ) need to be provived.'%(jj[0],self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \tMethod= %s \tLabe=l %s \tColumn= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))

            ## For 'slha' method
            self._OutSLHAVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'slha')]
            for jj in self._OutSLHAVar[ii]:
                if len(jj) < 6 :
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "SLHA" method, at least 6 items ( Name ID,  Input file ID,  Method, BLOCK/DECAY, Block name/PDG, Keys) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['BLOCK','DECAY']:
                    sf.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, the 4th item must be "BLOCK" or "DECAY". If you can to used other formats, please contact with the authors.'%(jj[0],self._ProgName) )
                sf.Info('  varID= %s \tfileID= %s \tMethod= %s \tB/D= %s \tName= %s \tKeys= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4],jj[5]))
            
    def setRunFlag(self, name):
        self._runflag = name
        sf.Info('Run this program if '+name)

    def getProgName(self):
        return self._ProgName
    def getCommand(self):
        return self._Command
    def getComPath(self):
        return self._ComPath
    def getInputFile(self):
        return self._InputFile
    def getInputVar(self):
        return self._InputVar
    def getOutputFile(self):
        return self._OutputFile
    def getOutputVar(self):
        return self._OutputVar
    def getRunFlag(self,par):
        try:
          return eval(self._runflag)
        except:
          sf.ErrorStop('The Run Flag "'+self._runflag+'" is wrong for program '+ self._ProgName)

    def WriteInputFile(self,par):
        if self._InFileID ==['']:
            return

        ## add for "math ..." in "Input variable" in [programX]
        sf.parseMath(par)

        ## self._InFileID is list of file ID in Input file in [programX]
        for ii in self._InFileID:
            ## For 'File' method
            #file_flag = False
            file_flag = True
            for jj in self._InFilVar[ii]:
                #file_flag = True
                if jj[3].lower()=='previous':
                    file_flag = False 
                ## "save" not support now
                elif jj[3].lower()=='save':
                    pass
                elif jj[3].lower()=='replace':
                    sf.Debug("For program",self._ProgName)
                    sf.Debug("Copied file",par[jj[0]])
                    sf.Debug("Copy file",self._InputFile[ii])
                    shutil.copy(par[jj[0]],self._InputFile[ii])
                else:
                    try:
                        open(self._InputFile[ii],'a').write( open(par[jj[0]].read()) )
                    except:
                        sf.ErrorStop('Can not open input file "%s" or "%s" in program "%s", by "%s".'%(self._InputFile[ii], par[jj[0]], self._ProgName, self._InFilVar[ii]))

            ## Open the input file
            if file_flag:
                try:
                    if len(self._InRepVar[ii])>0:
                        infile = open(self._InputFile[ii]+'.ESbackup','r').read()
                    else:
                        infile = open(self._InputFile[ii],'r').read()
                except:
                    sf.ErrorStop('Can not open the input file "%s" in program "%s".'%(self._InputFile[ii], self._ProgName))
            else:
                try:
                    infile = open(self._InputFile[ii],'r').read()
                except:
                    sf.ErrorStop('Can not open the input file "%s" in program "%s", which is obtained from previous program(s).'%(self._InputFile[ii],self._ProgName))
        
            ## For 'Replace' method
            for jj in self._InRepVar[ii]:
                #if file_flag:
                if True:
                    match = re.findall(r"\b%s\b"%jj[3],infile)
                    if len(match)==0:
                        sf.ErrorStop('For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s", which is obtained from previous program(s).'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                ## jj[3] is something being replaced and par is a dictionary and par[jj[0]] (value) will replace jj[3].
                ## "\b" will make sure ES_lam in ES_lamT would not be replaced
                infile = re.sub(r"\b%s\b"%jj[3],str(par[jj[0]]),infile)
            if len(self._InRepVar[ii])>0:
                open(self._InputFile[ii],'w').write(infile)
                
            ## inlines is a list of all lines
            ## invars is a list of list of words for all lines
            inlines = open(self._InputFile[ii]).readlines()
            ## new 20180425 liang
            #invar = [ss.split() for ss in inlines]
            invar = [re.split(r'[ \t,]+', ss.strip()) for ss in inlines]

            ## invarNotModify for giving beutiful file as same as possible compared to the original input file.
            ## new 20180425 liang
            #invarNotModify = [ss.split() for ss in inlines]
            invarNotModify = [re.split(r'[ \t,]+', ss.strip()) for ss in inlines]

            ## if return, only support "replace" method
            #return
    
            ## For 'Position' method
            ## self_InPosVar[ii] is a list of input variable with position method
            for jj in self._InPosVar[ii]:
                invar[jj[3]-1][jj[4]-1] = str(par[jj[0]])

            ## For 'Label' method
            for jj in self._InLabVar[ii]:
                labelinum = [xxi for xxi,xx in enumerate(inlines) if jj[3] in xx]
            
                for kk in labelinum:
                    invar[kk][jj[4]-1] = str(par[jj[0]])

            ## For 'SLHA' method
            for jj in self._InSLHAVar[ii]:
               
                blk = str(jj[4]).split()
                blk_flag = False
                #ks  = str(jj[5]).split()
                ks  = map(str, jj[5:])
                ks_flag  = False
                for kki, kk in enumerate( invar ):
                    if kk[0].startswith('#'): continue
                    if blk_flag:
                        if kk[0].upper() in ['BLOCK','DECAY']:
                            break
                        if len(kk) < len(jj)-4:
                            continue
                        if jj[3].upper() == 'BLOCK' and ''.join(ks) ==  ''.join(kk[0:len(ks)]):
                            ks_flag  = True
                            invar[kki][len(ks)]=str(par[jj[0]])
                        if jj[3].upper() == 'DECAY' and ''.join(ks) ==  ''.join(kk[1:len(ks)+1]):
                            ks_flag  = True
                            invar[kki][0]=str(par[jj[0]])
                    if jj[3].upper() == kk[0].upper() and ''.join(blk) == ''.join(kk[1:len(blk)+1]) :
                        blk_flag = True
                        if jj[3].upper() == 'DECAY' and jj[5] == 0:

                            invar[kki][2]=str(par[jj[0]])
                            ks_flag  = True
                            break
          
            ## input file with replacing value, line not changed is same as origin one in format.
            outlines=[]
            for xxi,xx in enumerate(inlines):
              if invar[xxi] == invarNotModify[xxi]:
                 outlines.append( inlines[xxi] )
              else:
                 ## keep format unchanged 20180512
                 patt=re.compile(r'[ \t,\n]+')
                 joinList=patt.findall(inlines[xxi])
                 #print joinList, invar[xxi]; raw_input()
                 newList=[]
                 if len(joinList) == len(invar[xxi]):
                     for yyi in range(len(invar[xxi])):
                         newList.append(invar[xxi][yyi]) 
                         newList.append(joinList[yyi])
                 elif len(joinList)-1 == len(invar[xxi]):
                     for yyi in range(len(invar[xxi])):
                         newList.append(joinList[yyi])
                         newList.append(invar[xxi][yyi])
                     newList.append(joinList[-1])
                 else:
                     sf.ErrorStop("Keep format unchanged Failed! Check src/mainfun.py Line 559.") 
                 #outlines.append( "  ".join(invar[xxi])+"\n" )
                 outlines.append( "".join(newList))
            open(self._InputFile[ii],'w').writelines(outlines)

    def RunProgram(self):
        cwd=self._ComPath
        timeout = 60*2   # if the program run more than 2 hour, it may be killed
        for cmd in self._Command:
          sf.Debug('Runing Program %s with command'%self._ProgName,cmd)
          Use_os_system = True
          if Use_os_system:
              ncwd = os.getcwd()
              os.chdir(cwd)
              os.system(cmd[0])
              os.chdir(ncwd)
          else:
            try:
                p=subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, cwd=cwd,
                         shell=True)
                start = time.time()
                timepassed = 0
                while p.poll() is None:
                     time.sleep(0.1)
                     now = time.time()
                     timepassed = int((now-start)/60)
                     if timepassed > timeout:
                         sf.WarningWait("Program %s has run more than 1 hour, It will be kiled ")
                         try:
                             p.terminate()
                         except Exception,e:
                             return
                         return
                (out, err) = p.communicate()
                if p.stdin:  p.stdin.close()
                if p.stdout: p.stdout.close()
                if p.stderr: p.stderr.close()
                try:
                    p.kill()
                except OSError:
                    pass
                sf.Debug('Program %s done in %s minutes, corresponding output:'%(self._ProgName,timepassed),out)
            except OSError as error:
                if cwd and not os.path.exists(cwd):
                    raise Exception('Directory %s doesn\'t exist. Impossible to run' % cwd)
                else:
                    error_text = 'Impossible to run \'%s\'\n' % cmd
                    error_text += '    in the directory \'%s\'\n' % cwd
                    error_text += 'Trying to recover from:\n'
                    error_text += '    %s\n' % str(error)
                    raise Exception(err)
#        print p.returncode
#
#        if p.returncode:
#            error_msg = 'The run \'%s\' fails with following message:\n' % cmd
#            error_msg += '    '+out.replace('\n','\n    ')+'\n'
#            error_msg += 'Please try to fix this issue and retry.\n'
#            error_msg += 'Help might be found at mail@mail.com\n'
#
#            raise Exception(error_msg)
#
#        return p.returncode

    def RemoveOutputFile(self):
        for ii in self._OutFileID:
            if not os.path.exists(self._OutputFile[ii]):
                sf.Debug('No output file for program %s'%self._ProgName)
                sf.Debug('output file name is %s'%self._OutputFile[ii])
                return False
            else:
                os.remove(self._OutputFile[ii])

    def ReadOutputFile(self,par,path):
        for ii in self._OutFileID:
            if not os.path.exists(self._OutputFile[ii]):
                sf.Debug('No output file for program %s'%self._ProgName)
                sf.Debug('output file name is %s'%self._OutputFile[ii])
                return False
            ## For 'File' method
            for jj in self._OutFileVar[ii]:
                ##marked 20180420 liang
                #if jj[3].lower() == 'save' :
                #    self._Count += 1
                #    path = os.path.join(path,"SavedFile")
                #    SavePath = os.path.join(path, os.path.basename(self._OutputFile[ii]))+"."+str(self._Count)
                #    shutil.copy(self._OutputFile[ii],SavePath)
                par[jj[0]] = self._OutputFile[ii]
        
            if len(self._OutPosVar[ii])+ len(self._OutLabelVar[ii]) + len(self._OutSLHAVar[ii])>0 :
                oulines = open(self._OutputFile[ii]).readlines()
                ## new 20180425 liang
                #ouvar = [ss.split() for ss in oulines]
                ouvar = [re.split(r'[ \t,]+', ss.strip()) for ss in oulines]
            
            ## For 'Position' method
            for jj in self._OutPosVar[ii]:
                try:
                    par[jj[0]] = float(ouvar[jj[3]-1][jj[4]-1])
                    sf.Debug('Output - %s='%jj[0],par[jj[0]])
                except:
                    sf.Info('Can not read the output var %s'%jj)
                    return False

            ## For 'Label' method
            for jj in self._OutLabelVar[ii]:
#                labelinum = [xxi for xx,xxi in enumerate(oulines) if jj[3] in xx]
                labeline = [xx for xx in oulines if jj[3] in xx]
                if len(labeline)>1:
                    sf.ErrorStop( 'For output variable "%s" in program "%s" with "Label" method, there is %d "%s" in output file "%s". Please choose other method.'%(jj[0],self._ProgName,len(labelinum),jj[3],self._OutputFile[ii]) )

                try:
                    ## new 20180425 liang
                    #par[jj[0]] = float(labeline[0].split()[jj[4]-1])
                    par[jj[0]] = float(re.split(r'[ \t,]+',labeline[0].strip())[jj[4]-1])
                    sf.Debug('Output - %s='%jj[0],par[jj[0]])
                except:
                    sf.Debug('Can not read the output var',jj[0])
                    return False
                        
            ## For 'SLHA' method
            for jj in self._OutSLHAVar[ii]:

                blk = str(jj[4]).split()
                blk_flag = False
                ks  = str(jj[5]).split()
                ks_flag  = False
                for kki, kk in enumerate( ouvar ):
                    if kk[0].startswith('#'): continue
                    if blk_flag:
                        if kk[0].upper() in ['BLOCK','DECAY']:
                            break
                        if len(kk) < len(jj)-4:
                            continue
                        if jj[3].upper() == 'BLOCK' and ''.join(ks) ==  ''.join(kk[0:len(ks)]):
                            ks_flag  = True
                            par[jj[0]] = float(ouvar[kki][len(ks)])
                        if jj[3].upper() == 'DECAY' and ''.join(ks) ==  ''.join(kk[1:len(ks)+1]):
                            ks_flag  = True
                            par[jj[0]] = float(ouvar[kki][0])
                    if jj[3].upper() == kk[0].upper() and ''.join(blk) == ''.join(kk[1:len(blk)+1]) :
                        blk_flag = True
                        if jj[3].upper() == 'DECAY' and jj[5] == 0:
                            if len(kk) < 3 :
                                sf.Debug('Can not read the output var',jj)
                                return False
                            else:
                                par[jj[0]]=float(ouvar[kki][3])
                                ks_flag  = True
                            break
                sf.Debug('Output - %s='%jj[0],par[jj[0]])
                if not ks_flag:
                    sf.Debug('Can not read the output var',jj)
                    return False

        return True

    def SetOutput(self,par):
        for ii in self._OutFileID:
            for jj in self._OutputVar[ii]:
                par[jj[0]] = float("nan")
        return True

    def Recover(self):
        for ii in self._InFileID:
            if (ii!= '') and os.path.isfile(self._InputFile[ii]+".ESbackup"):
                os.system("mv %s.ESbackup %s" %(self._InputFile[ii],self._InputFile[ii]))

    def __str__(self):
        return '\
## [program] #############################################################\n\
Name    : %(name)s\n\
Command :\n~~~~ %(command)s\n\
Command Path:\n~~~~ %(cpath)s\n\
Input File: \n~~~~ %(inputfile)s\n\
Input: \n~~~~ %(inputvar)s\n\
Output File: \n~~~~ %(outputfile)s\n\
Output: \n~~~~ %(outputvar)s'\
%{'name':self._ProgName, 'command':self._Command, 'cpath':self._ComPath,
'inputfile': self._InputFile, 'inputvar': self._InputVar,
'outputfile': self._OutputFile, 'outputvar': self._OutputVar}

    ## new function to use "math .." in [constrain]
    ## in order to add new variable in self.AllPar
    def setGaussian(self,var):
        var = sf.string2nestlist(var)
        sf.Info('Gaussian Constraint:')
        for ii in var:
            if len(ii) in [3]:
                pass
            elif len(ii) in [4,5]:
                if not ii[3] in ['symm','lower','upper']:
                    sf.ErrorStop( 'For the "Gaussian" constraint on "%s", the "Type" can only be "upper"/"lower", not "%s".'%(ii[0],ii[3]) )
            else:
                sf.ErrorStop( 'The "Gaussian" constraint on "%s" need 4 or 5 items( VarID, Mean, Deviation, Type [, Name] ).'%(ii[0]) )

            self.cgauvar[ii[0]] = sf.NaN

    ## new 20180430 liang
    def setFreeFormChi2(self,var):
        var = sf.string2nestlist(var)
        sf.Info('FreeFormChi2:')
        for ii in var:
            if (len(ii)==1 and ii[0] ) or len(ii)==2:
                pass
            else:
                sf.ErrorStop( 'The "FreeFormChi2" constraint on "%s" need 1 item or 2 items( VarID [, Name] ).'%(ii[0]) )

            self.cffchi2var[ii[0]] = sf.NaN

    ## for "Bound" in [programX]
    def setBound(self, boundvar):
        ## boundvar is a string of all content in "Bound" in [programX] of configure file
        self._BoundVar=sf.string2nestlist(boundvar)

        sf.Info('Bound condition in program %s:'%self._ProgName)
        for ii in self._BoundVar:
            if len(ii) <3:
                if ii[0] == '': ## Program can have zero input parameters
                    return
                sf.ErrorStop( 'The "Bound" in program "%s" must have at least 3 items.'%self._ProgName )
            elif len(ii) == 3:
                if ii[1] not in ["<=", ">=", ">", "<", "==", "!="]:
                    sf.ErrorStop( 'The second item "%s" in "Bound" in program "%s" must be "<=", ">=", ">", "<", "==", "!=" or a real number at 3 items.'%(ii[1], self._ProgName) )
                try:
                    float(ii[2]) 
                except: 
                    self.boundvar[ii[2]] = sf.NaN

                self.boundvar[ii[0]] = sf.NaN

                sf.Info('    NameID= %s \tOperator= %s \tLimit= %s'%(ii[0],ii[1],ii[2]))

            elif len(ii) in [4,5]:
                if ii[2].lower() not in ["upper", "lower"]:
                    sf.ErrorStop( 'The third item "%s" in "Bound" in program "%s" must be "lower" or "upper" at 4 items.'%(ii[2], self._ProgName) )
                if not (ii[3].startswith('/home') or ii[3].startswith('~')):
                    ii[3]=os.path.join(sf.CurrentPath, ii[3])
                try:
                    open(ii[3])
                except:
                    sf.ErrorStop('Can not find/open the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                try:
                    boundfile = numpy.loadtxt(ii[3])
                except:
                    sf.ErrorStop('Find string in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                try:
                    numpy.shape(boundfile)[1]
                except:
                    sf.ErrorStop('Only one row or column in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                if numpy.shape(boundfile)[1] < 2:
                    sf.ErrorStop('Less than two columns in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                sf.Info('    NameID= %s \tNameID= %s \tBound= %s \tBoundFile= %s'%(ii[0],ii[1],ii[2],ii[3]))


                ## new 20180429 liang
                if len(ii) == 4:
                    jj = ii + ['Bound_%s_%s_%s'%(ii[0], ii[1], ii[2])]
                else:
                    jj = ii

                self.boundvar[jj[4]] = sf.NaN

            else: 
                sf.ErrorStop( 'The "Bound" in program "%s" have at most 5 items.'%self._ProgName )

                self.boundvar[ii[0]] = sf.NaN
                self.boundvar[ii[1]] = sf.NaN


    ## for "Bound" in [programX]
    ## ReadBound() have survived conditions in SetBound()
    def ReadBound(self,par):

        ## return if no bound
        ## If no bound, self._BoundVar = [['']], self._BoundVar[0][0]=''.
        if not self._BoundVar:
            return True
        if not self._BoundVar[0][0]:
            return True

        ## new 20180429 liang
        for ii in self._BoundVar:
            if len(ii) in [4,5]:
                if len(ii) == 4:
                    jj = ii + ['Bound_%s_%s_%s'%(ii[0], ii[1], ii[2])]
                else:
                    jj = ii

                if not (ii[3].startswith('/home') or ii[3].startswith('~')):
                    ii[3]=os.path.join(sf.CurrentPath, ii[3])
                boundfile = numpy.loadtxt(ii[3])
                x=boundfile[:,0]
                y=boundfile[:,1]
                if par[ii[0]] < numpy.amin(x) or par[ii[0]] > numpy.amax(x):
                    sf.WarningNoWati('"%s" less(greater) than min(max) of the first column in limit file "%s" with method "%s" in "Bound" in program "%s".'%(ii[0], ii[3], ii[2], self._ProgName))
                    if ii[2].lower() == 'lower':
                        yinter = -1E100
                    else:
                        yinter = 1E100
                    sf.WarningNoWait('    So we set "%s=%e"'%(jj[4], yinter))
                else:
                    yinter = numpy.interp(par[ii[0]], x, y)
                par[jj[4]] = yinter 
                sf.Debug('"%s=%f" v.s. "%s=%f" compare to the %s limit "%s=%f" by interplotion in "Bound" for program %s'%(ii[0], par[ii[0]], ii[1], par[ii[1]], ii[2].lower(), ii[1], yinter, self._ProgName))


        ## add for "math ..." in "Bound" in [programX]
        sf.parseMath(par)

        phy=True
        for ii in self._BoundVar:
            if len(ii)== 3:
                sf.Debug('"%s=%f" compare to the limit "%s" in "Bound" for program %s'%(ii[0], par[ii[0]], ii[1:], self._ProgName))
                try:
                    float(ii[2])
                except:
                    phy = phy and eval("%f%s%f"%(par[ii[0]], ii[1], par[ii[2]]))
                else:
                    phy = phy and eval("%f%s%s"%(par[ii[0]], ii[1], ii[2]))

        return phy 


class EasyScanInput:
    def __init__(self):
        self._FileName   = 'test'
        self._PointNum   = 100
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
            self._FileName = os.path.join(sf.CurrentPath, name)
        if self._ScanMethod not in ['READ','PLOT']:
            # deal with the satuation that the result file path already exists.
            if os.path.exists(self._FileName):
                print("* The Result file [%s] already exists." % name )
                while True:
                    c = raw_input("Choose: (r)eplace, (b)ackup, (s)top\n")
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
                        print("Wrong input! Please type in one of ('r','b','s')")
            # generate the result file path
            os.mkdir(self._FileName)
            os.mkdir(os.path.join(self._FileName,'SavedFile'))
            if self._ScanMethod == 'MULTINEST':
                self.MNOutputFile = os.path.join(self._FileName, "MultiNestData/")
                os.mkdir(self.MNOutputFile)
            sf.Info('...............................................')
            sf.Info('...............................................')
            sf.Info('Result file name  = %s'%self._FileName)
        else:
            if not os.path.exists(self._FileName):
                sf.ErrorStop("The result file %s doesn't exist."%self._FileName)

            if self._ScanMethod in ['READ'] and not os.path.exists(os.path.join(self._FileName,'ScanInfINPUT.txt')):
                sf.ErrorStop("The information file %s/ScanInfINPUT.txt doesn't exist. Please copy ScanInf.txt as ScanInfINPUT.txt and check the first line in ScanInfINPUT.txt carefully."%self._FileName)
                    
            if self._ScanMethod in ['PLOT'] and not os.path.exists(os.path.join(self._FileName,'ScanInf.txt')):
                sf.ErrorStop("The information file %s/ScanInf.txt doesn't exist. Please give ScanInf.txt and check the first line in ScanInf.txt carefully."%self._FileName)

            if self._ScanMethod in ['READ', 'PLOT']:
                sf.Info("* Now you are in READ (PLOT) mode. Files in SavedFile and Figure (Figure) in\n* %s\n* would be deleted and other files in\n* %s\n* not be deleted. Please choose the way how to deal with files in SavedFile and Figure (Figure)." % (self._FileName, self._FileName) )
                while True:
                    c = raw_input("Choose: (d)elete, (b)ackup, (s)top\n")
                    if c == "s":
                        exit(1)
                    elif c == "d":
                        os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'Figure'))
                        if self._ScanMethod == 'READ':
                            os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'SavedFile'))
                        break
                    elif c == "b":
                        if not (os.path.exists(sf.CurrentPath+"/Backup")):
                            os.mkdir(sf.CurrentPath+"/Backup")
                        BackupTime = time.strftime("_%Y_%m_%d_%H_%M_%S", time.localtime())
                        BackupPath = os.path.join(sf.CurrentPath, 'Backup/'+name+BackupTime)
                        os.system(r"cp -r %s %s" %(self._FileName, BackupPath))
                        os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'Figure'))
                        if self._ScanMethod == 'READ':
                            os.system(r"find %s -type f -name '*' | xargs rm" %os.path.join(self._FileName,'SavedFile'))
                        break
                    else:
                        sf.Info("Wrong input! Please type in one of ('c', 's')")
            sf.Info('...............................................')
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
        self._FlagTuneR = True
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
                    continue
                elif self._ScanMethod == 'GRID':
                    if   lenii == 4 :
                        self.GridBin[ii[0]]=20
                        sf.WarningNoWait('As the scan method "Grid", the bins number of the parameter [%s] is not provided, which will be set to default value, %i, in this run.'%(ii[0],self.GridBin[ii[0]]) )
                    elif lenii == 5:
                        self.GridBin[ii[0]]=ii[4]
                    else :
                        sf.WarningWait('For the scan method you choosed, only 5 items ( ID, Prior distribution, Minimum, Maximum, Bins number ) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],lenii,lenii-5) )
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
                    continue
                sf.Info('  ID= %s\tPrior= %s\tMin= %f\tMax= %f'%(ii[0],ii[1],ii[2],ii[3]))

            if self._ScanMethod == 'OnePoint':
                if lenii < 2 :
                    sf.ErrorStop('For the scan method you choosed, the items of each input parameter should include at least 2 items ( ID, Value). But the paramter [%s] missed %i of them.'%(ii[0],2-lenii))
                elif lenii > 2:
                    sf.WarningWait('For the scan method you choosed, only 2 items ( ID, Value) are needed for each input parameter. But the parameter [%s] has %i items. The last %i will be ignored.'%(ii[0],linii,2-lenii))
                continue

    # not used for now
    def setProgID(self,progID):
        self._ProgID = progID
        print progID
    
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
                if jj not in self.AllPar.keys():
                    self.AllPar[jj] = prog[ii].invar[jj]
                    self.OutPar[jj] = prog[ii].invar[jj]
            for jj in prog[ii].boundvar:
                if jj not in self.AllPar.keys():
                    self.AllPar[jj] = prog[ii].boundvar[jj]
                    self.OutPar[jj] = prog[ii].boundvar[jj]
            for jj in prog[ii].cgauvar:
                if jj not in self.AllPar.keys():
                    self.AllPar[jj] = prog[ii].cgauvar[jj]
                    self.OutPar[jj] = prog[ii].cgauvar[jj]
            for jj in prog[ii].cffchi2var:
                if jj not in self.AllPar.keys():
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

    def getFlagTuneR(self):
        return self._FlagTuneR
    def getAccepRate(self):
        return self._AccepRate

    def ScreenPrint(self,par,loglike):
        self._Count +=1
        if self._Count%self._PrintNum == 0:
            print '------------ Num: %i ------------'%self._Count
            for ii in par:
                print ii+' : '+ str(par[ii])
            print 'LogLike = '+str(loglike)
    
    
    def __str__(self):
        return '\
## [scan] ################################################################\n\
Name   : %(name)s\n\
Ntot   : %(ntot)d\n\
Method : %(method)s\n\
Iseed  : %(iseed)d\n\
Debug  : %(debug)r\n\
Input  :\n~~~~ %(inputvar)s'\
%{'name':self._FileName,'ntot':self._PointNum, 'method':self._ScanMethod,
'iseed':self._RandomSeed,'debug': self._DebugFlag,}

############################################################
#################   Constraint    Class   ##################
############################################################
class constraint:
    def __init__(self):
        self._Gaussian=[]
        #self._Limit=[]
        self._FreeFormChi2=[]
        self.Chi2={'Chi2':sf.NaN}
    
    def setGaussian(self,var):
        var = sf.string2nestlist(var)
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

    ## new 20180430 liang
    def setFreeFormChi2(self,var):
        var = sf.string2nestlist(var)
        sf.Info('FreeFormChi2:')
        for ii in var:
            if (len(ii)==1 and ii[0] ) or len(ii)==2:
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

        ## new 20180430 liang
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


    def setPlotPar(self,path,ScanMethod):
        ## try this
        # self._data =  np.loadtxt(path)
       
        if ScanMethod in ['READ']:
            f_data = open(os.path.join(path,'ScanInfINPUT.txt'),'r')
        else:
            f_data = open(os.path.join(path,'ScanInf.txt'),'r')
        path   = map(str,f_data.readline().split())
        var    = {}
        while True:
            plot_line = f_data.readline()
            if not plot_line :
                break
            plot_line = map(str,plot_line.split())
            var["+".join(plot_line[:-1])] = int(plot_line[-1])
   
        self._path = os.path.join(path[0],'Figure')
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
            line_par = map(str,line.split())
            for ii in var:
                try:
                    self._data[ii].append(float( line_par[var[ii]] ))
                except:
                    sf.Debug('Skip parameter %s'%ii)

        ##new 20180418 liang
        if ScanMethod not in ['READ', 'MULTINEST', 'PLOT']:
            for ii in var:
                self._dataAllTry[ii] = []
            
            f_dataAllTry = open(os.path.join(path[0],'All_%s'%path[1]), 'r')
            while True:
                line = f_dataAllTry.readline()
                if not line :
                    break
                line_par = map(str,line.split())
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
                    map(float, self._data[par[jj]])
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

    def getPlot(self,ScanMethod):
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

            if ScanMethod not in ['READ', 'MULTINEST', 'PLOT']:
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


