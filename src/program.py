####################################################################
#    Class PROGRAM: contral External programs                      #
####################################################################
# External modules
import os,sys,signal
import re,shutil
import subprocess
import numpy
import time
# Internal modules
import auxfun as af

class PROGRAM:
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
        self._executor = True

        self._outputclean = True
        self._timelimit = 60

    def setProgName(self, name):
        self._ProgName=name
        af.Info('...............................................')
        af.Info('Program name    = %s'% self._ProgName)
    def setCommand(self, command):
        self._Command=af.string2nestlist(command)
        af.Info('Execute command = %s'% self._Command)
    def setComPath(self, cpath):
        if cpath.startswith('/home') or cpath.startswith('~'):
            self._ComPath=cpath
        else:
            self._ComPath=os.path.join(af.CurrentPath, cpath)
        if not os.path.exists(self._ComPath):
            af.ErrorStop('Command path "%s" do not exist.'%self._ComPath)
        af.Info('Command path    = %s'% self._ComPath)

    def setInputFile(self, inputfile):
        inputfile=af.string2nestlist(inputfile)
        self._InFileID = [x[0] for x in inputfile ]
        if self._InFileID != list(set(self._InFileID)):
            af.ErrorStop('Input file in program "%s" have same File ID.'%self._ProgName)
        af.Info('Input file      = ')
        for ii in inputfile:
            if len(ii) != 2:
                if ii[0] == '':
                    break
                af.ErrorStop('The input file of %s need two items (File ID, File path).'%self._ProgName)
            if not (ii[1].startswith('/home') or ii[1].startswith('~')):
                ii[1]=os.path.join(af.CurrentPath, ii[1])
            self._InputFile[ii[0]]=ii[1]
            af.Info('  fileID= %s \tFile= %s'%(ii[0],ii[1]))

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
                    af.ErrorStop( 'For input variable "%s" in program "%s" with "File" method, 4 items (Name, FileID, "File", Method) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['PREVIOUS', 'SAVE', 'REPLACE', 'ADD']:
                    af.ErrorStop( 'For input variable "%s" in program "%s" with "File" method, the 4th item must be "PREVIOUS", "SAVE", "REPLACE" or "ADD". If you can use other formats, please contact with the authors.'%(jj[0], self._ProgName) )
                if jj[3].upper() == "PREVIOUS":
                    file_flag = False
                af.Info( 'Becasue file (ID=%s) in program "%s" is obtained from previous program(s), check this input file is ignored.'%(jj[1], self._ProgName))
                af.Info('  Name= %s \tFileID= %s \t"File"= %s \tMethod %s'%(jj[0],jj[1],jj[2],jj[3]))
           
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
                    af.ErrorStop('For input/output variable "%s" in program "%s" with "Position" method, 5 items ( Name ID, Input file ID, Method, Line number, Column number ) need to be provived.'%(jj[0], self._ProgName) )
                af.Info('  varID= %s \t fileID= %s \t Method= %s \t Line num= %s \t Column num= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]) )
                
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

                    af.Debug('Position len(invar)',len(invar))
                    af.Debug('Position line num',jj[3])
                    if len(invar) < jj[3]:
                        af.ErrorStop('For input variable "%s" in program "%s" with "Position" method, the line number is larger than the number of lines of inoput file "%s". Please check your input file.'%(jj[0],self._ProgName,self._InputFile[ii]) )
                    af.Debug('Position len(invar[line num])',len(invar[jj[3]-1]))
                    af.Debug('Position Column number',jj[4])

                    if len(invar[jj[3]-1]) < jj[4]:
                        af.ErrorStop('For input variable "%s" in program "%s" with "Position" method, the column number is larger than the number of columns in line %s of input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )


    def checkVar_label(self, fileID, fileFlag):
            ## For 'Label' method
            ## save the input vars that use 'Label' method

            ii=fileID
            file_flag=fileFlag 

            self._InLabVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'label')]
            
            for jj in self._InLabVar[ii]:
                if len(jj) != 5 :
                    af.ErrorStop('For input variable "%s" in program "%s" with "Label" method, 5 items ( Name ID,  Input file ID,  Method,  Label name,  Input variable column number ) need to be provived.'%(jj[0],self._ProgName) )
                if int(jj[4]) - jj[4] != 0 or jj[4] == 0:
                    af.ErrorStop('For input variable "%s" in program "%s" with "Label" method, the fourth item Input variable column number need to be an integer and not zero.'%(jj[0],self._ProgName) )

                af.Info('  varID= %s \tfileID= %s \tMethod= %s \tLabel= %s \tColumn= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
                
                if file_flag:
                    ## inlines is a list of all lines
                    ## invars is a list of list of words for all lines
                    inlines = open(self._InputFile[ii]).readlines()
                    invar = [ss.split() for ss in inlines]

                    ## enumerate return object that could by use by for loop by parsing list where xxi is id and xx is value in list.
                    labelinum = [xxi for xxi,xx in enumerate(inlines) if jj[3] in xx]
                    if len(labelinum)==0:
                        af.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, there is no "%s" in the input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                    if len(labelinum)!=1:
                        af.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, there is more than one line with label "%s" in the input file "%s". Please check your input file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                    for kk in labelinum:
                        af.Debug('Labeled line',inlines[kk].strip('\n'))
                        if len(invar[kk]) < jj[4]:
                            af.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, the column number "%i" defined by user is larger than the number of columns "%i" in the corresponding labeled line "%s" of input file "%s".'%(jj[0], self._ProgName, jj[4], len(invar[labelinum[0]]), invar[labelinum[0]], self._InputFile[ii]) )
                        if jj[4] > 0 and invar[kk][int(jj[4]-1)] == jj[3]:
                            af.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, the extracting item with column ID "%i" is just the label "%s" in the corresponding labeled line "%s" of input file "%s".'%(jj[0], self._ProgName, jj[4], jj[3], invar[labelinum[0]], self._InputFile[ii]) )
                        if jj[4] < 0 and invar[kk][int(jj[4])] == jj[3]:
                            af.ErrorStop( 'For input variable "%s" in program "%s" with "Label" method, the ex    tracting item with column ID "%i" is just the label "%s" in the corresponding labeled line "%s" of input file     "%s".'%(jj[0], self._ProgName, jj[4], jj[3], invar[labelinum[0]], self._InputFile[ii]) )

    def checkVar_slha(self, fileID, fileFlag):
            ## For 'SLHA' method
            ## check the input vars that use 'SLHA' method

            ii=fileID
            file_flag=fileFlag 

            self._InSLHAVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'slha')]

            for jj in self._InSLHAVar[ii]:
                if len(jj) < 6 :
                    af.ErrorStop('For input variable "%s" in program "%s" with "SLHA" method, at least 6 items ( Name ID,  Input file ID,  Method, BLOCK/DECAY, Block name/PDG, Keys) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['BLOCK','DECAY']:
                    af.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, the 4th item must be "BLOCK" or "DECAY". If you can use other formats, please contact with the authors.'%(jj[0],self._ProgName) )

                af.Info('  varID= %s \t fileID= %s \t Method= %s \t Block/Decay= %s \tName= %s \tKeys= %s'%(jj[0], jj[1], jj[2], jj[3], jj[4], jj[5:]))

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
                    ks  = list(map(str, jj[5:]))
                    ks_flag  = False
                    for kk in invar:
                        if len(kk)==0: continue
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
                                af.Debug('SLHA match data line', kk)
                            if jj[3].upper() == 'DECAY' and ''.join(ks) ==  ''.join(kk[1:len(ks)+1]):
                                ks_flag  = True
                                af.Debug('SLHA match data line', kk)
                        if jj[3].upper() == kk[0].upper() and ''.join(blk) == ''.join(kk[1:len(blk)+1]) :
                            blk_flag = True
                            af.Debug('SLHA match line',kk)
                            if jj[3].upper() == 'DECAY' and jj[5] == 0:
                                if len(kk) < 3 :
                                    af.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, there are only %i column is the head line of "%s %s" in the input file.'%(jj[0],self._ProgName,len(kk),jj[3],jj[4],self._InputFile[ii]) )
                                else:
                                    af.Debug('SLHA match data line', kk)
                                    ks_flag  = True
                                break
                        
                    if not blk_flag:
                        af.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, can not find "%s %s" in the input file "%s".'%(jj[0], self._ProgName, jj[3], jj[4], self._InputFile[ii]) )

                    if not ks_flag:
                        af.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, can not find required line with key "%s" in "%s %s" of the input file "%s".'%(jj[0], self._ProgName, jj[5:], jj[3], jj[4], self._InputFile[ii]) )

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
                    af.ErrorStop('Can not find/open the input file "%s" in program "%s".'%(self._InputFile[ii],self._ProgName))

            self._InRepVar[ii] = [x for x in self._InputVar if (x[1] == ii) and (x[2].lower() == 'replace')]
           
            for jj in self._InRepVar[ii]:
                if len(jj) != 4 :
                    af.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, 4 items ( Name ID,  Input file ID,  Method, Name of replaced parameter ) need to be provived.'%(jj[0], self._ProgName) )
                af.Info('  varID= %s \tfileID= %s \tMethod= %s \tName= %s'%(jj[0],jj[1],jj[2],jj[3]))
               
                if file_flag:
                    ## re.findall used for returning a list of matching object.
                    match = re.findall(r"\b%s\b"%jj[3], infile)
                    if len(match)==0:
                        if not BackupFlag:
                            af.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s" and its ESbackup file.'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
                        else:
                            bk_match = re.findall(r"\b%s\b"%jj[3], infile_bk)
                            if len(bk_match)==0:
                                af.ErrorStop( 'For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s" and its ESbackup file.'%(jj[0], self._ProgName, jj[3], self._InputFile[ii]) )
                            else:
                                ## the input file is wrong, use the ESbackup file
                                infile = infile_bk
                                match = re.findall(r"\b%s\b"%jj[3], infile)
                                ## there is no backup file for the following par
                                BackupFlag = False
                                af.WarningNoWait( 'For input variable "%s" in program "%s" with "Replace" method, the input file "%s" does not contain "%s", but ESbackup file exists and contain it. In the following, ESbackup file will replace the content of the original input file.'%(jj[0],self._ProgName,self._InputFile[ii],jj[3]) )
                    if len(match)>1:
                        af.WarningNoWait( 'For input variable "%s" in program "%s" with "Replace" method, find %i "%s" in coressponding input file "%s". They will all be replaced by variable "%s" in the following program.'%(jj[0],self._ProgName,len(match),jj[3],self._InputFile[ii],jj[0]) )
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
        self._InputVar=af.string2nestlist(inputvar)
        
        for ii in self._InputVar:
            if len(ii) <4:
                if ii[0] == '': ## Program can have zero input parameters
                    return
                af.ErrorStop( 'The input variables in program "%s" must have at least 4 items (Name ID,  Input file ID,  Method,  Note).'%self._ProgName )
            ## self._InFileID is file ID.
            if not ii[1] in self._InFileID:
                af.ErrorStop( 'For input variable "%s" in program "%s", There is no input file with ID "%s"'%(ii[0],self._ProgName, ii[1]) )

            ## add for "math ..." in "Input variable" in [programX]
            self.invar[ii[0]] = af.NaN

        af.Info('Input variable = ')
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
        if outputfile=='':
            af.Debug('No OutFile in program %s'%self._ProgName)
            self._OutFileID = []
            return 
        outputfile=af.string2nestlist(outputfile)
        self._OutFileID = [x[0] for x in outputfile ]
        af.Debug('OutFileID',self._OutFileID)
        af.Info('Output file     = ')
        for ii in outputfile:
            if not (ii[1].startswith('/home') or ii[1].startswith('~')):
                ii[1] = os.path.join(af.CurrentPath, ii[1])
            self._OutputFile[ii[0]]=ii[1]
            af.Info('  ID= %s \tFile= %s'%(ii[0],ii[1]))

    def setOutputVar(self, outputvar):
        if len(self._OutFileID)==0:
            af.Debug('No OutFile for program %s'%self._ProgName)
            return 
        outputvar=af.string2nestlist(outputvar)
        for ii in outputvar:
            if ii[1] not in self._OutFileID:
                af.ErrorStop( 'ID for output variable "%s" in program "%s" is wrong.'%(ii[0],self._ProgName, ii[1]))
            if ii[2].upper() not in ['FILE', 'POSITION', 'LABEL', 'SLHA']:
                af.ErrorStop( 'For output variable "%s" in program "%s", EasyScan_HEP not supporting method="%s" now!'%(ii[0],self._ProgName, ii[2]) )
            
            self.outvar[ii[0]] = af.NaN
            
        af.Info('Output variable = ')
        for ii in self._OutFileID:
            self._OutputVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() != 'file')]

            ## For 'File' method
            self._OutFileVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'file')]
            if len(self._OutFileVar[ii])>1:
                af.ErrorStop( 'In program "%s", there is no need to use more than one vars to stand the output file "%s" where you use %s.'%(self._ProgName, self._OutputFile[ii],' '.join(self._OutFileVar[ii])) )
            for jj in self._OutFileVar[ii]:
                if len(jj) != 4 :
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "File" method, 4 items ( Name, FileID, "File", Method ) need to be provived.'%(jj[0],self._ProgName) )
                if jj[3].upper() != 'SAVE' :
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "File" method, 4 items ( Name, FileID, "File", Method ) need to be provived. Note Method=save only supporting now.'%(jj[0],self._ProgName) )
                af.Info('  varID= %s \tfileID= %s \t"File"= %s     Method=%s'%(jj[0],jj[1],jj[2],jj[3]))

            ## For 'Position' method
            self._OutPosVar[ii]  = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'position')]
            for jj in self._OutPosVar[ii]:
                if len(jj) != 5 :
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "Position" method, 5 items ( Name ID,  Input file ID,  Method, Line number,  Column number ) need to be provived.'%(jj[0],self._ProgName) )
                if int(jj[4]) - jj[4] != 0 or jj[4] == 0:
                    af.ErrorStop('For output variable "%s" in program "%s" with "Label" method, the fourth item Input variable column number need to be an integer and not zero.'%(jj[0], self._ProgName) )

                af.Info('  varID= %s \tfileID= %s \tMethod= %s \tLine num= %s \tColumn num= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
        
            ## For 'label' method
            self._OutLabelVar[ii]  = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'label')]
            for jj in self._OutLabelVar[ii]:
                if len(jj) != 5 :
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "Label" method, 5 items ( Name ID,  Input file ID,  Method, Label name,  Input variable column number ) need to be provived.'%(jj[0],self._ProgName) )
                af.Info('  varID= %s \tfileID= %s \tMethod= %s \tLabe=l %s \tColumn= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))

            ## For 'slha' method
            self._OutSLHAVar[ii] = [x for x in outputvar if (x[1] == ii) and (x[2].lower() == 'slha')]
            for jj in self._OutSLHAVar[ii]:
                if len(jj) < 6 :
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "SLHA" method, at least 6 items ( Name ID,  Input file ID,  Method, BLOCK/DECAY, Block name/PDG, Keys) need to be provived.'%(jj[0],self._ProgName) )
                if not jj[3].upper() in ['BLOCK','DECAY']:
                    af.ErrorStop( 'For input variable "%s" in program "%s" with "SLHA" method, the 4th item must be "BLOCK" or "DECAY". If you can to used other formats, please contact with the authors.'%(jj[0],self._ProgName) )
                af.Info('  varID= %s \tfileID= %s \tMethod= %s \tB/D= %s \tName= %s \tKeys= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4],jj[5]))
            
    def setExecutor(self, executor):
        if executor.lower() == 'os.system':
            self._executor = True
            af.Info('Use "%s" execute commands.'%executor)
        elif executor.lower() == 'subprocess.popen':
            self._executor = False
            af.Info('Use "%s" execute commands.'%executor)
        else:
            af.Info('The command executor for program "%s" should be either "os.system" or "subprocess.popen", not "%s".'%(self._ProgName,executor))
            self._executor = True
            af.WarningNoWait('Use "os.system" execute commands.')

    def setOutputClean(self, outputclean):
        if outputclean.lower()  in ['yes','y','t','true']:
            self._outputclean = True
            af.Info('Delete the output file of program "%s" before execute it. '%self._ProgName)
        elif outputclean.lower()  in ['no','n','f','false']:
            self._outputclean = False
            af.Info('Keep the output file of program "%s" before execute it. '%self._ProgName)      
        else:
            af.WarningNoWait('The item "Output clean" for program "%s" should be either "Yes" or "No", not "%s".'%(self._ProgName,outputclean))
            self._outputclean = True
            af.Info('Delete the output file of program "%s" before execute it. '%self._ProgName)

    def setTimeLimit(self, timelimit):
        self._executor = False
        self._timelimit = timelimit
        af.Info('Time limit = %i minutes.'%self._timelimit)


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

    def WriteInputFile(self,par):
        if self._InFileID ==['']:
            return
            
        af.parseMath(par)

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
                    af.Debug("For program",self._ProgName)
                    af.Debug("Copied file",par[jj[0]])
                    af.Debug("Copy file",self._InputFile[ii])
                    shutil.copy(par[jj[0]],self._InputFile[ii])
                else:
                    try:
                        open(self._InputFile[ii],'a').write( open(par[jj[0]].read()) )
                    except:
                        af.ErrorStop('Can not open input file "%s" or "%s" in program "%s", by "%s".'%(self._InputFile[ii], par[jj[0]], self._ProgName, self._InFilVar[ii]))

            ## Open the input file
            if file_flag:
                try:
                    if len(self._InRepVar[ii])>0:
                        infile = open(self._InputFile[ii]+'.ESbackup','r').read()
                    else:
                        infile = open(self._InputFile[ii],'r').read()
                except:
                    af.ErrorStop('Can not open the input file "%s" in program "%s".'%(self._InputFile[ii], self._ProgName))
            else:
                try:
                    infile = open(self._InputFile[ii],'r').read()
                except:
                    af.ErrorStop('Can not open the input file "%s" in program "%s", which is obtained from previous program(s).'%(self._InputFile[ii],self._ProgName))
        
            ## For 'Replace' method
            for jj in self._InRepVar[ii]:
                #if file_flag:
                if True:
                    match = re.findall(r"\b%s\b"%jj[3],infile)
                    if len(match)==0:
                        af.ErrorStop('For input variable "%s" in program "%s" with "Replace" method, can not find "%s" in coressponding input file "%s", which is obtained from previous program(s).'%(jj[0],self._ProgName,jj[3],self._InputFile[ii]) )
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
                    if jj[4] > 0:
                        invar[kk][int(jj[4]-1)] = str(par[jj[0]])
                    else:
                        invar[kk][int(jj[4])] = str(par[jj[0]])

            ## For 'SLHA' method
            for jj in self._InSLHAVar[ii]:
              
                blk = str(jj[4]).split()
                blk_flag = False
                ks  = list(map(str, jj[5:]))
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
                     af.ErrorStop("Keep format unchanged Failed! Check src/program.py at Line about 559.") 
                 #outlines.append( "  ".join(invar[xxi])+"\n" )
                 outlines.append( "".join(newList))
            open(self._InputFile[ii],'w').writelines(outlines)

    def RunProgram(self):
        af.Debug('Be about to run Program %s'%self._ProgName)
        cwd=self._ComPath
        # remove output file
        if self._outputclean:
            self.RemoveOutputFile()
        for cmd in self._Command:
          af.Debug('Runing Program %s with command'%self._ProgName,cmd)

          if self._executor:
              ncwd = os.getcwd()
              os.chdir(cwd)
              os.system(cmd[0])
              os.chdir(ncwd)
          else:
            try:
                t_begining = time.time()
                process=subprocess.Popen(cmd, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, cwd=cwd, preexec_fn=os.setpgrp,
                        shell=True)
                while process.poll() == None:
                    process.stdout.flush()  
                    output = process.stdout.readline()
                    if output:
                        print(output.strip())
                    seconds_passed = time.time() - t_begining
                    if seconds_passed > self._timelimit*60: 
                        os.killpg(os.getpgid(process.pid), signal.SIGTERM) 
                        af.WarningWait("Program %s has been running more than %f minutes, It will be killed."%(self._ProgName, self._timelimit))
                        return
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
            af.Debug('Remove remaining output file %s before running program %s'%(self._OutputFile[ii], self._ProgName))
            if not os.path.exists(self._OutputFile[ii]):
                af.Debug('No remaining output file %s for program %s'%(self._OutputFile[ii], self._ProgName))
                return False
            else:
                os.remove(self._OutputFile[ii])
                af.Debug('Successful remaining output file %s for program %s'%(self._OutputFile[ii], self._ProgName))

    def ReadOutputFile(self,par,path):
        for ii in self._OutFileID:
            if not os.path.exists(self._OutputFile[ii]):
                af.Debug('No output file "%s" for program %s'%(self._OutputFile[ii], self._ProgName))
                return False
            ## For 'File' method
            for jj in self._OutFileVar[ii]:
                par[jj[0]] = self._OutputFile[ii]
        
            if len(self._OutPosVar[ii])+ len(self._OutLabelVar[ii]) + len(self._OutSLHAVar[ii])>0 :
                oulines = open(self._OutputFile[ii]).readlines()
                ouvar = [re.split(r'[ \t,]+', ss.strip()) for ss in oulines]
            
            ## For 'Position' method
            for jj in self._OutPosVar[ii]:
                try:
                    par[jj[0]] = float(ouvar[jj[3]-1][jj[4]-1])
                    af.Debug('Output - %s='%jj[0],par[jj[0]])
                except:
                    af.Info('Can not read the output var %s'%jj)
                    return False

            ## For 'Label' method
            for jj in self._OutLabelVar[ii]:
                labeline = [xx for xx in oulines if re.search(str(jj[3]),xx)]
                if len(labeline)>1:
                    af.ErrorStop( 'For output variable "%s" in program "%s" with "Label" method, there is %d "%s" in output file "%s". Please choose other method.'%(jj[0],self._ProgName,len(labelinum),jj[3],self._OutputFile[ii]) )

                try:
                    ## new 20180425 liang
                    #par[jj[0]] = float(re.split(r'[ \t]+',labeline[0].strip())[int(jj[4]-1)])
                    if jj[4] > 0:
                        par[jj[0]] = re.split(r'[ \t]+',labeline[0].strip())[int(jj[4]-1)]
                    else:
                        par[jj[0]] = re.split(r'[ \t]+',labeline[0].strip())[int(jj[4])]
                    af.Debug('Output - %s='%jj[0],par[jj[0]])
                except:
                    af.Debug('Can not read the output var',jj[0])
                    return False
                        
            ## For 'SLHA' method
            for jj in self._OutSLHAVar[ii]:
                blk = str(jj[4]).split()
                blk_flag = False
                ks  = list(map(str, jj[5:]))
                ks_flag  = False
                for kki, kk in enumerate( ouvar ):
                    if len(kk)==0: continue
                    if kk[0].startswith('#'): continue
                    if blk_flag:
                        if kk[0].upper() in ['BLOCK','DECAY']:
                            break
                        if len(kk) < len(jj)-4:
                            continue
                        if jj[3].upper() == 'BLOCK' and ''.join(ks) ==  ''.join(kk[0:len(ks)]):
                            ks_flag  = True
                            par[jj[0]] = float(ouvar[kki][len(ks)])
                        if jj[3].upper() == 'DECAY' and ''.join(ks).replace('.0','') ==  ''.join(kk[1:len(ks)+1]):
                            ks_flag  = True
                            par[jj[0]] = float(ouvar[kki][0])
                    if jj[3].upper() == kk[0].upper() and ''.join(blk) == ''.join(kk[1:len(blk)+1]) :
                        blk_flag = True
                        if jj[3].upper() == 'DECAY' and jj[5] == 0:
                            if len(kk) < 3 :
                                af.Debug('Can not read the output var',jj)
                                return False
                            else:
                                par[jj[0]]=float(ouvar[kki][2])
                                ks_flag  = True
                            break
                af.Debug('Output - %s='%jj[0],par[jj[0]])
                if not ks_flag:
                    if jj[3].upper() == 'DECAY':
                        af.Debug('Can not read the output var',jj)
                        af.Debug('In DECAY mode, set it as zero!')
                        par[jj[0]]=0
                    else:
                        af.Debug('Can not read the output var',jj)
                        return False

        return True

    def Recover(self):
        for ii in self._InFileID:
            if (ii!= '') and os.path.isfile(self._InputFile[ii]+".ESbackup"):
                os.system("mv %s.ESbackup %s" %(self._InputFile[ii],self._InputFile[ii]))

    ## new function to use "math .." in [constrain]
    ## in order to add new variable in self.AllPar
    def setGaussian(self,var):
        var = af.string2nestlist(var)
        af.Info('Gaussian Constraint:')
        for ii in var:
            if len(ii) in [3]:
                pass
            elif len(ii) in [4,5]:
                if not ii[3].lower() in ['symm','lower','upper']:
                    af.ErrorStop( 'For the "Gaussian" constraint on "%s", the "Type" can only be "symm", "upper" or "lower", not "%s".'%(ii[0],ii[3]) )
            else:
                af.ErrorStop( 'The "Gaussian" constraint on "%s" need 4 or 5 items( VarID, Mean, Deviation, Type [, Name] ).'%(ii[0]) )

            self.cgauvar[ii[0]] = af.NaN

    ## new 20180430 liang
    def setFreeFormChi2(self,var):
        var = af.string2nestlist(var)
        af.Info('FreeFormChi2:')
        for ii in var:
            if (len(ii)==1 and ii[0] ) or len(ii)==2:
                pass
            else:
                af.ErrorStop( 'The "FreeFormChi2" constraint on "%s" need 1 item or 2 items( VarID [, Name] ).'%(ii[0]) )

            self.cffchi2var[ii[0]] = af.NaN

    ## for "Bound" in [programX]
    def setBound(self, boundvar):
        ## boundvar is a string of all content in "Bound" in [programX] of configure file
        self._BoundVar=af.string2nestlist(boundvar)

        af.Info('Bound condition in program %s:'%self._ProgName)
        for ii in self._BoundVar:
            if len(ii) <3:
                if ii[0] == '': ## Program can have zero input parameters
                    return
                af.ErrorStop( 'The "Bound" in program "%s" must have at least 3 items.'%self._ProgName )
            elif len(ii) == 3:
                if ii[1] not in ["<=", ">=", ">", "<", "==", "!="]:
                    af.ErrorStop( 'The second item "%s" in "Bound" in program "%s" must be "<=", ">=", ">", "<", "==", "!=" or a real number at 3 items.'%(ii[1], self._ProgName) )
                try:
                    float(ii[2]) 
                except: 
                    self.boundvar[ii[2]] = af.NaN

                self.boundvar[ii[0]] = af.NaN

                af.Info('    NameID= %s \tOperator= %s \tLimit= %s'%(ii[0],ii[1],ii[2]))

            elif len(ii) in [4,5]:
                if ii[2].lower() not in ["upper", "lower"]:
                    af.ErrorStop( 'The third item "%s" in "Bound" in program "%s" must be "lower" or "upper" at 4 items.'%(ii[2], self._ProgName) )
                if not (ii[3].startswith('/home') or ii[3].startswith('~')):
                    ii[3]=os.path.join(af.CurrentPath, ii[3])
                try:
                    open(ii[3])
                except:
                    af.ErrorStop('Can not find/open the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                try:
                    boundfile = numpy.loadtxt(ii[3])
                except:
                    af.ErrorStop('Find string in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                try:
                    numpy.shape(boundfile)[1]
                except:
                    af.ErrorStop('Only one row or column in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                if numpy.shape(boundfile)[1] < 2:
                    af.ErrorStop('Less than two columns in the limit file "%s" in "Bound" in program "%s".'%(ii[3],self._ProgName))
                af.Info('    NameID= %s \tNameID= %s \tBound= %s \tBoundFile= %s'%(ii[0],ii[1],ii[2],ii[3]))


                ## new 20180429 liang
                if len(ii) == 4:
                    jj = ii + ['Bound_%s_%s_%s'%(ii[0], ii[1], ii[2])]
                else:
                    jj = ii

                self.boundvar[jj[4]] = af.NaN

            else: 
                af.ErrorStop( 'The "Bound" in program "%s" have at most 5 items.'%self._ProgName )

    ## for "Bound" in [programX]
    ## ReadBound() have survived conditions in SetBound()
    def ReadBound(self,par):

        ## return if no bound
        ## If no bound, self._BoundVar = [['']], self._BoundVar[0][0]=''.
        if not self._BoundVar:
            return True
        if not self._BoundVar[0][0]:
            return True

        ## add for "math ..." in "Bound" in [programX]
        af.parseMath(par)

        ## new 20180429 liang
        phy=True
        for ii in self._BoundVar:
            if len(ii)== 3:
                af.Debug('"%s=%f" compare to the limit "%s" in "Bound" for program %s'%(ii[0], par[ii[0]], ii[1:], self._ProgName))
                try:
                    float(ii[2])
                except:
                    phy = phy and eval("%f%s%f"%(par[ii[0]], ii[1], par[ii[2]]))
                else:
                    phy = phy and eval("%f%s%s"%(par[ii[0]], ii[1], ii[2]))
            elif len(ii) in [4,5]:
                if len(ii) == 4:
                    jj = ii + ['Bound_%s_%s_%s'%(ii[0], ii[1], ii[2])]
                else:
                    jj = ii
                if not (ii[3].startswith('/home') or ii[3].startswith('~')):
                    ii[3]=os.path.join(af.CurrentPath, ii[3])
                boundfile = numpy.loadtxt(ii[3])
                x=boundfile[:,0]
                y=boundfile[:,1]
                if par[ii[0]] < numpy.amin(x) or par[ii[0]] > numpy.amax(x):
                    af.WarningNoWait('"%s" less(greater) than min(max) of the first column in limit file "%s" with method "%s" in "Bound" in program "%s".'%(ii[0], ii[3], ii[2], self._ProgName))
                    if ii[2].lower() == 'lower':
                        yinter = af.log_zero
                    elif ii[2].lower() == 'upper':
                        yinter = -1.0*af.log_zero
                    af.WarningNoWait('    So we set "%s=%e"'%(jj[4], yinter))
                else:
                    yinter = numpy.interp(par[ii[0]], x, y)
                par[jj[4]] = yinter 
                af.Debug('"x-axis: %s=%f, y-axis: %s=%f" compare to the %s limit "y-interplotion: %s=%f" by interplotion in "Bound" for program %s'%(ii[0], par[ii[0]], ii[1], par[ii[1]], ii[2].lower(), ii[1], yinter, self._ProgName))

                if ii[2].lower() == "upper":
                    phy = phy and eval("%f%s%s"%(par[ii[1]], '<=', par[jj[4]]))
                elif ii[2].lower() == "lower":
                    phy = phy and eval("%f%s%s"%(par[ii[1]], '>=', par[jj[4]]))

        return phy 

