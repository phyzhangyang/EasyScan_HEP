####################################################################
#    Readin configure file
####################################################################
# Internal modules
from program import PROGRAM
import auxfun as af
# External modules
import configparser
import sys,time 

# Check reduplicative variable name
def checkDuplicatedName(List, name, check_forbidden=True):
    for item in List:
      if check_forbidden:
        if item in af.forbidden_names:
          af.ErrorStop('%s name "%s" can be used as variable name.'%(name, item))
      counter = List.count(item)
      if counter>1:
        af.ErrorStop('%s name "%s" duplicates %i times in configure file.'%(name, item, counter))


def notFind(name):
  return 'Can not find "%s" in configure file. '%name

def notInteger(name):
  return '"%s" in configure file is not an integer.'%name

def takeDefault(name):
  return 'It will take default value, "%s".'%name

def ReadIn(Configfile, ES, Programs, Constraint, Ploter):
    cf=configparser.ConfigParser()
    cf.read(Configfile)

    # Read the basic scan parameters
    if not ('scan' in cf.sections()) :
        af.ErrorStop(notFind('[scan] section')) 
    try:
        ES.setScanMethod(cf.get('scan', 'Scan method'))
    except configparser.NoOptionError:
        af.ErrorStop(notFind('Scan method'))
    try:
        ES.setFolderName(cf.get('scan', 'Result folder name'))
    except configparser.NoOptionError:
        af.ErrorStop(notFind('Result folder name')) 


    # Read the plot information
    # This part is put here, because if ScanMethod = PLOT, the rest of sections can be ignore.
    plot_items = []
    try: 
        plot_items  = cf.options("plot")
    except configparser.NoSectionError:
        if ES.getScanMethod() == af._plot:
            af.ErrorStop(notFind('[plot] section')) 
    if 'histogram' in plot_items:
        Ploter.setHistogram(cf.get('plot', 'Histogram'))
    if 'scatter' in plot_items:
        Ploter.setScatter(cf.get('plot', 'Scatter'))
    if 'color' in plot_items:
        Ploter.setColor(cf.get('plot', 'Color'))
    if 'contour' in plot_items:
        Ploter.setContour(cf.get('plot', 'Contour'))
    checkDuplicatedName(Ploter._FigNames, "Figure", False)
    if ES.getScanMethod() == af._plot: return []

    # Back to read the basic scan parameters
    try:
        ES.setPointNum(cf.getint('scan', 'Number of points'))
        if ES.getScanMethod() in af._no_random:
            af.WarningNoWait('"Number of points" in configure file is not used.')
    except configparser.NoOptionError:
        if ES.getScanMethod() in af._no_random:
            pass
        else:
            af.WarningWait(notFind('Number of points')+takeDefault('10')) 
    except ValueError:
        af.ErrorStop(notInteger("Number of points"))
    try:
        ES.setRandomSeed(cf.getint('scan', 'Random seed'))
    except configparser.NoOptionError:
        if ES.getScanMethod() not in af._no_random:
           af.Info("Use current system time as random seed.")
    except ValueError:
        af.ErrorStop(notInteger("Random seed"))

    try:
        ES.setPrintNum(cf.getint('scan', 'Interval of print'))
    except configparser.NoOptionError:
        af.Info(notFind('Interval of print')+takeDefault('1'))
        ES.setPrintNum('1')
    except ValueError:
        af.WarningNoWait(notInteger("Interval of print")+takeDefault('1'))
        ES.setPrintNum('1')

    try:
        ES.setInputPar(cf.get('scan', 'Input parameters'))
    except configparser.NoOptionError:
        af.ErrorStop(notFind('Input parameters'))
    if ES.getScanMethod() in [af._onepoint, af._onepointbatch]:
        ES.InPar = ES.FixedPar
        ES.FixedPar = {}
    checkDuplicatedName(list(ES.InPar.keys()), "Input parameter")

    try:
        ES.setAccepRate(cf.get('scan', 'Acceptance rate'))
    except configparser.NoOptionError:
        if ES.getScanMethod() == af._mcmc:
            af.Info(notFind('Acceptance rate')+takeDefault('0.25')) 
        else:
            pass

    try:
        ES.setParallelThreads(cf.get('scan', 'Parallel threads'))
    except configparser.NoOptionError:
        af.Info('Parallel mode is off.')

    if ES.getParallelMode():
        try:
            ES.setParallelFolder(cf.get('scan', 'Parallel Folder'))
        except configparser.NoOptionError:
            af.ErrorStop(notFind('Parallel Folder'))
        
    
    

    ## sort programs by ID
    ProgID = [x for x in cf.sections() if x.startswith('program')]
    if len(ProgID)==0 :
        af.ErrorStop(notFind('[program] section'))
    for ii in ProgID:
        if not str.isdigit(ii[7:]):
            af.ErrorStop('The section name of [%s] is wrong'%ii)
    ProgID = sorted(ProgID, key=lambda x: int( x[7:] ) )


    ## Read the programs sections
    outputVarNames=[]   
    for ii in ProgID:
        items  = cf.options(ii)
        
        Programs[ii] = PROGRAM()
        if 'program name' in items:
            Programs[ii].setProgName(cf.get(ii, 'Program name'))
        else:
            Programs[ii].setProgName(ii)
        if 'execute command' in items:
            Programs[ii].setCommand(cf.get(ii, 'Execute command'))
        else:
            af.ErrorStop('No "Execute command" in "%s".'%ii)
        if 'command path' in items:
            Programs[ii].setComPath(cf.get(ii, 'Command path'))
        else:
            Programs[ii].setComPath('')
            af.WarningNoWait('Use current path as "Command path" for "%s".'%ii)
        
        if 'input file' in items and 'input variable' in items:
            Programs[ii].setInputFile(cf.get(ii, 'Input file'))
            Programs[ii].setInputVar(cf.get(ii, 'Input variable'))
        if 'output file' in items and 'output variable' in items:
            Programs[ii].setOutputFile(cf.get(ii, 'Output file'))
            Programs[ii].setOutputVar(cf.get(ii, 'Output variable'))
        
        for key, item in Programs[ii]._OutputVar.items():
            for subitem in item:
                outputVarNames.append(subitem[0])
        checkDuplicatedName(outputVarNames, "Output variable")
        
        # Additional optional commands
        try:
            Programs[ii].setExecutor(cf.get(ii, 'Command executor'))
        except:
            af.Info('Use "os.system" execute commands.')
        try:
            Programs[ii].setOutputClean(cf.get(ii, 'Clean output file'))
        except:
            af.Info('Delete output file(s) of %s before execute it. '%ii)
        try:
            Programs[ii].setBound(cf.get(ii, 'Bound'))
        except SystemExit:
            sys.exit(1)  
        except:
            af.Info('No Bound.')
        
        if Programs[ii]._executor: # 'os.system'
          try:
            Programs[ii].setTimeLimit(cf.getfloat(ii, 'Time limit in minute'))
            af.WarningNoWait('Time limit of %s is not used.'%ii)
          except:
            pass
        else: # 'subprocess.popen'
          try:
            Programs[ii].setTimeLimit(cf.getfloat(ii, 'Time limit in minute'))
          except: 
            af.Info('No time limit setting. Using defaut value 60 min.')
        
    ## Read the constraints
    constraint_items = []
    try:
        constraint_items  = cf.options("constraint")
        if len(constraint_items) == 0 and (ES.getScanMethod() not in af._no_like):
          af.ErrorStop('Section [constraint] is empty in the configure file.')
    except configparser.NoSectionError:
        if ES.getScanMethod() in af._no_like:
            pass 
        else:
            af.ErrorStop(notFind('[constraint] section'))
            
    af.Info('...............................................')
    af.Debug('Constraint items:',constraint_items)
    if 'gaussian' in constraint_items:
        Constraint.setGaussian(cf.get('constraint', 'Gaussian'))
        # In order to add "math .." in "Gaussian" to self.AllPar TODO
        if Programs:
            Programs[ProgID[0]].setGaussian(cf.get('constraint', 'Gaussian'))
    if 'freeformchi2' in constraint_items:
        Constraint.setFreeFormChi2(cf.get('constraint', 'FreeFormChi2'))
        if Programs:
            Programs[ProgID[0]].setFreeFormChi2(cf.get('constraint', 'FreeFormChi2')) 
    if ('gaussian' not in constraint_items) and ('freeformchi2' not in constraint_items) and (ES.getScanMethod() not in af._no_like):
        af.ErrorStop('No valid iterm in [constraint] section.')

    ES.setProgram(Programs)

    return ProgID
