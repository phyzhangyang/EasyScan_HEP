####################################################################
#    Read in the input configure file
####################################################################

# External modules.
import mainfun as mf
import init  as sf
import ConfigParser


def ReadIn(Configfile,ES,Programs,CS,Ploter):
    cf=ConfigParser.ConfigParser()
    cf.read(Configfile)

    ## Read the basic scan parameters
    try:
        ES.setScanMethod(cf.get('scan', 'Scan method'))
    except ConfigParser.NoOptionError:
        sf.WarningWait('Can not find "Scan method" in the input configure file, it will take the default value, "Random".')

    try:
        ES.setFileName(cf.get('scan', 'Result file name'))
    except ConfigParser.NoOptionError:
        sf.ErrorStop('Please provide "Result file name" in the input configure file.')

    ## If the scan method is 'plot', go directly to the 'plot' section
    ## Read the plot information
    plot_items = []
    try: 
        plot_items  = cf.options("plot")
    except ConfigParser.NoSectionError:
        if cf.get('scan', 'Scan method').upper() in ["PLOT"]:
            sf.ErrorStop('ConfigParser.NoSectionError: No section: [plot] in the configure file.')
        
    if 'histogram' in plot_items:
        Ploter.setHistogram(cf.get('plot', 'Histogram'))
    if 'scatter' in plot_items:
        Ploter.setScatter(cf.get('plot', 'Scatter'))
    if 'color' in plot_items:
        Ploter.setColor(cf.get('plot', 'Color'))
    if 'contour' in plot_items:
        Ploter.setContour(cf.get('plot', 'Contour'))
    # new 20180419 liang
    sf.checkItemInList(Ploter._FigNames)
    if ES.getScanMethod() == 'PLOT': return float('nan')

    ## Back to Read the basic scan parameters
    try:
        #ES.setPointNum(cf.getint('scan', 'Number of points'))
        ES.setPointNum(cf.getint('scan', 'Points number'))
	if cf.get('scan', 'Scan method').upper() in ["GRID"]:
            sf.WarningNoWait('"Points number" in the input configure file is not used in "GRID" scan mode.')
    except ConfigParser.NoOptionError:
	if cf.get('scan', 'Scan method').upper() in ["GRID"]:
            pass
        sf.WarningWait('Can not find "Points number" in the input configure file, it will take the default value, 100.')
    except ValueError:
        sf.ErrorStop('The "Points number" in the input configure file must be an integer.')

    try:
        ES.setRandomSeed(cf.getint('scan', 'Random seed'))
        sf.Info('"Random seed = %d" in the input configure file, using only in Random, MCMC or Multinest mode.'%cf.getint('scan', 'Random seed'))
    except ConfigParser.NoOptionError:
        #new 20180425 liang
        if cf.get('scan', 'Scan method').upper() in ["RANDOM","MCMC", "MULTINEST"]:
           sf.WarningWait('Can not find "Random seed" in the input configure file, it will take current system time as ramdom seed.')
    except ValueError:
        sf.ErrorStop('The "Random seed" in the input configure file must be an integer.')


    try:
        ES.setPrintNum(cf.getint('scan', 'Interval of print'))
    except ConfigParser.NoOptionError:
        sf.WarningNoWait('Can not find "Interval of print" in the input configure file, it will take the default value, 1.')
    except ValueError:
        sf.WarningNoWait('The "Interval of print" in the input configure file must be an integer, it will take the default value, 1.')

    try:
        ES.setInputPar(cf.get('scan', 'Input parameters'))
    except ConfigParser.NoOptionError:
        sf.ErrorStop('Can not find "Input parameters" in the input configure file.')

    try:
        ES.setAccepRate (cf.get('scan', 'Acceptance rate'))
    except ConfigParser.NoOptionError:
        pass

    ## Read the program(s) setting
    ProgID=filter(lambda x: x.startswith('program'), cf.sections())
    ProgID.sort(key=lambda x: int(filter(str.isdigit, x)) )

    # new 20180419 liang
    outputVarNames=[]   

    for ii in ProgID:
        items  = cf.options(ii)
        fullitems = ['program name', 'execute command', 'command path', 'input file', 'input variable', 'output file', 'output variable']
        for jj in fullitems:
            if not jj in items:
                sf.ErrorStop('For "%s", item "%s" missed.'%(ii,jj))
        Programs[ii] = mf.program()
        Programs[ii].setProgName(cf.get(ii, 'Program name'))
        Programs[ii].setCommand(cf.get(ii, 'Execute command'))
        Programs[ii].setComPath(cf.get(ii, 'Command path'))
        Programs[ii].setInputFile(cf.get(ii, 'Input file'))
        Programs[ii].setInputVar(cf.get(ii, 'Input variable'))
        Programs[ii].setOutputFile(cf.get(ii, 'Output file'))
        Programs[ii].setOutputVar(cf.get(ii, 'Output variable'))
        ##new 20180419 liang
        for key, item in Programs[ii]._OutputVar.items():
            for subitem in item:
                outputVarNames.append(subitem[0])
        sf.checkItemInList(outputVarNames)
        ## for "Bound" in [programX]
        ## This is not necessary for evry program -- Yang 05.08
        try:
            Programs[ii].setBound(cf.get(ii, 'Bound'))
        except:
            sf.Info('No Bound.')
        
        try:
            Programs[ii].setRunFlag(cf.get(ii, 'Run flag'))
        except:
            sf.Info('No Run flag.')

    ## Read the constraints
    # new 20180426 liang
    constraint_items = []
    try:
        # auto get lowered items by cf.options("constraint")
        constraint_items  = cf.options("constraint")
        if len(constraint_items) == 0:
            sf.ErrorStop('No item in the section [constraint] in the configure file!\n* For "GRID", "RANDOM" or "READ" scan mode, items could be calculated but not guide scanning. For "MCMC" and "MULTINEST" scan mode, items are calculated and guide scanning.')
    except ConfigParser.NoSectionError:
        #if cf.get('scan', 'Scan method').upper() in ["MCMC", "MULTINEST"]:
        sf.ErrorStop('ConfigParser.NoSectionError: No section: [constraint] in the configure file.')
            
    sf.Info('...............................................')
    sf.Info('...............................................')
    sf.Debug('constraint_items',constraint_items)
    if 'gaussian' in constraint_items:
        CS.setGaussian(cf.get('constraint', 'Gaussian'))
        ## In order to add "math .." in "Gaussian" to self.AllPar
        if Programs:
            Programs[ProgID[0]].setGaussian(cf.get('constraint', 'Gaussian'))
    
    ## marked 20180430 liang
    #if 'limit' in constraint_items:
    #    CS.setLimit(cf.get('constraint', 'Limit'))
    ## new 20180430 liang
    elif 'freeformchi2' in constraint_items:
        CS.setFreeFormChi2(cf.get('constraint', 'FreeFormChi2'))
        if Programs:
            Programs[ProgID[0]].setFreeFormChi2(cf.get('constraint', 'FreeFormChi2'))

    else:
        sf.ErrorStop('Only support "GAUSSIAN" and "FREEFORMCHI2" in [constraint] in the configure file.')
    #################################
    ## manage the vars into ES.AllPar
    #    ES.setProgID(ProgID)
    ES.setProgram(Programs)

    return ProgID
