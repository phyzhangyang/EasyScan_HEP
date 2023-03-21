####################################################################
#    Class PLOT: contral plot                                      #
####################################################################
# Internal modules
import auxfun as af
from scan_controller import CONTROLLER 
# External modules
import os
import time
import numpy

############################################################
#################   Plot    Class   ########################
############################################################

figconf={'figsize':(7,7), 'dpi':80}
labelconf={'fontsize':20}
legendconf={'fontsize':20}

class PLOTER():
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
        hist=af.string2nestlist(hist)
        jj=0
        for ii in hist:
            if len(ii)==1 :
                self._Histogram.append([ii[0],'Histogram_%s.png'%ii[0]])
            elif len(ii)==2 :
                self._Histogram.append([ii[0],'%s.png'%ii[1]])

            self._FigNames.append(self._Histogram[jj][1])
            jj+=1
  
    def setScatter(self, scatter):
        scatter = af.string2nestlist(scatter)
        jj=0
        for ii in scatter:
            if len(ii)==2 :
                self._Scatter.append([ii[0],ii[1],'Scatter_%s_%s.png'%(ii[0],ii[1])])
            elif len(ii)==3 :
                self._Scatter.append([ii[0],ii[1],'%s.png'%ii[2]])

            self._FigNames.append(self._Scatter[jj][2])
            jj+=1

    def setColor(self, color):
        color = af.string2nestlist(color)
        jj=0
        for ii in color:
            if len(ii)==3 :
                self._Color.append([ii[0],ii[1],ii[2],'Color_%s_%s_%s.png'%(ii[0],ii[1],ii[2])])
            elif len(ii)==4 :
                self._Color.append([ii[0],ii[1],ii[2],'%s.png'%ii[3]])

            self._FigNames.append(self._Color[jj][3])
            jj+=1

    def setContour(self, Contour):
        Contour = af.string2nestlist(Contour)
        jj=0
        for ii in Contour:
            if len(ii)==3 :
                self._Contour.append([ii[0],ii[1],ii[2],'Contour_%s_%s_%s.png'%(ii[0],ii[1],ii[2])])
            elif len(ii)==4 :
                self._Contour.append([ii[0],ii[1],ii[2],'%s.png'%ii[3]])

            self._FigNames.append(self._Contour[jj][3])
            jj+=1


    def setPlotPar(self, path, ScanMethod, postprocess=False, ScanFile="", Plot=True):
        try:
            import pandas
        except ImportError:
            af.ErrorStop("No pandas module. No plot will be generated.")

        # read result
        af.Info("Read data for ONEPOINTBATCH or PLOT ...")
        datafile = ""
        if ScanMethod not in af._post:
            # define datafile
            if ScanMethod == af._onepointbatch:
                if Plot:
                    datafile = os.path.join(path, af.ResultFile) 
                else:
                    datafile = ScanFile
            elif ScanMethod == af._multinest:
                datafile = os.path.join(path, af.ResultFile_MultiNest)
            else:
                datafile = os.path.join(path, af.ResultFile)
            # read data from file
            if ScanMethod == af._multinest:
                column_names = pandas.read_csv(os.path.join(path, af.ResultFile), header=0, index_col=False).columns.values.tolist()
                self._data = pandas.read_csv(datafile, header=None, names=column_names, delim_whitespace=True, index_col=False)
            else:
                self._data = pandas.read_csv(datafile, header=0, index_col=False, sep="[,\s\t]+", engine='python')
            # read specific data for mcmc
            if ScanMethod == af._mcmc:
                self._dataAllTry = pandas.read_csv(os.path.join(path, af.ResultFile_MCMC), header=0, index_col=False)
        else: # ScanMethod == plot or postprocess 
            ResultFile_name = af.ResultFile_post if postprocess else af.ResultFile
            if os.path.exists(os.path.join(path, af.ResultFile_MultiNest)):
                column_names = pandas.read_csv(os.path.join(path, ResultFile_name), header=0, index_col=False).columns.values.tolist()
                datafile = os.path.join(path, af.ResultFile_MultiNest)
                self._data = pandas.read_csv(datafile, header=None, names=column_names, delim_whitespace=True, index_col=False)
            else:
                datafile = os.path.join(path, ResultFile_name)
                self._data = pandas.read_csv(datafile, header=0, index_col=False)
        
        # check whether data is empty
        if self._data.shape[0] == 0:
            af.ErrorStop("Your data file %s is empty."%(datafile))
            return False 
        
        # make figure folder
        self._path = os.path.join(path,'Figures')
        if not os.path.exists(self._path):
            os.mkdir(self._path)
        else:
            __import__("shutil").rmtree(self._path) 
            os.mkdir(self._path)

        af.Info("Read data for ONEPOINTBATCH or PLOT successfully in the file %s"%datafile)

    def checkPar(self, par, num, section_name='plot', severe=False):                  
      for jj in range(num):
#        try:
#          if self._data[par[jj]].min() == self._data[par[jj]].max():
#            af.WarningNoWait("Parameter %s=%f is a cosntant number. No plot for it. "%(par[jj], self._data[par[jj]].min()))
#            return False 
#        except KeyError:
#          af.WarningNoWait("Parameter '%s' in [plot] section do not exist. No plot for it."%( par[jj] )  )
#          return False 
        if par[jj] not in self._data.columns.values.tolist():
            if severe:
                af.ErrorStop("The parameter named '%s' of [%s] of the configure file could not be found in headers from your data file."%(par[jj], section_name))
            af.WarningNoWait("The parameter named '%s' of [%s] of the configure file could not be found in headers from your data file."%(par[jj], section_name))
            return False 
        if not numpy.issubdtype(self._data[par[jj]].dtype, numpy.number):
            if severe:
                af.ErrorStop("Column '%s' from your data file involves a non-number."%(par[jj]))
            af.WarningNoWait("Column '%s' from your data file involves a non-number."%(par[jj]))
            return False 
      return True

    ##only for debugging
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
                paths.append(numpy.vstack(xy))
            contours.append(paths)
    
        return contours

    def getPlot(self,ScanMethod):
        try:
            import matplotlib
        except ImportError:
            af.ErrorStop("No matplotlib module. No plot will be generated.")
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        if len(self._Histogram) + len(self._Scatter) + len(self._Color) + len(self._Contour) == 0:
            return
        af.Info('Start to plot the result ... ')
        FigNames=[]
        for ii in self._Histogram :
            histconf={'bins':50, 'normed':False, 'facecolor':'green', 'alpha':0.7}
            af.Info('Generate histogram plot of parameter %s'%ii[0])
            if not self.checkPar(ii,1): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            subplot.hist(self._data[ii[0]], **histconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel('Count', **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[1]))

        for ii in self._Scatter :
            scatterconf={'s':50, 'marker':'o', 'edgecolors':'None', 'alpha':0.9}
            af.Info('Generate scatter plot of parameter %s and %s'%(ii[0],ii[1]))
            if not self.checkPar(ii,2): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            subplot.scatter(self._data[ii[0]],self._data[ii[1]],**scatterconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[2]))

            if ScanMethod == af._mcmc:
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
            colorconf={'edgecolors':'None','cmap':plt.get_cmap('winter')}
            af.Info('Generate color plot of parameter %s and %s with color %s'%(ii[0],ii[1],ii[2]))
            if not self.checkPar(ii,3): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)
            weigh = 50
            if "probability" in self._data.columns.values.tolist():
              weigh = self._data["probability"] * 100 / self._data["probability"].max() + 20
            elif "mult" in self._data.columns.values.tolist():
              weigh = self._data["mult"] * 100 / self._data["mult"].max() + 20
            sc1=subplot.scatter(self._data[ii[0]], self._data[ii[1]], c=self._data[ii[2]], s=weigh, **colorconf)
            cb1=plt.colorbar(sc1)
            cb1.set_label(ii[2], **labelconf)
            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[3]))

        for ii in self._Contour :
            try:
                from scipy.interpolate import griddata
            except ImportError:
                af.WarningNoWait("No scipy module. Contour plot will not be generated.")
                break
            af.Info('Generate contour plot of parameter %s and %s with contour %s'%(ii[0],ii[1],ii[2]))
            if not self.checkPar(ii,3): continue
            f=plt.figure(**figconf)
            subplot=f.add_subplot(111)

            x = self._data[ii[0]]
            X = numpy.linspace(min(x),max(x),100)
            y = self._data[ii[1]]
            Y = numpy.linspace(min(y),max(y),100)
            # z = [ numpy.log10(abs(u)) for u in self._data[ii[2]] ]  # log10
            z = self._data[ii[2]] 
            Z = griddata(x,y,z,X,Y,interp='linear')

            #C = subplot.contour(X,Y,Z, 3,linewidths=2)
            #plt.clabel(C, inline=True, fontsize=8)

            ## debug
            #Cpoint = self.get_contour_verts(C)  
            #numpy.savetxt(os.path.join(self._path, "contour_1_1.dat"),Cpoint[0][0])

            C = subplot.contourf(X,Y,Z,3, cmap=plt.cm.rainbow)
            
            cb1=plt.colorbar(C)
            cb1.set_label(ii[2], **labelconf)

            subplot.set_xlabel(ii[0], **labelconf)
            subplot.set_ylabel(ii[1], **labelconf)
            subplot.tick_params(which = 'both', direction = 'out')
            plt.savefig(os.path.join(self._path, ii[3]))


