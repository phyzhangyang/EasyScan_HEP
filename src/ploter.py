####################################################################
#    Class PLOT: contral plot                                      #
####################################################################
# External modules
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata
# Internal modules
import init as sf


############################################################
#################   Plot    Class   ########################
############################################################


histconf={'bins':50, 'density':False, 'facecolor':'green', 'alpha':0.7}
scatterconf={'s':50, 'marker':'o', 'edgecolors':'None', 'alpha':0.9}
colorconf={'s':50, 'edgecolors':'None','cmap':plt.get_cmap('winter')}
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


