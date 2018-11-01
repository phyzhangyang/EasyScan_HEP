####################################################################
#    Class CONSTRAINT: contral CONSTRAINT                          #
####################################################################
# External modules
# Internal modules
import init as sf

class CONSTRAINT:
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
