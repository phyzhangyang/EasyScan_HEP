####################################################################
#    Class CONSTRAINT: contral CONSTRAINT                          #
####################################################################
# Internal modules
import auxfun as af

class CONSTRAINT:
    def __init__(self):
        self._Gaussian=[]
        #self._Limit=[]
        self._FreeFormChi2=[]
        self.Chi2={}
        self.Chi2['Chi2'] = af.NaN 
    
    def setGaussian(self,var):
        var = af.string2nestlist(var)
        af.Info('Gaussian Constraint:')
        for ii in var:
            if len(ii) in [3]:
                jj = ii+['symm','Gaussian_%s'%ii[0]]
                self._Gaussian.append(jj)
                self.Chi2[jj[4]] = af.NaN 
                
                af.Info('    varID= %s\tMean= %e\tDeviation= %e\tType= %s\tName= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))
            elif len(ii) in [4,5]:
                if not ii[3].lower() in ['symm','lower','upper']:
                    af.ErrorStop( 'For the "Gaussian" constraint on "%s", the "Type" can only be "symm", "upper" or "lower", not "%s".'%(ii[0],ii[3]) )
                ## new 20180428 liang
                if len(ii) == 4: 
                    jj = ii+['Gaussian_%s'%ii[0]]
                else:
                    jj = ii
                self._Gaussian.append(jj)
                self.Chi2[jj[4]] = af.NaN 
                
                af.Info('    varID= %s\tMean= %e\tDeviation= %e\tType= %s\tName= %s'%(jj[0],jj[1],jj[2],jj[3],jj[4]))

            else:
                af.ErrorStop( 'The "Gaussian" constraint on "%s" need 4 or 5 items( ID, Mean, Deviation, Type [, Name] ).'%(ii[0]) )

        ## new 20180428 liang
        self.Chi2 = af.sortDic(self.Chi2) 

## marked 20180430 liang
#    def setLimit(self,var):
#        var = af.string2nestlist(var)
#        af.Info('Upper/Lower limit:')
#        for ii in var:
#            if len(ii) == 4:
#                self._Limit.append(ii)
#                af.Info('  varID(X)= %s\tvarID(Y)= %s\tConstraintFile= %s\tType= %s'%(ii[0],ii[1],ii[2],ii[3]))
#        ## add check the ConstraintFile exist
#        ## add check the ConstraintFile has two columns and >1 lines
#        ## very useful, you can simply let a>b
#            else:
#                af.ErrorStop( 'The "Limit" constraint on "(%s,%s)" need 4 items( X ID, Y ID, ConstraintFile, Type ).'%(ii[0],ii[1]) )

    ## new 20180430 liang
    def setFreeFormChi2(self,var):
        var = af.string2nestlist(var)
        af.Info('FreeFormChi2:')
        for ii in var:
            if len(ii) in [1,2]:
                if len(ii) == 1:
                    jj = ii + ['FreeFormChi2_%s'%ii[0]]
                else:
                    jj = ii
                self._FreeFormChi2.append(jj)
                af.Info('    varID= %s\tName= %s'%(jj[0], jj[1]))

                self.Chi2[jj[1]] = af.NaN
            else:
                af.ErrorStop( 'The "FreeFormChi2" constraint on "%s" need 1 item or 2 items( VarID [, Name] ).'%(ii[0]) )

        self.Chi2 = af.sortDic(self.Chi2)

    def getChisq(self,par):
        chisq = 0.0

        ## add for "math ..." in [constrain]
        af.parseMath(par)

        for ii in self._Gaussian:
            
            if not af.is_number(par[ii[0]]):
              return af.log_zero
              
            try:
                ii1 = eval(ii[1])
            except:
                ii1 = ii[1]
            try:
                ii2 = eval(ii[2])
            except:
                ii2 = ii[2]
            if ii[3].lower() == 'symm':
                ichisq = (ii1 - par[ii[0]])**2/ii2**2
            elif ii[3].lower() == 'upper':
                if par[ii[0]] <= ii1 :
                    ichisq = 0
                else:
                    ichisq = (ii1 - par[ii[0]])**2/ii2**2
            elif ii[3].lower() == 'lower':
                if par[ii[0]] >= ii1:
                    ichisq = 0
                else:
                    ichisq = (ii1 - par[ii[0]])**2/ii2**2

            ## new 20180428 liang
            self.Chi2[ii[4]]=ichisq
            af.Debug("Chi2: %s = %f"%(ii[4], ichisq))

            chisq += ichisq

        ## marked 20180430 liang
        #for ii in self._Limit:
        #    af.ErrorStop('Line limit constraint is not ready.')
        for ii in self._FreeFormChi2:
            ichisq = par[ii[0]]

            self.Chi2[ii[1]]=ichisq
            af.Debug("Chi2: %s = %f"%(ii[1], ichisq))

            chisq += ichisq

        ## new 20180428 liang
        self.Chi2['Chi2'] = chisq

        return chisq
