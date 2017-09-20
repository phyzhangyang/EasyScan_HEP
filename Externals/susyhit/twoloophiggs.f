c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  subroutines for the dominant two-loop radiative corrections to the
c  Higgs masses, consitently in DRbar scheme for SuSpect >= 2.3 interface
c  Routines written by P. Slavich
c  last modif: 05/12/2005 (J-L Kneur) 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      subroutine su_DSZHiggs(t,mg,T1,T2,st,ct,q,mu,tanb,v2,gs,
     $     OS,S11,S22,S12)

c     Two-loop O(a_t a_s) corrections to the CP-even Higgs mass matrix. 
c     Routine written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on G. Degrassi, P. Slavich and F. Zwirner, 
c     Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
c
c     Last update:  24/02/2004: mg is given as input instead of mg^2;
c                               value of pi corrected (10th digit);
c                               unused variables cleaned up.
c                   22/10/2002: gs is given as input.
c     
c
c     I/O PARAMETERS:
c     t = m_top^2, mg = m_gluino, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tanb = tan(beta), v2 = v^2, 
c     gs = strong coupling constant,
c     OS = renormalization scheme for 1-loop (0 = DRbar, 1 = On-Shell),
c     Sij = 2-loop corrections to the CP-even Higgs mass matrix elements.

      implicit none

      integer OS
      real*8 ht,gs,k,mt,pi,v2
      real*8 t,mg,T1,T2,st,ct,q,A,X,mu,tanb,sb,s2t,c2t
      real*8 F1,F2,F3,sF2,sF3
      real*8 DF1,DF2,DF3,DsF2,DsF3
      real*8 F2_s,sF2_A,sF3_A
      real*8 S11,S22,S12,osdr

c$$$      pi = 3.14159265897d0
      pi = 3.1415926535898d0

      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt    ! eq. (19) of DSZ
      A = X - mu/tanb           ! notice the sign convention for mu

      sb = dsin(datan(tanb))
      ht = dsqrt(2d0/v2)*mt/sb
      
      k = 4d0*gs**2/(16d0*Pi**2)**2 ! gs^2/(16 Pi^2)^2 CF Nc
      
      call strfuncs(t,mg,T1,T2,s2t,c2t,q,F1,F2,F3)
      call strsfuncs(mg,T1,T2,q,A,sF2,sF3)
      call strdfuncs(t,mg,T1,T2,s2t,c2t,q,A,X,
     $     DF1,DF2,DF3,DsF2,DsF3)

      osdr = 1d0*OS

      if(s2t.ne.0.and.A.ne.0) then
         
         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + osdr*DF3) ! eq. (25)
         
         S12 = .5d0 * ht**2 * mu * A  * s2t**2 * (F3 + sF3 +   ! eq. (26)
     $        osdr*(DF3 + DsF3)) + 
     $        ht**2 * mt * mu * s2t * (F2 + osdr*DF2)
         
         S22 = .5d0 * ht**2 * A**2 * s2t**2 * (F3 + 2d0*sF3 +  ! eq. (27)
     $        osdr*(DF3 + 2d0*DsF3)) + 
     $        2d0 * ht**2 * mt * A * s2t * (F2 + sF2 + 
     $        osdr*(DF2 + DsF2)) + 
     $        2d0 * ht**2 * mt**2 * (F1 + osdr*DF1)

c     some of the functions have poles in s2t=0 or in A=0. 
c     when necessary we consider the residues:
         
      elseif(s2t.eq.0.and.A.eq.0) then
         
         S11 = 0d0
         S12 = 0d0
         S22 = 2 * ht**2 * mt**2 * (F1 + osdr*DF1)

      elseif(s2t.eq.0.and.A.ne.0) then 

         call strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)     

         S11 = 0d0
         S12 = ht**2 * mt * mu * (F2_s + osdr*DF2)
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * A * (F2_s + osdr*DF2)

      elseif(s2t.ne.0.and.A.eq.0) then

         call strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)     

         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + osdr*DF3)
         S12 = .5d0 * ht**2 * mu * s2t**2 * (sF3_A + osdr*DsF3) + 
     $        ht**2 * mt * mu * s2t * (F2 + osdr*DF2)
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * s2t * (sF2_A + osdr*DsF2)
 
      endif

      S11 = k*S11
      S12 = k*S12
      S22 = k*S22

      return
      end

*
***********************************************************************
*
      
      subroutine strfuncs(t,mg,T1,T2,s2t,c2t,q,F1,F2,F3)

      implicit none
      real*8 t,mg,T1,T2,s2t,c2t,q,F1,F2,F3
      real*8 strF1ab,strF1c,strF2ab,strF2c,strF3ab,strF3c

      F1 = strF1ab(t,T1,T2,s2t,c2t,q) 
     $     + strF1c(t,mg,T1,s2t,q)
     $     + strF1c(t,mg,T2,-s2t,q)
      
      F2 = strF2ab(T1,T2,s2t,c2t,q) 
     $     + strF2c(t,mg,T1,T2,s2t,q)
     $     - strF2c(t,mg,T2,T1,-s2t,q)

      F3 = strF3ab(T1,T2,s2t,c2t,q) 
     $     + strF3c(t,mg,T1,T2,s2t,q)
     $     + strF3c(t,mg,T2,T1,-s2t,q)

      return
      end

*
*********************************************************************
*
            
      function strF1ab(t,T1,T2,s2t,c2t,q)

      implicit none
      real*8 t,T1,T2,s2t,c2t,q
      real*8 strF1ab

      strF1ab =                    ! eq. (32) 
     $     -6*(1d0-dLog(t/q))+5*dLog(T1*T2/t**2)+dLog(T1*T2/t**2)**2
     $     +8*dLog(t/q)**2-4*dLog(T1/q)**2-4*dLog(T2/q)**2
     $     -c2t**2*(2d0-dLog(T1/q)-dLog(T2/q)-dLog(T1/T2)**2)
     $     -s2t**2*(T1/T2*(1d0-dLog(T1/q))+T2/T1*(1d0-dLog(T2/q)))
      
      return
      end
      
*
*********************************************************************
*

      function strF1c(t,mg,T1,s2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,s2t,q
      real*8 strF1c,phi,del

      mt = dsqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF1c =                     ! eq. (A1)
     $     +4*(t+g-mg*mt*s2t)/T1*(1d0-dLog(g/q))
     $     +4*dLog(t/g) - 2*dLog(T1/g)
     $     +2d0/del*(4*g**2*dLog(T1/g)
     $     +(g**2-T1**2+t*(10*g+3*t+2*t*g/T1-2*t**2/T1))*dLog(t/g))
     $     +2*mg/mt*s2t*(dLog(T1/q)**2+2*dLog(t/q)*dLog(T1/q))
     $     +4*mg/mt*s2t/del*(g*(T1-t-g)
     $     *Log(T1/g)+t*(T1-3*g-2*t-(t*g-t**2)/T1)*dLog(t/g))
     $     +(4*g*(t+g-T1-2*mg*mt*s2t)/del
     $     -4*mg/mt*s2t)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*

      function strF2ab(T1,T2,s2t,c2t,q)
      
      implicit none
      real*8 T1,T2,s2t,c2t,q
      real*8 strF2ab
      
      strF2ab =                    ! eq. (33)
     $     5*dLog(T1/T2)-3*(dLog(T1/q)**2-dLog(T2/q)**2)
     $     +c2t**2*(5*dLog(T1/T2)
     $     -(T1+T2)/(T1-T2)*dLog(T1/T2)**2
     $     -2/(T1-T2)*(T1*dLog(T1/q)-T2*dLog(T2/q))*dLog(T1/T2))
     $     +s2t**2*(T1/T2*(1d0-dLog(T1/q))-T2/T1*(1d0-dLog(T2/q)))
      
      return
      end

*
*********************************************************************
*

      function strF2c(t,mg,T1,T2,s2t,q)

      implicit none
      real*8 t,g,mg,mt,T1,T2,s2t,q
      real*8 strF2c,phi,del

      mt = dsqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF2c =                     ! eq. (A2)
     $     4*(t+g)/T1-4*mg/mt*s2t/(T1-T2)*(3*T1-t*T2/T1)
     $     +2*mg/mt*s2t/(T1-T2)*(
     $     (4*t+5*T1+T2)*dLog(T1/q)-2*t*T2/T1*dLog(g/q))
     $     -4*(g+t)/T1*dLog(g/q) - 2*dLog(T1/g) 
     $     +2.d0/del*(2*g*(g+t-T1)*dLog(T1/g)
     $     +2*t*(3*g+2*t-T1+(g*t-t**2)/T1)*dLog(t/g))
     $     -4*mg*mt*s2t/del/T1*(2*g*T1*dLog(T1/g)-
     $     ((t-T1)**2-g*(t+T1))*dLog(t/g))
     $     -8*mg*mt/s2t/(T1-T2)*(dLog(T1/q)-dLog(t/q)*dLog(T1/q))
     $     -mg/mt*s2t/(T1-T2)*((T1+T2)*dLog(T1/q)**2
     $     +(10*t-2*g+T1+T2)*dLog(t/q)*dLog(T1/q)
     $     +(2*g-2*t+T1+T2)*dLog(T1/q)*dLog(g/q))
     $     +(8*g*t/del-8*mg*mt/s2t/(T1-T2)
     $     +2*s2t/mg/mt/(T1-T2)*(4*g*t-del)
     $     +s2t/mg/mt/del*(T1-g-t)**3)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*
      
      function strF3ab(T1,T2,s2t,c2t,q)

      implicit none
      real*8 T1,T2,s2t,c2t,q
      real*8 strF3ab

      strF3ab =                    ! eq. (34)
     $     (3+9*c2t**2)*(2d0-(T1+T2)/(T1-T2)*dLog(T1/T2))
     $     +4d0-(3d0+13*c2t**2)/(T1-T2)*(T1*dLog(T1/q)-T2*dLog(T2/q))
     $     +3*(T1+T2)/(T1-T2)*(dLog(T1/q)**2-dLog(T2/q)**2)
     $     -c2t**2*(4d0-((T1+T2)/(T1-T2))**2*dLog(T1/T2)**2
     $     -6*(T1+T2)/(T1-T2)**2
     $     *(T1*dLog(T1/q)-T2*dLog(T2/q))*dLog(T1/T2))
     $     -s2t**2*(T1/T2+T2/T1 + 2*dLog(T1*T2/q**2)
     $     -T1**2/T2/(T1-T2)*dLog(T1/q)
     $     +T2**2/T1/(T1-T2)*dLog(T2/q))

      return
      end
      
*
*********************************************************************
*

      function strF3c(t,mg,T1,T2,s2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,q
      real*8 strF3c,phi,del

      mt = dsqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF3c =                     ! eq. (A3)
     $     -4*T2/T1/(T1-T2)*(g+t)
     $     +4*mg*mt*s2t/(T1-T2)**2*(21*T1-T2**2/T1)
     $     +4d0/(T1-T2)*(g*T2/T1*dLog(g/q)-2*(t+g)*dLog(T1/q))
     $     -24*mg*mt*s2t/(T1-T2)**2*(3*T1+T2)*dLog(T1/q)
     $     +4*t/T1/del*(2*g*T1*dLog(T1/q)-g*(g-t+T1)*dLog(g/q)+
     $     (g*(T+T1)-(t-T1)**2)*dLog(t/q))
     $     -4*mg*mt*s2t/T1/del*(t*(g-t+T1)*dLog(t/q)
     $     -g*(g-t-T1)*dLog(g/q)+T1*(g+t-T1)*dLog(T1/q))
     $     +2*(2*g+2*t-T1-T2)/(T1-T2)*dLog(g*t/q**2)*dLog(T1/q)
     $     +12*mg*mt*s2t/(T1-T2)**2*(2*(g-t)*dLog(g/t)*dLog(T1/q)
     $     +(T1+T2)*dLog(t*g/q**2)*dLog(T1/q))
     $     +8*mg*mt/s2t/(T1-T2)**2*
     $     (-8*T1+2*(3*T1+T2)*dLog(T1/q)-2*(g-t)*dLog(g/t)*dLog(T1/q)
     $     -(T1+T2)*dLog(t*g/q**2)*dLog(T1/q))
     $     -((8/s2t-12*s2t)*mt/mg/(T1-T2)**2
     $     *(2*del+(g+t-T1)*(T1-T2))
     $     +(4*del+8*g*t)/g/(T1-T2)+2*(g+t-T1)/g
     $     -4*t*(g+t-T1- 2*mg*mt*s2t)/del)*phi(t,T1,g)

      return
      end

*
*********************************************************************
*

      subroutine strsfuncs(mg,T1,T2,q,A,sF2,sF3)
      
c     shift to the Fi functions due to the renormalization of A
      
      implicit none
      real*8 mg,T1,T2,q,A,sF2,sF3
      
      sF2  = mg/A *       ! eq. (35)
     $     2d0*(dlog(T2/q)**2 - dlog(T1/q)**2)
      
      sF3  = mg/A *       ! eq. (36)
     $     (8d0 - 2d0*(T1+T2)/(T1-T2)*(dlog(T2/q)**2 - dlog(T1/q)**2)
     $     + 8d0/(T1-T2)*(T2*dlog(T2/q) - T1*dlog(T1/q)))
      
      return
      end

*
*********************************************************************
*

      subroutine strresfuncs(t,mg,T1,T2,q,F2_s,sF2_A,sF3_A)
      
c     residues of some singular functions for s2t=0 and for A=0
      
      implicit none      
      real*8 t,g,mt,mg,T1,T2,q,sF2_A,sF3_A,F2_s,phi
      
      mt = dsqrt(t)
      g = mg**2

      F2_s =  -8*mg*mt/(T1-T2)*(
     $     (dLog(T1/q)-dLog(t/q)*dLog(T1/q)+phi(t,T1,g))-
     $     (dLog(T2/q)-dLog(t/q)*dLog(T2/q)+phi(t,T2,g)))
      
      sF2_A  = mg*
     $     2d0*(dlog(T2/q)**2 - dlog(T1/q)**2)
      
      sF3_A  = mg*
     $     (8d0 - 2d0*(T1+T2)/(T1-T2)*(dlog(T2/q)**2 - dlog(T1/q)**2)
     $     + 8d0/(T1-T2)*(T2*dlog(T2/q) - T1*dlog(T1/q)))

      return
      end
      
*     
***********************************************************************
*

      subroutine strdfuncs(t,mg,T1,T2,s2t,c2t,q,At,X,
     $     DF1,DF2,DF3,DsF2,DsF3)
      
c     shift of the parameters from DRbar to On-Shell scheme
 
      implicit none      
      real*8 t,g,mt,mg,T1,T2,s2t,c2t,q,At,X,myB0
      real*8 DF1,DF2,DF3,DsF2,DsF3
      real*8 msdr
      real*8 F1o,F2o,F3o,dm1,dm2,dmt,dAt,dth,ds2t

      msdr = -5d0
      mt = dSqrt(t)
      g = mg**2

      F1o = dLog(T1/q) + dLog(T2/q) - 2d0*dLog(t/q) ! eq. (31)
      F2o = dLog(T1/q) - dLog(T2/q) 
      F3o = 2d0 - (T1+T2)/(T1-T2)*(dLog(T1/q) - dLog(T2/q))

      dmt =                     ! eq. (B2)
     $     mt*(3*Log(t/q) + msdr + .5d0*(2*g/t*(dLog(g/q)-1d0)
     $      -T1/t*(dLog(T1/q)-1d0) - T2/t*(dLog(T2/q)-1d0)
     $      +(g+t-T1 - 2*s2t*mg*mt)/t*myB0(t,g,T1,q)
     $      +(g+t-T2 + 2*s2t*mg*mt)/t*myB0(t,g,T2,q)))

      dm1 =                     ! eq. (B3)
     $     T1*(3*dLog(T1/q) - 7d0 - c2t**2*(dLog(T1/q)-1d0)
     $      -s2t**2*T2/T1*(dLog(T2/q)-1d0) + 2*(
     $      g/T1*(dLog(g/q)-1d0) + t/T1*(dLog(t/q)-1d0)
     $      +(T1-g-t + 2*s2t*mg*mt)/T1*myB0(T1,t,g,q)))
      
      dm2 =                     ! eq. (B4)
     $     T2*(3*dLog(T2/q) - 7d0 - c2t**2*(dLog(T2/q)-1d0)
     $      -s2t**2*T1/T2*(dLog(T1/q)-1d0) + 2*(
     $      g/T2*(dLog(g/q)-1d0) + t/T2*(dLog(t/q)-1d0)
     $      +(T2-g-t - 2*s2t*mg*mt)/T2*myB0(T2,t,g,q)))

c$$$c     On-Shell theta-stop: asymmetric definition used in FeynHiggs
c$$$      dth = (4d0*mg*mt*c2t*myB0(T1,t,g,q) +
c$$$     $     c2t*s2t*(T2*(1d0-Log(T2/q))-T1*(1d0-Log(T1/q))))/(T1-T2)      

c     On-Shell theta-stop: eq. (B6)-(B7) of DSZ 
      dth = (4d0*mg*mt*c2t*(myB0(T1,t,g,q)+myB0(T2,t,g,q)) +
     $     2d0*c2t*s2t*(T2*(1d0-dLog(T2/q))-T1*(1d0-dLog(T1/q))))/
     $     2d0/(T1-T2)      

      ds2t = 2d0*c2t*dth

      dAt = ((dm1-dm2)/(T1-T2) + ds2t/s2t - dmt/mt)*X ! eq. (B8)

      DF1 = dm1/T1 + dm2/T2 - 4d0*dmt/mt + 4d0*dmt/mt*F1o ! eq. (37)
      DF2 = dm1/T1 - dm2/T2 + (3d0*dmt/mt + ds2t/s2t)*F2o ! eq. (38)
      DF3 = (2d0*T1*T2/(T1-T2)**2*dLog(T1/T2) - (T1+T2)/(T1-T2)) 
     $     *(dm1/T1-dm2/T2) + (2d0*dmt/mt + 2d0*ds2t/s2t)*F3o ! eq. (39)

      DsF2 = dAt/At * F2o       ! eq. (40)
      DsF3 = dAt/At * F3o       ! eq. (41)

c     residues of some singular functions for s2t=0 and for A=0

      if(s2t.eq.0d0) then
         DF2  = ds2t*F2o
         DsF2 = ds2t*X/At * F2o
      endif

      if(At.eq.0d0) then
         DsF2 = dAt * F2o
         DsF3 = dAt * F3o
      endif

      return
      end

*
**************************************************************************
**************************************************************************
*

      subroutine su_DSZodd(t,mg,T1,T2,st,ct,q,mu,tanb,v2,gs,DMA)

c     Two-loop O(a_t a_s) corrections to the CP-odd Higgs mass in the
c     DRbar scheme. Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on G. Degrassi, P. Slavich and F. Zwirner,
c     Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
c
c     Last update:  24/02/2004: mg is given as input instead of mg^2;
c                               value of pi corrected (10th digit);
c                               unused variables cleaned up;
c                               some formulae simplified.
c                   22/10/2002: gs is given as input.
c
c
c     I/O PARAMETERS:
c     t = m_top^2, mg = m_gluino, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tanb = tan(beta), v2 = v^2,
c     gs = strong coupling constant,
c     DMA = 2-loop corrections to the CP-odd Higgs mass.

      implicit none

      real*8 ht,gs,k,mt,pi,v2
      real*8 t,mg,T1,T2,st,ct,q,A,X,mu,tanb,sb,cb,s2t
      real*8 FA,FA_A,DMA

c$$$      pi = 3.14159265897d0
      pi = 3.1415926535898d0

      mt = dsqrt(t)

      s2t = 2d0*ct*st

      X = (T1-T2)*s2t/2d0/mt    ! eq. (19) of DSZ
      A = X - mu/tanb           ! notice the sign convention for mu
      
      sb = dsin(datan(tanb))
      cb = dcos(datan(tanb))

      ht = dsqrt(2d0/v2)*mt/sb
     
      k = 4d0*gs**2/(16d0*Pi**2)**2 ! gs^2/(16 Pi^2)^2 CF Nc

      call strfuncsodd(t,mg,T1,T2,s2t,q,A,FA,FA_A)

      if(A.ne.0d0) then

         DMA = ht**2*mu*A/(T1-T2)/sb/cb * FA ! eq. (C1)

      else
         
c     the function FA has poles in A=0. 
c     when necessary we consider the residues:
         
         DMA = ht**2*mu/(T1-T2)/sb/cb  * FA_A

      endif
      
      DMA = k*DMA

      return
      end

*
***********************************************************************
*
      
      subroutine strfuncsodd(t,mg,T1,T2,s2t,q,A,FA,FA_A)

      implicit none
      real*8 t,mg,T1,T2,s2t,q,A,FA,FA_A
      real*8 strFAab,strFAc,strresFAc

      FA = strFAab(T1,T2,s2t,q) 
     $     + strFAc(t,mg,T1,T2,s2t,A,q)
     $     - strFAc(t,mg,T2,T1,-s2t,A,q)

      FA_A = strresFAc(t,mg,T1,q) - strresFAc(t,mg,T2,q)
      
      return
      end

*
*********************************************************************
*
            
      function strFAab(T1,T2,s2t,q)

      implicit none
      real*8 T1,T2,s2t,q
      real*8 strFAab

      strFAab =                    ! eq. (C4)
     $     (8d0 - s2t**2*(2d0 - (T1+T2)/(T1-T2)*dLog(T1/T2)))
     $     *(T1*(1d0 - dLog(T1/q)) - T2*(1d0 - dLog(T2/q)))
     $     + 2d0*(T1*dLog(T1/q)**2 - T2*dLog(T2/q)**2)
     $     + 2d0/(T1-T2)*(T1*dLog(T1/q) - T2*dLog(T2/q))**2
      
      return
      end
      
*
*********************************************************************
*

      function strFAc(t,mg,T1,T2,s2t,A,q)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,A,q
      real*8 strFAc,phi,del

      mt = dSqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strFAc =                     ! eq. (C5)
     $     16d0*T1/(T1-T2)*mg*mt*s2t  
     $     -4d0*(g+t)*dLog(T1/q)
     $     - 2d0*mg/A*T1*(5d0- 4d0*dLog(T1/q))
     $     - 4d0*mg*mt*s2t*(3d0* T1 + T2)/(T1-T2)*dLog(T1/q)
     $     + 2d0*(g+t-T1)*dLog(t*g/q**2)*dLog(T1/q)
     $     + 2d0*T1*(1d0+mg/A)*dLog(g/q)*dLog(t/q) 
     $     - 2d0*mg/A*((g-t)*dLog(g/t)+ T1*dLog(t*g/q**2))*dLog(T1/q)
     $     + 2d0*mg*mt*s2t/(T1-T2)*
     $     (2d0*(g-t)*dLog(g/t) + (T1+T2)*dLog(t*g/q**2))*dLog(T1/q)
     $     - (4d0* t + 2d0* del/g*(1d0 + mg/A) 
     $     - 2d0*mt/mg*s2t/(T1-T2)*
     $     (2d0*del + (g+t-T1)*(T1-T2)))*phi(t,T1,g)

      return
      end

*
*********************************************************************
*

      function strresFAc(t,mg,T1,q)

c     residue of FA for A=0

      implicit none
      real*8 t,g,mg,T1,q
      real*8 strresFAc,phi,del

      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strresFAc =                 
     $     - 2d0*mg*T1*(5d0- 4d0*dLog(T1/q))
     $     + 2d0*mg*T1*dLog(g/q)*dLog(t/q) 
     $     - 2d0*mg*((g-t)*dLog(g/t)+ T1*dLog(t*g/q**2))*dLog(T1/q)
     $     - 2*del/mg*phi(t, T1, g)
      
      return
      end

*
*********************************************************************
*********************************************************************
*

      subroutine su_ewsb2loop(t,mg,T1,T2,st,ct,q,mu,tanb,vv,gs,
     $     S1,S2)

c     Two-loop O(a_t a_s) corrections to the Higgs tadpoles. 
c     Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on A. Dedes and P. Slavich, 
c     Nucl. Phys. B657 (2003) 333 [hep-ph/0212132].
c
c     Last update:  24/02/2004: mg is given as input instead of mg^2;
c                               value of pi corrected (10th digit);
c                               unused variables cleaned up.
c     Last update:  12/12/2002.
c
c     I/O PARAMETERS:
c     t = m_top^2, g = m_gluino^2, T1 = m_stop1^2, T2 = m_stop2^2,
c     st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tanb = tan(beta), vv = v^2,
c     gs = strong coupling constant
c     Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles.
c
c     Notice: we assume that the 1-loop part is computed in terms of 
c             running (DRbar) parameters, evaluated at the scale Q.

      implicit none

      real*8 gs,k,mt,vv
      real*8 t,mg,T1,T2,st,ct,q,A,X,mu,tanb,sb,cb,s2t,c2t,v1,v2
      real*8 F2l,G2l,S1,S2,pi

c$$$      pi = 3.14159265897d0
      pi = 3.1415926535898d0

      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt
      A = X - mu/tanb           ! notice the sign convention for mu

      sb = dsin(datan(tanb))
      cb = dcos(datan(tanb))
      v1 = sqrt(vv)*cb
      v2 = sqrt(vv)*sb

      k = 4d0*gs**2/(16d0*Pi**2)**2 ! gs^2/(16 Pi^2)^2 CF Nc
      
      call strfuncstad(t,mg,T1,T2,s2t,c2t,q,F2l,G2l)

      S1 = mt * mu/tanB * s2t * F2l
      S1 = S1/v1**2

      S2 = mt * A * s2t * F2l + 2d0 * mt**2 * G2l 
      S2 = S2/v2**2

      S1 = k*S1
      S2 = k*S2

      return
      end

*
***********************************************************************
*

      
      subroutine strfuncstad(t,mg,T1,T2,s2t,c2t,q,F2l,G2l)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,c2t,q
      real*8 F2l,G2l,strF2lc,strG2lc

      mt = sqrt(t)
      g = mg**2

      F2l = 4*mg*mt*(1+4*c2t**2)/s2t
     $     -(2*(T1-T2)+4*mg*mt/s2t)*dLog(g/q)*dLog(t/q)
     $     -2*(4-s2t**2)*(T1-T2)
     $     +(4*T1*T2-s2t**2*(T1+T2)**2)/(T1-T2)*dLog(T1/q)*dLog(T2/q)
     $     + strF2lc(t,mg,T1,T2,s2t,c2t,q)
     $     - strF2lc(t,mg,T2,T1,-s2t,c2t,q)

      G2l = 5*mg/mt*s2t*(T1-T2)-10*(T1+T2-2*t)-4*g
     $     + 12*t*(dLog(t/q)**2-2*dLog(t/q))
     $     +(4*g-s2t*mg/mt*(T1-T2))*dLog(g/q)*dLog(t/q)
     $     +s2t**2*(T1+T2)*dLog(T1/q)*dLog(T2/q)
     $     + strG2lc(t,mg,T1,T2,s2t,q)
     $     + strG2lc(t,mg,T2,T1,-s2t,q)
      
      return
      end

*
*********************************************************************
*

      function strF2lc(t,mg,T1,T2,s2t,c2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,c2t,q
      real*8 strF2lc,phi,del

      mt = dsqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strF2lc = (4*(g+t+2*T1)-s2t**2*(3*T1+T2)-4*s2t*mg*mt
     $     -16*mg*mt*T1*c2t**2/s2t/(T1-T2))*dLog(T1/q)       
     $     +T1/(T1-T2)*(s2t**2*(T1+T2)-2*(2*T1-T2))*dLog(T1/q)**2
     $     +2*(T1-g-t+mg*mt*s2t
     $     +2*c2t**2*mg*mt*T1/s2t/(T1-T2))*dLog(g*t/q**2)*dLog(T1/q)
     $     +4*mg*mt*c2t**2*(t-g)/s2t/(T1-T2)*dLog(t/g)*dLog(T1/q)
     $     +((2*del+4*g*t)/T1-2*mg*mt*s2t/T1*(g+t-T1)
     $     +4*c2t**2*mg*mt/T1/(T1-T2)/s2t*del)*phi(g,t,T1)

      return
      end


*
*********************************************************************
*

      function strG2lc(t,mg,T1,T2,s2t,q)

      implicit none
      real*8 t,g,mt,mg,T1,T2,s2t,q
      real*8 strG2lc,phi,del

      mt = sqrt(t)
      g = mg**2

      del = g**2 + t**2 + T1**2 - 2*(g*t + g*T1 + t*T1)
      
      strG2lc = (4*(g+t+2*T1)+s2t**2*(T1-T2)
     $     -4*mg/mt*s2t*(t+T1))*dLog(T1/q) 
     $     +(mg/mt*s2t*(5*t-g+T1)-2*(g+2*t))*dLog(t/q)*dLog(T1/q)
     $     +(mg/mt*s2t*(g-t+T1)-2*g)*dLog(g/q)*dLog(T1/q)
     $     -(2+s2t**2)*T1*dLog(T1/q)**2
     $     +(2*g/T1*(g+t-T1-2*mg*mt*s2t)
     $     +mg/mt*s2t*del/T1)*phi(g,t,T1)

      return
      end

*
**********************************************************************
**********************************************************************
*

c
c     Two-loop O(a_t^2 + at ab + ab^2) corrections to the Higgs masses 
c     and to the minimization conditions of the effective potential. 
c     Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127.
c
c     Last update:  25/09/2003: numerical value of Pi corrected (10th digit)
c                   17/09/2003: error corrected in F1q, F6q and FAq;
c                               optimization: calls to makederiv reduced to 2
c                   20/07/2003: routines optimized, bug for mu<0 corrected
c
c     I/O PARAMETERS:
c     t = m_top^2, b = m_bot^2, A0 = m_A^2, T1 = m_stop1^2, T2 = m_stop2^2,
c     B1 = m_sbot1^2, B2 = m_sbot2^2, st = sin(theta_stop), 
c     ct = cos(theta_stop), sb = sin(theta_sbot), cb = cos(theta_sbot),
c     q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanb = tan(beta), 
c     vv = v^2,
c     Sij = 2-loop corrections to the CP-even Higgs mass matrix elements,
c     Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles,
c     DMA = 2-loop corrections to the CP-odd Higgs mass.
c
c     Notice: we assume that the 1-loop part is computed in terms of 
c             running (DRbar) parameters, evaluated at the scale Q. The 
c             parameters in the bottom/sbottom sector should be computed 
c             in term of the "resummed" bottom Yukawa coupling.
c
c

      subroutine su_DDSHiggs(t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,
     $     vv,S11,S12,S22) 

      implicit none

      real*8 t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,S11,S12,S22
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,ht,hb,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp

      pi = 3.14159265358979d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu/tanb           
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = dsqrt(2d0/vv)*mt/sbe
      hb = dsqrt(2d0/vv)*mb/cbe

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncs(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      S11 = .5d0*ht**2*mu**2*s2t**2*F3t
     $     + 2d0*hb**2*mb**2*F1b + 2d0*hb**2*Ab*mb*s2b*F2b
     $     + .5d0*hb**2*Ab**2*s2b**2*F3b
     $     + 2d0*hb*ht*mb*mu*s2t*F4t + ht*hb*mu*Ab*s2t*s2b*F5

      S12 = ht**2*mu*mt*s2t*F2t + .5d0*ht**2*At*mu*s2t**2*F3t
     $     +hb**2*mu*mb*s2b*F2b + .5d0*hb**2*Ab*mu*s2b**2*F3b
     $     +ht*hb*mb*At*s2t*F4t + hb*ht*mt*Ab*s2b*F4b
     $     +.5d0*ht*hb*s2t*s2b*(At*Ab+mu**2)*F5
     $     +2d0*ht*hb*mt*mb*F6

      S22 = .5d0*hb**2*mu**2*s2b**2*F3b
     $     + 2d0*ht**2*mt**2*F1t + 2d0*ht**2*At*mt*s2t*F2t
     $     + .5d0*ht**2*At**2*s2t**2*F3t
     $     + 2d0*ht*hb*mt*mu*s2b*F4b + hb*ht*mu*At*s2b*s2t*F5

      S11 = k*S11
      S12 = k*S12
      S22 = k*S22

      return
      end


*
***********************************************************************
*
      
      subroutine su_DDSodd(t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,
     $     DMA) 

      implicit none

      real*8 t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,DMA
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,ht,hb,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp

      pi = 3.14159265358979d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu/tanb           
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = dsqrt(2d0/vv)*mt/sbe
      hb = dsqrt(2d0/vv)*mb/cbe

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncs(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      DMA = -(ht**2*mu*At/(T1-T2)*Ft + hb**2*mu*Ab/(B1-B2)*Fb
     $     + 2d0*ht*hb*FAp)/sbe/cbe

      DMA = k*DMA

      return
      end


*
***********************************************************************
*
      
      subroutine su_DDStad(t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,
     $     S1,S2) 

      implicit none

      real*8 t,b,A0,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp
      real*8 v1,v2,S1,S2

      pi = 3.14159265358979d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu/tanb           
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncs(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      v1 = Sqrt(vv)*cbe
      v2 = Sqrt(vv)*sbe

      S1 = mt*mu/tanb*s2t*Ft + mb*Ab*s2b*Fb + 2d0*mb**2*Gb
      S2 = mb*mu*tanb*s2b*Fb + mt*At*s2t*Ft + 2d0*mt**2*Gt
            
      S1 = k*S1/v1**2
      S2 = k*S2/v2**2

      return
      end

*
***********************************************************************
*
      
      subroutine makefuncs(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     $     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
     $     Ft,Fb,Gt,Gb,FAp)

      implicit none

      real*8 t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp
      
      real*8 D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     $     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     $     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     $     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     $     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      common/listderiv/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     $     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     $     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     $     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     $     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     $     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb

      real*8 D1b,DB1,DB2,Dc2b,DB1B1,DB2B2,Dbb,Dc2bc2b,DB1b,DB2b,DB1B2,
     $     Dbc2b,DB1c2b,DB2c2b,Dtc2b

      real*8 Xt,Xb,At,Ab

      call makederiv(b,t,A0,B1,B2,T1,T2,s2b,c2b,s2t,c2t,
     $     q,mu,vv,1d0/tanb)

      D1b = D1t
      DB1 = DT1
      DB2 = DT2
      Dc2b = Dc2t
      DB1B1 = DT1T1
      DB2B2 = DT2T2
      Dbb = Dtt
      Dc2bc2b = Dc2tc2t
      DB1b = DT1t
      DB2b = DT2t
      DB1B2 = DT1T2
      Dbc2b = Dtc2t
      DB1c2b = DT1c2t
      DB2c2b = DT2c2t
      Dtc2b = Dbc2t

      call makederiv(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)

      F1t = Dtt + DT1T1 + DT2T2 + 2d0*(DT1t + DT2t + DT1T2)
     $     +(Dcptpb + Dcptmptt + Dcptptb - 2d0*Dsptmpttsptptb)
     $     /4d0/t**2

      F2t = DT1T1 - DT2T2 + DT1t - DT2t
     $     -4d0*c2t**2/(T1-T2)*(Dtc2t + DT1c2t + DT2c2t)
     $     -(Dcptmptt - Dsptmpttsptptb)/s2t**2/t/(T1-T2)
      
      F3t = DT1T1 + DT2T2 - 2d0*DT1T2
     $     - 2d0/(T1-T2)*(DT1-DT2)
     $     + 16d0*c2t**2/(T1-T2)**2*(c2t**2*Dc2tc2t + 2d0*Dc2t)
     $     -8d0*c2t**2/(T1-T2)*(DT1c2t-DT2c2t)
     $     + 4d0/s2t**4/(T1-T2)**2*(Dcptmptt + Dcpbptt + Dcpttptb)

      F4t = DT1b + DT1B1 + DT1B2 - DT2b - DT2B1 - DT2B2
     $     -4d0*c2t**2/(T1-T2)*(DB1c2t + DB2c2t + Dbc2t)
     $     -(Dcpbptt + Dsptmpttspbmptb - Dspbmptbspbptt)
     $     /b/s2t**2/(T1-T2)

      F1b = Dbb + DB1B1 + DB2B2 + 2d0*(DB1b + DB2b + DB1B2)
     $     +(Dcptpb + Dcpbmptb + Dcpbptt - 2d0*Dspbmptbspbptt)
     $     /4d0/b**2

      F2b = DB1B1 - DB2B2 + DB1b - DB2b
     $     -4d0*c2b**2/(B1-B2)*(Dbc2b + DB1c2b + DB2c2b)
     $     -(Dcpbmptb - Dspbmptbspbptt)/s2b**2/b/(B1-B2)
      
      F3b = DB1B1 + DB2B2 - 2d0*DB1B2
     $     - 2d0/(B1-B2)*(DB1-DB2)
     $     + 16d0*c2b**2/(B1-B2)**2*(c2b**2*Dc2bc2b + 2d0*Dc2b)
     $     -8d0*c2b**2/(B1-B2)*(DB1c2b-DB2c2b)
     $     + 4d0/s2b**4/(B1-B2)**2*(Dcpbmptb + Dcptptb + Dcpttptb)

      F4b = DB1t + DT1B1 - DT1B2 - DB2t + DT2B1 - DT2B2
     $     -4d0*c2b**2/(B1-B2)*(DT1c2b + DT2c2b + Dtc2b)
     $     -(Dcptptb + Dsptmpttspbmptb - Dsptmpttsptptb)
     $     /t/s2b**2/(B1-B2)

      F5  = DT1B1 - DT1B2 - DT2B1 + DT2B2
     $     + 16d0*c2t**2*c2b**2/(T1-T2)/(B1-B2)*Dc2tc2b
     $     - 4d0*c2t**2/(T1-T2)*(DB1c2t-DB2c2t)
     $     - 4d0*c2b**2/(B1-B2)*(DT1c2b-DT2c2b)
     $     - 4d0/s2b**2/s2t**2/(T1-T2)/(B1-B2)*
     $     (Dcpttptb-Dsptmpttspbmptb + Dspbmptbspbptt + Dsptmpttsptptb)
      
      F6 = Dtb + DT1b + DT2b + DB1t + DB2t
     $     + DT1B1 + DT1B2 + DT2B1 + DT2B2
     $     -(Dcptpb - Dsptmpttspbmptb)/4d0/b/t 
      
      Ft = DT1 - DT2 - 4d0*c2t**2/(T1-T2)*Dc2t      
      
      Gt = D1t + DT1 + DT2

      Fb = DB1 - DB2 - 4d0*c2b**2/(B1-B2)*Dc2b
      
      Gb = D1b + DB1 + DB2

      Xt = (T1-T2)*s2t/2d0/sqrt(t)    
      Xb = (B1-B2)*s2b/2d0/sqrt(b)

      At = Xt - mu/tanb
      Ab = Xb - mu*tanb

      FAp = Dcptpb/4d0/Sqrt(b)/Sqrt(t)
     $     +4d0*(At*Ab - mu**2)**2*Sqrt(t)*Sqrt(b)
     $     /s2t**2/s2b**2/(T1-T2)**2/(B1-B2)**2*Dcpttptb
     $     +Sqrt(t)/Sqrt(b)/s2t**2/(T1-T2)**2
     $     *(At**2*Dcpbptt + mu**2/tanb**2*Dcptmptt)
     $     +Sqrt(b)/Sqrt(t)/s2b**2/(B1-B2)**2
     $     *(Ab**2*Dcptptb + mu**2*tanb**2*Dcpbmptb)
     $     -2d0*mu/s2t/s2b/(T1-T2)/(B1-B2)
     $     *(At*tanb*Dspbmptbspbptt + Ab/tanb*Dsptmpttsptptb
     $     + mu*Dsptmpttspbmptb)
      
      return
      end

*
***********************************************************************
*

      subroutine makederiv(t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     $     q,mu,vv,tanb)
      
      implicit real*8 (t)
      implicit character (a-s,u-z)

      real*8 t,b,A0,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,mt,mb,Xt,Xb,Yt,Yb,sbe,cbe,mu2,Nc
      real*8 Delt,phi,Li2
      external delt,phi,Li2

      real*8 D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     $     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     $     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     $     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     $     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb
      
      common/listderiv/D1t,DT1,DT2,Dc2t,DT1T1,DT2T2,
     $     Dtt,Dc2tc2t,DT1t,DT2t,DT1T2,
     $     Dtc2t,DT1c2t,DT2c2t,Dtb,DT1b,DT2b,DB1t,DB2t,DT1B1,DT2B1,
     $     DT1B2,DT2B2,Dbc2t,DB1c2t,DB2c2t,DT1c2b,DT2c2b,Dc2tc2b,
     $     Dcptpb,Dcpttptb,Dcpbptt,Dcptptb,Dcptmptt,Dcpbmptb,
     $     Dspbmptbspbptt,Dsptmpttsptptb,Dsptmpttspbmptb
      
      real*8  Logt,Logb,Logmu2,LogA0,LogT1,LogT2,LogB1,LogB2,
     $     phimu2tT1,phimu2tT2,phiB1tmu2,phiB2tmu2,phiT1bmu2,phiT2bmu2,
     $     phiA0T1T1,phiA0T2T2,phiA0T1T2,phiA0B1B1,phiA0B2B2,
     $     phiA0B1T1,phiA0B2T1,phiA0B1T2,phiA0B2T2,phiA0tt,phiA0bt,
     $     deltmu2tT1,deltmu2tT2,deltB1tmu2,deltB2tmu2,
     $     deltT1bmu2,deltT2bmu2,deltA0tt,deltA0bt,
     $     deltA0T1T1,deltA0T2T2,deltA0T1T2,deltA0B1B1,deltA0B2B2,
     $     deltA0B1T1,deltA0B2T1,deltA0B1T2,deltA0B2T2,
     $     Li2T1T2,Li2T1B1,Li2T1B2,Li2T2B1,Li2T2B2,Li2bt

      Nc = 3d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))
      
      ht = dsqrt(2d0/vv)*mt/sbe
      hb = dsqrt(2d0/vv)*mb/cbe
            
      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb           
      Yt  = Xt - mu/sbe/cbe
      Yb  = Xb - mu/sbe/cbe
      
      mu2 = mu**2

      Logt = dLog(t/q)
      Logb = dLog(b/q)
      Logmu2 = dLog(mu2/q)
      LogA0 = dLog(A0/q)
      LogT1 = dLog(T1/q)
      LogT2 = dLog(T2/q)
      LogB1 = dLog(B1/q)
      LogB2 = dLog(B2/q)
      phimu2tT1 = phi(mu2,t,T1)
      phimu2tT2 = phi(mu2,t,T2)
      phiB1tmu2 = phi(B1,t,mu2)
      phiB2tmu2 = phi(B2,t,mu2)
      phiT1bmu2 = phi(T1,b,mu2)
      phiT2bmu2 = phi(T2,b,mu2)
      phiA0T1T1 = phi(A0,T1,T1)
      phiA0T1T2 = phi(A0,T1,T2)
      phiA0T2T2 = phi(A0,T2,T2)
      phiA0B1B1 = phi(A0,B1,B1)
      phiA0B2B2 = phi(A0,B2,B2)
      phiA0B1T1 = phi(A0,B1,T1)
      phiA0B1T2 = phi(A0,B1,T2)
      phiA0B2T1 = phi(A0,B2,T1)
      phiA0B2T2 = phi(A0,B2,T2)
      phiA0tt = phi(A0,t,t)
      phiA0bt = phi(A0,b,t)
      deltmu2tT1 = delt(mu2,t,T1)
      deltmu2tT2 = delt(mu2,t,T2)
      deltB1tmu2 = delt(B1,t,mu2)
      deltB2tmu2 = delt(B2,t,mu2)
      deltT1bmu2 = delt(T1,b,mu2)
      deltT2bmu2 = delt(T2,b,mu2)
      deltA0T1T1 = delt(A0,T1,T1)
      deltA0T1T2 = delt(A0,T1,T2)
      deltA0T2T2 = delt(A0,T2,T2)
      deltA0B1B1 = delt(A0,B1,B1)
      deltA0B2B2 = delt(A0,B2,B2)
      deltA0B1T1 = delt(A0,B1,T1)
      deltA0B1T2 = delt(A0,B1,T2)
      deltA0B2T1 = delt(A0,B2,T1)
      deltA0B2T2 = delt(A0,B2,T2)
      deltA0tt = delt(A0,t,t)
      deltA0bt = delt(A0,b,t)
      Li2T1T2 = Li2(1d0-T1/T2)
      Li2T1B1 = Li2(1d0-T1/B1)
      Li2T2B1 = Li2(1d0-T2/B1)
      Li2T1B2 = Li2(1d0-T1/B2)
      Li2T2B2 = Li2(1d0-T2/B2)
      Li2bt = Li2(1d0-b/t)

	tmp1 = 2d0 - LogA0 - Logb + 2*Logt
	tmp2 = 2d0 - LogB1 - Logmu2 + 2*Logt
	tmp3 = 2d0 - LogB2 - Logmu2 + 2*Logt
	tmp4 = 2d0 - Logb - Logmu2 + 2*LogT1
	tmp5 = 2d0 - Logb - Logmu2 + 2*LogT2
	tmp6 = 0.25D0*hb**2/c2t - 0.25D0*ht**2/c2t
	tmp7 = -(0.25D0*hb**2/c2t) + 0.25D0*ht**2/c2t
	tmp8 = c2t**2 - 0.5D0*((-1d0 + Nc)*s2t**2)
	tmp9 = cbe**2*ht**2 + hb**2*sbe**2
	tmp10 = cbe**2*hb**2 + ht**2*sbe**2
        tmp11 = -(0.5D0*(cbe**2*ht**2*mb*mt*s2b*s2t)) - 
     -   0.5D0*(hb**2*mb*mt*s2b*s2t*sbe**2)
        tmp12 = 0.5D0*(cbe**2*ht**2*mb*mt*s2b*s2t) + 
     -   0.5D0*(hb**2*mb*mt*s2b*s2t*sbe**2)
        tmp13 = -(0.5D0*(cbe**2*hb**2*mb*mt*s2b*s2t)) - 
     -   0.5D0*(ht**2*mb*mt*s2b*s2t*sbe**2)
        tmp14 = 0.5D0*(cbe**2*hb**2*mb*mt*s2b*s2t) + 
     -   0.5D0*(ht**2*mb*mt*s2b*s2t*sbe**2)
	tmp15 = 2.d0/t**2 - (2*(-1d0 + Logt))/t**2
	tmp16 = 2.d0/t**2 - (2*Logt)/t**2
	tmp17 = 2.d0/T1**2 - (2*LogT1)/T1**2
	tmp18 = 2.d0/T2**2 - (2*LogT2)/T2**2
        tmp19 = -(0.0625D0*
     -      (B1*(-1d0 + LogB1)*(-1d0 + LogT1)*T1)/(c2b*c2t)) + 
     -   0.0625D0*(B2*(-1d0 + LogB2)*(-1d0 + LogT1)*T1)/(c2b*c2t) + 
     -   0.0625D0*(B1*(-1d0 + LogB1)*(-1d0 + LogT2)*T2)/(c2b*c2t) - 
     -   0.0625D0*(B2*(-1d0 + LogB2)*(-1d0 + LogT2)*T2)/(c2b*c2t)
        tmp20 = 0.25D0*((1d0 + c2b)*(1d0 + c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1d0 + c2b)*s2t*Xb)/mt
        tmp21 = 0.0625D0*b/(c2b*c2t) + 0.125D0*(mb*mt)/(s2b*s2t) + 
     -   0.0625D0*t/(c2b*c2t) - 0.125D0*(mb*Xb)/(c2t*s2b) + 
     -   0.125D0*(mt*Xb)/(c2b*s2t) - 0.0625D0*Xb**2/(c2b*c2t)
        tmp22 = -(0.0625D0*b/(c2b*c2t)) - 
     -   0.125D0*(mb*mt)/(s2b*s2t) - 0.0625D0*t/(c2b*c2t) + 
     -   0.125D0*(mb*Xb)/(c2t*s2b) - 0.125D0*(mt*Xb)/(c2b*s2t) + 
     -   0.0625D0*Xb**2/(c2b*c2t)
	tmp23 = 0.25D0*Xb/c2b - 0.25D0*Xt/c2b
	tmp24 = -(0.25D0*Xb/c2b) + 0.25D0*Xt/c2b
	tmp25 = 0.125D0*Xb/c2t**3 - 0.125D0*Xt/c2t**3
	tmp26 = -(0.125D0*Xb/c2t**3) + 0.125D0*Xt/c2t**3
	tmp27 = 0.25D0*Xb/c2t - 0.25D0*Xt/c2t
	tmp28 = -(0.25D0*Xb/c2t) + 0.25D0*Xt/c2t
        tmp29 = 0.25D0*((1d0 + c2b)*(1d0 + c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1d0 + c2b)*s2t*Xt)/mt
	tmp30 = (-2*mt*Xt)/s2t - Xt**2
	tmp31 = (2*mt*Xt)/s2t - Xt**2
        tmp32 = 0.0625D0*b/(c2b*c2t) + 0.125D0*(mb*mt)/(s2b*s2t) + 
     -   0.0625D0*t/(c2b*c2t) + 0.125D0*(mb*Xt)/(c2t*s2b) - 
     -   0.125D0*(mt*Xt)/(c2b*s2t) - 0.0625D0*Xt**2/(c2b*c2t)
        tmp33 = -(0.0625D0*b/(c2b*c2t)) - 
     -   0.125D0*(mb*mt)/(s2b*s2t) - 0.0625D0*t/(c2b*c2t) - 
     -   0.125D0*(mb*Xt)/(c2t*s2b) + 0.125D0*(mt*Xt)/(c2b*s2t) + 
     -   0.0625D0*Xt**2/(c2b*c2t)
	tmp34 = 4*t - 4*mt*s2t*Xt + s2t**2*Xt**2
	tmp35 = 4*t + 4*mt*s2t*Xt + s2t**2*Xt**2
        tmp36 = 0.25D0*((1d0 + c2b)*(1d0 + c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1d0 + c2b)*s2t*Yb)/mt
        tmp37 = 0.0625D0*b/(c2b*c2t) + 0.125D0*(mb*mt)/(s2b*s2t) + 
     -   0.0625D0*t/(c2b*c2t) - 0.125D0*(mb*Yb)/(c2t*s2b) + 
     -   0.125D0*(mt*Yb)/(c2b*s2t) - 0.0625D0*Yb**2/(c2b*c2t)
        tmp38 = -(0.0625D0*b/(c2b*c2t)) - 
     -   0.125D0*(mb*mt)/(s2b*s2t) - 0.0625D0*t/(c2b*c2t) + 
     -   0.125D0*(mb*Yb)/(c2t*s2b) - 0.125D0*(mt*Yb)/(c2b*s2t) + 
     -   0.0625D0*Yb**2/(c2b*c2t)
	tmp39 = 0.25D0*Yb/c2b - 0.25D0*Yt/c2b
	tmp40 = -(0.25D0*Yb/c2b) + 0.25D0*Yt/c2b
	tmp41 = 0.125D0*Yb/c2t**3 - 0.125D0*Yt/c2t**3
	tmp42 = -(0.125D0*Yb/c2t**3) + 0.125D0*Yt/c2t**3
	tmp43 = 0.25D0*Yb/c2t - 0.25D0*Yt/c2t
	tmp44 = -(0.25D0*Yb/c2t) + 0.25D0*Yt/c2t
        tmp45 = 0.25D0*((1d0 + c2b)*(1d0 + c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1d0 + c2b)*s2t*Yt)/mt
	tmp46 = (-2*mt*Yt)/s2t - Yt**2
	tmp47 = (2*mt*Yt)/s2t - Yt**2
        tmp48 = 0.0625D0*b/(c2b*c2t) + 0.125D0*(mb*mt)/(s2b*s2t) + 
     -   0.0625D0*t/(c2b*c2t) + 0.125D0*(mb*Yt)/(c2t*s2b) - 
     -   0.125D0*(mt*Yt)/(c2b*s2t) - 0.0625D0*Yt**2/(c2b*c2t)
        tmp49 = -(0.0625D0*b/(c2b*c2t)) - 
     -   0.125D0*(mb*mt)/(s2b*s2t) - 0.0625D0*t/(c2b*c2t) - 
     -   0.125D0*(mb*Yt)/(c2t*s2b) + 0.125D0*(mt*Yt)/(c2b*s2t) + 
     -   0.0625D0*Yt**2/(c2b*c2t)
	tmp50 = 4*t - 4*mt*s2t*Yt + s2t**2*Yt**2
	tmp51 = 4*t + 4*mt*s2t*Yt + s2t**2*Yt**2
	tmp52 = 0.0625D0*(1d0 -c2b)/c2t**3 -0.0625D0*(1d0 +c2b)/c2t**3
	tmp53 = -(0.0625D0*(1d0 - c2b)/c2t**3)
     $       + 0.0625D0*(1d0 + c2b)/c2t**3
	tmp54 = 0.125D0*(1d0 - c2b)/c2t - 0.125D0*(1d0 + c2b)/c2t
	tmp55 = -(0.125D0*(1d0 - c2b)/c2t) + 0.125D0*(1d0 + c2b)/c2t
	tmp56 = 0.25D0*((1d0 + c2b)*(1d0 - c2t))
     $       + 0.25D0*((1d0 - c2b)*(1d0 + c2t))
	tmp57 = 0.125D0*(1d0 - c2t)/c2b - 0.125D0*(1d0 + c2t)/c2b
	tmp58 = -(0.125D0*(1d0 - c2t)/c2b) + 0.125D0*(1d0 + c2t)/c2b
	tmp59 = 0.25D0*((1d0 - c2b)*(1d0 - c2t))
     $       + 0.25D0*((1d0 + c2b)*(1d0 + c2t))
	tmp60 = 0.5D0*((1d0 + c2b)*hb**2) +0.5D0*((1d0 -c2b)*ht**2)
	tmp61 = 0.5D0*((1d0 - c2b)*hb**2) +0.5D0*((1d0 +c2b)*ht**2)
	tmp62 = 0.5D0*((1d0 + c2t)*hb**2) + 0.5D0*((1d0 - c2t)*ht**2)
	tmp63 = 0.5D0*((1d0 - c2t)*hb**2) + 0.5D0*((1d0 + c2t)*ht**2)
	tmp64 = (A0 - b)/b + LogA0 - Logb - t/b
	tmp65 = (-LogA0 + Logt)*(A0 - t) + (-LogA0 + Logt)*t
	tmp66 = (A0 - b)*(-LogA0 + Logb) + (-LogA0 - Logb + 2*Logt)*t
        tmp67 = (-LogB1 + Logmu2)*(B1 - mu2) + 
     -   (-LogB1 - Logmu2 + 2*Logt)*t
        tmp68 = (-LogB2 + Logmu2)*(B2 - mu2) + 
     -   (-LogB2 - Logmu2 + 2*Logt)*t
	tmp69 = b*(-LogA0 + 2*Logb - Logt) + (LogA0 - Logt)*(-A0 + t)
	tmp70 = (-LogA0 + Logt)*t + (LogA0 - Logt)*(-A0 + t)
        tmp71 = B1*(2*LogB1 - Logmu2 - Logt) + 
     -   (Logmu2 - Logt)*(-mu2 + t)
        tmp72 = B2*(2*LogB2 - Logmu2 - Logt) + 
     -   (Logmu2 - Logt)*(-mu2 + t)
	tmp73 = Logmu2 - Logt - B1/t - (-mu2 + t)/t
	tmp74 = Logmu2 - Logt - B2/t - (-mu2 + t)/t
        tmp75 = 2*B1*LogB1 - 2*B2*LogB2 - 
     -   0.5D0*(deltB1tmu2*phiB1tmu2)/mu2 + 
     -   0.5D0*(deltB2tmu2*phiB2tmu2)/mu2 + 
     -   0.5D0*(Logmu2*Logt*(B1 - mu2 - t)) - 
     -   0.5D0*(Logmu2*Logt*(B2 - mu2 - t)) + 
     -   0.5D0*(LogB1*Logt*(-B1 + mu2 - t)) - 
     -   0.5D0*(LogB2*Logt*(-B2 + mu2 - t)) + 
     -   0.5D0*(LogB1*Logmu2*(-B1 - mu2 + t)) - 
     -   0.5D0*(LogB2*Logmu2*(-B2 - mu2 + t)) - 2.5D0*(B1 + mu2 + t) + 
     -   2.5D0*(B2 + mu2 + t)
        tmp76 = -5*b + 4*b*Logb + Logt**2*(b - t) - 5*t - 
     -   2*Li2bt*(-b + t) + Logt*(-2*b*Logb + 4*t)
	tmp77 = -6d0 +4*Logt +(2d0 -2*Logt)*Logt + (4*t - 2*Logt*t)/t
	tmp78 = -Logb + Logmu2 - (b - mu2)/b - T1/b
	tmp79 = 2*(-1d0 + LogT1)*T1 + 2*(-1d0 + LogT1)**2*T1
	tmp80 = (-LogA0 + LogT1)*(A0 - T1) + (-LogA0 + LogT1)*T1
        tmp81 = (A0 - B1)*(-LogA0 + LogB1) + 
     -   (-LogA0 - LogB1 + 2*LogT1)*T1
        tmp82 = (A0 - B2)*(-LogA0 + LogB2) + 
     -   (-LogA0 - LogB2 + 2*LogT1)*T1
        tmp83 = (-Logb + Logmu2)*(b - mu2) + 
     -   (-Logb - Logmu2 + 2*LogT1)*T1
        tmp84 = (Logmu2 - Logt)*(-mu2 + t) + 
     -   (-Logmu2 - Logt + 2*LogT1)*T1
        tmp85 = B1*(-LogA0 + 2*LogB1 - LogT1) + 
     -   (LogA0 - LogT1)*(-A0 + T1)
        tmp86 = B2*(-LogA0 + 2*LogB2 - LogT1) + 
     -   (LogA0 - LogT1)*(-A0 + T1)
	tmp87 = (-LogA0 + LogT1)*T1 + (LogA0 - LogT1)*(-A0 + T1)
        tmp88 = 2*A0*LogA0 + 2*B1*LogB1 + 2*LogT1*T1 + 
     -   0.5D0*(LogB1*LogT1*(A0 - B1 - T1)) + 
     -   0.5D0*(LogA0*LogT1*(-A0 + B1 - T1)) - 
     -   0.5D0*(deltA0B1T1*phiA0B1T1)/T1 + 
     -   0.5D0*(LogA0*LogB1*(-A0 - B1 + T1)) - 2.5D0*(A0 + B1 + T1)
        tmp89 = 2*A0*LogA0 + 2*B2*LogB2 + 2*LogT1*T1 + 
     -   0.5D0*(LogB2*LogT1*(A0 - B2 - T1)) + 
     -   0.5D0*(LogA0*LogT1*(-A0 + B2 - T1)) - 
     -   0.5D0*(deltA0B2T1*phiA0B2T1)/T1 + 
     -   0.5D0*(LogA0*LogB2*(-A0 - B2 + T1)) - 2.5D0*(A0 + B2 + T1)
        tmp90 = b*(2*Logb - Logmu2 - LogT1) + 
     -   (Logmu2 - LogT1)*(-mu2 + T1)
        tmp91 = (-Logmu2 + 2*Logt - LogT1)*t + 
     -   (Logmu2 - LogT1)*(-mu2 + T1)
        tmp92 = 2*b*Logb + 2*Logmu2*mu2 + 2*LogT1*T1 - 
     -   0.5D0*(deltT1bmu2*phiT1bmu2)/mu2 + 
     -   0.5D0*(Logmu2*LogT1*(b - mu2 - T1)) + 
     -   0.5D0*(Logb*LogT1*(-b + mu2 - T1)) + 
     -   0.5D0*(Logb*Logmu2*(-b - mu2 + T1)) - 2.5D0*(b + mu2 + T1)
        tmp93 = 2*A0*LogA0 - A0*LogA0*LogT1 + 4*LogT1*T1 + 
     -   0.5D0*(LogT1**2*(A0 - 2*T1)) - 
     -   0.5D0*(deltA0T1T1*phiA0T1T1)/T1 - 2.5D0*(A0 + 2*T1)
	tmp94 = -10*T1 + 4*LogT1*T1 + LogT1*(4*T1 - 2*LogT1*T1)
        tmp95 = -6d0 + 4*LogT1 + (2d0 - 2*LogT1)*LogT1 + 
     -   (4*T1 - 2*LogT1*T1)/T1
        tmp96 = (-LogA0 + 2*LogT1 - LogT2)*T1 + 
     -   (-LogA0 + LogT2)*(A0 - T2)
	tmp97 = -Logb + Logmu2 - (b - mu2)/b - T2/b
	tmp98 = 2*(-1d0 + LogT2)*T2 + 2*(-1d0 + LogT2)**2*T2
	tmp99 = (-LogA0 + LogT2)*(A0 - T2) + (-LogA0 + LogT2)*T2
        tmp100 = (A0 - B1)*(-LogA0 + LogB1) + 
     -   (-LogA0 - LogB1 + 2*LogT2)*T2
        tmp101 = (A0 - B2)*(-LogA0 + LogB2) + 
     -   (-LogA0 - LogB2 + 2*LogT2)*T2
        tmp102 = (-Logb + Logmu2)*(b - mu2) + 
     -   (-Logb - Logmu2 + 2*LogT2)*T2
        tmp103 = (Logmu2 - Logt)*(-mu2 + t) + 
     -   (-Logmu2 - Logt + 2*LogT2)*T2
        tmp104 = (LogA0 - LogT1)*(-A0 + T1) + 
     -   (-LogA0 - LogT1 + 2*LogT2)*T2
        tmp105 = B1*(-LogA0 + 2*LogB1 - LogT2) + 
     -   (LogA0 - LogT2)*(-A0 + T2)
        tmp106 = B2*(-LogA0 + 2*LogB2 - LogT2) + 
     -   (LogA0 - LogT2)*(-A0 + T2)
	tmp107 = (-LogA0 + LogT2)*T2 + (LogA0 - LogT2)*(-A0 + T2)
        tmp108 = 2*A0*LogA0 + 2*B1*LogB1 + 2*LogT2*T2 + 
     -   0.5D0*(LogB1*LogT2*(A0 - B1 - T2)) + 
     -   0.5D0*(LogA0*LogT2*(-A0 + B1 - T2)) - 
     -   0.5D0*(deltA0B1T2*phiA0B1T2)/T2 + 
     -   0.5D0*(LogA0*LogB1*(-A0 - B1 + T2)) - 2.5D0*(A0 + B1 + T2)
        tmp109 = 2*A0*LogA0 + 2*B2*LogB2 + 2*LogT2*T2 + 
     -   0.5D0*(LogB2*LogT2*(A0 - B2 - T2)) + 
     -   0.5D0*(LogA0*LogT2*(-A0 + B2 - T2)) - 
     -   0.5D0*(deltA0B2T2*phiA0B2T2)/T2 + 
     -   0.5D0*(LogA0*LogB2*(-A0 - B2 + T2)) - 2.5D0*(A0 + B2 + T2)
        tmp110 = b*(2*Logb - Logmu2 - LogT2) + 
     -   (Logmu2 - LogT2)*(-mu2 + T2)
        tmp111 = (-Logmu2 + 2*Logt - LogT2)*t + 
     -   (Logmu2 - LogT2)*(-mu2 + T2)
        tmp112 = 2*b*Logb + 2*Logmu2*mu2 + 2*LogT2*T2 - 
     -   0.5D0*(deltT2bmu2*phiT2bmu2)/mu2 + 
     -   0.5D0*(Logmu2*LogT2*(b - mu2 - T2)) + 
     -   0.5D0*(Logb*LogT2*(-b + mu2 - T2)) + 
     -   0.5D0*(Logb*Logmu2*(-b - mu2 + T2)) - 2.5D0*(b + mu2 + T2)
        tmp113 = 2*LogT1*T1 - 2*LogT2*T2 - 
     -   0.5D0*(deltT1bmu2*phiT1bmu2)/mu2 + 
     -   0.5D0*(deltT2bmu2*phiT2bmu2)/mu2 + 
     -   0.5D0*(Logmu2*LogT1*(b - mu2 - T1)) + 
     -   0.5D0*(Logb*LogT1*(-b + mu2 - T1)) + 
     -   0.5D0*(Logb*Logmu2*(-b - mu2 + T1)) - 2.5D0*(b + mu2 + T1) - 
     -   0.5D0*(Logmu2*LogT2*(b - mu2 - T2)) - 
     -   0.5D0*(Logb*LogT2*(-b + mu2 - T2)) - 
     -   0.5D0*(Logb*Logmu2*(-b - mu2 + T2)) + 2.5D0*(b + mu2 + T2)
        tmp114 = 2*A0*LogA0 - A0*LogA0*LogT2 + 4*LogT2*T2 + 
     -   0.5D0*(LogT2**2*(A0 - 2*T2)) - 
     -   0.5D0*(deltA0T2T2*phiA0T2T2)/T2 - 2.5D0*(A0 + 2*T2)
	tmp115 = -10*T2 + 4*LogT2*T2 + LogT2*(4*T2 - 2*LogT2*T2)
        tmp116 = -6d0 + 4*LogT2 + (2d0 - 2*LogT2)*LogT2 + 
     -   (4*T2 - 2*LogT2*T2)/T2
        tmp117 = 0.25D0*((1d0 - c2b)*(1d0 + c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1d0 - c2b)*s2t*Xb)/mt
        tmp118 = 0.25D0*((1d0 - c2b)*(1d0 - c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1d0 - c2b)*s2t*Xb)/mt
        tmp119 = 0.25D0*((1d0 + c2b)*(1d0 - c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1d0 + c2b)*s2t*Xb)/mt
        tmp120 = 0.25D0*(b*(1d0 + c2b)*(1d0 - c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 - c2b)*(1d0 + c2t)*t) + 
     -   0.5D0*((1d0 -c2t)*mb*s2b*Xb) -0.5D0*((1d0 - c2b)*mt*s2t*Xb) + 
     -   0.25D0*((1d0 - c2b)*(1d0 - c2t)*Xb**2)
        tmp121 = 0.25D0*(b*(1d0 - c2b)*(1d0 - c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 + c2b)*(1d0 + c2t)*t) - 
     -   0.5D0*((1d0 - c2t)*mb*s2b*Xb) - 0.5D0*((1d0 + c2b)*mt*s2t*Xb) + 
     -   0.25D0*((1d0 + c2b)*(1d0 - c2t)*Xb**2)
        tmp122 = -(0.125D0*(b*(1d0 + c2b))/c2t) + 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1d0 - c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Xb)/c2t + 0.25D0*((1d0 - c2b)*mt*Xb)/s2t - 
     -   0.125D0*((1d0 - c2b)*Xb**2)/c2t
        tmp123 = 0.125D0*(b*(1d0 + c2b))/c2t - 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1d0 - c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Xb)/c2t - 0.25D0*((1d0 - c2b)*mt*Xb)/s2t + 
     -   0.125D0*((1d0 - c2b)*Xb**2)/c2t
        tmp124 = -(0.125D0*(b*(1d0 - c2b))/c2t) - 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1d0 + c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Xb)/c2t + 0.25D0*((1d0 + c2b)*mt*Xb)/s2t - 
     -   0.125D0*((1d0 + c2b)*Xb**2)/c2t
        tmp125 = 0.125D0*(b*(1d0 - c2b))/c2t + 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1d0 + c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Xb)/c2t - 0.25D0*((1d0 + c2b)*mt*Xb)/s2t + 
     -   0.125D0*((1d0 + c2b)*Xb**2)/c2t
        tmp126 = 0.25D0*(b*(1d0 + c2b)*(1d0 + c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 - c2b)*(1d0 - c2t)*t) + 
     -   0.5D0*((1d0 + c2t)*mb*s2b*Xb) +0.5D0*((1d0 -c2b)*mt*s2t*Xb) + 
     -   0.25D0*((1d0 - c2b)*(1d0 + c2t)*Xb**2)
        tmp127 = 0.25D0*(b*(1d0 - c2b)*(1d0 + c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 + c2b)*(1d0 - c2t)*t) - 
     -   0.5D0*((1d0 + c2t)*mb*s2b*Xb) +0.5D0*((1d0 +c2b)*mt*s2t*Xb) + 
     -   0.25D0*((1d0 + c2b)*(1d0 + c2t)*Xb**2)
	tmp128 = 0.5D0*((1d0 + c2b)*Xb) + 0.5D0*((1d0 - c2b)*Xt)
	tmp129 = 0.5D0*((1d0 - c2b)*Xb) + 0.5D0*((1d0 + c2b)*Xt)
	tmp130 = 0.5D0*((1d0 + c2t)*Xb) + 0.5D0*((1d0 - c2t)*Xt)
	tmp131 = 0.5D0*((1d0 - c2t)*Xb) + 0.5D0*((1d0 + c2t)*Xt)
        tmp132 = 0.25D0*((1d0 - c2b)*(1d0 - c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1d0 - c2b)*s2t*Xt)/mt
        tmp133 = 0.25D0*((1d0 - c2b)*(1d0 + c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1d0 - c2b)*s2t*Xt)/mt
        tmp134 = 0.25D0*((1d0 + c2b)*(1d0 - c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1d0 + c2b)*s2t*Xt)/mt
        tmp135 = 0.25D0*(b*(1d0 + c2b)*(1d0 - c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 -c2b)*(1d0 + c2t)*t) - 
     -   0.5D0*((1d0 -c2t)*mb*s2b*Xt) + 0.5D0*((1d0 -c2b)*mt*s2t*Xt) + 
     -   0.25D0*((1d0 -c2b)*(1d0 - c2t)*Xt**2)
        tmp136 = 0.25D0*(b*(1d0 - c2b)*(1d0 - c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 + c2b)*(1d0 + c2t)*t) + 
     -   0.5D0*((1d0 - c2t)*mb*s2b*Xt) + 0.5D0*((1d0 + c2b)*mt*s2t*Xt) + 
     -   0.25D0*((1d0 + c2b)*(1d0 - c2t)*Xt**2)
        tmp137 = -(0.125D0*(b*(1d0 + c2b))/c2t) + 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1d0 - c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Xt)/c2t - 0.25D0*((1d0 - c2b)*mt*Xt)/s2t - 
     -   0.125D0*((1d0 - c2b)*Xt**2)/c2t
        tmp138 = 0.125D0*(b*(1d0 + c2b))/c2t - 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1d0 - c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Xt)/c2t + 0.25D0*((1d0 - c2b)*mt*Xt)/s2t + 
     -   0.125D0*((1d0 - c2b)*Xt**2)/c2t
        tmp139 = -(0.125D0*(b*(1d0 - c2b))/c2t) - 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1d0 + c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Xt)/c2t - 0.25D0*((1d0 + c2b)*mt*Xt)/s2t - 
     -   0.125D0*((1d0 + c2b)*Xt**2)/c2t
        tmp140 = 0.125D0*(b*(1d0 - c2b))/c2t + 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1d0 + c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Xt)/c2t + 0.25D0*((1d0 + c2b)*mt*Xt)/s2t + 
     -   0.125D0*((1d0 + c2b)*Xt**2)/c2t
        tmp141 = 0.25D0*(b*(1d0 + c2b)*(1d0 + c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1d0 - c2b)*(1d0 - c2t)*t) - 
     -   0.5D0*((1d0 +c2t)*mb*s2b*Xt) - 0.5D0*((1d0 -c2b)*mt*s2t*Xt) + 
     -   0.25D0*((1d0 - c2b)*(1d0 + c2t)*Xt**2)
        tmp142 = 0.25D0*(b*(1d0 - c2b)*(1d0 + c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 + c2b)*(1 - c2t)*t) + 
     -   0.5D0*((1 + c2t)*mb*s2b*Xt) - 0.5D0*((1 + c2b)*mt*s2t*Xt) + 
     -   0.25D0*((1 + c2b)*(1 + c2t)*Xt**2)
        tmp143 = 0.25D0*((1 - c2b)*(1 + c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1 - c2b)*s2t*Yb)/mt
        tmp144 = 0.25D0*((1 - c2b)*(1 - c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1 - c2b)*s2t*Yb)/mt
        tmp145 = 0.25D0*((1 + c2b)*(1 - c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1 + c2b)*s2t*Yb)/mt
        tmp146 = 0.25D0*(b*(1 + c2b)*(1 - c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 - c2b)*(1 + c2t)*t) + 
     -   0.5D0*((1 - c2t)*mb*s2b*Yb) - 0.5D0*((1 - c2b)*mt*s2t*Yb) + 
     -   0.25D0*((1 - c2b)*(1 - c2t)*Yb**2)
        tmp147 = 0.25D0*(b*(1 - c2b)*(1 - c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 + c2b)*(1 + c2t)*t) - 
     -   0.5D0*((1 - c2t)*mb*s2b*Yb) - 0.5D0*((1 + c2b)*mt*s2t*Yb) + 
     -   0.25D0*((1 + c2b)*(1 - c2t)*Yb**2)
        tmp148 = -(0.125D0*(b*(1 + c2b))/c2t) + 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1 - c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Yb)/c2t + 0.25D0*((1 - c2b)*mt*Yb)/s2t - 
     -   0.125D0*((1 - c2b)*Yb**2)/c2t
        tmp149 = 0.125D0*(b*(1 + c2b))/c2t - 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1 - c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Yb)/c2t - 0.25D0*((1 - c2b)*mt*Yb)/s2t + 
     -   0.125D0*((1 - c2b)*Yb**2)/c2t
        tmp150 = -(0.125D0*(b*(1 - c2b))/c2t) - 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1 + c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Yb)/c2t + 0.25D0*((1 + c2b)*mt*Yb)/s2t - 
     -   0.125D0*((1 + c2b)*Yb**2)/c2t
        tmp151 = 0.125D0*(b*(1 - c2b))/c2t + 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1 + c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Yb)/c2t - 0.25D0*((1 + c2b)*mt*Yb)/s2t + 
     -   0.125D0*((1 + c2b)*Yb**2)/c2t
        tmp152 = 0.25D0*(b*(1 + c2b)*(1 + c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 - c2b)*(1 - c2t)*t) + 
     -   0.5D0*((1 + c2t)*mb*s2b*Yb) + 0.5D0*((1 - c2b)*mt*s2t*Yb) + 
     -   0.25D0*((1 - c2b)*(1 + c2t)*Yb**2)
        tmp153 = 0.25D0*(b*(1 - c2b)*(1 + c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 + c2b)*(1 - c2t)*t) - 
     -   0.5D0*((1 + c2t)*mb*s2b*Yb) + 0.5D0*((1 + c2b)*mt*s2t*Yb) + 
     -   0.25D0*((1 + c2b)*(1 + c2t)*Yb**2)
	tmp154 = 0.5D0*((1 + c2b)*Yb) + 0.5D0*((1 - c2b)*Yt)
	tmp155 = 0.5D0*((1 - c2b)*Yb) + 0.5D0*((1 + c2b)*Yt)
	tmp156 = 0.5D0*((1 + c2t)*Yb) + 0.5D0*((1 - c2t)*Yt)
	tmp157 = 0.5D0*((1 - c2t)*Yb) + 0.5D0*((1 + c2t)*Yt)
        tmp158 = 0.25D0*((1 - c2b)*(1 - c2t)) + 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1 - c2b)*s2t*Yt)/mt
        tmp159 = 0.25D0*((1 - c2b)*(1 + c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt + 0.25D0*((1 - c2b)*s2t*Yt)/mt
        tmp160 = 0.25D0*((1 + c2b)*(1 - c2t)) - 
     -   0.25D0*(mb*s2b*s2t)/mt - 0.25D0*((1 + c2b)*s2t*Yt)/mt
        tmp161 = 0.25D0*(b*(1 + c2b)*(1 - c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 - c2b)*(1 + c2t)*t) - 
     -   0.5D0*((1 - c2t)*mb*s2b*Yt) + 0.5D0*((1 - c2b)*mt*s2t*Yt) + 
     -   0.25D0*((1 - c2b)*(1 - c2t)*Yt**2)
        tmp162 = 0.25D0*(b*(1 - c2b)*(1 - c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 + c2b)*(1 + c2t)*t) + 
     -   0.5D0*((1 - c2t)*mb*s2b*Yt) + 0.5D0*((1 + c2b)*mt*s2t*Yt) + 
     -   0.25D0*((1 + c2b)*(1 - c2t)*Yt**2)
        tmp163 = -(0.125D0*(b*(1 + c2b))/c2t) + 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1 - c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Yt)/c2t - 0.25D0*((1 - c2b)*mt*Yt)/s2t - 
     -   0.125D0*((1 - c2b)*Yt**2)/c2t
        tmp164 = 0.125D0*(b*(1 + c2b))/c2t - 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1 - c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Yt)/c2t + 0.25D0*((1 - c2b)*mt*Yt)/s2t + 
     -   0.125D0*((1 - c2b)*Yt**2)/c2t
        tmp165 = -(0.125D0*(b*(1 - c2b))/c2t) - 
     -   0.25D0*(mb*mt*s2b)/s2t + 0.125D0*((1 + c2b)*t)/c2t - 
     -   0.25D0*(mb*s2b*Yt)/c2t - 0.25D0*((1 + c2b)*mt*Yt)/s2t - 
     -   0.125D0*((1 + c2b)*Yt**2)/c2t
        tmp166 = 0.125D0*(b*(1 - c2b))/c2t + 
     -   0.25D0*(mb*mt*s2b)/s2t - 0.125D0*((1 + c2b)*t)/c2t + 
     -   0.25D0*(mb*s2b*Yt)/c2t + 0.25D0*((1 + c2b)*mt*Yt)/s2t + 
     -   0.125D0*((1 + c2b)*Yt**2)/c2t
        tmp167 = 0.25D0*(b*(1 + c2b)*(1 + c2t)) + 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 - c2b)*(1 - c2t)*t) - 
     -   0.5D0*((1 + c2t)*mb*s2b*Yt) - 0.5D0*((1 - c2b)*mt*s2t*Yt) + 
     -   0.25D0*((1 - c2b)*(1 + c2t)*Yt**2)
        tmp168 = 0.25D0*(b*(1 - c2b)*(1 + c2t)) - 
     -   0.5D0*(mb*mt*s2b*s2t) + 0.25D0*((1 + c2b)*(1 - c2t)*t) + 
     -   0.5D0*((1 + c2t)*mb*s2b*Yt) - 0.5D0*((1 + c2b)*mt*s2t*Yt) + 
     -   0.25D0*((1 + c2b)*(1 + c2t)*Yt**2)
	tmp169 = -Li2T1B1 - 0.5D0*(-LogB1 + LogT1)**2
	tmp170 = -Li2T1B2 - 0.5D0*(-LogB2 + LogT1)**2
	tmp171 = -Li2T1T2 - 0.5D0*(LogT1 - LogT2)**2
	tmp172 = -Li2T2B1 - 0.5D0*(-LogB1 + LogT2)**2
	tmp173 = -Li2T2B2 - 0.5D0*(-LogB2 + LogT2)**2
	tmp174 = -A0 + (A0 - t)**2/b - t
	tmp175 = -A0 + (A0 - t)**2/t - t
	tmp176 = -A0 + (A0 - T1)**2/B1 - T1
	tmp177 = -A0 + (A0 - T1)**2/B2 - T1
	tmp178 = -A0 + (A0 - T1)**2/T1 - T1
	tmp179 = -A0 - T1 + (A0 - T1)**2/T2
	tmp180 = -A0 + (A0 - T2)**2/B1 - T2
	tmp181 = -A0 + (A0 - T2)**2/B2 - T2
	tmp182 = -A0 + (A0 - T2)**2/T2 - T2
        tmp183 = cbe**2*hb**2*tmp119 + ht**2*sbe**2*tmp133 - 
     -   cbe*hb*ht*sbe*((mb*tmp56)/mt - 0.5D0*(s2b*s2t) - 
     -      0.5D0*(s2b*tmp130)/mt)
        tmp184 = cbe**2*hb**2*tmp118 + ht**2*sbe**2*tmp29 - 
     -   cbe*hb*ht*sbe*((mb*tmp59)/mt + 0.5D0*(s2b*s2t) + 
     -      0.5D0*(s2b*tmp130)/mt)
        tmp185 = cbe**2*hb**2*tmp127 + ht**2*sbe**2*tmp135 - 
     -   cbe*hb*ht*sbe*(mb*s2t*tmp128 - mt*s2b*tmp130 + 
     -      2*mb*mt*tmp56 - 0.5D0*(s2b*s2t*(b + t)) - 
     -      0.5D0*(s2b*s2t*Xb*Xt))
        tmp186 = cbe**2*hb**2*tmp126 + ht**2*sbe**2*tmp136 - 
     -   cbe*hb*ht*sbe*(mb*s2t*tmp129 + mt*s2b*tmp130 + 
     -      2*mb*mt*tmp59 + 0.5D0*(s2b*s2t*(b + t)) + 
     -      0.5D0*(s2b*s2t*Xb*Xt))
        tmp187 = ht**2*sbe**2*tmp132 + cbe**2*hb**2*tmp20 - 
     -   cbe*hb*ht*sbe*((mb*tmp59)/mt + 0.5D0*(s2b*s2t) - 
     -      0.5D0*(s2b*tmp131)/mt)
        tmp188 = cbe**2*hb**2*tmp117 + ht**2*sbe**2*tmp134 - 
     -   cbe*hb*ht*sbe*((mb*tmp56)/mt - 0.5D0*(s2b*s2t) + 
     -      0.5D0*(s2b*tmp131)/mt)
        tmp189 = cbe**2*hb**2*tmp121 + ht**2*sbe**2*tmp141 - 
     -   cbe*hb*ht*sbe*(-(mb*s2t*tmp128) - mt*s2b*tmp131 + 
     -      2*mb*mt*tmp59 + 0.5D0*(s2b*s2t*(b + t)) + 
     -      0.5D0*(s2b*s2t*Xb*Xt))
        tmp190 = cbe**2*hb**2*tmp120 + ht**2*sbe**2*tmp142 - 
     -   cbe*hb*ht*sbe*(-(mb*s2t*tmp129) + mt*s2b*tmp131 + 
     -      2*mb*mt*tmp56 - 0.5D0*(s2b*s2t*(b + t)) - 
     -      0.5D0*(s2b*s2t*Xb*Xt))
        tmp191 = hb**2*sbe**2*tmp145 + cbe**2*ht**2*tmp159 + 
     -   cbe*hb*ht*sbe*((mb*tmp56)/mt - 0.5D0*(s2b*s2t) - 
     -      0.5D0*(s2b*tmp156)/mt)
        tmp192 = hb**2*sbe**2*tmp144 + cbe**2*ht**2*tmp45 + 
     -   cbe*hb*ht*sbe*((mb*tmp59)/mt + 0.5D0*(s2b*s2t) + 
     -      0.5D0*(s2b*tmp156)/mt)
        tmp193 = hb**2*sbe**2*tmp153 + cbe**2*ht**2*tmp161 + 
     -   cbe*hb*ht*sbe*(mb*s2t*tmp154 - mt*s2b*tmp156 + 
     -      2*mb*mt*tmp56 - 0.5D0*(s2b*s2t*(b + t)) - 
     -      0.5D0*(s2b*s2t*Yb*Yt))
        tmp194 = hb**2*sbe**2*tmp152 + cbe**2*ht**2*tmp162 + 
     -   cbe*hb*ht*sbe*(mb*s2t*tmp155 + mt*s2b*tmp156 + 
     -      2*mb*mt*tmp59 + 0.5D0*(s2b*s2t*(b + t)) + 
     -      0.5D0*(s2b*s2t*Yb*Yt))
        tmp195 = cbe**2*ht**2*tmp158 + hb**2*sbe**2*tmp36 + 
     -   cbe*hb*ht*sbe*((mb*tmp59)/mt + 0.5D0*(s2b*s2t) - 
     -      0.5D0*(s2b*tmp157)/mt)
        tmp196 = hb**2*sbe**2*tmp143 + cbe**2*ht**2*tmp160 + 
     -   cbe*hb*ht*sbe*((mb*tmp56)/mt - 0.5D0*(s2b*s2t) + 
     -      0.5D0*(s2b*tmp157)/mt)
        tmp197 = hb**2*sbe**2*tmp147 + cbe**2*ht**2*tmp167 + 
     -   cbe*hb*ht*sbe*(-(mb*s2t*tmp154) - mt*s2b*tmp157 + 
     -      2*mb*mt*tmp59 + 0.5D0*(s2b*s2t*(b + t)) + 
     -      0.5D0*(s2b*s2t*Yb*Yt))
        tmp198 = hb**2*sbe**2*tmp146 + cbe**2*ht**2*tmp168 + 
     -   cbe*hb*ht*sbe*(-(mb*s2t*tmp155) + mt*s2b*tmp157 + 
     -      2*mb*mt*tmp56 - 0.5D0*(s2b*s2t*(b + t)) - 
     -      0.5D0*(s2b*s2t*Yb*Yt))
        tmp199 = cbe**2*hb**2*tmp125 + ht**2*sbe**2*tmp137 - 
     -   cbe*hb*ht*sbe*(-(mt*s2b*tmp27) + 2*mb*mt*tmp54 + 
     -      0.25D0*(s2b*(b + t))/s2t - 0.5D0*(mb*tmp128)/s2t + 
     -      0.25D0*(s2b*Xb*Xt)/s2t)
        tmp200 = cbe**2*hb**2*tmp123 + ht**2*sbe**2*tmp139 - 
     -   cbe*hb*ht*sbe*(mt*s2b*tmp27 + 2*mb*mt*tmp55 - 
     -      0.25D0*(s2b*(b + t))/s2t - 0.5D0*(mb*tmp129)/s2t - 
     -      0.25D0*(s2b*Xb*Xt)/s2t)
        tmp201 = cbe**2*hb**2*tmp124 + ht**2*sbe**2*tmp138 - 
     -   cbe*hb*ht*sbe*(-(mt*s2b*tmp28) + 2*mb*mt*tmp55 - 
     -      0.25D0*(s2b*(b + t))/s2t + 0.5D0*(mb*tmp128)/s2t - 
     -      0.25D0*(s2b*Xb*Xt)/s2t)
        tmp202 = cbe**2*hb**2*tmp122 + ht**2*sbe**2*tmp140 - 
     -   cbe*hb*ht*sbe*(mt*s2b*tmp28 + 2*mb*mt*tmp54 + 
     -      0.25D0*(s2b*(b + t))/s2t + 0.5D0*(mb*tmp129)/s2t + 
     -      0.25D0*(s2b*Xb*Xt)/s2t)
        tmp203 = hb**2*sbe**2*tmp151 + cbe**2*ht**2*tmp163 + 
     -   cbe*hb*ht*sbe*(-(mt*s2b*tmp43) + 2*mb*mt*tmp54 + 
     -      0.25D0*(s2b*(b + t))/s2t - 0.5D0*(mb*tmp154)/s2t + 
     -      0.25D0*(s2b*Yb*Yt)/s2t)
        tmp204 = hb**2*sbe**2*tmp149 + cbe**2*ht**2*tmp165 + 
     -   cbe*hb*ht*sbe*(mt*s2b*tmp43 + 2*mb*mt*tmp55 - 
     -      0.25D0*(s2b*(b + t))/s2t - 0.5D0*(mb*tmp155)/s2t - 
     -      0.25D0*(s2b*Yb*Yt)/s2t)
        tmp205 = hb**2*sbe**2*tmp150 + cbe**2*ht**2*tmp164 + 
     -   cbe*hb*ht*sbe*(-(mt*s2b*tmp44) + 2*mb*mt*tmp55 - 
     -      0.25D0*(s2b*(b + t))/s2t + 0.5D0*(mb*tmp154)/s2t - 
     -      0.25D0*(s2b*Yb*Yt)/s2t)
        tmp206 = hb**2*sbe**2*tmp148 + cbe**2*ht**2*tmp166 + 
     -   cbe*hb*ht*sbe*(mt*s2b*tmp44 + 2*mb*mt*tmp54 + 
     -      0.25D0*(s2b*(b + t))/s2t + 0.5D0*(mb*tmp155)/s2t + 
     -      0.25D0*(s2b*Yb*Yt)/s2t)
        tmp207 = (b**2*(Logb - Logt))/((1 - b/t)**2*t**4) + 
     -   b/((1 - b/t)*t**3) + (2*b*(Logb - Logt))/((1 - b/t)*t**3)
        tmp208 = (-2*b*(Logb - Logt))/((1 - b/t)*t**2) + 
     -   (-2 - 2*Logb)/t + (2*Logt)/t - 
     -   (2*(Logb - Logt))/((1 - b/t)*t) + 
     -   (2*b*(Logb - Logt)*(-b + t))/((1 - b/t)**2*t**3) + 
     -   (2*(-b + t))/((1 - b/t)*t**2) + 
     -   (2*(Logb - Logt)*(-b + t))/((1 - b/t)*t**2)
        tmp209 = -1 + 2*Li2bt + 4*Logb + (-2 - 2*Logb)*Logt + 
     -   Logt**2 - (2*(Logb - Logt)*(-b + t))/((1 - b/t)*t)
        tmp210 = -5 - 2*Li2bt + 4*Logt - Logt**2 + 
     -   (2*Logt*(b - t))/t + 
     -   (2*b*(Logb - Logt)*(-b + t))/((1 - b/t)*t**2) + 
     -   (-2*b*Logb + 4*t)/t
        tmp211 = -1 + 4*LogB1 + (-2 - 2*LogB1)*LogT1 + 
     -   LogT1**2 - (2*(LogB1 - LogT1)*(-B1 + T1))/
     -    ((1 - B1/T1)*T1) + 2*tmp169
        tmp212 = -1 + 4*LogB2 + (-2 - 2*LogB2)*LogT1 + 
     -   LogT1**2 - (2*(LogB2 - LogT1)*(-B2 + T1))/
     -    ((1 - B2/T1)*T1) + 2*tmp170
        tmp213 = -5 + 4*LogT1 - LogT1**2 + 
     -   (2*LogT1*(B1 - T1))/T1 + 
     -   (2*B1*(LogB1 - LogT1)*(-B1 + T1))/((1 - B1/T1)*T1**2) + 
     -   (-2*B1*LogB1 + 4*T1)/T1 - 2*tmp169
        tmp214 = -5 + 4*LogT1 - LogT1**2 + 
     -   (2*LogT1*(B2 - T1))/T1 + 
     -   (2*B2*(LogB2 - LogT1)*(-B2 + T1))/((1 - B2/T1)*T1**2) + 
     -   (-2*B2*LogB2 + 4*T1)/T1 - 2*tmp170
        tmp215 = -1 + 4*LogB1 + (-2 - 2*LogB1)*LogT2 + 
     -   LogT2**2 - (2*(LogB1 - LogT2)*(-B1 + T2))/
     -    ((1 - B1/T2)*T2) + 2*tmp172
        tmp216 = -1 + 4*LogB2 + (-2 - 2*LogB2)*LogT2 + 
     -   LogT2**2 - (2*(LogB2 - LogT2)*(-B2 + T2))/
     -    ((1 - B2/T2)*T2) + 2*tmp173
        tmp217 = -5 + 4*LogT2 - LogT2**2 + 
     -   (2*LogT2*(B1 - T2))/T2 + 
     -   (2*B1*(LogB1 - LogT2)*(-B1 + T2))/((1 - B1/T2)*T2**2) + 
     -   (-2*B1*LogB1 + 4*T2)/T2 - 2*tmp172
        tmp218 = -5 + 4*LogT2 - LogT2**2 + 
     -   (2*LogT2*(B2 - T2))/T2 + 
     -   (2*B2*(LogB2 - LogT2)*(-B2 + T2))/((1 - B2/T2)*T2**2) + 
     -   (-2*B2*LogB2 + 4*T2)/T2 - 2*tmp173
        tmp219 = -1 + LogT1**2 + LogT1*(-2 - 2*LogT2) + 
     -   4*LogT2 - (2*(-LogT1 + LogT2)*(T1 - T2))/
     -    (T1*(1 - T2/T1)) + 2*tmp171
        tmp220 = -5 + 4*LogT1 - LogT1**2 + 
     -   (2*LogT1*(-T1 + T2))/T1 + (4*T1 - 2*LogT2*T2)/T1 + 
     -   (2*(-LogT1 + LogT2)*(T1 - T2)*T2)/(T1**2*(1 - T2/T1)) - 
     -   2*tmp171
        tmp221 = -0.5D0 + 2*LogT1 - 
     -   (phiA0T1T2*(-A0 + T1 - T2))/T2 - 0.5D0*(LogA0*LogT1) + 
     -   0.5D0*(LogA0*LogT2) - 0.5D0*(LogT1*LogT2) + 
     -   0.5D0*(LogT2*(A0 - T1 - T2))/T1 + 
     -   0.5D0*(LogA0*(-A0 - T1 + T2))/T1 - 
     -   0.5D0*(deltA0T1T2*((phiA0T1T2*(A0 - T1 + T2))/deltA0T1T2 + 
     -         (T2*tmp96)/(deltA0T1T2*T1)))/T2
        tmp222 = -0.5D0 + 2*LogT2 - 
     -   (phiA0B1T2*(-A0 - B1 + T2))/T2 + 0.5D0*(LogA0*LogB1) - 
     -   0.5D0*(LogA0*LogT2) - 0.5D0*(LogB1*LogT2) + 
     -   0.5D0*(LogB1*(A0 - B1 - T2))/T2 + 
     -   0.5D0*(LogA0*(-A0 + B1 - T2))/T2 - 
     -   0.5D0*(deltA0B1T2*((B1*phiA0B1T2*(A0 + B1 - T2))/
     -          (deltA0B1T2*T2) + (B1*tmp100)/(deltA0B1T2*T2)))/B1
        tmp223 = -0.5D0 + 2*LogT2 - 
     -   (phiA0B2T2*(-A0 - B2 + T2))/T2 + 0.5D0*(LogA0*LogB2) - 
     -   0.5D0*(LogA0*LogT2) - 0.5D0*(LogB2*LogT2) + 
     -   0.5D0*(LogB2*(A0 - B2 - T2))/T2 + 
     -   0.5D0*(LogA0*(-A0 + B2 - T2))/T2 - 
     -   0.5D0*(deltA0B2T2*((B2*phiA0B2T2*(A0 + B2 - T2))/
     -          (deltA0B2T2*T2) + (B2*tmp101)/(deltA0B2T2*T2)))/B2
        tmp224 = -0.5D0 + 2*LogT2 - 
     -   (phiT2bmu2*(-b - mu2 + T2))/mu2 + 0.5D0*(Logb*Logmu2) - 
     -   0.5D0*(Logb*LogT2) - 0.5D0*(Logmu2*LogT2) + 
     -   0.5D0*(Logmu2*(b - mu2 - T2))/T2 + 
     -   0.5D0*(Logb*(-b + mu2 - T2))/T2 - 
     -   0.5D0*(deltT2bmu2*((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -         (mu2*tmp102)/(deltT2bmu2*T2)))/mu2
        tmp225 = 0.5D0 - 2*LogT2 + 
     -   (phiT2bmu2*(-b - mu2 + T2))/mu2 - 0.5D0*(Logb*Logmu2) + 
     -   0.5D0*(Logb*LogT2) + 0.5D0*(Logmu2*LogT2) - 
     -   0.5D0*(Logmu2*(b - mu2 - T2))/T2 - 
     -   0.5D0*(Logb*(-b + mu2 - T2))/T2 + 
     -   0.5D0*(deltT2bmu2*((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -         (mu2*tmp102)/(deltT2bmu2*T2)))/mu2
        tmp226 = -(phiT2bmu2/deltT2bmu2) - 
     -   (2*phiT2bmu2*(b + mu2 - T2)*(-b - mu2 + T2))/
     -    deltT2bmu2**2 - (mu2*tmp102)/(deltT2bmu2*T2**2) - 
     -   (2*mu2*(-b - mu2 + T2)*tmp102)/(deltT2bmu2**2*T2) + 
     -   ((b + mu2 - T2)*((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -        (mu2*tmp102)/(deltT2bmu2*T2)))/deltT2bmu2 + 
     -   (mu2*tmp5)/(deltT2bmu2*T2)
        tmp227 = phiT2bmu2/deltT2bmu2 - 
     -   (2*phiT2bmu2*(b - mu2 - T2)*(b + mu2 - T2))/
     -    deltT2bmu2**2 - (2*mu2*(b - mu2 - T2)*tmp102)/
     -    (deltT2bmu2**2*T2) + 
     -   ((b + mu2 - T2)*((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -        (mu2*tmp110)/(b*deltT2bmu2)))/deltT2bmu2 + 
     -   (mu2*tmp97)/(deltT2bmu2*T2)
        tmp228 = -0.5D0 + 2*Logt - 
     -   (phimu2tT2*(-mu2 + t - T2))/T2 - 0.5D0*(Logmu2*Logt) + 
     -   0.5D0*(Logmu2*LogT2) - 0.5D0*(Logt*LogT2) + 
     -   0.5D0*(LogT2*(mu2 - t - T2))/t + 
     -   0.5D0*(Logmu2*(-mu2 - t + T2))/t - 
     -   0.5D0*(deltmu2tT2*((mu2*phimu2tT2*(mu2 - t + T2))/
     -          (deltmu2tT2*T2) + (mu2*tmp111)/(deltmu2tT2*t)))/mu2
        tmp229 = -0.5D0 + 2*Logt - (phiA0bt*(-A0 - b + t))/t + 
     -   0.5D0*(LogA0*Logb) - 0.5D0*(LogA0*Logt) - 0.5D0*(Logb*Logt) + 
     -   0.5D0*(Logb*(A0 - b - t))/t + 0.5D0*(LogA0*(-A0 + b - t))/t - 
     -   0.5D0*(deltA0bt*((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -         (b*tmp66)/(deltA0bt*t)))/b
        tmp230 = -((b*phiA0bt)/(deltA0bt*t)) - 
     -   (2*b*phiA0bt*(A0 + b - t)*(-A0 - b + t))/
     -    (deltA0bt**2*t) + (b*tmp1)/(deltA0bt*t) - 
     -   (b*tmp66)/(deltA0bt*t**2) - 
     -   (2*b*(-A0 - b + t)*tmp66)/(deltA0bt**2*t) + 
     -   ((A0 + b - t)*((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -        (b*tmp66)/(deltA0bt*t)))/deltA0bt
        tmp231 = -0.5D0 + 2*Logt - 
     -   (phiB1tmu2*(-B1 - mu2 + t))/mu2 + 0.5D0*(LogB1*Logmu2) - 
     -   0.5D0*(LogB1*Logt) - 0.5D0*(Logmu2*Logt) + 
     -   0.5D0*(Logmu2*(B1 - mu2 - t))/t + 
     -   0.5D0*(LogB1*(-B1 + mu2 - t))/t - 
     -   0.5D0*(deltB1tmu2*((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -         (mu2*tmp67)/(deltB1tmu2*t)))/mu2
        tmp232 = -(phiB1tmu2/deltB1tmu2) - 
     -   (2*phiB1tmu2*(B1 + mu2 - t)*(-B1 - mu2 + t))/
     -    deltB1tmu2**2 + (mu2*tmp2)/(deltB1tmu2*t) - 
     -   (mu2*tmp67)/(deltB1tmu2*t**2) - 
     -   (2*mu2*(-B1 - mu2 + t)*tmp67)/(deltB1tmu2**2*t) + 
     -   ((B1 + mu2 - t)*((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -        (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2
        tmp233 = phiB1tmu2/deltB1tmu2 - 
     -   (2*phiB1tmu2*(-B1 - mu2 + t)*(-B1 + mu2 + t))/
     -    deltB1tmu2**2 + ((-B1 + mu2 + t)*
     -      ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -        (mu2*tmp67)/(deltB1tmu2*t)))/deltB1tmu2 - 
     -   (2*mu2*(-B1 - mu2 + t)*tmp71)/(B1*deltB1tmu2**2) + 
     -   (mu2*tmp73)/(B1*deltB1tmu2)
        tmp234 = -0.5D0 + 2*Logt - 
     -   (phiB2tmu2*(-B2 - mu2 + t))/mu2 + 0.5D0*(LogB2*Logmu2) - 
     -   0.5D0*(LogB2*Logt) - 0.5D0*(Logmu2*Logt) + 
     -   0.5D0*(Logmu2*(B2 - mu2 - t))/t + 
     -   0.5D0*(LogB2*(-B2 + mu2 - t))/t - 
     -   0.5D0*(deltB2tmu2*((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -         (mu2*tmp68)/(deltB2tmu2*t)))/mu2
        tmp235 = -((phiB1tmu2*(-B1 - mu2 + t))/mu2) + 
     -   (phiB2tmu2*(-B2 - mu2 + t))/mu2 + 0.5D0*(LogB1*Logmu2) - 
     -   0.5D0*(LogB2*Logmu2) - 0.5D0*(LogB1*Logt) + 
     -   0.5D0*(LogB2*Logt) + 0.5D0*(Logmu2*(B1 - mu2 - t))/t - 
     -   0.5D0*(Logmu2*(B2 - mu2 - t))/t + 
     -   0.5D0*(LogB1*(-B1 + mu2 - t))/t - 
     -   0.5D0*(LogB2*(-B2 + mu2 - t))/t - 
     -   0.5D0*(deltB1tmu2*((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -         (mu2*tmp67)/(deltB1tmu2*t)))/mu2 + 
     -   0.5D0*(deltB2tmu2*((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -         (mu2*tmp68)/(deltB2tmu2*t)))/mu2
        tmp236 = -(phiB2tmu2/deltB2tmu2) - 
     -   (2*phiB2tmu2*(B2 + mu2 - t)*(-B2 - mu2 + t))/
     -    deltB2tmu2**2 + (mu2*tmp3)/(deltB2tmu2*t) - 
     -   (mu2*tmp68)/(deltB2tmu2*t**2) - 
     -   (2*mu2*(-B2 - mu2 + t)*tmp68)/(deltB2tmu2**2*t) + 
     -   ((B2 + mu2 - t)*((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -        (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2
        tmp237 = phiB2tmu2/deltB2tmu2 - 
     -   (2*phiB2tmu2*(-B2 - mu2 + t)*(-B2 + mu2 + t))/
     -    deltB2tmu2**2 + ((-B2 + mu2 + t)*
     -      ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -        (mu2*tmp68)/(deltB2tmu2*t)))/deltB2tmu2 - 
     -   (2*mu2*(-B2 - mu2 + t)*tmp72)/(B2*deltB2tmu2**2) + 
     -   (mu2*tmp74)/(B2*deltB2tmu2)
        tmp238 = -0.5D0 + 2*LogT1 - 
     -   (phiA0B1T1*(-A0 - B1 + T1))/T1 + 0.5D0*(LogA0*LogB1) - 
     -   0.5D0*(LogA0*LogT1) - 0.5D0*(LogB1*LogT1) + 
     -   0.5D0*(LogB1*(A0 - B1 - T1))/T1 + 
     -   0.5D0*(LogA0*(-A0 + B1 - T1))/T1 - 
     -   0.5D0*(deltA0B1T1*((B1*phiA0B1T1*(A0 + B1 - T1))/
     -          (deltA0B1T1*T1) + (B1*tmp81)/(deltA0B1T1*T1)))/B1
        tmp239 = -0.5D0 + 2*LogT1 - 
     -   (phiA0B2T1*(-A0 - B2 + T1))/T1 + 0.5D0*(LogA0*LogB2) - 
     -   0.5D0*(LogA0*LogT1) - 0.5D0*(LogB2*LogT1) + 
     -   0.5D0*(LogB2*(A0 - B2 - T1))/T1 + 
     -   0.5D0*(LogA0*(-A0 + B2 - T1))/T1 - 
     -   0.5D0*(deltA0B2T1*((B2*phiA0B2T1*(A0 + B2 - T1))/
     -          (deltA0B2T1*T1) + (B2*tmp82)/(deltA0B2T1*T1)))/B2
        tmp240 = -0.5D0 + 2*LogT1 - 
     -   (phiT1bmu2*(-b - mu2 + T1))/mu2 + 0.5D0*(Logb*Logmu2) - 
     -   0.5D0*(Logb*LogT1) - 0.5D0*(Logmu2*LogT1) + 
     -   0.5D0*(Logmu2*(b - mu2 - T1))/T1 + 
     -   0.5D0*(Logb*(-b + mu2 - T1))/T1 - 
     -   0.5D0*(deltT1bmu2*((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -         (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
        tmp241 = -(phiT1bmu2/deltT1bmu2) - 
     -   (2*phiT1bmu2*(b + mu2 - T1)*(-b - mu2 + T1))/
     -    deltT1bmu2**2 + (mu2*tmp4)/(deltT1bmu2*T1) - 
     -   (mu2*tmp83)/(deltT1bmu2*T1**2) - 
     -   (2*mu2*(-b - mu2 + T1)*tmp83)/(deltT1bmu2**2*T1) + 
     -   ((b + mu2 - T1)*((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -        (mu2*tmp83)/(deltT1bmu2*T1)))/deltT1bmu2
        tmp242 = phiT1bmu2/deltT1bmu2 - 
     -   (2*phiT1bmu2*(b - mu2 - T1)*(b + mu2 - T1))/
     -    deltT1bmu2**2 + (mu2*tmp78)/(deltT1bmu2*T1) - 
     -   (2*mu2*(b - mu2 - T1)*tmp83)/(deltT1bmu2**2*T1) + 
     -   ((b + mu2 - T1)*((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -        (mu2*tmp90)/(b*deltT1bmu2)))/deltT1bmu2
        tmp243 = -0.5D0 + 2*Logt - 
     -   (phimu2tT1*(-mu2 + t - T1))/T1 - 0.5D0*(Logmu2*Logt) + 
     -   0.5D0*(Logmu2*LogT1) - 0.5D0*(Logt*LogT1) + 
     -   0.5D0*(LogT1*(mu2 - t - T1))/t + 
     -   0.5D0*(Logmu2*(-mu2 - t + T1))/t - 
     -   0.5D0*(deltmu2tT1*((mu2*phimu2tT1*(mu2 - t + T1))/
     -          (deltmu2tT1*T1) + (mu2*tmp91)/(deltmu2tT1*t)))/mu2
        tmp244 = 2/T1 - Logb/T1 - Logmu2/T1 - 
     -   0.5D0*(Logmu2*(b - mu2 - T1))/T1**2 - 
     -   0.5D0*(Logb*(-b + mu2 - T1))/T1**2 - 
     -   0.5D0*(2*phiT1bmu2 + deltT1bmu2*tmp241 + 
     -       4*(-b - mu2 + T1)*
     -        ((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -          (mu2*tmp83)/(deltT1bmu2*T1)))/mu2
        D1t = -(tmp109*tmp195) - tmp108*tmp196 - 
     -   2*hb*ht*mt*mu*s2b*tmp235 - (hb*ht*mu*s2b*tmp75)/mt + 
     -   tmp61*(-(B1*(-1 + LogB1)) + 2*B1*LogB1 - 
     -      B1*(-1 + LogB1)*(-1 + Logt) + (-1 + Logmu2)*mu2 + 
     -      2*Logmu2*mu2 + (-1 + Logmu2)*(-1 + Logt)*mu2 + 
     -      2*Logt*t - (B1 - mu2 - t)*tmp231 - 
     -      0.5D0*(deltB1tmu2*phiB1tmu2)/mu2 + 
     -      0.5D0*(Logmu2*Logt*(B1 - mu2 - t)) + 
     -      0.5D0*(LogB1*Logt*(-B1 + mu2 - t)) + 
     -      0.5D0*(LogB1*Logmu2*(-B1 - mu2 + t)) - 2.5D0*(B1 + mu2 + t)
     -      ) + tmp60*(-(B2*(-1 + LogB2)) + 2*B2*LogB2 - 
     -      B2*(-1 + LogB2)*(-1 + Logt) + (-1 + Logmu2)*mu2 + 
     -      2*Logmu2*mu2 + (-1 + Logmu2)*(-1 + Logt)*mu2 + 
     -      2*Logt*t - (B2 - mu2 - t)*tmp234 - 
     -      0.5D0*(deltB2tmu2*phiB2tmu2)/mu2 + 
     -      0.5D0*(Logmu2*Logt*(B2 - mu2 - t)) + 
     -      0.5D0*(LogB2*Logt*(-B2 + mu2 - t)) + 
     -      0.5D0*(LogB2*Logmu2*(-B2 - mu2 + t)) - 2.5D0*(B2 + mu2 + t)
     -      ) - 0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 
     -        5*T1 + LogT1*(-2*B2*LogB2 + 4*T1) - 
     -        2*(-B2 + T1)*tmp170)*tmp183) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      tmp184) - 0.5D0*((-5*B2 + 4*B2*LogB2 + 
     -        LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      tmp187) - 0.5D0*((-5*B1 + 4*B1*LogB1 + 
     -        LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      tmp188) - 4*cbe*hb*ht*mb*mt*sbe*
     -    (0.5D0 - 2*Logt + (phiA0bt*(-A0 - b + t))/t - 
     -      0.5D0*(LogA0*Logb) + 0.5D0*(LogA0*Logt) + 
     -      0.5D0*(Logb*Logt) - 0.5D0*(Logb*(A0 - b - t))/t - 
     -      0.5D0*(LogA0*(-A0 + b - t))/t + 0.5D0*tmp210 + 
     -      0.5D0*(deltA0bt*((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -            (b*tmp66)/(deltA0bt*t)))/b) - 
     -   (2*cbe*hb*ht*mb*sbe*
     -      (-2*A0*LogA0 - 2*b*Logb - 2*Logt*t - 
     -        0.5D0*(Logb*Logt*(A0 - b - t)) - 
     -        0.5D0*(LogA0*Logt*(-A0 + b - t)) + 
     -        0.5D0*(deltA0bt*phiA0bt)/t - 
     -        0.5D0*(LogA0*Logb*(-A0 - b + t)) + 2.5D0*(A0 + b + t) + 
     -        0.5D0*tmp76))/mt + 
     -   tmp10*(b*(-1 + Logb) + b*(-1 + Logb)*(-1 + Logt) + 
     -      0.5D0*((b + t)*tmp210) + 0.5D0*tmp76)
        D1t = D1t - tmp192*tmp88 - tmp191*tmp89 + 
     -   tmp9*(-(A0*(-1 + LogA0)) + 2*A0*LogA0 + b*(-1 + Logb) + 
     -      2*b*Logb - A0*(-1 + LogA0)*(-1 + Logt) + 
     -      b*(-1 + Logb)*(-1 + Logt) + 2*Logt*t - 
     -      (A0 - b - t)*tmp229 + 0.5D0*(Logb*Logt*(A0 - b - t)) + 
     -      0.5D0*(LogA0*Logt*(-A0 + b - t)) - 
     -      0.5D0*(deltA0bt*phiA0bt)/t + 
     -      0.5D0*(LogA0*Logb*(-A0 - b + t)) - 2.5D0*(A0 + b + t)) + 
     -   ht**2*(2*(-1 + Logmu2)*mu2 + 4*Logmu2*mu2 + 
     -      2*(-1 + Logmu2)*(-1 + Logt)*mu2 + 4*Logt*t - 
     -      (-1 + LogT1)*T1 - (-1 + Logt)*(-1 + LogT1)*T1 + 
     -      2*LogT1*T1 - (-1 + LogT2)*T2 - 
     -      (-1 + Logt)*(-1 + LogT2)*T2 + 2*LogT2*T2 - 
     -      (-mu2 - t + T2)*tmp228 - (-mu2 - t + T1)*tmp243 + 
     -      sbe**2*(-10*t + 2*(-1 + Logt)*t + 2*(-1 + Logt)**2*t + 
     -         4*Logt*t + Logt*(4*t - 2*Logt*t) + t*tmp77) + 
     -      0.5D0*(Logt*LogT1*(mu2 - t - T1)) + 
     -      0.5D0*(Logmu2*LogT1*(-mu2 + t - T1)) - 
     -      0.5D0*(deltmu2tT1*phimu2tT1)/T1 + 
     -      0.5D0*(Logmu2*Logt*(-mu2 - t + T1)) - 
     -      2.5D0*(mu2 + t + T1) + 0.5D0*(Logt*LogT2*(mu2 - t - T2)) + 
     -      0.5D0*(Logmu2*LogT2*(-mu2 + t - T2)) - 
     -      0.5D0*(deltmu2tT2*phimu2tT2)/T2 + 
     -      0.5D0*(Logmu2*Logt*(-mu2 - t + T2)) - 
     -      2.5D0*(mu2 + t + T2) + 
     -      cbe**2*(-2*A0*(-1 + LogA0) - 
     -         2*A0*(-1 + LogA0)*(-1 + Logt) + 2*(-1 + Logt)*t + 
     -         2*(-1 + Logt)**2*t + 
     -         2*(2*A0*LogA0 - A0*LogA0*Logt + 4*Logt*t + 
     -            0.5D0*(Logt**2*(A0 - 2*t)) - 
     -            0.5D0*(deltA0tt*phiA0tt)/t - 2.5D0*(A0 + 2*t)) - 
     -         (A0 - 2*t)*(-1 + 4*Logt - Logt**2 - (A0*LogA0)/t + 
     -            (2*A0*phiA0tt)/t + (Logt*(A0 - 2*t))/t + 
     -            0.5D0*(deltA0tt*phiA0tt)/t**2 - 
     -            0.5D0*(deltA0tt*
     -                ((A0*phiA0tt)/deltA0tt + 
     -                  (phiA0tt*tmp175)/deltA0tt + 
     -                  tmp65/deltA0tt + tmp70/deltA0tt))/t))) - 
     -   0.25D0*(ht**2*(cbe**2*tmp114*(4 - (2*s2t*Yt)/mt) + 
     -        cbe**2*tmp93*(4 + (2*s2t*Yt)/mt) + 
     -        0.5D0*(sbe**2*tmp115*(4 - (2*s2t*Xt)/mt)) + 
     -        0.5D0*(sbe**2*tmp94*(4 + (2*s2t*Xt)/mt))))
        DT1 = -(tmp194*tmp238) - tmp193*tmp239 - 
     -   2*hb*ht*mb*mu*s2t*tmp240 + 
     -   hb**2*(0.25D0*(B1*(1 - c2b)*(1 + c2t)*(-1 + LogB1)) + 
     -      0.25D0*(B2*(1 + c2b)*(1 + c2t)*(-1 + LogB2)) + 
     -      sbe**2*(0.5D0*(A0*(1 + c2t)*(-1 + LogA0)) + 
     -         0.5D0*(A0*(1 + c2t)*(-1 + LogA0)*(-1 + LogT1))) + 
     -      0.25D0*(B1*(1 - c2b)*(1 + c2t)*(-1 + LogB1)*
     -         (-1 + LogT1)) + 
     -      0.25D0*(B2*(1 + c2b)*(1 + c2t)*(-1 + LogB2)*(-1 + LogT1))
     -      ) + tmp62*(-(b*(-1 + Logb)) - 2*b*Logb - 
     -      b*(-1 + Logb)*(-1 + LogT1) - (-1 + Logmu2)*mu2 - 
     -      2*Logmu2*mu2 - (-1 + Logmu2)*(-1 + LogT1)*mu2 - 
     -      2*LogT1*T1 - (-b - mu2 + T1)*tmp240 + 
     -      0.5D0*(deltT1bmu2*phiT1bmu2)/mu2 - 
     -      0.5D0*(Logmu2*LogT1*(b - mu2 - T1)) - 
     -      0.5D0*(Logb*LogT1*(-b + mu2 - T1)) - 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T1)) + 2.5D0*(b + mu2 + T1))
     -     - 0.5D0*(tmp186*tmp213) - 0.5D0*(tmp185*tmp214) + 
     -   ht**2*((-1 + LogT2)*T2*tmp8 + 
     -      (-1 + LogT1)*(-1 + LogT2)*T2*tmp8 + 
     -      cbe**2*(A0*(-1 + LogA0)*(1 + 0.5D0*(1 - c2t)) + 
     -         A0*(-1 + LogA0)*(-1 + LogT1)*(1 + 0.5D0*(1 - c2t))) + 
     -      0.25D0*(B1*(1 + c2b)*(1 - c2t)*(-1 + LogB1)) + 
     -      0.25D0*(B2*(1 - c2b)*(1 - c2t)*(-1 + LogB2)) + 
     -      0.25D0*(B1*(1 + c2b)*(1 - c2t)*(-1 + LogB1)*
     -         (-1 + LogT1)) + 
     -      0.25D0*(B2*(1 - c2b)*(1 - c2t)*(-1 + LogB2)*
     -         (-1 + LogT1)) + 0.25D0*((1 + Nc)*s2t**2*tmp79)) + 
     -   ht**2*(-((-1 + Logmu2)*mu2) - 2*Logmu2*mu2 - 
     -      (-1 + Logmu2)*(-1 + LogT1)*mu2 - (-1 + Logt)*t - 
     -      2*Logt*t - (-1 + Logt)*(-1 + LogT1)*t - 2*LogT1*T1 - 
     -      0.5D0*(Logt*LogT1*(mu2 - t - T1)) - 
     -      0.5D0*(Logmu2*LogT1*(-mu2 + t - T1)) + 
     -      0.5D0*(deltmu2tT1*phimu2tT1)/T1 - 
     -      0.5D0*(Logmu2*Logt*(-mu2 - t + T1)) + 
     -      2.5D0*(mu2 + t + T1) - 
     -      (-mu2 - t + T1)*
     -       (-0.5D0 + 2*LogT1 - (phimu2tT1*(-mu2 - t + T1))/T1 + 
     -         0.5D0*(Logmu2*Logt) - 0.5D0*(Logmu2*LogT1) - 
     -         0.5D0*(Logt*LogT1) + 0.5D0*(Logt*(mu2 - t - T1))/T1 + 
     -         0.5D0*(Logmu2*(-mu2 + t - T1))/T1 - 
     -         0.5D0*(deltmu2tT1*
     -             ((mu2*phimu2tT1*(mu2 + t - T1))/
     -                (deltmu2tT1*T1) + (mu2*tmp84)/(deltmu2tT1*T1)
     -               ))/mu2))
        DT1 = DT1 - 0.25D0*
     -    (ht**2*((1 + c2t**2)*sbe**2*tmp220*Xt**2 + 
     -        2*(1 + c2t**2)*cbe**2*tmp221*Yt**2 + 
     -        cbe**2*tmp51*
     -         (-1 + 4*LogT1 - LogT1**2 - (A0*LogA0)/T1 + 
     -           (2*A0*phiA0T1T1)/T1 + (LogT1*(A0 - 2*T1))/T1 + 
     -           0.5D0*(deltA0T1T1*phiA0T1T1)/T1**2 - 
     -           0.5D0*(deltA0T1T1*
     -               ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                 (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                 tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/T1) + 
     -        0.5D0*(sbe**2*tmp35*tmp95)))
        DT2 = -(tmp198*tmp222) - tmp197*tmp223 - 
     -   2*hb*ht*mb*mu*s2t*tmp225 + 
     -   hb**2*(0.25D0*(B1*(1 - c2b)*(1 - c2t)*(-1 + LogB1)) + 
     -      0.25D0*(B2*(1 + c2b)*(1 - c2t)*(-1 + LogB2)) + 
     -      sbe**2*(0.5D0*(A0*(1 - c2t)*(-1 + LogA0)) + 
     -         0.5D0*(A0*(1 - c2t)*(-1 + LogA0)*(-1 + LogT2))) + 
     -      0.25D0*(B1*(1 - c2b)*(1 - c2t)*(-1 + LogB1)*
     -         (-1 + LogT2)) + 
     -      0.25D0*(B2*(1 + c2b)*(1 - c2t)*(-1 + LogB2)*(-1 + LogT2))
     -      ) + tmp63*(-(b*(-1 + Logb)) - 2*b*Logb - 
     -      b*(-1 + Logb)*(-1 + LogT2) - (-1 + Logmu2)*mu2 - 
     -      2*Logmu2*mu2 - (-1 + Logmu2)*(-1 + LogT2)*mu2 - 
     -      2*LogT2*T2 - (-b - mu2 + T2)*tmp224 + 
     -      0.5D0*(deltT2bmu2*phiT2bmu2)/mu2 - 
     -      0.5D0*(Logmu2*LogT2*(b - mu2 - T2)) - 
     -      0.5D0*(Logb*LogT2*(-b + mu2 - T2)) - 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T2)) + 2.5D0*(b + mu2 + T2))
     -     + ht**2*(-((-1 + Logmu2)*mu2) - 2*Logmu2*mu2 - 
     -      (-1 + Logmu2)*(-1 + LogT2)*mu2 - (-1 + Logt)*t - 
     -      2*Logt*t - (-1 + Logt)*(-1 + LogT2)*t - 2*LogT2*T2 - 
     -      0.5D0*(Logt*LogT2*(mu2 - t - T2)) - 
     -      0.5D0*(Logmu2*LogT2*(-mu2 + t - T2)) + 
     -      0.5D0*(deltmu2tT2*phimu2tT2)/T2 - 
     -      0.5D0*(Logmu2*Logt*(-mu2 - t + T2)) + 
     -      2.5D0*(mu2 + t + T2) - 
     -      (-mu2 - t + T2)*
     -       (-0.5D0 + 2*LogT2 - (phimu2tT2*(-mu2 - t + T2))/T2 + 
     -         0.5D0*(Logmu2*Logt) - 0.5D0*(Logmu2*LogT2) - 
     -         0.5D0*(Logt*LogT2) + 0.5D0*(Logt*(mu2 - t - T2))/T2 + 
     -         0.5D0*(Logmu2*(-mu2 + t - T2))/T2 - 
     -         0.5D0*(deltmu2tT2*
     -             ((mu2*phimu2tT2*(mu2 + t - T2))/
     -                (deltmu2tT2*T2) + 
     -               (mu2*tmp103)/(deltmu2tT2*T2)))/mu2)) - 
     -   0.5D0*(tmp190*tmp217) - 0.5D0*(tmp189*tmp218) + 
     -   ht**2*((-1 + LogT1)*T1*tmp8 + 
     -      (-1 + LogT1)*(-1 + LogT2)*T1*tmp8 + 
     -      cbe**2*(A0*(-1 + LogA0)*(1 + 0.5D0*(1 + c2t)) + 
     -         A0*(-1 + LogA0)*(-1 + LogT2)*(1 + 0.5D0*(1 + c2t))) + 
     -      0.25D0*(B1*(1 + c2b)*(1 + c2t)*(-1 + LogB1)) + 
     -      0.25D0*(B2*(1 - c2b)*(1 + c2t)*(-1 + LogB2)) + 
     -      0.25D0*(B1*(1 + c2b)*(1 + c2t)*(-1 + LogB1)*
     -         (-1 + LogT2)) + 
     -      0.25D0*(B2*(1 - c2b)*(1 + c2t)*(-1 + LogB2)*
     -         (-1 + LogT2)) + 0.25D0*((1 + Nc)*s2t**2*tmp98))
        DT2 = DT2 - 0.25D0*
     -    (ht**2*((1 + c2t**2)*sbe**2*tmp219*Xt**2 + 
     -        2*(1 + c2t**2)*cbe**2*Yt**2*
     -         (-0.5D0 + 2*LogT2 - (phiA0T1T2*(-A0 - T1 + T2))/T2 + 
     -           0.5D0*(LogA0*LogT1) - 0.5D0*(LogA0*LogT2) - 
     -           0.5D0*(LogT1*LogT2) + 
     -           0.5D0*(deltA0T1T2*phiA0T1T2)/T2**2 + 
     -           0.5D0*(LogT1*(A0 - T1 - T2))/T2 + 
     -           0.5D0*(LogA0*(-A0 + T1 - T2))/T2 - 
     -           0.5D0*(deltA0T1T2*
     -               (tmp104/deltA0T1T2 + 
     -                 (phiA0T1T2*tmp179)/deltA0T1T2))/T2) + 
     -        0.5D0*(sbe**2*tmp116*tmp34) + 
     -        cbe**2*tmp50*
     -         (-1 + 4*LogT2 - LogT2**2 - (A0*LogA0)/T2 + 
     -           (2*A0*phiA0T2T2)/T2 + (LogT2*(A0 - 2*T2))/T2 + 
     -           0.5D0*(deltA0T2T2*phiA0T2T2)/T2**2 - 
     -           0.5D0*(deltA0T2T2*
     -               ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                 tmp107/deltA0T2T2 + 
     -                 (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                 tmp99/deltA0T2T2))/T2)))
        Dc2t = (hb*ht*mb*mu*tmp113)/s2t - tmp109*tmp205 - 
     -   tmp108*tmp206 + (b*(-1 + Logb)*(-1 + Logmu2)*mu2 - 
     -      b*(-1 + Logb)*(-1 + LogT2)*T2 - 
     -      (-1 + Logmu2)*(-1 + LogT2)*mu2*T2 - 
     -      (-b - mu2 + T2)*tmp112)*tmp7 - tmp204*tmp88 - 
     -   tmp203*tmp89 + tmp6*
     -    (b*(-1 + Logb)*(-1 + Logmu2)*mu2 - 
     -      b*(-1 + Logb)*(-1 + LogT1)*T1 - 
     -      (-1 + Logmu2)*(-1 + LogT1)*mu2*T1 - 
     -      (-b - mu2 + T1)*tmp92) + 
     -   hb**2*(0.125D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT1)*T1)/
     -        c2t + 0.125D0*(B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT1)*
     -          T1)/c2t + sbe**2*
     -       (0.25D0*(A0*(-1 + LogA0)*(-1 + LogT1)*T1)/c2t - 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT2)*T2)/c2t) - 
     -      0.125D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/
     -        c2t - 0.125D0*(B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT2)*
     -          T2)/c2t) + 
     -   ht**2*((-1 + LogT1)*(-1 + LogT2)*T1*T2*
     -       (1 + 0.5D0*(-1 + Nc)) - 
     -      0.125D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT1)*T1)/
     -        c2t - 0.125D0*(B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT1)*
     -          T1)/c2t + cbe**2*
     -       (-(0.25D0*(A0*(-1 + LogA0)*(-1 + LogT1)*T1)/c2t) + 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT2)*T2)/c2t) + 
     -      0.125D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/
     -        c2t + 0.125D0*(B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT2)*
     -          T2)/c2t - 0.25D0*
     -       ((1 + Nc)*((-1 + LogT1)**2*T1**2 + 
     -           (-1 + LogT2)**2*T2**2))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      tmp199) - 0.5D0*((-5*B1 + 4*B1*LogB1 + 
     -        LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      tmp200) - 0.5D0*((-5*B2 + 4*B2*LogB2 + 
     -        LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      tmp201) - 0.5D0*((-5*B1 + 4*B1*LogB1 + 
     -        LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      tmp202)
        Dc2t = Dc2t - 0.25D0*
     -    (ht**2*(cbe**2*tmp114*tmp47 + cbe**2*tmp46*tmp93 + 
     -        sbe**2*(-5*T1 - 5*T2 + 4*LogT2*T2 + 
     -           LogT1**2*(-T1 + T2) + LogT1*(4*T1 - 2*LogT2*T2) - 
     -           2*(T1 - T2)*tmp171)*Xt**2 + 
     -        2*cbe**2*Yt**2*
     -         (2*A0*LogA0 + 2*LogT1*T1 + 2*LogT2*T2 + 
     -           0.5D0*(LogT1*LogT2*(A0 - T1 - T2)) + 
     -           0.5D0*(LogA0*LogT2*(-A0 + T1 - T2)) - 
     -           0.5D0*(deltA0T1T2*phiA0T1T2)/T2 + 
     -           0.5D0*(LogA0*LogT1*(-A0 - T1 + T2)) - 
     -           2.5D0*(A0 + T1 + T2)) + 0.5D0*(sbe**2*tmp115*tmp31) + 
     -        0.5D0*(sbe**2*tmp30*tmp94)))
        DT1T1 = -2*hb*ht*mb*mu*s2t*tmp244 + 
     -   hb**2*(0.25D0*(B1*(1 - c2b)*(1 + c2t)*(-1 + LogB1))/T1 + 
     -      0.25D0*(B2*(1 + c2b)*(1 + c2t)*(-1 + LogB2))/T1 + 
     -      0.5D0*(A0*(1 + c2t)*(-1 + LogA0)*sbe**2)/T1) + 
     -   ht**2*(((-1 + LogT2)*T2*tmp8)/T1 + 
     -      (A0*cbe**2*(-1 + LogA0)*(1 + 0.5D0*(1 - c2t)))/T1 + 
     -      0.25D0*(B1*(1 + c2b)*(1 - c2t)*(-1 + LogB1))/T1 + 
     -      0.25D0*(B2*(1 - c2b)*(1 - c2t)*(-1 + LogB2))/T1 + 
     -      0.25D0*((1 + Nc)*s2t**2*
     -         (8*(-1 + LogT1) + 2*(-1 + LogT1)**2 + 
     -           (2/T1**2 - (2*(-1 + LogT1))/T1**2)*T1**2))) - 
     -   0.5D0*(((4*B2*(LogB2 - LogT1))/((1 - B2/T1)*T1**2) + 8/T1 - 
     -        (4*LogT1)/T1 - 
     -        2*((B2**2*(LogB2 - LogT1))/((1 - B2/T1)**2*T1**4) + 
     -           B2/((1 - B2/T1)*T1**3) + 
     -           (2*B2*(LogB2 - LogT1))/((1 - B2/T1)*T1**3))*
     -         (-B2 + T1) - (-2*B2*LogB2 + 4*T1)/T1**2 + 
     -        (B2 - T1)*tmp17)*tmp185) - 
     -   0.5D0*(((4*B1*(LogB1 - LogT1))/((1 - B1/T1)*T1**2) + 8/T1 - 
     -        (4*LogT1)/T1 - 
     -        2*((B1**2*(LogB1 - LogT1))/((1 - B1/T1)**2*T1**4) + 
     -           B1/((1 - B1/T1)*T1**3) + 
     -           (2*B1*(LogB1 - LogT1))/((1 - B1/T1)*T1**3))*
     -         (-B1 + T1) - (-2*B1*LogB1 + 4*T1)/T1**2 + 
     -        (B1 - T1)*tmp17)*tmp186) - 
     -   tmp194*(2/T1 - LogA0/T1 - LogB1/T1 - 
     -      0.5D0*(LogB1*(A0 - B1 - T1))/T1**2 - 
     -      0.5D0*(LogA0*(-A0 + B1 - T1))/T1**2 - 
     -      0.5D0*((2*B1*phiA0B1T1)/T1 + 
     -          4*(-A0 - B1 + T1)*
     -           ((B1*phiA0B1T1*(A0 + B1 - T1))/(deltA0B1T1*T1) + 
     -             (B1*tmp81)/(deltA0B1T1*T1)) + 
     -          deltA0B1T1*
     -           ((B1*(2 - LogA0 - LogB1 + 2*LogT1))/
     -              (deltA0B1T1*T1) - 
     -             (B1*phiA0B1T1)/(deltA0B1T1*T1) - 
     -             (2*B1*phiA0B1T1*(A0 + B1 - T1)*(-A0 - B1 + T1))/
     -              (deltA0B1T1**2*T1) - 
     -             (B1*tmp81)/(deltA0B1T1*T1**2) - 
     -             (2*B1*(-A0 - B1 + T1)*tmp81)/
     -              (deltA0B1T1**2*T1) + 
     -             ((A0 + B1 - T1)*
     -                ((B1*phiA0B1T1*(A0 + B1 - T1))/
     -                   (deltA0B1T1*T1) + 
     -                  (B1*tmp81)/(deltA0B1T1*T1)))/deltA0B1T1))/
     -        B1)
        DT1T1 = DT1T1 - tmp193*
     -    (2/T1 - LogA0/T1 - LogB2/T1 - 
     -      0.5D0*(LogB2*(A0 - B2 - T1))/T1**2 - 
     -      0.5D0*(LogA0*(-A0 + B2 - T1))/T1**2 - 
     -      0.5D0*((2*B2*phiA0B2T1)/T1 + 
     -          4*(-A0 - B2 + T1)*
     -           ((B2*phiA0B2T1*(A0 + B2 - T1))/(deltA0B2T1*T1) + 
     -             (B2*tmp82)/(deltA0B2T1*T1)) + 
     -          deltA0B2T1*
     -           ((B2*(2 - LogA0 - LogB2 + 2*LogT1))/
     -              (deltA0B2T1*T1) - 
     -             (B2*phiA0B2T1)/(deltA0B2T1*T1) - 
     -             (2*B2*phiA0B2T1*(A0 + B2 - T1)*(-A0 - B2 + T1))/
     -              (deltA0B2T1**2*T1) - 
     -             (B2*tmp82)/(deltA0B2T1*T1**2) - 
     -             (2*B2*(-A0 - B2 + T1)*tmp82)/
     -              (deltA0B2T1**2*T1) + 
     -             ((A0 + B2 - T1)*
     -                ((B2*phiA0B2T1*(A0 + B2 - T1))/
     -                   (deltA0B2T1*T1) + 
     -                  (B2*tmp82)/(deltA0B2T1*T1)))/deltA0B2T1))/
     -        B2) + tmp62*(-((b*(-1 + Logb))/T1) - 
     -      ((-1 + Logmu2)*mu2)/T1 - (-b - mu2 + T1)*tmp244 + 
     -      2*(0.5D0 - 2*LogT1 + (phiT1bmu2*(-b - mu2 + T1))/mu2 - 
     -         0.5D0*(Logb*Logmu2) + 0.5D0*(Logb*LogT1) + 
     -         0.5D0*(Logmu2*LogT1) - 
     -         0.5D0*(Logmu2*(b - mu2 - T1))/T1 - 
     -         0.5D0*(Logb*(-b + mu2 - T1))/T1 + 
     -         0.5D0*(deltT1bmu2*
     -             ((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -               (mu2*tmp83)/(deltT1bmu2*T1)))/mu2))
        DT1T1 = DT1T1 + ht**2*
     -    (-(((-1 + Logmu2)*mu2)/T1) - ((-1 + Logt)*t)/T1 + 
     -      2*(0.5D0 - 2*LogT1 + (phimu2tT1*(-mu2 - t + T1))/T1 - 
     -         0.5D0*(Logmu2*Logt) + 0.5D0*(Logmu2*LogT1) + 
     -         0.5D0*(Logt*LogT1) - 0.5D0*(Logt*(mu2 - t - T1))/T1 - 
     -         0.5D0*(Logmu2*(-mu2 + t - T1))/T1 + 
     -         0.5D0*(deltmu2tT1*
     -             ((mu2*phimu2tT1*(mu2 + t - T1))/
     -                (deltmu2tT1*T1) + (mu2*tmp84)/(deltmu2tT1*T1)
     -               ))/mu2) - 
     -      (-mu2 - t + T1)*
     -       (2/T1 - Logmu2/T1 - Logt/T1 - 
     -         0.5D0*(Logt*(mu2 - t - T1))/T1**2 - 
     -         0.5D0*(Logmu2*(-mu2 + t - T1))/T1**2 - 
     -         0.5D0*((2*mu2*phimu2tT1)/T1 + 
     -             4*(-mu2 - t + T1)*
     -              ((mu2*phimu2tT1*(mu2 + t - T1))/
     -                 (deltmu2tT1*T1) + 
     -                (mu2*tmp84)/(deltmu2tT1*T1)) + 
     -             deltmu2tT1*
     -              (((2 - Logmu2 - Logt + 2*LogT1)*mu2)/
     -                 (deltmu2tT1*T1) - 
     -                (mu2*phimu2tT1)/(deltmu2tT1*T1) - 
     -                (2*mu2*phimu2tT1*(mu2 + t - T1)*
     -                   (-mu2 - t + T1))/(deltmu2tT1**2*T1) - 
     -                (mu2*tmp84)/(deltmu2tT1*T1**2) - 
     -                (2*mu2*(-mu2 - t + T1)*tmp84)/
     -                 (deltmu2tT1**2*T1) + 
     -                ((mu2 + t - T1)*
     -                   ((mu2*phimu2tT1*(mu2 + t - T1))/
     -                      (deltmu2tT1*T1) + 
     -                     (mu2*tmp84)/(deltmu2tT1*T1)))/deltmu2tT1
     -                ))/mu2))
        DT1T1 = DT1T1 - 0.25D0*
     -    (ht**2*((1 + c2t**2)*sbe**2*
     -         (8/T1 - (4*LogT1)/T1 - (4*T1 - 2*LogT2*T2)/T1**2 + 
     -           (4*(-LogT1 + LogT2)*T2)/(T1**2*(1 - T2/T1)) - 
     -           2*(T1 - T2)*
     -            (((-LogT1 + LogT2)*T2**2)/
     -               (T1**4*(1 - T2/T1)**2) + 
     -              T2/(T1**3*(1 - T2/T1)) + 
     -              (2*(-LogT1 + LogT2)*T2)/(T1**3*(1 - T2/T1))) + 
     -           (-T1 + T2)*tmp17)*Xt**2 + 
     -        0.5D0*(sbe**2*(4/T1 + (2*(2 - 2*LogT1))/T1 - 
     -             (2*LogT1)/T1 - (4*T1 - 2*LogT1*T1)/T1**2)*tmp35)
     -          + cbe**2*tmp51*
     -         (-((deltA0T1T1*phiA0T1T1)/T1**3) + 
     -           (A0*LogA0)/T1**2 + 4/T1 - (4*LogT1)/T1 + 
     -           (-4*A0*phiA0T1T1 + 
     -              deltA0T1T1*
     -               ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                 (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                 tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/T1**2+
     -             0.5D0*((A0 - 2*T1)*tmp17) - 
     -           0.5D0*(-8*A0*
     -                ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                  (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                  tmp80/deltA0T1T1 + tmp87/deltA0T1T1) + 
     -               deltA0T1T1*
     -                ((4*A0**2*phiA0T1T1)/deltA0T1T1**2 + 
     -                  (phiA0T1T1*
     -                     (-1 - (A0 - T1)**2/T1**2 - 
     -                       (2*(A0 - T1))/T1))/deltA0T1T1 + 
     -                  (1 + (A0 - T1)/T1)/deltA0T1T1 + 
     -                  (1 - (-A0 + T1)/T1)/deltA0T1T1 + 
     -                  (4*A0*phiA0T1T1*tmp178)/deltA0T1T1**2 + 
     -                  (4*A0*tmp80)/deltA0T1T1**2 + 
     -                  (4*A0*tmp87)/deltA0T1T1**2 + 
     -                  (A0*
     -                     ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                       (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                       tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/
     -                   deltA0T1T1 + 
     -                  (tmp178*
     -                     ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                       (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                       tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/
     -                   deltA0T1T1))/T1) + 
     -        2*(1 + c2t**2)*cbe**2*Yt**2*
     -         (2/T1 - LogA0/T1 - LogT2/T1 - 
     -           0.5D0*(LogT2*(A0 - T1 - T2))/T1**2 - 
     -           0.5D0*(LogA0*(-A0 - T1 + T2))/T1**2 - 
     -           0.5D0*(2*phiA0T1T2 + 
     -               4*(-A0 + T1 - T2)*
     -                ((phiA0T1T2*(A0 - T1 + T2))/deltA0T1T2 + 
     -                  (T2*tmp96)/(deltA0T1T2*T1)) + 
     -               deltA0T1T2*
     -                (-(phiA0T1T2/deltA0T1T2) + 
     -                  ((2 - LogA0 + 2*LogT1 - LogT2)*T2)/
     -                   (deltA0T1T2*T1) - 
     -                  (2*phiA0T1T2*(-A0 + T1 - T2)*
     -                     (A0 - T1 + T2))/deltA0T1T2**2 - 
     -                  (T2*tmp96)/(deltA0T1T2*T1**2) - 
     -                  (2*(-A0 + T1 - T2)*T2*tmp96)/
     -                   (deltA0T1T2**2*T1) + 
     -                  ((A0 - T1 + T2)*
     -                     ((phiA0T1T2*(A0 - T1 + T2))/
     -                      deltA0T1T2 + (T2*tmp96)/(deltA0T1T2*T1)
     -                       ))/deltA0T1T2))/T2)))
        DT2T2 = hb**2*(0.25D0*
     -       (B1*(1 - c2b)*(1 - c2t)*(-1 + LogB1))/T2 + 
     -      0.25D0*(B2*(1 + c2b)*(1 - c2t)*(-1 + LogB2))/T2 + 
     -      0.5D0*(A0*(1 - c2t)*(-1 + LogA0)*sbe**2)/T2) + 
     -   ht**2*(-(((-1 + Logmu2)*mu2)/T2) - ((-1 + Logt)*t)/T2 + 
     -      2*(0.5D0 - 2*LogT2 + (phimu2tT2*(-mu2 - t + T2))/T2 - 
     -         0.5D0*(Logmu2*Logt) + 0.5D0*(Logmu2*LogT2) + 
     -         0.5D0*(Logt*LogT2) - 0.5D0*(Logt*(mu2 - t - T2))/T2 - 
     -         0.5D0*(Logmu2*(-mu2 + t - T2))/T2 + 
     -         0.5D0*(deltmu2tT2*
     -             ((mu2*phimu2tT2*(mu2 + t - T2))/
     -                (deltmu2tT2*T2) + 
     -               (mu2*tmp103)/(deltmu2tT2*T2)))/mu2) - 
     -      (-mu2 - t + T2)*
     -       (2/T2 - Logmu2/T2 - Logt/T2 - 
     -         0.5D0*(Logt*(mu2 - t - T2))/T2**2 - 
     -         0.5D0*(Logmu2*(-mu2 + t - T2))/T2**2 - 
     -         0.5D0*((2*mu2*phimu2tT2)/T2 + 
     -             4*(-mu2 - t + T2)*
     -              ((mu2*phimu2tT2*(mu2 + t - T2))/
     -                 (deltmu2tT2*T2) + 
     -                (mu2*tmp103)/(deltmu2tT2*T2)) + 
     -             deltmu2tT2*
     -              (((2 - Logmu2 - Logt + 2*LogT2)*mu2)/
     -                 (deltmu2tT2*T2) - 
     -                (mu2*phimu2tT2)/(deltmu2tT2*T2) - 
     -                (2*mu2*phimu2tT2*(mu2 + t - T2)*
     -                   (-mu2 - t + T2))/(deltmu2tT2**2*T2) - 
     -                (mu2*tmp103)/(deltmu2tT2*T2**2) - 
     -                (2*mu2*(-mu2 - t + T2)*tmp103)/
     -                 (deltmu2tT2**2*T2) + 
     -                ((mu2 + t - T2)*
     -                   ((mu2*phimu2tT2*(mu2 + t - T2))/
     -                      (deltmu2tT2*T2) + 
     -                     (mu2*tmp103)/(deltmu2tT2*T2)))/
     -                 deltmu2tT2))/mu2)) - 
     -   0.5D0*(((4*B2*(LogB2 - LogT2))/((1 - B2/T2)*T2**2) + 8/T2 - 
     -        (4*LogT2)/T2 - 
     -        2*((B2**2*(LogB2 - LogT2))/((1 - B2/T2)**2*T2**4) + 
     -           B2/((1 - B2/T2)*T2**3) + 
     -           (2*B2*(LogB2 - LogT2))/((1 - B2/T2)*T2**3))*
     -         (-B2 + T2) - (-2*B2*LogB2 + 4*T2)/T2**2 + 
     -        (B2 - T2)*tmp18)*tmp189) - 
     -   0.5D0*(((4*B1*(LogB1 - LogT2))/((1 - B1/T2)*T2**2) + 8/T2 - 
     -        (4*LogT2)/T2 - 
     -        2*((B1**2*(LogB1 - LogT2))/((1 - B1/T2)**2*T2**4) + 
     -           B1/((1 - B1/T2)*T2**3) + 
     -           (2*B1*(LogB1 - LogT2))/((1 - B1/T2)*T2**3))*
     -         (-B1 + T2) - (-2*B1*LogB1 + 4*T2)/T2**2 + 
     -        (B1 - T2)*tmp18)*tmp190)
        DT2T2 = DT2T2 - tmp198*
     -    (2/T2 - LogA0/T2 - LogB1/T2 - 
     -      0.5D0*(LogB1*(A0 - B1 - T2))/T2**2 - 
     -      0.5D0*(LogA0*(-A0 + B1 - T2))/T2**2 - 
     -      0.5D0*((2*B1*phiA0B1T2)/T2 + 
     -          4*(-A0 - B1 + T2)*
     -           ((B1*phiA0B1T2*(A0 + B1 - T2))/(deltA0B1T2*T2) + 
     -             (B1*tmp100)/(deltA0B1T2*T2)) + 
     -          deltA0B1T2*
     -           ((B1*(2 - LogA0 - LogB1 + 2*LogT2))/
     -              (deltA0B1T2*T2) - 
     -             (B1*phiA0B1T2)/(deltA0B1T2*T2) - 
     -             (2*B1*phiA0B1T2*(A0 + B1 - T2)*(-A0 - B1 + T2))/
     -              (deltA0B1T2**2*T2) - 
     -             (B1*tmp100)/(deltA0B1T2*T2**2) - 
     -             (2*B1*(-A0 - B1 + T2)*tmp100)/
     -              (deltA0B1T2**2*T2) + 
     -             ((A0 + B1 - T2)*
     -                ((B1*phiA0B1T2*(A0 + B1 - T2))/
     -                   (deltA0B1T2*T2) + 
     -                  (B1*tmp100)/(deltA0B1T2*T2)))/deltA0B1T2))/
     -        B1) - tmp197*
     -    (2/T2 - LogA0/T2 - LogB2/T2 - 
     -      0.5D0*(LogB2*(A0 - B2 - T2))/T2**2 - 
     -      0.5D0*(LogA0*(-A0 + B2 - T2))/T2**2 - 
     -      0.5D0*((2*B2*phiA0B2T2)/T2 + 
     -          4*(-A0 - B2 + T2)*
     -           ((B2*phiA0B2T2*(A0 + B2 - T2))/(deltA0B2T2*T2) + 
     -             (B2*tmp101)/(deltA0B2T2*T2)) + 
     -          deltA0B2T2*
     -           ((B2*(2 - LogA0 - LogB2 + 2*LogT2))/
     -              (deltA0B2T2*T2) - 
     -             (B2*phiA0B2T2)/(deltA0B2T2*T2) - 
     -             (2*B2*phiA0B2T2*(A0 + B2 - T2)*(-A0 - B2 + T2))/
     -              (deltA0B2T2**2*T2) - 
     -             (B2*tmp101)/(deltA0B2T2*T2**2) - 
     -             (2*B2*(-A0 - B2 + T2)*tmp101)/
     -              (deltA0B2T2**2*T2) + 
     -             ((A0 + B2 - T2)*
     -                ((B2*phiA0B2T2*(A0 + B2 - T2))/
     -                   (deltA0B2T2*T2) + 
     -                  (B2*tmp101)/(deltA0B2T2*T2)))/deltA0B2T2))/
     -        B2) + tmp63*(-((b*(-1 + Logb))/T2) - 
     -      ((-1 + Logmu2)*mu2)/T2 + 2*tmp225 - 
     -      (-b - mu2 + T2)*
     -       (2/T2 - Logb/T2 - Logmu2/T2 - 
     -         0.5D0*(Logmu2*(b - mu2 - T2))/T2**2 - 
     -         0.5D0*(Logb*(-b + mu2 - T2))/T2**2 - 
     -         0.5D0*(2*phiT2bmu2 + 
     -             4*(-b - mu2 + T2)*
     -              ((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -                (mu2*tmp102)/(deltT2bmu2*T2)) + 
     -             deltT2bmu2*tmp226)/mu2)) - 
     -   2*hb*ht*mb*mu*s2t*
     -    (-2/T2 + Logb/T2 + Logmu2/T2 + 
     -      0.5D0*(Logmu2*(b - mu2 - T2))/T2**2 + 
     -      0.5D0*(Logb*(-b + mu2 - T2))/T2**2 + 
     -      0.5D0*(2*phiT2bmu2 + 
     -          4*(-b - mu2 + T2)*
     -           ((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -             (mu2*tmp102)/(deltT2bmu2*T2)) + 
     -          deltT2bmu2*tmp226)/mu2)
        DT2T2 = DT2T2 + ht**2*
     -    (((-1 + LogT1)*T1*tmp8)/T2 + 
     -      (A0*cbe**2*(-1 + LogA0)*(1 + 0.5D0*(1 + c2t)))/T2 + 
     -      0.25D0*(B1*(1 + c2b)*(1 + c2t)*(-1 + LogB1))/T2 + 
     -      0.25D0*(B2*(1 - c2b)*(1 + c2t)*(-1 + LogB2))/T2 + 
     -      0.25D0*((1 + Nc)*s2t**2*
     -         (8*(-1 + LogT2) + 2*(-1 + LogT2)**2 + 
     -           (2/T2**2 - (2*(-1 + LogT2))/T2**2)*T2**2)))
        DT2T2 = DT2T2 - 0.25D0*
     -    (ht**2*((1 + c2t**2)*sbe**2*
     -         (4/T2 - (2*LogT1)/T2 + 
     -           (4*(-LogT1 + LogT2))/(T1*(1 - T2/T1)) - 
     -           2*(T1 - T2)*
     -            ((-LogT1 + LogT2)/(T1**2*(1 - T2/T1)**2) + 
     -              1/(T1*T2*(1 - T2/T1))))*Xt**2 + 
     -        2*(1 + c2t**2)*cbe**2*Yt**2*
     -         (-((deltA0T1T2*phiA0T1T2)/T2**3) + 2/T2 - 
     -           LogA0/T2 - LogT1/T2 + 
     -           (2*phiA0T1T2*(-A0 - T1 + T2) + 
     -              deltA0T1T2*
     -               (tmp104/deltA0T1T2 + 
     -                 (phiA0T1T2*tmp179)/deltA0T1T2))/T2**2 - 
     -           0.5D0*(LogT1*(A0 - T1 - T2))/T2**2 - 
     -           0.5D0*(LogA0*(-A0 + T1 - T2))/T2**2 - 
     -           0.5D0*(2*phiA0T1T2 + 
     -               4*(-A0 - T1 + T2)*
     -                (tmp104/deltA0T1T2 + 
     -                  (phiA0T1T2*tmp179)/deltA0T1T2) + 
     -               deltA0T1T2*
     -                ((2 - LogA0 - LogT1 + 2*LogT2)/deltA0T1T2 - 
     -                  (phiA0T1T2*(A0 - T1)**2)/
     -                   (deltA0T1T2*T2**2) - 
     -                  (2*(-A0 - T1 + T2)*tmp104)/deltA0T1T2**2 - 
     -                  (2*phiA0T1T2*(-A0 - T1 + T2)*tmp179)/
     -                   deltA0T1T2**2 + 
     -                  (tmp179*
     -                     (tmp104/deltA0T1T2 + 
     -                       (phiA0T1T2*tmp179)/deltA0T1T2))/
     -                   deltA0T1T2))/T2) + 
     -        0.5D0*(sbe**2*(4/T2 + (2*(2 - 2*LogT2))/T2 - 
     -             (2*LogT2)/T2 - (4*T2 - 2*LogT2*T2)/T2**2)*tmp34)
     -          + cbe**2*tmp50*
     -         (-((deltA0T2T2*phiA0T2T2)/T2**3) + 
     -           (A0*LogA0)/T2**2 + 4/T2 - (4*LogT2)/T2 + 
     -           (-4*A0*phiA0T2T2 + 
     -              deltA0T2T2*
     -               ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                 tmp107/deltA0T2T2 + 
     -                 (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                 tmp99/deltA0T2T2))/T2**2 + 
     -           0.5D0*((A0 - 2*T2)*tmp18) - 
     -           0.5D0*(-8*A0*
     -                ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                  tmp107/deltA0T2T2 + 
     -                  (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                  tmp99/deltA0T2T2) + 
     -               deltA0T2T2*
     -                ((4*A0**2*phiA0T2T2)/deltA0T2T2**2 + 
     -                  (phiA0T2T2*
     -                     (-1 - (A0 - T2)**2/T2**2 - 
     -                       (2*(A0 - T2))/T2))/deltA0T2T2 + 
     -                  (1 + (A0 - T2)/T2)/deltA0T2T2 + 
     -                  (1 - (-A0 + T2)/T2)/deltA0T2T2 + 
     -                  (4*A0*tmp107)/deltA0T2T2**2 + 
     -                  (4*A0*phiA0T2T2*tmp182)/deltA0T2T2**2 + 
     -                  (4*A0*tmp99)/deltA0T2T2**2 + 
     -                  (A0*
     -                     ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                       tmp107/deltA0T2T2 + 
     -                       (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                       tmp99/deltA0T2T2))/deltA0T2T2 + 
     -                  (tmp182*
     -                     ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                       tmp107/deltA0T2T2 + 
     -                       (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                       tmp99/deltA0T2T2))/deltA0T2T2))/T2)))
        tmp245 = (2*(-1 + Logmu2)*mu2)/t - ((-1 + LogT1)*T1)/t - 
     -   ((-1 + LogT2)*T2)/t + 2*tmp228 + 2*tmp243 + 
     -   sbe**2*(8*(-1 + Logt) + 2*(-1 + Logt)**2 + 
     -      t*(4/t + (2*(2 - 2*Logt))/t - (2*Logt)/t - 
     -         (4*t - 2*Logt*t)/t**2) + t**2*tmp15 + 2*tmp77) - 
     -   (-mu2 - t + T2)*(2/t - Logmu2/t - LogT2/t - 
     -      0.5D0*(LogT2*(mu2 - t - T2))/t**2 - 
     -      0.5D0*(Logmu2*(-mu2 - t + T2))/t**2 - 
     -      0.5D0*((2*mu2*phimu2tT2)/T2 + 
     -          4*(-mu2 + t - T2)*
     -           ((mu2*phimu2tT2*(mu2 - t + T2))/(deltmu2tT2*T2) + 
     -             (mu2*tmp111)/(deltmu2tT2*t)) + 
     -          deltmu2tT2*
     -           (((2 - Logmu2 + 2*Logt - LogT2)*mu2)/
     -              (deltmu2tT2*t) - 
     -             (mu2*phimu2tT2)/(deltmu2tT2*T2) - 
     -             (2*mu2*phimu2tT2*(-mu2 + t - T2)*
     -                (mu2 - t + T2))/(deltmu2tT2**2*T2) - 
     -             (mu2*tmp111)/(deltmu2tT2*t**2) - 
     -             (2*mu2*(-mu2 + t - T2)*tmp111)/
     -              (deltmu2tT2**2*t) + 
     -             ((mu2 - t + T2)*
     -                ((mu2*phimu2tT2*(mu2 - t + T2))/
     -                   (deltmu2tT2*T2) + 
     -                  (mu2*tmp111)/(deltmu2tT2*t)))/deltmu2tT2))/
     -        mu2) + cbe**2*
     -    (8*(-1 + Logt) + 2*(-1 + Logt)**2 - 
     -      (2*A0*(-1 + LogA0))/t + t**2*tmp15 + 
     -      4*(-1 + 4*Logt - Logt**2 - (A0*LogA0)/t + 
     -         (2*A0*phiA0tt)/t + (Logt*(A0 - 2*t))/t + 
     -         0.5D0*(deltA0tt*phiA0tt)/t**2 - 
     -         0.5D0*(deltA0tt*
     -             ((A0*phiA0tt)/deltA0tt + 
     -               (phiA0tt*tmp175)/deltA0tt + tmp65/deltA0tt + 
     -               tmp70/deltA0tt))/t) - 
     -      (A0 - 2*t)*(-((deltA0tt*phiA0tt)/t**3) + 
     -         (A0*LogA0)/t**2 + 4/t - (4*Logt)/t + 
     -         (-4*A0*phiA0tt + 
     -            deltA0tt*
     -             ((A0*phiA0tt)/deltA0tt + 
     -               (phiA0tt*tmp175)/deltA0tt + tmp65/deltA0tt + 
     -               tmp70/deltA0tt))/t**2 + 
     -         0.5D0*((A0 - 2*t)*tmp16) - 
     -         0.5D0*(-8*A0*((A0*phiA0tt)/deltA0tt + 
     -                (phiA0tt*tmp175)/deltA0tt + tmp65/deltA0tt + 
     -                tmp70/deltA0tt) + 
     -             deltA0tt*
     -              ((4*A0**2*phiA0tt)/deltA0tt**2 + 
     -                (phiA0tt*
     -                   (-1 - (A0 - t)**2/t**2 - (2*(A0 - t))/t))/
     -                 deltA0tt + (1 + (A0 - t)/t)/deltA0tt + 
     -                (1 - (-A0 + t)/t)/deltA0tt + 
     -                (4*A0*phiA0tt*tmp175)/deltA0tt**2 + 
     -                (4*A0*tmp65)/deltA0tt**2 + 
     -                (4*A0*tmp70)/deltA0tt**2 + 
     -                (A0*((A0*phiA0tt)/deltA0tt + 
     -                     (phiA0tt*tmp175)/deltA0tt + 
     -                     tmp65/deltA0tt + tmp70/deltA0tt))/
     -                 deltA0tt + 
     -                (tmp175*
     -                   ((A0*phiA0tt)/deltA0tt + 
     -                     (phiA0tt*tmp175)/deltA0tt + 
     -                     tmp65/deltA0tt + tmp70/deltA0tt))/
     -                 deltA0tt))/t))
        tmp245 = tmp245 - 
     -   (-mu2 - t + T1)*(2/t - Logmu2/t - LogT1/t - 
     -      0.5D0*(LogT1*(mu2 - t - T1))/t**2 - 
     -      0.5D0*(Logmu2*(-mu2 - t + T1))/t**2 - 
     -      0.5D0*((2*mu2*phimu2tT1)/T1 + 
     -          4*(-mu2 + t - T1)*
     -           ((mu2*phimu2tT1*(mu2 - t + T1))/(deltmu2tT1*T1) + 
     -             (mu2*tmp91)/(deltmu2tT1*t)) + 
     -          deltmu2tT1*
     -           (((2 - Logmu2 + 2*Logt - LogT1)*mu2)/
     -              (deltmu2tT1*t) - 
     -             (mu2*phimu2tT1)/(deltmu2tT1*T1) - 
     -             (2*mu2*phimu2tT1*(-mu2 + t - T1)*
     -                (mu2 - t + T1))/(deltmu2tT1**2*T1) - 
     -             (mu2*tmp91)/(deltmu2tT1*t**2) - 
     -             (2*mu2*(-mu2 + t - T1)*tmp91)/
     -              (deltmu2tT1**2*t) + 
     -             ((mu2 - t + T1)*
     -                ((mu2*phimu2tT1*(mu2 - t + T1))/
     -                   (deltmu2tT1*T1) + 
     -                  (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))/
     -        mu2)
        Dtt = ht**2*tmp245 + 
     -   tmp10*(-5 - 2*Li2bt + 4*Logt - Logt**2 + 
     -      (b*(-1 + Logb))/t + (2*Logt*(b - t))/t + 
     -      (2*b*(Logb - Logt)*(-b + t))/((1 - b/t)*t**2) + 
     -      (-2*b*Logb + 4*t)/t + 
     -      0.5D0*((b + t)*((4*b*(Logb - Logt))/((1 - b/t)*t**2) + 
     -           8/t - (4*Logt)/t - (-2*b*Logb + 4*t)/t**2 + 
     -           (b - t)*tmp16 - 2*(-b + t)*tmp207))) + 
     -   tmp61*(-((B1*(-1 + LogB1))/t) + ((-1 + Logmu2)*mu2)/t + 
     -      2*tmp231 - (B1 - mu2 - t)*
     -       (2/t - LogB1/t - Logmu2/t - 
     -         0.5D0*(Logmu2*(B1 - mu2 - t))/t**2 - 
     -         0.5D0*(LogB1*(-B1 + mu2 - t))/t**2 - 
     -         0.5D0*(2*phiB1tmu2 + deltB1tmu2*tmp232 + 
     -             4*(-B1 - mu2 + t)*
     -              ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -                (mu2*tmp67)/(deltB1tmu2*t)))/mu2)) + 
     -   tmp60*(-((B2*(-1 + LogB2))/t) + ((-1 + Logmu2)*mu2)/t + 
     -      2*tmp234 - (B2 - mu2 - t)*
     -       (2/t - LogB2/t - Logmu2/t - 
     -         0.5D0*(Logmu2*(B2 - mu2 - t))/t**2 - 
     -         0.5D0*(LogB2*(-B2 + mu2 - t))/t**2 - 
     -         0.5D0*(2*phiB2tmu2 + deltB2tmu2*tmp236 + 
     -             4*(-B2 - mu2 + t)*
     -              ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -                (mu2*tmp68)/(deltB2tmu2*t)))/mu2)) - 
     -   2*hb*ht*mu*s2b*(tmp235/mt + 
     -      mt*(-(LogB1/t) + LogB2/t - 
     -         0.5D0*(Logmu2*(B1 - mu2 - t))/t**2 + 
     -         0.5D0*(Logmu2*(B2 - mu2 - t))/t**2 - 
     -         0.5D0*(LogB1*(-B1 + mu2 - t))/t**2 + 
     -         0.5D0*(LogB2*(-B2 + mu2 - t))/t**2 - 
     -         0.5D0*(2*phiB1tmu2 + deltB1tmu2*tmp232 + 
     -             4*(-B1 - mu2 + t)*
     -              ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -                (mu2*tmp67)/(deltB1tmu2*t)))/mu2 + 
     -         0.5D0*(2*phiB2tmu2 + deltB2tmu2*tmp236 + 
     -             4*(-B2 - mu2 + t)*
     -              ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -                (mu2*tmp68)/(deltB2tmu2*t)))/mu2) - 
     -      0.25D0*tmp75/t**1.5D0)
        Dtt = Dtt + tmp9*
     -    (-((A0*(-1 + LogA0))/t) + (b*(-1 + Logb))/t + 2*tmp229 - 
     -      (A0 - b - t)*(2/t - LogA0/t - Logb/t - 
     -         0.5D0*(Logb*(A0 - b - t))/t**2 - 
     -         0.5D0*(LogA0*(-A0 + b - t))/t**2 - 
     -         0.5D0*((2*b*phiA0bt)/t + deltA0bt*tmp230 + 
     -             4*(-A0 - b + t)*
     -              ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -                (b*tmp66)/(deltA0bt*t)))/b)) - 
     -   4*cbe*hb*ht*mb*sbe*
     -    ((0.5D0 - 2*Logt + (phiA0bt*(-A0 - b + t))/t - 
     -         0.5D0*(LogA0*Logb) + 0.5D0*(LogA0*Logt) + 
     -         0.5D0*(Logb*Logt) - 0.5D0*(Logb*(A0 - b - t))/t - 
     -         0.5D0*(LogA0*(-A0 + b - t))/t + 0.5D0*tmp210 + 
     -         0.5D0*(deltA0bt*
     -             ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -               (b*tmp66)/(deltA0bt*t)))/b)/mt + 
     -      mt*(-2/t + LogA0/t + Logb/t + 
     -         0.5D0*(Logb*(A0 - b - t))/t**2 + 
     -         0.5D0*(LogA0*(-A0 + b - t))/t**2 + 
     -         0.5D0*((4*b*(Logb - Logt))/((1 - b/t)*t**2) + 8/t - 
     -            (4*Logt)/t - (-2*b*Logb + 4*t)/t**2 + 
     -            (b - t)*tmp16 - 2*(-b + t)*tmp207) + 
     -         0.5D0*((2*b*phiA0bt)/t + deltA0bt*tmp230 + 
     -             4*(-A0 - b + t)*
     -              ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -                (b*tmp66)/(deltA0bt*t)))/b) - 
     -      0.25D0*(-2*A0*LogA0 - 2*b*Logb - 2*Logt*t - 
     -          0.5D0*(Logb*Logt*(A0 - b - t)) - 
     -          0.5D0*(LogA0*Logt*(-A0 + b - t)) + 
     -          0.5D0*(deltA0bt*phiA0bt)/t - 
     -          0.5D0*(LogA0*Logb*(-A0 - b + t)) + 
     -          2.5D0*(A0 + b + t) + 0.5D0*tmp76)/t**1.5D0) + 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (cbe*hb*ht*sbe*
     -         (0.25D0*(s2b*tmp130)/t**1.5D0
     $       - 0.5D0*(mb*tmp56)/t**1.5D0)-
     -          cbe**2*hb**2*
     -         (0.125D0*(mb*s2b*s2t)/t**1.5D0 - 
     -           0.125D0*((1 + c2b)*s2t*Xb)/t**1.5D0) - 
     -        ht**2*sbe**2*
     -         (0.125D0*(mb*s2b*s2t)/t**1.5D0 - 
     -           0.125D0*((1 - c2b)*s2t*Xt)/t**1.5D0)))
        Dtt = Dtt + tmp89*
     -    (-(cbe*hb*ht*sbe*
     -         (0.25D0*(s2b*tmp156)/t**1.5D0
     $       - 0.5D0*(mb*tmp56)/t**1.5D0))-
     -        hb**2*sbe**2*
     -       (0.125D0*(mb*s2b*s2t)/t**1.5D0 - 
     -         0.125D0*((1 + c2b)*s2t*Yb)/t**1.5D0) - 
     -      cbe**2*ht**2*(0.125D0*(mb*s2b*s2t)/t**1.5D0 - 
     -         0.125D0*((1 - c2b)*s2t*Yt)/t**1.5D0)) + 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (cbe*hb*ht*sbe*
     -         (0.25D0*(s2b*tmp131)/t**1.5D0
     $       - 0.5D0*(mb*tmp59)/t**1.5D0)-
     -          cbe**2*hb**2*
     -         (-(0.125D0*(mb*s2b*s2t)/t**1.5D0) + 
     -           0.125D0*((1 + c2b)*s2t*Xb)/t**1.5D0) - 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(mb*s2b*s2t)/t**1.5D0) + 
     -           0.125D0*((1 - c2b)*s2t*Xt)/t**1.5D0))) + 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (cbe*hb*ht*sbe*
     -         (-(0.25D0*(s2b*tmp130)/t**1.5D0) - 
     -           0.5D0*(mb*tmp59)/t**1.5D0) - 
     -        cbe**2*hb**2*
     -         (-(0.125D0*(mb*s2b*s2t)/t**1.5D0) - 
     -           0.125D0*((1 - c2b)*s2t*Xb)/t**1.5D0) - 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(mb*s2b*s2t)/t**1.5D0) - 
     -           0.125D0*((1 + c2b)*s2t*Xt)/t**1.5D0))) + 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (cbe*hb*ht*sbe*
     -         (-(0.25D0*(s2b*tmp131)/t**1.5D0) - 
     -           0.5D0*(mb*tmp56)/t**1.5D0) - 
     -        cbe**2*hb**2*
     -         (0.125D0*(mb*s2b*s2t)/t**1.5D0 + 
     -           0.125D0*((1 - c2b)*s2t*Xb)/t**1.5D0) - 
     -        ht**2*sbe**2*
     -         (0.125D0*(mb*s2b*s2t)/t**1.5D0 + 
     -           0.125D0*((1 + c2b)*s2t*Xt)/t**1.5D0))) - 
     -   0.25D0*(ht**2*((cbe**2*s2t*tmp114*Yt)/t**1.5D0 - 
     -        (cbe**2*s2t*tmp93*Yt)/t**1.5D0 + 
     -        0.5D0*(s2t*sbe**2*tmp115*Xt)/t**1.5D0 - 
     -        0.5D0*(s2t*sbe**2*tmp94*Xt)/t**1.5D0))
        Dtt = Dtt + tmp109*
     -    (-(cbe*hb*ht*sbe*
     -         (0.25D0*(s2b*tmp157)/t**1.5D0
     $       - 0.5D0*(mb*tmp59)/t**1.5D0))-
     -        hb**2*sbe**2*
     -       (-(0.125D0*(mb*s2b*s2t)/t**1.5D0) + 
     -         0.125D0*((1 + c2b)*s2t*Yb)/t**1.5D0) - 
     -      cbe**2*ht**2*(-(0.125D0*(mb*s2b*s2t)/t**1.5D0) + 
     -         0.125D0*((1 - c2b)*s2t*Yt)/t**1.5D0)) + 
     -   tmp88*(-(cbe*hb*ht*sbe*
     -         (-(0.25D0*(s2b*tmp156)/t**1.5D0) - 
     -           0.5D0*(mb*tmp59)/t**1.5D0)) - 
     -      hb**2*sbe**2*(-(0.125D0*(mb*s2b*s2t)/t**1.5D0) - 
     -         0.125D0*((1 - c2b)*s2t*Yb)/t**1.5D0) - 
     -      cbe**2*ht**2*(-(0.125D0*(mb*s2b*s2t)/t**1.5D0) - 
     -         0.125D0*((1 + c2b)*s2t*Yt)/t**1.5D0)) + 
     -   tmp108*(-(cbe*hb*ht*sbe*
     -         (-(0.25D0*(s2b*tmp157)/t**1.5D0) - 
     -           0.5D0*(mb*tmp56)/t**1.5D0)) - 
     -      hb**2*sbe**2*(0.125D0*(mb*s2b*s2t)/t**1.5D0 + 
     -         0.125D0*((1 - c2b)*s2t*Yb)/t**1.5D0) - 
     -      cbe**2*ht**2*(0.125D0*(mb*s2b*s2t)/t**1.5D0 + 
     -         0.125D0*((1 + c2b)*s2t*Yt)/t**1.5D0))
        Dc2tc2t = (b*(-1 + Logb)*(-1 + Logmu2)*mu2 - 
     -      b*(-1 + Logb)*(-1 + LogT2)*T2 - 
     -      (-1 + Logmu2)*(-1 + LogT2)*mu2*T2 - 
     -      (-b - mu2 + T2)*tmp112)*
     -    (0.125D0*hb**2/c2t**3 - 0.125D0*ht**2/c2t**3) + 
     -   (b*(-1 + Logb)*(-1 + Logmu2)*mu2 - 
     -      b*(-1 + Logb)*(-1 + LogT1)*T1 - 
     -      (-1 + Logmu2)*(-1 + LogT1)*mu2*T1 - 
     -      (-b - mu2 + T1)*tmp92)*
     -    (-(0.125D0*hb**2/c2t**3) + 0.125D0*ht**2/c2t**3) + 
     -   ht**2*(0.0625D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT1)*T1)/
     -        c2t**3 + 0.0625D0*
     -       (B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT1)*T1)/c2t**3 + 
     -      cbe**2*(0.125D0*(A0*(-1 + LogA0)*(-1 + LogT1)*T1)/
     -           c2t**3 - 0.125D0*
     -          (A0*(-1 + LogA0)*(-1 + LogT2)*T2)/c2t**3) - 
     -      0.0625D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/
     -        c2t**3 - 0.0625D0*
     -       (B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT2)*T2)/c2t**3) + 
     -   hb**2*(-(0.0625D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT1)*
     -            T1)/c2t**3) - 
     -      0.0625D0*(B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT1)*T1)/
     -        c2t**3 + sbe**2*
     -       (-(0.125D0*(A0*(-1 + LogA0)*(-1 + LogT1)*T1)/c2t**3) + 
     -         0.125D0*(A0*(-1 + LogA0)*(-1 + LogT2)*T2)/c2t**3) + 
     -      0.0625D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/
     -        c2t**3 + 0.0625D0*
     -       (B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT2)*T2)/c2t**3) + 
     -   0.5D0*(hb*ht*mb*mu*tmp113)/s2t**3 + 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(cbe**2*hb**2*
     -           (0.0625D0*(b*(1 - c2b))/c2t**3 - 
     -             0.125D0*(mb*mt*s2b)/s2t**3 - 
     -             0.0625D0*((1 + c2b)*t)/c2t**3 - 
     -             0.125D0*(mb*s2b*Xb)/c2t**3 + 
     -             0.125D0*((1 + c2b)*mt*Xb)/s2t**3 + 
     -             0.0625D0*((1 + c2b)*Xb**2)/c2t**3)) + 
     -        cbe*hb*ht*sbe*
     -         (-(mt*s2b*tmp25) + 2*mb*mt*tmp52 - 
     -           0.125D0*(s2b*(b + t))/s2t**3 + 
     -           0.25D0*(mb*tmp128)/s2t**3 - 0.125D0*(s2b*Xb*Xt)/s2t**3
     -           ) - ht**2*sbe**2*
     -         (-(0.0625D0*(b*(1 + c2b))/c2t**3) - 
     -           0.125D0*(mb*mt*s2b)/s2t**3 + 
     -           0.0625D0*((1 - c2b)*t)/c2t**3 + 
     -           0.125D0*(mb*s2b*Xt)/c2t**3 + 
     -           0.125D0*((1 - c2b)*mt*Xt)/s2t**3 - 
     -           0.0625D0*((1 - c2b)*Xt**2)/c2t**3)))
        Dc2tc2t = Dc2tc2t + 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (-(cbe**2*hb**2*
     -           (-(0.0625D0*(b*(1 - c2b))/c2t**3) + 
     -             0.125D0*(mb*mt*s2b)/s2t**3 + 
     -             0.0625D0*((1 + c2b)*t)/c2t**3 + 
     -             0.125D0*(mb*s2b*Xb)/c2t**3 - 
     -             0.125D0*((1 + c2b)*mt*Xb)/s2t**3 - 
     -             0.0625D0*((1 + c2b)*Xb**2)/c2t**3)) + 
     -        cbe*hb*ht*sbe*
     -         (-(mt*s2b*tmp26) + 2*mb*mt*tmp53 + 
     -           0.125D0*(s2b*(b + t))/s2t**3 - 
     -           0.25D0*(mb*tmp128)/s2t**3 + 0.125D0*(s2b*Xb*Xt)/s2t**3
     -           ) - ht**2*sbe**2*
     -         (0.0625D0*(b*(1 + c2b))/c2t**3 + 
     -           0.125D0*(mb*mt*s2b)/s2t**3 - 
     -           0.0625D0*((1 - c2b)*t)/c2t**3 - 
     -           0.125D0*(mb*s2b*Xt)/c2t**3 - 
     -           0.125D0*((1 - c2b)*mt*Xt)/s2t**3 + 
     -           0.0625D0*((1 - c2b)*Xt**2)/c2t**3))) + 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (-(cbe**2*hb**2*
     -           (0.0625D0*(b*(1 + c2b))/c2t**3 + 
     -             0.125D0*(mb*mt*s2b)/s2t**3 - 
     -             0.0625D0*((1 - c2b)*t)/c2t**3 + 
     -             0.125D0*(mb*s2b*Xb)/c2t**3 + 
     -             0.125D0*((1 - c2b)*mt*Xb)/s2t**3 + 
     -             0.0625D0*((1 - c2b)*Xb**2)/c2t**3)) + 
     -        cbe*hb*ht*sbe*
     -         (mt*s2b*tmp25 + 2*mb*mt*tmp53 + 
     -           0.125D0*(s2b*(b + t))/s2t**3 + 
     -           0.25D0*(mb*tmp129)/s2t**3 + 0.125D0*(s2b*Xb*Xt)/s2t**3
     -           ) - ht**2*sbe**2*
     -         (-(0.0625D0*(b*(1 - c2b))/c2t**3) + 
     -           0.125D0*(mb*mt*s2b)/s2t**3 + 
     -           0.0625D0*((1 + c2b)*t)/c2t**3 - 
     -           0.125D0*(mb*s2b*Xt)/c2t**3 + 
     -           0.125D0*((1 + c2b)*mt*Xt)/s2t**3 - 
     -           0.0625D0*((1 + c2b)*Xt**2)/c2t**3)))
        Dc2tc2t = Dc2tc2t + 
     -   tmp109*(-(hb**2*sbe**2*
     -         (0.0625D0*(b*(1 - c2b))/c2t**3 - 
     -           0.125D0*(mb*mt*s2b)/s2t**3 - 
     -           0.0625D0*((1 + c2b)*t)/c2t**3 - 
     -           0.125D0*(mb*s2b*Yb)/c2t**3 + 
     -           0.125D0*((1 + c2b)*mt*Yb)/s2t**3 + 
     -           0.0625D0*((1 + c2b)*Yb**2)/c2t**3)) - 
     -      cbe*hb*ht*sbe*(-(mt*s2b*tmp41) + 2*mb*mt*tmp52 - 
     -         0.125D0*(s2b*(b + t))/s2t**3 + 
     -         0.25D0*(mb*tmp154)/s2t**3 - 0.125D0*(s2b*Yb*Yt)/s2t**3)-
     -        cbe**2*ht**2*
     -       (-(0.0625D0*(b*(1 + c2b))/c2t**3) - 
     -         0.125D0*(mb*mt*s2b)/s2t**3 + 
     -         0.0625D0*((1 - c2b)*t)/c2t**3 + 
     -         0.125D0*(mb*s2b*Yt)/c2t**3 + 
     -         0.125D0*((1 - c2b)*mt*Yt)/s2t**3 - 
     -         0.0625D0*((1 - c2b)*Yt**2)/c2t**3)) - 
     -   0.25D0*(ht**2*((cbe**2*mt*tmp114*Yt)/s2t**3 - 
     -        (cbe**2*mt*tmp93*Yt)/s2t**3 + 
     -        0.5D0*(mt*sbe**2*tmp115*Xt)/s2t**3 - 
     -        0.5D0*(mt*sbe**2*tmp94*Xt)/s2t**3)) + 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (-(cbe**2*hb**2*
     -           (-(0.0625D0*(b*(1 + c2b))/c2t**3) - 
     -             0.125D0*(mb*mt*s2b)/s2t**3 + 
     -             0.0625D0*((1 - c2b)*t)/c2t**3 - 
     -             0.125D0*(mb*s2b*Xb)/c2t**3 - 
     -             0.125D0*((1 - c2b)*mt*Xb)/s2t**3 - 
     -             0.0625D0*((1 - c2b)*Xb**2)/c2t**3)) + 
     -        cbe*hb*ht*sbe*
     -         (mt*s2b*tmp26 + 2*mb*mt*tmp52 - 
     -           0.125D0*(s2b*(b + t))/s2t**3 - 
     -           0.25D0*(mb*tmp129)/s2t**3 - 0.125D0*(s2b*Xb*Xt)/s2t**3
     -           ) - ht**2*sbe**2*
     -         (0.0625D0*(b*(1 - c2b))/c2t**3 - 
     -           0.125D0*(mb*mt*s2b)/s2t**3 - 
     -           0.0625D0*((1 + c2b)*t)/c2t**3 + 
     -           0.125D0*(mb*s2b*Xt)/c2t**3 - 
     -           0.125D0*((1 + c2b)*mt*Xt)/s2t**3 + 
     -           0.0625D0*((1 + c2b)*Xt**2)/c2t**3)))
        Dc2tc2t = Dc2tc2t + 
     -   tmp89*(-(hb**2*sbe**2*
     -         (-(0.0625D0*(b*(1 - c2b))/c2t**3) + 
     -           0.125D0*(mb*mt*s2b)/s2t**3 + 
     -           0.0625D0*((1 + c2b)*t)/c2t**3 + 
     -           0.125D0*(mb*s2b*Yb)/c2t**3 - 
     -           0.125D0*((1 + c2b)*mt*Yb)/s2t**3 - 
     -           0.0625D0*((1 + c2b)*Yb**2)/c2t**3)) - 
     -      cbe*hb*ht*sbe*(-(mt*s2b*tmp42) + 2*mb*mt*tmp53 + 
     -         0.125D0*(s2b*(b + t))/s2t**3 - 
     -         0.25D0*(mb*tmp154)/s2t**3 + 0.125D0*(s2b*Yb*Yt)/s2t**3)-
     -        cbe**2*ht**2*
     -       (0.0625D0*(b*(1 + c2b))/c2t**3 + 
     -         0.125D0*(mb*mt*s2b)/s2t**3 - 
     -         0.0625D0*((1 - c2b)*t)/c2t**3 - 
     -         0.125D0*(mb*s2b*Yt)/c2t**3 - 
     -         0.125D0*((1 - c2b)*mt*Yt)/s2t**3 + 
     -         0.0625D0*((1 - c2b)*Yt**2)/c2t**3)) + 
     -   tmp108*(-(hb**2*sbe**2*
     -         (0.0625D0*(b*(1 + c2b))/c2t**3 + 
     -           0.125D0*(mb*mt*s2b)/s2t**3 - 
     -           0.0625D0*((1 - c2b)*t)/c2t**3 + 
     -           0.125D0*(mb*s2b*Yb)/c2t**3 + 
     -           0.125D0*((1 - c2b)*mt*Yb)/s2t**3 + 
     -           0.0625D0*((1 - c2b)*Yb**2)/c2t**3)) - 
     -      cbe*hb*ht*sbe*(mt*s2b*tmp41 + 2*mb*mt*tmp53 + 
     -         0.125D0*(s2b*(b + t))/s2t**3 + 
     -         0.25D0*(mb*tmp155)/s2t**3 + 0.125D0*(s2b*Yb*Yt)/s2t**3)-
     -        cbe**2*ht**2*
     -       (-(0.0625D0*(b*(1 - c2b))/c2t**3) + 
     -         0.125D0*(mb*mt*s2b)/s2t**3 + 
     -         0.0625D0*((1 + c2b)*t)/c2t**3 - 
     -         0.125D0*(mb*s2b*Yt)/c2t**3 + 
     -         0.125D0*((1 + c2b)*mt*Yt)/s2t**3 - 
     -         0.0625D0*((1 + c2b)*Yt**2)/c2t**3)) + 
     -   tmp88*(-(hb**2*sbe**2*
     -         (-(0.0625D0*(b*(1 + c2b))/c2t**3) - 
     -           0.125D0*(mb*mt*s2b)/s2t**3 + 
     -           0.0625D0*((1 - c2b)*t)/c2t**3 - 
     -           0.125D0*(mb*s2b*Yb)/c2t**3 - 
     -           0.125D0*((1 - c2b)*mt*Yb)/s2t**3 - 
     -           0.0625D0*((1 - c2b)*Yb**2)/c2t**3)) - 
     -      cbe*hb*ht*sbe*(mt*s2b*tmp42 + 2*mb*mt*tmp52 - 
     -         0.125D0*(s2b*(b + t))/s2t**3 - 
     -         0.25D0*(mb*tmp155)/s2t**3 - 0.125D0*(s2b*Yb*Yt)/s2t**3)-
     -        cbe**2*ht**2*
     -       (0.0625D0*(b*(1 - c2b))/c2t**3 - 
     -         0.125D0*(mb*mt*s2b)/s2t**3 - 
     -         0.0625D0*((1 + c2b)*t)/c2t**3 + 
     -         0.125D0*(mb*s2b*Yt)/c2t**3 - 
     -         0.125D0*((1 + c2b)*mt*Yt)/s2t**3 + 
     -         0.0625D0*((1 + c2b)*Yt**2)/c2t**3))
        DT1t = -(tmp192*tmp238) - tmp191*tmp239 - 
     -   0.5D0*(tmp184*tmp213) - 0.5D0*(tmp183*tmp214) + 
     -   ht**2*(1 - 3*Logt + Logmu2*Logt - 
     -      (-1 + Logt)*(-1 + LogT1) + LogT1 - Logmu2*LogT1 + 
     -      (phimu2tT1*(-mu2 + t - T1))/T1 - 
     -      (phimu2tT1*(-mu2 - t + T1))/T1 - 
     -      0.5D0*(LogT1*(mu2 - t - T1))/t + 
     -      0.5D0*(Logt*(mu2 - t - T1))/T1 + 
     -      0.5D0*(Logmu2*(-mu2 + t - T1))/T1 - 
     -      0.5D0*(Logmu2*(-mu2 - t + T1))/t - 
     -      0.5D0*(deltmu2tT1*
     -          ((mu2*phimu2tT1*(mu2 + t - T1))/(deltmu2tT1*T1) + 
     -            (mu2*tmp84)/(deltmu2tT1*T1)))/mu2 + 
     -      0.5D0*(deltmu2tT1*
     -          ((mu2*phimu2tT1*(mu2 - t + T1))/(deltmu2tT1*T1) + 
     -            (mu2*tmp91)/(deltmu2tT1*t)))/mu2 - 
     -      (-mu2 - t + T1)*
     -       (phimu2tT1/T1 - 
     -         ((-mu2 + t - T1)*
     -            ((mu2*phimu2tT1*(mu2 + t - T1))/
     -               (deltmu2tT1*T1) + (mu2*tmp84)/(deltmu2tT1*T1))
     -            )/mu2 - ((-mu2 - t + T1)*
     -            ((mu2*phimu2tT1*(mu2 - t + T1))/
     -               (deltmu2tT1*T1) + (mu2*tmp91)/(deltmu2tT1*t)))
     -           /mu2 + 0.5D0*Logmu2/t - 0.5D0*LogT1/t + 
     -         0.5D0*Logmu2/T1 - 0.5D0*Logt/T1 + 
     -         0.5D0*(mu2 - t - T1)/(t*T1) - 
     -         0.5D0*(deltmu2tT1*
     -             ((mu2*phimu2tT1)/(deltmu2tT1*T1) - 
     -               (2*mu2*phimu2tT1*(-mu2 + t - T1)*
     -                  (mu2 + t - T1))/(deltmu2tT1**2*T1) + 
     -               (mu2*(Logmu2 - Logt - (-mu2 + t)/t - T1/t))/
     -                (deltmu2tT1*T1) - 
     -               (2*mu2*(-mu2 + t - T1)*tmp84)/
     -                (deltmu2tT1**2*T1) + 
     -               ((mu2 + t - T1)*
     -                  ((mu2*phimu2tT1*(mu2 - t + T1))/
     -                     (deltmu2tT1*T1) + 
     -                    (mu2*tmp91)/(deltmu2tT1*t)))/deltmu2tT1))
     -            /mu2)) - 
     -   0.25D0*(ht**2*(cbe**2*(4 + (2*s2t*Yt)/mt)*
     -         (-1 + 4*LogT1 - LogT1**2 - (A0*LogA0)/T1 + 
     -           (2*A0*phiA0T1T1)/T1 + (LogT1*(A0 - 2*T1))/T1 + 
     -           0.5D0*(deltA0T1T1*phiA0T1T1)/T1**2 - 
     -           0.5D0*(deltA0T1T1*
     -               ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                 (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                 tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/T1) + 
     -        0.5D0*(sbe**2*tmp95*(4 + (2*s2t*Xt)/mt))))
        DT2t = -(tmp196*tmp222) - tmp195*tmp223 + 
     -   ht**2*(1 - 3*Logt + Logmu2*Logt - 
     -      (-1 + Logt)*(-1 + LogT2) + LogT2 - Logmu2*LogT2 + 
     -      (phimu2tT2*(-mu2 + t - T2))/T2 - 
     -      (phimu2tT2*(-mu2 - t + T2))/T2 - 
     -      0.5D0*(LogT2*(mu2 - t - T2))/t + 
     -      0.5D0*(Logt*(mu2 - t - T2))/T2 + 
     -      0.5D0*(Logmu2*(-mu2 + t - T2))/T2 - 
     -      0.5D0*(Logmu2*(-mu2 - t + T2))/t - 
     -      0.5D0*(deltmu2tT2*
     -          ((mu2*phimu2tT2*(mu2 + t - T2))/(deltmu2tT2*T2) + 
     -            (mu2*tmp103)/(deltmu2tT2*T2)))/mu2 + 
     -      0.5D0*(deltmu2tT2*
     -          ((mu2*phimu2tT2*(mu2 - t + T2))/(deltmu2tT2*T2) + 
     -            (mu2*tmp111)/(deltmu2tT2*t)))/mu2 - 
     -      (-mu2 - t + T2)*
     -       (phimu2tT2/T2 - 
     -         ((-mu2 + t - T2)*
     -            ((mu2*phimu2tT2*(mu2 + t - T2))/
     -               (deltmu2tT2*T2) + (mu2*tmp103)/(deltmu2tT2*T2)
     -              ))/mu2 - 
     -         ((-mu2 - t + T2)*
     -            ((mu2*phimu2tT2*(mu2 - t + T2))/
     -               (deltmu2tT2*T2) + (mu2*tmp111)/(deltmu2tT2*t))
     -            )/mu2 + 0.5D0*Logmu2/t - 0.5D0*LogT2/t + 
     -         0.5D0*Logmu2/T2 - 0.5D0*Logt/T2 + 
     -         0.5D0*(mu2 - t - T2)/(t*T2) - 
     -         0.5D0*(deltmu2tT2*
     -             ((mu2*phimu2tT2)/(deltmu2tT2*T2) - 
     -               (2*mu2*phimu2tT2*(-mu2 + t - T2)*
     -                  (mu2 + t - T2))/(deltmu2tT2**2*T2) + 
     -               (mu2*(Logmu2 - Logt - (-mu2 + t)/t - T2/t))/
     -                (deltmu2tT2*T2) - 
     -               (2*mu2*(-mu2 + t - T2)*tmp103)/
     -                (deltmu2tT2**2*T2) + 
     -               ((mu2 + t - T2)*
     -                  ((mu2*phimu2tT2*(mu2 - t + T2))/
     -                     (deltmu2tT2*T2) + 
     -                    (mu2*tmp111)/(deltmu2tT2*t)))/deltmu2tT2)
     -             )/mu2)) - 0.5D0*(tmp188*tmp217) - 
     -   0.5D0*(tmp187*tmp218) - 
     -   0.25D0*(ht**2*(cbe**2*(4 - (2*s2t*Yt)/mt)*
     -         (-1 + 4*LogT2 - LogT2**2 - (A0*LogA0)/T2 + 
     -           (2*A0*phiA0T2T2)/T2 + (LogT2*(A0 - 2*T2))/T2 + 
     -           0.5D0*(deltA0T2T2*phiA0T2T2)/T2**2 - 
     -           0.5D0*(deltA0T2T2*
     -               ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                 tmp107/deltA0T2T2 + 
     -                 (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                 tmp99/deltA0T2T2))/T2) + 
     -        0.5D0*(sbe**2*tmp116*(4 - (2*s2t*Xt)/mt))))
        DT1T2 = ht**2*(c2t**2 + (-1 + LogT1)*tmp8 + 
     -      (-1 + LogT2)*tmp8 + (-1 + LogT1)*(-1 + LogT2)*tmp8 - 
     -      0.5D0*((-1 + Nc)*s2t**2)) - 
     -   0.25D0*(ht**2*((1 + c2t**2)*sbe**2*
     -         ((2*LogT1)/T1 + (-2 - 2*LogT2)/T1 + 
     -           (2*(-LogT1 + LogT2)*(T1 - T2)*T2)/
     -            (T1**3*(1 - T2/T1)**2) - 
     -           (2*(-LogT1 + LogT2))/(T1*(1 - T2/T1)) + 
     -           (2*(T1 - T2))/(T1**2*(1 - T2/T1)) + 
     -           (2*(-LogT1 + LogT2)*(T1 - T2))/
     -            (T1**2*(1 - T2/T1)) - 
     -           (2*(-LogT1 + LogT2)*T2)/(T1**2*(1 - T2/T1)))*Xt**2
     -          + 2*(1 + c2t**2)*cbe**2*Yt**2*
     -         ((phiA0T1T2*(-A0 + T1 - T2))/T2**2 + phiA0T1T2/T2 - 
     -           ((-A0 + T1 - T2)*
     -              (tmp104/deltA0T1T2 + 
     -                (phiA0T1T2*tmp179)/deltA0T1T2))/T2 - 
     -           ((-A0 - T1 + T2)*
     -              ((phiA0T1T2*(A0 - T1 + T2))/deltA0T1T2 + 
     -                (T2*tmp96)/(deltA0T1T2*T1)))/T2 + 
     -           0.5D0*LogA0/T1 - 0.5D0*LogT2/T1 + 0.5D0*LogA0/T2 - 
     -           0.5D0*LogT1/T2 + 0.5D0*(A0 - T1 - T2)/(T1*T2) + 
     -           0.5D0*(deltA0T1T2*
     -               ((phiA0T1T2*(A0 - T1 + T2))/deltA0T1T2 + 
     -                 (T2*tmp96)/(deltA0T1T2*T1)))/T2**2 - 
     -           0.5D0*(deltA0T1T2*
     -               (phiA0T1T2/deltA0T1T2 + 
     -                 ((LogA0 - LogT2 - T1/T2 + (A0 - T2)/T2)*T2)/
     -                  (deltA0T1T2*T1) - 
     -                 (2*phiA0T1T2*(-A0 - T1 + T2)*
     -                    (A0 - T1 + T2))/deltA0T1T2**2 + 
     -                 ((A0 - T1 + T2)*
     -                    (tmp104/deltA0T1T2 + 
     -                      (phiA0T1T2*tmp179)/deltA0T1T2))/
     -                  deltA0T1T2 + tmp96/(deltA0T1T2*T1) - 
     -                 (2*T2*(-A0 - T1 + T2)*tmp96)/
     -                  (deltA0T1T2**2*T1)))/T2)))
        Dtc2t = -(0.5D0*((-5*B2 + 4*B2*LogB2 + 
     -          LogT1**2*(B2 - T1) - 5*T1 + 
     -          LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -        (-(cbe*hb*ht*sbe*
     -             ((mb*tmp54)/mt + 0.25D0*s2b/s2t - 
     -               0.5D0*(s2b*tmp27)/mt)) + 
     -          cbe**2*hb**2*
     -           (-(0.125D0*(1 + c2b)/c2t) + 
     -             0.125D0*(mb*s2b)/(mt*s2t) - 
     -             0.125D0*((1 + c2b)*Xb)/(mt*s2t)) + 
     -          ht**2*sbe**2*
     -           (0.125D0*(1 - c2b)/c2t + 0.125D0*(mb*s2b)/(mt*s2t) - 
     -             0.125D0*((1 - c2b)*Xt)/(mt*s2t))))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mb*tmp55)/mt - 0.25D0*s2b/s2t - 
     -             0.5D0*(s2b*tmp28)/mt)) + 
     -        cbe**2*hb**2*
     -         (0.125D0*(1 + c2b)/c2t - 0.125D0*(mb*s2b)/(mt*s2t) + 
     -           0.125D0*((1 + c2b)*Xb)/(mt*s2t)) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mb*s2b)/(mt*s2t) + 
     -           0.125D0*((1 - c2b)*Xt)/(mt*s2t)))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mb*tmp55)/mt - 0.25D0*s2b/s2t + 
     -             0.5D0*(s2b*tmp27)/mt)) + 
     -        cbe**2*hb**2*
     -         (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mb*s2b)/(mt*s2t) - 
     -           0.125D0*((1 - c2b)*Xb)/(mt*s2t)) + 
     -        ht**2*sbe**2*
     -         (0.125D0*(1 + c2b)/c2t - 0.125D0*(mb*s2b)/(mt*s2t) - 
     -           0.125D0*((1 + c2b)*Xt)/(mt*s2t)))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mb*tmp54)/mt + 0.25D0*s2b/s2t + 
     -             0.5D0*(s2b*tmp28)/mt)) + 
     -        cbe**2*hb**2*
     -         (0.125D0*(1 - c2b)/c2t + 0.125D0*(mb*s2b)/(mt*s2t) + 
     -           0.125D0*((1 - c2b)*Xb)/(mt*s2t)) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(1 + c2b)/c2t) + 0.125D0*(mb*s2b)/(mt*s2t) + 
     -           0.125D0*((1 + c2b)*Xt)/(mt*s2t))))
        Dtc2t = Dtc2t - tmp89*
     -    (cbe*hb*ht*sbe*((mb*tmp54)/mt + 0.25D0*s2b/s2t - 
     -         0.5D0*(s2b*tmp43)/mt) + 
     -      hb**2*sbe**2*(-(0.125D0*(1 + c2b)/c2t) + 
     -         0.125D0*(mb*s2b)/(mt*s2t) - 
     -         0.125D0*((1 + c2b)*Yb)/(mt*s2t)) + 
     -      cbe**2*ht**2*(0.125D0*(1 - c2b)/c2t + 
     -         0.125D0*(mb*s2b)/(mt*s2t) - 
     -         0.125D0*((1 - c2b)*Yt)/(mt*s2t))) - 
     -   tmp109*(cbe*hb*ht*sbe*
     -       ((mb*tmp55)/mt - 0.25D0*s2b/s2t - 0.5D0*(s2b*tmp44)/mt) + 
     -      hb**2*sbe**2*(0.125D0*(1 + c2b)/c2t - 
     -         0.125D0*(mb*s2b)/(mt*s2t) + 
     -         0.125D0*((1 + c2b)*Yb)/(mt*s2t)) + 
     -      cbe**2*ht**2*(-(0.125D0*(1 - c2b)/c2t) - 
     -         0.125D0*(mb*s2b)/(mt*s2t) + 
     -         0.125D0*((1 - c2b)*Yt)/(mt*s2t))) - 
     -   tmp88*(cbe*hb*ht*sbe*
     -       ((mb*tmp55)/mt - 0.25D0*s2b/s2t + 0.5D0*(s2b*tmp43)/mt) + 
     -      hb**2*sbe**2*(-(0.125D0*(1 - c2b)/c2t) - 
     -         0.125D0*(mb*s2b)/(mt*s2t) - 
     -         0.125D0*((1 - c2b)*Yb)/(mt*s2t)) + 
     -      cbe**2*ht**2*(0.125D0*(1 + c2b)/c2t - 
     -         0.125D0*(mb*s2b)/(mt*s2t) - 
     -         0.125D0*((1 + c2b)*Yt)/(mt*s2t))) - 
     -   tmp108*(cbe*hb*ht*sbe*
     -       ((mb*tmp54)/mt + 0.25D0*s2b/s2t + 0.5D0*(s2b*tmp44)/mt) + 
     -      hb**2*sbe**2*(0.125D0*(1 - c2b)/c2t + 
     -         0.125D0*(mb*s2b)/(mt*s2t) + 
     -         0.125D0*((1 - c2b)*Yb)/(mt*s2t)) + 
     -      cbe**2*ht**2*(-(0.125D0*(1 + c2b)/c2t) + 
     -         0.125D0*(mb*s2b)/(mt*s2t) + 
     -         0.125D0*((1 + c2b)*Yt)/(mt*s2t))) - 
     -   0.25D0*(ht**2*((cbe**2*tmp114*Yt)/(mt*s2t) - 
     -        (cbe**2*tmp93*Yt)/(mt*s2t) + 
     -        0.5D0*(sbe**2*tmp115*Xt)/(mt*s2t) - 
     -        0.5D0*(sbe**2*tmp94*Xt)/(mt*s2t)))
        DT1c2t = -(tmp204*tmp238) - tmp203*tmp239 + 
     -   (hb*ht*mb*mu*tmp240)/s2t + 
     -   hb**2*(0.125D0*(B1*(1 - c2b)*(-1 + LogB1))/c2t + 
     -      0.125D0*(B2*(1 + c2b)*(-1 + LogB2))/c2t + 
     -      sbe**2*(0.25D0*(A0*(-1 + LogA0))/c2t + 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT1))/c2t) + 
     -      0.125D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT1))/c2t + 
     -      0.125D0*(B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT1))/c2t) + 
     -   tmp6*(-(b*(-1 + Logb)) - 2*b*Logb - 
     -      b*(-1 + Logb)*(-1 + LogT1) - (-1 + Logmu2)*mu2 - 
     -      2*Logmu2*mu2 - (-1 + Logmu2)*(-1 + LogT1)*mu2 - 
     -      2*LogT1*T1 - (-b - mu2 + T1)*tmp240 + 
     -      0.5D0*(deltT1bmu2*phiT1bmu2)/mu2 - 
     -      0.5D0*(Logmu2*LogT1*(b - mu2 - T1)) - 
     -      0.5D0*(Logb*LogT1*(-b + mu2 - T1)) - 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T1)) + 2.5D0*(b + mu2 + T1))
     -     - 0.5D0*(tmp200*tmp213) - 0.5D0*(tmp199*tmp214) + 
     -   ht**2*(-(0.125D0*(B1*(1 + c2b)*(-1 + LogB1))/c2t) - 
     -      0.125D0*(B2*(1 - c2b)*(-1 + LogB2))/c2t + 
     -      cbe**2*(-(0.25D0*(A0*(-1 + LogA0))/c2t) - 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT1))/c2t) - 
     -      0.125D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT1))/c2t - 
     -      0.125D0*(B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT1))/c2t + 
     -      (-1 + LogT2)*T2*(1 + 0.5D0*(-1 + Nc)) + 
     -      (-1 + LogT1)*(-1 + LogT2)*T2*(1 + 0.5D0*(-1 + Nc)) - 
     -      0.25D0*((1 + Nc)*tmp79)) - 
     -   0.25D0*(ht**2*(sbe**2*tmp220*Xt**2 + 
     -        2*cbe**2*tmp221*Yt**2 + 
     -        cbe**2*tmp46*
     -         (-1 + 4*LogT1 - LogT1**2 - (A0*LogA0)/T1 + 
     -           (2*A0*phiA0T1T1)/T1 + (LogT1*(A0 - 2*T1))/T1 + 
     -           0.5D0*(deltA0T1T1*phiA0T1T1)/T1**2 - 
     -           0.5D0*(deltA0T1T1*
     -               ((A0*phiA0T1T1)/deltA0T1T1 + 
     -                 (phiA0T1T1*tmp178)/deltA0T1T1 + 
     -                 tmp80/deltA0T1T1 + tmp87/deltA0T1T1))/T1) + 
     -        0.5D0*(sbe**2*tmp30*tmp95)))
        DT2c2t = -(tmp206*tmp222) - tmp205*tmp223 + 
     -   (hb*ht*mb*mu*tmp225)/s2t + 
     -   hb**2*(-(0.125D0*(B1*(1 - c2b)*(-1 + LogB1))/c2t) - 
     -      0.125D0*(B2*(1 + c2b)*(-1 + LogB2))/c2t + 
     -      sbe**2*(-(0.25D0*(A0*(-1 + LogA0))/c2t) - 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT2))/c2t) - 
     -      0.125D0*(B1*(1 - c2b)*(-1 + LogB1)*(-1 + LogT2))/c2t - 
     -      0.125D0*(B2*(1 + c2b)*(-1 + LogB2)*(-1 + LogT2))/c2t) + 
     -   tmp7*(-(b*(-1 + Logb)) - 2*b*Logb - 
     -      b*(-1 + Logb)*(-1 + LogT2) - (-1 + Logmu2)*mu2 - 
     -      2*Logmu2*mu2 - (-1 + Logmu2)*(-1 + LogT2)*mu2 - 
     -      2*LogT2*T2 - (-b - mu2 + T2)*tmp224 + 
     -      0.5D0*(deltT2bmu2*phiT2bmu2)/mu2 - 
     -      0.5D0*(Logmu2*LogT2*(b - mu2 - T2)) - 
     -      0.5D0*(Logb*LogT2*(-b + mu2 - T2)) - 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T2)) + 2.5D0*(b + mu2 + T2))
     -     - 0.5D0*(tmp202*tmp217) - 0.5D0*(tmp201*tmp218) + 
     -   ht**2*(0.125D0*(B1*(1 + c2b)*(-1 + LogB1))/c2t + 
     -      0.125D0*(B2*(1 - c2b)*(-1 + LogB2))/c2t + 
     -      cbe**2*(0.25D0*(A0*(-1 + LogA0))/c2t + 
     -         0.25D0*(A0*(-1 + LogA0)*(-1 + LogT2))/c2t) + 
     -      0.125D0*(B1*(1 + c2b)*(-1 + LogB1)*(-1 + LogT2))/c2t + 
     -      0.125D0*(B2*(1 - c2b)*(-1 + LogB2)*(-1 + LogT2))/c2t + 
     -      (-1 + LogT1)*T1*(1 + 0.5D0*(-1 + Nc)) + 
     -      (-1 + LogT1)*(-1 + LogT2)*T1*(1 + 0.5D0*(-1 + Nc)) - 
     -      0.25D0*((1 + Nc)*tmp98)) - 
     -   0.25D0*(ht**2*(sbe**2*tmp219*Xt**2 + 
     -        2*cbe**2*Yt**2*
     -         (-0.5D0 + 2*LogT2 - (phiA0T1T2*(-A0 - T1 + T2))/T2 + 
     -           0.5D0*(LogA0*LogT1) - 0.5D0*(LogA0*LogT2) - 
     -           0.5D0*(LogT1*LogT2) + 
     -           0.5D0*(deltA0T1T2*phiA0T1T2)/T2**2 + 
     -           0.5D0*(LogT1*(A0 - T1 - T2))/T2 + 
     -           0.5D0*(LogA0*(-A0 + T1 - T2))/T2 - 
     -           0.5D0*(deltA0T1T2*
     -               (tmp104/deltA0T1T2 + 
     -                 (phiA0T1T2*tmp179)/deltA0T1T2))/T2) + 
     -        0.5D0*(sbe**2*tmp116*tmp31) + 
     -        cbe**2*tmp47*
     -         (-1 + 4*LogT2 - LogT2**2 - (A0*LogA0)/T2 + 
     -           (2*A0*phiA0T2T2)/T2 + (LogT2*(A0 - 2*T2))/T2 + 
     -           0.5D0*(deltA0T2T2*phiA0T2T2)/T2**2 - 
     -           0.5D0*(deltA0T2T2*
     -               ((A0*phiA0T2T2)/deltA0T2T2 + 
     -                 tmp107/deltA0T2T2 + 
     -                 (phiA0T2T2*tmp182)/deltA0T2T2 + 
     -                 tmp99/deltA0T2T2))/T2)))
        Dtb = tmp10*(-1 + Logb + (-1 + Logb)*(-1 + Logt) + 
     -      Logt + 0.5D0*((b + t)*tmp208)
     $       + 0.5D0*tmp209 + 0.5D0*tmp210)-
     -     tmp108*(-(0.125D0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)) - 
     -      0.125D0*(hb**2*s2b*s2t*sbe**2)/(mb*mt) + 
     -      0.5D0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)) - 
     -   tmp109*(0.125D0*(cbe**2*ht**2*s2b*s2t)/(mb*mt) + 
     -      0.125D0*(hb**2*s2b*s2t*sbe**2)/(mb*mt) + 
     -      0.5D0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)) - 
     -   (2*cbe*hb*ht*mt*sbe*
     -      (0.5D0 - 2*Logt + (phiA0bt*(-A0 - b + t))/t - 
     -        0.5D0*(LogA0*Logb) + 0.5D0*(LogA0*Logt) + 
     -        0.5D0*(Logb*Logt) - 0.5D0*(Logb*(A0 - b - t))/t - 
     -        0.5D0*(LogA0*(-A0 + b - t))/t + 0.5D0*tmp210 + 
     -        0.5D0*(deltA0bt*
     -            ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -              (b*tmp66)/(deltA0bt*t)))/b))/mb - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (-(0.125D0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)) - 
     -        0.125D0*(ht**2*s2b*s2t*sbe**2)/(mb*mt) - 
     -        0.5D0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (-(0.125D0*(cbe**2*hb**2*s2b*s2t)/(mb*mt)) - 
     -        0.125D0*(ht**2*s2b*s2t*sbe**2)/(mb*mt) - 
     -        0.5D0*(cbe*hb*ht*sbe*tmp56)/(mb*mt))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (0.125D0*(cbe**2*hb**2*s2b*s2t)/(mb*mt) + 
     -        0.125D0*(ht**2*s2b*s2t*sbe**2)/(mb*mt) - 
     -        0.5D0*(cbe*hb*ht*sbe*tmp59)/(mb*mt))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (0.125D0*(cbe**2*hb**2*s2b*s2t)/(mb*mt) + 
     -        0.125D0*(ht**2*s2b*s2t*sbe**2)/(mb*mt) - 
     -        0.5D0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)))
        Dtb = Dtb - tmp89*
     -    (-(0.125D0*(cbe**2*ht**2*s2b*s2t)/(mb*mt)) - 
     -      0.125D0*(hb**2*s2b*s2t*sbe**2)/(mb*mt) + 
     -      0.5D0*(cbe*hb*ht*sbe*tmp56)/(mb*mt)) - 
     -   tmp88*(0.125D0*(cbe**2*ht**2*s2b*s2t)/(mb*mt) + 
     -      0.125D0*(hb**2*s2b*s2t*sbe**2)/(mb*mt) + 
     -      0.5D0*(cbe*hb*ht*sbe*tmp59)/(mb*mt)) - 
     -   (2*cbe*hb*ht*mb*sbe*
     -      (0.5D0 - 2*Logb + (phiA0bt*(-A0 + b - t))/t + 
     -        0.5D0*(LogA0*Logb) - 0.5D0*(LogA0*Logt) + 
     -        0.5D0*(Logb*Logt) - 0.5D0*(Logt*(A0 - b - t))/b - 
     -        0.5D0*(deltA0bt*phiA0bt)/(b*t) - 
     -        0.5D0*(LogA0*(-A0 - b + t))/b + 0.5D0*tmp209 + 
     -        0.5D0*(deltA0bt*
     -            ((b*phiA0bt*tmp174)/(deltA0bt*t) + 
     -              tmp69/deltA0bt))/b))/mt - 
     -   4*cbe*hb*ht*mb*mt*sbe*
     -    (-(phiA0bt/t) - (phiA0bt*(-A0 - b + t))/(b*t) + 
     -      ((-A0 + b - t)*
     -         ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -           (b*tmp66)/(deltA0bt*t)))/b + 
     -      ((-A0 - b + t)*
     -         ((b*phiA0bt*tmp174)/(deltA0bt*t) + tmp69/deltA0bt))/
     -       b - 0.5D0*LogA0/b + 0.5D0*Logt/b - 0.5D0*LogA0/t + 
     -      0.5D0*Logb/t - 0.5D0*(A0 - b - t)/(b*t) + 0.5D0*tmp208 - 
     -      0.5D0*(deltA0bt*((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -            (b*tmp66)/(deltA0bt*t)))/b**2 + 
     -      0.5D0*(deltA0bt*((b*phiA0bt)/(deltA0bt*t) - 
     -            (2*b*phiA0bt*(-A0 + b - t)*(A0 + b - t))/
     -             (deltA0bt**2*t) + (b*tmp64)/(deltA0bt*t) + 
     -            tmp66/(deltA0bt*t) - 
     -            (2*b*(-A0 + b - t)*tmp66)/(deltA0bt**2*t) + 
     -            ((A0 + b - t)*
     -               ((b*phiA0bt*tmp174)/(deltA0bt*t) + 
     -                 tmp69/deltA0bt))/deltA0bt))/b) - 
     -   (cbe*hb*ht*sbe*(-2*A0*LogA0 - 2*b*Logb - 2*Logt*t - 
     -        0.5D0*(Logb*Logt*(A0 - b - t)) - 
     -        0.5D0*(LogA0*Logt*(-A0 + b - t)) + 
     -        0.5D0*(deltA0bt*phiA0bt)/t - 
     -        0.5D0*(LogA0*Logb*(-A0 - b + t)) + 2.5D0*(A0 + b + t) + 
     -        0.5D0*tmp76))/(mb*mt)
        Dtb = Dtb + tmp9*
     -    (-2 + 3*Logb + (-1 + Logb)*(-1 + Logt) + 3*Logt - 
     -      Logb*Logt - (phiA0bt*(-A0 + b - t))/t - 
     -      (phiA0bt*(-A0 - b + t))/t + 
     -      0.5D0*(Logt*(A0 - b - t))/b + 
     -      0.5D0*(deltA0bt*phiA0bt)/(b*t) + 
     -      0.5D0*(Logb*(A0 - b - t))/t + 
     -      0.5D0*(LogA0*(-A0 + b - t))/t + 
     -      0.5D0*(LogA0*(-A0 - b + t))/b - 
     -      0.5D0*(deltA0bt*((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -            (b*tmp66)/(deltA0bt*t)))/b - 
     -      0.5D0*(deltA0bt*((b*phiA0bt*tmp174)/(deltA0bt*t) + 
     -            tmp69/deltA0bt))/b - 
     -      (A0 - b - t)*(phiA0bt/t + 
     -         (phiA0bt*(-A0 - b + t))/(b*t) - 
     -         ((-A0 + b - t)*
     -            ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -              (b*tmp66)/(deltA0bt*t)))/b - 
     -         ((-A0 - b + t)*
     -            ((b*phiA0bt*tmp174)/(deltA0bt*t) + 
     -              tmp69/deltA0bt))/b + 0.5D0*LogA0/b - 
     -         0.5D0*Logt/b + 0.5D0*LogA0/t - 0.5D0*Logb/t + 
     -         0.5D0*(A0 - b - t)/(b*t) + 
     -         0.5D0*(deltA0bt*
     -             ((b*phiA0bt*(A0 + b - t))/(deltA0bt*t) + 
     -               (b*tmp66)/(deltA0bt*t)))/b**2 - 
     -         0.5D0*(deltA0bt*
     -             ((b*phiA0bt)/(deltA0bt*t) - 
     -               (2*b*phiA0bt*(-A0 + b - t)*(A0 + b - t))/
     -                (deltA0bt**2*t) + (b*tmp64)/(deltA0bt*t) + 
     -               tmp66/(deltA0bt*t) - 
     -               (2*b*(-A0 + b - t)*tmp66)/(deltA0bt**2*t) + 
     -               ((A0 + b - t)*
     -                  ((b*phiA0bt*tmp174)/(deltA0bt*t) + 
     -                    tmp69/deltA0bt))/deltA0bt))/b))
        DT1b = -((hb*ht*mu*s2t*tmp240)/mb) - 
     -   2*hb*ht*mb*mu*s2t*
     -    (phiT1bmu2/mu2 - 
     -      ((b - mu2 - T1)*
     -         ((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -           (mu2*tmp83)/(deltT1bmu2*T1)))/mu2 - 
     -      ((-b - mu2 + T1)*
     -         ((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -           (mu2*tmp90)/(b*deltT1bmu2)))/mu2 + 0.5D0*Logmu2/b - 
     -      0.5D0*LogT1/b - 0.5D0*Logb/T1 + 0.5D0*Logmu2/T1 + 
     -      0.5D0*(-b + mu2 - T1)/(b*T1) - 
     -      0.5D0*(deltT1bmu2*tmp242)/mu2) + 
     -   tmp62*(1 - 3*Logb + Logb*Logmu2 - 
     -      (-1 + Logb)*(-1 + LogT1) + LogT1 - Logmu2*LogT1 + 
     -      (phiT1bmu2*(b - mu2 - T1))/mu2 - 
     -      (phiT1bmu2*(-b - mu2 + T1))/mu2 - 
     -      0.5D0*(LogT1*(-b + mu2 - T1))/b + 
     -      0.5D0*(Logmu2*(b - mu2 - T1))/T1 + 
     -      0.5D0*(Logb*(-b + mu2 - T1))/T1 - 
     -      0.5D0*(Logmu2*(-b - mu2 + T1))/b - 
     -      (-b - mu2 + T1)*
     -       (phiT1bmu2/mu2 - 
     -         ((b - mu2 - T1)*
     -            ((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -              (mu2*tmp83)/(deltT1bmu2*T1)))/mu2 - 
     -         ((-b - mu2 + T1)*
     -            ((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -              (mu2*tmp90)/(b*deltT1bmu2)))/mu2 + 
     -         0.5D0*Logmu2/b - 0.5D0*LogT1/b - 0.5D0*Logb/T1 + 
     -         0.5D0*Logmu2/T1 + 0.5D0*(-b + mu2 - T1)/(b*T1) - 
     -         0.5D0*(deltT1bmu2*tmp242)/mu2) - 
     -      0.5D0*(deltT1bmu2*
     -          ((phiT1bmu2*(b + mu2 - T1))/deltT1bmu2 + 
     -            (mu2*tmp83)/(deltT1bmu2*T1)))/mu2 + 
     -      0.5D0*(deltT1bmu2*
     -          ((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -            (mu2*tmp90)/(b*deltT1bmu2)))/mu2) - 
     -   0.5D0*(tmp214*(-(cbe*hb*ht*sbe*
     -           ((mt*tmp56)/mb - 0.5D0*(s2b*s2t) + 
     -             0.5D0*(s2t*tmp128)/mb)) + 
     -        cbe**2*hb**2*
     -         (0.25D0*((1 - c2b)*(1 + c2t)) - 
     -           0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 + c2t)*s2b*Xb)/mb)
     -          + ht**2*sbe**2*
     -         (0.25D0*((1 + c2b)*(1 - c2t)) - 
     -           0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 - c2t)*s2b*Xt)/mb)
     -        ))
        DT1b = DT1b - tmp239*
     -    (cbe*hb*ht*sbe*((mt*tmp56)/mb - 0.5D0*(s2b*s2t) + 
     -         0.5D0*(s2t*tmp154)/mb) + 
     -      hb**2*sbe**2*(0.25D0*((1 - c2b)*(1 + c2t)) - 
     -         0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 + c2t)*s2b*Yb)/mb)+
     -        cbe**2*ht**2*
     -       (0.25D0*((1 + c2b)*(1 - c2t)) - 0.25D0*(mt*s2b*s2t)/mb - 
     -         0.25D0*((1 - c2t)*s2b*Yt)/mb)) - 
     -   tmp238*(cbe*hb*ht*sbe*
     -       ((mt*tmp59)/mb + 0.5D0*(s2b*s2t) + 0.5D0*(s2t*tmp155)/mb)+
     -        hb**2*sbe**2*
     -       (0.25D0*((1 + c2b)*(1 + c2t)) + 0.25D0*(mt*s2b*s2t)/mb + 
     -         0.25D0*((1 + c2t)*s2b*Yb)/mb) + 
     -      cbe**2*ht**2*(0.25D0*((1 - c2b)*(1 - c2t)) + 
     -         0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 - c2t)*s2b*Yt)/mb))-
     -     0.5D0*(tmp213*(-(cbe*hb*ht*sbe*
     -           ((mt*tmp59)/mb + 0.5D0*(s2b*s2t) + 
     -             0.5D0*(s2t*tmp129)/mb)) + 
     -        cbe**2*hb**2*
     -         (0.25D0*((1 + c2b)*(1 + c2t)) + 
     -           0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 + c2t)*s2b*Xb)/mb)
     -          + ht**2*sbe**2*
     -         (0.25D0*((1 - c2b)*(1 - c2t)) + 
     -           0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 - c2t)*s2b*Xt)/mb)
     -        ))
        DT2b = -((hb*ht*mu*s2t*tmp225)/mb) + 
     -   tmp63*(1 - 3*Logb + Logb*Logmu2 - 
     -      (-1 + Logb)*(-1 + LogT2) + LogT2 - Logmu2*LogT2 + 
     -      (phiT2bmu2*(b - mu2 - T2))/mu2 - 
     -      (phiT2bmu2*(-b - mu2 + T2))/mu2 - 
     -      0.5D0*(LogT2*(-b + mu2 - T2))/b + 
     -      0.5D0*(Logmu2*(b - mu2 - T2))/T2 + 
     -      0.5D0*(Logb*(-b + mu2 - T2))/T2 - 
     -      0.5D0*(Logmu2*(-b - mu2 + T2))/b - 
     -      0.5D0*(deltT2bmu2*
     -          ((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -            (mu2*tmp102)/(deltT2bmu2*T2)))/mu2 + 
     -      0.5D0*(deltT2bmu2*
     -          ((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -            (mu2*tmp110)/(b*deltT2bmu2)))/mu2 - 
     -      (-b - mu2 + T2)*
     -       (phiT2bmu2/mu2 - 
     -         ((b - mu2 - T2)*
     -            ((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -              (mu2*tmp102)/(deltT2bmu2*T2)))/mu2 - 
     -         ((-b - mu2 + T2)*
     -            ((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -              (mu2*tmp110)/(b*deltT2bmu2)))/mu2 + 
     -         0.5D0*Logmu2/b - 0.5D0*LogT2/b - 0.5D0*Logb/T2 + 
     -         0.5D0*Logmu2/T2 + 0.5D0*(-b + mu2 - T2)/(b*T2) - 
     -         0.5D0*(deltT2bmu2*tmp227)/mu2)) - 
     -   2*hb*ht*mb*mu*s2t*
     -    (-(phiT2bmu2/mu2) + 
     -      ((b - mu2 - T2)*
     -         ((phiT2bmu2*(b + mu2 - T2))/deltT2bmu2 + 
     -           (mu2*tmp102)/(deltT2bmu2*T2)))/mu2 + 
     -      ((-b - mu2 + T2)*
     -         ((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -           (mu2*tmp110)/(b*deltT2bmu2)))/mu2 - 
     -      0.5D0*Logmu2/b + 0.5D0*LogT2/b + 0.5D0*Logb/T2 - 
     -      0.5D0*Logmu2/T2 - 0.5D0*(-b + mu2 - T2)/(b*T2) + 
     -      0.5D0*(deltT2bmu2*tmp227)/mu2) - 
     -   0.5D0*(tmp218*(-(cbe*hb*ht*sbe*
     -           ((mt*tmp59)/mb + 0.5D0*(s2b*s2t) - 
     -             0.5D0*(s2t*tmp128)/mb)) + 
     -        cbe**2*hb**2*
     -         (0.25D0*((1 - c2b)*(1 - c2t)) + 
     -           0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 - c2t)*s2b*Xb)/mb)
     -          + ht**2*sbe**2*
     -         (0.25D0*((1 + c2b)*(1 + c2t)) + 
     -           0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 + c2t)*s2b*Xt)/mb)
     -        ))
        DT2b = DT2b - tmp223*
     -    (cbe*hb*ht*sbe*((mt*tmp59)/mb + 0.5D0*(s2b*s2t) - 
     -         0.5D0*(s2t*tmp154)/mb) + 
     -      hb**2*sbe**2*(0.25D0*((1 - c2b)*(1 - c2t)) + 
     -         0.25D0*(mt*s2b*s2t)/mb - 0.25D0*((1 - c2t)*s2b*Yb)/mb)+
     -        cbe**2*ht**2*
     -       (0.25D0*((1 + c2b)*(1 + c2t)) + 0.25D0*(mt*s2b*s2t)/mb - 
     -         0.25D0*((1 + c2t)*s2b*Yt)/mb)) - 
     -   tmp222*(cbe*hb*ht*sbe*
     -       ((mt*tmp56)/mb - 0.5D0*(s2b*s2t) - 0.5D0*(s2t*tmp155)/mb)+
     -        hb**2*sbe**2*
     -       (0.25D0*((1 + c2b)*(1 - c2t)) - 0.25D0*(mt*s2b*s2t)/mb + 
     -         0.25D0*((1 - c2t)*s2b*Yb)/mb) + 
     -      cbe**2*ht**2*(0.25D0*((1 - c2b)*(1 + c2t)) - 
     -         0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 + c2t)*s2b*Yt)/mb))-
     -     0.5D0*(tmp217*(-(cbe*hb*ht*sbe*
     -           ((mt*tmp56)/mb - 0.5D0*(s2b*s2t) - 
     -             0.5D0*(s2t*tmp129)/mb)) + 
     -        cbe**2*hb**2*
     -         (0.25D0*((1 + c2b)*(1 - c2t)) - 
     -           0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 - c2t)*s2b*Xb)/mb)
     -          + ht**2*sbe**2*
     -         (0.25D0*((1 - c2b)*(1 + c2t)) - 
     -           0.25D0*(mt*s2b*s2t)/mb + 0.25D0*((1 + c2t)*s2b*Xt)/mb)
     -        ))
        DB1t = -(tmp196*(-0.5D0 + 2*LogB1 - 
     -        (phiA0B1T2*(-A0 + B1 - T2))/T2 - 0.5D0*(LogA0*LogB1) + 
     -        0.5D0*(LogA0*LogT2) - 0.5D0*(LogB1*LogT2) + 
     -        0.5D0*(LogT2*(A0 - B1 - T2))/B1 + 
     -        0.5D0*(deltA0B1T2*phiA0B1T2)/(B1*T2) + 
     -        0.5D0*(LogA0*(-A0 - B1 + T2))/B1 - 
     -        0.5D0*(deltA0B1T2*
     -            (tmp105/deltA0B1T2 + 
     -              (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1)) - 
     -   0.5D0*(tmp184*tmp211) - 0.5D0*(tmp188*tmp215) - 
     -   2*hb*ht*mt*mu*s2b*
     -    (phiB1tmu2/mu2 - 
     -      ((B1 - mu2 - t)*
     -         ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -           (mu2*tmp67)/(deltB1tmu2*t)))/mu2 - 
     -      ((-B1 - mu2 + t)*
     -         ((phiB1tmu2*(-B1 + mu2 + t))/deltB1tmu2 + 
     -           (mu2*tmp71)/(B1*deltB1tmu2)))/mu2 + 
     -      0.5D0*Logmu2/B1 - 0.5D0*Logt/B1 - 0.5D0*LogB1/t + 
     -      0.5D0*Logmu2/t + 0.5D0*(-B1 + mu2 - t)/(B1*t) - 
     -      0.5D0*(deltB1tmu2*tmp233)/mu2) - 
     -   (hb*ht*mu*s2b*(-0.5D0 + 2*LogB1 - 
     -        (phiB1tmu2*(B1 - mu2 - t))/mu2 - 
     -        0.5D0*(LogB1*Logmu2) - 0.5D0*(LogB1*Logt) + 
     -        0.5D0*(Logmu2*Logt) + 0.5D0*(Logt*(-B1 + mu2 - t))/B1 + 
     -        0.5D0*(Logmu2*(-B1 - mu2 + t))/B1 - 
     -        0.5D0*(deltB1tmu2*
     -            ((phiB1tmu2*(-B1 + mu2 + t))/deltB1tmu2 + 
     -              (mu2*tmp71)/(B1*deltB1tmu2)))/mu2))/mt + 
     -   tmp61*(1 + LogB1 - LogB1*Logmu2 - 
     -      (-1 + LogB1)*(-1 + Logt) - 3*Logt + Logmu2*Logt - 
     -      (phiB1tmu2*(B1 - mu2 - t))/mu2 + 
     -      (phiB1tmu2*(-B1 - mu2 + t))/mu2 + 
     -      0.5D0*(Logt*(-B1 + mu2 - t))/B1 - 
     -      0.5D0*(Logmu2*(B1 - mu2 - t))/t - 
     -      0.5D0*(LogB1*(-B1 + mu2 - t))/t + 
     -      0.5D0*(Logmu2*(-B1 - mu2 + t))/B1 - 
     -      (B1 - mu2 - t)*
     -       (phiB1tmu2/mu2 - 
     -         ((B1 - mu2 - t)*
     -            ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -              (mu2*tmp67)/(deltB1tmu2*t)))/mu2 - 
     -         ((-B1 - mu2 + t)*
     -            ((phiB1tmu2*(-B1 + mu2 + t))/deltB1tmu2 + 
     -              (mu2*tmp71)/(B1*deltB1tmu2)))/mu2 + 
     -         0.5D0*Logmu2/B1 - 0.5D0*Logt/B1 - 0.5D0*LogB1/t + 
     -         0.5D0*Logmu2/t + 0.5D0*(-B1 + mu2 - t)/(B1*t) - 
     -         0.5D0*(deltB1tmu2*tmp233)/mu2) + 
     -      0.5D0*(deltB1tmu2*
     -          ((phiB1tmu2*(B1 + mu2 - t))/deltB1tmu2 + 
     -            (mu2*tmp67)/(deltB1tmu2*t)))/mu2 - 
     -      0.5D0*(deltB1tmu2*
     -          ((phiB1tmu2*(-B1 + mu2 + t))/deltB1tmu2 + 
     -            (mu2*tmp71)/(B1*deltB1tmu2)))/mu2)
        DB1t = DB1t - tmp192*
     -    (-0.5D0 + 2*LogB1 - (phiA0B1T1*(-A0 + B1 - T1))/T1 - 
     -      0.5D0*(LogA0*LogB1) + 0.5D0*(LogA0*LogT1) - 
     -      0.5D0*(LogB1*LogT1) + 0.5D0*(LogT1*(A0 - B1 - T1))/B1 + 
     -      0.5D0*(deltA0B1T1*phiA0B1T1)/(B1*T1) + 
     -      0.5D0*(LogA0*(-A0 - B1 + T1))/B1 - 
     -      0.5D0*(deltA0B1T1*
     -          ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1) + 
     -            tmp85/deltA0B1T1))/B1)
        DB2t = -(tmp195*(-0.5D0 + 2*LogB2 - 
     -        (phiA0B2T2*(-A0 + B2 - T2))/T2 - 0.5D0*(LogA0*LogB2) + 
     -        0.5D0*(LogA0*LogT2) - 0.5D0*(LogB2*LogT2) + 
     -        0.5D0*(LogT2*(A0 - B2 - T2))/B2 + 
     -        0.5D0*(deltA0B2T2*phiA0B2T2)/(B2*T2) + 
     -        0.5D0*(LogA0*(-A0 - B2 + T2))/B2 - 
     -        0.5D0*(deltA0B2T2*
     -            (tmp106/deltA0B2T2 + 
     -              (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2)) - 
     -   0.5D0*(tmp183*tmp212) - 0.5D0*(tmp187*tmp216) - 
     -   2*hb*ht*mt*mu*s2b*
     -    (-(phiB2tmu2/mu2) + 
     -      ((B2 - mu2 - t)*
     -         ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -           (mu2*tmp68)/(deltB2tmu2*t)))/mu2 + 
     -      ((-B2 - mu2 + t)*
     -         ((phiB2tmu2*(-B2 + mu2 + t))/deltB2tmu2 + 
     -           (mu2*tmp72)/(B2*deltB2tmu2)))/mu2 - 
     -      0.5D0*Logmu2/B2 + 0.5D0*Logt/B2 + 0.5D0*LogB2/t - 
     -      0.5D0*Logmu2/t - 0.5D0*(-B2 + mu2 - t)/(B2*t) + 
     -      0.5D0*(deltB2tmu2*tmp237)/mu2) + 
     -   tmp60*(1 + LogB2 - LogB2*Logmu2 - 
     -      (-1 + LogB2)*(-1 + Logt) - 3*Logt + Logmu2*Logt - 
     -      (phiB2tmu2*(B2 - mu2 - t))/mu2 + 
     -      (phiB2tmu2*(-B2 - mu2 + t))/mu2 + 
     -      0.5D0*(Logt*(-B2 + mu2 - t))/B2 - 
     -      0.5D0*(Logmu2*(B2 - mu2 - t))/t - 
     -      0.5D0*(LogB2*(-B2 + mu2 - t))/t + 
     -      0.5D0*(Logmu2*(-B2 - mu2 + t))/B2 - 
     -      (B2 - mu2 - t)*
     -       (phiB2tmu2/mu2 - 
     -         ((B2 - mu2 - t)*
     -            ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -              (mu2*tmp68)/(deltB2tmu2*t)))/mu2 - 
     -         ((-B2 - mu2 + t)*
     -            ((phiB2tmu2*(-B2 + mu2 + t))/deltB2tmu2 + 
     -              (mu2*tmp72)/(B2*deltB2tmu2)))/mu2 + 
     -         0.5D0*Logmu2/B2 - 0.5D0*Logt/B2 - 0.5D0*LogB2/t + 
     -         0.5D0*Logmu2/t + 0.5D0*(-B2 + mu2 - t)/(B2*t) - 
     -         0.5D0*(deltB2tmu2*tmp237)/mu2) + 
     -      0.5D0*(deltB2tmu2*
     -          ((phiB2tmu2*(B2 + mu2 - t))/deltB2tmu2 + 
     -            (mu2*tmp68)/(deltB2tmu2*t)))/mu2 - 
     -      0.5D0*(deltB2tmu2*
     -          ((phiB2tmu2*(-B2 + mu2 + t))/deltB2tmu2 + 
     -            (mu2*tmp72)/(B2*deltB2tmu2)))/mu2) - 
     -   (hb*ht*mu*s2b*(0.5D0 - 2*LogB2 + 
     -        (phiB2tmu2*(B2 - mu2 - t))/mu2 + 
     -        0.5D0*(LogB2*Logmu2) + 0.5D0*(LogB2*Logt) - 
     -        0.5D0*(Logmu2*Logt) - 0.5D0*(Logt*(-B2 + mu2 - t))/B2 - 
     -        0.5D0*(Logmu2*(-B2 - mu2 + t))/B2 + 
     -        0.5D0*(deltB2tmu2*
     -            ((phiB2tmu2*(-B2 + mu2 + t))/deltB2tmu2 + 
     -              (mu2*tmp72)/(B2*deltB2tmu2)))/mu2))/mt
        DB2t = DB2t - tmp191*
     -    (-0.5D0 + 2*LogB2 - (phiA0B2T1*(-A0 + B2 - T1))/T1 - 
     -      0.5D0*(LogA0*LogB2) + 0.5D0*(LogA0*LogT1) - 
     -      0.5D0*(LogB2*LogT1) + 0.5D0*(LogT1*(A0 - B2 - T1))/B2 + 
     -      0.5D0*(deltA0B2T1*phiA0B2T1)/(B2*T1) + 
     -      0.5D0*(LogA0*(-A0 - B2 + T1))/B2 - 
     -      0.5D0*(deltA0B2T1*
     -          ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1) + 
     -            tmp86/deltA0B2T1))/B2)
        DT1B1 = ht**2*(0.25D0*((1 + c2b)*(1 - c2t)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogB1)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogT1)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogB1)*(-1 + LogT1)))+
     -     hb**2*(0.25D0*((1 - c2b)*(1 + c2t)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogB1)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogT1)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogB1)*(-1 + LogT1)))-
     -     0.5D0*(((-2*B1*(LogB1 - LogT1))/((1 - B1/T1)*T1**2) + 
     -        (-2 - 2*LogB1)/T1 + (2*LogT1)/T1 - 
     -        (2*(LogB1 - LogT1))/((1 - B1/T1)*T1) + 
     -        (2*B1*(LogB1 - LogT1)*(-B1 + T1))/
     -         ((1 - B1/T1)**2*T1**3) + 
     -        (2*(-B1 + T1))/((1 - B1/T1)*T1**2) + 
     -        (2*(LogB1 - LogT1)*(-B1 + T1))/((1 - B1/T1)*T1**2))*
     -      tmp186) - tmp194*
     -    (phiA0B1T1/T1 + (phiA0B1T1*(-A0 - B1 + T1))/(B1*T1) - 
     -      ((-A0 + B1 - T1)*
     -         ((B1*phiA0B1T1*(A0 + B1 - T1))/(deltA0B1T1*T1) + 
     -           (B1*tmp81)/(deltA0B1T1*T1)))/B1 - 
     -      ((-A0 - B1 + T1)*
     -         ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1) + 
     -           tmp85/deltA0B1T1))/B1 + 0.5D0*LogA0/B1 - 
     -      0.5D0*LogT1/B1 + 0.5D0*LogA0/T1 - 0.5D0*LogB1/T1 + 
     -      0.5D0*(A0 - B1 - T1)/(B1*T1) + 
     -      0.5D0*(deltA0B1T1*
     -          ((B1*phiA0B1T1*(A0 + B1 - T1))/(deltA0B1T1*T1) + 
     -            (B1*tmp81)/(deltA0B1T1*T1)))/B1**2 - 
     -      0.5D0*(deltA0B1T1*
     -          ((B1*phiA0B1T1)/(deltA0B1T1*T1) - 
     -            (2*B1*phiA0B1T1*(-A0 + B1 - T1)*(A0 + B1 - T1))/
     -             (deltA0B1T1**2*T1) + 
     -            (B1*((A0 - B1)/B1 + LogA0 - LogB1 - T1/B1))/
     -             (deltA0B1T1*T1) + tmp81/(deltA0B1T1*T1) - 
     -            (2*B1*(-A0 + B1 - T1)*tmp81)/
     -             (deltA0B1T1**2*T1) + 
     -            ((A0 + B1 - T1)*
     -               ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1) + 
     -                 tmp85/deltA0B1T1))/deltA0B1T1))/B1)
        DT2B1 = hb**2*(0.25D0*((1 - c2b)*(1 - c2t)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogB1)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogT2)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogB1)*(-1 + LogT2)))+
     -     ht**2*(0.25D0*((1 + c2b)*(1 + c2t)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogB1)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogT2)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogB1)*(-1 + LogT2)))-
     -     tmp198*(phiA0B1T2/T2 + 
     -      (phiA0B1T2*(-A0 - B1 + T2))/(B1*T2) - 
     -      ((-A0 + B1 - T2)*
     -         ((B1*phiA0B1T2*(A0 + B1 - T2))/(deltA0B1T2*T2) + 
     -           (B1*tmp100)/(deltA0B1T2*T2)))/B1 - 
     -      ((-A0 - B1 + T2)*
     -         (tmp105/deltA0B1T2 + 
     -           (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1 + 
     -      0.5D0*LogA0/B1 - 0.5D0*LogT2/B1 + 0.5D0*LogA0/T2 - 
     -      0.5D0*LogB1/T2 + 0.5D0*(A0 - B1 - T2)/(B1*T2) + 
     -      0.5D0*(deltA0B1T2*
     -          ((B1*phiA0B1T2*(A0 + B1 - T2))/(deltA0B1T2*T2) + 
     -            (B1*tmp100)/(deltA0B1T2*T2)))/B1**2 - 
     -      0.5D0*(deltA0B1T2*
     -          ((B1*phiA0B1T2)/(deltA0B1T2*T2) - 
     -            (2*B1*phiA0B1T2*(-A0 + B1 - T2)*(A0 + B1 - T2))/
     -             (deltA0B1T2**2*T2) + 
     -            (B1*((A0 - B1)/B1 + LogA0 - LogB1 - T2/B1))/
     -             (deltA0B1T2*T2) + tmp100/(deltA0B1T2*T2) - 
     -            (2*B1*(-A0 + B1 - T2)*tmp100)/
     -             (deltA0B1T2**2*T2) + 
     -            ((A0 + B1 - T2)*
     -               (tmp105/deltA0B1T2 + 
     -                 (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/
     -             deltA0B1T2))/B1) - 
     -   0.5D0*(((-2*B1*(LogB1 - LogT2))/((1 - B1/T2)*T2**2) + 
     -        (-2 - 2*LogB1)/T2 + (2*LogT2)/T2 - 
     -        (2*(LogB1 - LogT2))/((1 - B1/T2)*T2) + 
     -        (2*B1*(LogB1 - LogT2)*(-B1 + T2))/
     -         ((1 - B1/T2)**2*T2**3) + 
     -        (2*(-B1 + T2))/((1 - B1/T2)*T2**2) + 
     -        (2*(LogB1 - LogT2)*(-B1 + T2))/((1 - B1/T2)*T2**2))*
     -      tmp190)
        DT1B2 = ht**2*(0.25D0*((1 - c2b)*(1 - c2t)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogB2)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogT1)) + 
     -      0.25D0*((1 - c2b)*(1 - c2t)*(-1 + LogB2)*(-1 + LogT1)))+
     -     hb**2*(0.25D0*((1 + c2b)*(1 + c2t)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogB2)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogT1)) + 
     -      0.25D0*((1 + c2b)*(1 + c2t)*(-1 + LogB2)*(-1 + LogT1)))-
     -     0.5D0*(((-2*B2*(LogB2 - LogT1))/((1 - B2/T1)*T1**2) + 
     -        (-2 - 2*LogB2)/T1 + (2*LogT1)/T1 - 
     -        (2*(LogB2 - LogT1))/((1 - B2/T1)*T1) + 
     -        (2*B2*(LogB2 - LogT1)*(-B2 + T1))/
     -         ((1 - B2/T1)**2*T1**3) + 
     -        (2*(-B2 + T1))/((1 - B2/T1)*T1**2) + 
     -        (2*(LogB2 - LogT1)*(-B2 + T1))/((1 - B2/T1)*T1**2))*
     -      tmp185) - tmp193*
     -    (phiA0B2T1/T1 + (phiA0B2T1*(-A0 - B2 + T1))/(B2*T1) - 
     -      ((-A0 + B2 - T1)*
     -         ((B2*phiA0B2T1*(A0 + B2 - T1))/(deltA0B2T1*T1) + 
     -           (B2*tmp82)/(deltA0B2T1*T1)))/B2 - 
     -      ((-A0 - B2 + T1)*
     -         ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1) + 
     -           tmp86/deltA0B2T1))/B2 + 0.5D0*LogA0/B2 - 
     -      0.5D0*LogT1/B2 + 0.5D0*LogA0/T1 - 0.5D0*LogB2/T1 + 
     -      0.5D0*(A0 - B2 - T1)/(B2*T1) + 
     -      0.5D0*(deltA0B2T1*
     -          ((B2*phiA0B2T1*(A0 + B2 - T1))/(deltA0B2T1*T1) + 
     -            (B2*tmp82)/(deltA0B2T1*T1)))/B2**2 - 
     -      0.5D0*(deltA0B2T1*
     -          ((B2*phiA0B2T1)/(deltA0B2T1*T1) - 
     -            (2*B2*phiA0B2T1*(-A0 + B2 - T1)*(A0 + B2 - T1))/
     -             (deltA0B2T1**2*T1) + 
     -            (B2*((A0 - B2)/B2 + LogA0 - LogB2 - T1/B2))/
     -             (deltA0B2T1*T1) + tmp82/(deltA0B2T1*T1) - 
     -            (2*B2*(-A0 + B2 - T1)*tmp82)/
     -             (deltA0B2T1**2*T1) + 
     -            ((A0 + B2 - T1)*
     -               ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1) + 
     -                 tmp86/deltA0B2T1))/deltA0B2T1))/B2)
        DT2B2 = hb**2*(0.25D0*((1 + c2b)*(1 - c2t)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogB2)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogT2)) + 
     -      0.25D0*((1 + c2b)*(1 - c2t)*(-1 + LogB2)*(-1 + LogT2)))+
     -     ht**2*(0.25D0*((1 - c2b)*(1 + c2t)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogB2)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogT2)) + 
     -      0.25D0*((1 - c2b)*(1 + c2t)*(-1 + LogB2)*(-1 + LogT2)))-
     -     tmp197*(phiA0B2T2/T2 + 
     -      (phiA0B2T2*(-A0 - B2 + T2))/(B2*T2) - 
     -      ((-A0 + B2 - T2)*
     -         ((B2*phiA0B2T2*(A0 + B2 - T2))/(deltA0B2T2*T2) + 
     -           (B2*tmp101)/(deltA0B2T2*T2)))/B2 - 
     -      ((-A0 - B2 + T2)*
     -         (tmp106/deltA0B2T2 + 
     -           (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2 + 
     -      0.5D0*LogA0/B2 - 0.5D0*LogT2/B2 + 0.5D0*LogA0/T2 - 
     -      0.5D0*LogB2/T2 + 0.5D0*(A0 - B2 - T2)/(B2*T2) + 
     -      0.5D0*(deltA0B2T2*
     -          ((B2*phiA0B2T2*(A0 + B2 - T2))/(deltA0B2T2*T2) + 
     -            (B2*tmp101)/(deltA0B2T2*T2)))/B2**2 - 
     -      0.5D0*(deltA0B2T2*
     -          ((B2*phiA0B2T2)/(deltA0B2T2*T2) - 
     -            (2*B2*phiA0B2T2*(-A0 + B2 - T2)*(A0 + B2 - T2))/
     -             (deltA0B2T2**2*T2) + 
     -            (B2*((A0 - B2)/B2 + LogA0 - LogB2 - T2/B2))/
     -             (deltA0B2T2*T2) + tmp101/(deltA0B2T2*T2) - 
     -            (2*B2*(-A0 + B2 - T2)*tmp101)/
     -             (deltA0B2T2**2*T2) + 
     -            ((A0 + B2 - T2)*
     -               (tmp106/deltA0B2T2 + 
     -                 (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/
     -             deltA0B2T2))/B2) - 
     -   0.5D0*(((-2*B2*(LogB2 - LogT2))/((1 - B2/T2)*T2**2) + 
     -        (-2 - 2*LogB2)/T2 + (2*LogT2)/T2 - 
     -        (2*(LogB2 - LogT2))/((1 - B2/T2)*T2) + 
     -        (2*B2*(LogB2 - LogT2)*(-B2 + T2))/
     -         ((1 - B2/T2)**2*T2**3) + 
     -        (2*(-B2 + T2))/((1 - B2/T2)*T2**2) + 
     -        (2*(LogB2 - LogT2)*(-B2 + T2))/((1 - B2/T2)*T2**2))*
     -      tmp189)
        Dbc2t = tmp7*(2*b*Logb + (-1 + Logmu2)*mu2 + 
     -      (-1 + Logb)*(-1 + Logmu2)*mu2 + 2*Logmu2*mu2 - 
     -      (-1 + LogT2)*T2 - (-1 + Logb)*(-1 + LogT2)*T2 + 
     -      2*LogT2*T2 - 0.5D0*(deltT2bmu2*phiT2bmu2)/mu2 + 
     -      0.5D0*(Logmu2*LogT2*(b - mu2 - T2)) + 
     -      0.5D0*(Logb*LogT2*(-b + mu2 - T2)) + 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T2)) - 
     -      2.5D0*(b + mu2 + T2) - 
     -      (-b - mu2 + T2)*
     -       (-0.5D0 + 2*Logb - (phiT2bmu2*(b - mu2 - T2))/mu2 - 
     -         0.5D0*(Logb*Logmu2) - 0.5D0*(Logb*LogT2) + 
     -         0.5D0*(Logmu2*LogT2) + 
     -         0.5D0*(LogT2*(-b + mu2 - T2))/b + 
     -         0.5D0*(Logmu2*(-b - mu2 + T2))/b - 
     -         0.5D0*(deltT2bmu2*
     -             ((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -               (mu2*tmp110)/(b*deltT2bmu2)))/mu2)) + 
     -   0.5D0*(hb*ht*mu*tmp113)/(mb*s2t) + 
     -   tmp6*(2*b*Logb + (-1 + Logmu2)*mu2 + 
     -      (-1 + Logb)*(-1 + Logmu2)*mu2 + 2*Logmu2*mu2 - 
     -      (-1 + LogT1)*T1 - (-1 + Logb)*(-1 + LogT1)*T1 + 
     -      2*LogT1*T1 - 0.5D0*(deltT1bmu2*phiT1bmu2)/mu2 + 
     -      0.5D0*(Logmu2*LogT1*(b - mu2 - T1)) + 
     -      0.5D0*(Logb*LogT1*(-b + mu2 - T1)) + 
     -      0.5D0*(Logb*Logmu2*(-b - mu2 + T1)) - 
     -      2.5D0*(b + mu2 + T1) - 
     -      (-b - mu2 + T1)*
     -       (-0.5D0 + 2*Logb - (phiT1bmu2*(b - mu2 - T1))/mu2 - 
     -         0.5D0*(Logb*Logmu2) - 0.5D0*(Logb*LogT1) + 
     -         0.5D0*(Logmu2*LogT1) + 
     -         0.5D0*(LogT1*(-b + mu2 - T1))/b + 
     -         0.5D0*(Logmu2*(-b - mu2 + T1))/b - 
     -         0.5D0*(deltT1bmu2*
     -             ((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -               (mu2*tmp90)/(b*deltT1bmu2)))/mu2)) + 
     -   (hb*ht*mb*mu*(-((phiT1bmu2*(b - mu2 - T1))/mu2) + 
     -        (phiT2bmu2*(b - mu2 - T2))/mu2 - 0.5D0*(Logb*LogT1) + 
     -        0.5D0*(Logmu2*LogT1) + 0.5D0*(Logb*LogT2) - 
     -        0.5D0*(Logmu2*LogT2) + 0.5D0*(LogT1*(-b + mu2 - T1))/b + 
     -        0.5D0*(Logmu2*(-b - mu2 + T1))/b - 
     -        0.5D0*(LogT2*(-b + mu2 - T2))/b - 
     -        0.5D0*(Logmu2*(-b - mu2 + T2))/b + 
     -        0.5D0*(deltT2bmu2*
     -            ((phiT2bmu2*(-b + mu2 + T2))/deltT2bmu2 + 
     -              (mu2*tmp110)/(b*deltT2bmu2)))/mu2 - 
     -        0.5D0*(deltT1bmu2*
     -            ((phiT1bmu2*(-b + mu2 + T1))/deltT1bmu2 + 
     -              (mu2*tmp90)/(b*deltT1bmu2)))/mu2))/s2t
        Dbc2t = Dbc2t - 0.5D0*
     -    ((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mt*tmp55)/mb - 0.25D0*s2b/s2t - 
     -             0.25D0*tmp129/(mb*s2t))) + 
     -        cbe**2*hb**2*
     -         (0.125D0*(1 + c2b)/c2t - 0.125D0*(mt*s2b)/(mb*s2t) + 
     -           0.125D0*(s2b*Xb)/(c2t*mb)) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mt*s2b)/(mb*s2t) - 
     -           0.125D0*(s2b*Xt)/(c2t*mb)))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mt*tmp55)/mb - 0.25D0*s2b/s2t + 
     -             0.25D0*tmp128/(mb*s2t))) + 
     -        cbe**2*hb**2*
     -         (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mt*s2b)/(mb*s2t) + 
     -           0.125D0*(s2b*Xb)/(c2t*mb)) + 
     -        ht**2*sbe**2*
     -         (0.125D0*(1 + c2b)/c2t - 0.125D0*(mt*s2b)/(mb*s2t) - 
     -           0.125D0*(s2b*Xt)/(c2t*mb)))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mt*tmp54)/mb + 0.25D0*s2b/s2t + 
     -             0.25D0*tmp129/(mb*s2t))) + 
     -        cbe**2*hb**2*
     -         (-(0.125D0*(1 + c2b)/c2t) + 0.125D0*(mt*s2b)/(mb*s2t) - 
     -           0.125D0*(s2b*Xb)/(c2t*mb)) + 
     -        ht**2*sbe**2*
     -         (0.125D0*(1 - c2b)/c2t + 0.125D0*(mt*s2b)/(mb*s2t) + 
     -           0.125D0*(s2b*Xt)/(c2t*mb)))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (-(cbe*hb*ht*sbe*
     -           ((mt*tmp54)/mb + 0.25D0*s2b/s2t - 
     -             0.25D0*tmp128/(mb*s2t))) + 
     -        cbe**2*hb**2*
     -         (0.125D0*(1 - c2b)/c2t + 0.125D0*(mt*s2b)/(mb*s2t) - 
     -           0.125D0*(s2b*Xb)/(c2t*mb)) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(1 + c2b)/c2t) + 0.125D0*(mt*s2b)/(mb*s2t) + 
     -           0.125D0*(s2b*Xt)/(c2t*mb))))
        Dbc2t = Dbc2t - tmp88*
     -    (cbe*hb*ht*sbe*((mt*tmp55)/mb - 0.25D0*s2b/s2t - 
     -         0.25D0*tmp155/(mb*s2t)) + 
     -      hb**2*sbe**2*(0.125D0*(1 + c2b)/c2t - 
     -         0.125D0*(mt*s2b)/(mb*s2t) + 0.125D0*(s2b*Yb)/(c2t*mb))+
     -        cbe**2*ht**2*
     -       (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mt*s2b)/(mb*s2t) - 
     -         0.125D0*(s2b*Yt)/(c2t*mb))) - 
     -   tmp109*(cbe*hb*ht*sbe*
     -       ((mt*tmp55)/mb - 0.25D0*s2b/s2t + 0.25D0*tmp154/(mb*s2t))+
     -        hb**2*sbe**2*
     -       (-(0.125D0*(1 - c2b)/c2t) - 0.125D0*(mt*s2b)/(mb*s2t) + 
     -         0.125D0*(s2b*Yb)/(c2t*mb)) + 
     -      cbe**2*ht**2*(0.125D0*(1 + c2b)/c2t - 
     -         0.125D0*(mt*s2b)/(mb*s2t) - 0.125D0*(s2b*Yt)/(c2t*mb)))-
     -     tmp108*(cbe*hb*ht*sbe*
     -       ((mt*tmp54)/mb + 0.25D0*s2b/s2t + 0.25D0*tmp155/(mb*s2t))+
     -        hb**2*sbe**2*
     -       (-(0.125D0*(1 + c2b)/c2t) + 0.125D0*(mt*s2b)/(mb*s2t) - 
     -         0.125D0*(s2b*Yb)/(c2t*mb)) + 
     -      cbe**2*ht**2*(0.125D0*(1 - c2b)/c2t + 
     -         0.125D0*(mt*s2b)/(mb*s2t) + 0.125D0*(s2b*Yt)/(c2t*mb)))-
     -     tmp89*(cbe*hb*ht*sbe*
     -       ((mt*tmp54)/mb + 0.25D0*s2b/s2t - 0.25D0*tmp154/(mb*s2t))+
     -        hb**2*sbe**2*
     -       (0.125D0*(1 - c2b)/c2t + 0.125D0*(mt*s2b)/(mb*s2t) - 
     -         0.125D0*(s2b*Yb)/(c2t*mb)) + 
     -      cbe**2*ht**2*(-(0.125D0*(1 + c2b)/c2t) + 
     -         0.125D0*(mt*s2b)/(mb*s2t) + 0.125D0*(s2b*Yt)/(c2t*mb)))
        DB1c2t = hb**2*(0.125D0*((1 - c2b)*(-1 + LogT1)*T1)/c2t + 
     -      0.125D0*((1 - c2b)*(-1 + LogB1)*(-1 + LogT1)*T1)/c2t - 
     -      0.125D0*((1 - c2b)*(-1 + LogT2)*T2)/c2t - 
     -      0.125D0*((1 - c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/c2t) + 
     -   ht**2*(-(0.125D0*((1 + c2b)*(-1 + LogT1)*T1)/c2t) - 
     -      0.125D0*((1 + c2b)*(-1 + LogB1)*(-1 + LogT1)*T1)/c2t + 
     -      0.125D0*((1 + c2b)*(-1 + LogT2)*T2)/c2t + 
     -      0.125D0*((1 + c2b)*(-1 + LogB1)*(-1 + LogT2)*T2)/c2t) - 
     -   tmp206*(-0.5D0 + 2*LogB1 - (phiA0B1T2*(-A0 + B1 - T2))/T2 - 
     -      0.5D0*(LogA0*LogB1) + 0.5D0*(LogA0*LogT2) - 
     -      0.5D0*(LogB1*LogT2) + 0.5D0*(LogT2*(A0 - B1 - T2))/B1 + 
     -      0.5D0*(deltA0B1T2*phiA0B1T2)/(B1*T2) + 
     -      0.5D0*(LogA0*(-A0 - B1 + T2))/B1 - 
     -      0.5D0*(deltA0B1T2*
     -          (tmp105/deltA0B1T2 + 
     -            (B1*phiA0B1T2*tmp180)/(deltA0B1T2*T2)))/B1) - 
     -   0.5D0*(tmp200*tmp211) - 0.5D0*(tmp202*tmp215) - 
     -   tmp204*(-0.5D0 + 2*LogB1 - (phiA0B1T1*(-A0 + B1 - T1))/T1 - 
     -      0.5D0*(LogA0*LogB1) + 0.5D0*(LogA0*LogT1) - 
     -      0.5D0*(LogB1*LogT1) + 0.5D0*(LogT1*(A0 - B1 - T1))/B1 + 
     -      0.5D0*(deltA0B1T1*phiA0B1T1)/(B1*T1) + 
     -      0.5D0*(LogA0*(-A0 - B1 + T1))/B1 - 
     -      0.5D0*(deltA0B1T1*
     -          ((B1*phiA0B1T1*tmp176)/(deltA0B1T1*T1) + 
     -            tmp85/deltA0B1T1))/B1)
        DB2c2t = ht**2*(-(0.125D0*
     -         ((1 - c2b)*(-1 + LogT1)*T1)/c2t) - 
     -      0.125D0*((1 - c2b)*(-1 + LogB2)*(-1 + LogT1)*T1)/c2t + 
     -      0.125D0*((1 - c2b)*(-1 + LogT2)*T2)/c2t + 
     -      0.125D0*((1 - c2b)*(-1 + LogB2)*(-1 + LogT2)*T2)/c2t) + 
     -   hb**2*(0.125D0*((1 + c2b)*(-1 + LogT1)*T1)/c2t + 
     -      0.125D0*((1 + c2b)*(-1 + LogB2)*(-1 + LogT1)*T1)/c2t - 
     -      0.125D0*((1 + c2b)*(-1 + LogT2)*T2)/c2t - 
     -      0.125D0*((1 + c2b)*(-1 + LogB2)*(-1 + LogT2)*T2)/c2t) - 
     -   tmp205*(-0.5D0 + 2*LogB2 - (phiA0B2T2*(-A0 + B2 - T2))/T2 - 
     -      0.5D0*(LogA0*LogB2) + 0.5D0*(LogA0*LogT2) - 
     -      0.5D0*(LogB2*LogT2) + 0.5D0*(LogT2*(A0 - B2 - T2))/B2 + 
     -      0.5D0*(deltA0B2T2*phiA0B2T2)/(B2*T2) + 
     -      0.5D0*(LogA0*(-A0 - B2 + T2))/B2 - 
     -      0.5D0*(deltA0B2T2*
     -          (tmp106/deltA0B2T2 + 
     -            (B2*phiA0B2T2*tmp181)/(deltA0B2T2*T2)))/B2) - 
     -   0.5D0*(tmp199*tmp212) - 0.5D0*(tmp201*tmp216) - 
     -   tmp203*(-0.5D0 + 2*LogB2 - (phiA0B2T1*(-A0 + B2 - T1))/T1 - 
     -      0.5D0*(LogA0*LogB2) + 0.5D0*(LogA0*LogT1) - 
     -      0.5D0*(LogB2*LogT1) + 0.5D0*(LogT1*(A0 - B2 - T1))/B2 + 
     -      0.5D0*(deltA0B2T1*phiA0B2T1)/(B2*T1) + 
     -      0.5D0*(LogA0*(-A0 - B2 + T1))/B2 - 
     -      0.5D0*(deltA0B2T1*
     -          ((B2*phiA0B2T1*tmp177)/(deltA0B2T1*T1) + 
     -            tmp86/deltA0B2T1))/B2)
        DT1c2b = ht**2*(0.125D0*(B1*(1 - c2t)*(-1 + LogB1))/c2b - 
     -      0.125D0*(B2*(1 - c2t)*(-1 + LogB2))/c2b + 
     -      0.125D0*(B1*(1 - c2t)*(-1 + LogB1)*(-1 + LogT1))/c2b - 
     -      0.125D0*(B2*(1 - c2t)*(-1 + LogB2)*(-1 + LogT1))/c2b) + 
     -   hb**2*(-(0.125D0*(B1*(1 + c2t)*(-1 + LogB1))/c2b) + 
     -      0.125D0*(B2*(1 + c2t)*(-1 + LogB2))/c2b - 
     -      0.125D0*(B1*(1 + c2t)*(-1 + LogB1)*(-1 + LogT1))/c2b + 
     -      0.125D0*(B2*(1 + c2t)*(-1 + LogB2)*(-1 + LogT1))/c2b) - 
     -   0.5D0*(tmp214*(cbe**2*hb**2*
     -         (-(0.125D0*(b*(1 + c2t))/c2b) + 
     -           0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 - c2t)*t)/c2b + 
     -           0.25D0*((1 + c2t)*mb*Xb)/s2b + 
     -           0.25D0*(mt*s2t*Xb)/c2b + 0.125D0*((1 + c2t)*Xb**2)/c2b
     -           ) - cbe*hb*ht*sbe*
     -         (mb*s2t*tmp23 + 2*mb*mt*tmp57 + 
     -           0.25D0*(s2t*(b + t))/s2b + 0.5D0*(mt*tmp130)/s2b + 
     -           0.25D0*(s2t*Xb*Xt)/s2b) + 
     -        ht**2*sbe**2*
     -         (0.125D0*(b*(1 - c2t))/c2b + 0.25D0*(mb*mt*s2t)/s2b - 
     -           0.125D0*((1 + c2t)*t)/c2b + 
     -           0.25D0*((1 - c2t)*mb*Xt)/s2b - 
     -           0.25D0*(mt*s2t*Xt)/c2b - 0.125D0*((1 - c2t)*Xt**2)/c2b
     -           ))) - 0.5D0*
     -    (tmp213*(cbe**2*hb**2*
     -         (0.125D0*(b*(1 + c2t))/c2b - 0.25D0*(mb*mt*s2t)/s2b - 
     -           0.125D0*((1 - c2t)*t)/c2b - 
     -           0.25D0*((1 + c2t)*mb*Xb)/s2b - 
     -           0.25D0*(mt*s2t*Xb)/c2b - 0.125D0*((1 + c2t)*Xb**2)/c2b
     -           ) - cbe*hb*ht*sbe*
     -         (mb*s2t*tmp24 + 2*mb*mt*tmp58 - 
     -           0.25D0*(s2t*(b + t))/s2b - 0.5D0*(mt*tmp130)/s2b - 
     -           0.25D0*(s2t*Xb*Xt)/s2b) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(b*(1 - c2t))/c2b) - 
     -           0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 + c2t)*t)/c2b - 
     -           0.25D0*((1 - c2t)*mb*Xt)/s2b + 
     -           0.25D0*(mt*s2t*Xt)/c2b + 0.125D0*((1 - c2t)*Xt**2)/c2b
     -           )))
        DT1c2b = DT1c2b - 
     -   tmp239*(hb**2*sbe**2*
     -       (-(0.125D0*(b*(1 + c2t))/c2b) + 0.25D0*(mb*mt*s2t)/s2b + 
     -         0.125D0*((1 - c2t)*t)/c2b + 
     -         0.25D0*((1 + c2t)*mb*Yb)/s2b + 0.25D0*(mt*s2t*Yb)/c2b + 
     -         0.125D0*((1 + c2t)*Yb**2)/c2b) + 
     -      cbe*hb*ht*sbe*(mb*s2t*tmp39 + 2*mb*mt*tmp57 + 
     -         0.25D0*(s2t*(b + t))/s2b + 0.5D0*(mt*tmp156)/s2b + 
     -         0.25D0*(s2t*Yb*Yt)/s2b) + 
     -      cbe**2*ht**2*(0.125D0*(b*(1 - c2t))/c2b + 
     -         0.25D0*(mb*mt*s2t)/s2b - 0.125D0*((1 + c2t)*t)/c2b + 
     -         0.25D0*((1 - c2t)*mb*Yt)/s2b - 0.25D0*(mt*s2t*Yt)/c2b - 
     -         0.125D0*((1 - c2t)*Yt**2)/c2b)) - 
     -   tmp238*(hb**2*sbe**2*
     -       (0.125D0*(b*(1 + c2t))/c2b - 0.25D0*(mb*mt*s2t)/s2b - 
     -         0.125D0*((1 - c2t)*t)/c2b - 
     -         0.25D0*((1 + c2t)*mb*Yb)/s2b - 0.25D0*(mt*s2t*Yb)/c2b - 
     -         0.125D0*((1 + c2t)*Yb**2)/c2b) + 
     -      cbe*hb*ht*sbe*(mb*s2t*tmp40 + 2*mb*mt*tmp58 - 
     -         0.25D0*(s2t*(b + t))/s2b - 0.5D0*(mt*tmp156)/s2b - 
     -         0.25D0*(s2t*Yb*Yt)/s2b) + 
     -      cbe**2*ht**2*(-(0.125D0*(b*(1 - c2t))/c2b) - 
     -         0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 + c2t)*t)/c2b - 
     -         0.25D0*((1 - c2t)*mb*Yt)/s2b + 0.25D0*(mt*s2t*Yt)/c2b + 
     -         0.125D0*((1 - c2t)*Yt**2)/c2b))
        DT2c2b = hb**2*(-(0.125D0*
     -         (B1*(1 - c2t)*(-1 + LogB1))/c2b) + 
     -      0.125D0*(B2*(1 - c2t)*(-1 + LogB2))/c2b - 
     -      0.125D0*(B1*(1 - c2t)*(-1 + LogB1)*(-1 + LogT2))/c2b + 
     -      0.125D0*(B2*(1 - c2t)*(-1 + LogB2)*(-1 + LogT2))/c2b) + 
     -   ht**2*(0.125D0*(B1*(1 + c2t)*(-1 + LogB1))/c2b - 
     -      0.125D0*(B2*(1 + c2t)*(-1 + LogB2))/c2b + 
     -      0.125D0*(B1*(1 + c2t)*(-1 + LogB1)*(-1 + LogT2))/c2b - 
     -      0.125D0*(B2*(1 + c2t)*(-1 + LogB2)*(-1 + LogT2))/c2b) - 
     -   0.5D0*(tmp218*(cbe**2*hb**2*
     -         (-(0.125D0*(b*(1 - c2t))/c2b) - 
     -           0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 + c2t)*t)/c2b + 
     -           0.25D0*((1 - c2t)*mb*Xb)/s2b - 
     -           0.25D0*(mt*s2t*Xb)/c2b + 0.125D0*((1 - c2t)*Xb**2)/c2b
     -           ) - cbe*hb*ht*sbe*
     -         (-(mb*s2t*tmp23) + 2*mb*mt*tmp58 - 
     -           0.25D0*(s2t*(b + t))/s2b + 0.5D0*(mt*tmp131)/s2b - 
     -           0.25D0*(s2t*Xb*Xt)/s2b) + 
     -        ht**2*sbe**2*
     -         (0.125D0*(b*(1 + c2t))/c2b - 0.25D0*(mb*mt*s2t)/s2b - 
     -           0.125D0*((1 - c2t)*t)/c2b + 
     -           0.25D0*((1 + c2t)*mb*Xt)/s2b + 
     -           0.25D0*(mt*s2t*Xt)/c2b - 0.125D0*((1 + c2t)*Xt**2)/c2b
     -           ))) - 0.5D0*
     -    (tmp217*(cbe**2*hb**2*
     -         (0.125D0*(b*(1 - c2t))/c2b + 0.25D0*(mb*mt*s2t)/s2b - 
     -           0.125D0*((1 + c2t)*t)/c2b - 
     -           0.25D0*((1 - c2t)*mb*Xb)/s2b + 
     -           0.25D0*(mt*s2t*Xb)/c2b - 0.125D0*((1 - c2t)*Xb**2)/c2b
     -           ) - cbe*hb*ht*sbe*
     -         (-(mb*s2t*tmp24) + 2*mb*mt*tmp57 + 
     -           0.25D0*(s2t*(b + t))/s2b - 0.5D0*(mt*tmp131)/s2b + 
     -           0.25D0*(s2t*Xb*Xt)/s2b) + 
     -        ht**2*sbe**2*
     -         (-(0.125D0*(b*(1 + c2t))/c2b) + 
     -           0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 - c2t)*t)/c2b - 
     -           0.25D0*((1 + c2t)*mb*Xt)/s2b - 
     -           0.25D0*(mt*s2t*Xt)/c2b + 0.125D0*((1 + c2t)*Xt**2)/c2b
     -           )))
        DT2c2b = DT2c2b - 
     -   tmp223*(hb**2*sbe**2*
     -       (-(0.125D0*(b*(1 - c2t))/c2b) - 0.25D0*(mb*mt*s2t)/s2b + 
     -         0.125D0*((1 + c2t)*t)/c2b + 
     -         0.25D0*((1 - c2t)*mb*Yb)/s2b - 0.25D0*(mt*s2t*Yb)/c2b + 
     -         0.125D0*((1 - c2t)*Yb**2)/c2b) + 
     -      cbe*hb*ht*sbe*(-(mb*s2t*tmp39) + 2*mb*mt*tmp58 - 
     -         0.25D0*(s2t*(b + t))/s2b + 0.5D0*(mt*tmp157)/s2b - 
     -         0.25D0*(s2t*Yb*Yt)/s2b) + 
     -      cbe**2*ht**2*(0.125D0*(b*(1 + c2t))/c2b - 
     -         0.25D0*(mb*mt*s2t)/s2b - 0.125D0*((1 - c2t)*t)/c2b + 
     -         0.25D0*((1 + c2t)*mb*Yt)/s2b + 0.25D0*(mt*s2t*Yt)/c2b - 
     -         0.125D0*((1 + c2t)*Yt**2)/c2b)) - 
     -   tmp222*(hb**2*sbe**2*
     -       (0.125D0*(b*(1 - c2t))/c2b + 0.25D0*(mb*mt*s2t)/s2b - 
     -         0.125D0*((1 + c2t)*t)/c2b - 
     -         0.25D0*((1 - c2t)*mb*Yb)/s2b + 0.25D0*(mt*s2t*Yb)/c2b - 
     -         0.125D0*((1 - c2t)*Yb**2)/c2b) + 
     -      cbe*hb*ht*sbe*(-(mb*s2t*tmp40) + 2*mb*mt*tmp57 + 
     -         0.25D0*(s2t*(b + t))/s2b - 0.5D0*(mt*tmp157)/s2b + 
     -         0.25D0*(s2t*Yb*Yt)/s2b) + 
     -      cbe**2*ht**2*(-(0.125D0*(b*(1 + c2t))/c2b) + 
     -         0.25D0*(mb*mt*s2t)/s2b + 0.125D0*((1 - c2t)*t)/c2b - 
     -         0.25D0*((1 + c2t)*mb*Yt)/s2b - 0.25D0*(mt*s2t*Yt)/c2b + 
     -         0.125D0*((1 + c2t)*Yt**2)/c2b))
        Dc2tc2b = hb**2*tmp19 + ht**2*tmp19 - 
     -   tmp89*(hb**2*sbe**2*tmp38 + cbe**2*ht**2*tmp49 + 
     -      cbe*hb*ht*sbe*(-(0.25D0*(mb*mt)/(c2b*c2t)) - 
     -         0.125D0*(b + t)/(s2b*s2t) - 0.5D0*(mb*tmp39)/s2t + 
     -         0.5D0*(mt*tmp43)/s2b - 0.125D0*(Yb*Yt)/(s2b*s2t))) - 
     -   tmp108*(hb**2*sbe**2*tmp38 + cbe**2*ht**2*tmp49 + 
     -      cbe*hb*ht*sbe*(-(0.25D0*(mb*mt)/(c2b*c2t)) - 
     -         0.125D0*(b + t)/(s2b*s2t) + 0.5D0*(mb*tmp40)/s2t - 
     -         0.5D0*(mt*tmp44)/s2b - 0.125D0*(Yb*Yt)/(s2b*s2t))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (cbe**2*hb**2*tmp22 + ht**2*sbe**2*tmp33 - 
     -        cbe*hb*ht*sbe*
     -         (-(0.25D0*(mb*mt)/(c2b*c2t)) - 
     -           0.125D0*(b + t)/(s2b*s2t) - 0.5D0*(mb*tmp23)/s2t + 
     -           0.5D0*(mt*tmp27)/s2b - 0.125D0*(Xb*Xt)/(s2b*s2t)))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (cbe**2*hb**2*tmp22 + ht**2*sbe**2*tmp33 - 
     -        cbe*hb*ht*sbe*
     -         (-(0.25D0*(mb*mt)/(c2b*c2t)) - 
     -           0.125D0*(b + t)/(s2b*s2t) + 0.5D0*(mb*tmp24)/s2t - 
     -           0.5D0*(mt*tmp28)/s2b - 0.125D0*(Xb*Xt)/(s2b*s2t)))) - 
     -   0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (cbe**2*hb**2*tmp21 + ht**2*sbe**2*tmp32 - 
     -        cbe*hb*ht*sbe*
     -         (0.25D0*(mb*mt)/(c2b*c2t) + 0.125D0*(b + t)/(s2b*s2t) - 
     -           0.5D0*(mb*tmp24)/s2t - 0.5D0*(mt*tmp27)/s2b + 
     -           0.125D0*(Xb*Xt)/(s2b*s2t)))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (cbe**2*hb**2*tmp21 + ht**2*sbe**2*tmp32 - 
     -        cbe*hb*ht*sbe*
     -         (0.25D0*(mb*mt)/(c2b*c2t) + 0.125D0*(b + t)/(s2b*s2t) + 
     -           0.5D0*(mb*tmp23)/s2t + 0.5D0*(mt*tmp28)/s2b + 
     -           0.125D0*(Xb*Xt)/(s2b*s2t))))
        Dc2tc2b = Dc2tc2b - 
     -   tmp88*(hb**2*sbe**2*tmp37 + cbe**2*ht**2*tmp48 + 
     -      cbe*hb*ht*sbe*(0.25D0*(mb*mt)/(c2b*c2t) + 
     -         0.125D0*(b + t)/(s2b*s2t) - 0.5D0*(mb*tmp40)/s2t - 
     -         0.5D0*(mt*tmp43)/s2b + 0.125D0*(Yb*Yt)/(s2b*s2t))) - 
     -   tmp109*(hb**2*sbe**2*tmp37 + cbe**2*ht**2*tmp48 + 
     -      cbe*hb*ht*sbe*(0.25D0*(mb*mt)/(c2b*c2t) + 
     -         0.125D0*(b + t)/(s2b*s2t) + 0.5D0*(mb*tmp39)/s2t + 
     -         0.5D0*(mt*tmp44)/s2b + 0.125D0*(Yb*Yt)/(s2b*s2t)))
        Dcptpb = -2*cbe*hb*ht*mb*mt*sbe*tmp108*tmp56 + 
     -   cbe*hb*ht*mb*mt*sbe*
     -    (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -      LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*tmp56
     -     + cbe*hb*ht*mb*mt*sbe*
     -    (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -      LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*tmp56
     -     - 2*cbe*hb*ht*mb*mt*sbe*tmp109*tmp59 + 
     -   cbe*hb*ht*mb*mt*sbe*
     -    (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -      LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*tmp59
     -     + cbe*hb*ht*mb*mt*sbe*
     -    (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -      LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*tmp59
     -     - 2*cbe*hb*ht*mb*mt*sbe*tmp59*tmp88 - 
     -   2*cbe*hb*ht*mb*mt*sbe*tmp56*tmp89 - 
     -   4*cbe*hb*ht*mb*mt*sbe*
     -    (-2*A0*LogA0 - 2*b*Logb - 2*Logt*t - 
     -      0.5D0*(Logb*Logt*(A0 - b - t)) - 
     -      0.5D0*(LogA0*Logt*(-A0 + b - t)) + 
     -      0.5D0*(deltA0bt*phiA0bt)/t - 
     -      0.5D0*(LogA0*Logb*(-A0 - b + t)) + 2.5D0*(A0 + b + t) + 
     -      0.5D0*tmp76)
        Dcpttptb = 0.25D0*(cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*Xb*
     -      Xt) - 0.25D0*(cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*Xb*
     -      Xt) - 0.25D0*(cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*Xb*
     -      Xt) + 0.25D0*(cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*Xb*
     -      Xt) + 0.5D0*(cbe*hb*ht*s2b*s2t*sbe*tmp108*Yb*Yt) - 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*tmp109*Yb*Yt) - 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*tmp88*Yb*Yt) + 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*tmp89*Yb*Yt)
        Dcpbptt = -2*hb*ht*mb*mu*s2t*tmp113 - 
     -   cbe*hb*ht*sbe*tmp89*(mb*s2t*tmp154 - 0.5D0*(b*s2b*s2t)) - 
     -   cbe*hb*ht*sbe*tmp108*
     -    (-(mb*s2t*tmp155) - 0.5D0*(b*s2b*s2t)) - 
     -   cbe*hb*ht*sbe*tmp109*
     -    (-(mb*s2t*tmp154) + 0.5D0*(b*s2b*s2t)) - 
     -   cbe*hb*ht*sbe*tmp88*(mb*s2t*tmp155 + 0.5D0*(b*s2b*s2t)) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (mb*s2t*tmp128 - 0.5D0*(b*s2b*s2t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (-(mb*s2t*tmp129) - 0.5D0*(b*s2b*s2t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(mb*s2t*tmp128) + 0.5D0*(b*s2b*s2t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (mb*s2t*tmp129 + 0.5D0*(b*s2b*s2t)))
        Dcptptb = -2*hb*ht*mt*mu*s2b*tmp75 - 
     -   cbe*hb*ht*sbe*tmp89*
     -    (-(mt*s2b*tmp156) - 0.5D0*(s2b*s2t*t)) - 
     -   cbe*hb*ht*sbe*tmp108*(mt*s2b*tmp157 - 0.5D0*(s2b*s2t*t)) - 
     -   cbe*hb*ht*sbe*tmp88*(mt*s2b*tmp156 + 0.5D0*(s2b*s2t*t)) - 
     -   cbe*hb*ht*sbe*tmp109*
     -    (-(mt*s2b*tmp157) + 0.5D0*(s2b*s2t*t)) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (-(mt*s2b*tmp130) - 0.5D0*(s2b*s2t*t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (mt*s2b*tmp131 - 0.5D0*(s2b*s2t*t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)*
     -      (mt*s2b*tmp130 + 0.5D0*(s2b*s2t*t))) + 
     -   0.5D0*(cbe*hb*ht*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(mt*s2b*tmp131) + 0.5D0*(s2b*s2t*t)))
        Dcptmptt = -(tmp109*
     -      (0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t) + 
     -        hb**2*sbe**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 + c2b)*mt*s2t*Yb)) + 
     -        cbe**2*ht**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 - c2b)*mt*s2t*Yt))))-
     -     tmp89*(-(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t)) + 
     -      hb**2*sbe**2*(-(0.5D0*(mb*mt*s2b*s2t)) + 
     -         0.5D0*((1 + c2b)*mt*s2t*Yb)) + 
     -      cbe**2*ht**2*(-(0.5D0*(mb*mt*s2b*s2t)) + 
     -         0.5D0*((1 - c2b)*mt*s2t*Yt))) - 
     -   tmp108*(-(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t)) + 
     -      hb**2*sbe**2*(-(0.5D0*(mb*mt*s2b*s2t)) - 
     -         0.5D0*((1 - c2b)*mt*s2t*Yb)) + 
     -      cbe**2*ht**2*(-(0.5D0*(mb*mt*s2b*s2t)) - 
     -         0.5D0*((1 + c2b)*mt*s2t*Yt))) - 
     -   0.25D0*(ht**2*(-2*mt*s2t*sbe**2*tmp115*Xt + 
     -        2*mt*s2t*sbe**2*tmp94*Xt - 
     -        4*cbe**2*mt*s2t*tmp114*Yt + 4*cbe**2*mt*s2t*tmp93*Yt)
     -      ) - 0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 
     -        5*T2 + LogT2*(-2*B2*LogB2 + 4*T2) - 
     -        2*(-B2 + T2)*tmp173)*
     -      (-(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t)) + 
     -        cbe**2*hb**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 + c2b)*mt*s2t*Xb)) + 
     -        ht**2*sbe**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 - c2b)*mt*s2t*Xt))))-
     -     0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t) + 
     -        cbe**2*hb**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) + 0.5D0*((1 + c2b)*mt*s2t*Xb))
     -          + ht**2*sbe**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) + 0.5D0*((1 - c2b)*mt*s2t*Xt))
     -        )) - 0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 
     -        5*T2 + LogT2*(-2*B1*LogB1 + 4*T2) - 
     -        2*(-B1 + T2)*tmp172)*
     -      (0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t) + 
     -        cbe**2*hb**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 - c2b)*mt*s2t*Xb))
     -          + ht**2*sbe**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 + c2b)*mt*s2t*Xt))
     -        )) - 0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 
     -        5*T1 + LogT1*(-2*B1*LogB1 + 4*T1) - 
     -        2*(-B1 + T1)*tmp169)*
     -      (-(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t)) + 
     -        cbe**2*hb**2*
     -         (0.5D0*(mb*mt*s2b*s2t) + 0.5D0*((1 - c2b)*mt*s2t*Xb)) + 
     -        ht**2*sbe**2*
     -         (0.5D0*(mb*mt*s2b*s2t) + 0.5D0*((1 + c2b)*mt*s2t*Xt))))
        Dcptmptt = Dcptmptt - 
     -   tmp88*(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t) + 
     -      hb**2*sbe**2*(0.5D0*(mb*mt*s2b*s2t) + 
     -         0.5D0*((1 - c2b)*mt*s2t*Yb)) + 
     -      cbe**2*ht**2*(0.5D0*(mb*mt*s2b*s2t) + 
     -         0.5D0*((1 + c2b)*mt*s2t*Yt)))
        Dcpbmptb = -(tmp89*
     -      (-(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe)) + 
     -        hb**2*sbe**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 + c2t)*mb*s2b*Yb))
     -          + cbe**2*ht**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 - c2t)*mb*s2b*Yt))
     -        )) - 0.25D0*(hb**2*
     -      (2*cbe**2*(-10*B1 + 4*B1*LogB1 + 
     -           LogB1*(4*B1 - 2*B1*LogB1))*mb*s2b*Xb - 
     -        2*cbe**2*(-10*B2 + 4*B2*LogB2 + 
     -           LogB2*(4*B2 - 2*B2*LogB2))*mb*s2b*Xb + 
     -        4*mb*s2b*sbe**2*Yb*
     -         (2*A0*LogA0 + 4*B1*LogB1 - A0*LogA0*LogB1 - 
     -           2.5D0*(A0 + 2*B1) + 0.5D0*((A0 - 2*B1)*LogB1**2) - 
     -           0.5D0*(deltA0B1B1*phiA0B1B1)/B1) - 
     -        4*mb*s2b*sbe**2*Yb*
     -         (2*A0*LogA0 + 4*B2*LogB2 - A0*LogA0*LogB2 - 
     -           2.5D0*(A0 + 2*B2) + 0.5D0*((A0 - 2*B2)*LogB2**2) - 
     -           0.5D0*(deltA0B2B2*phiA0B2B2)/B2))) - 
     -   0.5D0*((-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)*
     -      (0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe) + 
     -        cbe**2*hb**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 + c2t)*mb*s2b*Xb))
     -          + ht**2*sbe**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) - 0.5D0*((1 - c2t)*mb*s2b*Xt))
     -        )) - 0.5D0*((-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 
     -        5*T1 + LogT1*(-2*B1*LogB1 + 4*T1) - 
     -        2*(-B1 + T1)*tmp169)*
     -      (-(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe)) + 
     -        cbe**2*hb**2*
     -         (0.5D0*(mb*mt*s2b*s2t) + 0.5D0*((1 + c2t)*mb*s2b*Xb)) + 
     -        ht**2*sbe**2*
     -         (0.5D0*(mb*mt*s2b*s2t) + 0.5D0*((1 - c2t)*mb*s2b*Xt))))-
     -     0.5D0*((-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)*
     -      (-(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe)) + 
     -        cbe**2*hb**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 - c2t)*mb*s2b*Xb)) + 
     -        ht**2*sbe**2*
     -         (0.5D0*(mb*mt*s2b*s2t) - 0.5D0*((1 + c2t)*mb*s2b*Xt))))-
     -     0.5D0*((-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)*
     -      (0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe) + 
     -        cbe**2*hb**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) + 0.5D0*((1 - c2t)*mb*s2b*Xb))
     -          + ht**2*sbe**2*
     -         (-(0.5D0*(mb*mt*s2b*s2t)) + 0.5D0*((1 + c2t)*mb*s2b*Xt))
     -        ))
        Dcpbmptb = Dcpbmptb - 
     -   tmp88*(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe) + 
     -      hb**2*sbe**2*(0.5D0*(mb*mt*s2b*s2t) + 
     -         0.5D0*((1 + c2t)*mb*s2b*Yb)) + 
     -      cbe**2*ht**2*(0.5D0*(mb*mt*s2b*s2t) + 
     -         0.5D0*((1 - c2t)*mb*s2b*Yt))) - 
     -   tmp109*(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe) + 
     -      hb**2*sbe**2*(0.5D0*(mb*mt*s2b*s2t) - 
     -         0.5D0*((1 - c2t)*mb*s2b*Yb)) + 
     -      cbe**2*ht**2*(0.5D0*(mb*mt*s2b*s2t) - 
     -         0.5D0*((1 + c2t)*mb*s2b*Yt))) - 
     -   tmp108*(-(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe)) + 
     -      hb**2*sbe**2*(-(0.5D0*(mb*mt*s2b*s2t)) + 
     -         0.5D0*((1 - c2t)*mb*s2b*Yb)) + 
     -      cbe**2*ht**2*(-(0.5D0*(mb*mt*s2b*s2t)) + 
     -         0.5D0*((1 + c2t)*mb*s2b*Yt)))
        Dspbmptbspbptt = 
     -  -(0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp108)) + 
     -   0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp109) - 
     -   0.25D0*(b*cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)) + 
     -   0.25D0*(b*cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)) + 
     -   0.25D0*(b*cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)) - 
     -   0.25D0*(b*cbe*hb*ht*s2b*s2t*sbe*
     -      (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)) + 
     -   0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp88) - 
     -   0.5D0*(b*cbe*hb*ht*s2b*s2t*sbe*tmp89)
        Dsptmpttsptptb = 
     -  -(0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp108)) + 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp109) - 
     -   0.25D0*(cbe*hb*ht*s2b*s2t*sbe*t*
     -      (-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)) + 
     -   0.25D0*(cbe*hb*ht*s2b*s2t*sbe*t*
     -      (-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 5*T1 + 
     -        LogT1*(-2*B2*LogB2 + 4*T1) - 2*(-B2 + T1)*tmp170)) + 
     -   0.25D0*(cbe*hb*ht*s2b*s2t*sbe*t*
     -      (-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 5*T2 + 
     -        LogT2*(-2*B1*LogB1 + 4*T2) - 2*(-B1 + T2)*tmp172)) - 
     -   0.25D0*(cbe*hb*ht*s2b*s2t*sbe*t*
     -      (-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 5*T2 + 
     -        LogT2*(-2*B2*LogB2 + 4*T2) - 2*(-B2 + T2)*tmp173)) + 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp88) - 
     -   0.5D0*(cbe*hb*ht*s2b*s2t*sbe*t*tmp89)
        Dsptmpttspbmptb = 
     -  -(tmp108*tmp11) - tmp109*tmp12 - tmp12*tmp88 - 
     -   tmp11*tmp89 - 0.5D0*
     -    (tmp14*(-5*B1 + 4*B1*LogB1 + LogT1**2*(B1 - T1) - 5*T1 + 
     -        LogT1*(-2*B1*LogB1 + 4*T1) - 2*(-B1 + T1)*tmp169)) - 
     -   0.5D0*(tmp13*(-5*B2 + 4*B2*LogB2 + LogT1**2*(B2 - T1) - 
     -        5*T1 + LogT1*(-2*B2*LogB2 + 4*T1) - 
     -        2*(-B2 + T1)*tmp170)) - 
     -   0.5D0*(tmp13*(-5*B1 + 4*B1*LogB1 + LogT2**2*(B1 - T2) - 
     -        5*T2 + LogT2*(-2*B1*LogB1 + 4*T2) - 
     -        2*(-B1 + T2)*tmp172)) - 
     -   0.5D0*(tmp14*(-5*B2 + 4*B2*LogB2 + LogT2**2*(B2 - T2) - 
     -        5*T2 + LogT2*(-2*B2*LogB2 + 4*T2) - 
     -        2*(-B2 + T2)*tmp173))

      return
      end

*
***********************************************************************
***********************************************************************
*
*      
*     SOME AUXILIARY FUNCTIONS CALLED BY THE DSZ, BDSZ & DDS ROUTINES
*
*     Last updates:
*
*                 24/11/2003: 1) phi expansion around lambda=0 removed
*
*                 06/08/2003: 1) phi expansion around x/z=0,y/z=0
*
*                 13/05/2003: 1) new routine for the complex dilogarithm
*                             2) phi expansion around lambda=0
*                             3) more Passarino-Veltman functions
*
***********************************************************************
*

      real*8 function myAA(m,q)      
      real*8 m,q

      if(m.ne.0d0) then
         myAA = m*(1d0-dLog(m/q))
      else
         myAA = 0d0
      endif

      return
      end

*
***********************************************************************
*


      real*8 function myB0(q,m1,m2,mu2) 

c     from Degrassi and Sirlin, Phys. Rev. D46 (1992) 3104.
      
      real*8 q,m1,m2,Omega,mu2

      if(q.eq.0d0) then

         if(m1.eq.0d0.and.m2.ne.0d0) then
            myB0 = 1d0-dLog(m2/mu2)
         elseif(m1.ne.0d0.and.m2.eq.0d0) then
            myB0 = 1d0-dLog(m1/mu2)
         elseif(m1.eq.m2) then
            myB0 = -dLog(m1/mu2)
         else
            myB0 = 1d0 - dLog(m2/mu2) + m1/(m1-m2)*dLog(m2/m1)
         endif
         
      else

         if(m1.eq.0d0.and.m2.ne.0d0) then
            
            if(m2.ne.q) then
          myB0 = -(dLog(m2/mu2)-2d0-(m2/q-1d0)*dLog(dabs(1d0-q/m2))) 
            else 
               myB0 = -(dLog(m2/mu2) - 2d0)
            endif
            
         elseif(m2.eq.0d0.and.m1.ne.0d0) then
            
            if(m1.ne.q) then
         myB0 = -(dLog(m1/mu2)-2d0-(m1/q-1d0)*dLog(dabs(1d0-q/m1))) 
            else
               myB0 = -(dLog(m1/mu2) - 2d0)
            endif
            
         elseif(m2.eq.0d0.and.m1.eq.0d0) then
            
            myB0 = -(dLog(q/mu2) - 2d0) ! cut the imaginary part (I Pi)
            
         else
            
            myB0 = -( dlog(q/mu2)-2.d0 + 
     1           1.d0/2.d0*( 1.d0 + (m1/q-m2/q))*dlog(m1/q) +
     2           1.d0/2.d0*( 1.d0 - (m1/q-m2/q))*dlog(m2/q) +
     3           2.d0*Omega(m1/q,m2/q))
            
         endif
         
      endif

      return
      end
      
c     function Omega(a,b) contained in myB0
      real*8 function Omega(a,b)
      real*8 a,b,cbig
      Cbig = (a+b)/2.d0 - (a-b)**2.d0/4.d0 -1.d0/4.d0
      if(Cbig.gt.0.d0) then
         Omega = dsqrt(Cbig)*
     1        (datan((1.d0 + a - b)/(2.d0*dsqrt(Cbig))) +
     2        datan((1.d0 - a + b)/(2.d0*dsqrt(Cbig))) )
      elseif(Cbig.lt.0d0) then
         Cbig = - Cbig
         Omega = 1.d0/2.d0*dsqrt(Cbig)*
     1        dlog((a/2.d0 +b/2.d0 -1.d0/2.d0 -dsqrt(Cbig))/
     2        (a/2.d0 + b/2.d0 -1.d0/2.d0 + dsqrt(Cbig)))
      else
         Omega = 0         
      endif

      return
      end

*
**********************************************************************
*
      
      real*8 function myB1(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0
      
      if(p.eq.0d0) then
         myB1 = (1d0-Log(m2/q)+m1**2/(m1-m2)**2*Log(m2/m1)
     $        +(m1+m2)/(m1-m2)/2d0)/2d0
      else
         myB1 = (myAA(m2,q)-myAA(m1,q)+(p+m1-m2)*myB0(p,m1,m2,q))/2d0/p
      endif
      
      return
      end

*
**********************************************************************
*
      
      real*8 function myG(p,m1,m2,q)

      implicit none

      real*8 p,m1,m2,q
      real*8 myAA,myB0

      myG = (p-m1-m2)*myB0(p,m1,m2,q) - myAA(m1,q) - myAA(m2,q)

      return
      end
      
*
**********************************************************************
*

      function phi(x,y,z)

c     from Davydychev and Tausk, Nucl. Phys. B397 (1993) 23

      implicit none
      real*8 x,y,z,phi,pphi,myphi
      
      if(x.le.z.and.y.le.z) then
         pphi = myphi(x,y,z)
      elseif(z.le.x.and.y.le.x) then
         pphi = z/x*myphi(z,y,x)
      elseif(z.le.y.and.x.le.y) then
         pphi = z/y*myphi(z,x,y)
      endif

      phi = pphi
      end
      
      function myphi(x,y,z)
      
      implicit none

      real*8 x,y,z,myphi
      real*8 u,v
      real*8 Pi,Li2
      complex*16 clam,cxp,cxm,CLI2,ccphi
        logical su_isnan
        real*8 flip,flop
      Pi = 3.14159265358979d0

c     auxiliary variables

      u = x/z
      v = y/z
      
      if(u.le.1d-8) then
         
         if(v.ne.1d0) then
            myphi = (log(u)*log(v)+2d0*Li2(1d0-v))/(1d0-v)
         else
            myphi = 2d0-log(u)
         endif

      elseif(v.le.1d-8) then

         if(u.ne.1d0) then
            myphi = (log(v)*log(u)+2d0*Li2(1d0-u))/(1d0-u)
         else
            myphi = 2d0-log(v)
         endif

      else
         
         if((1d0-u-v)**2.ge.4d0*u*v) then         
            clam = dcmplx(dsqrt((1d0-u-v)**2 - 4d0*u*v),0d0)
         else
            clam = dcmplx(0d0,dsqrt(4d0*u*v - (1d0-u-v)**2))
         endif
         
         cxp = (1d0+(u-v)-clam)/2d0
         cxm = (1d0-(u-v)-clam)/2d0
         
c     phi function from eq. (A4)
            
         ccphi = (2d0*log(cxp)*log(cxm) - log(u)*log(v) - 
     &        2d0*(CLI2(cxp) + CLI2(cxm)) + Pi**2/3d0)/clam
         myphi = dreal(ccphi)
                     
      endif
         if(su_isNaN(myphi)) then            !jlk test 
         flop=flip
         endif
      return
      end


*
***********************************************************************
*

      function delt(x,y,z)

      implicit none

      real*8 x,y,z,delt

      delt = x**2 + y**2 + z**2 - 2d0*(x*y + x*z + y*z)

      return
      end

*
***********************************************************************
*

      function Li2(x)

      implicit none

      complex*16 CLI2,z
      real*8 x,Li2

      z = DCMPLX(x,0d0)
      Li2 = dreal(CLI2(z))

      return
      end

*
***********************************************************************
*

      COMPLEX*16 FUNCTION CLI2(Z)

c     just call the CSPEN routine
      
      COMPLEX*16 Z,CSPEN

      CLI2 = CSPEN(Z)

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        FUNCTION CSPEN(ZZ)

        COMPLEX*16 CSPEN,W,SUM,ZZ,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...

      B(1)=0.1666666666666666666666666667d0
      B(2)=-0.0333333333333333333333333333d0
      B(3)=0.0238095238095238095238095238d0
      B(4)=-0.0333333333333333333333333333d0
      B(5)=0.0757575757575757575757575758d0
      B(6)=-0.2531135531135531135531135531d0
      B(7)=1.1666666666666666666666666667d0
      B(8)=-7.09215686274509804d0
      B(9)=54.97117794486215539d0

c      write(*,*) 'z:',z
      Z =ZZ*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c      write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
        CSPEN=-CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
        RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
        CSPEN=1.64493406684822643D0
c        write(*,*) 'cspen:', cspen
        RETURN
      END IF
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CDLOG(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 2
c      write(*,*) 'u:',u
c      write(*,*) 'sum:',sum
      DO 1 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2
      SUM=SUM+U*B(K)
 1    CONTINUE
 2    CSPEN=SUM
c        write(*,*) 'cspen:', cspen
      RETURN
10    W=-CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 12

      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12
      SUM=SUM+U*B(K)
11    CONTINUE
12    CSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
c        write(*,*) 'cspen:', cspen
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 22
      DO 21 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22
      SUM=SUM+U*B(K)
21    CONTINUE
22    CSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
      RETURN
30    W=CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 32
      DO 31 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32
      SUM=SUM+U*B(K)
31    CONTINUE
32    CSPEN=SUM+3.28986813369645287D0
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cccccccccccccccccc  Now add the last programs by Pietro on the tau's
c
c     Two-loop O(a_tau*a_bot) corrections to the Higgs masses 
c     and to the minimization conditions of the effective potential. 
c     Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on unpublished work of P. Slavich.
c
c     Last update:  19/06/2003.
c
c     I/O PARAMETERS:
c     t = m_tau^2, b = m_bot^2, T1 = m_stau1^2, T2 = m_stau2^2,
c     B1 = m_sbot1^2, B2 = m_sbot2^2, st = sin(theta_stau), 
c     ct = cos(theta_stau), sb = sin(theta_sbot), cb = cos(theta_sbot),
c     q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanb = tan(beta), 
c     vv = v^2,
c     Sij = 2-loop corrections to the CP-even Higgs mass matrix elements,
c     Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles,
c     DMA = 2-loop corrections to the CP-odd Higgs mass.
c
c     Notice: we assume that the 1-loop part is computed in terms of 
c             running (DRbar) parameters, evaluated at the scale Q. The 
c             parameters in the bottom/sbottom (tau/stau) sector should be
c             computed in term of the "resummed" bottom (tau) Yukawa coupling.
c
c

      subroutine SU_taubot(t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,
     $     S11,S12,S22) 

      implicit none

      real*8 t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,S11,S12,S22
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,ht,hb,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp

      pi = 3.14159265897d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu*tanb         ! notice the sign convention for mu
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = dsqrt(2d0/vv)*mt/cbe
      hb = dsqrt(2d0/vv)*mb/cbe

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncstau(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      S11 = 
     $     2d0*hb**2*mb**2*F1b + 2d0*hb**2*Ab*mb*s2b*F2b
     $     + .5d0*hb**2*Ab**2*s2b**2*F3b
     $     + 2d0*ht**2*mt**2*F1t + 2d0*ht**2*At*mt*s2t*F2t
     $     + .5d0*ht**2*At**2*s2t**2*F3t
     $     + 2d0*hb*ht*mb*At*s2t*F4t + 2d0*ht*hb*mt*Ab*s2b*F4b
     $     + ht*hb*At*Ab*s2t*s2b*F5 + 4d0*ht*hb*mt*mb*F6

      S12 = ht**2*mu*mt*s2t*F2t + .5d0*ht**2*At*mu*s2t**2*F3t
     $     +hb**2*mu*mb*s2b*F2b + .5d0*hb**2*Ab*mu*s2b**2*F3b
     $     +ht*hb*mb*mu*s2t*F4t + hb*ht*mt*mu*s2b*F4b
     $     +.5d0*ht*hb*s2t*s2b*mu*(At+Ab)*F5
     $     +2d0*ht*hb*mt*mb*F6

      S22 = .5d0*hb**2*mu**2*s2b**2*F3b+ .5d0*ht**2*mu**2*s2t**2*F3t
     $      + hb*ht*mu**2*s2b*s2t*F5

      S11 = k*S11
      S12 = k*S12
      S22 = k*S22

      return
      end


*
***********************************************************************
*
      
      subroutine SU_taubotodd(t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,
     $     DMA) 

      implicit none

      real*8 t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,DMA
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,ht,hb,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp

      pi = 3.14159265897d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu*tanb           
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = dsqrt(2d0/vv)*mt/cbe
      hb = dsqrt(2d0/vv)*mb/cbe

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncstau(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      DMA = -(ht**2*mu*At/(T1-T2)*Ft + hb**2*mu*Ab/(B1-B2)*Fb
     $     + 2d0*ht*hb*FAp)/sbe/cbe

      DMA = k*DMA

      return
      end


*
***********************************************************************
*
      
      subroutine SU_taubottad(t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv,
     $     S1,S2) 

      implicit none

      real*8 t,b,T1,T2,B1,B2,st,ct,sb,cb,q,mu,tanb,vv
      real*8 c2t,s2t,c2b,s2b,At,Ab,Xt,Xb,mt,mb,cbe,sbe,ht,hb,pi,k
      real*8 F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp
      real*8 v1,v2,S1,S2

      pi = 3.14159265897d0

      mt = dsqrt(t)
      mb = dsqrt(b)

      s2t = 2d0*ct*st
      s2b = 2d0*cb*sb
      c2t = ct**2 - st**2
      c2b = cb**2 - sb**2

      Xt = (T1-T2)*s2t/2d0/mt    
      Xb = (B1-B2)*s2b/2d0/mb    
      At = Xt - mu*tanb           
      Ab = Xb - mu*tanb           

      sbe = dsin(datan(tanb))
      cbe = dcos(datan(tanb))

      ht = dsqrt(2d0/vv)*mt/cbe
      hb = dsqrt(2d0/vv)*mb/cbe

      k = 3d0/(16d0*Pi**2)**2 
 
      call makefuncstau(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp)

      v1 = Sqrt(vv)*cbe
      v2 = Sqrt(vv)*sbe

      S1 = mt*At*s2t*Ft + mb*Ab*s2b*Fb + 2d0*mt**2*Gt + 2d0*mb**2*Gb
      S2 = mt*mu*tanb*s2t*Ft + mb*mu*tanb*s2b*Fb
            
      S1 = k*S1/v1**2
      S2 = k*S2/v2**2

      return
      end

*
***********************************************************************
*
      
      subroutine makefuncstau(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,
     $     q,mu,vv,tanb,F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,
     $     Ft,Fb,Gt,Gb,FAp)

      implicit none

      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb,
     $     F1t,F2t,F3t,F4t,F1b,F2b,F3b,F4b,F5,F6,Ft,Fb,Gt,Gb,FAp
      
      real*8 tauF1q,tauF2q,tauF3q,tauF4q,tauF5q,tauF6q,
     $     tauFq,tauGq,tauFAq

      F1t = tauF1q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      F2t = tauF2q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      F3t = tauF3q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      F4t = tauF4q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)

      F1b = tauF1q(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb) 
      F2b = tauF2q(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb) 
      F3b = tauF3q(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb) 
      F4b = tauF4q(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb) 

      F5 = tauF5q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      F6 = tauF6q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)

      Ft = tauFq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      Gt = tauGq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
     
      Fb = tauFq(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb)
      Gb = tauGq(b,t,B1,B2,T1,T2,s2b,c2b,s2t,c2t,q,mu,vv,tanb)
      
      FAp = tauFAq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)

      return
      end

*
***********************************************************************
*
      
      function tauF1q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)

      implicit none
      
      real*8 tauF1q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))

      tauF1q = -.5d0*ht*hb*s2t*s2b*(T1-T2)/T1/T2*pippob

      return
      end

*
***********************************************************************
*
      
      function tauF2q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
  
      implicit none
   
      real*8 tauF2q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))

      tauF2q = .5d0*ht*hb*s2b/s2t/T1/T2/(T1-T2)*pippob*
     $     (s2t**2*(T1**2-T2**2) + 2d0*c2t**2*T1*T2*Log(T1/T2))
     
      return
      end

*
***********************************************************************
*
      
      function tauF3q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauF3q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
            
      tauF3q = -.5d0*ht*hb*s2b/s2t/T1/T2/(T1-T2)**2*pippob*
     $     ((T1-T2)*(s2t**2*(T1+T2)**2 - 8d0*c2t**2*T1*T2)
     $     -2d0*(3d0*s2t**2 - 2d0)*T1*T2*(T1+T2)*Log(T1/T2))

      return
      end

*
***********************************************************************
*
      
      function tauF4q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauF4q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
                       
      tauF4q = .5d0*ht*hb*s2b/s2t/(T1-T2)*Log(B1/B2)*
     $     (2d0*c2t**2*pippot + s2t**2*(T1-T2)*Log(T1*T2/q**2))
      
      return
      end

*
***********************************************************************
*
      
      function tauF5q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauF5q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
            
      tauF5q = .5d0*hb*ht/s2b/s2t* 
     $     (-4d0*pippob*pippot/(B1-B2)/(T1-T2)*(1d0-c2b**2*c2t**2)
     $     + 2d0*pippot*s2b**2*c2t**2/(T1-T2)*Log(B1*B2/q**2)
     $     + 2d0*pippob*s2t**2*c2b**2/(B1-B2)*Log(T1*T2/q**2)
     $     + s2t**2*s2b**2*(Log(B1/q)*Log(T1/q) + Log(B1/q)*Log(T2/q)
     $     + Log(B2/q)*Log(T1/q) + Log(B2/q)*Log(T2/q)))
      
      return
      end

*
***********************************************************************
*
      
      function tauF6q(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauF6q
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
            
      tauF6q = .5d0*ht*hb*s2t*s2b*Log(B1/B2)*Log(T1/T2)
      
      return
      end

*
***********************************************************************
*
      
      function tauFq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauFq
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
            
      tauFq = .5d0*ht*hb*s2b/s2t/(T1-T2)*pippob*
     $     (2d0*c2t**2*pippot + s2t**2*(T1-T2)*Log(T1*T2/q**2))
      
      return
      end

*
***********************************************************************
*
      
      function tauGq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauGq
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))
            
      tauGq = .5d0*hb*ht*s2b*s2t*pippob*Log(T1/T2)
      
      return
      end

*
***********************************************************************
*
      
      function tauFAq(t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb)
      
      implicit none

      real*8 tauFAq
      real*8 t,b,T1,T2,B1,B2,s2t,c2t,s2b,c2b,q,mu,vv,tanb
      real*8 Xt,Xb,At,Ab
      real*8 ht,hb,pippob,pippot

      pippob = B1*(Log(B1/q)-1d0) - B2*(Log(B2/q)-1d0)
      pippot = T1*(Log(T1/q)-1d0) - T2*(Log(T2/q)-1d0)
      
      ht = dsqrt(2d0*t/vv)/cos(atan(tanb))
      hb = dsqrt(2d0*b/vv)/cos(atan(tanb))

      Xt = (T1-T2)*s2t/2d0/sqrt(t)    
      Xb = (B1-B2)*s2b/2d0/sqrt(b)

      At = Xt - mu*tanb
      Ab = Xb - mu*tanb
      
      tauFAq = 2d0*ht*hb*Sqrt(t)*Sqrt(b)*tanb*(Ab-At)**2*
     $     mu**2/s2b/s2t/(T1-T2)**2/(B1-B2)**2*pippot*pippob
      
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SU_tausqHiggs(t,A0,BL,T1,T2,st,ct,q,mu,tb,v2,OS,
     $     S11,S22,S12)

c     Two-loop O(a_tau^2) corrections to the CP-even Higgs mass matrix. 
c     Routine written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on A. Brignole, G. Degrassi, P. Slavich and F. Zwirner, 
c     hep-ph/0112177 with appropriate substitutions for top -> tau.
c
c     Last update:  19/06/2003.
c
c
c     I/O PARAMETERS:
c     t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
c     st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tb = tan(beta), v2 = v^2, 
c     OS = renormalization scheme for 1-loop (0 = DRbar, 1 = On-Shell),
c     Sij = 2-loop corrections to the CP-even Higgs mass matrix elements.

      implicit none

      integer OS
      real*8 ht,k,mt,pi,v2,tb
      real*8 t,mu2,A0,BL,T1,T2,st,ct,q,A,X,mu,tanb,sb,cb,s2t,c2t
      real*8 F1,F2,F3,dmuF2,dmuF3,dAtF2,dAtF3,DM12,DM22
      real*8 DF1,DF2,DF3,DdmuF2,DdmuF3,DdAtF2,DdAtF3,F2_s
      real*8 S11,S22,S12,osdr,DMom,ShiftB,ShiftB2,ShiftB3,sw

      pi = 3.14159265897d0

      tanb = 1d0/tb           ! swap H1 <-> H2 !!!

      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt    
      A = X - mu/tanb           

      sb = dsin(datan(tanb))
      cb = dcos(datan(tanb))
      
      ht = dsqrt(2d0/v2)*mt/sb

c$$$      k = 3d0*ht**2/(16d0*Pi**2)**2 
      k = ht**2/(16d0*Pi**2)**2 ! remove color factor !!!
      
      call taufuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,F1,F2,F3)
      call tausfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,A,ht,
     $     dmuF2,dmuF3,dAtF2,dAtF3,DM12,DM22)
      call taudfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,A,v2,
     $     DF1,DF2,DF3,DdmuF2,DdmuF3,DdAtF2,DdAtF3)

      osdr = 1d0*OS
 
      if(s2t.ne.0d0.and.A.ne.0d0) then

         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + 2d0*dmuF3 +
     $        osdr*(DF3 + 2d0*DdmuF3))
         
         S12 = .5d0 * ht**2 * mu * A  * s2t**2 * (F3 + dmuF3 + dAtF3 +
     $        osdr*(DF3 + DdmuF3 + DdAtF3)) + 
     $        ht**2 * mt * mu * s2t * (F2 + dmuF2 +
     $        osdr*(DF2 + DdmuF2))
         
         S22 = .5d0 * ht**2 * A**2 * s2t**2 * (F3 + 2d0*dAtF3 + 
     $        osdr*(DF3 + 2d0*DdAtF3)) + 
     $        2d0 * ht**2 * mt * A * s2t * (F2 + dAtF2 + 
     $        osdr*(DF2 + DdAtF2)) + 
     $        2d0 * ht**2 * mt**2 * (F1 + osdr*DF1)
         
c     some of the functions have poles in s2t=0 or in A=0. 
c     when necessary we consider the residues:
         
      elseif(s2t.eq.0d0.and.A.eq.0d0) then
         
         S11 = 0d0
         S12 = 0d0
         S22 = 2 * ht**2 * mt**2 * (F1 + osdr*DF1)

      elseif(s2t.eq.0d0.and.A.ne.0d0) then 

         call tauresfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,F2_s)

         S11 = 0d0
         S12 = ht**2 * mt * mu * (F2_s + osdr*DF2)
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * A * (F2_s + osdr*(DF2 + DdAtF2))

      elseif(s2t.ne.0d0.and.A.eq.0d0) then

         S11 = .5d0 * ht**2 * mu**2 * s2t**2 * (F3 + 2d0*dmuF3 +
     $        osdr*(DF3 + 2d0*DdmuF3))
         S12 = .5d0 * ht**2 * mu * s2t**2 * osdr*DdAtF3 +
     $        ht**2 * mt * mu * s2t * (F2 + dmuF2 +
     $        osdr*(DF2 + DdmuF2))
         S22 = 2d0 * ht**2 * mt**2 * (F1 + osdr*DF1) +
     $        2d0 * ht**2 * mt * s2t * osdr*DdAtF2
 
      endif
         
      S11 = k*S11
      S12 = k*(S12 + DM12)
      S22 = k*(S22 + DM22)

      sw = S11                  ! swap H1 <-> H2 !!!
      S11 = S22
      S22 = sw

      return
      end

*
***********************************************************************
*

      subroutine taufuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,F1,F2,F3)

      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,F1,F2,F3
      real*8 tauF1ab,tauF1c,tauF2ab,tauF2c,tauF3ab,tauF3c

      F1 = tauF1ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu) 
     $     + tauF1c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
     $     + tauF1c(t,A0,BL,T2,T1,-s2t,-c2t,cb,sb,q,mu)
      
      F2 = tauF2ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu) 
     $     + tauF2c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
     $     - tauF2c(t,A0,BL,T2,T1,-s2t,-c2t,cb,sb,q,mu)

      F3 = tauF3ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu) 
     $     + tauF3c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
     $     + tauF3c(t,A0,BL,T2,T1,-s2t,-c2t,cb,sb,q,mu)

      return
      end

*
*********************************************************************
*
            
      function tauF1ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi
      real*8 tauF1ab

      mu2 = mu**2
      if(mu2.eq.0d0) mu2 = 1d-10

      tauF1ab = 
     $     (2.*BL*mu2/t*(mu2-BL+t)/delt(BL,mu2,t)
     $     +(BL+mu2-t)/t)*phi(BL,mu2,t)
     $     +2.*A0*cb**2*(A0-6.*t)/(A0-4.*t)/t*phi(A0,t,t)
     $     -2.*cb**2*Li2(1.-A0/t)
     $     -2. -(2.+sb**2)/3.*Pi**2
     $     + Log(t/q)*(
     $     (4.*(BL-mu2-10.*t)*t+A0*
     $     (mu2-BL+(6.+4.*sb**2)*t))/(A0-4.*t)/t
     $     +1./delt(BL,mu2,t)*((BL-mu2)**3/t
     $     +(2.*mu2**2+2.*BL*mu2+5.*BL*t+mu2*t-4.*BL**2-2.*t**2)))
     $     +Log(A0/q)*(4.*A0*cb**2/(A0-4.*t))
     $     +Log(BL/q)*(-BL/t+BL*(-BL+mu2+t)**2/t/delt(BL,mu2,t))
     $     +Log(mu2/q)*(mu2/t+mu2*(t**2-(BL-mu2)**2)/t/delt(BL,mu2,t))
     $     +Nc*(Log(t/q)**2-Log(T1/q)**2/2.-Log(T2/q)**2/2.)
     $     +Log(BL/q)*Log(mu2/q)-Log(BL/q)*Log(t/q)-Log(mu2/q)*Log(t/q)
     $     -3.*Log(t/q)**2-2.*Log(T1/q)**2-2.*Log(T2/q)**2

      return
      end
      
*
*********************************************************************
*

      function tauF1c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi,Xt,Yt,st2,ct2
      real*8 tauF1c

      mu2 = mu**2
      if(mu2.eq.0d0) mu2 = 1d-10

      Xt  = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt  = s2t*(T1-T2)/2d0/Sqrt(t) - mu/sb/cb
      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0
 
      tauF1c =
     $     2.*mu2**2*(mu2+t-T1)/T1/delt(T1,mu2,t)*phi(mu2,t,T1)
     $     -A0**2*(1+c2t**2)*cb**2*Yt**2
     $     /2./T2/delt(A0,T1,T2)*phi(A0,T1,T2)
     $     -cb**2*A0/T1*((2.*Sqrt(t)+s2t*Yt)/Sqrt(t)
     $     +(2.*Sqrt(t)+s2t*Yt)**2/2./(A0-4.*T1))*phi(A0,T1,T1)
     $     -cb**2/T1*phi(A0,BL,T1)*(
     $     2.*A0*BL*(ct2*t+s2t*Yt*Sqrt(t)+st2*Yt**2)/delt(A0,T1,BL)
     $     +(A0+BL-T1)*((1.+c2t)*Sqrt(t)+s2t*Yt)/2./Sqrt(t))
     $     +sb**2*(1+c2t+s2t*Xt/Sqrt(t))*Li2(1-BL/T1)
     $     +(1-c2t)*Li2(1-mu2/T1)
     $     +s2t**2*(T1-T2)**2/4./T1/T2*Nc
     $     +(1./3.+s2t*(sb**2*(Xt-Yt)+Yt)/4./Sqrt(t))*Pi**2
     $     -3.*s2t*(sb**2*Xt+cb**2*Yt)*(2.*t-T1)/2./Sqrt(t)/T1
     $     -(sb**2*Xt**2+cb**2*Yt**2)/4./T1/T2*
     $     ((1.+c2t**2)*T1+(5.-2.*c2t-c2t**2)*T2)
     $     +(3.-c2t)*mu2/T1-(3.-c2t)*cb**2*A0/2./T1
     $     -(1.+c2t)*t/2./T1-(1.-c2t)*BL/2./T1
     $     -(1.+c2t**2)*(T1**2+T2**2)/4./T1/T2
     $     +5./2.+c2t/2.-s2t**2/2.
     $     +Log(t/q)*(1.+mu2/t-t/T1-T1/t
     $     +(-(mu2-T1)**3/t + 4.*mu2**2 + 5.*mu2*t + 2.*t**2
     $     +mu2**2*t/T1-t**3/T1-2.*mu2*T1-2.*T1**2)/delt(T1,mu2,t))
     $     +Log(mu2/q)*(mu2*((-2.+c2t)*t+T1)/t/T1
     $     -mu2*((T1-t)**3+2.*mu2*(t-T1)*T1+mu2**2*(t+T1))
     $     /t/T1/delt(mu2,t,T1))
     $     +Log(BL/q)*((1-c2t)/2.*BL/T1
     $     -cb**2*BL*(A0-BL+T1)/T1/delt(A0,T1,BL)*
     $     (ct2*t+s2t*Yt*Sqrt(t)+st2*Yt**2)
     $     +sb**2*BL*(ct2*t+s2t*Xt*Sqrt(t)+st2*Xt**2)/(BL-T1)/T1)
     $     +Log(A0/q)*((3.-c2t)*A0*cb**2/2./T1
     $     +A0*cb**2*(2.*Sqrt(t)+s2t*Yt)**2/2./(A0-4.*T1)/T1
     $     +A0*(1+c2t**2)*cb**2*(A0*(T1+T2)-(T1-T2)**2)*Yt**2
     $     /4./T1/T2/delt(A0,T1,T2)
     $     +A0*cb**2*(A0-BL-T1)/T1/delt(A0,BL,T1)*
     $     (ct2*t+s2t*Yt*Sqrt(t)+st2*Yt**2))
     $     +Log(T2/q)*((1+c2t**2)*T2/4./T1-Nc*s2t**2*T2/4./T1
     $     -cb**2*(1+c2t**2)*Yt**2/4./T2
     $     +sb**2*(1+c2t**2)*Xt**2/4./T1
     $     +cb**2*(1+c2t**2)*Yt**2/4./T1/T2/delt(A0,T1,T2)*
     $     (A0**2*T1-A0*(2.*T1**2+5*T1*T2+T2**2)+(T1-T2)**2*(T1+T2)))
     $     +Log(T1/q)*(cb**2*(ct2*t+s2t*Yt*Sqrt(t)+st2*Yt**2)*
     $     (1./T1-(A0+BL-T1)/delt(A0,T1,BL))
     $     +cb**2*(1+c2t**2)*Yt**2/4./T1/T2/delt(A0,T1,T2)*
     $     (A0**2*T2-A0*(T1**2+5*T1*T2+2.*T2**2)+(T1-T2)**2*(T1+T2))
     $     +1./delt(T1,mu2,t)*((mu2-T1)**2*T1/t-(mu2-t)**2*(mu2+t)/T1
     $     +6.*mu2**2+4.*mu2*t+2.*t**2-3.*mu2*T1-2.*T1**2)
     $     +cb**2*(A0-8.*T1)/2./T1/(A0-4.*T1)*(2*Sqrt(t)+s2t*Yt)**2
     $     +sb**2*(BL-2.*T1)/T1/(BL-T1)*
     $     (ct2*t+s2t*Xt*Sqrt(t)+st2*Xt**2)
     $     -(1-c2t)*(mu2-2.*T1)/2./T1-s2t**2*Nc*(T1-2.*T2)/4./T2
     $     +sb**2*c2t*(1-c2t)*Xt**2/2./T1
     $     +sb**2*(1+c2t**2)*Xt**2/4./T2
     $     +sb**2*s2t*(t-3.*T1)/Sqrt(t)/T1*Xt
     $     -3.*cb**2*s2t*(t+T1)/Sqrt(t)/T1*Yt
     $     -cb**2*(5.-2.*c2t-c2t**2)*Yt**2/4./T1
     $     -(3.+c2t-8.*sb**2)*t/2./T1
     $     +(3.-c2t)/2./T1*mu2+(1+c2t**2)*T1/4./T2
     $     -T1/t+(-14+c2t-c2t**2)/2.)
     $     +(Nc+1.)*s2t**2/4.*
     $     (3.*Log(T1/q)**2-2.*Log(T1/q)*Log(T2/q)-Log(T2/q)**2)
     $     +s2t*(6.*sb**2*Xt+5*cb**2*Yt)/2./Sqrt(t)*Log(T1/q)**2
     $     +(9.+sb**2-cb**2*c2t)/2.*Log(T1/q)**2
     $     +Log(T1/q)*Log(T2/q)-2.*Log(t/q)*Log(T1/q)
     $     +cb**2*(1.+c2t+s2t*Yt/Sqrt(t))/2.*(Log(A0/q)*Log(T1/q)
     $     +Log(BL/q)*Log(T1/q)-Log(A0/q)*Log(BL/q))

      return
      end

*
*********************************************************************
*

      function tauF2ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
      
      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi
      real*8 tauF2ab
      
      tauF2ab = -(3.+Nc)/2.*(Log(T1/q)**2-Log(T2/q)**2)

      return
      end

*
*********************************************************************
*

      function tauF2c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi,ct2,st2,Xt,Yt,At
      real*8 tauF2c

      mu2 = mu**2
      if(mu2.eq.0d0) mu2 = 1d-10

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb
      At = sb**2*Xt+cb**2*Yt

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0

      tauF2c = 4*mu2**2*t/T1/delt(mu2,t,T1)*phi(mu2,t,T1)
     $ +(A0*c2t**2*Yt**2/(T1-T2)/T2 
     $ +(1+c2t**2)/2.*(T1/T2-1)*Yt**2*A0/delt(A0,T1,T2))
     $ *cb**2*phi(A0,T1,T2)
     $ -(2*A0*c2t**2*Sqrt(t)*Yt/s2t/T1/(T1-T2)
     $ +s2t*(A0*Yt/2./Sqrt(t)/T1+Yt*A0*(A0-4*T1)/2./Sqrt(t)/T1/(T1-T2))
     $ + A0*(2*Sqrt(t)+s2t*Yt)**2/2./T1/(A0-4*T1)
     $ + A0/T1+A0*c2t**2*Yt**2/T1/(T1-T2))*cb**2*phi(A0,T1,T1)
     $ -(2*A0*BL/T1*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)/delt(A0,BL,T1)
     $ +c2t**2*(A0+BL-T1)*Yt*Sqrt(t)/s2t/T1/(T1-T2)
     $ +s2t*((A0+BL-T1)*Yt/4./Sqrt(t)/T1
     $ +delt(A0,BL,T1)*Yt/2./Sqrt(t)/T1/(T1-T2))
     $ +ct2*(A0+BL-T1)/2./T1+c2t*Yt**2*(A0+BL-T1)/2./T1/(T1-T2)
     $ +c2t*delt(A0,BL,T1)/2./T1/(T1-T2)
     $ -c2t*t*(A0+BL-T1)/2./T1/(T1-T2))*cb**2*phi(A0,BL,T1)
     $ +(s2t*(BL-T1)*Xt/Sqrt(t)/(T1-T2)+s2t*Xt/2./Sqrt(t)+ct2
     $ -c2t*(t+T1-BL)/(T1-T2)+c2t*Xt*(2*c2t*Sqrt(t)+s2t*Xt)
     $ /s2t/(T1-T2))*sb**2*Li2(1-BL/T1)
     $ +(1-c2t-2*c2t*(mu2-T1)/(T1-T2))*Li2(1-mu2/T1)
     $ -21*s2t*T1*At/2./(T1-T2)/Sqrt(t)+c2t*T1/2./(T1-T2)
     $ +3*s2t*At/4./Sqrt(t)-3*s2t*Sqrt(t)*At/T1
     $ +c2t*(2*BL+2*A0*cb**2-4*mu2-2*t+T1
     $ +2*(sb**2*Xt**2+cb**2*Yt**2))/4./T1
     $ +(-5+c2t**2)/4./T1*(sb**2*Xt**2+cb**2*Yt**2)
     $ -(2*BL+6*A0*cb**2-12*mu2+2*t+(1+c2t**2-Nc*s2t**2)*T2)/4./T1
     $ +((1+c2t**2-Nc*s2t**2)*T1
     $ +(1+c2t**2)*(sb**2*Xt**2+cb**2*Yt**2))/4./T2
     $ -t/T1*Log(t/q)+(2*mu2*(t-T1)*T1+mu2**2*(t+T1)-(t-T1)**3)
     $ /T1/delt(T1,mu2,t)*Log(t*T1/q**2)-Log(mu2/q)*(2-c2t)*mu2/T1
     $ +mu2*((t-T1)**2-mu2**2)/T1/delt(T1,mu2,t)*Log(mu2*T1/q**2)
     $ +Log(BL/q)*(st2*BL/T1+sb**2*BL/T1/(BL-T1)*
     $ (ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2))
     $ -cb**2*BL*(A0-BL+T1)*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)
     $ /T1/delt(A0,T1,BL)*Log(BL*T1/q**2)
     $ +Log(A0/q)*((3-c2t)*cb**2*A0/2./T1
     $ +A0*cb**2*(2*Sqrt(t)+s2t*Yt)**2/2./(A0-4*T1)/T1)
     $ -cb**2*A0*(1+c2t**2)*(A0-T1-T2)*(T1-T2)*Yt**2
     $ /4./T1/T2/delt(A0,T1,T2)*Log(A0*T1/q**2)
     $ +cb**2*A0*(A0-BL-T1)*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)
     $ /T1/delt(A0,T1,BL)*Log(A0*T1/q**2)
     $ +Log(T2/q)*(c2t**2*(1+Nc)*(T1+T2)/(T1-T2)
     $ +c2t**2*(sb**2*Xt**2+cb**2*Yt**2)/(T1-T2)
     $ -(1+c2t**2)*sb**2*(T1-T2)*Xt**2/4./T1/T2
     $ +(1+c2t**2-s2t**2*Nc)*T2/4./T1
     $ +(1+c2t**2)*(sb**2*Xt**2+cb**2*Yt**2)/4./T2)
     $ -(1+c2t**2)*cb**2*Yt**2*(A0**2*T1+(T1-T2)**3
     $ -A0*(2*T1**2+3*T1*T2-T2**2))
     $ /4./T1/T2/delt(A0,T1,T2)*Log(T2*T1/q**2)
     $ +Log(T1/q)*(4*mu2*(mu2+t-T1)/delt(T1,mu2,t)
     $ +2*(1+c2t**2)*cb**2*Yt**2*A0*(A0-T1-3*T2)/4./T2/delt(A0,T1,T2)
     $ -2*cb**2*(A0+BL-T1)*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)
     $ /delt(A0,T1,BL)+cb**2*(2*Sqrt(t)+s2t*Yt)**2/2./T1
     $ -2*cb**2*(2*Sqrt(t)+s2t*Yt)**2/(A0-4*T1)
     $ +sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/T1
     $ -sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/(BL-T1)
     $ -st2*mu2/T1+1-c2t-6*c2t**2*Sqrt(t)*At/s2t/(T1-T2)
     $ +9*s2t*T1*At/Sqrt(t)/(T1-T2)
     $ -c2t*(1+c2t)*(sb**2*Xt**2+cb**2*Yt**2)/(T1-T2)
     $ -(1+c2t**2)*sb**2*(T1-T2)*Xt**2/4./T1/T2
     $ -c2t*(BL+A0*cb**2-2*mu2-t+(1+3*c2t*(1+Nc))*T1
     $ -c2t*(1+Nc)*T2)/(T1-T2)
     $ -3*s2t*At/2./Sqrt(t)+s2t*Sqrt(t)*(sb**2*Xt-3*cb**2*Yt)/T1
     $ -(3-2*c2t+c2t**2)/4./T1*(sb**2*Xt**2+cb**2*Yt**2)
     $ +s2t**2/2./T1*(sb**2*Xt**2-cb**2*Yt**2)
     $ -6+c2t+(Nc+1)*s2t**2/2.+((3-c2t)*mu2-(c2t+8*cb**2-5)*t)/2./T1
     $ -(2*c2t**2+(1-Nc)*s2t**2)*T1/4./T2)
     $ +Log(T2/q)**2*(3*c2t**2*Sqrt(t)*At/2./s2t/(T1-T2)
     $ +c2t*(1-2*c2t)*(sb**2*Xt**2+cb**2*Yt**2)/4./(T1-T2)+3/8.
     $ +c2t*(2*BL+2*A0*cb**2-4*mu2-2*t+(1-2*c2t*(1+Nc))*T1
     $ +(1-6*c2t*(1+Nc))*T2)/8./(T1-T2))
     $ +Log(T1/q)**2*(3*s2t/2./Sqrt(t)*At
     $ +s2t*cb**2*A0*Yt/(T1-T2)/Sqrt(t)
     $ +s2t*(BL-6*T1)*At/2./(T1-T2)/Sqrt(t)
     $ +9*c2t**2*Sqrt(t)/2./s2t*At/(T1-T2)
     $ +3*c2t*(1+2*c2t)*(sb**2*Xt**2+cb**2*Yt**2)/4./(T1-T2)
     $ +c2t*(6*BL+6*A0*cb**2-12*mu2-6*t+(7+26*c2t*(1+Nc))*T1
     $ -(1+2*c2t*(1+Nc))*T2)/8./(T1-T2)+25/8.-c2t/2.+s2t**2*(1+Nc))
     $ -(s2t*(2*A0+2*BL-T1-T2)*Yt/4./Sqrt(t)/(T1-T2)
     $ +(1+c2t)/4.+c2t**2*Sqrt(t)*Yt/s2t/(T1-T2)
     $ +c2t*(A0+BL-t-T1+Yt**2)/2./(T1-T2))*cb**2*Log(A0/T1)*Log(BL/T1)
     $ -cb**2*s2t*2*A0*Yt/(T1-T2)/Sqrt(t)*Log(A0/q)*Log(T1/q)
     $ -s2t*BL*At/(T1-T2)/Sqrt(t)*Log(BL/q)*Log(T1/q)
     $ -(c2t**2*(1+Nc)*(T1+T2)/(T1-T2)
     $ +c2t**2*(sb**2*Xt**2+cb**2*Yt**2)/(T1-T2))*Log(T1/q)*Log(T2/q)


      return
      end

*
*********************************************************************
*
      
      function tauF3ab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi
      real*8 tauF3ab

      tauF3ab = (2.+Nc)/2.*(2.-Log(T1/q)-Log(T2/q))
     $     *(2.-(T1+T2)/(T1-T2)*Log(T1/T2))

      return
      end
      
*
*********************************************************************
*

      function tauF3c(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi,ct2,st2,Xt,Yt,At
      real*8 tauF3c

      mu2 = mu**2
      if(mu2.eq.0d0) mu2 = 1d-10

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb
      At = sb**2*Xt+cb**2*Yt

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0

      tauF3c =
     $ (2*mu2*t*(mu2+t-T1)/T1/delt(T1,mu2,t)
     $ -(4*mu2*t+2*delt(T1,mu2,t))/T1/(T1-T2)
     $ -(mu2+t-T1)/T1)*phi(mu2,t,T1)
     $ +(A0*(1+c2t**2)*(A0-2*(T1+T2))*Yt**2/2./T2/delt(A0,T1,T2)
     $ +4*A0*c2t**2*(A0-2*(T1+T2))*Yt**2/T2/(T1-T2)**2
     $ -(1-3*c2t**2)*Yt**2/2./T2)*cb**2*phi(A0,T1,T2)
     $ +(-A0*(2*Sqrt(t)+s2t*Yt)**2/2./T1/(A0-4*T1)
     $ +A0*(2*Sqrt(t)+s2t*Yt)**2/2./T1/(T1-T2)
     $ -2*A0*c2t**2*Yt**2*(2*A0-7*T1-T2)/T1/(T1-T2)**2
     $ +2*A0*(1-3*c2t**2)*Sqrt(t)*(A0-4*T1)*Yt/s2t/T1/(T1-T2)**2
     $ -4*A0*c2t**2*Sqrt(t)*Yt/s2t/T1/(T1-T2))*cb**2*phi(A0,T1,T1)
     $ +(-2*A0*BL*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)/T1/delt(A0,T1,BL)
     $ +2*(1-3*c2t**2)*Sqrt(t)*Yt*delt(A0,T1,BL)/s2t/T1/(T1-T2)**2
     $ -2*c2t**2*Sqrt(t)*Yt*(A0+BL-T1)/s2t/T1/(T1-T2)
     $ +3*c2t*delt(A0,T1,BL)*(t-Yt**2)/T1/(T1-T2)**2
     $ +(A0+BL-T1)*(c2t*(t-Yt**2)+(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2))
     $ /T1/(T1-T2))*cb**2*phi(A0,BL,T1)
     $ -(1-3*c2t**2)*sb**2*Xt**2/(T1-T2)*Li2(1-T2/T1)
     $ +(1-c2t+2*(1-3*c2t)*(mu2-T1)/(T1-T2)
     $ -6*c2t*(mu2-T1)**2/(T1-T2)**2)*Li2(1-mu2/T1)
     $ +(-4*(1-3*c2t**2)*sb**2*Sqrt(t)*(BL-T1)*Xt/s2t/(T1-T2)**2
     $ +4*c2t**2*sb**2*Sqrt(t)*Xt/s2t/(T1-T2)
     $ -2*sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/(T1-T2)
     $ -2*sb**2*c2t*(3*BL-2*T1-T2)*(t-Xt**2)/(T1-T2)**2)*Li2(1-BL/T1)
     $ +2*Sqrt(t)*At*(3*(7-18*c2t**2)*T1-9*c2t**2*T2)/s2t/(T1-T2)**2
     $ -3*s2t*Sqrt(t)*At*(4*T1-T2)/T1/(T1-T2)+2*c2t*At**2/(T1-T2)
     $ -21*c2t*T1*(sb**2*Xt**2+cb**2*Yt**2)/(T1-T2)**2
     $ +5*c2t/2./(T1-T2)*(sb**2*Xt**2+cb**2*Yt**2)
     $ -(6-2*c2t)/4./T1*(sb**2*Xt**2+cb**2*Yt**2)
     $ +3*mu2/T1+c2t*mu2*(T1+T2)/T1/(T1-T2)+18*c2t*mu2*T1/(T1-T2)**2
     $ +3*cb**2*A0*c2t/(T1-T2)-cb**2*A0*(3-c2t)/2./T1
     $ -12*(cb**2*A0+BL)*T1*c2t/(T1-T2)**2
     $ +3*c2t*BL/(T1-T2)-(1-c2t)*BL/2./T1
     $ +15*c2t*T1*t/(T1-T2)**2-3*c2t*t/2./(T1-T2)-(1+c2t)*t/2./T1
     $ -s2t**2/2.*(1+Nc)-c2t**2*(1+Nc)*(T1-T2)**2/4./T1/T2
     $ -9*(1+Nc)*T1*c2t**2/(T1-T2)-3*T1*(2*T1+3*T2)/(T1-T2)**2*c2t
     $ -(1-Nc)*(T1-T2)**2/4./T1/T2-(28+6*Nc)*T1/2./(T1-T2)
     $ +Log(t/q)*(-t/T1)
     $ -Log(t*T1/q**2)*t*((t-T1)**2-mu2**2)/T1/delt(T1,mu2,t)
     $ +Log(mu2/q)*(-6*mu2*(3*mu2-T2)*c2t/(T1-T2)**2-(2-c2t)*mu2/T1)
     $ -Log(mu2*T1/q**2)*mu2*(T1**2+mu2**2-t**2-2*mu2*T1)
     $ /T1/delt(T1,mu2,t)+Log(BL/q)*(12*BL*c2t*T1/(T1-T2)**2
     $ +sb**2*BL*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/T1/(BL-T1)
     $ +(1-c2t)*BL/2./T1-3*c2t*BL/(T1-T2))
     $ -Log(BL*T1/q**2)*cb**2*BL*(A0-BL+T1)/T1/delt(A0,T1,BL)
     $ *(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)
     $ +Log(A0/q)*A0*cb**2*((2*Sqrt(t)+s2t*Yt)**2/2./(A0-4*T1)/T1
     $ +12*c2t*T1/(T1-T2)**2+(3-c2t)/2./T1-3*c2t/(T1-T2))
     $ -Log(A0*T1/q**2)*A0*cb**2*
     $ (Yt**2*(1+c2t**2)/4./T1/T2/delt(A0,T1,T2)
     $ *((T1+T2)**2-A0*(T1+T2)+4*T1*T2)-(A0-BL-T1)/T1/delt(A0,T1,BL)
     $ *(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2))
     $ +Log(T2/q)*(12*c2t**2*T2*Sqrt(t)*At/s2t/(T1-T2)**2
     $ -6*s2t*Sqrt(t)*T2*At/(T1-T2)**2
     $ -3*(1-c2t)*cb**2*A0*T2/(T1-T2)**2-(1-3*c2t)*T2*BL/(T1-T2)**2
     $ +6*(1-c2t)*T2*mu2/(T1-T2)**2-(1+3*c2t)*T2*t/(T1-T2)**2
     $ +(sb**2*Xt**2+cb**2*Yt**2)*(-(2-3*c2t+23*c2t**2)*T2/(T1-T2)**2
     $ -(1+5*c2t**2)/2./(T1-T2)-(1+c2t**2)/4./T2)
     $ +(1+c2t**2)*sb**2*Xt**2*(T1**2-4*T1*T2-T2**2)/4./T1/T2/(T1-T2)
     $ -(1+Nc)*c2t**2/4./(T1-T2)**2*(9*T1**2+32*T1*T2+19*T2**2)
     $ -(1+Nc)*(T1-T2)/4./T1*c2t**2+T2*(2*T1+T2)/(T1-T2)**2*c2t
     $ -T2*((7+3*Nc)*T1+(1-Nc)*T2)/2./(T1-T2)**2
     $ -((1-Nc)*T1-(3-Nc)*T2)/2./(T1-T2)+(1-Nc)*T2/4./T1)
     $ +Log(T2*T1/q**2)*(1+c2t**2)*cb**2*Yt**2/4./T1/T2/delt(A0,T1,T2)
     $ *(A0**2*T1+(T1-T2)**3-2*T2*(T1**2-T2**2)
     $ -A0*(2*T1**2+T1*T2+T2**2))
     $ +Log(T1/q)*(2*(mu2+t-T1)**2/delt(T1,mu2,t)
     $ -cb**2*(1+c2t**2)*Yt**2/2./T2/delt(A0,T1,T2)
     $ *(A0**2-4*T2*(T1-T2)-A0*(T1+3*T2))
     $ -2*cb**2*(A0+BL-T1)*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)
     $ /delt(A0,T1,BL)
     $ +cb**2*(2*Sqrt(t)+s2t*Yt)**2/2./T1
     $ -2*cb**2*(2*Sqrt(t)+s2t*Yt)**2/(A0-4*T1)
     $ +sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/T1
     $ -sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)/(BL-T1)
     $ -(1-c2t)/2.*mu2/T1+1-c2t
     $ +12*Sqrt(t)*At*((-3+7*c2t**2)*T1+c2t**2*T2)/s2t/(T1-T2)**2
     $ +6*Sqrt(t)*(2*T1-T2)*s2t*At/(T1-T2)**2+sb**2*s2t**2*Xt**2/T1
     $ +s2t*Sqrt(t)/T1*(sb**2*Xt-3*cb**2*Yt)
     $ +(1+c2t**2)*sb**2*Xt**2*(T1**2+4*T1*T2-T2**2)/4./T1/T2/(T1-T2)
     $ +(sb**2*Xt**2+cb**2*Yt**2)*((2+15*c2t+23*c2t**2)*T1/(T1-T2)**2
     $ +(5-6*c2t-5*c2t**2)/2./(T1-T2)-(5-2*c2t-c2t**2)/4./T1)
     $ +3*A0*cb**2*(2*(1+c2t)*T1-(1-c2t)*T2)/(T1-T2)**2
     $ +BL*(2*(1+3*c2t)*T1-(1-3*c2t)*T2)/(T1-T2)**2
     $ +t*(2*(1-6*c2t)*T1-(1+3*c2t)*T2)/(T1-T2)**2
     $ +(5-c2t-8*cb**2)*t/2./T1+(3-c2t)*mu2/2./T1
     $ -2*(3+12*c2t)*mu2*T1/(T1-T2)**2-6*(1-c2t)*mu2/(T1-T2)
     $ -9/2.+s2t**2/2.*(1+Nc)+3*c2t/2.+c2t*T1*(4*T1+11*T2)/(T1-T2)**2
     $ +(1+Nc)*(65*T1**2+4*T1*T2-9*T2**2)*c2t**2/4./(T1-T2)**2
     $ +(1+Nc)*(T1-T2)*c2t**2/4./T2+(1-Nc)*T1/4./T2
     $ +T1*((7+3*Nc)*T1+(1-Nc)*T2)/2./(T1-T2)**2
     $ +(5*(5+Nc)*T1+(1-Nc)*T2)/2./(T1-T2))
     $ +(2*Log(mu2/q)**2-Log(mu2/T1)**2)*3*c2t*mu2**2/(T1-T2)**2
     $ +(2*Log(BL/q)**2-Log(BL/T1)**2)*
     $ (2*BL*(1-3*c2t**2)*Sqrt(t)*At/s2t/(T1-T2)**2-BL/2./(T1-T2)
     $ -3*BL*c2t*(T1+T2-2*t+2*sb**2*Xt**2+2*cb**2*Yt**2)/2./(T1-T2)**2)
     $ +(2*Log(A0/q)**2-Log(A0/T1)**2)*cb**2*
     $ (4*A0*(1-3*c2t**2)*Sqrt(t)*Yt/s2t/(T1-T2)**2-3*A0/2./(T1-T2)
     $ -3*A0*c2t*(T1+T2-2*t+2*Yt**2)/2./(T1-T2)**2)
     $ +Log(t/T1)*Log(mu2/T1)*(T1+T2-2*t-2*mu2)/(T1-T2)
     $ +Log(T2/T1)*Log(A0/T1)*Yt**2*cb**2*
     $ (8*A0*c2t**2+(1-3*c2t**2)*(T1-T2))/2./(T1-T2)**2
     $ +Log(BL/T1)*Log(A0/T1)*cb**2*(2*Sqrt(t)*Yt/(T1-T2)**2/s2t
     $ *((A0+BL)*(1-3*c2t**2)-(1-2*c2t**2)*T1+c2t**2*T2)
     $ +s2t*Sqrt(t)*Yt/(T1-T2)+(t+Yt**2)/2./(T1-T2)
     $ +3*c2t*(2*A0+2*BL-T1-T2)*(t-Yt**2)/2./(T1-T2)**2)
     $ +Log(T1/q)*Log(T2/q)*(2*c2t**2*(1+Nc)*(T1+T2)**2/(T1-T2)**2
     $ +((1+5*c2t**2)*T1
     $ -(1-11*c2t**2)*T2)/2./(T1-T2)**2*(sb**2*Xt**2+cb**2*Yt**2))
     $ +Log(T1/q)**2*(9*c2t**2*Sqrt(t)*At/s2t/(T1-T2)
     $ +6*Sqrt(t)*T1*(2-5*c2t**2)*At/s2t/(T1-T2)**2
     $ -3*Sqrt(t)*s2t*(2*T1-T2)*At/(T1-T2)**2
     $ -((5+5*c2t+13*c2t**2)*T1-(3-4*c2t-10*c2t**2)*T2)
     $ /2./(T1-T2)**2*(sb**2*Xt**2+cb**2*Yt**2)
     $ +mu2*((6+5*c2t)*T1-(3-4*c2t)*T2)/(T1-T2)**2
     $ +9/4.-c2t+5/4.*s2t**2*(1+Nc)
     $ +c2t*(t-BL-A0*cb**2)*(5*T1+4*T2)/2./(T1-T2)**2
     $ -(t+BL+3*A0*cb**2)*(2*T1-T2)/2./(T1-T2)**2
     $ +c2t/4.-9*T1*T2*c2t/2./(T1-T2)**2
     $ -c2t**2*(1+Nc)*(8*T1**2+16*T1*T2-T2**2)/2./(T1-T2)**2
     $ -T1*((14+5*Nc)*T1-(10+4*Nc)*T2)/2./(T1-T2)**2)
     $ +Log(T2/q)**2*(-3*c2t**2*Sqrt(t)*(T1+T2)*At/s2t/(T1-T2)**2
     $ +3*s2t*Sqrt(t)*T2*At/(T1-T2)**2
     $ +(2*T2-c2t*(T1+2*T2)+c2t**2*(2*T1+5*T2))
     $ /2./(T1-T2)**2*(sb**2*Xt**2+cb**2*Yt**2)
     $ +T2*(3*cb**2*A0+BL-6*mu2+t)/2./(T1-T2)**2
     $ -c2t*(cb**2*A0+BL-mu2-t)*(T1+2*T2)/2./(T1-T2)**2
     $ +c2t*mu2*(T1+2*T2)/2./(T1-T2)**2-3/4.+s2t**2/4.*(1+Nc)
     $ +c2t**2*((T1+T2)**2+3*T2**2)/2./(T1-T2)**2*(1+Nc)
     $ -c2t*((T1+T2)**2+2*T1*T2)/4./(T1-T2)**2
     $ +T2*(2*(1+Nc)*T1+(2-Nc)*T2)/2./(T1-T2)**2)

      return
      end

*
*********************************************************************
*

      subroutine tauresfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,F2_s)

      implicit none
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 delt,Li2,phi,ct2,st2,Xt,Yt,At
      real*8 F2_s

      mu2 = mu**2

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb
      At = sb**2*Xt+cb**2*Yt

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0
      
      F2_s =
     $     2*c2t**2*sb**2*Sqrt(t)*Xt/(T1-T2)*Li2(1-BL/T1)
     $     -2*c2t**2*sb**2*Sqrt(t)*Xt/(T1-T2)*Li2(1-BL/T2)
     $     -c2t**2*cb**2*Sqrt(t)*Yt/(T1-T2)*Log(A0/T1)*Log(BL/T1)
     $     +c2t**2*cb**2*Sqrt(t)*Yt/(T1-T2)*Log(A0/T2)*Log(BL/T2)
     $     -6*At*c2t**2*Sqrt(t)/(T1-T2)*Log(T1/T2)
     $     +3*c2t**2*At*Sqrt(t)/(T1-T2)*Log(T1/q)**2
     $     -3*c2t**2*At*Sqrt(t)/(T1-T2)*Log(T2/q)**2
     $     -c2t**2*cb**2*Yt*Sqrt(t)*
     $     (A0+BL-T1)/T1/(T1-T2)*phi(A0,BL,T1)
     $     +c2t**2*cb**2*Yt*Sqrt(t)*
     $     (A0+BL-T2)/T2/(T1-T2)*phi(A0,BL,T2)
     $     -2*A0*c2t**2*cb**2*Yt*Sqrt(t)
     $     /T1/(T1-T2)*phi(A0,T1,T1)
     $     +2*A0*c2t**2*cb**2*Yt*Sqrt(t)
     $     /T2/(T1-T2)*phi(A0,T2,T2)

      return
      end

*
*********************************************************************
*

      subroutine tausfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,At,ht,
     $     dmuF2,dmuF3,dAtF2,dAtF3,DM12,DM22)
      
      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,At,ht
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 dmuF2,dmuF3,dAtF2,dAtF3,DM12,DM22


      dmuF2 = -Nc/4.*(Log(T1/q)**2-Log(T2/q)**2)

      dmuF3 = -Nc/4.*(Log(T1/q)+Log(T2/q)-2.)
     $     *(2.-(T1+T2)/(T1-T2)*Log(T1/T2))

      dAtF2 = -(3.+Nc)/2.*(Log(T1/q)**2-Log(T2/q)**2)

      dAtF3 = -(3.+Nc)/2.*(Log(T1/q)+Log(T2/q)-2.)
     $     *(2.-(T1+T2)/(T1-T2)*Log(T1/T2))

      DM12 = ht**2*Nc*s2t*mu*Sqrt(t)/4.
     $     *(Log(T1/q)**2-Log(T2/q)**2)
     $     +ht**2*Nc*s2t**2*mu*At/8.*(Log(T1/q)+Log(T2/q)-2.)
     $     *(2.-(T1+T2)/(T1-T2)*Log(T1/T2))
      
      DM22 = ht**2*Nc*t*(Log(T1/q)**2+Log(T2/q)**2-2.*Log(t/q)**2)
     $     +ht**2*Nc*s2t*At*Sqrt(t)*(Log(T1/q)**2-Log(T2/q)**2)
     $     +ht**2*Nc*s2t**2*At**2/4.*(Log(T1/q)+Log(T2/q)-2.)
     $     *(2.-(T1+T2)/(T1-T2)*Log(T1/T2))
      
      return
      end

*     
***********************************************************************
*


      subroutine taudfuncs(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,A,v2,
     $     DF1,DF2,DF3,DdmuF2,DdmuF3,DdAtF2,DdAtF3)
      
c     shift of the parameters from DRbar to On-Shell scheme
 
      implicit none      
      real*8 t,mu2,A0,BL,T1,T2,s2t,c2t,cb,sb,q,A,mu,Xt,Yt
      real*8 myB0,myAA
      real*8 DF1,DF2,DF3,DdmuF2,DdmuF3,DdAtF2,DdAtF3
      real*8 Pi/3.141592654/,mt,ct2,st2,v2,v22,Nc/1d0/ ! color factor !!!
      real*8 F1o,F2o,F3o,dm1,dm2,dmt,dAt,dth,ds2t,dv2,dv22,dmu,dcotb
      real*8 pi12_1,pi12_2
      
      mu2 = mu**2

      Xt = A + mu*cb/sb
      Yt = A - mu*sb/cb

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0

      v22 = v2*sb**2
      
      mt = Sqrt(t)
      
      F1o = Log(T1/q) + Log(T2/q) - 2d0*Log(t/q)
      F2o = Log(T1/q) - Log(T2/q) 
      F3o = 2d0 - (T1+T2)/(T1-T2)*(Log(T1/q) - Log(T2/q))

c     counterterms:

      dv2 = v22 *Nc/2d0* (2d0 *Log(t/q) - 1d0 - BL/t + 
     $     T1/t * ct2 * (2d0* BL/(BL-T1) * Log(BL/T1) - 1d0)+
     $     T2/t * st2 * (2d0* BL/(BL-T2) * Log(BL/T2) - 1d0))

      dcotb = 0d0               ! beta DRbar
      dv22 = v22/v2*dv2

      dmt = Sqrt(t)/2d0*((1d0-5d0/2d0*sb**2)*Log(t/q)
     $     +3d0*cb**2/2d0*A0/t*(1d0-Log(A0/q))
     $     -3d0/2d0*mu2/t*(1d0-Log(mu2/q))
     $     +5d0*sb**2-1d0+cb**2/2d0*(1d0-A0/t)*myB0(t,0d0,A0,q)
     $     +cb**2*(2d0-A0/t)*myB0(t,t,A0,q)
     $     +(T1/t*(1d0-Log(T1/q))+(t-T1+mu2)/t*myB0(t,mu2,T1,q)
     $     +T2/t*(1d0-Log(T2/q))+(t-T2+mu2)/t*myB0(t,mu2,T2,q)
     $     +BL/t*(1d0-Log(BL/q))+(t-BL+mu2)/t*myB0(t,mu2,BL,q))/2d0)

      dmu = 0d0                 ! mu DRbar
      
      dm1 = ((T1-t-mu2)*myB0(T1,t,mu2,q) - myAA(t,q)
     $     + st2*(T1-mu2)*myB0(T1,0d0,mu2,q) - (1d0+st2)*myAA(mu2,q)
     $     + cb**2*(1d0+st2)* myAA(A0,q) + st2* myAA(BL,q) +
     $     (c2t**2-(Nc-1d0)/2d0*s2t**2)*myAA(T2,q)
     $     + (Nc+1d0)/2d0* s2t**2* myAA(T1,q) + 1d0/2d0*  
     $     (sb**2* (2d0 *Sqrt(t)+ s2t* Xt)**2* myB0(T1,T1,0d0,q) +
     $     cb**2*(2d0* Sqrt(t)+ s2t* Yt)**2* myB0(T1,T1,A0,q) + 
     $     sb**2*(1d0+c2t**2)*Xt**2*myB0(T1,T2,0d0,q) + 
     $     cb**2*(1d0+c2t**2)*Yt**2*myB0(T1,T2,A0,q)) +
     $     sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2)*myB0(T1,BL,0d0,q) +
     $     cb**2*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)*myB0(T1,BL,A0,q))
      
      dm2 = ((T2-t-mu2)*myB0(T2,t,mu2,q) - myAA(t,q)
     $     + ct2*(T2-mu2)*myB0(T2,0d0,mu2,q) - (1d0+ct2)*myAA(mu2,q)
     $     + cb**2* (1d0+ct2)*myAA(A0,q) + ct2* myAA(BL,q) +
     $     (c2t**2 - (Nc-1d0)/2d0* s2t**2)*myAA(T1,q)
     $     +(Nc+1d0)/2d0* s2t**2* myAA(T2,q) + 1d0/2d0* 
     $     (sb**2* (2d0* Sqrt(t)- s2t* Xt)**2* myB0(T2,T2,0d0,q) +
     $     cb**2*(2d0*Sqrt(t)- s2t*Yt)**2*myB0(T2,T2,A0,q) + 
     $     sb**2*(1d0+c2t**2)* Xt**2* myB0(T2,T1,0d0,q) + 
     $     cb**2*(1d0+c2t**2)* Yt**2* myB0(T2,T1,A0,q)) +
     $     sb**2*(st2*t-s2t*Sqrt(t)*Xt+ct2*Xt**2)*myB0(T2,BL,0d0,q) +
     $     cb**2*(st2*t-s2t*Sqrt(t)*Yt+ct2*Yt**2)*myB0(T2,BL,A0,q))

      pi12_1 =1d0/2d0*(s2t*((T1-mu2)*myB0(T1,0d0,mu2,q)-myAA(mu2,q))
     $     + s2t* cb**2* myAA(A0,q) + s2t* myAA(BL,q) +
     $     (Nc+1d0)* c2t* s2t* (myAA(T1,q) - myAA(T2,q)) +
     $     sb**2*c2t*Xt*(2d0* Sqrt(t) + s2t* Xt)* myB0(T1,T1,0d0,q) +
     $     cb**2*c2t*Yt*(2d0* Sqrt(t) + s2t* Yt)* myB0(T1,T1,A0,q) +
     $     sb**2*c2t*Xt*(2d0* Sqrt(t) - s2t* Xt)* myB0(T1,T2,0d0,q) +
     $     cb**2*c2t*Yt*(2d0* Sqrt(t) - s2t* Yt)* myB0(T1,T2,A0,q) -
     $     sb**2*(s2t*t-2d0*c2t*Sqrt(t)*Xt-s2t*Xt**2)
     $     *myB0(T1,BL,0d0,q) -
     $     cb**2*(s2t*t-2d0*c2t*Sqrt(t)*Yt-s2t*Yt**2)
     $     *myB0(T1,BL,A0,q))

      pi12_2 =1d0/2d0*(s2t*((T2-mu2)*myB0(T2,0d0,mu2,q)-myAA(mu2,q))
     $     + s2t* cb**2* myAA(A0,q) + s2t* myAA(BL,q) +
     $     (Nc+1d0)* c2t* s2t* (myAA(T1,q) - myAA(T2,q)) +
     $     sb**2*c2t*Xt*(2d0* Sqrt(t) + s2t* Xt)* myB0(T2,T1,0d0,q) +
     $     cb**2*c2t*Yt*(2d0* Sqrt(t) + s2t* Yt)* myB0(T2,T1,A0,q) +
     $     sb**2*c2t*Xt*(2d0* Sqrt(t) - s2t* Xt)* myB0(T2,T2,0d0,q) +
     $     cb**2*c2t*Yt*(2d0* Sqrt(t) - s2t* Yt)* myB0(T2,T2,A0,q) -
     $     sb**2*(s2t*t-2d0*c2t*Sqrt(t)*Xt-s2t*Xt**2)
     $     *myB0(T2,BL,0d0,q) -
     $     cb**2*(s2t*t-2d0*c2t*Sqrt(t)*Yt-s2t*Yt**2)
     $     *myB0(T2,BL,A0,q))

      dth = (pi12_1 + pi12_2)/2d0/(T1-T2)

      ds2t = 2d0*c2t*dth

      dAt = ((dm1-dm2)/(T1-T2) + ds2t/s2t - dmt/mt)*Xt
     $     - mu * dcotb - dmu * cb/sb

      DF1 = dm1/T1 + dm2/T2 - 4d0*dmt/mt + (4d0*dmt/mt - dv22/v22)*F1o
      DF2 = dm1/T1 - dm2/T2 + (3d0*dmt/mt - dv22/v22 + ds2t/s2t)*F2o
      DF3 = (2d0*T1*T2/(T1-T2)**2*Log(T1/T2) - (T1+T2)/(T1-T2)) 
     $     *(dm1/T1-dm2/T2) + (2d0*dmt/mt-dv22/v22+2d0*ds2t/s2t)*F3o

      DdmuF2 = dmu/mu * F2o       
      DdmuF3 = dmu/mu * F3o       

      DdAtF2 = dAt/A * F2o       
      DdAtF3 = dAt/A * F3o       

c     residues of some singular functions for s2t=0 and for A=0

      if(s2t.eq.0d0) then         
         DF2 = ds2t*F2o
         DdAtF2 = ds2t*Xt/A
      endif

      if(mu.eq.0d0) then
         DdmuF2 = dmu * F2o       
         DdmuF3 = dmu * F3o       
      endif

      if(A.eq.0d0) then
         DdAtF2 = dAt * F2o       
         DdAtF3 = dAt * F3o       
      endif
      
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine SU_tausqodd(t,A0,BL,T1,T2,st,ct,q,mu,tb,v2,DMA)

c     Two-loop O(a_tau^2) corrections to the CP-odd Higgs mass in the
c     DRbar scheme. Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on A. Brignole, G. Degrassi, P. Slavich and F. Zwirner, 
c     hep-ph/0112177 with appropriate substitutions for top -> tau.
c
c     Last update:  19/06/2003.
c
c
c     I/O PARAMETERS:
c     t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
c     st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tb = tan(beta), v2 = v^2, 
c     DMA = 2-loop corrections to the CP-odd Higgs mass.

      implicit none

      real*8 ht,k,mt,pi,v2,tb
      real*8 t,mu2,A0,BL,T1,T2,st,ct,q,A,X,mu,tanb,sb,cb,s2t,c2t
      real*8 FA,FA_A,DMA

      pi = 3.14159265897d0

      tanb = 1d0/tb           ! swap H1 <-> H2 !!!

      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt  
      A = X - mu/tanb         

      sb = dsin(datan(tanb))
      cb = dcos(datan(tanb))

      ht = dsqrt(2d0/v2)*mt/sb
      
c$$$      k = 3d0*ht**2/(16d0*Pi**2)**2 
      k = ht**2/(16d0*Pi**2)**2 ! remove color factor !!!
    
      call taufuncsodd(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,FA,FA_A)

      if(A.ne.0d0) then

         DMA = ht**2*mu*A/(T1-T2)/sb/cb * FA

      else

c     the function FA has poles in A=0. 
c     when necessary we consider the residues:

         DMA = ht**2*mu/(T1-T2)/sb/cb * FA_A

      endif

      DMA = k*DMA

      return
      end

*
***********************************************************************
*
      
      subroutine taufuncsodd(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu,FA,FA_A)

      implicit none

      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 FA,FA_A,tauFAab,tauFAc,tauresFAc

      FA = tauFAab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu) 
     $     + tauFAc(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
     $     - tauFAc(t,A0,BL,T2,T1,-s2t,-c2t,cb,sb,q,mu)
      
      FA_A = tauresFAc(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)
     $     - tauresFAc(t,A0,BL,T2,T1,-s2t,-c2t,cb,sb,q,mu)
 
      return
      end

*
*********************************************************************
*
            
      function tauFAab(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 tauFAab


      tauFAab = (5.+2.*Nc)*(T1*(1.-Log(T1/q)+Log(T1/q)**2/2)-
     $     T2*(1.-Log(T2/q)+Log(T2/q)**2/2))
      
      return
      end
      
*
*********************************************************************
*

      function tauFAc(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Nc/1d0/ ! color factor !!!
      real*8 mu2,Xt,Yt,At,st2,ct2
      real*8 tauFAc,phi,delt,Li2

      mu2 = mu**2
      if(mu2.eq.0d0) mu2 = 1d-10

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb
      At = sb**2*Xt+cb**2*Yt

      ct2 = (1d0+c2t)/2d0
      st2 = (1d0-c2t)/2d0

      tauFAc =  -(delt(mu2,t,T1)+2*t*mu2)/T1*phi(mu2,t,T1)
     $ +cb**2*Yt**2*(A0*c2t**2*(A0-2*(T1+T2))/2/(T1-T2)/T2 
     $ -s2t**2*(T1-T2)/4/T2)*phi(A0,T1,T2)
     $ +cb**2*(-A0*Sqrt(t)*Yt*(A0-4*T1)
     $ *(1./2/At/T1/Sqrt(t)-s2t/T1/(T1-T2))
     $ -A0*c2t**2*Yt**2*(A0-4*T1)/2/T1/(T1-T2)
     $ +A0*(2*Sqrt(t)+s2t*Yt)**2/4/T1)*phi(A0,T1,T1)
     $ +cb**2*(-Sqrt(t)*Yt*delt(A0,BL,T1)*
     $ (1./2/At/T1/Sqrt(t)-s2t/T1/(T1-T2))
     $ +c2t*delt(A0,BL,T1)*(t-Yt**2)/2/T1/(T1-T2)
     $ +(A0+BL-T1)*(ct2*t+s2t*Sqrt(t)*Yt+st2*Yt**2)/2/T1)*phi(A0,BL,T1)
     $ -s2t**2*sb**2*Xt**2/2*Li2(1-T2/T1)
     $ +(mu2-T1)*(1-c2t*(1+(mu2-T1)/(T1-T2)))*Li2(1-mu2/T1)
     $ +(2*sb**2*Sqrt(t)*Xt*(BL-T1)
     $ *(1./2/At/Sqrt(t)-s2t/(T1-T2))
     $ -c2t*sb**2*(t-Xt**2)*(BL-T1)/(T1-T2)
     $ -sb**2*(ct2*t+s2t*Sqrt(t)*Xt+st2*Xt**2))*Li2(1-BL/T1)
     $ -21*Sqrt(t)*T1*(1./2/Sqrt(t)-s2t*At/(T1-T2))
     $ -Nc*c2t**2*T1-9*s2t*Sqrt(t)*At/2+c2t/2*BL
     $ -c2t*T1/(T1-T2)*(A0*cb**2+2*BL-2*t
     $ -3*mu2+2*sb**2*Xt**2+2*cb**2*Yt**2)
     $ -T1/2-T1*(2*T1+5*T2)/2/(T1-T2)*c2t-T1*c2t**2-(3+Nc)*(T1-T2)
     $ +mu2*c2t*T2/(T1-T2)*Log(mu2/q)
     $ +BL*c2t*(3*T1+T2)/2/(T1-T2)*Log(BL/q)
     $ +A0*cb**2*c2t*(3*T1+T2)/2/(T1-T2)*Log(A0/q)
     $ +Log(T2/q)*(-3*s2t*Sqrt(t)*T2*At/(T1-T2)
     $ +c2t*T2/2/(T1-T2)*(BL+cb**2*A0-2*mu2-t+T1+
     $ sb**2*Xt**2+cb**2*Yt**2)
     $ +((s2t**2-2)*T1+3*(3*s2t**2-4)*T2)/4/(T1-T2)
     $ *(sb**2*Xt**2+cb**2*Yt**2)
     $ +Nc*(T1+T2)/4/(T1-T2)*(s2t**2*(T1+2*T2)-4*T2)
     $ +s2t**2/4/(T1-T2)*(T1**2+3*T1*T2+2*T2**2)
     $ -(T1**2-3*T1*T2+9*T2**2)/2/(T1-T2)-(3+Nc)*T2
     $ -T2*(BL+3*A0*cb**2-6*mu2+t)/2/(T1-T2))
     $ +Log(T1/q)*(18*Sqrt(t)*T1*(1./2/Sqrt(t)-s2t*At/(T1-T2))
     $ +3*s2t*Sqrt(t)*(2*T1-T2)/(T1-T2)*At
     $ +((9+8*c2t+9*c2t**2)*T1+(-5+2*c2t+c2t**2)*T2)/
     $ 4/(T1-T2)*(sb**2*Xt**2+cb**2*Yt**2)
     $ +Nc*(4*(3-2*s2t**2)*T1**2
     $ +(-4+s2t**2)*T1*T2+s2t**2*T2**2)/4/(T1-T2)
     $ +c2t**2*(8*T1**2-T1*T2-T2**2)/4/(T1-T2)
     $ +cb**2*A0*(2*(3+c2t)*T1-(3-c2t)*T2)/2/(T1-T2)
     $ +BL*(2*(1+c2t)*T1-(1-c2t)*T2)/2/(T1-T2)
     $ -mu2*(3*(2+c2t)*T1-(3-c2t)*T2)/(T1-T2)
     $ +t*(2*(1-2*c2t)*T1-(1+c2t)*T2)/2/(T1-T2)
     $ +c2t*T1*(T1+4*T2)/2/(T1-T2)+(3+Nc)*T1
     $ +(20*T1**2-11*T1*T2-T2**2)/4/(T1-T2))
     $ +T1*Log(mu2/q)*Log(t/q)+(mu2+t-T1)*Log(T1/q)*Log(t/q)
     $ +(mu2+t-T1+c2t*mu2**2/(T1-T2))*Log(T1/q)*Log(mu2/q)
     $ +Log(BL/q)*Log(T1/q)*
     $ (cb**2*(A0-BL-T1)*Yt/At/2-BL*sb**2*Xt/At
     $ +2*BL*Sqrt(t)*s2t*At/(T1-T2)
     $ -s2t*cb**2*Sqrt(t)*Yt*(A0+BL-T1)/(T1-T2)
     $ -cb**2*Sqrt(t)*Yt*s2t/2
     $ -c2t*cb**2*(A0+BL-T1)*(t-Yt**2)/2/(T1-T2)
     $ +c2t*BL*(t-T1-sb**2*Xt**2-cb**2*Yt**2)/(T1-T2)
     $ -(1-c2t)/2*BL-cb**2*((1+c2t)*t+(1-c2t)*Yt**2)/4)
     $ +cb**2*Log(A0/q)*Log(T1/q)*(
     $ -Sqrt(t)*(3*A0-BL+T1)*Yt*(1./2/At/Sqrt(t)-s2t/(T1-T2))
     $ -Sqrt(t)*Yt*s2t/2+c2t**2*Yt**2*(T1-T2-2*A0)/4/(T1-T2)
     $ -(6*A0+t+2*Yt**2)/4-c2t*(2*BL-2*A0-T1-T2)*(t-Yt**2)/4/(T1-T2)
     $ -c2t*A0*(T1+T2)/2/(T1-T2))
     $ +cb**2*Log(A0/q)*Log(BL/q)*
     $ (Sqrt(t)*T1*Yt*(1./2/At/Sqrt(t)-s2t/(T1-T2))
     $ +s2t*Sqrt(t)*Yt/2-c2t*(t-Yt**2)*(T1+T2)/2/(T1-T2))
     $ +cb**2*(A0*c2t**2*Yt**2/2/(T1-T2)
     $ +s2t**2*Yt**2/4)*Log(A0/q)*Log(T2/q)
     $ +Log(T2/q)**2*(3*s2t*Sqrt(t)*T2*At/2/(T1-T2)+
     $ (2-c2t+c2t**2)*T2/4/(T1-T2)*(sb**2*Xt**2+cb**2*Yt**2)
     $ +Nc*(1+c2t**2)*T2**2/4/(T1-T2)-s2t**2*T2**2/4/(T1-T2)
     $ +BL*T2*(1-c2t)/4/(T1-T2)-(3-c2t)*mu2*T2/2/(T1-T2)
     $ +(3-c2t)*A0*T2*cb**2/4/(T1-T2)+(1+c2t)*t*T2/4/(T1-T2)
     $ +T2*(9*T2-(4+c2t)*T1)/4/(T1-T2)+(3+Nc)*T2/2)
     $ +Log(T1/q)*Log(T2/q)*(
     $ Nc*(4*T1*T2-s2t**2*(T1+T2)**2)/4/(T1-T2)+(1+c2t**2)/4*(T1-T2)
     $ +cb**2*c2t**2*Yt**2*(T1+T2-A0)/2/(T1-T2)
     $ +c2t**2*T2*(T1+sb**2*Xt**2)/(T1-T2))
     $ +Log(T1/q)**2*(
     $ cb**2*Sqrt(t)*Yt*(A0-BL+T1)*(1./2/At/Sqrt(t)-s2t/(T1-T2))
     $ +Sqrt(t)*(BL-6*T1)*(1./2/Sqrt(t)-s2t*At/(T1-T2))
     $ +cb**2*s2t*Sqrt(t)*Yt/2-3*At*Sqrt(t)*(2*T1-T2)*s2t/2/(T1-T2)
     $ +Nc*T1*(3*(-2+s2t**2)*T1+2*s2t**2*T2)/4/(T1-T2)
     $ -sb**2*Xt**2/4/(T1-T2)*((5+2*c2t+3*c2t**2)*T1+
     $ (-3+c2t+2*c2t**2)*T2-2*BL*c2t)
     $ -cb**2*Yt**2/4/(T1-T2)*((3+c2t+4*c2t**2)*T1
     $ -s2t**2*T2-2*A0*c2t**2)
     $ -(3+c2t)*cb**2*A0*T1/4/(T1-T2)-T1*BL*(1+c2t)/4/(T1-T2)
     $ -t*c2t*sb**2*BL/2/(T1-T2)
     $ +mu2/2/(T1-T2)*(2*(2+c2t)*T1-(1-c2t)*T2)-c2t*mu2**2/2/(T1-T2)
     $ +t/4/(T1-T2)*((-5+c2t-(1-c2t)*sb**2)*T1+(4+(1+c2t)*sb**2)*T2)
     $ -T1*((4+3*c2t**2)*T1+c2t*(3+2*c2t)*T2)/4/(T1-T2)-(3+Nc)*T1/2)

      return
      end

*
*********************************************************************
*

      function tauresFAc(t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu)

      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,cb,sb,q,mu
      real*8 Pi/3.141592654/,Xt,Yt
      real*8 tauresFAc,phi,delt,Li2

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb

      tauresFAc =
     $     -cb**2*A0*Yt*(A0-4*T1)/2/T1*phi(A0,T1,T1)
     $     -cb**2*Yt*delt(A0,BL,T1)/2/T1*phi(A0,BL,T1)
     $     +sb**2*Xt*(BL-T1)*Li2(1-BL/T1)
     $     +Log(BL/q)*Log(T1/q)*(cb**2*(A0-BL-T1)*Yt/2-BL*sb**2*Xt)
     $     -cb**2*Log(A0/q)*Log(T1/q)*(3*A0-BL+T1)*Yt/2
     $     +cb**2*Log(A0/q)*Log(BL/q)*T1*Yt/2
     $     +Log(T1/q)**2*cb**2*Yt*(A0-BL+T1)/2
      
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine SU_tausqtad(t,A0,BL,T1,T2,st,ct,q,mu,tb,vv,
     $     S1,S2)

c     Two-loop O(a_t^2) corrections to the Higgs tadpoles. 
c     Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
c     Based on A. Dedes and P. Slavich, hep-ph/0212132 
c     with appropriate substitutions for top -> tau.
c
c     Last update:  19/06/2003.
c
c     I/O PARAMETERS:
c     t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
c     st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
c     mu = Higgs mixing parameter, tb = tan(beta), vv = v^2,
c     Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles.
c
c     Notice: we assume that the 1-loop part is computed in terms of 
c             running (DRbar) parameters, evaluated at the scale Q.

      implicit none

      real*8 ht,k,mt,pi/3.14159265897d0/,vv,v1,v2,tb
      real*8 t,mu2,A0,BL,T1,T2,st,ct,q,A,X,mu,tanb,sb,cb,s2t,c2t
      real*8 F2l,G2l,S1,S2,sw

      tanb = 1d0/tb           ! swap H1 <-> H2 !!!
      
      mt = dsqrt(t)

      s2t = 2d0*ct*st
      c2t = ct**2 - st**2

      X = (T1-T2)*s2t/2d0/mt
      A = X - mu/tanb           ! notice the sign convention for mu

      sb = dsin(datan(tanb))
      cb = dcos(datan(tanb))
      v1 = sqrt(vv)*cb
      v2 = sqrt(vv)*sb

      ht = dsqrt(2d0/vv)*mt/sb
      
c$$$      k = 3d0*ht**2/(16d0*Pi**2)**2 
      k = ht**2/(16d0*Pi**2)**2 ! remove color factor !!!
 

      call tautadfuncs(t,A0,BL,T1,T2,s2t,c2t,q,mu,sb,cb,F2l,G2l)

      S1 = mt * mu/tanb * s2t * F2l
      S1 = S1/v1**2

      S2 = mt * A * s2t * F2l + 2d0 * mt**2 * G2l
      S2 = S2/v2**2

      S1 = k*S1
      S2 = k*S2

      sw = S1                  ! swap H1 <-> H2 !!!
      S1 = S2
      S2 = sw

      return
      end

*
***********************************************************************
*

      
      subroutine tautadfuncs(t,A0,BL,T1,T2,s2t,c2t,q,mu,sb,cb,F2l,G2l)

      implicit none
      real*8 t,A0,BL,T1,T2,s2t,c2t,q,mu,mu2,sb,cb,sb2,cb2
      real*8 F2l,G2l,phi,delt,Li2,Xt,Yt,At
      real*8 Pi/3.14159265897d0/,Nc/1d0/ ! color factor !!!

      sb2 = sb**2
      cb2 = cb**2
      mu2 = mu**2

      if(mu2.eq.0d0) mu2 = 1d-10

      Xt = s2t*(T1-T2)/2d0/Sqrt(t)
      Yt = Xt - mu/cb/sb
      At = sb2*Xt+cb2*Yt
      
      F2l = (delt(mu2,t,T1)+2*mu2*t)/T1*phi(mu2,t,T1)
     $     -(delt(mu2,t,T2)+2*mu2*t)/T2*phi(mu2,t,T2)
     $     +A0*cb2*(2*Sqrt(t)+s2t*Yt)/4/s2t/T1/(T1-T2)*
     $     (-2*s2t*Sqrt(t)*(T1-T2)+2*(A0-4*T1)*Yt
     $     +s2t**2*(7*T1+T2-2*A0)*Yt)*phi(A0,T1,T1)
     $     -A0*cb2*(2*Sqrt(t)-s2t*Yt)/4/s2t/T2/(T1-T2)*
     $     (-2*s2t*Sqrt(t)*(T1-T2)+2*(A0-4*T2)*Yt
     $     +s2t**2*(7*T2+T1-2*A0)*Yt)*phi(A0,T2,T2)
     $     +cb2*Yt**2/2/(T1-T2)/T2*(s2t**2*(T1-T2)**2-2*A0**2*c2t**2
     $     +A0*4*c2t**2*(T1+T2))*phi(A0,T1,T2)
     $     +(-cb2*delt(A0,BL,T1)/2/s2t/T1/(T1-T2)
     $     *(c2t*s2t*(t-Yt**2)-2*c2t**2*Sqrt(t)*Yt)
     $     -cb2/4/T1*(A0+BL-T1)*((1+c2t)*t+2*s2t*Sqrt(t)*Yt
     $     +Yt**2*(1-c2t)))*phi(A0,BL,T1)
     $     -(-cb2*delt(A0,BL,T2)/2/s2t/T2/(T1-T2)
     $     *(c2t*s2t*(t-Yt**2)-2*c2t**2*Sqrt(t)*Yt)
     $     -cb2/4/T2*(A0+BL-T2)*((1-c2t)*t-2*s2t*Sqrt(t)*Yt
     $     +Yt**2*(1+c2t)))*phi(A0,BL,T2)
     $     +(mu2-T1)*(T2-T1+c2t*(mu2-T2))/(T1-T2)*Li2(1-mu2/T1)
     $     +(mu2-T2)*(T1-T2-c2t*(mu2-T1))/(T1-T2)*Li2(1-mu2/T2)
     $     +(-sb2*(BL-T1)*(c2t*s2t*(t-Xt**2)-2*c2t**2*Sqrt(t)*Xt)
     $     /s2t/(T1-T2)
     $     -sb2/2*((1+c2t)*t+2*s2t*Sqrt(t)*Xt+(1-c2t)*Xt**2))
     $     *Li2(1-T1/BL)
     $     -(-sb2*(BL-T2)*(c2t*s2t*(t-Xt**2)-2*c2t**2*Sqrt(t)*Xt)
     $     /s2t/(T1-T2)
     $     -sb2/2*((1-c2t)*t-2*s2t*Sqrt(t)*Xt+(1+c2t)*Xt**2))
     $     *Li2(1-T2/BL)
     $     +s2t**2*sb2*Xt**2/2*(Li2(1-T2/T1)-Li2(1-T1/T2))	       
     $     +(-((-2+s2t**2)*T1+s2t**2*T2)/2/(T1-T2)
     $     *(sb2*Xt**2+cb2*Yt**2+2*Sqrt(t)*At/s2t)
     $     -(2*A0*c2t**2+(T1-T2)*s2t**2)/4/(T1-T2)*cb2*Yt**2
     $     -A0*cb2*c2t**2*Sqrt(t)*Yt/s2t/(T1-T2)
     $     -(Nc+1)*T1*((-3+s2t**2)*T1+(1+s2t**2)*T2)/2/(T1-T2)-Nc/2*T1
     $     +c2t*(mu2-T1)*(mu2-T2)/2/(T1-T2)-mu2/2+t)*Log(T1/q)**2
     $     -(((-2+s2t**2)*T2+s2t**2*T1)/2/(T1-T2)
     $     *(sb2*Xt**2+cb2*Yt**2-2*Sqrt(t)*At/s2t)
     $     +(2*A0*c2t**2-(T1-T2)*s2t**2)/4/(T1-T2)*cb2*Yt**2
     $     -A0*cb2*c2t**2*Sqrt(t)*Yt/s2t/(T1-T2)
     $     +(Nc+1)*T2*((-3+s2t**2)*T2+(1+s2t**2)*T1)/2/(T1-T2)-Nc/2*T2
     $     +c2t*(mu2-T1)*(mu2-T2)/2/(T1-T2)-mu2/2+t)*Log(T2/q)**2
     $     +(-c2t**2*(T1+T2)/(T1-T2)*(sb2*Xt**2+cb2*Yt**2)
     $     +c2t**2*cb2*A0*Yt**2/(T1-T2)-2*(1+Nc)*c2t**2*T1*T2/(T1-T2)
     $     +((1+Nc)*s2t**2-2)/2*(T1-T2))*Log(T1/q)*Log(T2/q)
     $     +(T1-mu2-t)*Log(t/q)*Log(T1/q)-(T2-mu2-t)*Log(t/q)*Log(T2/q)
     $     +((T1-mu2-t)-c2t*mu2**2/(T1-T2))*Log(mu2/q)*Log(T1/q)
     $     -((T2-mu2-t)-c2t*mu2**2/(T1-T2))*Log(mu2/q)*Log(T2/q)
     $     -cb2*Sqrt(t)*Yt/s2t*Log(A0/q)*Log(BL/q)
     $     -(T1-T2)*Log(mu2/q)*Log(t/q)-sb2*Sqrt(t)*Xt/s2t*Log(BL/q)**2
     $     +(-cb2*Sqrt(t)*Yt/2/s2t/(T1-T2)*
     $     (2*BL*c2t**2-6*A0*c2t**2+(-2+s2t**2)*T1+s2t**2*T2)
     $     +c2t*cb2/4/(T1-T2)*((2*BL-T1-T2)*(t-Yt**2)
     $     +2*A0*(T1+T2+Yt**2-t))
     $     +cb2/4*(6*A0+t+(1+2*s2t**2)*Yt**2))*Log(A0/q)*Log(T1/q)
     $     -(-cb2*Sqrt(t)*Yt/2/s2t/(T1-T2)*
     $     (2*BL*c2t**2-6*A0*c2t**2+(-2+s2t**2)*T2+s2t**2*T1)
     $     +c2t*cb2/4/(T1-T2)*((2*BL-T1-T2)*(t-Yt**2)
     $     +2*A0*(T1+T2+Yt**2-t))
     $     +cb2/4*(6*A0+t+(1+2*s2t**2)*Yt**2))*Log(A0/q)*Log(T2/q)
     $     +(cb2*Sqrt(t)*Yt/2/s2t/(T1-T2)*
     $     (2*BL*c2t**2-2*A0*c2t**2+(-2+s2t**2)*T1+s2t**2*T2)
     $     -Sqrt(t)*At/s2t/(T1-T2)*((-2+s2t**2)*T1+s2t**2*T2)
     $     +c2t/4/(T1-T2)*(2*A0*cb2*(t-Yt**2)+(T1+T2)*
     $     (t*(-2+cb2)+2*sb2*Xt**2+cb2*Yt**2)
     $     +2*BL*(T1+T2-cb2*(t-Yt**2)))
     $     +(2*BL-(-2+cb2)*t+2*sb2*Xt**2+cb2*Yt**2)/4)
     $     *Log(BL/q)*Log(T1/q)
     $     -(cb2*Sqrt(t)*Yt/2/s2t/(T1-T2)*
     $     (2*BL*c2t**2-2*A0*c2t**2+(-2+s2t**2)*T2+s2t**2*T1)
     $     -Sqrt(t)*At/s2t/(T1-T2)*((-2+s2t**2)*T2+s2t**2*T1)
     $     +c2t/4/(T1-T2)*(2*A0*cb2*(t-Yt**2)+(T1+T2)*
     $     (t*(-2+cb2)+2*sb2*Xt**2+cb2*Yt**2)
     $     +2*BL*(T1+T2-cb2*(t-Yt**2)))
     $     +(2*BL-(-2+cb2)*t+2*sb2*Xt**2+cb2*Yt**2)/4)
     $     *Log(BL/q)*Log(T2/q)
     $     +(3*Sqrt(t)*At/s2t/(T1-T2)*((-4+3*s2t**2)*T1+s2t**2*T2)
     $     -(3*(1+c2t)*T1+(-3+c2t)*T2)/2/(T1-T2)*(sb2*Xt**2+cb2*Yt**2)
     $     -c2t/2/(T1-T2)*((T1+T2)*(BL+A0*cb2)
     $     -T1*(4*mu2+3*t-T1)-T2*(2*mu2+t-3*T1))
     $     -BL/2-3./2.*A0*cb2+3*mu2-t/2
     $     +Nc/2*((-4+3*s2t**2)*T1+s2t**2*T2)
     $     +((-13+3*s2t**2)*T1+(-2+s2t**2)*T2)/2)*Log(T1/q)
     $     -(3*Sqrt(t)*At/s2t/(T1-T2)*((-4+3*s2t**2)*T2+s2t**2*T1)
     $     +(3*(1-c2t)*T2+(-3-c2t)*T1)/2/(T1-T2)*(sb2*Xt**2+cb2*Yt**2)
     $     -c2t/2/(T1-T2)*((T1+T2)*(BL+A0*cb2)
     $     -T2*(4*mu2+3*t-T2)-T1*(2*mu2+t-3*T2))
     $     -BL/2-3./2.*A0*cb2+3*mu2-t/2
     $     +Nc/2*((-4+3*s2t**2)*T2+s2t**2*T1)
     $     +((-13+3*s2t**2)*T2+(-2+s2t**2)*T1)/2)*Log(T2/q)
     $     +c2t*(mu2*Log(mu2/q)-BL*Log(BL/q)-A0*cb2*Log(A0/q))
     $     +3*At*Sqrt(t)/s2t*(5-4*s2t**2)
     $     +2*c2t*(sb2*Xt**2+cb2*Yt**2)+c2t*(BL+A0*cb2-3*mu2-2*t+T1+T2)
     $     -(2*Nc*(-1+s2t**2)-11+2*s2t**2)/2*(T1-T2)

      G2l =
     $     -2*A0*cb2*(A0-3*t)/t*phi(A0,t,t)
     $     +(2*BL+mu2-BL**2/t+BL*mu2/t-t)*phi(BL,mu2,t)
     $     +mu2/T1*(mu2+t-T1)*phi(mu2,t,T1)
     $     +mu2/T2*(mu2+t-T2)*phi(mu2,t,T2)
     $     +A0*cb2*(2*Sqrt(t)+s2t*Yt)*(A0-2*t-4*T1-s2t*Sqrt(t)*Yt)
     $     /4/Sqrt(t)/T1*phi(A0,T1,T1)
     $     +A0*cb2*(2*Sqrt(t)-s2t*Yt)*(A0-2*t-4*T2+s2t*Sqrt(t)*Yt)
     $     /4/Sqrt(t)/T2*phi(A0,T2,T2)
     $     +A0*cb2*(-2+s2t**2)*Yt**2/2/T2*phi(A0,T1,T2)
     $     +((s2t*Yt/Sqrt(t)+1+c2t)*delt(A0,BL,T1)
     $     -(2*s2t*Sqrt(t)*Yt+(1-c2t)*Yt**2+(1+c2t)*t)
     $     *(A0+BL-T1))*cb2/4/T1*phi(A0,BL,T1)
     $     +((-s2t*Yt/Sqrt(t)+1-c2t)*delt(A0,BL,T2)
     $     -(-2*s2t*Sqrt(t)*Yt+(1+c2t)*Yt**2+(1-c2t)*t)
     $     *(A0+BL-T2))*cb2/4/T2*phi(A0,BL,T2)
     $     +2*cb2*(t-A0)*Li2(1-t/A0)
     $     +(-1+c2t)*(mu2-T1)*Li2(1-mu2/T1)
     $     +(-1-c2t)*(mu2-T2)*Li2(1-mu2/T2)
     $     +(-sb2*s2t*(2*t+T1-BL)*Xt/2/Sqrt(t)-sb2/2*(1-c2t)*Xt**2
     $     +sb2/2*(1+c2t)*(BL-t-T1))*Li2(1-T1/BL)
     $     +(sb2*s2t*(2*t+T2-BL)*Xt/2/Sqrt(t)-sb2/2*(1+c2t)*Xt**2
     $     +sb2/2*(1-c2t)*(BL-t-T2))*Li2(1-T2/BL)
     $     -sb2*(s2t**2-2)*Xt**2/2*(Li2(1-T1/T2)+Li2(1-T2/T1))
     $     +cb2*(t-A0)*Log(A0/q)**2-2*cb2*(t+2*A0)*Log(A0/q)*Log(t/q)
     $     +(cb2*(t+2*A0)-3*t)*Log(t/q)**2-(BL-t)*Log(BL/q)*Log(mu2/q)
     $     +(BL+2*mu2-t)*Log(mu2/q)*Log(t/q)-(BL+t)*Log(BL/q)*Log(t/q)
     $     +(-sb2*s2t*(T1-T2)*Xt/4/Sqrt(t)-sb2/2*Xt**2
     $     +sb2/4*(2*BL-2*t-(1+c2t)*T1-(1-c2t)*T2))*Log(BL/q)**2
     $     +(-cb2*s2t*(T1-T2)*Yt/4/Sqrt(t)-cb2/2*(Yt**2-A0)
     $     +cb2/4*(2*BL-2*t-(1+c2t)*T1-(1-c2t)*T2))*Log(A0/q)*Log(BL/q)
     $     +(s2t*(2*t+T1)/2/Sqrt(t)*At-cb2*s2t*A0*Yt/4/Sqrt(t)
     $     +sb2/2*Xt**2+cb2*s2t**2/4*Yt**2-A0*cb2/2-mu2/2*(1-c2t)
     $     +t+(3-c2t+(1+Nc)*s2t**2)/2*T1)*Log(T1/q)**2
     $     +(-s2t*(2*t+T2)/2/Sqrt(t)*At+cb2*s2t*A0*Yt/4/Sqrt(t)
     $     +sb2/2*Xt**2+cb2*s2t**2/4*Yt**2-A0*cb2/2-mu2/2*(1+c2t)
     $     +t+(3+c2t+(1+Nc)*s2t**2)/2*T2)*Log(T2/q)**2
     $     +(cb2*(2-s2t**2)*Yt**2/2+(2-(1+Nc)*s2t**2)/2*(T1+T2))
     $     *Log(T1/q)*Log(T2/q)
     $     -(mu2+2*t)*Log(t/q)*Log(T1/q)-mu2*Log(mu2/q)*Log(T1/q)
     $     -(mu2+2*t)*Log(t/q)*Log(T2/q)-mu2*Log(mu2/q)*Log(T2/q)
     $     +(cb2*s2t*(3*A0-BL+2*t+T1)*Yt/4/Sqrt(t)+cb2/4*(1-c2t)*Yt**2
     $     -cb2/4*((1+c2t)*(BL-t-T1)-(11-c2t)*A0))*Log(A0/q)*Log(T1/q)
     $     +(-cb2*s2t*(3*A0-BL+2*t+T2)*Yt/4/Sqrt(t)+cb2/4*(1+c2t)*Yt**2
     $     -cb2/4*((1-c2t)*(BL-t-T2)-(11+c2t)*A0))*Log(A0/q)*Log(T2/q)
     $     +(cb2*s2t*Yt/4/Sqrt(t)*(BL-A0+2*t+T1)
     $     +sb2*s2t*Xt/2/Sqrt(t)*(2*t+T1)+sb2/2*(1-c2t)*Xt**2
     $     +cb2/4*(1-c2t)*Yt**2-cb2/4*A0*(1+c2t)
     $     +BL/4*(2*(1-c2t)+cb2*(1+c2t))
     $     -(1+c2t)*(-2+cb2)/4*(t+T1))*Log(BL/q)*Log(T1/q)
     $     +(-cb2*s2t*Yt/4/Sqrt(t)*(BL-A0+2*t+T2)
     $     -sb2*s2t*Xt/2/Sqrt(t)*(2*t+T2)+sb2/2*(1+c2t)*Xt**2
     $     +cb2/4*(1+c2t)*Yt**2-cb2/4*A0*(1-c2t)
     $     +BL/4*(2*(1+c2t)+cb2*(1-c2t))
     $     -(1-c2t)*(-2+cb2)/4*(t+T2))*Log(BL/q)*Log(T2/q)
     $     +(-3*s2t*(t+T1)/Sqrt(t)*At
     $     -(3-c2t)/2*(sb2*Xt**2+cb2*Yt**2)-(1-c2t)/2*BL
     $     -(3-c2t)/2*cb2*A0+(3-c2t)*mu2-(1+c2t)/2*t
     $     +(c2t-15-(1+Nc)*s2t**2)/2*T1
     $     +(-2+(1+Nc)*s2t**2)/2*T2)*Log(T1/q)
     $     +(3*s2t*(t+T2)/Sqrt(t)*At
     $     -(3+c2t)/2*(sb2*Xt**2+cb2*Yt**2)-(1+c2t)/2*BL
     $     -(3+c2t)/2*cb2*A0+(3+c2t)*mu2-(1-c2t)/2*t
     $     +(-c2t-15-(1+Nc)*s2t**2)/2*T2
     $     +(-2+(1+Nc)*s2t**2)/2*T1)*Log(T2/q)
     $     -Pi**2*sb2*t/3
     $     +15*s2t*(T1-T2)/4/Sqrt(t)*At+3./2.*(sb2*Xt**2+cb2*Yt**2)
     $     +18*t*Log(t/q)+BL/2+3./2.*A0*cb2-3*mu2
     $     +35./4.*(T1+T2-2*t)-c2t/4*(T1-T2)

      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%  End of routines of Pietro %%%%%%%%%%%%%%%%%%%%%%%%%
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MDR,MTOP,AU,AD,MU,MCHI
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR/MDR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C     MCHI IS THE HEAVIEST CHARGINO MASS. 
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.

C     OUTPUT: MH,HM,MCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C     WHERE MHP AND HPM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Program based on the work by M. Carena, M. Quiros
c       and C.E.M. Wagner, "Effective potential methods and
c       the Higgs mass spectrum in the MSSM", Nucl. Phys.
c       B461 (1996) 407. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE SUBH_TLH(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
     *                 MHP,HMP,MCH,SA,CA,TANBA,MGLU)

      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/su_PARAM/GF,ALPH,AMZ,AMW
      COMMON/HSELF_HDEC/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

      MCHI = MCHI0
      TANBA = TANB
      TANBT = TANB
      
      PI = 4*DATAN(1D0)
      amtau= 1.777d0
      ammuon= 0.106d0
      amc=1.4d0
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)
      MB      = RUNM(MTOP,5)
      RMTOP   = RUNM(MTOP,6)

      TQ = LOG((MQ**2+MTOP**2)/MTOP**2)
      TU = LOG((MUR**2 + MTOP**2)/MTOP**2)
      TD = LOG((MD**2 + MTOP**2)/MTOP**2)
      SINB = TANB/DSQRT(1.D0 + TANB**2)
      COSB = SINB/TANB

      IF(MA.GT.MTOP)
     *       TANBA = TANB*(1.D0-3.D0/32.D0/PI**2*
     *       (RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *       DLOG(MA**2/MTOP**2))
      IF(MA.LT.MTOP.OR.MA.EQ.MTOP) TANBT = TANBA

      SINB = TANBT/DSQRT(1.D0 + TANBT**2)
      COSB = 1.D0/DSQRT(1.D0 + TANBT**2)
      COS2B = (TANBT**2 - 1.D0)/(TANBT**2 + 1.D0)
      G1 = DSQRT(ALPHA1*4.D0*PI)
      G2 = DSQRT(ALPHA2*4.D0*PI)
      G3 = DSQRT(ALPHA3*4.D0*PI)
      HU = RMTOP/V/SINB
      HD =  MB/V/COSB
C

      IF(MQ.GT.MUR) TP = TQ - TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TP = TU - TQ
      IF(MQ.GT.MUR) TDP = TU
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) TDP = TQ
      IF(MQ.GT.MD) TPD = TQ - TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TPD = TD - TQ
      IF(MQ.GT.MD) TDPD = TD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) TDPD = TQ

      IF(MQ.GT.MD) DLAMBDA1 = 6./96./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA1 = 3./32./PI**2*
     * HD**2*(G1**2/3.+G2**2)*TPD

      IF(MQ.GT.MUR) DLAMBDA2 =12./96./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA2 = 3./32./PI**2*
     * HU**2*(-G1**2/3.+G2**2)*TP

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  dlambdap1 and dlambdap2 are the new log corrections due to
c  the presence of the gluino mass. They are in general very small,
c  and only present if there is a hierarchy of masses between the
c  two stops.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        dlambdap2 = 0
        tglu = log(mglu**2/mtop**2)

        if(mglu.lt.mur.or.mglu.lt.mq) then
        if(mq.gt.mur.and.mglu.gt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tglu**2)
        endif

        if(mq.gt.mur.and.mglu.lt.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
        endif

        if(mq.gt.mur.and.mglu.eq.mur) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tq**2-tu**2)
        endif

        if(mur.gt.mq.and.mglu.gt.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tglu**2)
        endif

        if(mur.gt.mq.and.mglu.lt.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
        endif

        if(mur.gt.mq.and.mglu.eq.mq) then
        dlambdap2 = -4./(16.*pi**2)**2*hu**4*(tu**2-tq**2)
        endif
        endif

      DLAMBDA3 = 0.
      DLAMBDA4 = 0.

      IF(MQ.GT.MD) DLAMBDA3 = -1./32./PI**2*G1**2*HD**2*TPD
      IF(MQ.LT.MD.OR.MQ.EQ.MD) DLAMBDA3 = 3./64./PI**2*HD**2*
     *(G2**2-G1**2/3.)*TPD
      
      IF(MQ.GT.MUR) DLAMBDA3 = DLAMBDA3 - 
     *1./16./PI**2*G1**2*HU**2*TP
      IF(MQ.LT.MUR.OR.MQ.EQ.MUR) DLAMBDA3 = DLAMBDA3 + 
     * 3./64./PI**2*HU**2*(G2**2+G1**2/3.)*TP

      IF(MQ.LT.MUR) DLAMBDA4 = -3./32./PI**2*G2**2*HU**2*TP
      IF(MQ.LT.MD) DLAMBDA4 = DLAMBDA4 - 3./32./PI**2*G2**2*
     *                        HD**2*TPD
C
      LAMBDA1 = ((G1**2 + G2**2)/4.)*
     *(1.-3.*HD**2*(TPD + TDPD)/8./PI**2)
     *+(3.*HD**4./16./PI**2) *TPD*(1.   
     *+ (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * (TPD + 2.*TDPD)/16./PI**2) 
     *+(3.*HD**4./8./PI**2) *TDPD*(1.  + (3.*HD**2/2. + HU**2/2.       
     *- 8.*G3**2) * TDPD/16./PI**2) + DLAMBDA1 
C
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*
     *(TP + TDP)/8./PI**2)
     *+(3.*HU**4./16./PI**2) *TP*(1.   
     *+ (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * (TP + 2.*TDP)/16./PI**2) 
     *+(3.*HU**4./8./PI**2) *TDP*(1. + (3.*HU**2/2. + HD**2/2.       
     *- 8.*G3**2) * TDP/16./PI**2) + DLAMBDA2  + DLAMBDAP2
C
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2)*(TP + TDP)/16./PI**2 -3.*
     *(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA3 
C
      LAMBDA4 = (- G2**2/2.)*(1.
     *-3.*(HU**2)*(TP + TDP)/16./PI**2
     *-3.*(HD**2)*(TPD + TDPD)/16./PI**2) +DLAMBDA4
C     
	LAMBDA5 = 0.
	LAMBDA6 = 0.
	LAMBDA7 = 0.

C
C     THIS IS THE CONTRIBUTION FROM LIGHT CHARGINOS/NEUTRALINOS
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  	 MSSUSY=DSQRT(0.5D0*(MQ**2+MUR**2)+MTOP**2)
	IF(MCHI.GT.MSSUSY)GOTO 3790
	IF(MCHI.LT.MTOP) MCHI=MTOP
	TCHAR=LOG(MSSUSY**2/MCHI**2)
	DELTAL12=(9./64./PI**2*G2**4+5./192./PI**2*G1**4)*TCHAR
	DELTAL3P4=(3./64./PI**2*G2**4+7./192./PI**2*G1**4
     *       +4./32/PI**2*G1**2*G2**2)*TCHAR
	DELTAM112=2.*DELTAL12*V**2*COSB**2
	DELTAM222=2.*DELTAL12*V**2*SINB**2
	DELTAM122=2.*DELTAL3P4*V**2*SINB*COSB
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
        DLAM1 = DELTAM112/2.D0/V**2/COSB**2
        DLAM2 = DELTAM222/2.D0/V**2/SINB**2
        DLAM3 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(G1**2-G2**2)/(G1**2+G2**2)
        DLAM4 = DELTAM122/2.D0/V**2/SINB/COSB
     .        *(2*G2**2)/(G1**2+G2**2)
        LAMBDA1 = LAMBDA1+DLAM1
        LAMBDA2 = LAMBDA2+DLAM2
        LAMBDA3 = LAMBDA3+DLAM3
        LAMBDA4 = LAMBDA4+DLAM4
C--END OF EXTENSION
 3790	CONTINUE
CCCCCCCCCCCCCCC    END OF CHARGINOS AND NEUTRALINOS  CCCCCCCCCCCC 


C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
      CALL GFUN_TLH(MA,TANBA,MQ,MUR,MD,MTOP,AU,AD,MU,MGLU,
     *                 DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7)

      LAMBDA1 = LAMBDA1+DLAM1
      LAMBDA2 = LAMBDA2+DLAM2
      LAMBDA3 = LAMBDA3+DLAM3
      LAMBDA4 = LAMBDA4+DLAM4
      LAMBDA5 = LAMBDA5+DLAM5
      LAMBDA6 = LAMBDA6+DLAM6
      LAMBDA7 = LAMBDA7+DLAM7
      
      M2(1,1) = 2.*V**2*(LAMBDA1*COSB**2+2.*LAMBDA6*
     *COSB*SINB + LAMBDA5*SINB**2) + MA**2*SINB**2
      M2(2,2) = 2.*V**2*(LAMBDA5*COSB**2+2.*LAMBDA7*
     *COSB*SINB + LAMBDA2*SINB**2) + MA**2*COSB**2
      M2(1,2) = 2.*V**2*(LAMBDA6*COSB**2+(LAMBDA3+LAMBDA4)*
     *COSB*SINB + LAMBDA7*SINB**2) - MA**2*SINB*COSB
      M2(2,1) = M2(1,2)

      M2P(1,1) = M2(1,1)
      M2P(2,2) = M2(2,2)
      M2P(1,2) = M2(1,2)
      M2P(2,1) = M2(2,1)

C--END OF EXTENSION

      TRM2P  = M2P(1,1) + M2P(2,2)
      DETM2P = M2P(1,1)*M2P(2,2) - M2P(1,2)*M2P(2,1)

      MH2P = (TRM2P - DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
      HM2P = (TRM2P + DSQRT(TRM2P**2 - 4.D0* DETM2P))/2.D0
C !!!!!!!!!!!!!!!!!!!
      MCH2=MA**2+(LAMBDA5-LAMBDA4)*V**2
C !!!!!!!!!!!!!!!!!!!
      MCH=DSQRT(MCH2)
      HMP = DSQRT(HM2P) 
      IF(MH2P.LT.0.)GOTO 5555
      MHP = DSQRT(MH2P) 
C
      SIN2ALPHA = 2.*M2P(1,2)/DSQRT(TRM2P**2-4.D0*DETM2P)
      COS2ALPHA = (M2P(1,1)-M2P(2,2))/DSQRT(TRM2P**2-4.D0*DETM2P)
      IF(COS2ALPHA.GT.0.) ALPHA = DASIN(SIN2ALPHA)/2.D0
      IF(COS2ALPHA.LT.0.) ALPHA = -PI/2.D0-DASIN(SIN2ALPHA)/2.D0
      SA = DSIN(ALPHA)
      CA = DCOS(ALPHA)  
      SQBMA = (SINB*CA - COSB*SA)**2

5555  RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCC NON DEGENERATE STOP/SBOTTOM EFFECTS CCCCCCCCC
C
        SUBROUTINE GFUN_TLH(MA,TANB,MQ,MUR,MD,MTOP,AT,AB,MU,MGLU,
     *                     DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        DIMENSION VH(2,2),VH1(2,2),VH2(2,2),
     *            VH3T(2,2),VH3B(2,2),AL(2,2)
        COMMON/su_PARAM/GF,ALPH,AMZ,AMW
        G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y)
        amtau = 1.777d0
        ammuon = 0.106d0
        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MQ.GT.MUR) MST = MQ
        IF(MUR.GT.MQ.OR.MUR.EQ.MQ) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MQ.GT.MD) MSB = MQ
	IF(MD.GT.MQ.OR.MD.EQ.MQ) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        RMTOP   = RUNM(MTOP,6)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        AL(1,1) = AL1
        AL(1,2) = (AL2+AL1)/2.
        AL(2,1) = (AL2+AL1)/2.
        AL(2,2) = AL2

	IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
	ELSE
	VI =  V
	H1I = VI*COSB
	H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
	END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

      CALL DELMB_TLH(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)

        IF(STOP22.LT.0.) GOTO 4237
        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
        MBOT2 = DSQRT(MBOT4)
        mb = mb/(1+deltamb)

        VH1(1,1) = 1./TANBST
        VH1(2,1) = -1.
        VH1(1,2) = -1.
        VH1(2,2) = TANBST
        VH2(1,1) = TANBST
        VH2(1,2) = -1.
        VH2(2,1) = -1.
        VH2(2,2) = 1./TANBST

C CCCCCCCCCCCCCCCCCCCCCCCCCCC  D-terms CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	STW=SW

	F1T=(MQ2-MUR2)/(STOP12-STOP22)*(.5-4./3.*STW)*
     *         LOG(STOP1/STOP2)
     *        +(.5-2./3.*STW)*LOG(STOP1*STOP2/(MQ2+MTOP2))
     *        + 2./3.*STW*LOG(STOP1*STOP2/(MUR2+MTOP2))

	F1B=(MQ2-MD2)/(SBOT12-SBOT22)*(-.5+2./3.*STW)*
     *        LOG(SBOT1/SBOT2)
     *        +(-.5+1./3.*STW)*LOG(SBOT1*SBOT2/(MQ2+MBOT2))
     *        - 1./3.*STW*LOG(SBOT1*SBOT2/(MD2+MBOT2))

	F2T=1/(STOP12-STOP22)*
     *         (-.5*LOG(STOP12/STOP22)
     *        +(4./3.*STW-.5)*(MQ2-MUR2)/(STOP12-STOP22)*
     *         G(STOP12,STOP22))

	F2B=1/(SBOT12-SBOT22)*
     *         (.5*LOG(SBOT12/SBOT22)
     *        +(-2./3.*STW+.5)*(MQ2-MD2)/(SBOT12-SBOT22)*
     *        G(SBOT12,SBOT22))

C*************************************************************
C
C--EXTENSION OF CARENA ET AL.: TRAFO MASS MATRIX -> LAMBDA_I
C
C TRAFOS APPROXIMATE -> EXACT:
C
C (i)  1/M_{SUSY}^2 -> LOG(M1^2/M2^2) / (M1^2-M2^2)
C
C (ii) 1/M_{SUSY}^4 -> -6 G(M1^2,M2^2) / (M1^2-M2^2)^2
C
C Then use results of Phys. Lett. B355 (1995) 209 in order to
C obtain the results for lambda_1 - lambda_7 according to
C Nucl. Phys. B461 (1996) 407. Perform a full evolution from
C M_SUSY -> m_t for lambdas (anomalous dimensions, v_i).
C
C - ht^2*hb^2 terms neglected in lambda_3,4 (according to
C   Nucl. Phys. B461 (1996) 407)
C
C*************************************************************

        DLAM1T = MTOP4/(SINBT**4)*(MU**2/(STOP1**2
     *    -STOP2**2))**2*G(STOP12,STOP22)
     *  - MZ**2*MTOP2*MU**2/TANBST**2*F2T/COSBT**2

        DLAM1B = MBOT4/(COSBB**4)*(LOG(SBOT1**2*SBOT2**2/
     *    (MQ2+MBOT2)/(MD2+MBOT2))
     *    + 2*AB**2/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2))
     *  + MBOT4/(COSBB**4)*(AB**2/
     *    (SBOT1**2-SBOT2**2))**2*G(SBOT12,SBOT22)
     *  + MZ**2*(2*MBOT2*F1B-MBOT2*AB**2*F2B)/COSBB**2

        DLAM2T = MTOP4/(SINBT**4)*(LOG(STOP1**2*STOP2**2/
     *    (MQ2+MTOP2)/(MUR2+MTOP2))
     *  + 2*AT**2/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2))
     *  + MTOP4/(SINBT**4)*(AT**2/
     *    (STOP1**2-STOP2**2))**2*G(STOP12,STOP22)
     *  + MZ**2*(-2*MTOP2*F1T+MTOP2*AT**2*F2T)/SINBT**2
 
        DLAM2B = MBOT4/(COSBB**4)*MU**4/(SBOT1**2
     *    -SBOT2**2)**2*G(SBOT12,SBOT22)
     *    + MZ**2*MBOT2*MU**2*TANBSB**2*F2B/SINBB**2
 
        DLAM3T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  + MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM3B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM4T = MTOP4/(SINBT**4)*
     *    MU**2/(STOP1**2-STOP2**2)*(LOG(STOP1**2/STOP2**2)/2.D0
     *  + AT**2/(STOP1**2-STOP2**2)*G(STOP12,STOP22))
     *  + MZ**2*(MTOP2/TANBST*F1T-MTOP2*(AT**2-MU**2)/TANBST/2.*F2T)
     *    /SINBT/COSBT/2
c    *  - MTOP2*MBOT2/(SINBT**2*COSBB**2)*(
c    *    LOG(STOP1**2*STOP2**2/(MQ2+MTOP2)/(MUR2+MTOP2))
c    *  + LOG(SBOT1**2*SBOT2**2/(MQ2+MBOT2)/(MD2+MBOT2))
c    *  + ((AT+AB)**2/2-MU**2)*(
c    *      1.D0/(STOP1**2-SBOT1**2)*LOG(STOP1**2/SBOT1**2)
c    *    + 1.D0/(STOP2**2-SBOT2**2)*LOG(STOP2**2/SBOT2**2))
c    *  - (MU**2-AT*AB)**2*(
c    *    - 1.D0/(STOP1**2-SBOT1**2)**2*G(STOP12,SBOT12)
c    *    - 1.D0/(STOP2**2-SBOT2**2)**2*G(STOP22,SBOT22)))

        DLAM4B = MBOT4/(COSBB**4)*MU**2/(SBOT1**2-SBOT2**2)*(
     *    LOG(SBOT1**2/SBOT2**2)/2.D0
     *  + AB**2/(SBOT1**2-SBOT2**2)*G(SBOT12,SBOT22))
     *  + MZ**2*(-MBOT2*TANBSB*F1B+MBOT2*(AB**2-MU**2)*TANBSB/2.*F2B)
     *    /SINBB/COSBB/2

        DLAM5T = MTOP4/(SINBT**4)*
     *    (MU**2*AT**2)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)

        DLAM5B = MBOT4/(COSBB**4)*
     *    (MU**2*AB**2)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)

        DLAM6T = MTOP4/(SINBT**4)*
     *    (-MU**3*AT)/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22)
     *  + MZ**2*MTOP2*MU*AT/TANBST*F2T/(2*SINBT*COSBT)

        DLAM6B = MBOT4/(COSBB**4)*MU*AB*
     *    (-1.D0/(SBOT1**2-SBOT2**2)*LOG(SBOT1**2/SBOT2**2)
     *    -AB**2/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22))
     *  - MZ**2*(-MBOT2*AB*MU*TANBSB*F2B)/(2*SINBB*COSBB)

        DLAM7T = MTOP4/(SINBT**4)*MU*AT*
     *    (-1.D0/(STOP1**2-STOP2**2)*LOG(STOP1**2/STOP2**2)
     *    -AT**2/(STOP1**2-STOP2**2)**2*G(STOP12,STOP22))
     *  - MZ**2*MTOP2*AT*MU/TANBST*F2T/(2*SINBT*COSBT)

        DLAM7B = MBOT4/(COSBB**4)*
     *    (-MU**3*AB)/(SBOT1**2-SBOT2**2)**2*G(SBOT12,SBOT22)
     *    - MZ**2*MBOT2*MU*AB*TANBSB*F2B/(2*SINBB*COSBB)

       TQ = LOG((MQ2 + MTOP2)/MTOP2)
       TU = LOG((MUR2+MTOP2)/MTOP2)
       TQD = LOG((MQ2 + MB**2)/MB**2)
       TD = LOG((MD2+MB**2)/MB**2)

        FACT = 3.D0/(16.D0*PI**2*(H1T**2+H2T**2)**2)
        FACB = 3.D0/(16.D0*PI**2*(H1B**2+H2B**2)**2)

        DLAM1 = FACT*DLAM1T*(1.-AL1*TT) + FACB*DLAM1B*(1.-AL1*TB)

        DLAM2 = FACT*DLAM2T*(1.-AL2*TT) + FACB*DLAM2B*(1.-AL2*TB)

        DLAM3 = FACT*DLAM3T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM3B*(1.-(AL1+AL2)/2*TB)

        DLAM4 = FACT*DLAM4T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM4B*(1.-(AL1+AL2)/2*TB)

        DLAM5 = FACT*DLAM5T*(1.-(AL1+AL2)/2*TT)
     *        + FACB*DLAM5B*(1.-(AL1+AL2)/2*TB)

        DLAM6 = FACT*DLAM6T*(1.-(3*AL1+AL2)/4*TT)
     *        + FACB*DLAM6B*(1.-(3*AL1+AL2)/4*TB)

        DLAM7 = FACT*DLAM7T*(1.-(AL1+3*AL2)/4*TT)
     *        + FACB*DLAM7B*(1.-(AL1+3*AL2)/4*TB)

        FACTOR = 1.D0
        DLAM1 = DLAM1 * FACTOR
        DLAM2 = DLAM2 * FACTOR
        DLAM3 = DLAM3 * FACTOR
        DLAM4 = DLAM4 * FACTOR
        DLAM5 = DLAM5 * FACTOR
        DLAM6 = DLAM6 * FACTOR
        DLAM7 = DLAM7 * FACTOR

C--END OF EXTENSION

        GOTO 4236
 4237   CONTINUE

        DLAM1 = -1.D+15
        DLAM2 = -1.D+15
        DLAM3 = -1.D+15
        DLAM4 = -1.D+15
        DLAM5 = -1.D+15
        DLAM6 = -1.D+15
        DLAM7 = -1.D+15

4236    RETURN
        END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       End of program from M. Carena, M. Quiros and C.E.M. Wagner.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION T_TLH(X,Y,Z)
      implicit real*8(a-h,l,m,o-z)
      if(x.eq.y) x = x - 0.00001
      if(x.eq.z) x = x - 0.00002
      if(y.eq.z) y = y - 0.00003
c       write(*,*) 'xyz',x,y,z
      T_TLH = (X**2*Y**2*log(X**2/Y**2) + X**2*Z**2*log(Z**2/X**2)
     * + Y**2*Z**2*log(Y**2/Z**2))/((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      return
      end

      SUBROUTINE DELMB_TLH(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        COMMON/su_PARAM/GF,ALPH,AMZ,AMW
        amtau=1.77d0
        ammuon=0.106d0
        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      RMTOP = RUNM(MTOP,6)
      MB = RUNM(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS(AMZ,2)
      ALPHA3  = ALPHAS(MTOP,2)

      G1 = DSQRT(ALPHA1*4.*PI)
      G2 = DSQRT(ALPHA2*4.*PI)
      G3 = DSQRT(ALPHA3*4.*PI)
      
        IF(MQ.GT.MUR) MST = MQ
        IF(MUR.GT.MQ.OR.MUR.EQ.MQ) MST = MUR
        MSUSYT = DSQRT(MST**2  + MTOP**2)

	IF(MQ.GT.MD) MSB = MQ
	IF(MD.GT.MQ.OR.MD.EQ.MQ) MSB = MD
	MSUSYB = DSQRT(MSB**2 + MB**2)

	TT = LOG(MSUSYT**2/MTOP**2)
	TB = LOG(MSUSYB**2/MTOP**2)

        HT = RMTOP/V/SINB
        HTST = RMTOP/V
        HB =  MB/V/COSB
        G32 = ALPHA3*4.*PI

        BT2 = -(8.*G32 - 9.*HT**2/2. - HB**2/2.)/(4.*PI)**2
	BB2 = -(8.*G32 - 9.*HB**2/2. - HT**2/2.)/(4.*PI)**2
        AL2 = 3./8./PI**2*HT**2
        BT2ST = -(8.*G32 - 9.*HTST**2/2.)/(4.*PI)**2
        ALST = 3./8./PI**2*HTST**2
        AL1 = 3./8./PI**2*HB**2

        IF(MA.GT.MTOP) THEN
        VI = V*(1. + 3./32./PI**2*HTST**2*LOG(MTOP**2/MA**2))
        H1I = VI*COSBA
        H2I = VI*SINBA
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MA**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MA**2/MSUSYB**2))**.25
        ELSE
        VI =  V
        H1I = VI*COSB
        H2I = VI*SINB
        H1T = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYT**2))**.25
        H2T = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYT**2))**.25
        H1B = H1I*(1.+3./8./PI**2*HB**2*LOG(MTOP**2/MSUSYB**2))**.25
        H2B = H2I*(1.+3./8./PI**2*HT**2*LOG(MTOP**2/MSUSYB**2))**.25
        END IF

        TANBST = H2T/H1T
        SINBT = TANBST/(1.+TANBST**2)**.5
        COSBT = SINBT/TANBST

        TANBSB = H2B/H1B
        SINBB = TANBSB/(1.+TANBSB**2)**.5
        COSBB = SINBB/TANBSB

        deltamt = 0
        deltamb = 0

        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MQ2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MQ2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MQ2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MQ2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        IF(SBOT22.LT.0.) GOTO 4237

        STOP1 = STOP12**.5
        STOP2 = STOP22**.5
        SBOT1 = SBOT12**.5
        SBOT2 = SBOT22**.5

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Here is the definition of deltamb and deltamt, which
c     are the vertex corrections to the bottom and top quark
c     mass, keeping the dominant QCD and top Yukawa coupling
c     induced corrections.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        deltamb = -2*alpha3/3./pi*mglu*(ab-mu*tanb)*
     *  T_TLH(sbot1,sbot2,mglu)
     *  + ht**2/(4.*pi)**2*(at-mu/tanb)*mu*tanb*
     *  T_TLH(stop1,stop2,mu)


        deltamt = -2.*alpha3/3./pi*(at-mu/tanb)*mglu*
     *  T_TLH(stop1,stop2,mglu)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Here the new values of the top and bottom quark masses at
c   the scale MS are defined, to be used in the effective
c   potential approximation. They are just the old ones, but
c   including the finite corrections deltamt and deltamb.
c   The deltamb corrections can become large and are resummed
c   to all orders, as suggested in the two recent works by M. Carena,
c   S. Mrenna and C.E.M. Wagner, as well as in the work by M. Carena,
c   D. Garcia, U. Nierste and C.E.M. Wagner, to appear. The top
c   quark mass corrections are small and are kept in the perturbative
c   formulation. The function T(X,Y,Z) is necessary for the calculation.
c   the entries are masses and NOT their squares !
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        mtop4 = rmtop**4.*(1.+2.*bt2*tt- al2*tt - 4.*deltamt)
c     * /(1.+deltamt)**4.
        mbot4 = mb**4.*(1.+2.*bb2*tb - al1*tb)
     * /(1.+deltamb)**4.
        MTOP2 = DSQRT(MTOP4)
	MBOT2 = DSQRT(MBOT4)

        STOP12 = (MQ2 + MUR2)*.5 + MTOP2 
     *   +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2)
     *   +(((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *   MQ2 - MUR2)**2*0.25 + MTOP2*(AT-MU/TANBST)**2)**.5

        STOP22 = (MQ2 + MUR2)*.5 + MTOP2 
     *  +1./8.*(G2**2+G1**2)*(H1T**2-H2T**2) 
     *   - (((G2**2-5.*G1**2/3.)/4.*(H1T**2-H2T**2) +
     *  MQ2 - MUR2)**2*0.25 
     *  + MTOP2*(AT-MU/TANBST)**2)**.5

        IF(STOP22.LT.0.) GOTO 4237

        SBOT12 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *  + (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *  MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

        SBOT22 = (MQ2 + MD2)*.5  
     *   - 1./8.*(G2**2+G1**2)*(H1B**2-H2B**2)
     *   - (((G1**2/3.-G2**2)/4.*(H1B**2-H2B**2) +
     *   MQ2 - MD2)**2*0.25 + MBOT2*(AB-MU*TANBSB)**2)**.5

4237    RETURN
        END
c---------------------------------------------------------
      COMPLEX*16 FUNCTION F0_TLH(M1,M2,QSQ)
      IMPLICIT REAL*8 (A-H,M,O-Z)
      COMPLEX*16 CD,CR,CQ2,IEPS,CBET,CXX
      M1SQ = M1*M1
      M2SQ = M2*M2
      AQSQ = DABS(QSQ)
      IEPS = DCMPLX(1.D0,1.D-12)
      CQ2 = QSQ*IEPS
      CD = (M1SQ-M2SQ)/CQ2
      CR = CDSQRT((1+CD)**2 - 4*M1SQ/CQ2)
      IF(QSQ.EQ.0.D0) THEN
       F0_TLH = 0.D0
      ELSE
       IF(M1.EQ.M2) THEN
        F0_TLH = -2.D0 + CR*CDLOG(-(1+CR)/(1-CR))
       ELSE
        CBET = CDSQRT(1-4*M1*M2/(CQ2 - (M1-M2)**2))
        CXX = (CBET-1)/(CBET+1)
        F0_TLH = -1 + ((QSQ+M2SQ-M1SQ)/2/QSQ - M2SQ/(M2SQ-M1SQ))
     .                                           *DLOG(M2SQ/M1SQ)
     .     - (QSQ-(M1-M2)**2)/QSQ*CBET*CDLOG(CXX)
       ENDIF
      ENDIF
      RETURN
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c    End of two loop corrections to the neutral Higgs boson masses 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%