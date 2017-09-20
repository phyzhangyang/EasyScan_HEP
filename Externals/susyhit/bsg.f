CC ----------------------------------------------------------------------------
CC  !last changes: J-L Kneur, March 01, 2005 
C ------------------------------------------------------------------------------
C New version of the bsg routine (January 10, 2001) 
C Many new things ...
C Original file name from P.Gambino : susyNLO.f
C ---------------- Beginning of code --------------------------------------------
C       La routine matching calcola i coeff di Wilson per un dato modello.
C       La routine bsg calcola il branching ratio a partire dai 
C       Coeff di Wilson. version 1/6/98 bug corretto
C      our ew corrections are implemented only for mh=100 7/2000

      subroutine su_bsg(as,t,bc,rcb,aeinv,scw,scb,scl,vkm,bsl,deltp,io,
     +   c70,c71,c80,c81,ee,bbox,br)
C       Input:
C       as          alpha_s (M_Z)
C       t           top mass [enters only in the calculation of alpha_s(scw)]
C       bc          pole bottom mass - pole charm mass
C       rcb         pole charm mass / pole bottom mass
C       aeinv       1/alpha(muw)
C       scw         matching scale
C       scb         low-energy scale
C       scl         scale for semileptonic
C       vkm         |V_ts V_tb / V_cb |^2
C       bsl         semileptonic BR
C       deltp       cut on photon energy
C       io          LO (io=0) NLO (io=1)
C       c70,c71,c80,
C       c81,ee      matching conditions
C       bbox        ratio box/box_SM for CKM matrix elemts
C       Output:
C       br          BR(b->s gamma) 
        implicit double precision (a-h,o-z)
        dimension aa(8),eee(8),ff(8),gg(8),hh(8),clead(8),
     . fnum(8,8),bb(8)
        common/alp/asc,tc,pi,zm,ioc
        common/fun/z
        common/su_bsllsign/iwrongsign
        external func1,func2,func3,func4

C       PARAMETRI DI INPUT
        wm=80.419d0
        zm=91.187d0
        rbs=50.d0
        delt=deltp
        hlam=0.12d0    ! lambda_2 (HQET parameter)
        ityp=1    ! calcola f_ij tramite integrali
c       ityp=2    ! valori numerici di f_ij per delt=0.99 o 0.9 e sqrt(z)=0.29

c EW corrections for mh=100 ??
        iew=1

C       ALTRI PARAMETRI
        asc=as
        ioc=io
        tc=t
        z=rcb**2
        ae=1.d0/aeinv
        pi=4.d0*datan(1.d0)
        pi2=pi**2

        if(c70.eq.0) then
           br=0
           return
        endif

C       MAGIC NUMBERS
        aa(1)=14./23.
        aa(2)=16./23.
        aa(3)=6./23.
        aa(4)=-12./23.
        aa(5)=0.4086
        aa(6)=-0.4230
        aa(7)=-0.8994
        aa(8)=0.1456
        eee(1)=4661194./816831.
        eee(2)=-8516./2217.
        eee(3)=0.
        eee(4)=0.
        eee(5)=-1.9043
        eee(6)=-0.1008
        eee(7)=0.1216
        eee(8)=0.0183
        ff(1)=-17.3023
        ff(2)=8.5027
        ff(3)=4.5508
        ff(4)=0.7519
        ff(5)=2.0040
        ff(6)=0.7476
        ff(7)=-0.5385
        ff(8)=0.0914
        gg(1)=14.8088
        gg(2)=-10.8090
        gg(3)=-0.8740
        gg(4)=0.4218
        gg(5)=-2.9347
        gg(6)=0.3971
        gg(7)=0.1600
        gg(8)=0.0225
        hh(1)=626126./272277.
        hh(2)=-56281./51730.
        hh(3)=-3./7.
        hh(4)=-1./14.
        hh(5)=-0.6494
        hh(6)=-0.0380
        hh(7)=-0.0186
        hh(8)=-0.0057
        bb(1)=0.5784
        bb(2)=-0.3921
        bb(3)=-0.1429
        bb(4)=0.04762
        bb(5)=-0.1275
        bb(6)=0.03174
        bb(7)=0.007811
        bb(8)=-0.003101

C       CALCOLO ALPHA
        asscw=asf(scw)
        asscb=asf(scb)
        eta=asscw/asscb

C       COEFF WILSON A SCB
        c10b=-eta**(-12./23.)+eta**(6./23.)
        c20b=(eta**(-12./23.)+2*eta**(6./23.))/3.
        sumh=0.
        do i=1,8
         sumh=sumh+hh(i)*eta**aa(i)
        enddo


        c70b=eta**(16./23.)*c70+8./3.*(eta**(14./23.)-eta**(16./23.))*
     +   c80+sumh
c
        if(c70b .lt. 0.d0) then
          iwrongsign = 0       !Sign as in SM
        else
          iwrongsign = 1       !Sign has been flipped: trouble w/ b --> sll!
        endif

c addition of sign of C7 = sign of SM constraint (equiv b->s l+l- sign):
c electroweak corrections only for the calculation of c70b
        if(iew.eq.1)then
           c70t= c70 *(0.974d0) !-0.005d0)
           c80t= c80 *(0.993d0) ! -0.005d0)
c    c70(mb) including ew corrections; last term from ew from mixings
      c70bew=eta**(16./23.)*c70t+8./3.*(eta**(14./23.)-eta**(16./23.))*
     +   c80t+sumh + 0.00054d0
        else
        c70bew=c70b
        endif
c       print*,"c70bew",c70bew,c70b
        
        c80b=(c80+313063./363036.)*eta**(14./23.)-0.9135*eta**0.4086+
     +   0.0873*eta**(-0.4230)-0.0571*eta**(-0.8994)+0.0209*eta**0.1456

        if(io.eq.1)then
         sumg=0.
         do i=1,8
c old     sumg=sumg+(eee(i)*eta*ee+
c     +   ff(i)+gg(i)*eta)*eta**aa(i)
         sumg=sumg+(eee(i)*eta*(ee+4.d0/3.d0*Log(scw/wm))+
     +    ff(i)+(gg(i)+12*bb(i)*Log(scw/wm))*eta)*eta**aa(i)
         enddo
         c71b=eta**(39./23)*c71+8./3*(eta**(37./23)-eta**(39./23))*c81
         c71b=c71b+(297664./14283.*eta**(16./23.)-7164416./357075.*
     +    eta**(14./23.)+256868./14283.*eta**(37./23.)-6698884./357075.*
     +    eta**(39./23.))*c80
        c71b=c71b+37208./4761.*(eta**(39./23.)-eta**(16./23.))*c70+sumg
        endif
        clead(1)=c10b
        clead(2)=c20b
        clead(3)=-eta**(-12./23.)/27.+2*eta**(6./23.)/63.-0.0659*
     +   eta**0.4086+0.0595*eta**(-0.423)-0.0218*eta**(-0.8994)+
     +   0.0335*eta**0.1456
        clead(4)=-eta**(-12./23.)/9.+eta**(6./23.)/21.+0.0237*
     +   eta**0.4086-0.0173*eta**(-0.423)-0.1336*eta**(-0.8994)-
     +   0.0316*eta**0.1456
        clead(5)=eta**(-12./23.)/108.-eta**(6./23.)/126.+0.0094*
     +   eta**0.4086-0.01*eta**(-0.423)+0.0010*eta**(-0.8994)-
     +   0.0017*eta**0.1456
        clead(6)=-eta**(-12./23.)/36.-eta**(6./23.)/84.+0.0108*
     +   eta**0.4086+0.0163*eta**(-0.423)+0.0103*eta**(-0.8994)+
     +   0.0023*eta**0.1456
        clead(7)=c70b
        clead(8)=c80b

C       CORREZIONI HQET
        bm2=(bc/(1-sqrt(z)))**2
        bm=sqrt(bm2)   ! bottom mass
        cm2=z*bm2      ! charm squared mass
        fkin=1-8*z+8*z**3-z**4-12*z**2*log(z)
        hqb=-6*hlam*(1-(1-z)**4/fkin)
        hqc=  hlam*(c10b/54.d0-c20b/9.d0)/c70b
                        !corrected according to NPB version 
        hqet=1+hqb/bm2+hqc/cm2

c       print*,c71b
C       CALCOLO LEADING 
c     vkm from input is the result of a SM fit here we change to other models
        box= 1.d0+(1.04884-vkm)*(1.d0-1.d0/dsqrt(bbox))/vkm
c        box=1.d0   !in case we don't want that 

c     Czarnecki-Marciano QED correction
        alpha0= 1.d0/137.036
        cqed0=alpha0/ae*(1- 2.d0*ae*log(zm/bm)/pi)*
     -   (1 -ae/pi*(104.d0/243.d0 - 8.d0*c70/9.d0)*log(wm/bm)/abs(c70b))

        br=bsl*vkm*6*ae*c70b**2/(pi*fkin)*hqet*box
        
c        brsenza=bsl*vkm*6*ae/(pi*fkin)*hqet*box
        if(io.eq.0)then
           br= br*cqed0
c          print*, cqed0
        endif

C       CORREZIONE NLO
        if(io.eq.1)then
        cb=-8.*asf(scb)/(3*pi) !-18.6d0*(asf(scb)/pi)**2
c dalla massa polo del b

        z2=z*z
        z3=z2*z
        z4=z3*z
        zs=sqrt(z)
        zl=log(z)
C       **** ho portato 17./ nella riga sotto
        hz=-(1-z2)*(25./4.-239./3.*z+25./4.*z2)+z*zl*(20+90*z-4./3.*z2+
     +   17./3.*z3)+z2*zl**2*(36+z2)+(1-z2)*(17-64*z+17*z2)/3.*log(1-z)
        hz=hz-4*(1+30*z2+z4)*zl*log(1-z)-(1+16*z2+z4)*(6*sp(z)-pi2)-32*
     +   zs**3*(1+z)*(pi2-4*sp(zs)+4*sp(-zs)-2*zl*log((1-zs)/(1+zs)))
        csl=2*asf(scl)/(3*pi)*hz/fkin        ! dal semileptonico

        rie3=1.202057
        r2=-833+144*pi2*zs**3+(1728-180*pi2-1296*rie3+(1296-324*pi2)*zl+
     +   108*zl**2+36*zl**3)*z+(648+72*pi2+(432-216*pi2)*zl+36*zl**3)*z2
        r2=2./243.*(r2+(-54-84*pi2+1092*zl-756*zl**2)*z3)
        r1=-r2/6.
        r7=-10./3.-8*pi2/9.
        r8=-4./27.*(-33+2*pi2)
        geff1=-208./243.
        geff2=416./81.
        geff7=32./3.
        geff8=-32./9.
        cd=   c10b*(r1+geff1*log(bm/scb))+c20b*(r2+geff2*log(bm/scb))
        cd=cd+c70b*(r7+geff7*log(bm/scb))+c80b*(r8+geff8*log(bm/scb))
c       print*,c71b,cd
        cd=2.d0*asscb/(4*pi)*c70b*(c71b+cd)      ! dall'ampiezza D

        deltg=log(delt)
        if(ityp.eq.1)then
        delt2=delt**2
        delt3=delt2*delt
        deltg1=log(1-delt)
        fff77=(10*delt+delt2-2./3.*delt3+delt*(delt-4)*deltg)/3.
c       old formula
c       fff88=-2*log(rbs)*(delt2+2*delt+2*deltg1)+2*sp(1-delt)-pi2/3.
c       fff88=(fff88-delt*(1+2*delt)*deltg+8*deltg1-2./3.*delt3+3*delt2
c     +   +7*delt)/27.
c       fff78=8./9.*(sp(1-delt)-pi2/6.-delt*deltg+(11*delt-3*delt2+
c     +   delt3/3.)/4.)
c       new and correct formula
        fff88=-2*log(rbs)*(delt2+2*delt+4*deltg1)+4*sp(1-delt)-2*pi2/3.
        fff88=(fff88-delt*(2+delt)*deltg+8*deltg1-2./3.*delt3+3*delt2
     +   +7*delt)/27.
        fff78=8./9.*(sp(1-delt)-pi2/6.-delt*deltg+(9*delt-delt2+
     +   delt3/3.)/4.)
        eps=1.d-4
        est0=0.d0
        est1=(1-delt)/z
        est2=1./z
        fff22=16*z/27.*(delt*gauss1(func1,est0,est1,eps)+
     +   gauss1(func2,est1,est2,eps))   
        fff27=-8*z2/9.*(delt*gauss1(func3,est0,est1,eps)+
     +   gauss1(func4,est1,est2,eps))
        else if(ityp.eq.2)then
           if(abs(delt-0.99d0).lt.0.0001d0)then
        fff77=3.42106099217725
c old numbers
c       fff78=0.394065893946580
c       fff88=0.668299012321822
c new numbers (correct) for delta=0.99
        fff78=0.389665893956536
        fff88=3.216166804631025
        fff22=3.402419197348582d-002
        fff27=1.817503700368764d-002
        else if(abs(delt-0.9d0).lt.0.0001d0)then
        fff77=3.20598527956
        fff78=0.38734061186
        fff88=1.31742266135
        fff22=3.40367951120922d-002
        fff27=1.81788563733362d-002
        endif
        endif
        fff28=-fff27/3.
        fff11=fff22/36.
        fff12=-fff22/3.
        fff17=-fff27/6.
        fff18=-fff28/6.
        fnum(1,1)=fff11
        fnum(1,2)=fff12
        fnum(1,3)=-0.003493067823644294
        fnum(1,4)=0.0005821779706073824
        fnum(1,5)= -0.04589256102467917
        fnum(1,6)=-0.05997538991269699
        fnum(1,7)=fff17
        fnum(1,8)=fff18
        fnum(2,2)=fff22
        fnum(2,3)=0.02095840694186576
        fnum(2,4)=-0.003493067823644294
        fnum(2,5)=0.2753553661480751
        fnum(2,6)=0.359852339476182
        fnum(2,7)=fff27
        fnum(2,8)=fff28
        fnum(3,3)=0.01399990588537613
        fnum(3,4)=-0.004666635295125378
        fnum(3,5)=0.3276651145509532
        fnum(3,6)=0.06660649932577872
        fnum(3,7)=0.04208256917382194
        fnum(3,8)=-0.01402752305794064
        fnum(4,4)=0.008791470008620995
        fnum(4,5)=-0.05461085242515887
        fnum(4,6)=0.1569505914595811
        fnum(4,7)=-0.007013761528970321
        fnum(4,8)=0.002337920509656775
        fnum(5,5)=1.93694064002221
        fnum(5,6)=0.9505876392301583
        fnum(5,7)=0.5925919999999998
        fnum(5,8)=-0.1975306666666666
        fnum(6,6)=0.9305162079553912
        fnum(6,7)=0.0001919555740097084
        fnum(6,8)=-0.00006398519133657409
        fnum(7,7)=fff77
        fnum(7,8)=fff78
        fnum(8,8)=fff88

        a1=(exp(-asscb*deltg*(7+2*deltg)/(3*pi))-1)*c70b**2
        a2=0.d0
        do j=1,8
        do i=1,j
        a2=a2+fnum(i,j)*clead(i)*clead(j)
        enddo
        enddo
        a2=asscb/pi*a2
        ca=a1+a2                     !dal bremsstrahlung


c QED correction after eq 11 of Kagan-Neubert
        c7em= (32d0/75d0*eta**(-9d0/23d0) -40d0/69d0*eta**(-7d0/23d0)+    
     &   88./575.*eta**(16d0/23d0))*c70+ 
     &   (-32./575.*eta**(-9d0/23d0)+32./1449.*eta**(-7d0/23d0)+
     &   640./1449.*eta**(14.d0/23.d0)-704./1725.*eta**(16d0/23d0))*c80
     &   -190./8073.*eta**(-35d0/23d0)-359./3105.*eta**(-17./23.)+
     &   4276./121095.*eta**(-12./23.)+350531./1009125.*eta**(-9./23.)+
     &   2./4347.*eta**(-7./23.)-5956/15525.*eta**(6./23.)+
     &   38380./169533.*eta**(14./23.)-748./8625.*eta**(16./23.)

c here we can decide where to put the QED correction: using cqed or in the 
c overall QED factor qed1 below
        cqed=2.d0*ae/asscb*c7em*c70b
c       TRUNCATED VERSION
c       deff=sqrt(abs((1+cb+csl+(cd+ca+cqed)/c70bew**2)*c70bew**2))

c        NB here the term |D|**2 in eq.16 was expanded in powers of alphas
c        this may cause problems when c70b is small therefore here we put 
c        also the option without truncation
C       NON TRUNCATED VERSION
      deff=sqrt(abs((c70bew+cd/2.d0/c70b)**2+ca+cqed+(cb+csl)*c70b**2))

c - for the same reason we need to recalculate hqc
        hqc=  hlam*(c10b/54.d0-c20b/9.d0)*c70b/deff**2 
                        !corrected according to NPB version 
        hqetnlo=1+hqb/bm2+hqc/cm2
        
        br=br*deff**2/c70b**2*hqetnlo/hqet

c semileptonic qed corr
        cqed1= (1- 2.d0*ae*log(zm/bm)/pi) ! -
c     &  2.d0*ae/asscb*c7em/deff)
c     -    (1-0.55*ae/pi/deff*(104.d0/243.d0 -8.d0*c70/9.d0)*log(wm/bm))

        br= br*alpha0/ae*cqed1
        endif
        return
        end
        
C***************************************************************

        double precision function gre(t)
C       PARTE REALE DI G(T)
        implicit double precision (a-h,o-z)
        if(t.lt.4.d0)then
        gre=-2*(datan(sqrt(t/(4-t))))**2
        else
        pi=4.*datan(1.d0)
        pi2=pi**2
        gre=-pi2/2+2*(log((sqrt(t)+sqrt(t-4))/2.))**2
        endif
        return
        end

C*************************************************************************

        double precision function gim(t)
C       PARTE IMMAGINARIA DI G(T)
        implicit double precision (a-h,o-z)
        if(t.lt.4.d0)then
        gim=0.
        else
        pi=4.*datan(1.d0)
        gim=-2*pi*(log((sqrt(t)+sqrt(t-4))/2.))
        endif
        return
        end

C*************************************************************************

        double precision function func1(t)
        implicit double precision (a-h,o-z)
        common/fun/z
        func1=(1-z*t)*((gre(t)/t+.5)**2+(gim(t)/t)**2)
        return
        end

C*************************************************************************

        double precision function func2(t)
        implicit double precision (a-h,o-z)
        common/fun/z
        func2=(1-z*t)**2*((gre(t)/t+.5)**2+(gim(t)/t)**2)
        return
        end

C*************************************************************************

        double precision function func3(t)
        implicit double precision (a-h,o-z)
        common/fun/z
        func3=gre(t)+t/2.
        return
        end

C*************************************************************************

        double precision function func4(t)
        implicit double precision (a-h,o-z)
        common/fun/z
        func4=(1-z*t)*(gre(t)+t/2.)
        return
        end

C*************************************************************************

      double precision FUNCTION GAUSS1(F,A,B,EPS)
        implicit double precision (a-h,o-z)
                                                 
      REAL*8 W(12),X(12),A,B,EPS,DELTA,CONST,AA,BB,Y,C1,C2,S8,
     1                 S16,U,F
                              
      DATA CONST /1.0E-25/
                          
      DATA W
     1       / 0.10122 85362 90376 25915 25313 543,
     2         0.22238 10344 53374 47054 43559 944,
     3         0.31370 66458 77887 28733 79622 020,
     4         0.36268 37833 78361 98296 51504 493,
     5         0.02715 24594 11754 09485 17805 725,
     6         0.06225 35239 38647 89286 28438 370,
     7         0.09515 85116 82492 78480 99251 076,
     8         0.12462 89712 55533 87205 24762 822,
     9         0.14959 59888 16576 73208 15017 305,
     A         0.16915 65193 95002 53818 93120 790,
     B         0.18260 34150 44923 58886 67636 680,
     C         0.18945 06104 55068 49628 53967 232 /
                                                      
      DATA X
     1       / 0.96028 98564 97536 23168 35608 686,
     2         0.79666 64774 13626 73959 15539 365,
     3         0.52553 24099 16328 98581 77390 492,
     4         0.18343 46424 95649 80493 94761 424,
     5         0.98940 09349 91649 93259 61541 735,
     6         0.94457 50230 73232 57607 79884 155,
     7         0.86563 12023 87831 74388 04678 977,
     8         0.75540 44083 55003 03389 51011 948,
     9         0.61787 62444 02643 74844 66717 640,
     A         0.45801 67776 57227 38634 24194 430,
     B         0.28160 35507 79258 91323 04605 015,
     C         0.09501 25098 37637 44018 53193 354 /
                                                      
      DELTA=CONST*ABS(A-B)
      GAUSS1=0.
      AA=A
    5 Y=B-AA
      IF(ABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5*(AA+BB)
      C2=C1-AA
      S8=0.
      S16=0.
      DO 1 I = 1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I = 5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(ABS(S16-S8) .GT. EPS*(1.0+ABS(S16))) GO TO 4
      GAUSS1=GAUSS1+S16
      AA=BB
      GO TO 5
    4 Y=0.5*Y
      IF(ABS(Y) .GT. DELTA) GO TO 2
                                    
      PRINT 7
      GAUSS1=0.
      RETURN
    7 FORMAT(1X,37HGAUSS1 ... TOO HIGH ACCURACY REQUIRED)
      END                                                

C       *****************************************************************

        subroutine matching(imod,io,nlosusy,ihv,scw,as,t,h,tanb,
     +   stopl,stoph,xstop,sb1,sb2,xsbot,sqk,gl,aat,abo,amu,
     +   chm,u,v,c70,c80,c71,c81,ee,bbox,ierr)
C       Input:
C       imod        SM (0), Two Higgs (1), SUSY (2)
C       io          LO (io=0) NLO (io=1)
C       nlosusy     susy coeff at LO (=0) or NLO (=1)
C       ihv         if nlosusy=1: 1= heavy sqk/gluino case (new paper), 
C                                 0= light stop_R case (old paper)
C       scw         matching scale
C       as          alpha_s (M_Z)
C       t           pole top mass
C       h           Higgs mass 
C       tanb        tan(beta)
C       stopl       lighter stop mass
C       stoph       heavier stop mass
C       xstop       stop mixing angle (in units of PI)
C       sb1,sb2     sbottom masses 
C       xsbot       sbottom mixing angle (in units of PI)
C       sqk         common left-squark mass
C       gl         mgluino
C       aat        At soft parameter
C       abo        Ab soft parameter
C       amu         MU parameter
C       ch_m         lightest chargino mass
C       u,v         chargino U,V matrices
C       Output:
C       c70         C7 (LO)
C       c71         C7 (NLO)
C       c80         C8 (LO)
C       c81         C8 (NLO)
C       ee          E(x)    
C       ierr        OK (0), manca soluzione per M (1,2)
C       R           R_delta = Delta_susy/Delta_SM per bbbar boxes
        implicit double precision (a-h,o-z)
        dimension u(2,2),v(2,2),chm(2),tmix(2,2),acc(2,2),accd(2),
     +  str(2),xsqc(2),xstc(2,2),yy(2),uno(2),due(2),dues(2),aacc(2,2),
     +  sbr(2),bmix(2,2)
        double precision nem(4), nn(4,4)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad

C       PARAMETRI DI INPUT
        wm=80.419d0
        zm=91.1867d0

C       ALTRI PARAMETRI
        pi=4.*datan(1.d0)
        pi2=pi**2
        sq2=dsqrt(2.d0)
        r23=2.d0/3.d0
        sinb=tanb/dsqrt(1+tanb**2)
        cosb=sinb/tanb
        asc=as
        tc=t
        stoplc=stopl
        stophc=stoph
        xstopc=xstop
        xsbotc=xsbot
        amuc= amu
        sqkc=sqk
        glc=gl
        ioc=io


C       SCELTA MODELLO HIGGS
        au=1./tanb
        ad=-tanb         ! Modello II
C       ad=1./tanb       ! Modello I

C       MASSE RUNNING A SCW
c        topr= runt(scw)
        topr = runm(scw,6)
c       print*, topr,tc,runt(tc)
        xt= (topr/wm)**2
        if(imod.ge.1) then
          yt=(topr/h)**2
        end if
        if(imod.ge.2) then
        str(1)=stoph         ! no running mass for stops
        str(2)=stopl
c
        tmix(1,1)=dcos(xstop*pi)
        tmix(1,2)=dsin(xstop*pi)
        tmix(2,1)=-tmix(1,2)
        tmix(2,2)=tmix(1,1)

cgg add sbottom
        sbr(1)=sb2         ! no running mass for sbottoms
        sbr(2)=sb1
        bmix(1,1)=dcos(xsbot*pi)
        bmix(1,2)=dsin(xsbot*pi)
        bmix(2,1)=-bmix(1,2)
        bmix(2,2)=bmix(1,1)

c - now the chargino routine is called outside
c       call chargino(tanb,ch_m,amu,amg,chm,u,v,ierr)
c       if(ierr.ne.0) then
c          stop 'error chargino'
c          WRITE(6,*)'error chargino'
c          return
c          goto 3
c       endif

c this is for the susy limit only
c       ierr=0
c       u(1,1)=1./sq2
c       u(1,2)=1./sq2
c       u(2,1)=-1./sq2
c       u(2,2)=1./sq2
c       v(1,1)=1./sq2
c       v(1,2)=-1./sq2
c       v(2,1)=1./sq2
c       v(2,2)=1./sq2
c       chm(1)=wm
c       chm(2)=wm

        ierr= 0

        do i=1,2
        xsqc(i)=(sqk/chm(i))**2
        do k=1,2
        xstc(k,i)=(str(k)/chm(i))**2
        enddo
        enddo
        endif   

C       COEFF WILSON A SCW (LO)
        c70=ff1(xt) 
        c80=fg1(xt) 

c ==> NNB electroweak corrections only for the calculation of c70b
c be careful that it is not used both here and in the routine bsg
c          c70= c70 *(0.974d0) !-0.005d0)
c          c80= c80 *(0.993d0) ! -0.005d0)

        if(imod.ge.1)then
         c70=c70+au**2/3.*ff1(yt)-au*ad*ff2(yt)
         c80=c80+au**2/3.*fg1(yt)-au*ad*fg2(yt)
        endif

        if(imod.ge.2)then

         sint= tmix(1,2)
         cost= tmix(1,1)

         sbmx= bmix(1,2)
         cbmx= bmix(1,1)

         do j=1,2
         accd(j)=u(j,2)*wm/(sq2*cosb*chm(j))
         yy(j)= v(j,2)*tmix(2,2)*topr/(sq2*wm*sinb)
         do k=1,2
c       this is t_jk, t_1j=acc(j,1); t_2j=-acc(j,2)
         acc(j,k)=v(j,1)*tmix(k,1)-v(j,2)*tmix(k,2)*topr/(sq2*wm*sinb)
         enddo
         enddo

         c70s=0.d0
         c80s=0.d0
         do j=1,2
c NB for the limit of heavy squarks but stop2 comment next4 lines + do k=2,2
        c70s=c70s+r23*v(j,1)**2*(wm/sqk)**2*ff1(xsqc(j))
        c80s=c80s+r23*v(j,1)**2*(wm/sqk)**2*fg1(xsqc(j))
        c70s=c70s+accd(j)*v(j,1)*ff3(xsqc(j))
        c80s=c80s+accd(j)*v(j,1)*fg3(xsqc(j))

         do k=1,2  ! k is the stop index
        c70s=c70s-r23*acc(j,k)**2*(wm/str(k))**2*ff1(xstc(k,j))
        c80s=c80s-r23*acc(j,k)**2*(wm/str(k))**2*fg1(xstc(k,j))
        c70s=c70s-accd(j)*acc(j,k)*tmix(k,1)*ff3(xstc(k,j))
        c80s=c80s-accd(j)*acc(j,k)*tmix(k,1)*fg3(xstc(k,j))
         enddo

c      large tan beta terms
c        uno(j)=u(j,2)*v(j,1)*wm/(sq2*cosb*chm(j))
c         due(j)=u(j,2)*v(j,2)*topr/(2.d0*cosb*sinb*chm(j))*sint*cost

c        c70s= c70s+uno(j)*(ff3(xsqc(j))-
c     &    cost**2*ff3(xstc(1,j))-sint**2*ff3(xstc(2,j)))
c        c70s= c70s+due(j)*(ff3(xstc(1,j))-ff3(xstc(2,j)))
c        c80s=c80s+ uno(j)*(fg3(xsqc(j))-
c     &    cost**2*fg3(xstc(1,j))-sint**2*fg3(xstc(2,j)))
c        c80s=c80s+ due(j)*(fg3(xstc(1,j))-fg3(xstc(2,j)))

         enddo

         c70=c70+c70s
         c80=c80+c80s
        endif


C       COEFF WILSON A SCW (NLO)
        if(io.eq.1)then
Cgg - Standard Model part
         c71=ffsm(xt)+ffsmlog(xt)*2*log(scw/wm)
         c81=fgsm(xt)+fgsmlog(xt)*2*log(scw/wm)
         ee=esm(xt)
cgg      print*,"sm ",c71,c81
         if(imod.ge.1)then
Cgg - Charged Higgs two-doublets part
          c71=c71+ffh(yt)+ffhlog(yt)*2*log(scw/h)
          c81=c81+fgh(yt)+fghlog(yt)*2*log(scw/h)
          ee=ee+au**2*eh(yt)
         endif

         if(imod.ge.2 .and. nlosusy.eq.1 )then
Cgg - SUSY part: depends on which approximation is valid ...
          c71s=0.d0
          c81s=0.d0
          ees=0.d0
          do j=1,2
            ees = ees + acc(j,2)**2*(wm/str(2))**2*echi(xstc(2,j))
          enddo

          if( ihv.EQ.0 ) then   ! This is the old light stop_R case

           do j=1,2
C         ees = ees + acc(j,2)**2*(wm/str(2))**2*echi(xstc(2,j))

c  notice that in the following sqk is used for all the squarks of 
c  light flavors
c  the functions g7,8chi1,2 are not exactly the same as in the paper
c   but they are combinations 
       g7chi1= -8.d0/9.d0*ff1(xstc(2,j))*(-1.d0+6.d0*log(glc/chm(j))+
     -        2.d0*delta1(sqk) ) -
     -        4/9.d0 *echi(xstc(2,j)) + xstc(2,j)*ffs5(xstc(2,j)) +
     -        32.*(fg1(xstc(2,j))-3.*ff1(xstc(2,j)))*2*
     -        log(scw/chm(j))/27.d0     
          g8chi1= -8.d0/9.d0*fg1(xstc(2,j))*(-1.d0 +
     -        2.*delta1(sqk)+6.d0*log(glc/chm(j)))-
     -        1/6.d0 *echi(xstc(2,j)) + xstc(2,j)*fgs5(xstc(2,j)) -
     -        28.*fg1(xstc(2,j))*2*log(scw/chm(j))/9.d0 
          g7chi2= -4.d0/3.d0*ff3(xstc(2,j))*(2.*delta1(sqk)+
     -        delta2(sqk,sqk)-2.d0) - ffs1(xstc(2,j)) +
     -        16.*(fg3(xstc(2,j)) -3.*ff3(xstc(2,j)))*2*
     -        log(scw/chm(j))/9.d0      
          g8chi2= -4.d0/3.d0*fg3(xstc(2,j))*(2.*delta1(sqk)+
     -        delta2(sqk,sqk)-2.d0)  - fgs1(xstc(2,j)) -
     -        14.*fg3(xstc(2,j))*2*log(scw/chm(j))/3.d0 
          g7chi3= 8.d0/9.d0*ff1(xstc(2,j))*(12.*log(scw/glc)+
     -          2.d0*deltat2(sqk)-2.d0)
          g8chi3= 8.d0/9.d0*fg1(xstc(2,j))*(12.*log(scw/glc)+
     -          2.d0*deltat2(sqk)-2.d0)
          g7chi4= 4.d0/3.d0*ff3(xstc(2,j))*(6.*log(scw/glc)+
     -          deltat2(sqk)-1.d0)
          g8chi4= 4.d0/3.d0*fg3(xstc(2,j))*(6.*log(scw/glc)+
     -          deltat2(sqk)-1.d0)

          c71s= c71s + acc(j,2)**2*(wm/str(2))**2*g7chi1  +      
     -          accd(j)*tmix(2,1)*acc(j,2)*g7chi2 +
     -          acc(j,2)*yy(j)*(wm/str(2))**2*g7chi3 +
     -          accd(j)*yy(j)*tmix(2,1)*g7chi4
cgg       print*,"c71ss j",j, delta2(sqk,sqk)

          c81s= c81s +acc(j,2)**2*(wm/str(2))**2*g8chi1 +
     -          accd(j)*tmix(2,1)*acc(j,2)*g8chi2+
     -          acc(j,2)*yy(j)*(wm/str(2))**2*g8chi3 +
     -          accd(j)*yy(j)*tmix(2,1)*g8chi4  

          enddo
c       additional contributions from shifts of SM and 2HDM
c       NB in ..sw we use ww(sqk,0)=1/2 for sqk=mgl
          cf= 4.d0/3.d0
          c71sw= cf*(tmix(1,1)**2*ww(sqk,stophc)+
     -               tmix(1,2)**2/2.d0)*gg7(xt)
c         print*,c71sw
          c81sw= cf*(tmix(1,1)**2*ww(sqk,stophc)+
     -               tmix(1,2)**2/2.d0)*gg8(xt)
          c71sp= cf*(2*utd(sqk)*ff1(xt)/3.d0  
     -            -(utd(sqk)+udd(sqk,sqk))*ff2(xt))
c         print*,c71sp
          c81sp= cf*(2*utd(sqk)*fg1(xt)/3.d0
     -            -(utd(sqk)+udd(sqk,sqk))*fg2(xt))
          c71sh= cf*(2*au**2/3.d0*ftd(sqk)*ff1(yt) 
     -            +(ftd(sqk)+fdd(sqk,sqk))*ff2(yt))
c         print*,c71sh
          c81sh= cf*(2*au**2/3.d0*ftd(sqk)*fg1(yt) 
     -            +(ftd(sqk)+fdd(sqk,sqk))*fg2(yt))
c if you want to have only LO susy coeff comment the following 3 lines
          c71=c71 +c71s+c71sw+c71sp+c71sh
          c81=c81 +c81s+c81sw+c81sp+c81sh

          else if( ihv.EQ.1 ) then ! This is the new heavy colored particles case 

cc LARGE LOGS ONLY (ex large tanb scenario)
c >>>> NNB Msusy is set HERE!!! 
C         amsusy= 1200.d0 ! msusy, to be set equal to average?
          amsusy= sqk ! set equal to average squark mass

c  NB in dmb H2=-1/2 this is true only for degenerate sbottoms!!
          dmb=asf(amsusy)/pi/3.d0*amuc*tanb/glc

         
c  notation as in the new paper
          xx1=(sb1/glc)**2
          xx2=(sb2/glc)**2
          u1=(stophc/glc)**2
          u2=(stoplc/glc)**2
          y11=(stophc/chm(1))**2
          y12=(stophc/chm(2))**2
          y21=(stoplc/chm(1))**2
          y22=(stoplc/chm(2))**2
          yukt= 0.65d0*runt(amsusy)/wm/sq2 !g(mz)=0.65

          epsb= -asf(amsusy)*2.d0/3.d0/pi*amuc/glc*h2(xx1,xx2)

          epsbew=-(yukt/4/pi)**2*aat/chm(1)*u(1,2)*v(1,2)*h2(y11,y21) - 
     &    (yukt/4/pi)**2*aat/chm(2)*u(2,2)*v(2,2)*h2(y12,y22)  
c         print*,"epsbew",epsbew
          epsbew=0.

          isb=1
          if( isb.eq.0 ) then
            epsbp= -asf(amsusy)*2.d0/3.d0/pi*amuc/glc*
     &        (cost**2*h2(u1,xx2)+sint**2*h2(u2,xx2))
          else if( isb.eq.1 ) then
            epsbp= -asf(amsusy)*2.d0/3.d0/pi*amuc/glc*
     &        (cost**2*(h2(u1,xx2)*cbmx**2+h2(u1,xx1)*sbmx**2)+
     &         sint**2*(h2(u2,xx2)*cbmx**2+h2(u2,xx1)*sbmx**2))
          endif         

          if( isb.eq.0 ) then
            epst= -asf(amsusy)*2.d0/3.d0/pi*amuc/glc*
     &        (cost**2*h2(u2,xx1)+sint**2*h2(u1,xx1))
          else if( isb.eq.1 ) then
            epst= -asf(amsusy)*2.d0/3.d0/pi*amuc/glc*
     &        (cost**2*(h2(u2,xx1)*cbmx**2+h2(u2,xx2)*sbmx**2)+
     &         sint**2*(h2(u1,xx1)*cbmx**2+h2(u1,xx2)*sbmx**2))
          endif         
c         print*,"x12",epsb,epsbp,epst,sb1,sb2,aat

          akk= 1/(1+(epsb+epsbew)*tanb)
          eta= asf(amsusy)/asf(scw)

c first redefine the couplings; 

         do j=1,2
         do k=1,2
         aacc(j,k)=v(j,1)*tmix(k,1)-
     &          v(j,2)*tmix(k,2)*runt(amsusy)/(sq2*wm*sinb)
         enddo
         enddo

          c7imp=0.d0
          c8imp=0.d0

          do j=1,2
c == first we write the improved LO coeff
        c7imp=c7imp+r23*v(j,1)**2*(wm/sqk)**2*ff1(xsqc(j))
        c8imp=c8imp+r23*v(j,1)**2*(wm/sqk)**2*fg1(xsqc(j))
        c7imp=c7imp+akk*accd(j)*v(j,1)*ff3(xsqc(j))
        c8imp=c8imp+akk*accd(j)*v(j,1)*fg3(xsqc(j))
         do k=1,2  ! k is the stop index
        c7imp=c7imp-r23*aacc(j,k)**2*(wm/str(k))**2*ff1(xstc(k,j))
        c8imp=c8imp-r23*aacc(j,k)**2*(wm/str(k))**2*fg1(xstc(k,j))
        c7imp=c7imp-akk*accd(j)*aacc(j,k)*tmix(k,1)*ff3(xstc(k,j))
        c8imp=c8imp-akk*accd(j)*aacc(j,k)*tmix(k,1)*fg3(xstc(k,j))
         enddo
         enddo

c  only large tan beta
c        c7imp= c7imp+uno(j)*(ff3(xsqc(j))-
c     &    cost**2*ff3(xstc(1,j))-sint**2*ff3(xstc(2,j)))
c        c7imp= c7imp+dues(j)*(ff3(xstc(1,j))-ff3(xstc(2,j)))    
c        c8imp= c8imp+uno(j)*(fg3(xsqc(j))-
c     &    cost**2*fg3(xstc(1,j))-sint**2*fg3(xstc(2,j)))
c        c8imp=c8imp+ dues(j)*(fg3(xstc(1,j))-fg3(xstc(2,j)))
c        enddo 

cgg      print*,"c7imp",c7imp,c8imp,runt(amsusy)

c ==   First, the LL from the running of the ops
c         c7ll= -32.d0/3.d0*(c70s-c80s/3.d0)*log(msy/scw)
c         c8ll= -28.d0/3.d0*c80s*log(msy/scw)
c == the same but resummed

         c7ll=eta**(16.d0/21.d0)*c7imp+8.d0/3.d0*(eta**(14.d0/21.d0)-
     -    eta**(16.d0/21.d0))*c8imp !-c70s
c        c7ll=4*pi/asf(scw)*c7ll
         c8ll=eta**(14.d0/21.d0)*c8imp !-c80s
c        c8ll=4*pi/asf(scw)*c8ll

c        print*,"c7ll",c7ll,c8ll,c70s,c80s,akk

         c7sg=tanb*(epsb-epsbp)*akk*ff2(xt)/(asf(amsusy)/4.d0/pi)
         c8sg=tanb*(epsb-epsbp)*akk*fg2(xt)/(asf(amsusy)/4.d0/pi)
        
c         c7sg= -8.d0/3.d0*amuc*tanb/glc*ff2(xt) *
c     -   (-0.5d0 -cost**2*h2((stophc/glc)**2,0.99d0)-
c     -    sint**2*h2((stoplc/glc)**2,0.99d0))    
c         c8sg= -8.d0/3.d0*amuc*tanb/glc*fg2(xt)*
c     -   (-0.5d0 -cost**2*h2((stophc/glc)**2,0.99d0)-
c     -    sint**2*h2((stoplc/glc)**2,0.99d0))
cgg       print*,"c7sg",c7sg,c8sg,akk

          c7sh= 8.d0/3.d0*amuc*tanb/glc*ff2(yt)*
     -    (-0.5d0 +cost**2*h2((stophc/glc)**2,0.99d0)+
     -    sint**2*h2((stoplc/glc)**2,0.99d0))     
          c8sh= 8.d0/3.d0*amuc*tanb/glc*fg2(yt)*
     -    (-0.5d0 +cost**2*h2((stophc/glc)**2,0.99d0)+
     -    sint**2*h2((stoplc/glc)**2,0.99d0))
                  
cgg       print*,c7sg,dmb,c7ll

c if you want to have only LO susy coeff 
C         if(nlosusy.eq.1)then
             c71s=c7sg          !+ c7ll+c7sh
             c81s=c8sg          !+c8ll+c8sh
             c70=c70 +c7ll-c70s
             c80=c80 +c8ll-c80s
             c71=c71 +c71s !+c71sw+c71sp+c71sh these are only for light stop2
             c81=c81 +c81s !+c81sw+c81sp+c81sh
C         endif

        endif ! switch on ihv

             ee= ee + ees
         endif
        endif


C       the ratio R_delta for the B-Bbar mixing and epsilon

        deltabox=aab(xt)                !SM contribution
c       print*,"d",deltabox
        if(imod.ge.1)then       !2HDM contribution
         deltabox=deltabox + au**4*yt*gboxh(yt)/4.d0  + 
     -       2*au**2*xt*(fpbox(xt,xt/yt)+ gpbox(xt,xt/yt)/4.d0)
        endif

        if(imod.ge.2)then       ! here we include only the light stop 
           deltaboxs=0.d0
           do ibox=1,2
              do jbox=1,2
                 deltaboxs=deltaboxs+  acc(jbox,2)**2*acc(ibox,2)**2*
     -            (wm**2/chm(jbox)**2)/xt*       
     -           gpbox(xstc(2,jbox), (chm(ibox)**2/chm(jbox)**2))
              enddo
           enddo
           deltabox=deltabox+deltaboxs
c               print*,deltaboxs
        endif

        bbox= deltabox/aab(xt)

 3      return
        end


C***************************************************************

        double precision function runt(x)
C       CALCOLA M_top^{running} (X)
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        if(x.le.tc)then
         fn=5.
        else
         fn=6.
        endif
        b0=11.-2*fn/3.
        b1=102.-38.*fn/3.
        g0=8.
        g1=404./3.-40*fn/9.
        asx=asf(x)
        ast=asf(tc)
        rrr=tc*(asx/ast)**(g0/(2*b0))
        if(ioc.eq.0)then
         runt=rrr
        else
         corr=1+ast*g0/(4*pi*2*b0)*(-b1/b0+g1/g0)*(asx/ast-1)
c this is the relation between mpole/mt(mpole)
c        pol= 1+4*ast/(3*pi) !+10.9d0*(ast/pi)**2
c this is the relation between mpole/mtrun(mtrun)
         pol=1+4*ast/(3*pi)+8.243d0*(ast/pi)**2
c this is the 1loop order result 
c        pol=1+4*ast/(3*pi)
         runt=rrr*corr/pol
        endif
        return
        end

C***************************************************************
c we do not use stop running mass
c       double precision function runst(x)
C       CALCOLA M_stop^{running} (X)
c       implicit double precision (a-h,o-z)
c       common/alp/asc,tc,pi,zm,ioc
c       common/st/stoplc,stophc
c       runst=0.
c       return
c       end
C***************************************************************

        double precision function asf(x)
C       CALCOLA ALPHA_S (X)
c NB it uses 5 flavors for x<mt and 6 if x>=mt
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        fn=5.
        b0=11.-2*fn/3.
        b1=102.-38.*fn/3.
        vvv=1-b0*asc/(2*pi)*log(zm/x)
        if(ioc.eq.0)then
         asf=asc/vvv   
        else
         asf=asc/vvv*(1-b1/b0*asc/(4*pi*vvv)*log(vvv))
        endif
        if(x.gt.tc)then
         vvv=1-b0*asc/(2*pi)*log(zm/tc)
        if(ioc.eq.0)then
         ast=asc/vvv
        else
         ast=asc/vvv*(1-b1/b0*asc/(4*pi*vvv)*log(vvv))
        endif
        b0t=b0-2./3.
        b1t=b1-38./3.
        vvv=1-b0t*ast/(2*pi)*log(tc/x)
        if(ioc.eq.0)then
         asf=ast/vvv
        else
         asf=ast/vvv*(1-b1t/b0t*ast/(4*pi*vvv)*log(vvv))
        endif
        endif
        return
        end

C***************************************************************

        double precision function sp(x)
C       FUNZIONE DILOG
        implicit double precision (a-h,o-z)
        external f
        pi=4.*datan(1.d0)
        if(x.ge.-1.d0.and.x.le.0.5d0)then
         sp=f(x)
        else if(x.gt.0.5d0.and.x.lt.1.d0)then
         sp=-f(1.d0-x)+pi**2/6.-dlog(x)*dlog(1.d0-x)
        else if(x.lt.-1.d0)then
         sp=-f(1.d0/x)-pi**2/6-.5d0*(dlog(-x))**2
        else
         write(6,*)'error in dilog',x
        endif
        return
        end

        double precision function f(x)
        implicit double precision (a-h,o-z)
        dimension b(12)
        z=-dlog(1.d0-x)
        b(1)=-.5d0
        b(2)=1.d0/6.
        b(3)=0.
        b(4)=-1.d0/30.
        b(5)=0.d0
        b(6)=1.d0/42.
        b(7)=0.d0
        b(8)=-1.d0/30.
        b(9)=0.d0
        b(10)=5.d0/66.
        b(11)=0.d0
        b(12)=-691.d0/2730.
        cc=z
        sum=z   
        do i=1,12
        cc=cc*z/(i+1.d0)
        sum=sum+b(i)*cc
        enddo
        f=sum
        end

C*************************************************************************

        double precision function ff1(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        d=1./(x-1)**3
        dg=log(x)/(x-1)**4
        ff1=x*(7-5*x-8*x**2)*d/24.+x**2*(3*x-2)*dg/4.
        else
        ff1=-5./48.
        endif
        return
        end

C*************************************************************************

        double precision function fg1(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
        d=1./(x-1)**3
        dg=log(x)/(x-1)**4
        fg1=x*(2+5*x-x**2)*d/8.-3*x**2*dg/4.
        else
        fg1=-1./16.
        endif
        return
        end

C*************************************************************************

        double precision function ff2(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
        d=1./(x-1)**2
        dg=log(x)/(x-1)**3
        ff2=x*(3-5*x)*d/12.+x*(3*x-2)*dg/6.
        else
        ff2=-7./36.
        endif
        return
        end

C*************************************************************************

        double precision function fg2(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
        d=1./(x-1)**2
        dg=log(x)/(x-1)**3
        fg2=x*(3-x)*d/4.-x*dg/2.
        else
        fg2=-1./6.
        endif
        return
        end

C*************************************************************************

        double precision function ff3(x)
        implicit double precision (a-h,o-z)
        ff3=2.d0/3.d0*(1-1.d0/x)*ff1(x)+ff2(x)+23.d0/36.d0
        return
        end

C*************************************************************************

        double precision function fg3(x)
        implicit double precision (a-h,o-z)
        fg3=2.d0/3.d0*(1-1.d0/x)*fg1(x)+fg2(x)+1.d0/3.d0
        return
        end

C*************************************************************************

        double precision function gg7(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
        d=1./(x-1)**3
        dg=log(x)/(x-1)**4
        gg7=-(23-67*x+50*x**2)*d/36.d0 +x*(2-7*x+6*x**2)*dg/6.d0
        else
        gg7=3./8.d0
        endif
        return
        end

C*************************************************************************

        double precision function gg8(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        d=1./(x-1)**3
        dg=log(x)/(x-1)**4
        gg8=-(4-5*x-5*x**2)*d/12.d0 +x*(1-2*x)*dg/2.d0
        else
        gg8=1./8.d0
        endif
        return
        end

C*************************************************************************

        double precision function esm(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        d=1./(x-1)**3
        dg=log(x)/(x-1)**4
        esm=x*(-18+11*x+x**2)*d/12.+x**2*(15-16*x+4*x**2)*dg/6.
        esm=esm-2./3.*log(x)
        else
        esm=43./72.
        endif
        return
        end

C*************************************************************************

        double precision function eh(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        d=1./(x-1.)**3
        dg=log(x)/(x-1.d0)**4
        eh=x*(16-29*x+7*x**2)*d/36.d0+x*(3*x-2)*dg/6.d0
        else
        eh=1./8.d0
        endif
        return
        end
C*************************************************************************
c NB there is a misprint in Gabrielli-Giudice a factor2
c moreover this corresponds to E(1/x) in GG notation
        double precision function echi(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        d=1./(x-1.d0)**3
        dg=log(x)/(x-1.d0)**4
        echi= x*(11.d0 - 7.d0*x + 2.d0*x**2)*d/18.d0 - x*dg/3.d0
        else
        echi=1./12.d0
        endif
        return
        end

C*************************************************************************

        double precision function ffsm(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        x2=x*x
        x3=x2*x
        x4=x3*x
        x5=x4*x
        d4=(x-1)**4
        d5=d4*(x-1)
        spx=sp(1-1./x)
        xl=log(x)
        xl2=xl**2
        ffsm=(-16*x4-122*x3+80*x2-8*x)*spx/(9*d4)+(6*x4+46*x3-28*x2)*xl2
     +   /(3*d5)+(-102*x5-588*x4-2262*x3+3244*x2-1364*x+208)*xl/(81*d5)
        ffsm=ffsm+(1646*x4+12205*x3-10740*x2+2509*x-436)/(486*d4)
        else
        ffsm=-3451./9720.
        endif
        return
        end

C*************************************************************************

        double precision function fgsm(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        x2=x*x
        x3=x2*x
        x4=x3*x
        x5=x4*x
        d4=(x-1)**4
        d5=d4*(x-1)
        spx=sp(1-1./x)
        xl=log(x)
        xl2=xl**2
        fgsm=(-4*x4+40*x3+41*x2+x)*spx/(6*d4)+(-17*x3-31*x2)*xl2/(2*
     +   d5)+(-210*x5+1086*x4+4893*x3+2857*x2-1994*x+280)*xl/(216*d5)
        fgsm=fgsm+(737*x4-14102*x3-28209*x2+610*x-508)/(1296*d4)
        else
        fgsm=-9821./12960.
        endif
        return
        end

C*************************************************************************

        double precision function ffsmlog(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        diff=(-x*(18*x**2+30*x-24)*log(x)+47*x**3-63*x**2+9*x+7)/24.
     +   /(x-1)**5
        ffsmlog=8*x*diff+16./3.*ff1(x)-16./9.*fg1(x)+208./81.
        else
        ffsmlog=1369./810.
        endif
        return
        end

C*************************************************************************

        double precision function fgsmlog(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
        diff=(6*x*(1+x)*log(x)-x**3-9*x**2+9*x+1)/4.d0
     +   /(x-1)**5
        fgsmlog=8*x*diff+14./3.d0*fg1(x)+35./27.d0
        else
        fgsmlog=869./1080.d0
        endif
        return
        end

C*************************************************************************

        double precision function ffh(wth) !this includes -4/9 E_H
        implicit double precision (a-h,o-z)
        common/bs_hig/au,ad
        th= wth
        if( abs(wth-1d0).lt.1d-3 ) th= 1d0+sign(1d0,(wth-1d0))*1d-3
        ffh = Ad*Au*(4.d0*(7.d0*th - 13.d0*th**2 + 2.d0*th**3)/
     -      (3.d0*(-1.d0 + th)**3) + 
     -     16.d0*(-3.d0*th + 7.d0*th**2 - 2.d0*th**3)*
     -     sp(1.d0 - 1.d0/th)/(9.d0*(-1.d0 + th)**3) + 
     -     8.d0*(-3.d0*th - th**2+12.d0*th**3 - 2.d0*th**4)*Log(th)/
     -     (9.d0*(-1.d0 + th)**4) + 
     -   4.d0*(8.d0*th - 14.d0*th**2 - 3.d0*th**3)*Log(th)**2/
     -    (9.d0*(-1.d0 + th)**4)) + 
     -  Au**2*((893.d0*th - 5706.d0*th**2+7785.d0*th**3- 1244*th**4)/
     -      (486.d0*(-1.d0 + th)**4) + 
     -     2.d0*(18.d0*th**2 - 37.d0*th**3 + 8.d0*th**4)*
     -     sp(1.d0 - 1.d0/th)/(9.d0*(-1.d0 + th)**4) + 
     -     2.d0*(-56.d0*th+266.d0*th**2-183.d0*th**3- 192.d0*th**4 + 
     -     21.d0*th**5)*Log(th)/(81.d0*(-1.d0 + th)**5) + 
     -     2.d0*(-14.d0*th**2 + 23.d0*th**3 + 3.d0*th**4)*Log(th)**2/
     -     (9.d0*(-1.d0 + th)**5)) 
        return
        end

C*************************************************************************

        double precision function fgh(wth)
        implicit double precision (a-h,o-z)
        common/bs_hig/au,ad
        th= wth
        if( abs(wth-1d0).lt.1d-3 ) th= 1d0+sign(1d0,(wth-1d0))*1d-3
        fgh = Ad*Au*((143.d0*th - 44.d0*th**2 + 29.d0*th**3)/
     -      (8.d0*(-1.d0 + th)**3) + 
     -   (-36.d0*th + 25.d0*th**2 - 17.d0*th**3)*sp(1.d0 - 1.d0/th)/
     -     (6.d0*(-1.d0 + th)**3) + 
     -   (-3.d0*th - 187.d0*th**2 + 12.d0*th**3 - 14.d0*th**4)*Log(th)/
     -    (12.d0*(-1.d0 + th)**4) + (19.d0*th + 17.d0*th**2)*Log(th)**2/
     -    (3.d0*(-1.d0 + th)**4) )
     -      + Au**2*((1226.d0*th - 18423.d0*th**2 + 7866.d0*th**3 - 
     -      4493.d0*th**4)/(1296.d0*(-1.d0 + th)**4) + 
     -     (30.d0*th**2 - 17.d0*th**3 + 13.d0*th**4)*
     -     sp(1.d0 - 1.d0/th)/(6.d0*(-1.d0 + th)**4) + 
     -     (-238.d0*th + 847.d0*th**2 + 1335.d0*th**3 + 318.d0*th**4 + 
     -     42.d0*th**5)*Log(th)/(216.d0*(-1.d0 + th)**5) - 
     -     (31.d0*th**2 + 17.d0*th**3)*Log(th)**2.d0/
     -    (6.d0*(-1.d0 + th)**5)) 
        return
        end

C*************************************************************************

        double precision function ffhlog(wth)
        implicit double precision (a-h,o-z)
        common/bs_hig/au,ad
        th= wth
        if( abs(wth-1d0).lt.1d-3 ) th= 1d0+sign(1d0,(wth-1d0))*1d-3
        ffhlog=-(     (Ad*Au*
     -      (2.d0*(-21.d0*th + 47.d0*th**2 - 8.d0*th**3)/
     -     (9.d0*(-1.d0 + th)**3) + 
     -    4.d0*(8.d0*th - 14.d0*th**2 - 3.d0*th**3)*Log(th)/
     -    (9.d0*(-1.d0 + th)**4)) + 
     - Au**2*((31.d0*th + 18.d0*th**2 - 135.d0*th**3 + 14.d0*th**4)/
     -    (27.d0*(-1.d0 + th)**4) + 
     -   2.d0*(-14.d0*th**2 + 23.d0*th**3 + 3.d0*th**4)*Log(th)/
     -   (9.d0*(-1.d0 + th)**5)))   )
        return
        end

C*************************************************************************

        double precision function fghlog(wth)
        implicit double precision (a-h,o-z)
        common/bs_hig/au,ad
        th= wth
        if( abs(wth-1d0).lt.1d-3 ) th= 1d0+sign(1d0,(wth-1d0))*1d-3
        fghlog = -( (Ad*Au*
     -  ((-81.d0*th + 16.d0*th**2 - 7.d0*th**3)/(6.d0*(-1.d0 + th)**3) + 
     -        (19.d0*th + 17.d0*th**2)*Log(th)/(3.d0*(-1.d0 + th)**4)) + 
     -     Au**2*((38.d0*th + 261.d0*th**2 - 18.d0*th**3 + 7.d0*th**4)/
     -         (36.d0*(-1.d0 + th)**4) - 
     -  (31.d0*th**2 + 17.d0*th**3)*Log(th)/(6.d0*(-1.d0 + th)**5))) )
        return
        end
C*************************************************************************

        double precision function ffs1(wth)
        implicit double precision (a-h,o-z)
        x= wth
        if( abs(wth-1d0).lt.1d-3 ) x= 1d0+sign(1d0,(wth-1d0))*1d-3
        ffs1 =   (-4.d0*(-5.d0 + 3.d0*x))/(9.d0*(-1.d0 + x)**2) + 
     -   (16.d0*(3.d0-7.d0*x)*x*sp(1.d0 - 1.d0/x))/
     -   (9.d0*(-1.d0 + x)**3) + 
     -   (4.d0*(4.d0-30.d0*x + 40.d0*x**2)*Log(x))/
     -   ( 9.d0*(-1.d0 + x)**3) + 
     -   (16.d0 *(1.d0 - 3.d0*x)*x* Log(x)**2) / 
     -          (9.d0 *(-1.d0 + x )**3)
        return
        end

C*************************************************************************

        double precision function ffs5(wth)
        implicit double precision (a-h,o-z)
        x= wth
        if( abs(wth-1d0).lt.1d-3 ) x= 1d0+sign(1d0,(wth-1d0))*1d-3
        ffs5 =(-85.d0 + 347.d0*x - 526.d0*x**2)/(243.d0*(-1.d0 +x)**3)+
     -   (4.d0*x*(-8.d0 + 13.d0*x + 6.d0*x**2)*sp(1.d0 - 
     -     1.d0/x))/(9.d0*(-1.d0 + x)**4) + 
     -   (4.d0*(-20.d0 + 126.d0*x - 144.d0*x**2 - 
     -    39.d0*x**3)*Log(x))/(81.d0*(-1.d0 + x)**4) + 
     -   (2.d0*x*(-10.d0 + 21.*x)*Log(x)**2)/ (9.d0*(-1.d0 + x)**4)
        return
        end

C*************************************************************************

        double precision function fgs1(wth)
        implicit double precision (a-h,o-z)
        x= wth
        if( abs(wth-1d0).lt.1d-3 ) x= 1d0+sign(1d0,(wth-1d0))*1d-3
        fgs1 =  (61.d0 - 39.d0*x)/(12.d0 *(-1.d0 + x)**2) + 
     -   (4.d0*x*(3.d0 + 4.d0*x)*sp(1.d0 - 1.d0/x))/
     -    (3.*(-1.d0 + x)**3) + 
     -   ((7.d0 - 60.d0*x - 14.d0*x**2)*Log(x)) / 
     -   (6.d0 *(-1.d0 + x)**3) + 
     -   (14.d0*x*Log(x)**2)/(3.d0*(-1.d0 + x)**3)
        return
        end

C*************************************************************************

        double precision function fgs5(wth)
        implicit double precision (a-h,o-z)
        x= wth
        if( abs(wth-1d0).lt.1d-3 ) x= 1d0+sign(1d0,(wth-1d0))*1d-3
        fgs5 = (-1210.d0 + 437.d0*x+1427.d0*x**2)/
     -   (648.d0*(-1.d0+x)**3) - 
     - (x*(49.d0+46.d0*x + 9.d0*x**2)*sp(1.d0 - 1.d0/x))/
     -    (12.d0 *(-1.d0 + x)**4) + 
     -   ((-85.d0 + 603.d0*x + 387.d0*x**2 - 
     -     78.d0*x**3)*Log(x))/(108.d0*(-1.d0 + x)**4) - 
     -   (13.d0*x*Log(x)**2)/(3.d0*(-1.d0 + x)**4)
        return
        end

C*************************************************************************
        double precision function h1(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
           h1=1.d0/(1-x) + (2*x - x**2)*Log(x)/(1-x)**2
        else
           h1=-1/2.d0
        endif
        return
        end
C*************************************************************************
        double precision function h3(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
           h3=2.d0*x*Log(x)/(1-x)
        else
           h3=-2.d0
        endif
        return
        end
C*************************************************************************
        double precision function h4(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-2)then
           h4=2.d0*Log(x)/(1-x)
        else
           h4=-2.d0
        endif
        return
        end
C*************************************************************************
        double precision function h2(x,y)
        implicit double precision (a-h,o-z)
        if(abs(x-y).gt.1.d-2)then
           if(abs(x-1.).lt.1.d-2)then
              h2=1.d0/(-1 + y) - y*Log(y)/(-1 + y)**2
           else 
              if(abs(y-1.).lt.1.d-2)then
                 h2=1.d0/(-1 + x) - x*Log(x)/(-1 + x)**2
              else      
                 h2= x*Log(x)/((1 - x)*(x - y)) + 
     -                y*Log(y)/((1 - y)*(-x + y))
              endif 
           endif  
        else
           if(abs(x-1.).lt.1.d-2)then   
              h2=-1/2.d0                
           else
              h2=1.d0/(1 - x) + Log(x)/(-1 + x)**2
           endif
        endif      
        return
        end
C*************************************************************************

        double precision function delta1(sq) !checked
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        x= (sq/glc)**2
        sst= dsin(pi*xstopc) 
        cst= dcos(pi*xstopc) 
           delta1 = -3.d0/4.d0 - h1(x)/2.d0 !+ 2.d0*sst*cst*tc/glc
c       the last term has been neglected according to the paper
           return
        end

C*************************************************************************

        double precision function deltat2(sqd1) ! checked
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        t1= (stophc/glc)**2
        d1= (sqd1/glc)**2

        tlr =dsin(xstopc*pi)*dcos(xstopc*pi)*(stophc**2-stoplc**2)
        deltat2= 3.d0+ h3(d1)+ h1(t1)/2.d0 - (tlr/tc)*h4(t1)/glc 
        return
        end

C*************************************************************************

        double precision function delta2(sqd1,sqd2) !  checked
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad
      COMMON/SU_atri3/dal,dau,dad
        d1= (sqd1/glc)**2
        d2= (sqd2/glc)**2
        tb= 1.d0/au
        tlr =dsin(xstopc*pi)*dcos(xstopc*pi)*(stophc**2-stoplc**2)
c        aad= tlr/tc + amuc/tb  ! corrected with true Ad now (jlk):
        aad = dad
   
        ctt =dcos(xstopc*pi)/dsin(xstopc*pi)

        delta2= h3(d2)*(1-ctt*tc*glc/sqd2**2) +5./2.d0 + 
     -    (h1(d1)+h1(d2))/2.d0 - 2.d0*(aad-amuc*tb)*h2(d1,d2)/glc
c left over from a test?
c       delta2= - 2.d0*(-amuc*tb)*h2(d1,d2)/glc
        return
        end

C*************************************************************************

        double precision function ww(sq1,sq2) !checked
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        x= (sq1/glc)**2
        y= (sq2/glc)**2
        if(abs(x-y).gt.1.d-3)then
           if(abs(x-1.).lt.1.d-3)then
              ww=(1 + 5*y)/(2*(1 - y)) + y*(2 + y)*Log(y)/(-1 + y)**2
           else 
              if(abs(y-1.).lt.1.d-3)then
                 ww=(1 + 5*x)/(2*(1 - x)) + x*(2 + x)*Log(x)/(-1 + x)**2
              else      
           ww=        (x + y - 2*x*y)/((-1 + x)*(-1 + y)) + 
     -   (x**3 - 2*x*y + x**2*y)*Log(x)/((-1 + x)**2*(x - y)) + 
     -   (2*x*y - x*y**2 - y**3)*Log(y)/((x - y)*(-1 + y)**2)
              endif 
           endif  
        else
           ww=0.d0              
        endif      
        return
        end

C*************************************************************************
c    in the following function we have used h1(0)=1 and h2(1,1)=-1/2
        double precision function utd(sqd1) ! checked 
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad
        t1= (stophc/glc)**2
        d1= (sqd1/glc)**2
        sst= dsin(pi*xstopc) 
        cst= dcos(pi*xstopc) 

        tlr= sst*cst*(stophc**2-stoplc**2)
      utd=-h1(d1)/2.d0+(cst**2*h1(t1)+sst**2)/2.d0-(tlr/tc)*h4(t1)/glc +
     -    (tlr/tc )*(cst**2*h4(d1)-sst**2)/glc
c     -    (tlr/tc  +2.d0*amuc*au)*(cst**2*h4(d1)-sst**2)/glc
        return
        end

c for the checks h1[x_]:= 1/(1-x) + (2x-x^2) Log[x]/(1-x)^2
c h4[x_]:= 2 Log[x]/(1-x) 
C*************************************************************************

        double precision function udd(sqd1,sqd2) ! to be checked 
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad
      COMMON/SU_atri3/dal,dau,dad
        t1= (stophc/glc)**2
        d1= (sqd1/glc)**2
        d2= (sqd2/glc)**2
        tb= 1.d0/au
        sst= dsin(pi*xstopc) 
        cst= dcos(pi*xstopc) 
c  here we set Ad= At   ! corrected with true Ad now (jlk):
        tlr =sst*cst*(stophc**2-stoplc**2)
c        aad= tlr/tc + amuc/tb
        aad = dad
c   
        udd= h1(d1)/2.d0 - (cst**2*h1(t1)+sst**2)/2.d0
     -      -2.d0*(aad-amuc*tb)/glc*h2(d1,d2) +
     -     2.d0*(aad-amuc*tb)/glc*(cst**2*h2(t1,d2)+sst**2*h4(d2)/2.d0)
c    -     2.d0*(aad+amuc*tb)/glc*(cst**2*h2(t1,d2)+sst**2*h4(d2)/2.d0)
        return
        end

C*************************************************************************
c in the paper this is H_td
        double precision function ftd(sqd1) ! checked 
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad
        t1= (stophc/glc)**2
        d1= (sqd1/glc)**2
        sst= dsin(pi*xstopc) 
        cst= dcos(pi*xstopc) 
        tlr =dsin(xstopc*pi)*dcos(xstopc*pi)*(stophc**2-stoplc**2)

        ftd= -h1(d1)/2.d0 + (cst**2*h1(t1)+sst**2)/2.d0- 
     . (tlr/tc)*h4(t1)/glc+
     -    (tlr/tc +amuc*(au+1.d0/au))*(cst**2*h4(d1)-sst**2)/glc
c     -    (tlr/tc +amuc*(au-1.d0/au))*(cst**2*h4(d1)-sst**2)/glc
        return
        end
C*************************************************************************
c in the paper this is H_d
        double precision function fdd(sqd1,sqd2) ! to be checked 
        implicit double precision (a-h,o-z)
        common/alp/asc,tc,pi,zm,ioc
        common/st/stoplc,stophc,glc,xstopc,amuc
        common/bs_hig/au,ad
      COMMON/SU_atri3/dal,dau,dad
        t1= (stophc/glc)**2
        d1= (sqd1/glc)**2
        d2= (sqd2/glc)**2
        sst= dsin(pi*xstopc) 
        cst= dcos(pi*xstopc) 
        tb= 1.d0/au
c  here we set Ad= At  ! corrected with true Ad coupling now (jlk)
        tlr =dsin(xstopc*pi)*dcos(xstopc*pi)*(stophc**2-stoplc**2)
c        aad= tlr/tc + amuc/tb
        aad = dad
c
        fdd= h1(d1)/2.d0 - (cst**2*h1(t1)+sst**2)/2.d0
     -      -2.d0*(aad-amuc*tb)/glc*h2(d1,d2) +
     -     2.d0*(aad+amuc/tb)/glc*(cst**2*h2(t1,d2)+sst**2*h4(d2)/2.d0)
c     -     2.d0*(aad-amuc/tb)/glc*(cst**2*h2(t1,d2)+sst**2*h4(d2)/2.d0)
        return
        end

C*************************************************************************
c  here there are functions which are useful for the derivation of
c  the ckm elements from fit
        double precision function aab(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
           aab =(-4 +15*x -12*x**2 + x**3 + 6*x**2*Log(x))/(4.d0*(-1 + 
     -          x)**3)
        else
           aab= 3.d0/4.d0
        endif
        return
        end

C*************************************************************************
        double precision function gboxh(x)
        implicit double precision (a-h,o-z)
        if(abs(x-1.).gt.1.d-3)then
           gboxh =(-1 + x**2 - 2*x*Log(x))/(-1 + x)**3
        else
           gboxh= 1.d0/3.d0
        endif
        return
        end

C*************************************************************************
        double precision function fpbox(x,y)
        implicit double precision (a-h,o-z)
        if(abs(x-y).gt.1.d-3)then
           fpbox = (-1 + x - Log(x))/((-1 + x)**2*(-x + y)) + 
     -   (x*Log(x)/(-1 + x) + y*Log(y)/(1 - y))/(x - y)**2
        else
           fpbox=(-1 + x**2 - 2*x*Log(x))/(2*(-1 + x)**3*x)
        endif
        return
        end

C*************************************************************************
        double precision function gpbox(x,y)
        implicit double precision (a-h,o-z)
        if(abs(x-y).gt.1.d-3)then
           if(abs(y-1).gt.1.d-3)then
           gpbox = (3 - 4*x + x**2 + 4*x*Log(x) - 2*x**2*Log(x))/
     -    (2.d0*(-1 + x)**2*(-x + y)) - 
     -   (3*(-x + y)/2.d0 + x**2*Log(x)/(-1 + x) + y**2*Log(y)/(1 - y))/
     -    (x - y)**2
        else 
           gpbox= (-1 + x**2 - 2*x*Log(x))/(-1 + x)**3
        endif
        else
           if(abs(x-1).gt.1.d-3)then
              gpbox=(3/2 - 2*x + x**2/2 + Log(x))/(-1 + x)**3
           else
              gpbox=1.d0/3.d0
           endif
        endif
        return
        end
C*************************************************************************
C*************************************************************************
C========================================================================

        SUBROUTINE CHARGINO(TGB,CH_M,AMU,AMG,CHM,U,V,IERR)
C========================================================================
C  Input (all masses in GeV):
C    TGB       vev's ratio
C    CH_M      lightest chargino mass
C    AMU       MU parameter
C  Output:
C    AMG       gaugino mass
C    CHM(2)    chargino masses ( CHM(1) > CHM(2) )
C    U,V (2,2) chargino diagonalization matrices
C    IERR      1 (no solution for MU) 2 (divergent solution for MU)
C========================================================================
C       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        implicit double precision (a-h,o-z)
        DIMENSION U(2,2),V(2,2),CHM(2)
        WM=80.419D0
        IERR=0
        PI=4.*DATAN(1.D0)
        SQRT2=DSQRT(2.D0)
        EPS=1.D-2
C========================================================================
C       Definition of chargino's parameters
C========================================================================
        CBE=1./DSQRT(1.+TGB**2)
        SBE=TGB*CBE
        S2BE=2.*SBE*CBE
        C2BE=CBE**2-SBE**2

C*********** Chargino masses   *******************************************
        AMUQ=AMU**2
        CH_M2=CH_M**2
        WM2=WM**2
        XX=AMUQ-CH_M2
       DSQRT_MU=XX*(XX+2.*WM2)+(WM2*S2BE)**2
       IF(DSQRT_MU.LT.0.D+0) THEN
        IERR=1
        RETURN
       ENDIF
        IF(DABS(XX).LT.EPS**6) THEN
          IERR=2
          RETURN
        ENDIF
        IF(DABS(XX).LT.EPS) THEN
          AMG1=(WM2*S2BE-2.*AMUQ/S2BE)/(2.*AMU)
          AMG2=2.*AMU*WM2*S2BE/XX
        ELSE
          AMG1=(WM2*AMU*S2BE+CH_M*DSQRT(DSQRT_MU))/XX
          AMG2=(WM2*AMU*S2BE-CH_M*DSQRT(DSQRT_MU))/XX
        ENDIF
c NB here we can CHOOSE either the lighter or the heavier gaugino mass
c according to whether we prefer a light or heavy heavier chargino
c       however, choosing the max it is more likely that it will be positive
c as M2 should be (see notes)
        AMG=MAX(AMG1,AMG2)
c       AMG=MIN(AMG1,AMG2)      
        IF(AMG.LT.0.D0)THEN
         IERR=1
         RETURN
        ENDIF
        AMGQ=AMG**2
       CHM(2)=DSQRT((AMGQ+AMUQ+2.*WM2-DSQRT((AMGQ-AMUQ)**2+
     &  4.*WM2*(WM2*C2BE**2+AMGQ+AMUQ+2.*AMG*AMU*S2BE)))/2.)
       DIFF_CH=ABS(2.*(CHM(2)-CH_M)/(CHM(2)+CH_M))
       IF(DIFF_CH.GT.EPS**2) THEN
        AMG=MIN(AMG1,AMG2)
        IF(AMG.LT.0.D0)THEN
         IERR=1
         RETURN
        ENDIF
        AMGQ=AMG**2
       CHM(2)=DSQRT((AMGQ+AMUQ+2.*WM2-DSQRT((AMGQ-AMUQ)**2+
     &  4.*WM2*(WM2*C2BE**2+AMGQ+AMUQ+2.*AMG*AMU*S2BE)))/2.)
        DIFF_CH=ABS(2.*(CHM(2)-CH_M)/(CHM(2)+CH_M))
           IF(DIFF_CH.GT.EPS**2) THEN
            IERR=1
            RETURN
           ENDIF
       ENDIF
       CHM(1)=DSQRT((AMGQ+AMUQ+2.*WM2+DSQRT((AMGQ-AMUQ)**2+
     &  4.*WM2*(WM2*C2BE**2+AMGQ+AMUQ+2.*AMG*AMU*S2BE)))/2.)

C************* U,V matrices *****************************************
       HHH=(AMG*AMU-WM2*S2BE)
       IF(HHH.GT.0.D+0) THEN
       EEP=1.
       ELSE
       EEP=-1.
       ENDIF

        D_P=(AMG+AMU)**2+2.*WM2*(1.-S2BE)
        D_M=(AMG-AMU)**2+2.*WM2*(1.+S2BE)
       AA_P=(CHM(1)+EEP*CHM(2))*(AMG+AMU)/D_P
       AA_M=(CHM(1)-EEP*CHM(2))*(AMG-AMU)/D_M
       BB_P=WM*SQRT2*(CHM(1)+EEP*CHM(2))*(SBE-CBE)/D_P
       BB_M=WM*SQRT2*(CHM(1)-EEP*CHM(2))*(SBE+CBE)/D_M

       SQR=DSQRT((AA_P+AA_M)**2+(BB_P+BB_M)**2)
       CTETA=SQR/2.
       STETA=CTETA*(AA_P-AA_M)/(BB_P+BB_M)
       CPHI=(AA_P+AA_M)/SQR
       SPHI=(BB_P+BB_M)/SQR

       V(1,1)=CPHI
       V(1,2)=SPHI
       V(2,1)=-SPHI*EEP
       V(2,2)=CPHI*EEP
       U(1,1)=CTETA
       U(1,2)=STETA
       U(2,1)=-STETA
       U(2,2)=CTETA
C       ------------------------------------------------
C       Check
        CCK1=U(1,1)*AMG+U(1,2)*SQRT2*WM*CBE
        CCK2=U(2,1)*AMG+U(2,2)*SQRT2*WM*CBE
        CBK1=U(1,2)*AMU+U(1,1)*SQRT2*WM*SBE
        CBK2=U(2,2)*AMU+U(2,1)*SQRT2*WM*SBE
        C1=CCK1*V(1,1)+CBK1*V(1,2)-CHM(1)
        C2=CCK1*V(2,1)+CBK1*V(2,2)
        C3=CCK2*V(1,1)+CBK2*V(1,2)
        C4=CCK2*V(2,1)+CBK2*V(2,2)-CHM(2)
        IF(DABS(C1).GT.EPS)WRITE(6,*)'check failed c1=',c1
        IF(DABS(C2).GT.EPS)WRITE(6,*)'check failed c2=',c2
        IF(DABS(C3).GT.EPS)WRITE(6,*)'check failed c3=',c3
        IF(DABS(C4).GT.EPS)WRITE(6,*)'check failed c4=',c4
        RETURN
        END