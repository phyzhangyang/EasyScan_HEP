c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c  The program SuSpect: calculating the MSSM (or related models) Spectrum 
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: Abdelhak Djouadi, Jean-Loïc Kneur and Gilbert Moultaka
c  (main coding and maintenance: J-L Kneur)
c              (LPTA, CNRS & Universite de Montpellier II). 
c  VERSION 2.41  
c  Last modifications : August 1, 2008 by the SuSpect authors
c  Last modifications : November 20, 2008 by M.M.Muehlleitner.
c  The reference to be used for the program is: hep-ph/0211331 
c                 (published in Comput. Phys. Commun. 176:426) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  This code calculates the supersymmetric + Higgs particle spectrum in the
c  unconstrained Minimal Supersymmetric Standard Model (MSSM), as well as 
c  constrained models such as the minimal Supergravity (mSUGRA), the 
c  gauge mediated SUSY (GMSB) and anomaly mediated SUSY (AMSB) breaking 
c  models. All important features are included:
c  - Renormalization Group evolution between low and high energy scales.
c  - Consistent implementation of radiative electroweak symmetry breaking. 
c  - Calculation of the physical masses with radiative corrections.
c  (including e.g. all leading 2-loop corrections to the Higgs masses
c   in the DRbar scheme -courtesy of Pietro Slavich)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  For the users manual, updated information, changes, maintenance, see 
c  Home page: http://w3.lpta.univ-montp2.fr/~kneur/Suspect/
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
         SUBROUTINE SUSPECT2(iknowl,input,ichoice,errmess)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cc  This is the MAIN routine of the program, to be used as it is or to be
c  called by any other routine (such as SuSpect2_call.f),as discussed below.  
c  The routine has the following four important input control parameters:
c
c  IKNOWL: which sets the degree of control on the various parts of 
c          the algorithm. It has two possible values:
c  -- IKNOWL=0: blind use of the program, no warning and other messages.
c  default values are set for the control parameters and the program 
c  gives just the results from the physical input. 
c  -- IKNOWL=1: 
c  some warning/error messages are collected in the SuSpect2.out file
c  (this is the recommended choice).  
c
c  INPUT: sets input (and output) control, it offers now 4 possibilities:
c  =0: model and option input parameters ONLY read in file suspect2.in.
c  (output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c  =1: define yourself IN suspect2_call.f all relevant input and parameters.
c  (i.e. NO reading of input files): 
c  see example of  input in accompanying file suspect2_call.f
c  Maybe more convenient e.g. for scan over the model parameter space.
c  (output generated in both suspect2.out and SLHA format suspect2_lha.out) 
c  =2: same as input =0 but read SLHA format input file: suspect2_lha.in
c (it writes also all output in the SLHA format file: suspect2_lha.out)
c  =11: same as input=1, but NO output file(s) suspect*.out generated
c   
c  ICHOICE: initializes the various options for the models to be considered, 
c           the degree of accuracy to be required, the features to be 
c           included, etc. There are 10 possible choice at present and the
c           options are described in more details in the input file:
c !! NB: ICHOICE(..) superseded if using SLHA input file mode: in this
c   case we follow SLHA standards and conventions (now adapted to ver. 2).  
c  -- ICHOICE(1): Choice of the model to be considered.
c  -- ICHOICE(2): choice of perturbative order (1 or 2 loop) of the RGEs. 
c  -- ICHOICE(3): To impose or not the GUT scale. 
c  -- ICHOICE(4): Accuracy of the RGEs.
c  -- ICHOICE(5): To impose or not the radiative EWSB. 
c  -- ICHOICE(6): Chose different input in the Higgs sector (MA,MU,Mhu,Mhd)
c  -- ICHOICE(7): Choice of radiative corrections to (s)particles masses. 
c  -- ICHOICE(8): Choice of the EWSB scale
c  -- ICHOICE(9): Accuracy of the physical spectrum calculation. 
c  -- ICHOICE(10): different options for calculation of Higgs boson masses.
c  -- ICHOICE(11):  running/pole H masses used in loops.
c
c  ERRMESS: which provides a useful set of warning/error message flags,
c           that are automatically written in the output file SUSPECT2.OUT:
c  -- ERRMESS(i)= 0: Everything is fine.
c  -- ERRMESS(1)=-1: tachyon 3rd gen. sfermion from RGE
c  -- ERRMESS(2)=-1: tachyon 1,2 gen. sfermion from RGE
c  -- ERRMESS(3)=-1: tachyon A    (maybe temporary: see final mass) 
c  -- ERRMESS(4)=-1: tachyon 3rd gen. sfermion from mixing
c  -- ERRMESS(5)=-1: mu inconsistent (or unstable) after many iterations
c  -- ERRMESS(6)=-1: non-convergent mu from EWSB 
c  -- ERRMESS(7)=-1: EWSB maybe inconsistent  (!but RG-improved only check)
c  -- ERRMESS(8)=-1: V_Higgs maybe UFB or CCB (!but RG-improved only check)
c  -- ERRMESS(9)=-1: Higgs boson masses are NaN 
c  -- ERRMESS(10)=-1: RGE problems (non-pert and/or Landau poles)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  The program starts here:
c==========================
c
      implicit real*8(a-h,m,o-z)
      real*8 nf,nl,nq
      logical su_isNaN
      parameter(ni=87,nout=88,ninlha=85,noutlha=86)
      parameter(n=31)
      dimension ichoice(11),errmess(10),imod(1:2)
      dimension y(n),ysave(n),ygut(n),yewsb(n),ysav2(n)
      dimension U(2,2),VV(2,2),Z(4,4),dxmn(4)
      dimension gcen(2,2),gctb(2,2),glee(2,2),gltt(2,2),
     .          glbb(2,2),ghee(2,2),ghtt(2,2),ghbb(2,2)
      dimension ac1(2,2),ac2(2,2),ac3(2,2),an1(4,4),an2(4,4),an3(4,4),
     .          acnl(2,4),acnr(2,4)
      dimension gmn(4),xmn(4),gmc(2),gmst(2),msb(2),gmsl(2),
     .          gmsu(2),gmsd(2),gmse(2),gmsn(4)
      dimension bsgchm(2), ubsg(2,2),vbsg(2,2)
c
c     ****************************************************************
c     INPUT parameters for interface with other codes: 
c    
c NB: the parameters defined in the commons below in the INPUT/OUTPUT
c section are necessary (and sufficient in most situations!) for 
c interface with other codes.
c     ****************************************************************
c   "Standard model" INPUT parameters (couplings and fermion masses):
      COMMON/SU_smpar/dalfinv,dsw2,dalphas,dmt,dmb,dmc,dmtau
c !!NB: dmt,dmtau are pole masses but dmb is mb(mb)_MSbar 
c  RG evolution scale parameters EWSB scale, high and low RGE ends):
      COMMON/SU_rgscal/dqewsb,dehigh,delow
c   MSSM parameters of the scalar sector:     
      COMMON/SU_mssmhpar/dmhu2,dmhd2,dma,dmu
c   The U(1), SU(2), SU(3) soft SUSY-breaking gaugino masses:
      COMMON/SU_mssmgpar/dm1,dm2,dm3 
c   The soft-SUSY breaking slepton mass terms (3d and then 1/2 gen.): 
      COMMON/SU_mssmslep/dmsl,dmtaur,dmel,dmer
c   The soft-SUSY breaking squark mass terms (3d and then 1/2 gen.):
      COMMON/SU_mssmsqua/dmsq,dmtr,dmbr,dmuq,dmur,dmdr
c   The soft-SUSY breaking trilinear couplings (3d and then 1/2 gen.):
      COMMON/SU_atri3/dal,dau,dad
      COMMON/SU_atri12/dal1,dau1,dad1
c
c   GUT scale MSSM parameters output:     
      COMMON/SU_mssmhgut/mhu2gut,mhd2gut,magut,mugut
      COMMON/SU_mssmggut/m1gut,m2gut,m3gut 
      COMMON/SU_mssmslgut/mslgut,mtaurgut,melgut,mergut
      COMMON/SU_mssmsqgut/msqgut,mtrgut,mbrgut,muqgut,murgut,mdrgut
      COMMON/SU_A3gut/algut,augut,adgut
      COMMON/SU_A12gut/al1gut,au1gut,ad1gut
c  tan(beta) and sign(mu)
      COMMON/SU_radewsb/sgnmu0,tgbeta
c  mSUGRA case input parameters:
      COMMON/SU_msugra/m0,mhalf,a0
c  GMSB case input parameters:
      COMMON/SU_gmsb/mgmmess,mgmsusy,nl,nq
c  AMSB case input parameters:
      COMMON/SU_amsb/m32,am0,cq,cu,cd,cl,ce,chu,chd
c     ****************************************************************
c     COMMONS for OUTPUT masses and mixing angles:
c
c !! However some of the INPUT parameters above can also be OUTPUT
c at the end of the run: typically the soft terms like mu,mhu2, etc ..
c     ****************************************************************
      COMMON/SU_outhiggs/dml,dmh,dmch,alfa
c  light, heavy, charged Higgs masses, neutral (h,H) mix angle alpha 
      COMMON/SU_outginos/dmc1,dmc2,dmn1,dmn2,dmn3,dmn4,mgluino
c   charginos 1,2 masses, neutralinos 1-4 masses, gluino mass 
      COMMON/SU_outsqu/dmst1,dmst2,dmsu1,dmsu2
c  stop 1,2 and sup 1,2 = scharm 1,2 masses
      COMMON/SU_outsqd/dmsb1,dmsb2,dmsd1,dmsd2
c  sbottom 1,2 and sdown 1,2 = sstrange 1,2 masses
      COMMON/SU_outslep/dmsl1,dmsl2,dmse1,dmse2,dmsn1,dmsntau
c  stau 1,2 ; selectron (=smuon) 1,2; sneut_e,mu, sneut_tau masses
      COMMON/SU_outmix/thet,theb,thel
c  stop, sbottom, stau mixing angles
c low-energy contrained parameter values: rho-1, g_mu-2, Br(b->s gamma):
      COMMON/SU_lowen/crho,gmuon,brsg
c     ****************************************************************
c     COMMONs INTERNAL to the routine 
c
c ("internal" means that normally the user does not have to care about 
c any parameters defined by the commons etc below: in particular none 
c of these commons below should be necessary for interface with other codes)
c     ****************************************************************  
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
      COMMON/SU_stepwi/wistep,h1,kpole
      COMMON/SU_stegut/ifirst,jfirst,ygut
      COMMON/SU_errsf/sterr,sberr,stauerr,stnuerr
      COMMON/SU_qcdflag/nnlo,idrflag
      COMMON/SU_hflag/ihflag
      COMMON/SU_tachyrc/tachsqrc
      COMMON/SU_good/iflop
      COMMON/SU_sthresh/rmtop,susym,egut
      COMMON/SU_gunif/kunif
      COMMON/SU_param/gf,alpha,mz,mw
      COMMON/SU_cte/nf,cpi,zm,wm,tbeta
      COMMON/SU_als/xlambda,mc0,mb0,mt0,n0
      COMMON/SU_fmasses/mtau,mbpole,mtpole
      COMMON/SU_runmasses/mtaurun,mbrun,mtrun
      COMMON/SU_yuka/ytau,yb,ytop
      COMMON/SU_yukaewsb/ytauewsb,ybewsb,ytewsb,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_tbewsb/vuewsb,vdewsb 
       common/su_allewsb/yewsb
      COMMON/SU_treesfer/msbtr1,msbtr2,msttr1,msttr2
      common/su_mixflip/istflip,isbflip
      COMMON/SU_hmass/ma,ml,mh,mch,marun
      COMMON/SU_break/msl,mtaur,msq,mtr,mbr,al,au,ad,
     .            mu,m1,m2,m3
      COMMON/SU_break2/mel,mer,muq,mur,mdr
      COMMON/SU_smass/gmn,xmn,gmc,gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn
      COMMON/SU_hcoup/bcoup,a,gat,gab,glt,glb,ght,ghb,ghvv,glvv
      COMMON/SU_HMIX/BETA,Adum
      COMMON/SU_cplhsf/gcen,gctb,glee,gltt,glbb,ghee,ghtt,ghbb,
     .                 gatt,gabb,gaee
      COMMON/SU_cplhino/ac1,ac2,ac3,an1,an2,an3,acnl,acnr
      COMMON/SU_cteloop/vu,vd,atop,ab,atau,rmllt,rmllb,rmlltau,
     . rmrrt,rmrrb,rmrrtau
      COMMON/SU_soft/rmtaur,rml,rmbr,rmtr,rmq
      COMMON/SU_cpl/g12,g22,sw2
      COMMON/SU_sgnm123/sgnm1,sgnm2,sgnm3
      COMMON/SU_matino/u,vv,z,dxmn
      common/su_MAinput/piaa,tadba,D2MA,kMAflag  
c    (!! add for MA input case)
      common/su_errma/errma2z   ! added for ma^2(mz) <0 control
      common/su_mbmb/mbmb,imbmb   ! added for mb(mb) input
      common/su_nonpert/inonpert  ! added for non-pert problems
      COMMON/run_p/pizzp       
      COMMON/pietro/mApole,mCHpole  
      common/su_polemz/ipolemz
      COMMON/SU_renscale/scale
      common/SU_ftune/czmu,czbmu,ctmu,ctbmu
c
c     ****************************************************************
       external su_deriv1,su_deriv2,su_rkqc
       sgn(x)=dsign(x,x)/dabs(x)
c     ****************************************************************
c
c  Here comes the initialisation and reading part:
c  ==============================================
c  reinitialize various control parameters + other parameters:
           do ierr=1,10
           errmess(ierr)=0.d0
           enddo
          errnogo=0.d0
          irge=0
          iflop=0
          irpb=0
          tachsqrc=0.d0
          icount=0
          iremember=0
           do irg=1,31
           y(irg)=0.d0
           enddo
c       
       ml=0.d0
       ytauewsb=0.d0
       ybewsb=0.d0
       ytewsb=0.d0
       pizz_mz=0.d0
       pizzp =0.d0  
       inorc=0     
       inonpert=0  ! added for non-pert pbs control
       bup=0.d0  
       sterr=0.d0
       sberr=0.d0
       stauerr=0.d0
       stnuerr=0.d0
       errma2z=0.d0
c%%%  further reinitializations added
      alsewsb=0.d0
      g2ewsb=0.d0
      g1ewsb=0.d0
      vuewsb=0.d0
      vdewsb=0.d0 
c
      mh=0.d0
      mch=0.d0
      marun=0.d0
       mapole = 1.d0  ! initialization at 1st call (value later superseded)
c
      piaa=0.d0
      tadba=0.d0
       D2MA=0.d0
       kMAflag=0  
       imbmb  =0  
c
       dmc1=0.d0
       dmc2=0.d0
       dmn1=0.d0
       dmn2=0.d0
       dmn3=0.d0
       dmn4=0.d0
       mgluino=0.d0
c
       dmst1=0.d0
       dmst2=0.d0 
       dmsu1=0.d0 
       dmsu2=0.d0
       dmsb1=0.d0
       dmsb2=0.d0
       dmsd1=0.d0 
       dmsd2=0.d0
       dmsl1=0.d0
       dmsl2=0.d0
       dmse1=0.d0 
       dmse2=0.d0 
       dmsn1=0.d0 
       dmsntau=0.d0
c
       thet=0.d0
       theb=0.d0 
       thel=0.d0
c
       dml=0.d0
       dmh=0.d0
       dmch=0.d0
       alfa=0.d0 
c%%%%%%
c open OUTPUT file:
         if(input.ne.11) then
      OPEN(NOUT,FILE='suspect2.out',status='unknown')
         endif
c  Read input:
c  Physical input parameters:
      if(input.eq.0) then
c  Read all relevant physical parameters from the input file:)
      OPEN(NI,FILE='suspect2.in',status='unknown')
      do i=1,10
      read(ni,*)
      enddo
      read(ni,*) ichoice(1)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(2)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(3)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(4)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(5)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(6)
c
      do i=1,5
      read(ni,*)
      enddo
      read(ni,*) ichoice(7)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) ichoice(8)
c
      read(ni,*)
      read(ni,*)
      read(ni,*) ichoice(9)
c
      do i=1,5
      read(ni,*)
      enddo
      read(ni,*) ichoice(10)
c
      do i=1,5
      read(ni,*)
      enddo
      read(ni,*) ichoice(11)
c
      do i=1,3
      read(ni,*)
      enddo
      read(ni,*) alfinv, alphas, mt, mbmb, mtau   !now mb(mb)
c
      read(ni,*)
      read(ni,*)
      read(ni,*) Ehigh,  qewsb
c
      read(ni,*)
      read(ni,*)
      read(ni,*)
c   Minimal supergravity input 
      if(ichoice(1).eq.10) then
      read(ni,*) rm0,  rmhalf,  a0,  tgbeta,  sgnmu0
      m0=rm0
      mhalf=rmhalf
      else if(ichoice(1).eq.11) then
c  Gauge Mediated Supersymmetry Breaking  input: 
      do i=1,4
      read(ni,*)
      enddo
      read(ni,*) mgmmess, mgmsusy, tgbeta,  sgnmu0, nl, nq
      else if(ichoice(1).eq.12) then
c  Anomaly Mediated Supersymmetry Breaking  input: 
      do i=1,8
      read(ni,*)
      enddo
      read(ni,*) m32,am0,tgbeta,sgnmu0,cq,cu,cd,cl,ce,chu,chd
      else
c  i.e. non-universal arbitrary input case
      do i=1,12
      read(ni,*)
      enddo
      read(ni,*) mhu2,  mhd2, tgbeta, sgnmu0
      read(ni,*)
      read(ni,*) m1,  m2,  m3
      read(ni,*)
      read(ni,*) msl,  mtaur, msq,  mtr,  mbr
      read(ni,*) 
      read(ni,*) mel, mer, muq,  mur, mdr
      read(ni,*) 
      read(ni,*) al, au, ad, al1, au1, ad1
      read(ni,*)
      read(ni,*)
      read(ni,*) ma, mu
      mu = sgnmu0*dabs(mu)     ! add to avoid inconsistent user's input
      endif
      close(ni)
       if(gf.eq.0d0) gf = 1.16639d-5  ! only if not already defined
       if(mz.eq.0d0) mz = 91.187d0  ! only if not already defined
c  
c simple renaming for SLHA output purpose:
      dalfinv = alfinv 
      dalphas = alphas
      dmt = mt
      dmb = mbmb
      dmc =mc
      dmtau =mtau
      dqewsb= qewsb
      dehigh = ehigh
      dMHU2     = mhu2
      dMHD2     = mhd2
      dM1       = m1
      dM2       = m2
      dM3       = m3
      dMSL      = msl
      dMTAUR    = mtaur
      dMSQ      = msq
      dMTR      = mtr
      dMBR      = mbr
      dMEL      = mel
      dMER      = mer
      dMUQ      = muq
      dMUR      = mur
      dMDR      = mdr
      dAL       = al
      dAU       = au
      dAD       = ad
      dAL1      = al1
      dAU1      = au1
      dAD1      = ad1
      dMU       = mu
      dma = ma
      call SU_read_leshouches(input,ninlha,ichoice,imod)
c 
       else if(input.eq.2) then 
c this is the case where input file is read in SLHA format:
        open(ninlha,FILE='suspect2_lha.in',status='unknown')
        call SU_read_leshouches(input,ninlha,ichoice,imod)  
        close(ninlha)
      endif 
c      
      if(ichoice(1).eq.10.or.ichoice(1).eq.1) then
      igut_in=1
      else
      igut_in=0
      endif
      if(input.ne.2.and.igut_in.eq.0) goto 8765 ! all cases =/= SLHA input file
c this is the case when input is read from suspect2.in or from a calling
c file for mSUGRA        
c nb this series of IF is due to the various ways in which input may be defined
      if(input.eq.0) then
      m0=rm0
      mhalf=rmhalf
      endif
      if(input.ne.2.and.ichoice(1).eq.10) then
      al0=a0
      ad0=a0
      au0=a0
      al10=a0
      ad10=a0
      au10=a0
      mhu2=m0**2
      mhd2=m0**2
      mtaur0=m0
      msl0=m0
      mbr0=m0
      mtr0=m0
      msq0=m0
      mer0=m0
      mel0=m0
      mdr0=m0
      mur0=m0
      muq0=m0
      m10=mhalf
      m20=mhalf
      m30=mhalf
      endif
 8765 continue
c      endif
c
      if(input.ne.0) then  
      alfinv   =dalfinv
cc      sw2      =dsw2
      alphas   =dalphas
c
      mt       =dmt
      mbmb       = dmb      ! mbmb is mb(mb)_MSbar input
      mtau       = dmtau
c
      qewsb    = dqewsb
      ehigh    = dehigh
c
      if(ichoice(1).eq.10) then
      rm0      = m0
      else if(ichoice(1).eq.12) then
      rm0 = am0
      endif
      rmhalf   = mhalf
c (nb parameters a0 and sgnmu already defined via common)
c
      if(input.eq.2.or.input.eq.11.or.ichoice(1).eq.1) then
        if(input.eq.11.and.ichoice(1).eq.10) goto 8766
      mhu2     = dmhu2
      mhd2     = dmhd2
c
      m1       = dm1
      m2       = dm2
      m3       = dm3
c
      msl      = dmsl
      mtaur    = dmtaur
      msq      = dmsq
      mtr      = dmtr
      mbr      = dmbr
      mel      = dmel
      mer      = dmer
      muq      = dmuq
      mur      = dmur
      mdr      = dmdr
c
      al       = dal
      au       = dau
      ad       = dad
c 
      al1       = dal1
      au1       = dau1
      ad1       = dad1
c
      ma       = dma
      mu       = dmu      
 8766 continue
      endif
        if(sgnmu0.ne.1d0.and.sgnmu0.ne.-1d0) then
           if(mu.ne.0d0) then 
	   sgnmu0 = mu/dabs(mu)
           else 
	   sgnmu0=1d0
	   endif
	else
      mu = sgnmu0*dabs(mu)     !added to avoid inconsistent user's input
        endif
      if(input.eq.1) call SU_read_leshouches(input,ninlha,ichoice,imod)  
      endif
c
c  ! added reinitialization of mhu2,mhd2 for scans:
      if(ichoice(1).ne.2.and.ichoice(6).eq.0) then
      mhu2=1.d4
      mhd2=1.d4
c (nb ichoice(6)=0 -> MA,MU input thus these mhu2,mhd2 values are
c irrelevant but only initialized for convergence of iteration control)
      endif
      ihflag= ichoice(10)
      ihrcsave=ihflag
      ipolemz=ichoice(11)
      tgbet0 = tgbeta
      beta_z = datan(tgbet0)  
        if(ichoice(1).eq.12) rm0 = am0 
c  blind use: assign protection default values to control parameters:
        if(ichoice(2).eq.0) ichoice(2)=11
c  (i.e. 1-loop RGE at first run by default, if 2-loop not chosen) 
c:essai        if(ichoice(3).eq.0.and.ichoice(1).ne.2) ichoice(3)=1
c (i.e. GUT scale always consitently recalculated as g1(gut)=g2(gut)
        if(ichoice(4).ne.1.and.ichoice(4).ne.2) ichoice(4)=1
        iaccsave=ichoice(4)
c (i.e. protections against wrong or undefined rge accuracy setup)
        if(ichoice(5).ne.1) ichoice(5)=1
c maggie changed 2/12/2008
c (i.e. always EWSB)
        if(input.eq.0.and.ichoice(1).eq.10) ichoice(6)=1  ! added 29/11/08
c (protection against wrong use of MA input if mSUGRA in old input file)
c end maggie changed 2/12/2008
        if(ichoice(8).ne.0) ichoice(8) = 1
c (i.e to be sure that not choosing default EWSB scale 
c  =(m_t_L*m_t_R)^(1/2) is on purpose) 
c choose frozen scale in RGE parameters:
        kpole = ichoice(8)
c for some NO RGE purposes:
        inorge = ichoice(1)
c for special case where MA(pole) is really input:
        kmaflag=ichoice(6)
c choose susy R.C. options:
        if(ichoice(7).eq.1) then
c only mt,mb,mtau susy R.C.
        isfrc = 0
        else if(ichoice(7).eq.2) then
c  mt,mb,mtau  + (all) squarks + (all) gaugino susy RC:
        isfrc = 1
        endif
c
c   optimize number of long (RG+ Full spectrum) iterations
      irgmax = 50
      irgsave=irgmax 
c
      if (ichoice(1) .le. 2 ) then
c one could have arbitrary m1,m2,m3 signs
      sgnm1  = m1/dabs(m1)
      sgnm2  = m2/dabs(m2)
      sgnm3  = m3/dabs(m3)
      else if(ichoice(1) .eq. 10 .or. ichoice(1).eq.11) then
c msugra (or gmsb) case
      sgnm1  = 1.d0
      sgnm2  = 1.d0
      sgnm3  = 1.d0
      else if(ichoice(1) .eq. 12) then
c amsb case
      sgnm1  = 1.d0
      sgnm2  = 1.d0
      sgnm3  = -1.d0
      endif
c
c  rename input parameters for SLHA  and other input choice purpose
cc      if(input.eq.2.or.input.eq.0) then
      if(input.ne.2.and.ichoice(1).eq.10) goto 9876 
      mhu20 = mhu2
      mhd20 = mhd2
c
      msl0=msl
      mtaur0=mtaur
      msq0=msq
      mtr0=mtr
      mbr0=mbr
c 
      mel0=mel
      mer0=mer
      muq0=muq
      mur0=mur
      mdr0=mdr
c
      al0=al
      au0=au
      ad0=ad
      al10=al1
      au10=au1
      ad10=ad1
      mu0=mu
c
      m10=m1
      m20=m2
      m30=m3
 9876 continue
cc      endif
      beta = datan(tgbeta)
c
c some other basic parameter definitions
       pi = 4*datan(1.d0)
       cpi = 1.d0/(16*pi*pi)
       if(gf.eq.0d0) gf = 1.16639d-5  ! only if not already defined
       sw2 = .2221d0   ! only starting value! sw2_DRbar calculated below
       fermi=gf
       if(mz.eq.0d0) mz = 91.187d0  ! only if not already defined
       zm = mz
c  guess starting point for susym , elow, ehigh scales:
       elow = mz
       if(ichoice(3).ne.0) ehigh = 1.d17
       if(ichoice(1).eq.10.or.ichoice(1).eq.12) then
       susym = .5*(rm0+rmhalf) +mz
       else if(ichoice(1).eq.11) then
       susym = mz 
       rm0 = susym
c NB this "rm0" is not m0, only used at 1rst iter to guess mu(0),b(0) 
c
       else if(ichoice(1).eq.1) then
       rm0= (msl+mtaur+msq+mtr+mbr+mel+mer+muq+mur+mdr)/10
       susym = .5*(rm0+(m1+m2+m3)/3) +mz
       endif 
       gut=ehigh
       kunif=ichoice(3)
       wistep= 1.d2
       nf = 6.d0
c
c  mw, sw^2 msbar at Z scale (values may be changed):
       cw2= 1.d0-sw2
       sw=dsqrt(sw2)
       cw =dsqrt(cw2)
       mw = mz*cw
       wm = mw
       mc=1.40d0
       rmtau=mtau
       rmtau2=rmtau**2
c
c  Some saving 
      mc0=mc
c      mb0=mb
      mb0=4.9d0    ! value just used for very first initialization
      mb=mb0
      mt0=mt
      mbpole=mb
      mtpole=mt
      mtaurun=mtau
      mbrun=mb0
      mtrun=mt0
c (initial values only! at 1st iter mrun = mpole)

c passing from alpha_S(MZ) MSbar to alpha_S(MZ) DRbar:      
      alphas0=alphas
      g32=4*pi*alphas0/(1.d0-(1.d0/2)*alphas0/(2*pi) )
      alphas=g32/4/pi
c (NB value in fact used at first RG only, does not include SUSY etc R.C.)
c
c passing from alpha(MZ) MSbar to alpha(MZ) DRbar:      
       alpha =1.d0/(alfinv -1.d0/pi/6)   
c (NB value used at first RG iteration only, does not include SUSY etc R.C.)
       e2=4*pi*alpha
       sw20=sw2
       cw20=1.d0-sw20
       g12= e2/cw20
       g22=e2/sw20
       g120=g12
       g220=g22
c
      acc=1.d-8
      nloop=2 
      nnlo = 1
      idrflag =0
      xlambda=xitla(nloop,alphas0,acc)   ! alphas0 is MSbar
      n0 = 5
      CALL alsini(acc)
      imbmb=0   !   (just reinitialization) 
      rmbms=runm(mz,5)
      mb=mbpole
      mb0=mb
c rmbms is mb running(MZ) in MSbar scheme
      mc=runm(mz,4)
c
c Now defining running quark masses in DRbar at Z scale:
       rmb = rmbms*(1.d0-alphas/(3*pi) -23*alphas**2/(72*pi**2) 
     . +3*g22/(128*pi**2) +13*g12/(1152*pi**2))
c rmb is mb running(MZ) in DRbar scheme (what is mostly used after)
c
c      xlambda=xitla(nloop,alphas0,acc)
c      CALL alsini(acc)
c
       rmb2=rmb**2
       rmtop = runm(mt,6)
       rmt2=rmtop**2
      if(ichoice(1).eq.2.and.ichoice(6).eq.1) then
      iremember=1
      ichoice(1)=0   ! a trick to simplify the bottom-up case
      endif
c
c     ****************************************************************
c     Long iteration (on RGE + spectrum once defined) starts here:
c     ****************************************************************
c
 44    irge=irge+1
c ! reinitialize at each RGE iter Higgs RC choice (1 or 2 loop):
       ihflag=ihrcsave
c reinitialize  error messages until last iteration:
       if(irge.le.irgmax) then
       do i=1,10
       errmess(i) =0.d0
       enddo 
       endif
c
       tbeta=tgbet0
c 
c calculating s^2_W_DRbar(MZ), g1_DRbar(MZ), g2_DRbar(MZ) incl. SUSY R.C.:
       if(irge.ge.2) then
c (because at first call no susy physical masses etc are defined)
c first need to compute PIzz(Q=mz), PIww(Q=mz):
       scale = mz
       call SU_PIXX(sw2,dsqrt(g22),dsqrt(g12),tbeta,pizz,piww,piww0
     $      ,0d0)               ! PiXX with pole mt 
       pizz_mz=pizz
c   Now the more complete calculation of g1,g2,sw2 (MZ) in DRbar:
      call su_runningcp(alphas0,mt,rmtop,m3z,tbeta,pizz,piww,piww0,
     .      alphadr,alphas,sw2) 
       e2=4*pi*alphadr
       cw2=1.d0-sw2
c!!!following redef of sw etc added
       sw=dsqrt(sw2)
       cw =dsqrt(cw2)
       mw = sqrt((mz**2+pizz)*cw2 - piww) 
       wm = mw
       g12= e2/cw2
       g22=e2/sw2
       g32=4*pi*alphas         
       endif

c - higgs vev at Z scale: tbeta = vu/vd
c (NB in our normalization MZ = (g12+g22)/2*(vu2+vd2), and there
c  are no factors of sqrt(2) in the Higgs doublet components
c (cf Ramond et al PRD49(1994) 4882)

       if(irge.eq.1) then       
          pizz =0.d0
       else          
          call SU_PIXX(sw2,dsqrt(g22),dsqrt(g12),tbeta,pizz,piww,piww0
     $         ,rmtop)          ! PiZZ with running mt 
          pizz_mz=pizz
       endif
c
      if(su_isNaN(pizz).or.mz**2+pizz.le.0.d0) then
c !!! protections added
c non-pert or NaN pb, uses tree-level values temporarily:
      pizz = 0.d0
      if(irge.eq.irgmax) inonpert=-1    
      endif
c
       vd2 = 2*(mz**2+pizz)/(g12+g22)/(1.d0+tbeta**2)
       cbeta= 1.d0/dsqrt(1.d0+tbeta**2)
       sbeta=tbeta*cbeta
       vu2 = vd2*tbeta**2
       vd= dsqrt(vd2)
       vu= dsqrt(vu2)
       v= dsqrt(vu2+vd2)
       vd_mz=vd                 
       vu_mz=vu
c
c  defining Yukawa couplings at Z scale:
       if(irge.eq.1) then
       y(4) = mtau/vd
c QCD corrections to mt(mz) (yt(mz)=y(6)) in DRbar including Logs:
        mtlog = dlog((mt/mz)**2)
        delmt = alphas/pi*(5.d0/3 -mtlog)
     . +alphas**2*(0.876d0 -0.384*mtlog +0.038*mtlog**2) 
c       
       y(6) = mt/vu*(1.d0-delmt)
       y(5) = rmb/vd
       endif
c - higgs vev at Z scale: y(7)=Ln vu, y(8)=Ln vd
       y(7) = .5*log(vu2)
       y(8) = .5*log(vd2)
c      1st stage: evolution of gauge + yukawa cpl from Mz to GUT:
c ! for irge=1 (iter. 1) yukawa's determined from QCD corrections only
       y(1) = 5.d0*g12/3.d0
c (i.e usual SU(5) normalisation of g1)
       y(2) = g22
       y(3) = g32
c set RGE accuracy choices (3 different)
       if(ichoice(4).eq.0) then
       h1=.2d0
       eps=1.d-3
       else if(ichoice(4).eq.1) then
       h1=.06d0
       eps=1.d-3
       else if(ichoice(4).eq.2) then
       h1=.01d0
       if(ichoice(1).eq.0.or.ichoice(1).eq.2) h1=.005d0  
c     ! more precise rge for pmssm
       eps=2.d-5
       endif
c
         if(ichoice(3).ne.0.and.irge.eq.1) ehigh =1.d17
c note ehigh = 1.e17 will be superseded 
c by true unification scale (where y(1)=y(2) by def.)):
c
       if(ichoice(1).eq.0) then
c Case where only mass spectrum at EWSB scale is calculated:
c it is then assumed that all MSSM parameters are defined at EWSB scale,
c except tanbeta(mz). The EWSB scale is an input arbitrarily chosen, and
c the only RGE performed is to calculate the gauge+yukawa +vevs from their
c input values at mz scale to their consistent values at EWSB scale.
c
         if(ichoice(8).eq.0) then
       if(qewsb.eq.0.d0) qewsb = 1.05*zm
c (protections in case of badly chosen ewsb scale input in this case)
         else
           if(irge.eq.1) then
           qewsb=dsqrt(msq*mtr)
           else
           qewsb=dsqrt(msttr1*msttr2)
           endif
          if(qewsb.lt.mz) qewsb=mz+1.d-1   !! added protection 
         endif
       x1 = dlog(zm)
       x2 = dlog(qewsb)
       else
c means all other cases where RGE is performed from mz to GUT scales
      if(ichoice(8).eq.1) then
        if(irge.eq.1) then
        qewsb= mz
        else
        qewsb = dsqrt(msttr1*msttr2)
        if(qewsb.lt.mz) qewsb=mz+1.d-1    !! added protection
        endif  
      endif 
       x1 = dlog(zm)
       x2 = dlog(1d20)       
c!essai       if(ichoice(1).eq.2.and.ehigh.ne.0d0) x2=dlog(ehigh)
       if(ichoice(3).eq.0.and.ehigh.ne.0d0) x2=dlog(ehigh)
       endif
c
        ifirst=0
        jfirst=0
        scale = qewsb

c  first step: run from mz to high scale with initial conditions
c g_i(mz), yukawa_i, find GUT scale etc.
c
       if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv1,su_rkqc)
       else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv2,su_rkqc)
       endif
c protection against RGE num. pbs (Landau poles, non-perturbativity):
       if(iflop.eq.1) then
       errmess(10)=-1.d0
       goto 801
       endif
       if(ichoice(1).eq.0) then
       g1ewsb = dsqrt(3*y(1)/5)
       g2ewsb = dsqrt(y(2))
       alsewsb = y(3)/4/pi
       ytauewsb=y(4)
       ybewsb= y(5)
       ytewsb= y(6)       
       vuewsb=dexp(y(7))
       vdewsb=dexp(y(8))
       tbeta= vuewsb/vdewsb
       goto 880
       endif

c        
c (exact) gauge (g1=g2) unif. if required:
 882   if(egut.ne.0.d0.and.ichoice(3).ne.0) then
       ehigh=egut
         do irg=1,31
         y(irg)=ygut(irg)
         enddo
       y(2)=y(1)
       endif
c
       do i=1,8
       ysave(i)=y(i)
       end do
       vu = dexp(y(7))
       vd = dexp(y(8)) 
       mtaugut=vd*y(4)
       mbgut = vd*y(5)
       mtgut=vu*y(6)
c
       if(ichoice(1).eq.2.and.irge.eq.irgmax) then
       dmhu2 =y(12)
       dmhd2 =y(13)
       dmtaur = dsqrt(y(14))
       dmsL = dsqrt(y(15))
       dmbr =dsqrt(y(16))
       dmtr =dsqrt(y(17))
       dmsQ =dsqrt(y(18))
       dmer =dsqrt(y(24))
       dmel =dsqrt(y(25))
       dmdr =dsqrt(y(26))
       dmur =dsqrt(y(27))
       dmuQ =dsqrt(y(28))
       dal =y(9)
       dad =y(10)
       dau =y(11)
       dal1 =y(29)
       dad1 =y(30)
       dau1 =y(31)
c
       dB = y(19)
       dmu = sgnmu0*dexp(y(23))
       dM1=sgnM1*dexp(y(20))
       dM2=sgnM2*dexp(y(21))
       dM3=sgnM3*dexp(y(22))
       goto 801
       endif
c******************************************************************
c      2d stage: evolution from HIGH (GUT) scale down to low energy
c******************************************************************
c
c Now taking input rmu0,B0 values (!only guess initialization values)
       if(ichoice(6).eq.0.and.irge.gt.1) then
       mhu2 = y(12)
       mhd2=y(13)
c i.e. for MA_pole,MU(EWSB) input: in this case mhu2,mhd2(GUT) should not 
c  be reinitialized (except at first RGE iteration where they are undefined)
       endif
c guess mu(GUT) value at first time run (later superseeded by EWSB MU) 
c this apply in particular in mSUGRA or non-univ cases:
       if(rm0.eq.0d0) rm0 = 1d-4   ! protection
       if(irge.eq.1) rmu0 = 1.1*rm0
c
   7   rmu0=sgnmu0*dabs(rmu0)
c also guess value for b(GUT):
       b0 = 2*rm0
c set up boundary conditions at GUT scale:
c   yukawa coupling (eventual) unification at GUT scale:
       if(ichoice(3).ge.2) then
       y(5)=y(4)
       ysave(5)=y(5)
       endif
c  (tau- b unification)
       if(ichoice(3).eq.3) then
       y(6)=y(5)
       ysave(6)=y(6)
       endif
c  (tau-b-top unification): ! caution then tanbeta is constrained!
c (!! NOT YET operative  !!)
c
c - Higgs initial vev at GUT scale: fixed from evolution
c   from Z scale (see above)
  2    icount=icount+1
c icount is a counter for things to be done only at first iteration
       errhu=1.d3
       errhd=1.d3
       ifix =0
c errhu,hd,ifix are convergence control parameters for ichoice(6)=0
c if on different high-energy input (msugra, amsb,gmsb,..) starts here:
c
  77     if(ichoice(1).eq.1.or.ichoice(1).eq.10) then
c  Unconstrained MSSM: general case; ! now includes SUGRA case with 
c  universality as special case in same algorithm
       y(9)=al0
       y(10)=ad0
       y(11)=au0
       y(29)=al10
       y(30)=ad10
       y(31)=au10
       y(12)=mhu2
       y(13)=mhd2
       y(14)=mtaur0**2
       y(15)=msl0**2
       y(16)=mbr0**2
       y(17)=mtr0**2
       y(18)=msq0**2
       if(irge.eq.1) y(19)=b0
       if(irge.eq.1) y(23)=dlog(dabs(rmu0))
       y(24)=mer0**2
       y(25)=mel0**2
       y(26)=mdr0**2
       y(27)=mur0**2
       y(28)=muq0**2
       y(20)=dlog(dabs(m10))
       y(21)=dlog(dabs(m20))
       y(22)=dlog(dabs(m30))
c                 
       else if(ichoice(1).eq.2.and.irge.eq.1) then
c Bottom-up RGE case with soft(EWSB) input: initialize reasonable GUT
c scale (guess) values at firt RGE iteration to catch convergence
       do j= 9,11
       y(j) = 0d0
       end do
       do j= 29,31
       y(j) = 0d0
       end do
       do kk=12,18
       y(kk) = 100d0
       end do
       do kk=24,28
       y(kk) = 100d0
       end do
       y(19) = b0
       do l=20,22
       y(l) = dlog(100d0)
       end do
       y(23) = 0d0
c           
       else if(ichoice(1).eq.12) then
c AMSB case:  m3/2, m0,c_q, etc (coeffs of m0),sgn(mu0) input)
      CALL SU_AMSBSUB(rm0,m32,cq,cu,cd,cl,ce,chu,chd,y(1),y(2),y(3),
     . y(4),y(5),y(6),y(9),y(10),y(11),y(29),y(30),y(31),y(12),y(13),
     . y(14),y(15),y(16),y(17),y(18),y(24),y(25),y(26),y(27),y(28),
     . m10,m20,m30)
c
       y(20)=dlog(dabs(m10))
       y(21)=dlog(dabs(m20))
       y(22)=dlog(dabs(m30))
c remaining needed parameters:
       if(irge.eq.1) y(19)=b0
       if(irge.eq.1) y(23) = dlog(dabs(rmu0)) 
c forces radiative EW breaking (if was not chosen before:)
       ichoice(5)=1
       endif
c=======================================
 883    x1= dlog(ehigh)
c
c  Generic end of running scale: determined "consistently" by default:
c  - at MZ scale for gauge+yukawas couplings, that are defined at MZ,
c  - at EWSB scale (!by default = sqrt(mst_L*mst_R)).
c  For all others RG parameters 
c
       if(scale.eq.0d0) scale =mz+1d-1    ! protection against undefined
       xewsb=dlog(scale)
       x2=dlog(mz)
       h1=-h1
 3     issb=0
       istab=0
       ifirst=0
       jfirst=0
           if(ichoice(1).ne.11) then
c  RGE is made in two steps from Gut scale to EWSB; then MZ
       if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,x1,xewsb,eps,h1,1.d-8,nok,nbad,
     . su_deriv1,su_rkqc)
       else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,x1,xewsb,eps,h1,1.d-8,nok,nbad,
     . su_deriv2,su_rkqc)
       endif
           else
c this is GMSB case: input are mgmmess,mgmgsusy,nl,nq, sgn(mu) and tbeta)
c but intermediate (messenger) scale: mgmmess for RGE of soft terms
c
       if(irge.eq.1) y(19) = B0
       if(irge.eq.1) y(23) = dlog(dabs(rmu0)) 
        xint = dlog(mgmmess)
       if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,x1,xint,eps,h1,1.d-8,nok,nbad,su_deriv1,
     .              su_rkqc)
       else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,x1,xint,eps,h1,1.d-8,nok,nbad,su_deriv2,
     .               su_rkqc)
       endif

c - Now input soft susy-breaking terms at messenger scale:
        CALL SU_GMSBSUB(mgmmess,mgmsusy,nl,nq, y(1),y(2),y(3),
     . y(9),y(10),y(11),y(29),y(30),y(31),y(12),y(13),y(14),y(15),y(16),
     . y(17),y(18),y(24),y(25),y(26),y(27),y(28),m10,m20,m30)
c
       y(20)=dlog(dabs(M10))
       y(21)=dlog(dabs(m20))
       y(22)=dlog(dabs(m30))
       do i=29,31
       y(i) = 0.d0
       enddo

c next RGE down to EWSB scale: forces as usual radiative EW breaking:
       ichoice(5)=1
       if(ichoice(2).eq.11) then
c$$$      CALL SU_ODEINT(y,n,xint,x2,eps,h1,1.d-8,nok,nbad,su_deriv1,
      CALL SU_ODEINT(y,n,xint,xewsb,eps,h1,1.d-8,nok,nbad,su_deriv1,
     .              su_rkqc)
       else if(ichoice(2).eq.21) then
c$$$      CALL SU_ODEINT(y,n,xint,x2,eps,h1,1.d-8,nok,nbad,su_deriv2,
      CALL SU_ODEINT(y,n,xint,xewsb,eps,h1,1.d-8,nok,nbad,su_deriv2,
     .              su_rkqc)
       endif
           endif
c (last endif = end of non-univ/mSUGRA/GMSB/AMSB distinctions)
c    
c protection against big troubles in RGE (e.g. Landau poles):
       if(iflop.eq.1) then
       errmess(10)=-1.d0
       goto 801
       endif
c    
c check for problems (non-perturbativity or/and Landau poles) in RGE:
       do i=1,31
       if(su_isnan(y(i))) then
       errmess(10)=-1.d0
       endif
       end do
       if(errmess(10).eq.-1.d0) then
       goto 801
       endif
       if(ichoice(1).eq.2) then
c  ! new algorithm for EWSB soft terms input with bottom-up RGE
c  Unconstrained MSSM: general case: note in this case al0 etc are
c supposed to be soft terms input values at low EWSB input scale
       y(9)=al0
       y(10)=ad0
       y(11)=au0
       y(29)=al10
       y(30)=ad10
       y(31)=au10
       if(ichoice(6).eq.1.or.irge.eq.1) then
       y(12)=mhu2
       y(13)=mhd2
       endif
c (nb otherwise it means that MA_pole and MU(EWSB) are input)        
       y(14)=mtaur0**2
       y(15)=msl0**2
       y(16)=mbr0**2
       y(17)=mtr0**2
       y(18)=msq0**2
       if(irge.eq.1) y(19)=b0
       if(irge.eq.1) y(23)=dlog(dabs(rmu0))
       y(24)=mer0**2
       y(25)=mel0**2
       y(26)=mdr0**2
       y(27)=mur0**2
       y(28)=muq0**2
       y(20)=dlog(dabs(m10))
       y(21)=dlog(dabs(m20))
       y(22)=dlog(dabs(m30))
c      

       endif
       vu=dexp(y(7))
       vd=dexp(y(8))
c
c saving all rge parameters at ewsb scale:
       do ip=1,31
       yewsb(ip)=y(ip)
       enddo
c saving also Yukawas and others at EWSB scale:
 886   ytauewsb=y(4)
       ybewsb=y(5)
       ytewsb=y(6)
       alsewsb=y(3)/(4*pi)
       g2ewsb= dsqrt(y(2))
       g1ewsb= dsqrt(3*y(1)/5)
       vuewsb = dexp(y(7))
       vdewsb=dexp(y(8))
       tbeta=dexp(y(7))/dexp(y(8))
       atau =y(9)
       ab=y(10)
       atop=y(11)
       al1 =y(29)
       ad1=y(30)
       au1=y(31)       
       rmhu2 = y(12)
       rmhd2 = y(13)
c  !! change (after 10 iter) of standard fixed point algorithm:
c  mhu_new = mhu_ewsb -> mhu_new = (1-c)*mhu_old + c*mhu_ewsb, c=0.3
c  this trick may help recovering convergence if on the boarder:
c  also increasing RGE accuracy in this case:
           if(irge.ge.10) then
       rmhu2= .7d0*rmhu2old +.3d0*rmhu2
       ichoice(4)=2 
c  (i.e. also increasing RGE accuracy in this case)
           endif
c
       do kk=14,18
       if(y(kk).lt.0.d0) then  
       if(irge.eq.irgmax) errmess(1)=-1.d0   
            if(iknowl.eq.2) then
       write(*,'(a)') 'Bad input: one m^2(3rd gen. sf) <0 from RGE '
       write(*,'(a)') 'maybe temporary due to iteration. wait and see'
            endif
c
       endif
       enddo
       do kk=24,28
       if(y(kk).lt.0.d0) then
       if(irge.eq.irgmax) errmess(2)=-1.d0  
             if(iknowl.eq.2) then
       write(*,'(a)') 'Bad input: one m^2(1,2 gen. sf) <0 from RGE '
       write(*,'(a)') 'maybe temporary due to iteration. wait and see'
             endif
       endif
       enddo
c 
       if(errmess(1).eq.-1.d0.or.errmess(2).eq.-1.d0) then
       goto 801
       endif
c
     
       rmtaur = dsqrt(y(14))
       rmL = dsqrt(y(15))
       rmbr=dsqrt(y(16))
       rmtr=dsqrt(y(17))
       rmQ =dsqrt(y(18))
       rmer =dsqrt(y(24))
       rmel=dsqrt(y(25))
       rmdr=dsqrt(y(26))
       rmur=dsqrt(y(27))
       rmuQ=dsqrt(y(28))
c!! modif (temporary, until final conv) rescue in case tachyon RGE sf 
       if(irge.lt.irgmax) then    ! protections against NaN 
       if(y(14).lt.0.d0) rmtaur=1.d0
       if(y(15).lt.0.d0) rmL=1.d0
       if(y(16).lt.0.d0) rmbr=1.d0
       if(y(17).lt.0.d0) rmtr=1.d0
       if(y(18).lt.0.d0) rmQ=1.d0
c
       if(y(24).lt.0.d0) rmer=1.d0
       if(y(25).lt.0.d0) rmel=1.d0
       if(y(26).lt.0.d0) rmdr=1.d0
       if(y(27).lt.0.d0) rmur=1.d0
       if(y(28).lt.0.d0) rmuQ=1.d0
       else
       if(y(14).lt.0.d0) errmess(1)=-1d0
       if(y(15).lt.0.d0)  errmess(1)=-1d0
       if(y(16).lt.0.d0)  errmess(1)=-1d0
       if(y(17).lt.0.d0)  errmess(1)=-1d0
       if(y(18).lt.0.d0)  errmess(1)=-1d0
c
       if(y(24).lt.0.d0) errmess(2)=-1d0
       if(y(25).lt.0.d0) errmess(2)=-1d0
       if(y(26).lt.0.d0) errmess(2)=-1d0
       if(y(27).lt.0.d0) errmess(2)=-1d0
       if(y(28).lt.0.d0) errmess(2)=-1d0
c
       if(errmess(1).eq.-1d0.or.errmess(2).eq.-1.d0) goto 801
       endif
c
       B = y(19)
       if(ichoice(1).eq.2.and.ichoice(6).eq.0) then    ! if MU(ewsb) input
       rmu0=mu   
       y(23)=dlog(dabs(mu))
       endif 
       rmu = sgn(rmu0)*dexp(y(23))
       if(ichoice(6).eq.0) rmu=mu   ! if MU(ewsb) input
       if(ichoice(5).eq.1) then
       Bold  =B
       rmuold =1.d0
c
       else
c  means no radiative EW required
       continue
       endif
       rmino1=sgnM1*dexp(y(20))
       rmino2=sgnM2*dexp(y(21))
       rmino3=sgnM3*dexp(y(22))
c
c  interface with Higgs mass spectrum calculations:
      ihflag=ichoice(10)
c 
      msl=rml
      mtaur=rmtaur
      msq=rmq
      mtr=rmtr
      mbr=rmbr
c
      mel=rmel
      mer=rmer
      muq=rmuq
      mur=rmur
      mdr=rmdr
c
      al=atau
      au=atop
      ad=ab
      mu=rmu
c
c 
      m1=rmino1
      m2=rmino2
      m3=rmino3

c  tan beta at the relevant (ewsb) scale: tbeta
c (extracted from COMMON/SU_tbewsb/vuewsb,vdewsb )
      tbeta= vuewsb/vdewsb
      beta = datan(tbeta)
c
c set EWSB scale IF not defaut value (used in RGE, V_eff, and susy R.C):
      if(ichoice(8).eq.0) then
      scale= qewsb 
      ewsb2= qewsb**2
      else if(ichoice(8).eq.1) then
c Default EWSB scale:
         if(msttr1.ne.0.d0.and.irge.gt.1) then
      scale= dsqrt(msttr1*msttr2)
         else
      scale = dsqrt(msq*mtr)
         endif
      if(scale.lt.mz) then   !!  added protections 
      scale=mz+1.d-1
        if(scale.lt.dsqrt(msq*mtr)) then
        scale=   dsqrt(msq*mtr)
        qewsb= scale
        endif
      endif
      ewsb2= scale**2
      endif
c
c    Gaugino masses:
       CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
c
       dmc1=gmc(1)
       dmc2=gmc(2)
       dmn1=xmn(1)
       dmn2=xmn(2)
       dmn3=xmn(3)
       dmn4=xmn(4)
c       
c******************************************************************
c- Set up the conditions for radiative sym. break. and stability:
c*******************************************************************
       cbeta=1.d0/dsqrt(1.d0+tbeta*tbeta)
       sbeta = tbeta*cbeta
       c2beta = cbeta*cbeta-sbeta*sbeta
       wm2=wm*wm
       zm2=zm*zm
c
                 if(ichoice(6).eq.0) then
c========================================
c  input is MA_pole!, MU(EWSB). Consistent M^2_Hu, M^2_d from EWSB
c  with iteration. 
 66      ifix=ifix+1
       inonpert=0
       if(ichoice(1).eq.0) then
       ewsb2= qewsb**2
       ytewsb = rmtop/vu
       endif
c   Gaugino masses
        if(ichoice(7).eq.2.and.irge.eq.irgmax) inorc =1
        CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
c
c  sfermion masses
       CALL SU_SFERMION(msq,mtr,mbr,msl,mtaur,al,au,ad,mu,
     .                 gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn)
         if(sterr.eq.-1.d0.or.sberr.eq.-1.d0.or.stauerr.eq.-1.d0.
     . or.stnuerr.eq.-1.d0) then
c means there is really the tachyonic sfermion mass problem:
c can't even calculate Higgs spectrum etc, so has to stop really.
         errmess(4)=-1.d0
         goto 801
         endif
c
c   Higgs masses
c the user's Ma input value is used in Higgs spectrum calc.:
        ma2 = ma**2
        mapole2=ma2     !! in that case the input is really MA_pole
       if(ewsb2.lt.mz**2) ewsb2=qewsb**2
      CALL SU_SUSYCP(tgbeta)      
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
      alfa= a
c
c  call one-loop effective potential corrections to Mh^2_1,2:
c  dVdvd2, dVdvu2 are d(V_eff)/d(vd^2) and d(V_eff)/d(vu^2) which
c  add corrections to m^2_Phid (rmhd2) and m^2_Phiu (rmhu2) resp. 
       rmtaur = mtaur
       rml = msl
       rmbr= mbr
       rmtr= mtr
       rmq = msq
       atau= al
       ab= ad
       atop = au
       rmst12= msttr1**2
       rmst22= msttr2**2
       rmsb12= msbtr1**2
       rmsb22= msbtr2**2
       rmstau12=gmsl(1)**2
       rmstau22=gmsl(2)**2
       dmsu1=gmsu(1)
       dmsu2=gmsu(2)
       dmsd1=gmsd(1)
       dmsd2=gmsd(2)
       dmse1=gmse(1)
       dmse2=gmse(2)
       dmsn1=gmsn(1)
       dmsntau=gmsn(3)
       dmc1=gmc(1)
       dmc2=gmc(2)
       dmn1=xmn(1)
       dmn2=xmn(2)
       dmn3=xmn(3)
       dmn4=xmn(4)
       rmu=mu
       ewsb2 = scale**2
c   
       if(ytewsb.eq.0.d0) ytewsb=rmtop/vu
       if(ytauewsb.eq.0.d0) ytauewsb=rmtau/vu
       if(ybewsb.eq.0.d0) ybewsb= rmb/vu
       CALL SU_VLOOP2(ewsb2,MU,AU,AD,AL,dVdvd2,dVdvu2,pizz)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
          if(ifix.eq.1) then
          dVdvd2=0.d0
          dVdvu2=0.d0
          rmhu2old=0.d0
          rmhd2old=0.d0
          endif
c
c   
      sb2=dsin(beta)**2
      cb2=dcos(beta)**2
      mzdr2= mz**2+pizz
      madr2= mapole2 +piaa -tadba -D2MA
      rmhu2 = (cb2-sb2)*mzdr2/2 +cb2*madr2 -mu**2   -dVdvu2 
      rmhd2 = (sb2-cb2)*mzdr2/2 +sb2*madr2 -mu**2   -dVdvd2 
c
c (Note: -dVdvi2: to get "tree-level" values of M^2_Hu, M^2_Hd, 
c  thus without V_eff loop corrections)
      dmhu2=rmhu2   
      dmhd2=rmhd2
      B=(rmhd2+rmhu2+dVdvd2+dVdvu2+2*rmu**2)*dsin(2*beta)/(2*rmu)
c - to be compared to previous M^2_Hu,Hd values:
       errhuold=errhu
       errhdold=errhd
       errhu = (rmhu2-rmhu2old)/rmhu2
       errhd = (rmhd2-rmhd2old)/rmhd2
       errstop=1d-2
       if(ifix.gt.1) errstop=1d-4
       if(dabs(errhu).gt.errstop.and.dabs(errhd).gt.errstop) then
       rmhu2old=rmhu2
       rmhd2old=rmhd2
       goto 66
       endif
c
       y(12) = rmhu2 
       y(13) = rmhd2
c
 8     continue
       if(ichoice(6).eq.0) then
c stop long (RGE) iterations on spectrum when xx % accuracy reached:
c (usually needs ~ 3-4 iterations). NB conv. test is made on MA_run(EWSB)
        if(irge.eq.1) madr2old=0.d0
        if(ichoice(9).le.1) then                 
        if(dabs(1.d0-madr2old/madr2).lt.2d-2) irgmax=irge
        else                                     
        if(dabs(1.d0-madr2old/madr2).lt.2d-4) irgmax=irge
        endif
        madr2old=madr2    
        endif
c
c========================== Now comes ichoice(6).neq.0 i.e:
c      input parameters M_Hu, M_Hd
c      consistent mu, B from EWSB conditions
                        else
c===========================                        
c stop long iterations on spectrum when xx % accuracy reached:
c (usually needs ~ 3-4 iterations)
        if(ichoice(1).ne.2) then
        if(irge.eq.1) rmhu2old=0.d0
        if(ichoice(9).le.1) then                 ! 1% accuracy
        if(dabs(1.d0-rmhu2old/rmhu2).lt.2d-2) irgmax=irge
        else                                     ! 0.01% accuracy
        if(dabs(1.d0-rmhu2old/rmhu2).lt.2d-4) irgmax=irge
        endif
        rmhu2old=rmhu2 
        if(irge.eq.irgsave) errmess(5)=-1.d0
        endif
c --- algorithm to find a consistent mu with V_eff corrections:
       errB=1.d5
       errmu=1.d5
       if(ichoice(5).eq.1) then
c i.e. we want EWSB to determine mu and B
       dVdvd2=0.d0
       dVdvu2=0.d0
       ifix=0
  80   ifix=ifix+1
       inonpert=0
       mu=rmu
       CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
c
       dmc1=gmc(1)
       dmc2=gmc(2)
       dmn1=xmn(1)
       dmn2=xmn(2)
       dmn3=xmn(3)
       dmn4=xmn(4)
c 
c equation for MA_run:
        if(irge.eq.1) pizz=0.d0
        ma2 =(rmhu2+dVdvu2-rmhd2-dVdvd2)/dcos(2*beta)-zm**2-pizz 
        dmhu2=rmhu2
        dmhd2=rmhd2
c     
        if(ma2.ge.0.d0) then
        MA=dsqrt(ma2)
        masave=ma
        errmess(3)=0.d0
        else
c Allows for temporary MA^2 < 0 (before EWSB converges)
c and attempt to retrieve a correct MA via a correct MU etc.
c Gives approximate MA_run(ewsb) values just so that calculation 
c (EWSB iteration) can proceed for a while:
        ma=1.1d0
        if(.NOT.(su_isNaN(mapole)).and.mapole.ne.0.d0) ma=mapole   
        masave=ma
       CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)       
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
        endif
c
c   Now Calculate sfermion masses and mixing angle:
c
       CALL SU_SFERMION(msq,mtr,mbr,msl,mtaur,al,au,ad,mu,
     .                  gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn)
         if(sterr.eq.-1.d0.or.sberr.eq.-1.d0.or.stauerr.eq.-1.d0.
     . or.stnuerr.eq.-1.d0) then
c means there is really the tachyonic sfermion mass problem:
c can't even calculate Higgs spectrum etc, so has to stop really.
         errmess(4)=-1.d0
             if(iknowl.eq.2) then
         write(*,'(a)')' CAUTION: m^2_sf < 0! . Has been changed to 0 '
             endif
         goto 801
         endif
         if(tachsqrc.eq.-1.d0) then
          errmess(4)=-1.d0
         goto 801
        endif
c Otherwise (= no tachyonic sfermions) calculate Higgs mass spectrum:
       CALL SU_SUSYCP(tbeta) 
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
c protection against NAN Higgs that could occurs despite previous protec.
       if(su_isNAN(ml).or.su_isNAN(mH).or.su_isNAN(MCH)) then
       errmess(9)=-1.d0
       goto 801
       endif
       if(ml.eq.0.d0.or.ml.gt.1.d10.or.mH.gt.1.d10) then
       if(irge.eq.irgmax) then
       errmess(9)=-1.d0
       goto 801
       endif
       endif
c
       rmst12= msttr1**2
       rmst22= msttr2**2
       rmsb12= msbtr1**2
       rmsb22= msbtr2**2
       rmstau12=gmsl(1)**2
       rmstau22=gmsl(2)**2
       dmsu1=gmsu(1)
       dmsu2=gmsu(2)
       dmsd1=gmsd(1)
       dmsd2=gmsd(2)
       dmse1=gmse(1)
       dmse2=gmse(2)
       dmsn1=gmsn(1)
       dmsntau=gmsn(3)
       alfa= a
c  call one-loop effective potential corrections to Mh^2_1,2:
c 
c  dVdvd2, dVdvu2 are d(V_eff)/d(vd^2) and d(V_eff)/d(vu^2) which
c  add corrections to m^2_Hd (rmhd2) and m^2_Hu (rmhu2) resp. 
c   
       if(ytewsb.eq.0.d0) ytewsb=rmtop/vu
       if(ytauewsb.eq.0.d0) ytauewsb=rmtau/vu
       if(ybewsb.eq.0.d0) ybewsb= rmb/vu
       CALL SU_VLOOP2(ewsb2,MU,AU,AD,AL,dVdvd2,dVdvu2,pizz)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
c
              if(su_isNAN(dvdvd2).or.su_isNaN(dvdvu2)) then
              if(irge.eq.irgmax.and.ifix.ne.1) then
              errmess(3)=-1.d0
              goto 801
              else
c Maybe due to uncorrect spectrum at 1rst iter., give it a chance
              if(su_isNAN(dvdvd2)) dvdvd2=0.d0
              if(su_isNAN(dvdvu2)) dvdvu2=0.d0
              endif  
              endif

c Now the radiative breaking conditions DEFINE true mu(mz):
c 
c Tree-level EWSB conditions as (first time!) MU guess:
       if(ifix.eq.1) then
       rmu2 =(rmhd2-(rmhu2)*tbeta**2)/(tbeta**2-1.d0)
     .        -(zm**2+pizz)/2.d0 
       else
       rmu2 =(rmhd2+dVdvd2-(rmhu2+dVdvu2)*tbeta**2)/(tbeta**2-1.d0)
     .        -(zm**2+pizz)/2.d0
       endif
       if(rmu2.le.0.d0) then
            if(iknowl.eq.2) then
       write(*,'(a)')'Warning: MU^2(EWSB) <0 (may be temporary) '
             endif    
       if(irge.eq.irgmax.and.ifix.ge.5) then 
c Consider the MU^2 < 0 problem irremediable:
       errmess(6)=-1.d0
       goto 801
       else
c Take approximate MU "RG" =f(MA,M_Hu,Mhd) to attempt to retrieve
c EWSB convergence:
         if(ma**2-rmhu2 -rmhd2.gt.0.d0) then
       rmu= sgn(rmu0)*dsqrt((ma**2-rmhu2 -rmhd2)/2)
         else
c take arbitrary small MU to attempt to retrieve EWSB convergence:
        rmu = sgn(rmu0)*10d0
         endif
       endif
c       rmu= rmu/2.d0
      else
       rmu =sgn(rmu0)*dsqrt(rmu2)
c  !! added: change (after 10 iter) of standard fixed point algorithm:
c   mu_new = mu_ewsb -> mu_new = (1-c)* mu_old + c*mu_ewsb, c=0.3
c   to try recovering convergence if on the boarder:
       if(ifix.ge.10) rmu= .7d0*rmuold +.3d0*rmu 
       MU=rmu
       endif
c      
c - ..and true B(EWSB):
c  tree-level EWSB conditions as first time MU guess:
       if(ifix.eq.1) then
       B = (rmhd2 +rmhu2 +2*rmu**2)*sbeta*cbeta/rmu 
       else
       B = (rmhu2+dVdvu2 +rmhd2+dVdvd2 +2*rmu**2)*sbeta*cbeta/rmu  
       endif
c
c - to be compared to evolved mu values:
       errmuold=errmu
       errmu= (rmu-rmuold)/rmuold

       if(dabs(errmu).lt.5.d-5.and.ma2.gt.0.d0.and.rmu2.gt.0.d0
     & .or.ifix.eq.20) then
c i.e. considers as unconvergent MU from EWSB either:
c   -inaccurate (> 1e-4) convergence;
c   - more than 5 tolerated iterations IF MA^2 was in fact <0, 
c so that convergence is around fake MA,MU
c since MA was articifially forced = MZ temporarily in that case
c  
       goto 81
       else
       if(ma2.le.0.d0.and.ifix.eq.5) goto 81
c !!added to get out if really unconvergent EWSB:
       if(dabs(errmu).gt.dabs(errmuold).and.ifix.gt.15) then  
         if(irge.eq.irgmax) then
       errmess(6)=-1.d0
       goto 801
         else
       goto 81
         endif
       endif
c
       rmuold=rmu
       goto 80
       endif
c  ( end of the iterative loop on consistent MU,B )
c
  81    continue 
                     endif
c (previous endif = end of the choice M_Hu,MHd or MA,MU input)
       if(ma2.le.0.d0.and.ifix.eq.5.and.irge.eq.irgmax) then
       errmess(6)=-1.d0
       errmess(3)=-1.d0
          if(iknowl.eq.2) then
       write(*,'(a)')' consistent EWSB unconvergent below 1d-4' 
          endif
       endif
       if(ifix.eq.20.and.irge.eq.irgmax) then
       errmess(6)=-1.d0
       errmess(3)=-1.d0
          if(iknowl.eq.2) then
       write(*,'(a)')' consistent EWSB unconvergent below 1d-4' 
          endif
       endif
c
 88    if(ichoice(1).eq.1.and.ifix.eq.20) then
       errmess(6)=-1.d0
          if(iknowl.eq.2) then
       write(*,'(a)')' consistent EWSB unconvergent below 1d-4' 
          endif
       endif
c
       endif
c    control SSB V stability (naive RG improved checks of UFB/CCB):
      r1= rmhd2 +dvdvd2 +mu**2
      r2= rmhu2 +dvdvu2 +mu**2
      r3= B*mu
      test1= r1*r2-r3*r3
      test2 = ma2 +2*r3
      test3 = ma2 -2*r3
      if(ichoice(5).eq.1) then
      if(test1.ge.0.d0.and.irge.eq.irgmax) then
      errmess(7)=-1.d0
            if(iknowl.eq.2) then
      write(*,'(a)')' Warning!: EW Sym. Break may be not realized '
            endif
      endif
      if(test2.lt.0.d0.or.test3.lt.0.d0.and.irge.eq.irgmax) then
      errmess(8)=-1.d0
           if(iknowl.eq.2) then
      write(*,'(a)')' Warning: Potential maybe unbounded from below  '
            endif
      endif
c CCB (simplest!) constraints, checked at EWSB scale:
            if(irge.eq.irgmax) then
          ccbt= atop**2-3*(msq**2 +mtr**2 +rmhu2 +rmu**2)
          ccbb= ab**2-3*(msq**2 +mbr**2 +rmhd2 +rmu**2)
          ccbtau= atau**2-3*(msl**2 +mtaur**2 +rmhd2 +rmu**2)
          ccbu= au1**2-3*(muq**2 +mur**2 +rmhu2 +rmu**2)
          ccbd= ad1**2-3*(muq**2 +mdr**2 +rmhd2 +rmu**2)
          ccbl= al1**2-3*(mel**2 +mer**2 +rmhd2 +rmu**2)
          if(ccbt.gt.0.d0.or.ccbb.gt.0.d0.or.ccbtau.gt.0.d0) then
c ! these are points which do not pass those naive CCB constraints
          errmess(8)=-1.d0
          endif
          if(ccbu.gt.0.d0.or.ccbd.gt.0.d0.or.ccbl.gt.0.d0) then
          errmess(8)=-1.d0
          endif
            endif
       else
c Means no radiative EW required
c Now B = y(19) and mu =exp(y(23)) are determined from EW breaking
c (however not radiative breaking in this case)
       rmu2 =(rmhd2+dVdvd2-(rmhu2+dVdvu2)*tgbeta**2)/(tgbeta**2-1.d0)
     .        -zm**2/2.d0
       if(rmu2.le.0.d0) then
            if(iknowl.eq.2) then
       write(*,'(a)')' Warning: initial rmu0(HIGH) inconsistent. '
       write(*,'(a)')' has been changed '
            endif
       rmu0  = rmu0/2
       do i=1,8
       y(i)=ysave(i)
       end do
       x2=x1
       h1=-h1
       goto 7
       else
       continue 
       endif
c      
       rmu =sgn(rmu0)*dsqrt(rmu2)
c - .. B(mz):
       B = (rmhd2+dVdvd2+rmhu2+dVdvu2+2*rmu2)*sbeta*cbeta/rmu 
c    control of SSB and V stability scales:
      r1= rmhd2 +dVdvd2 +rmu2
      r2= rmhu2 +dVdvu2 +rmu2
      r3= B*rmu
c
      test1= r1*r2-r3*r3
      test2= r1+r2+2*r3
      test3=r1+r2-2*r3
       if(test1.gt.0.d0) then
       errmess(7)=-1.d0
             if(iknowl.eq.2) then
       write(*,'(a)')'Warning: m^2(Hu),m^2(Hd) inconsistent with EWSB' 
       write(*,13) rmhu2,rmhd2
       write(*,'(a)')' have been changed '
             endif       
       mhu2= 1.5*mhu2
       mhd2= mhd2
       do i=1,8
       y(i)=ysave(i)
       end do
       x2=x1
       h1=-h1
       goto 7 
       endif
       if(test2.lt.0.d0.or.test3.lt.0.d0) then
       errmess(8)=-1.d0
            if(iknowl.eq.2) then
       write(*,'(a)')' Warning: Potential unbounded from below! ' 
       write(*,'(a)')' m^2(Hu),m^2(Hd) values been changed ' 
             endif
       mhu2 = 1.5*mhu2
       mhd2  =mhd2
       do i=1,8
       y(i)=ysave(i)
       end do
       x2=x1
       h1=-h1
       goto 7 
       endif
       endif
       if(ichoice(5).ne.1) then
       ma = dsqrt(rmhu2 +rmhd2 +2*rmu**2 )
c  calculate Higgs mass spectrum
       if(ewsb2.lt.mz**2) ewsb2 = qewsb**2
       if(ytewsb.eq.0.d0) ytewsb=rmtop/vu
       CALL SU_SUSYCP(tgbeta)       
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
c
c  calculate sfermion masses and mixing angle:

       CALL SU_SFERMION(msq,mtr,mbr,msl,mtaur,al,au,ad,mu,
     .                  gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn)
        endif
c
c Special case of unconstrained MSSM with low-en input:
c =====================================================

 880   if(ichoice(1).eq.0) then
c  case of the pMSSM (unconstrained MSSM, low-en input) 
c
                                  if(ichoice(6).eq.0) then
c   Input parameter of the pMSSM  is MA_pole , MU(EWSB)
c stop long iterations on spectrum when xx % accuracy reached:
c (usually needs ~ 3-4 iterations)
        if(irge.eq.1) mhu2old=0.d0
           if(irge.lt.irgmax) then
        if(ichoice(9).le.1) then
        if(dabs(1.d0-mhu2old/mhu2).lt.2d-2) irgmax=irge  ! 1% accuracy
        else
        if(dabs(1.d0-mhu2old/mhu2).lt.2d-4) irgmax=irge ! .01%
        endif
           endif
        mhu2old=mhu2 
        if(irge.eq.irgsave) errmess(5)=-1.d0
c   Gaugino masses
  881   beta =datan(tbeta)
        if(ichoice(7).eq.2.and.irge.eq.irgmax) inorc =1
        CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
c
c  sfermion masses
       CALL SU_SFERMION(msq,mtr,mbr,msl,mtaur,al,au,ad,mu,
     .                 gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn)
         if(sterr.eq.-1.d0.or.sberr.eq.-1.d0.or.stauerr.eq.-1.d0.
     . or.stnuerr.eq.-1.d0) then
c means there is really the tachyonic sfermion mass problem:
c can't even calculate Higgs spectrum etc, so has to stop really.
         errmess(4)=-1.d0
         goto 801
         endif
c
c   Higgs masses
        ma2 = ma**2
        mapole2=ma2     !! in that case the input is really Ma_pole
       if(ewsb2.lt.mz**2) ewsb2=qewsb**2
      CALL SU_SUSYCP(tgbeta)      
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
      alfa= a
c Check of EWSB in this parametrization:
c Note we include in the EWSB consistency relations all the 
c V_eff contributions +loop: indeed, it is consistent with the fact
c that all Higgs masses are calculated with 1- +2-loop contributions:
c  
       rmtaur = mtaur
       rml = msl
       rmbr= mbr
       rmtr= mtr
       rmq = msq
       atau= al
       ab= ad
       atop = au
c
       rmst12= msttr1**2
       rmst22= msttr2**2
       rmsb12= msbtr1**2
       rmsb22= msbtr2**2
       rmstau12=gmsl(1)**2
       rmstau22=gmsl(2)**2
       rmt2=mtrun**2
       rmtop=mtrun
       rmb2=mbrun**2
       rmtau2= mtaurun**2
       dmsu1=gmsu(1)
       dmsu2=gmsu(2)
       dmsd1=gmsd(1)
       dmsd2=gmsd(2)
       dmse1=gmse(1)
       dmse2=gmse(2)
       dmsn1=gmsn(1)
       dmsntau=gmsn(3)
       alfa = a
       dmc1=gmc(1)
       dmc2=gmc(2)
       dmn1=xmn(1)
       dmn2=xmn(2)
       dmn3=xmn(3)
       dmn4=xmn(4)
       rmu=mu
       ewsb2 = scale**2
       CALL SU_VLOOP2(ewsb2,MU,AU,AD,AL,dVdvd2,dVdvu2,pizz)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
c   
      sb2=dsin(beta)**2
      cb2=dcos(beta)**2
      mzdr2= mz**2+pizz
      madr2= mapole2 +piaa -tadba -D2MA
      rmhu2 = (cb2-sb2)*mzdr2/2 +cb2*madr2 -mu**2   -dVdvu2 
      rmhd2 = (sb2-cb2)*mzdr2/2 +sb2*madr2 -mu**2   -dVdvd2 
c
c (Note: -dVdvi2: to get "tree-level" values of M^2_Hu, M^2_Hd, 
c  thus without V_eff loop corrections)
      dmhu2=rmhu2   
      dmhd2=rmhd2
      B=(rmhd2+rmhu2+dVdvd2+dVdvu2+2*rmu**2)*dsin(2*beta)/(2*rmu)
c
c  Control of SSB and V stability scales:
      r1= rmhd2 +dVdvd2  +rmu**2
      r2= rmhu2 +dVdvu2  +rmu**2
      r3= B*rmu
      test1= r1*r2-r3*r3
      test2= r1+r2+2*r3
      test3=r1+r2-2*r3
c
      mhu2=rmhu2
      mhd2=rmhd2
         if(ichoice(1).eq.0.or.ichoice(1).eq.2) then
        rmino1=m1
        rmino2=m2
        rmino3=m3
         sgnm1 = m1/dabs(m1)
         sgnm2 = m2/dabs(m2)
         sgnm3 = m3/dabs(m3)
          else
       rmino1=sgnM1*dexp(y(20))
       rmino2=sgnM2*dexp(y(21))
       rmino3=sgnM3*dexp(y(22))
         endif
        if(test1.ge.0.d0) then
      errmess(7)=-1.d0
        endif
        if(test2.lt.0.d0.or.test3.lt.0.d0) then
       errmess(8)=-1.d0
        endif
      if(ichoice(1).eq.2.and.irge.eq.irgmax) goto 801
c       
c========================== Endif of ichoice(1) pMSSM
                       else
c==========================
c   Input parameter of pMSSM  is MHd2,MHu2:    
      ichoice(5) = 1
       rmhd2=mhd2
       rmhu2=mhu2
        if(irge.eq.1) then
        pizz=0.d0
        dvdvd2=0.d0
        dvdvu2=0.d0
        endif 
       ewsb2=qewsb**2
       if(irge.ge.2) then
       CALL SU_VLOOP2(ewsb2,MU,AU,AD,AL,dVdvd2,dVdvu2,pizz)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
c
       endif
        ma2 =(rmhu2+dVdvu2-rmhd2-dVdvd2)/dcos(2*beta)-zm**2-pizz
c
c --- Algorithm to find a consistent MU with V_eff corrections:
c --- the radiative breaking conditions DEFINE true mu(mz):
c 
       tgbeta=tbeta
       rmu2 =(rmhd2+dVdvd2-(rmhu2+dVdvu2)*tbeta**2)/(tbeta**2-1.d0)
     .        -(zm**2+pizz)/2.d0
           if(rmu2.le.0.d0) then
               if(iknowl.eq.2) then
       write(*,'(a)')' CAUTION: initial M^2_Hu,Hd inconsistent'
       write(*,'(a)')' their values were changed so that mu^2 >=0! ' 
c
               endif
c  find the minimal values of M^2_Hu,Hd to guarantee mu^2 >0,MA>0:
      rmhu2 = (1.d-6+mz**2/2)*(1.d0-tgbeta**2)/(1.d0+tgbeta**2) +
     . (ma**2-2*1.d-6)/(1.d0+tgbeta**2)    
      rmhd2 = -rmhu2    
      rmu = sgnmu0*1.d-6
      rmu2=rmu**2
            else
c rmu^2 >0 from the input
      rmu = sgnmu0*dsqrt(rmu2)
      rmu2=rmu**2
            endif
c stop long iterations on spectrum when xx % accuracy reached:
c (usually needs ~ 3-4 iterations)
        if(irge.eq.1) rmu2old=0.d0
        if(ichoice(9).le.1) then
        if(dabs(1.d0-rmu2old/rmu2).lt.0.02d0) irgmax=irge
        else
        if(dabs(1.d0-rmu2old/rmu2).lt.0.0002d0) irgmax=irge  !.002->.0002
        endif
        rmu2old=rmu2 
        if(irge.eq.irgsave) errmess(5)=-1.d0
c      
c - .. and true B(mz):
       b = (rmhd2+dvdvd2 +rmhu2+dVdvu2 +2*rmu2)*sbeta*cbeta/rmu 
       mu=rmu
       CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then
       errmess(10)=-1.d0
       goto 801
       endif
c  calculate sfermion masses and mixing angle:
       CALL SU_SFERMION(msq,mtr,mbr,msl,mtaur,al,au,ad,mu,
     .               gmst,msb,gmsl,gmsu,gmsd,gmse,gmsn)
         if(sterr.eq.-1.d0.or.sberr.eq.-1.d0.or.stauerr.eq.-1.d0.
     . or.stnuerr.eq.-1.d0) then
c means there is really the tachyonic sfermion mass problem:
c can't even calculate Higgs spectrum etc, so has to stop really.
         errmess(4)=-1.d0
         goto 801
         endif
c      
       if(ma2.ge.0.d0) then
       ma = dsqrt(ma2 )
       else
       ma = 1.d-6
       errmess(3)=-1.d0
       endif
       if(ewsb2.lt.mz**2) ewsb2=qewsb**2
       CALL SU_SUSYCP(tgbeta)
       if(inonpert.eq.-1.and.irge.eq.irgmax) then 
       errmess(10) =-1.d0
       goto 801
       endif
       alfa = a       
c
c
c    control of SSB and V stability scales:
      r1= rmhd2 +dVdvd2 +rmu2
      r2= rmhu2 +dVdvu2 +rmu2
      r3= B*rmu
      test1= r1*r2-r3*r3
      test2= r1+r2+2*r3
      test3=r1+r2-2*r3
  
            if(test1.ge.0.d0) then
           errmess(7)=-1.d0
            endif
            if(test2.lt.0.d0.or.test3.lt.0.d0) then
       errmess(8)=-1.d0
            endif
         if(ichoice(1).eq.0.or.ichoice(1).eq.2) then
        rmino1=m1
        rmino2=m2
        rmino3=m3
         sgnm1 = m1/dabs(m1)
         sgnm2 = m2/dabs(m2)
         sgnm3 = m3/dabs(m3)
         else
       rmino1=sgnM1*dexp(y(20))
       rmino2=sgnM2*dexp(y(21))
       rmino3=sgnM3*dexp(y(22))
         endif
c
        endif
        endif
c===================================
       if(irge.eq.irgmax.and.iremember.eq.1) ichoice(1)=2 
c (a trick to simplify the case of bottom-up evol with Mhu,Mhd input)       
       if(ichoice(1).eq.2.and.irge.eq.irgmax) then
c
c unconstrained MSSM runned up to High scale 
c  Now the final run from Q_ewsb to Q_final:
       x1 = dlog(qewsb)
       x2=dlog(ehigh)
       h1= dsign(h1,x2-x1)
c
       if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv1,
     .              su_rkqc)
       else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv2,
     .              su_rkqc)
       endif
       goto 882
       endif
c
c     ****************************************************************
c      SUSY radiative corrections to tau,b,t and sparticle masses:
c     ****************************************************************
c recovering all rge parameter values at mz scale:
 884      if(ichoice(1).ne.0.and.irge.ge.2) then
       y(19)=b
       y(23)=dlog(dabs(mu))
       xewsb=dlog(qewsb)    !! added to be consistent with new
c                           !! protections for tachyonic sfermions
        if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,xewsb,x2,eps,h1,1.d-8,nok,nbad,
     . su_deriv1,su_rkqc)
        else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,xewsb,x2,eps,h1,1.d-8,nok,nbad,
     . su_deriv2,su_rkqc)
        endif

        vu = vu_mz               
        vd = vd_mz

       rmtau=y(4)*vd
       rmb = y(5)*vd
       rmtop =y(6)*vu
       mtaurun = rmtau
       mbrun = rmb
       mtrun = rmtop
       else if(ichoice(1).eq.0) then 
       y(9)=al
       y(10)=ad
       y(11)=au
       y(29)=al1
       y(30)=ad1
       y(31)=au1
       y(12)=mhu2
       y(13)=mhd2
       y(14)=mtaur**2
       y(15)=msl**2
       y(16)=mbr**2
       y(17)=mtr**2
       y(18)=msq**2
       y(19)=b
       y(23)=dlog(dabs(mu))
       y(24)=mer**2
       y(25)=mel**2
       y(26)=mdr**2
       y(27)=mur**2
       y(28)=muq**2
       y(20)=dlog(dabs(m1))
       y(21)=dlog(dabs(m2))
       y(22)=dlog(dabs(m3))
        x1=dlog(qewsb)
c 
        x2=dlog(mz)
        h1=-h1
        if(ichoice(2).eq.11) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv1,su_rkqc)
        else if(ichoice(2).eq.21) then
      CALL SU_ODEINT(y,n,x1,x2,eps,h1,1.d-8,nok,nbad,su_deriv2,su_rkqc)
        endif
        vu = vu_mz               
        vd = vd_mz

       rmtau=y(4)*vd
       rmb = y(5)*vd
       rmtop =y(6)*vu
       mtaurun = rmtau
       mbrun = rmb
       mtrun = rmtop
       endif

       if(ichoice(7).eq.2) then
c
c====== Incorporating leading susy RC to gluino mass:  
        CALL SU_GINOCR(alsewsb,m3, mb0,mt0, delgino)
       mgluino = dabs(m3)/(1.d0 -delgino/dabs(m3))
       else
       mgluino= dabs(m3)
       endif
       mglu=mgluino
c
         if(ichoice(7).ge.1) then
c======  Incorporating mb,mt,mtau corrections:
c first redefining all needed soft etc parameters now at mz scale:
       alz=y(9)
       adz=y(10)
       auz=y(11)
       mtaurz=dsqrt(y(14))
       mslz=dsqrt(y(15))
       mbrz=dsqrt(y(16))
       mtrz=dsqrt(y(17))
       msqz=dsqrt(y(18))
       merz=dsqrt(y(24))
       melz=dsqrt(y(25))
       mdrz=dsqrt(y(26))
       murz=dsqrt(y(27))
       muqz=dsqrt(y(28))
c!! modif (temporary, until final conv) rescue in case tachyon RGE sf 
       if(irge.lt.irgmax) then   ! protections against NaN
       if(y(14).lt.0.d0) mtaurz=1.d0
       if(y(15).lt.0.d0) msLz=1.d0
       if(y(16).lt.0.d0) mbrz=1.d0
       if(y(17).lt.0.d0) mtrz= 1.d0
       if(y(18).lt.0.d0) msQz=1.d0
c
       if(y(24).lt.0.d0) merz=1.d0
       if(y(25).lt.0.d0) melz=1.d0
       if(y(26).lt.0.d0) mdrz=1.d0
       if(y(27).lt.0.d0) murz=1.d0
       if(y(28).lt.0.d0) muQz=1.d0
       else 
       if(y(14).lt.0.d0) errmess(1)=-1d0
       if(y(15).lt.0.d0)  errmess(1)=-1d0
       if(y(16).lt.0.d0)  errmess(1)=-1d0
       if(y(17).lt.0.d0)  errmess(1)=-1d0
       if(y(18).lt.0.d0)  errmess(1)=-1d0
c
       if(y(24).lt.0.d0) errmess(2)=-1d0
       if(y(25).lt.0.d0) errmess(2)=-1d0
       if(y(26).lt.0.d0) errmess(2)=-1d0
       if(y(27).lt.0.d0) errmess(2)=-1d0
       if(y(28).lt.0.d0) errmess(2)=-1d0
c
       if(errmess(1).eq.-1d0.or.errmess(2).eq.-1.d0) goto 801
       endif
c
       mu_mz=sgnmu0*dexp(y(23))
       B_mz = y(19)      
       m1z=sgnm1*dexp(y(20))
       m2z=sgnm2*dexp(y(21))
       m3z=sgnm3*dexp(y(22))
      if(irge.eq.1) then
       mtausave = rmtau
       mbsave = rmb
       mtsave= rmtop
       endif
c calculating all sfermion parameters at mz scale:
      call SU_SFBPMZ(pizz_mz,msqz,mtrz,mbrz,mslz,mtaurz,muqz,murz,mdrz,
     . melz,merz,alz,auz,adz,mu_mz,B_mz,tgbet0,rmtau,rmb,rmtop)
         if(sterr.eq.-1.d0.or.sberr.eq.-1.d0.or.stauerr.eq.-1.d0.
     . or.stnuerr.eq.-1.d0) then
c means there is really the tachyonic sfermion mass problem at Q=MZ
         errmess(4)=-1.d0
       if(errma2z.eq.-1.d0) then
c stop/ put error flag: ma^2(mz)<0 at last iter, considered irremediable
       errmess(3)=errma2z
       endif
         goto 801
         endif
       if(errma2z.eq.-1.d0) then
c stop/ put error flag: ma^2(mz)<0 at last iter, considered irremediable
       errmess(3)=errma2z
       goto 801
       endif
       CALL SU_BMSUSYCr(alphas,mb,rmtop,rmb,y(6),tgbet0,m2z
     .      ,m3z,msqz,mtrz,mbrz,auz,adz,mu_mz,  delmb) 
c Now susy RC to tau and top  masses:
       msntau_mz = dsqrt(msLz**2+0.5d0*(mz**2+pizz_mz)*dcos(2*beta_z)) 
       if(su_isNaN(msntau_mz)) msntau_mz = 1d0   ! protection 
       CALL SU_TAUMSCR(tgbet0,mu_mz,m2z,msntau_mz,  delmtau) ! changed 
c
       CALL SU_TOPMSCR(alphas,mt,mb0,rmtop,rmb,y(6),y(5),tgbet0,
     .                  m3z,msqz,mtrz,mbrz, auz,adz,mu_mz, delmtop)
c   
c  NB: SUSY RC to quark masses redefines their respective yukawas
c (we assume the top, b, tau pole masses do not change, within exp.acc.)
       if(irge.lt.irgmax) then
c  redefining running mtau,mb,mtop masses and Yuk. cplgs at Z scale:
c  modif in mb resummations (since 2.11 version): 
c  for t,b we have generically: M(pole) = M(run,Q) * (1 +CR_QCD(Q)+CR_SUSY(Q) )
c  from which we want to extract e.g. Mb(run,MZ). 
c 1) NO resummation for mtop: (mt = mt_pole,delmtop = CR_QCD(mt)+CR_SUSY(mt) 
c i.e. delmtop contains all corrections): 
       rmtop = mtpole*(1.d0 +delmtop)
c similarly for mtau:
       rmtau= mtau*(1.d0 +delmtau)
c 2) Now for mb: note that in eqs. below: rmb is mb(run,MZ)(QCD+SUSY); 
c  delmb =  CR_SUSY(MZ)only, as CR_QCD(MZ) is already taken into account before
c  Also resummation is made for mb which may be relevant for large tb
c
       rmb = mbsave/(1.d0 +delmb)
c

       y(4) = rmtau/vd
       y(5) = rmb/vd
       y(6) = rmtop/vu
       endif
c===========================
c  Now this will redefine Yukawas at high scale as well: 
       mtaurun = rmtau
       mbrun = rmb
       mtrun = rmtop
              if(irge.lt.irgmax) then
c saving some parameters:
      dmhu2=rmhu2
      dmhd2=rmhd2
      dm1=rmino1
      dm2=rmino2
      dm3=rmino3
      dtgbeta=tgbeta
      dma=ma
      dml=ml
      dmh=mh
      dmch=mch
      dmc1=gmc(1)
      dmc2=gmc(2)
      dmn1=xmn(1)
      dmn2=xmn(2)
      dmn3=xmn(3)
      dmn4=xmn(4)
      dmst1=gmst(1)
      dmst2=gmst(2)
      dmsu1=gmsu(1)
      dmsu2=gmsu(2)
      dmsb1=msb(1)
      dmsb2=msb(2)
      dmsd1=gmsd(1)
      dmsd2=gmsd(2)
      dmsl1=gmsl(1)
      dmsl2=gmsl(2)
      dmse1=gmse(1)
      dmse2=gmse(2)
      dmsn1=gmsn(1)
      dmsntau=gmsn(3)
c 
      dMSL      = msl
      dMTAUR    = mtaur
      dMSQ      = msq
      dMTR      = mtr
      dMBR      = mbr
      dMEL      = mel
      dMER      = mer
      dMUQ      = muq
      dMUR      = mur
      dMDR      = mdr
      dAL       = al
      dAU       = au
      dAD       = ad
      dAL1      = al1
      dAU1      = au1
      dAD1      = ad1
      dMA       = ma
      dMU       = mu      
c 
              goto 44 
              endif
       else
c==========
c Means that no RC are required
        mtcr=mt
        mbcr=mb
        mtaucr=mtau
         endif
c last thing: calculating now the R.C to chargino, neutralino masses:
       if(ichoice(7).eq.2) inorc=1    
       CALL SU_GAUGINO(mu,m1,m2,m3,beta,a,gmc,gmn,xmn)
      dmc1=gmc(1)
      dmc2=gmc(2)
      dmn1=xmn(1)
      dmn2=xmn(2)
      dmn3=xmn(3)
      dmn4=xmn(4)
c    ****************************************************************
c    Now comes the writing in the outputs part. 
c     ****************************************************************
 801  continue
c
c additional theoretical and experimental limits checks (g-2 etc)
       errnogo =errmess(4)+errmess(9)+errmess(10)
       if(errnogo.eq.0.d0) then
c 1) the Rho parameter (SU(2)_custodial breaking at loop-level): 
        crho=0.d0
        call su_delrho(mt,gmst,msb,gmsl,gmsn(3),thet,theb,thel,crho) 
c
c%%%%%%%%%%%%%%%%%%%%%%%%%
c  2) g_mu -2 SM + SUSY contributions:
       call su_gminus2(mel,mer,al1,mu,tgbeta,u,vv,z,dxmn,
     . gmc(1),gmc(2), gmuon)
c  3) What follow is for interface with b-> s gamma calculation:
        imod_bs=2 
        io_bs= 1 
        bsdeltp=0.9d0
        bsvkm=0.95d0    
        bsl=0.105d0
c (re)define st,sb mixing to match b->s gamma routine conventions:
c = flip angles def so that m_sf_1 > m_sf_2 (bsg conventions)
        bsthet= (thet -pi/2)/pi
        bstheb= (theb -pi/2)/pi
c        
        xsuh = min(gmst(2),mgluino,gmsu(1),gmsd(1))
        xsul = max(gmst(1),gmc(1))
        xsvl = min(gmst(1),msb(1),gmsu(1),gmsd(1),mgluino)        
        if(xsvl.ge.400.D0) then
        inlosusy =1 
        ihv = 1
        else if(dabs(bsthet).lt.0.1d0.and.xsuh.gt.2*xsul) then
        inlosusy =1  
        ihv = 0
        else
        inlosusy =0
        ihv = 0
        endif       
        bsgchm(1)=gmc(2)
        bsgchm(2)=gmc(1) 
        bsgflag=0.d0
      call chargino(tgbeta,gmc(1),mu,mmm2,bsgchm,ubsg,vbsg,ierr)
      call matching(imod_bs,io_bs,inlosusy,ihv,mw,alphas0,mt,mch,tgbeta,
     .  gmst(1),gmst(2),bsthet,msb(1),msb(2),bstheb,gmsd(1),
     .  mgluino,Au,Ad,rmu,bsgchm,
     .    ubsg,vbsg,c70,c80,c71,c81,ee,Rbox,ierr)      
      call su_bsg(alphas0,mt,mbpole-mc0,mc0/mbpole,alfinv,mw,rmb,rmb,
     . bsvkm,bsl,bsdeltp,io_bs,c70,c71,c80,c81,ee,Rbox,brsg)           
c
c 4) calculating some fine-tuning parameters for info
       call su_finetune(mu,tgbeta,rmhd2,rmhu2, 
     . czmu,czbmu,ctmu,ctbmu)
      endif
c%%%%%%%%%%%%%%%%%%%%
c saving final soft etc parameters and output masses:
c    special case:
      if(ichoice(1).eq.2) then
      rmhu2=y(12)
      rmhd2=y(13)
      endif
      if(ichoice(1).ne.2) then
      dmhu2=rmhu2
      dmhd2=rmhd2
      dm1=rmino1
      dm2=rmino2
      dm3=rmino3
      dtgbeta=tgbeta
      dMSL      = msl
      dMTAUR    = mtaur
      dMSQ      = msq
      dMTR      = mtr
      dMBR      = mbr
      dMEL      = mel
      dMER      = mer
      dMUQ      = muq
      dMUR      = mur
      dMDR      = mdr
      dAL       = al
      dAU       = au
      dAD       = ad
      dAL1      = al1
      dAU1      = au1
      dAD1      = ad1
      dMU       = mu 
      endif
c     
      dma=ma
      dml=ml
      dmh=mh
      dmch=mch
      dmc1=gmc(1)
      dmc2=gmc(2)
      dmn1=xmn(1)
      dmn2=xmn(2)
      dmn3=xmn(3)
      dmn4=xmn(4)
      dmst1=gmst(1)
      dmst2=gmst(2)
      dmsu1=gmsu(1)
      dmsu2=gmsu(2)
      dmsb1=msb(1)
      dmsb2=msb(2)
      dmsd1=gmsd(1)
      dmsd2=gmsd(2)
      dmsl1=gmsl(1)
      dmsl2=gmsl(2)
      dmse1=gmse(1)
      dmse2=gmse(2)
      dmsn1=gmsn(1)
      dmsntau=gmsn(3)
c
              if(input.ne.11) then
c writing output in the SLHA format
                open(noutlha,file='slhaspectrum.in',status='unknown')
                call su_lhaout(noutlha,ichoice,errmess,imod)
                close(noutlha)
C ************  SUSPECT OUTPUT WRITING (in SUSPECT2.out)
         if(errmess(1).eq.-1.d0 .or. 
     .      errmess(2).eq.-1.d0  .or.
     .      errmess(4).eq.-1.d0  .or.
     .      errmess(6).eq.-1.d0  .or.
     .      errmess(9).eq.-1.d0  .or.
     .      errmess(10).eq.-1.d0) then
      write(nout,'(a)')'CAUTION UNRELIABLE OUTPUT! check errmess below'
         endif
        if(ichoice(1).eq.10) then  
        write(nout,'(a)')'             SUSPECT2.4 OUTPUT: MSUGRA CASE'
        write(nout,'(a)')'             ------------------------------'
        write(nout,'(a)')
        else if(ichoice(1).eq.11) then 
        write(nout,'(a)')'             SUSPECT2.4 OUTPUT: GMSB CASE'
        write(nout,'(a)')'             ----------------------------'
        write(nout,'(a)')
        else if(ichoice(1).eq.12) then 
        write(nout,'(a)')'             SUSPECT2.4 OUTPUT: AMSB CASE'
        write(nout,'(a)')'             ----------------------------'
        write(nout,'(a)')
        else
        write(nout,'(a)')'             SUSPECT2.4 OUTPUT: pMSSM CASE'
        write(nout,'(a)')'             -----------------------------'
       endif
        if(ichoice(1).eq.0) then        
       write(nout,'(a)')'Spectrum calculation only at low (EWSB) energy
     . scale'
        write(nout,'(a)')'             -----------------------------'
       endif
        if(ichoice(1).eq.2) then        
       write(nout,'(a)')' Bottom-up: RGE from low (EWSB) to GUT energy
     . scale'
        write(nout,'(a)')'             -----------------------------'
       endif
c         
       write(nout,'(a)')'Input values:'
       write(nout,'(a)')'-------------'
        if(ichoice(1).eq.10) then
       write(nout,578)'m_0','m_1/2','A_0','tan(beta)','sign(mu)'
            write(nout,102) rm0,rmhalf,A0,tgbet0,sgnmu0
       write(nout,'(a)')
        else   if(ichoice(1).eq.11) then
       write(nout,579)'M_mess','M_susy','nl','nq','tan(beta)','sign(mu)'
            write(nout,108) mgmmess,mgmsusy, nl,nq, tgbet0,sgnmu0
       write(nout,'(a)')
        else   if(ichoice(1).eq.12) then
       write(nout,580)'M_3/2','m_0','tan(beta)','sign(mu)'
       write(nout,109)m32,am0,tgbet0,sgnmu0
       write(nout,'(a)')
       write(nout,5800)'cQ ','cuR','cdR','cL ','ceR','cHu','cHd'
       write(nout,1080)cq,cu,cd,cl,ce,chu,chd
       write(nout,'(a)')
        endif
       write(nout,581)'M_top','mb_mb','M_tau','1/alpha','sw**2(M_Z)',
     . 'alpha_S' 
             write(nout,1040) mt,mbmb,mtau,alfinv,sw20,alphas0 
       write(nout,'(a)') 
        if(ichoice(1).ne.0)then
       write(nout,582)'M_GUT','M_EWSB','E_LOW','(input or ouput scales)'
        if(ichoice(3).eq.0) then
       write(nout,105) ehigh,dsqrt(ewsb2),elow
       write(nout,'(a)')
        else if(ichoice(3).eq.1) then
       write(nout,105) egut,dsqrt(ewsb2),elow
       write(nout,'(a)')
       endif
       endif
        if(ichoice(1).eq.1) then        
       write(nout,'(a)')'Input non-universal soft terms at M_GUT'
       write(nout,'(a)')'---------------------------------------'
        endif
        if(ichoice(1).eq.0.or.ichoice(1).eq.2) then        
       write(nout,'(a)')'Input non-universal soft terms at M_EWSB'
       write(nout,'(a)')'----------------------------------------'
       endif
        if(ichoice(1).eq.0.or.ichoice(1).eq.1) then        
        if(ichoice(6).eq.0) then
       write(nout,5840)'Q_EWSB',' mu   ','M_A   ','tan(beta)','sign(mu)'
       write(nout,102) qewsb,mu0,MA,tbeta,sgnmu0
        else if(ichoice(6).eq.1) then
       write(nout,5840)'Q_EWSB','M^2_Hu','M^2_Hd','tan(beta)','sign(mu)'
       write(nout,102) qewsb,mhu20,mhd20,tbeta,sgnmu0
       write(nout,'(a)')
        endif
c      
       write(nout,585)'M_1','M_2','M_3'
       write(nout,105) m10,m20,m30
       write(nout,'(a)')
c        
       write(nout,586) 'm_eR','m_eL','m_dR','m_uR','m_qL' 
       write(nout,102) mer0,mel0,mdr0,mur0,muq0
       write(nout,'(a)')
c      
       write(nout,587)'m_tauR','m_tauL','m_bR','m_tR','m_QL'
       write(nout,102) mtaur0,msl0,mbr0,mtr0,msq0
       write(nout,'(a)')
c  
       write(nout,588)'Atau','Abottom','Atop','Al','Ad','Au'
       write(nout,104) al0,ad0,au0,al10,ad10,au10
       write(nout,'(a)')       
c  
         endif
        if(ichoice(1).eq.1.or.ichoice(1).ge.10) then
       write(nout,'(a)')
     $          'Fermion masses and gauge couplings: Q=HIGH/EWSB'
       write(nout,'(a)')'---------------------------------------------'
       write(nout,583)'M_top','M_bot','M_tau','g1','g2','g3'
       write(nout,104) mtgut,mbgut,mtaugut,
     $      sqrt(ysave(1)),sqrt(ysave(2)),sqrt(ysave(3)) 
       write(nout,104) ytewsb*vuewsb, ybewsb*vdewsb,ytauewsb*vdewsb,
     .      sqrt(5./3.)*g1ewsb,g2ewsb,sqrt(4*pi*alsewsb) 
        write(nout,'(a)')
        else
       write(nout,'(a)')'Fermion masses and gauge couplings: Q=EWSB'
       write(nout,'(a)')'------------------------------------------'
       write(nout,583)'M_top','M_bot','M_tau','g1','g2','g3'
       write(nout,104) ytewsb*vuewsb, ybewsb*vdewsb,ytauewsb*vdewsb,
     .      sqrt(5./3.)*g1ewsb,g2ewsb,sqrt(4*pi*alsewsb) 
        endif
c
         if(ichoice(1).ne.0) then
       write(nout,'(a)')'mu parameter and soft terms at M_EWSB:'
       write(nout,'(a)')'--------------------------------------'
       write(nout,5841)'mu','B','M^2_Hu','M^2_Hd'
       write(nout,1010)rmu,B,rmhu2,rmhd2
       write(nout,'(a)')
       write(nout,585)'M_1','M_2','M_3'
       write(nout,105) m1,m2,m3
       write(nout,'(a)')
       write(nout,587)'m_tauR','m_tauL','m_bR','m_tR','m_QL'
       write(nout,102) rmtaur,rml,rmbr,rmtr,rmq
       write(nout,'(a)')
       write(nout,586) 'm_eR','m_eL','m_dR','m_uR','m_qL' 
       write(nout,102) rmer,rmel,rmdr,rmur,rmuq
       write(nout,'(a)')
       write(nout,588)'Atau','Abottom','Atop','Al','Ad','Au'
       write(nout,104) al,ad,au,al1,ad1,au1
        endif
         if(ichoice(1).eq.2) then
       write(nout,'(a)')'mu parameter and soft terms at M_GUT:'
       write(nout,'(a)')'--------------------------------------'
       write(nout,5841)'mu','B','M^2_Hu','M^2_Hd'
       write(nout,1010) mugut,Bgut,mhu2gut,mhd2gut
       write(nout,'(a)')
       write(nout,585)'M_1','M_2','M_3'
       write(nout,105) m1gut,m2gut,m3gut
       write(nout,'(a)')
       write(nout,587)'m_tauR','m_tauL','m_bR','m_tR','m_QL'
       write(nout,102) mtaurgut,mslgut,mbrgut,mtrgut,msqgut
       write(nout,'(a)')
       write(nout,586) 'm_eR','m_eL','m_dR','m_uR','m_qL' 
       write(nout,102) mergut,melgut,mdrgut,murgut,muqgut
       write(nout,'(a)')
       write(nout,588)'Atau','Abottom','Atop','Al','Ad','Au'
       write(nout,104) algut,adgut,augut,al1gut,ad1gut,au1gut
        endif
       write(nout,'(a)')
       write(nout,'(a)')'Mass matrices and mixing angles:'
       write(nout,'(a)')'--------------------------------'
       write(nout,596)'tan(beta)','alpha_(h,H)'
       write(nout,103) tbeta,alfa
       write(nout,'(a)')
       write(nout,597)'thet_tau','thet_b','thet_t'
       write(nout,105) thel,theb,thet
       write(nout,'(a)')
       write(nout,598)'Z(i,j)'
       write(nout,1015) Z(1,1),Z(1,2),Z(1,3),Z(1,4)
       write(nout,1015) Z(2,1),Z(2,2),Z(2,3),Z(2,4)
       write(nout,1015) Z(3,1),Z(3,2),Z(3,3),Z(3,4)
       write(nout,1015) Z(4,1),Z(4,2),Z(4,3),Z(4,4)
       write(nout,'(a)')
       write(nout,600)'U(i,j)','V(i,j)'
       write(nout,1015) U(1,1),U(1,2),VV(1,1),VV(1,2) 
       write(nout,1015) U(2,1),U(2,2),VV(2,1),VV(2,2) 
       write(nout,'(a)')
c
       write(nout,'(a)')'Final Higgs and SUSY particle masses: '
       write(nout,'(a)')'------------------------------------- '
        if(ma2.gt.0.d0) then
       write(nout,589)'h  ','H','A','H+'
       write(nout,111) ml, mh, ma, mch
        else
       write(nout,'(a)')
       write(nout,'(a)')'MA**2 <0! NO further Higgs masses calculated'
        endif
       write(nout,'(a)')
       write(nout,590)'chi+_1','chi+_2','chi0_1','chi0_2','chi0_3',
     . 'chi0_4' 
       write(nout,104) gmc(1),gmc(2),xmn(1),xmn(2),xmn(3),xmn(4)
       write(nout,'(a)')
       write(nout,5820)'gluino'
       write(nout,106) mgluino
       write(nout,'(a)')
       write(nout,591)'stop_1','stop_2','sup_1','sup_2'
       write(nout,101) gmst(1),gmst(2),gmsu(1),gmsu(2)
       write(nout,'(a)')
       write(nout,592)'sbot_1','sbot_2','sdown_1','sdown_2'
       write(nout,101) msb(1),msb(2),gmsd(1),gmsd(2)                  
       write(nout,'(a)')
       write(nout,593)'stau_1','stau_2','snutau','selec_1','selec_2',
     .'snuelec'
       write(nout,104) gmsl(1),gmsl(2),gmsn(3),gmse(1),gmse(2),gmsn(1)
       write(nout,'(a)')
        write(nout,'(a)')'Low-energy/LEP precision parameter values:'
       write(nout,597)'Delta_rho','g_mu -2','Br(b->s gamma)'
       write(nout,105) crho,gmuon,brsg
      write(nout,'(a)')'Fine-tuning values for info: fine-tuned if >>1'
      write(nout,'(a)')'dmZ^2/mZ^2(mu^2) dmZ^2/mZ^2(B.mu) dmt/mt(mu)
     . dmt/mt(B.mu)' 
      write(nout,101) czmu,czbmu,ctmu,ctbmu
       write(nout,'(a)')
 1000  if(iknowl.ne.0) then
        write(nout,'(a)')'Warning/Error Flags: errmess(1)-(10):'
        write(nout,'(a)')'-------------------------------------'
        write(nout,595) (errmess(ierr),ierr=1,10)
        write(nout,'(a)')'---------------------------------'
        write(nout,'(a)')'errmess(i)= 0: Everything is fine.'
        write(nout,'(a)')'errmess(1)=-1: tachyon 3rd gen. sfermion from 
     .RGE'
        write(nout,'(a)')'errmess(2)=-1: tachyon 1,2 gen. sfermion from 
     .RGE'
        write(nout,'(a)')'errmess(3)=-1: tachyon A    (maybe temporary: 
     .see final mass) '
        write(nout,'(a)')'errmess(4)=-1: tachyon 3rd gen. sfermion from
     .mixing'
        write(nout,'(a)')'errmess(5)=-1: mu unstable after many iter'
        write(nout,'(a)')'errmess(6)=-1: non-convergent mu from EWSB '
        write(nout,'(a)')'errmess(7)=-1: EWSB maybe inconsistent 
     .(!but RG-improved only check)' 
        write(nout,'(a)')'errmess(8)=-1: V_eff maybe UFB or CCB 
     .(!but RG-improved only check)' 
        write(nout,'(a)')'errmess(9)=-1: Higgs boson masses are NaN '
        write(nout,'(a)')'errmess(10)=-1: RGE problems (non-pert and/or
     .Landau poles)'
      if(errmess(1).eq.-1.d0) then
      write(nout,'(a)') 'Bad input: one m^2(3rd gen. sf) <0 from RGE '
      write(nout,'(a)') 'maybe artefact of algorithm, see final result'
      endif
      if(errmess(2).eq.-1.d0) then
      write(nout,'(a)') 'Bad input: one m^2(1,2 gen. sf) <0 from RGE '
      write(nout,'(a)') 'maybe artefact of algorithm, see final result'
      endif
       if(errmess(1).eq.-1.d0.or.errmess(2).eq.-1.d0) then
      write(nout,'(a)')' Tachyonic RGE: UNRELIABLE OUTPUT! '
       goto 900
        endif
      if(errmess(3).eq.-1.d0) then
      write(nout,'(a)') 'Warning:  MA^2(Q) <0 at a scale MZ<Q<EWSB ! '
      write(nout,'(a)') 'check final results '
      endif
      if(errmess(4).eq.-1.d0) then
      write(nout,'(a)') 'STOP: one tachyonic m^2(3rd gen. sf) <0 '
      write(nout,'(a)') 'UNRELIABLE OUTPUT! '
      goto 900
      endif
      if(errmess(5).eq.-1.d0) then
      write(nout,'(a)')' Warning: MU unstable after many iter'
      endif
      if(errmess(6).eq.-1.d0) then
      write(nout,'(a)')'WARNING: EWSB unconvergent after 20 iter.' 
      endif
      if(errmess(7).eq.-1.d0) then
      write(nout,'(a)')'EW Sym. Break. may be not realized '
      write(nout,'(a)')'(however from naive tree-level analysis) '
      endif
      if(errmess(8).eq.-1.d0) then
      write(nout,'(a)')'Potential may be unbounded from below '
      write(nout,'(a)')'(however from naive tree-level analysis) '
      endif
      if(errmess(9).eq.-1.d0) then
      write(nout,'(a)')' PROBLEM: some Higgs masses are NaN! '
      endif
      if(errmess(10).eq.-1.d0) then
      write(nout,'(a)')'STOP: non-pert. R.C., or Landau pole in RGE!'
      goto 900
      endif
      endif
c
c saving some parameters for output:
       mhu2=rmhu2
       mhd2=rmhd2
       m1=rmino1
       m2=rmino2
       m3=rmino3
      dtgbeta=tgbeta
c 
 578    format(2x,a3,8x,a5,6x,a3,8x,a9,2x,a8,3x)
 579    format(2x,a6,5x,a6,5x,a2,9x,a2,8x,a10,2x,a8)
 580    format(2x,a5,6x,a3,8x,a9,2x,a8,3x,a3,4x)
 581    format(2x,a5,6x,a5,6x,a5,6x,a7,4x,a10,3x,a7,5x)
 582    format(2x,a5,5x,a7,5x,a5,5x,a25,6x)
 5820   format(2x,a6,5x,a7,5x,a5,5x,a25,6x)
 583    format(2x,a5,6x,a5,6x,a5,5x,a5,6x,a5,6x,a5,6x)
 584    format(2x,a6,5x,a6,5x,a2,9x,a2,9x,a7,1x)
 5841   format(2x,a2,9x,a1,10x,a6,7x,a6,7x,a7,1x)
 5840   format(2x,a6,3x,a6,6x,a6,5x,a9,2x,a8,6x)
 585    format(2x,a3,8x,a3,8x,a3,9x)
 586    format(2x,a4,7x,a4,7x,a4,7x,a4,7x,a4,7x)
 587    format(2x,a6,5x,a6,5x,a4,7x,a4,7x,a4,7x)
 588    format(2x,a4,7x,a7,4x,a4,7x,a2,9x,a2,9x,a2)
 589    format(2x,a1,11x,a1,11x,a1,11x,a2,3x)
 590    format(2x,a6,5x,a6,5x,a6,5x,a6,5x,a6,5x,a6,5x)
 591    format(2x,a6,5x,a6,5x,a5,5x,a6,6x)
 592    format(2x,a6,5x,a6,5x,a7,4x,a7,4x)
 593    format(2x,a6,5x,a6,5x,a6,5x,a7,4x,a7,4x,a7)
 5800   format(2x,7(a3,5x))
 1080   format(2x,7(g8.3))
 596    format(a11,2x,a11)
 597    format(1x,a9,2x,a8,3x,a14)
 598    format(2x,a6,4x,a6,4x,a6)
 595    format(2x,10(f3.0))
 600    format(2x,a6,18x,a6)
 105    format(3(g11.4))
 101    format(4(g11.4))
 111    format(4(g12.5))
 1015   format(1x,4(g11.4,1x))
 1016   format(2(g11.4,1x),2x,2(g11.4,1x))
 1010   format(3(g11.4),2x,g11.4)
 102    format(5(g11.4))
 103    format(2(g11.4))
 1040   format(3(g11.4),1x,g11.5,g11.4,2x,g11.4)
 104    format(6(g11.4))
 106    format(1(g11.4))
 108    format(2(g11.4),1x,g11.4,3x,g11.4,3x,2(g11.4))
 109    format(4(g11.4))
 13     format(1x,'M(phi_u)^2, M(phi_d)^2 ',g14.6,1x,g14.6)
 777    format(3(g14.6))
c
c============== Closing the file and saving
      write(*,'(a)')' RUN TERMINATED : OUTPUT in suspect2.out'
      CLOSE(NOUT)      
                 endif
c saving some parameters for loops etc
      mhu2=rmhu2
      mhd2=rmhd2
c!!added  for rge accur:
      if(iaccsave.ge.1) ichoice(4)=iaccsave
c  final interface for all soft terms:
c!!added if NaN higgs occured in scan, to reinitialize for next point:
      if(errmess(9).eq.-1.d0) then
      ml =0.d0
      mh=0.d0
      ma=0.d0
      mch=0.d0
      endif
      dMHU2     = mhu2
      dMHD2     = mhd2
      dM1       = m1
      dM2       = m2
      dM3       = m3
      dMSL      = msl
      dMTAUR    = mtaur
      dMSQ      = msq
      dMTR      = mtr
      dMBR      = mbr
      dMEL      = mel
      dMER      = mer
      dMUQ      = muq
      dMUR      = mur
      dMDR      = mdr
      dAL       = al
      dAU       = au
      dAD       = ad
      dAL1      = al1
      dAU1      = au1
      dAD1      = ad1
      dMU       = mu
      if(ichoice(1).ge.10) dMA = ma   ! if added if MA NOT input
 900   continue      
       end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++End of the main program SuSpect2+++++++++++++++
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      LOGICAL FUNCTION SU_IsNaN(a)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
c test if a is NaN. Comparison of NaN with any number gives FALSE
c=================================================================
      real*8 a
      SU_IsNaN = .not. (a.gt.0d0.or.a.lt.0d0.or.a.eq.0d0)
      end
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c               Here come the various subroutines for the models
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     
            SUBROUTINE SU_AMSBSUB(m0,m32,cq,cu,cd,cl,ce,chu,chd,
     .    g12,g22,g32,ytau,yb,yt,al,ad,au,al1,au1,ad1,mhu2,mhd2,
     .    mtaur2,msl2,mbr2,mtr2,msq2,mer2,mel2,mdr2,mur2,muq2,m1,m2,m3)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the initial conditions at initial scale (where the RGE starts)
c  in the general AMSB model (i.e. including a soft-SUSY breaking scalar mass
c  m_0 with a different weight c_i for every Higgs and sfermion scalar mass). 
c     The input parameters at the initial scale are:
c  m32: the gravitino mass,
c  m0 : the soft-SUSY breaking scalar mass term,
c  cq,cu,cd,cl,ce,chu,chd: weights  of m0 for the different soft terms, 
c  (for the original AMSB model: c_i=0 and usual minimal AMSB model: c_i=1), 
c  g12,g22,g23: gauge couplings squared, 
c  ytau,yb,yt : third generation Yukawa gauge couplings.
c     The ouputs at the initial scale are:
c  m1,m2,m3: gaugino mass terms, 
c  au,ad,al,au1,ad1,al1: 3d and 1st/2d generation trilinear couplings,
c  mhu2,mhd2,mtaur2,msl2,mbr2,mtr2,msq2,mer2,mel2,mdr2,mur2,muq2, 
c  Higgs and sfermion soft mass terms squared.   
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      real*8 nf
      nf = 6.d0
      pi = 4*datan(1.d0)
      cpi =1.d0/(16*pi**2)
      ytau2=ytau**2
      yb2=yb**2
      ytop2=yt**2
c
c  gauginos 
      betg1 = 33.d0/5*g12 
c adding 2-loop gauge+yukawa in boundary conditions:
     . +cpi*g12*((19*nf/15+9.d0/25)*g12
     . +(3*nf/5+9.d0/5)
     . *g22+44*nf/15*g32-26*ytop2/5-14*yb2/5-18*ytau2/5) 
c
      betg2 = g22
     . +cpi*g22*((nf/5+3.d0/5)*g12
     . +(7*nf-17.d0)*g22
     . +4*nf*g32 -6*ytop2 -6*yb2-2*ytau2 )   
c
      betg3 = -3*g32
     . +cpi*g32*(11*nf/30*g12+3*nf/2*g22
     . +(34*nf/3-54.d0)*g32 -4*ytop2 -4*yb2 )  
c      m1 = 33.d0/5*cpi*g12*m32 
c      m2 =         cpi*g22*m32
c      m3 = -3*     cpi*g32*m32
      m1 = cpi*betg1*m32 
      m2 = cpi*betg2*m32
      m3 = cpi*betg3*m32
c
c soft scalar masses and couplings 
      byt = yt*(-13*g12/15 -3*g22 -16*g32/3 +6*ytop2 +yb2)
     . +cpi*yt*(-22*ytop2*ytop2-5*yb2*yb2-5*yb2*ytop2-yb2*ytau2
     . +(6*g12/5+6*g22+16*g32)*ytop2+2*g12/5*yb2 
     . +(13*nf/15+403.d0/450)*g12**2 +(3*nf-21.d0/2)*g22**2
     . +(16*nf/3-304.d0/9)*g32**2 +g12*g22+136.d0/45
     . *g12*g32 +8*g22*g32 ) 
c
      byb = yb*(-7*g12/15 -3*g22 -16*g32/3 +ytop2 +6*yb2 +ytau2)
     . +cpi*yb*(-22*yb2*yb2-5*ytop2*ytop2-5*yb2*ytop2-3*yb2*ytau2 
     . -3*ytau2*ytau2 
     . +4*g12/5*ytop2+(2*g12/5+6*g22+16*g32)*yb2
     . +6*g12/5*ytau2 
     . +(7*nf/15+7.d0/18)*g12**2+(3*nf-21.d0/2)*g22**2
     . +(16*nf/3-304.d0/9)*g32**2 +g12*g22+8*g12*g32/9
     . +8*g22*g32 )  
c
      byta = ytau*(-9*g12/5 -3*g22 +3*yb2 +4*ytau2)
     . +cpi*ytau*(-10*ytau2*ytau2-9*yb2*yb2 -9*yb2*ytau2-3*yb2*ytop2
     . +(6*g22+6*g12/5)*ytau2 +(-2*g12/5+16*g32)*yb2
     . +(9*nf/5+27.d0/10)*g12**2+(3*nf-21.d0/2)*g22**2
     . +9*g12*g22/5 )   

c
c trilinear A terms (3rd generation)
      au = -byt/yt*cpi*m32
      ad = -byb/yb*cpi*m32
      al = -byta/ytau*cpi*m32 
c
c trilinear A terms (1st and 2d generations)
      dyovery4 =  ytau2 +3*yb2  -9*g12/5  -3*g22
     . +cpi*(
     . -3*ytau2*ytau2-9*yb2*yb2 -3*yb2*ytop2
     . +6*g12/5*ytau2 +(-2*g12/5+16*g32)*yb2
     . +(9*nf/5+27.d0/10)*g12**2 +(3*nf-21.d0/2)*g22**2
     . +9*g12*g22/5 )   


      dyovery5 =  3*yb2+ytau2
     . -7*g12/15 -3*g22 -16*g32/3
     . +cpi*( -9*yb2*yb2 -3*yb2*ytop2 -3*ytau2*ytau2 
     . +(-2*g12/5+16*g32)*yb2 +6*g12/5*ytau2 
     . +(7*nf/15+7.d0/18)*g12**2 +(3*nf-21.d0/2)*g22**2
     . +(16*nf/3-304.d0/9)*g32**2 +g12*g22+8*g12*g32/9
     . +8*g22*g32  )


      dyovery6 =  3*ytop2
     . -13*g12/15 -3*g22 -16*g32/3
     . +cpi*( -9*ytop2*ytop2-3*yb2*ytop2
     . +(4*g12/5 +16*g32)*ytop2
     . + (13*nf/15+403.d0/450)*g12**2+(3*nf-21.d0/2)*g22**2
     . +(16*nf/3-304.d0/9)*g32**2 +g12*g22+136.d0/45*g12*g32
     . +8*g22*g32 )
c
c      au1 = -(-13*g12/15 -3*g22 -16*g32/3)*cpi*m32
c      ad1 = -(-7*g12/15 -3*g22 -16*g32/3)*cpi*m32
c      al1 = -(-9*g12/5 -3*g22 )*cpi*m32 
c
      au1 = -dyovery6*cpi*m32
      ad1 = -dyovery5*cpi*m32
      al1 = -dyovery4*cpi*m32 
c 3rd generation fermion masses
      msl2 =(-99*g12**2/50-3*g22**2/2+ytau*byta)*cpi**2*m32**2
     .     +cl*m0**2
      mtaur2 = (-198*g12**2/25 +2*ytau*byta)*cpi**2*m32**2
     .     +ce*m0**2
      msq2 = (-11*g12**2/50-3*g22**2/2 +8*g32**2+
     .     yt*byt +yb*byb)*cpi**2*m32**2 
     .    +cq*m0**2
      mtr2 = (-88*g12**2/25 +8*g32**2 +2*yt*byt)*cpi**2*m32**2 
     .     +cu*m0**2
      mbr2 = (-22*g12**2/25 +8*g32**2 +2*yb*byb)*cpi**2*m32**2 
     .     +cd*m0**2
c
c 1rst and 2d generations (degenerate with 3rd):
      mel2 = (-99*g12**2/50-3*g22**2/2)*cpi**2*m32**2 
     .     +cl*m0**2
      mer2 = (-198*g12**2/25 )*cpi**2*m32**2 
     .     +ce*m0**2
      muq2 =  (-11*g12**2/50-3*g22**2/2 +8*g32**2)*cpi**2*m32**2 
     .     +cq*m0**2
      mur2 = (-88*g12**2/25 +8*g32**2 )*cpi**2*m32**2 
     .     +cu*m0**2
      mdr2  = (-22*g12**2/25 +8*g32**2)*cpi**2*m32**2 
     .      +cd*m0**2
c
c Higgs mass term^2:
      mhu2 = (-99*g12**2/50-3*g22**2/2 +3*yt*byt)*cpi**2*m32**2 
     .     +chu*m0**2
      mhd2 = (-99*g12**2/50-3*g22**2/2 +3*yb*byb +ytau*byta)*
     .     cpi**2*m32**2 
     .     +chd*m0**2
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            SUBROUTINE SU_GMSBSUB(mgmmess,mgmsusy,nl,nq, g12,g22,g32,
     . al,ad,au,al1,ad1,au1,mhu2,mhd2,mtaur2,msl2,mbr2,mtr2,msq2,mer2,
     . mel2,mdr2,mur2,muq2,m1,m2,m3)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the  GMSB model initial conditions at the messenger scale 
c  (where the RGE starts). The input at the messenger scale areMgmmess:
c  mgmmess,mgmsusy: messenger and SUSY-breaking scales,
c  nl, nq number of lepton/ quark messengers (in minimal GMSB, nl=nq=1), 
c  g12,g22,g23: gauge couplings squared.
c  The output parameters at the messenger scale are: 
c  m1,m2,m3, the gaugino masses,
c  au,ad,al,au1,ad1,al1: the trilinear sfermion couplings,  
c  mhu2,mhd2,mtaur2,msl2,mbr2,mtr2,msq2,mer2,mel2,mdr2,mur2,muq2: 
c  Higgs and sfermion soft mass terms squared.
c  The routine needs to evaluate a Spence function which is supplied: 
c            REAL*8 FUNCTION SU_PLI2(x)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit real*8(a-h,m,o-z)
      real*8 nl,nq
c
       COMMON/SU_bernou/b0,b1,b(12)
       COMMON/SU_facto/f(20)
c===== defining one-loop functions and related material:
c
      f1(x)= ( (1.d0+x)*dLog(1.d0+x) +(1.d0 -x)*dLog(1.d0-x))/x**2
c
      f2(x)= (1.d0 +x)/x**2*(dLog(1.d0+x) -2*SU_PLi2(x/(1.d0 +x)) 
     . + 1.d0/2*SU_Pli2( 2*x/(1.d0+x)) ) +
     .   (1.d0 - x)/x**2 *( dLog(1.d0-x) -2*SU_PLi2(-x/(1.d0 -x)) 
     . + 1.d0/2*SU_PLi2( -2*x/(1.d0-x)) )
c
      pi = 4*datan(1.d0)
c ---  computation of factorial for Spence fction:
c
       f(1) = dble(1)
       do 1 k = 2,19
       f(k) = f(k-1)*k
  1    continue
c
c ---  computation of the first relevant Bernouilli numbers:
c
       b0 = dble(1)
       b1 = dble(-1)/dble(2)
       b(1) = dble(1)/dble(6)
       b(2) = dble(-1)/dble(30)
       b(3) = dble(1)/dble(42)
       b(4) = dble(-1)/dble(30)
       b(5) = dble(5)/dble(66)
       b(6) = dble(-691)/dble(2730)
       b(7) = dble(7)/dble(6)
       b(8) = dble(-3617)/dble(510)
       b(9) = dble(43867)/dble(798)
c
c---- Def: mgmSUSY= F/S; mgmmess=lambdaS, x=F/(lambdaS^2)=mgmSUSY/mgmmess  
c     where lambda is the coupling of messengers to the goldstino
c     in the superpotential  gaugino masses: 
       x=  mgmSUSY/mgmmess  

c nq "quark" messengers, nl "lepton" messengers: under SU(3)*SU(2)*U(1) 
c we choose:
c   q ~ ( 3,   1,   -1/3)
c   l ~  (1,   2,    1/2) *)
c
c DC=Sum of Dynkin(messengers) times Casimir(scalar partners) convoluted by 
c the corresponding gauge couplings 
c minimal GMSB: DC for the various MSSM scalar partners *)
c
       yq=1.d0/6 
       yu=-2.d0/3
       yd=1.d0/3
       yl=-1.d0/2
       yr=1.d0 
       yh1=-1.d0/2
       yh2= 1.d0/2
c
       al1= g12/(16*pi**2)
       al2= g22/(16*pi**2)
       al3= g32/(16*pi**2)
      dcq= 4.d0/3*nq*al3**2+ 3.d0/4*nl*al2**2 +
     .     3.d0/5*yq**2/5*(2*nq +3*nl)*al1**2
      dcu= 4.d0/3*nq *al3**2 + 3.d0/5*yu**2/5*(2*nq + 3*nl)*al1**2
      dcd= 4.d0/3 *nq *al3**2 + 3.d0/5* yd**2/5*(2*nq + 3*nl)*al1**2
      dcl= 3.d0/4 *nl* al2**2 + 3.d0/5*yl**2/5 *(2 *nq + 3 *nl)*al1**2
      dce=  3.d0/5* yr**2/5 *(2* nq + 3 *nl)* al1**2 
      dchd= 3.d0/4*nl* al2**2 +3.d0/5* yh1**2/5*(2* nq + 3 *nl)*al1**2 
      dchu= 3.d0/4*nl* al2**2 +3.d0/5* yh2**2/5*(2* nq + 3* nl)*al1**2
c
c DC_i are group factor combinations for soft terms mi^2
c
c  gauginos 
      if(x.eq.1.d0) then
      f1x = 1.38629
      else
      f1x = f1(x)
      endif
      m1 = g12/(16*pi*pi)* mgmsusy/5 *(2* nq + 3* nl)* f1x 
      m2 = g22/(16*pi*pi)* mgmsusy * nl * f1x   
      m3 = g32/(16*pi*pi) *mgmsusy* nq *f1x
c
c  soft scalar masses and couplings 
      au  =0.d0
      ad  =0.D0 
      al  =0.d0 
      au1 =0.d0
      ad1 =0.D0 
      al1 =0.d0       
c 3rd generation
      if(x.eq.1.d0) then
      f2x =  0.702266
      else
      f2x = f2(x)
      endif
      msl2 = 2* mgmsusy**2* dcl* f2x 
      mtaur2 = 2* mgmsusy**2* dce* f2x
      msq2 =  2* mgmsusy**2* dcq* f2x
      mtr2 = 2* mgmsusy**2* dcu* f2x
      mbr2  = 2* mgmsusy**2* dcd* f2x
c 1rst and 2d generations (degenerate with 3rd):
      mel2 = 2* mgmsusy**2* dcl* f2x 
      mer2 = 2* mgmsusy**2* dce* f2x
      muq2 =  2* mgmsusy**2* dcq* f2x
      mur2 = 2* mgmsusy**2* dcu* f2x
      mdr2  = 2* mgmsusy**2* dcd* f2x
c Higgs mass term^2:
      mhd2 = 2* mgmsusy**2* dchd* f2x
      mhu2 = 2 *mgmsusy**2* dchu* f2x
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c  Here is defined the Spence function that is needed:      
c  ===================================================
        REAL*8 FUNCTION SU_PLI2(X)
        implicit real*8(a-h,o-z)
        COMMON/SU_bernou/b0,b1,b(12)
        COMMON/SU_facto/f(20)
        pi = 4*datan(1.d0)
        rz = 0.d0
        if(x.le.0.5d0) then
        z = -dlog(1.d0 - x)
        do 1 i = 1,9
        ii = 2*i+1
        rz = rz + b(i)/f(ii)*z**ii
  1     continue
        sp = b0*z + b1*z*z/2.d0 + rz
        else
        z = -dlog(x)
        do 2 i = 1,9
        ii = 2*i+1
        rz = rz + b(i)/f(ii)*z**ii
  2     continue
        sp = -(b0*z + b1*z*z/2.d0 + rz) + pi*pi/dble(6)
     .       -dlog(x)*dlog(1.d0-x)
        endif
        SU_PLI2 = sp
        end
c  
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   ++++++++++++++++++++++ End of the routines for models  ++++++++++++++
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c         Here come the subroutines for the radiative corrections to 
c electroweak parameters (s^2_W etc) in Drbar scheme:
c This routine is for the running gauge couplings and sin^2theta_W a la BPMZ. 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE su_RUNNINGCP(alphas,mt,rmt,rmg,tbeta,pizz,piww,piww0,
     . alphadr,alphasdr,sw2dr)
      implicit real*8 (a-h,m,o-z)
      dimension u(2,2),v(2,2),z(4,4),dxmn(4)
      COMMON/su_PARAM/GF,ALPH,MZ,MW
      COMMON/SU_bpmz/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     .      msta1,msta2,msb1,msb2,mst1,mst2,thet,theb,thel
      COMMON/SU_higgsrunz/ml,mh,ma,mch,alfa 
      COMMON/su_OUTGINOS/mc1,mc2,mn1,mn2,mn3,mn4,mgluino  
      COMMON/su_matino/U,V,Z,dxmn
      COMMON/SU_cpl/g12,g22,sw2  ! to have correct sw2
c Define first the functions to be used in the calculation.
      SU_frhol(r)=19.d0-33.d0/2*r+43.d0/12*r**2+7.d0/120*r**3
     .  -pi*dsqrt(r)*(4.d0-3.d0/2*r+3.d0/32*r**2+r**3/256)
     .  -pi**2*(2.d0-2*r+r**2/2)-dlog(r)*(3*r-r**2/2)      
      SU_frhoh(r)=dlog(r)**2*(3.d0/2-9.d0/r-15.d0/r**2-48.d0/r**3
     .                       -168.d0/r**4-612.d0/r**5)
     . -dlog(r)*(27.d0/2+4.d0/r-125.d0/4/r**2-558.d0/5/r**3
     .           -8307.d0/20/r**4-109321.d0/70/r**5)
     . +pi**2*(1.d0-4.d0/r-5.d0/r**2-16.d0/r**3-56.d0/r**4-204.d0/r**5)
     . +49.d0/4+2.d0/3/r+1613.d0/48/r**2+87.57d0/r**3
     . +341959.d0/1200/r**4+9737663.d0/9800/r**5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c First thing to calculate: the running of alpha, including only the
c contributions of SUSY particles and the one of the top quark. The 
c commons are as those of the HLOOP routine. 
c Note that alphamz and alphasmz are the input values read by suspect  i.e. 
c alpha(MZ) and alphas(MZ) in the MSbar scheme. 
c-------------
c
      pi=4*datan(1.d0)
      beta=datan(tbeta)
      ct= dcos(thet)
      st= dsin(thet)
      cb= dcos(theb)
      sb= dsin(theb)
      dalpha=  -alph/(2*pi)*( !!      -7*dlog(mw/mz)            
     .      +16.d0/9*dlog(mt/mz) +dlog(mch/mz)/3
     .      +4.d0/9*(2*dlog(msu1/mz)+2*dlog(msu2/mz)
     .                + dlog(mst1/mz)+dlog(mst2/mz)          ) 
     .      +1.d0/9*(2*dlog(msd1/mz)+2*dlog(msd2/mz)
     .                + dlog(msb1/mz)+dlog(msb2/mz)          )
     .      +1.d0/3*(2*dlog(mse1/mz)+2*dlog(mse2/mz)
     .                + dlog(msta1/mz)+dlog(msta2/mz)        )
     .      +4.d0/3*(dlog(mc1/mz)+dlog(mc2/mz) ) ) 

       alphadr=alph/(1.d0-dalpha)
c
      dalphas=alphas/(2*pi)*(1.d0/2              
     .     -2.d0/3*dlog(rmt/mz) 
     .     -2*dlog(abs(rmg)/mz) 
     .      -1.d0/6*(2*dlog(msu1/mz)+2*dlog(msu2/mz)
     .                + dlog(mst1/mz)+dlog(mst2/mz)  
c---- ramona removed double + in following line on 11/9/14         
     .                + 2*dlog(msd1/mz)+2*dlog(msd2/mz)
     .                + dlog(msb1/mz)+dlog(msb2/mz)          ) )
      alphasdr=alphas/(1.d0-dalphas)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Now we calculate sin^2theta_W in the DR-scheme. IN fact it turns out that
c we don't need input value measured at LEP or for the W mass in principle.
c However, to define the tree level quantities we need the value of M_W
c ----------
      cw2 = 1.d0 -sw2
c We first start with the calculation of the rho parameter at two-loop
c---- Two--loop QCD corrections to the t/b contributions. 

      drho2f=alphadr/4/pi/sw2*alphas/pi*(-2.145*mt**2/mw**2
     .  +1.262*dlog(mt/mz)-2.24d0-0.85*mz**2/mt**2)     
c---- Two loop Higgs contribution:    
      drhol=SU_frhol(ml/mt)
      if(ma/mt.le.1.9d0)then 
      drhoh=SU_frhol(mh/mt)
      drhoa=SU_frhol(mh/mt)
      else 
      drhoh=SU_frhoh(mh/mt)
      drhoa=SU_frhoh(mh/mt)
      endif      
      drho2h=(3*GF*mt**2/8/pi**2/dsqrt(2.d0))**2/3*(
     .       dcos(alfa)**2/dsin(beta)**2*drhol
     .     + dsin(alfa)**2/dsin(beta)**2*drhoh
     .     - 1.d0/tbeta**2*drhoa)   
c     
c    two--loop QCD corrections to the stop/sbottom contributions. 
      drho2s= 3*GF/(8*pi**2*dsqrt(2.d0))*2*alphas/3/pi*(1.d0+pi**2/3)*
     .     ( ct**2*(cb**2*SU_frho(mst1,msb1)+sb**2*SU_frho(mst1,msb2))+
     .       st**2*(cb**2*SU_frho(mst2,msb1)+sb**2*SU_frho(mst2,msb2))-
     .   ct**2*st**2*SU_frho(mst1,mst2)-cb**2*sb**2*SU_frho(msb1,msb2))

      drho2=drho2h+drho2f+drho2s
      drho=(pizz/mz**2-piww/mw**2)/(1.d0+pizz/mz**2) +drho2
      rhohat=1.d0/(1.d0-drho)
c------------------------------------------
c Then we calculate deltar, first evaluating the higher orders
c       
      dr2f=alphadr/4/pi/sw2/cw2*alphas/pi*(2.145d0*mt**2/mz**2
     .  +0.575d0*dlog(mt/mz)-0.224d0-0.144d0*mz**2/mt**2)
      drvb=rhohat*alphadr/4/pi/sw2*(6.d0+dlog(cw2)/sw2*(7.d0/2
     .     -5.d0/2*sw2-sw2*(5.d0-3.d0/2) ) )
c  (NB! drvb is the non-universal pure SM contribution: its contribution
c   is about 0.01 thus non-negligible in principle)
      dr1loop=rhohat*piww0/mw**2-pizz/mz**2+drvb
      dr2h= -drho2h*rhohat*(1.d0-dr1loop) 
      deltar=dr1loop +dr2h +dr2f         
c Now calculate sin^2theta_W by solving the usual equation:
       deter=alphadr*pi/dsqrt(2.d0)/mz**2/Gf/(1.d0-deltar)
       sw2dr=(1.d0-dsqrt(1.d0-4*deter))/2
       cw2dr=1.d0-sw2dr
       end
c------------------------------------------------------------
      real*8 function SU_frho(x,y)
      implicit real*8 (a-h,m,o-z)
      if(x.eq.y) then
      su_frho=0.d0
      else
      SU_frho = x+y-2*x*y/(x-y)*dlog(x/y)
      endif
      end
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SU_PIXX(sw2,g,g1,tbeta,pizz,piww,piww0,rmtop) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
c  Calculates essentially the PIzz(Q),PIww(Q) at Q=MZ and Q=0 for calculating
c  s^2_w, g1, g2 (MZ) in DRbar scheme
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B0,SU_BH,SU_BT22
      dimension u(2,2),v(2,2),z(4,4),dxmn(4),gmn(4),gmc(2)
      COMMON/su_PARAM/GF,ALPH,MZ,MW
      COMMON/SU_higgsrunz/ml,mh,ma,mch,alfa 
      COMMON/SU_outginos/mc1,mc2,mn1,mn2,mn3,mn4,mgluino
      COMMON/SU_matino/u,v,z,dxmn
      COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_fmasses/mtau,mb,mt
      COMMON/SU_renscale/scale     
      COMMON/SU_bpmz/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     .      msta1,msta2,msb1,msb2,mst1,mst2,thet,theb,thel

c     
c       
c basic parameters and definitions used:
       sq2=dsqrt(2.d0)
       pi = 4*datan(1.d0)
       cw2 = 1.d0-sw2
       sw = dsqrt(sw2)
       cw = dsqrt(cw2)
       cwm2 =1.d0/cw2
c defining other parameters for the Higgs mass calculation
       B=datan(tbeta)
       beta=B
       mup=1.d-2
       mdo=1.d-2
       me=0.5d-3
       mmu=0.106d0
       ms=0.190d0
       mcq=1.40d0
       eps=1.d-2
       eps0=eps**2
c  : 
       gmn(1)=dabs(dxmn(1))
       gmn(2)=dabs(dxmn(2))
       gmn(3)=dabs(dxmn(3))
       gmn(4)=dabs(dxmn(4))
       gmc(1)=mc1
       gmc(2)=mc2
c
       ct=dcos(thet)
       st=dsin(thet)
       cb=dcos(theb)
       sb=dsin(theb)
       cta=dcos(thel)
       sta=dsin(thel)
c
       cbeta2=1.d0/(1.d0+tbeta**2)
       cbet= dsqrt(cbeta2)
       sbet=dsqrt(1.d0-cbeta2)
       c2b =2*cbeta2-1.d0
c 
       sal=dsin(alfa)
       cal=dcos(alfa)
       s2a = 2*sal*cal 
       s2al=s2a

       mtsave = mt              
       if(rmtop.ne.0d0) mt = rmtop
c
c-----------------------------------------------------------------
c                 Z boson self-energy at q**2=mz**2 and q**2=0
c-----------------------------------------------------------------
       qsz=mz**2
c
      pizzf = 3*( (.5d0-2*sw2/3)**2+(2*sw2/3)**2)
     .*(dble(SU_BH(qsz,mt,mt))+dble(SU_BH(qsz,mcq,mcq))
     .       +dble(SU_BH(qsz,mup,mup)))
     .  + 3*((-.5d0+sw2/3)**2+(-sw2/3)**2)
     .*(dble(SU_BH(qsz,mb,mb))+dble(SU_BH(qsz,ms,ms))
     .       +dble(SU_BH(qsz,mdo,mdo)))
     .  + ((-.5d0+sw2)**2+(-sw2)**2)*(dble(SU_BH(qsz,me,me))
     .    +dble(SU_BH(qsz,mmu,mmu))+dble(SU_BH(qsz,mtau,mtau)))
     .  + .5d0**2*3*dble(SU_BH(qsz,eps,eps))
     .  -12*(.5d0-2*sw2/3)*(2*sw2/3)
     .  *(mt**2*dble(SU_B0(qsz,mt,mt))+mcq**2*dble(SU_B0(qsz,mcq,mcq))
     .  +mup**2*dble(SU_B0(qsz,mup,mup))) 
     .  -12*(-.5d0+sw2/3)*(-sw2/3)
     .  *(mb**2*dble(SU_B0(qsz,mb,mb))+ms**2*dble(SU_B0(qsz,ms,ms))
     .  +mdo**2*dble(SU_B0(qsz,mdo,mdo))) 
     .  -4*(-.5d0+sw2)*(-sw2)*(me**2*dble(SU_B0(qsz,me,me))+mmu**2
     .  *dble(SU_B0(qsz,mmu,mmu))+mtau**2*dble(SU_B0(qsz,mtau,mtau)))
c     
      pizzb = -2*cw**4*(2*qsz+mw**2-mz**2*sw**4/cw**2)
     . *dble(SU_B0(qsz,mw,mw))
     . -(8*cw**4+(cw2-sw2)**2)*dble(SU_BT22(qsz,mw,mw))
c
      pizzh0=-dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml))
      pizzhS = -dsin(beta-alfa)**2*(dble(SU_BT22(qsz,ma,mh))
     .  + dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsz,mz,mh))
     .  + dble(SU_BT22(qsz,ma,ml))-mz**2*dble(SU_B0(qsz,mz,mh)) )
     .  -(cw**2-sw**2)**2*dble(SU_BT22(qsz,mch,mch))
     .  -pizzh0
c
       pizzsu= -12*( (.5d0-2*sw2/3)*dcos(thet)**2
     .-(2*sw2/3)*dsin(thet)**2 )**2*dble(SU_BT22(qsz,mst1,mst1))
     .         -12*(-(.5d0-2*sw2/3)*dsin(thet)**2
     .+(2*sw2/3)*dcos(thet)**2 )**2*dble(SU_BT22(qsz,mst2,mst2))
     .      -24*( (.5d0)*dsin(thet)*dcos(thet) )**2
     .  *dble(SU_BT22(qsz,mst1,mst2))
     .    -24*(.5d0-2*sw2/3)**2*dble(SU_BT22(qsz,msu1,msu1))
     .    -24*(+2*sw2/3)**2*dble(SU_BT22(qsz,msu2,msu2))
c
      pizzsd= -12*( (-.5d0+sw2/3)*dcos(theb)**2
     .-(-sw2/3)*dsin(theb)**2)**2*dble(SU_BT22(qsz,msb1,msb1))
     .       -12*( -(-.5d0+sw2/3)*dsin(theb)**2
     .+(-sw2/3)*dcos(theb)**2)**2*dble(SU_BT22(qsz,msb2,msb2))
     .      -24*((-0.5d0)*dsin(theb)*dcos(theb))**2
     .  *dble(SU_BT22(qsz,msb1,msb2))
     .    -24*(-.5d0+sw2/3)**2*dble(SU_BT22(qsz,msd1,msd1))
     .    -24*(-sw2/3)**2*dble(SU_BT22(qsz,msd2,msd2))
c
      pizzsl=-4*( (-.5d0+sw2)*dcos(thel)**2
     .- (-sw2)*dsin(thel)**2 )**2*dble(SU_BT22(qsz,msta1,msta1))
     .       -4*( -(-.5d0+sw2)*dsin(thel)**2
     .  +(-sw2)*dcos(thel)**2 )**2*dble(SU_BT22(qsz,msta2,msta2))
     .      -8*((-.5d0)*dsin(thel)*dcos(thel))**2
     .  *dble(SU_BT22(qsz,msta1,msta2))
     .      -8*(-.5d0+sw2)**2*dble(SU_BT22(qsz,mse1,mse1))
     .       -8*(-sw2)**2*dble(SU_BT22(qsz,mse2,mse2))
c     .       -12*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
c correction msn1-> msntau (jlk)
     .       -8*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
     .       -4*(.5d0)**2*dble(SU_BT22(qsz,msntau,msntau))
c
      pizzs=pizzsl+pizzsd+pizzsu
c
      pizzn=0.d0
      do  i=1,4
      do  j=1,4
      pizzn = pizzn + 1.d0/4*(Z(i,3)*Z(j,3) -Z(i,4)*Z(j,4))**2*
     .       (dble(SU_BH(qsz,gmn(i),gmn(j)))
     .       -2*dxmn(i)*dxmn(j)*dble(SU_B0(qsz,gmn(i),gmn(j))) )
      enddo
      enddo
c
      pizzc=0.d0
      do i=1,2
      do j=1,2
      pizzc = pizzc +1.d0/4*( 
     .( ( 2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2) )**2+
     .  ( 2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2) )**2 )
     .            *dble(SU_BH(qsz,gmc(i),gmc(j))) 
     .     +4*(2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2))*
     .        (2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2))*
     .            gmc(i)*gmc(j)*dble(SU_B0(qsz,gmc(i),gmc(j))) )
      enddo
      enddo
c
c Sum of the susy contributions for pizz and final pizz(MZ**2) 
      pizzsm=alph/4/pi/sw2/cw2*(pizzf+pizzb+pizzh0)
      pizzsusy=alph/4/pi/sw2/cw2*
     .        (pizzhS+pizzs+pizzn+pizzc)
      pizz=pizzsm+pizzsusy
c
c----------------------------------------------------------------
      qsz=eps 
c
      pizzf0 = 3*( (.5d0-2*sw2/3)**2+(2*sw2/3)**2)
     .*(dble(SU_BH(qsz,mt,mt))+dble(SU_BH(qsz,mcq,mcq))
     .  +dble(SU_BH(qsz,mup,mup)))
     .  + 3*((-.5d0+sw2/3)**2+(-sw2/3)**2)
     .*(dble(SU_BH(qsz,mb,mb))+dble(SU_BH(qsz,ms,ms))
     .  +dble(SU_BH(qsz,mdo,mdo)))
     .  + ((-.5d0+sw2)**2+(-sw2)**2)*(dble(SU_BH(qsz,me,me))
     .    +dble(SU_BH(qsz,mmu,mmu))+dble(SU_BH(qsz,mtau,mtau)))
     .  + .5d0**2*3*dble(SU_BH(qsz,eps,eps))
     .  -12*(.5d0-2*sw2/3)*(2*sw2/3)
     .  *(mt**2*dble(SU_B0(qsz,mt,mt))+mcq**2*dble(SU_B0(qsz,mcq,mcq))
     .  +mup**2*dble(SU_B0(qsz,mup,mup))) 
     .  -12*(-.5d0+sw2/3)*(-sw2/3)
     .  *(mb**2*dble(SU_B0(qsz,mb,mb))+ms**2*dble(SU_B0(qsz,ms,ms))
     .  +mdo**2*dble(SU_B0(qsz,mdo,mdo))) 
     .  -4*(-.5d0+sw2)*(-sw2)*(me**2*dble(SU_B0(qsz,me,me))+mmu**2
     .  *dble(SU_B0(qsz,mmu,mmu))+mtau**2*dble(SU_B0(qsz,mtau,mtau)))
c
      pizzb0 = -2*cw**4*(2*qsz+mw**2-mz**2*sw**4/cw**2)
     . *dble(SU_B0(qsz,mw,mw))
     . -(8*cw**4+(cw2-sw2)**2)*dble(SU_BT22(qsz,mw,mw))
c
      pizzh00=-dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml))
      pizzhS0 = -dsin(beta-alfa)**2*(dble(SU_BT22(qsz,ma,mh))
     .  + dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsz,mz,mh))
     .  + dble(SU_BT22(qsz,ma,ml))-mz**2*dble(SU_B0(qsz,mz,mh)) )
     .  -(cw**2-sw**2)**2*dble(SU_BT22(qsz,mch,mch))
     .  -pizzh00
c
       pizzsu0= -12*( (.5d0-2*sw2/3)*dcos(thet)**2
     .-(2*sw2/3)*dsin(thet)**2 )**2*dble(SU_BT22(qsz,mst1,mst1))
     .         -12*(-(.5d0-2*sw2/3)*dsin(thet)**2
     .+(2*sw2/3)*dcos(thet)**2 )**2*dble(SU_BT22(qsz,mst2,mst2))
     .      -24*( (.5d0)*dsin(thet)*dcos(thet) )**2
     .  *dble(SU_BT22(qsz,mst1,mst2))
     .    -24*(.5d0-2*sw2/3)**2*dble(SU_BT22(qsz,msu1,msu1))
     .    -24*(+2*sw2/3)**2*dble(SU_BT22(qsz,msu2,msu2))
c
      pizzsd0= -12*( (-.5d0+sw2/3)*dcos(theb)**2
     .-(-sw2/3)*dsin(theb)**2)**2*dble(SU_BT22(qsz,msb1,msb1))
     .       -12*( -(-.5d0+sw2/3)*dsin(theb)**2
     .+(-sw2/3)*dcos(theb)**2)**2*dble(SU_BT22(qsz,msb2,msb2))
     .      -24*((-0.5d0)*dsin(theb)*dcos(theb))**2
     .  *dble(SU_BT22(qsz,msb1,msb2))
     .    -24*(-.5d0+sw2/3)**2*dble(SU_BT22(qsz,msd1,msd1))
     .    -24*(-sw2/3)**2*dble(SU_BT22(qsz,msd2,msd2))
c
      pizzsl0=-4*( (-.5d0+sw2)*dcos(thel)**2
     .- (-sw2)*dsin(thel)**2 )**2*dble(SU_BT22(qsz,msta1,msta1))
     .       -4*( -(-.5d0+sw2)*dsin(thel)**2
     .  +(-sw2)*dcos(thel)**2 )**2*dble(SU_BT22(qsz,msta2,msta2))
     .      -8*((-.5d0)*dsin(thel)*dcos(thel))**2
     .  *dble(SU_BT22(qsz,msta1,msta2))
     .      -8*(-.5d0+sw2)**2*dble(SU_BT22(qsz,mse1,mse1))
     .       -8*(-sw2)**2*dble(SU_BT22(qsz,mse2,mse2))
c     .       -12*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
c correction msn1-> msntau (jlk)
     .       -8*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
     .       -4*(.5d0)**2*dble(SU_BT22(qsz,msntau,msntau))
c
      pizzs0=pizzsl0+pizzsd0+pizzsu0
c
      pizzn0=0.d0
      do  i=1,4
      do  j=1,4
      pizzn0 = pizzn0 + 1.d0/4*(Z(i,3)*Z(j,3) -Z(i,4)*Z(j,4))**2*
     .       (dble(SU_BH(qsz,gmn(i),gmn(j)))
     .       -2*dxmn(i)*dxmn(j)*dble(SU_B0(qsz,gmn(i),gmn(j))) )
      enddo
      enddo
c
      pizzc0=0.d0
      do i=1,2
      do j=1,2
      pizzc0 = pizzc0 +1.d0/4*( 
     .( ( 2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2) )**2+
     .  ( 2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2) )**2 )
     .            *dble(SU_BH(qsz,gmc(i),gmc(j))) 
     .     +4*(2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2))*
     .        (2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2))*
     .            gmc(i)*gmc(j)*dble(SU_B0(qsz,gmc(i),gmc(j))) )
      enddo
      enddo
c
c Sum of the susy contributions for pizz and final pizz(MZ**2) 
      pizzsm0=alph/4/pi/sw2/cw2*(pizzf0+pizzb0+pizzh00)
      pizzsusy0=alph/4/pi/sw2/cw2*
     .        (pizzhS0+pizzs0+pizzn0+pizzc0)
      pizz0=pizzsm0+pizzsusy0
c
c-----------------------------------------------------------------
c                W boson self-energy at q**2=mw**2 and q**2=0
c-----------------------------------------------------------------
      qsw=mw**2
c
      piwwf=3.d0/2*(dble(SU_BH(qsw,mt,mb))+dble(SU_BH(qsw,mcq,ms))
     . +dble(SU_BH(qsw,mup,mdo)))+0.5d0*(dble(SU_BH(qsw,me,eps))
     . +dble(SU_BH(qsw,mmu,eps))+dble(SU_BH(qsw,mtau,eps)))
c
      piwwb=-(1.d0+8*cw**2)*dble(SU_BT22(qsw,mz,mw))-sw**2*(
     . 8*dble(SU_BT22(qsw,mw,eps))+4*qsw*dble(SU_B0(qsw,mw,eps)))
     . -((4*qsw+mz**2+mw**2)*cw**2-mz**2*sw**4)
     . *dble(SU_B0(qsw,mz,mw))
c
      piwwh0=-   dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw))         
      piwwhS = -dsin(beta-alfa)**2*(dble(SU_BT22(qsw,mh,mch))
     .  + dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsw,ml,mch))
     .  + dble(SU_BT22(qsw,mh,mw))-mw**2*dble(SU_B0(qsw,mh,mw)) )
     .  -dble(SU_BT22(qsw,ma,mch)) -piwwh0
c
      piwws =-2*3*( 2*dble(SU_BT22(qsw,msu1,msd1)) 
     .+dcos(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst1,msb1))
     .+dcos(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst1,msb2))
     .+dsin(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst2,msb1))
     .+dsin(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst2,msb2)) )
     .       -2*(  2*dble(SU_BT22(qsw,msn1,mse1)) 
c     . + dcos(thel)**2*dble(SU_BT22(qsw,msn1,msta1))
c     . + dsin(thel)**2*dble(SU_BT22(qsw,msn1,msta2)) )
c correction msn1 -> msntau (jlk)
     . + dcos(thel)**2*dble(SU_BT22(qsw,msntau,msta1))
     . + dsin(thel)**2*dble(SU_BT22(qsw,msntau,msta2)) )
c 
      piwwnc=0.d0
       do i=1,4
       do j=1,2
       piwwnc= piwwnc +
     . ( (-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)**2+
     .   (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)**2 )*
     . dble(SU_BH(qsw,gmn(i),gmc(j))) 
     . + 4*(-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)*
     .        (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)*
     . dxmn(i)*gmc(j)*dble(SU_B0(qsw,gmn(i),gmc(j)))
       enddo
       enddo   
c
c Sum of the susy contributions for piww and final piww(Mw**2)  
      piwwsm=alph/4/pi/sw2*(piwwf+piwwb+piwwh0)
      piwwsusy=alph/4/pi/sw2*(piwwhS+piwws+piwwnc)
      piww=piwwsm+piwwsusy
c
c-----------------------------------------------------------------
      qsw=eps
c
      piwwf0=3.d0/2*(dble(SU_BH(qsw,mt,mb))+dble(SU_BH(qsw,mcq,ms))
     . +dble(SU_BH(qsw,mup,mdo)))+0.5d0*(dble(SU_BH(qsw,me,eps))
     . +dble(SU_BH(qsw,mmu,eps))+dble(SU_BH(qsw,mtau,eps)))
c
      piwwb0=-(1.d0+8*cw**2)*dble(SU_BT22(qsw,mz,mw))-sw**2*(
     . 8*dble(SU_BT22(qsw,mw,eps))+4*qsw*dble(SU_B0(qsw,mw,eps)))
     . -((4*qsw+mz**2+mw**2)*cw**2-mz**2*sw**4)
     . *dble(SU_B0(qsw,mz,mw))
c
      piwwh00 = -dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw))         
      piwwhS0 = -dsin(beta-alfa)**2*(dble(SU_BT22(qsw,mh,mch))
     .  + dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsw,ml,mch))
     .  + dble(SU_BT22(qsw,mh,mw))-mw**2*dble(SU_B0(qsw,mh,mw)) )
     .  -dble(SU_BT22(qsw,ma,mch)) -piwwh00
c
      piwws0 =-2*3*( 2*dble(SU_BT22(qsw,msu1,msd1)) 
     .+dcos(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst1,msb1))
     .+dcos(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst1,msb2))
     .+dsin(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst2,msb1))
     .+dsin(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst2,msb2)) )
     .       -2*(  2*dble(SU_BT22(qsw,msn1,mse1)) 
c     . + dcos(thel)**2*dble(SU_BT22(qsw,msn1,msta1))
c     . + dsin(thel)**2*dble(SU_BT22(qsw,msn1,msta2)) )
c correction msn1 -> msntau (jlk)
     . + dcos(thel)**2*dble(SU_BT22(qsw,msntau,msta1))
     . + dsin(thel)**2*dble(SU_BT22(qsw,msntau,msta2)) )
c 
      piwwnc0=0.d0
       do i=1,4
       do j=1,2
       piwwnc0= piwwnc0 +
     . ( (-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)**2+
     .   (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)**2 )*
     . dble(SU_BH(qsw,gmn(i),gmc(j))) 
     . + 4*(-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)*
     .        (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)*
     . dxmn(i)*gmc(j)*dble(SU_B0(qsw,gmn(i),gmc(j)))
       enddo
       enddo   
c
c Sum of the susy contributions for piww and final piww(0)  
      piwwsm0=alph/4/pi/sw2*(piwwf0+piwwb0+piwwh00)
      piwwsusy0=alph/4/pi/sw2*(piwwhS0+piwws0+piwwnc0)
      piww0=piwwsm0+piwwsusy0
c
      mt = mtsave              
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
c         Here come the subroutines for the radiative corrections to 
c         the third generation fermion masses: mb,mt and mtau. 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  The following three routines are for the evaluation of the (SUSY) radiative 
c  corrections to the generation fermion masses with analytical expressions
c  from the paper of Pierce, Bagger, Matchev, Zhang, hep-ph/9606211 (PBMZ). 
c  They will need to evaluate the one--loop real (A) and two-loop complex  
c  (B0 and B1) Passarino-Veltman functions which are supplied:
c      REAL*8 FUNCTION SU_A(m)
c      COMPLEX*16 FUNCTION SU_B0(qsq,m1,m2)
c      COMPLEX*16 FUNCTION SU_B1(s,mi,mj)
c  (the arguments are the internal masses and the momentum transfer squared).
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SU_BMSUSYCR(alphas,mb,rmt,rmb,yt,tbeta,m2,mgluino,
     .             mql,mur,mdr,at,ab,mu, delmb)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the SUSY radiative corrections to the bottom mass including
c  the SUSY QCD corrections (the standard ones are calculated with RUNM)
c  and the dominant electroweak corrections due to the Yukawa couplings.
c  The input are respectively: the strong coupling constant, pole b mass,
c  the running top and bottom masses, the top Yukawa coupling, tan(beta), 
c  the SU(2) gaugino mass, the gluino mass, the 3d generation squark mass 
c  terms, the 3d generation trilinear couplings and the parameter mu. 
c  The output delmtop is the SUSY radiative correction to the bottom mass.
c  These corrections are then re-summed in the main routine.  
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B0
      COMMON/SU_param/gf,alph,mz,mw
      COMMON/SU_cpl/g12,g22,sw2    ! added to get correct sw2
      COMMON/SU_bpew/msu1e,msu2e,msd1e,msd2e,
     . mse1e,mse2e,msn1e,msntaue,
     . msta1e,msta2e,msb1e,msb2e,mst1e,mst2e,
     . thete,thebe,thele
      COMMON/SU_bpmz/msu1,msu2,msd1,msd2,
     . mse1,mse2,msn1,msntau,
     . msta1,msta2,msb1,msb2,mst1,mst2,
     . thet,theb,thel
      COMMON/SU_renscale/scale
      B11(x)= .5d0*(.5d0+1.d0/(1.d0-x)+dlog(x)/(1.d0-x)**2)
      B12(x)= .5d0*(.5d0+1.d0/(1.d0-x)+dlog(x)/(1.d0-x)**2-dlog(x))
c   Fix ren. scale (used in B0, B1) which should be M_Z in m_b
       scalsave=scale
       scale= mz
       pi=4*datan(1.d0)
       b = datan(tbeta)
c
       ct2=dcos(thet)**2
       st2=dsin(thet)**2
c
       x1= (msb1/mgluino)**2
       x2= (msb2/mgluino)**2
       mm1 = max(msb1,mgluino)
       mm2 = max(msb2,mgluino)
c
        if(x1.ge.1.d0) then
       r1= b11(x1) -dlog(mm1**2/scale**2)/2
        else
       r1=b12(x1) -dlog(mm1**2/scale**2)/2
        endif
        if(x2.ge.1.d0) then
       r2= b11(x2) -dlog(mm2**2/scale**2)/2
        else
       r2=b12(x2) -dlog(mm2**2/scale**2)/2
       endif
c
       ginosq = -alphas/pi/3*(r1+r2 -mgluino/rmb*dsin(2*theb)*
     . (dble(SU_B0(mb**2,mgluino,msb1))
     . -dble(SU_B0(mb**2,mgluino,msb2)) ))
c
      cinost = -yt**2/pi**2/16*mu* tbeta/2.d0 *dsin(2*thet)/rmt * 
     . (dble(SU_B0(mb**2,mu,mst1))-dble(SU_B0(mb**2,mu,mst2)) )   
     . -g22/16/pi**2 *mu*m2*tbeta/(mu**2-m2**2)*  
     . (ct2*dble(SU_B0(mb**2,m2,mst1))+st2*dble(SU_B0(mb**2,m2,mst2)) 
     . -ct2*dble(SU_B0(mb**2,mu,mst1)) -st2*dble(SU_B0(mb**2,mu,mst2))) 
c
       delmb= ginosq +cinost
c
       scale=scalsave
       end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_TOPMSCR(alphas,mt,mb,rmt,rmb,yt,yb,tbeta,
     .            mgl,mql,mur,mdr,at,ab,mu, delmtop)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the radiative corrections to the top quark mass including the
c  standard and SUSY QCD corrections (the standard corrections are also 
c  calculable with RUNM) and the electroweak corrections including the 
c  contributions of gauge bosons, Higgs bosons, charginos and neutralinos.
c  The input are respectively: the strong coupling constant, the pole masses,
c  running masses and Yukawa couplings of the top and bottom quarks, tan(beta),
c  the 3d generation squark mass terms and trilinear couplings and mu. 
c  The output delmtop is the radiative correction to the top quark mass. 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B0,SU_B1
      dimension u(2,2),v(2,2),z(4,4),dxmn(4),mc(2),
     . antr(4),bntr(4),antl(4),bntl(4),ant1(4),bnt1(4),ant2(4),bnt2(4),
     . actbl(2),bctbl(2),actbr(2),bctbr(2),actb1(2),bctb1(2),
     . actb2(2),bctb2(2)
      COMMON/SU_param/gf,alph,mz,mw

      COMMON/SU_higgsrunz/ml,mh,ma,mhp,alfa 

      COMMON/SU_outginos/mc1,mc2,mn1,mn2,mn3,mn4,mgluino
      COMMON/SU_matino/u,v,z,dxmn
      COMMON/SU_cpl/g12,g22,sw2
      COMMON/SU_bpmz/msu1,msu2,msd1,msd2,
     . mse1,mse2,msn1,msntau,
     . msta1,msta2,msb1,msb2,mst1,mst2,
     . thet,theb,thel
      COMMON/SU_renscale/scale     
c
c basic parameters and definitions used:
       pi=4*datan(1.d0)
       b = datan(tbeta)
c       sw2= 1.d0-(mw/mz)**2
       sw=dsqrt(sw2)
       cw=dsqrt(1.d0-sw2)
c       e=dsqrt(4*pi*alph)
c       g2=e/sw
        g2=dsqrt(g22)
        e =g2*sw
       cbeta2=1.d0/(1.d0+tbeta**2)
       cbet= dsqrt(cbeta2)
       sbet=dsqrt(1.d0-cbeta2)
       sa=dsin(alfa)
       ca=dcos(alfa)
c
c Fix ren. scale (used in B0, B1): MZ here.
      scalsave=scale
      scale= mz    
c
       ct=dcos(thet)
       st=dsin(thet)
       cb=dcos(theb)
       sb=dsin(theb)
c
c Defining couplings:  
       gtL= .5d0 -(2.d0/3)*sw**2
       gtR= (2.d0/3)*sw**2
c
c Various contributions:        
c
        dqcd = 16*pi*alphas/3*
     . (dble(SU_B1(mt**2,mgl,mst1))+dble(SU_B1(mt**2,mgl,mst2))
c$$$     .  -mgluino/mt*dsin(2*thet)*
     .  -mgl/mt*dsin(2*thet)*
     .(dble(SU_B0(mt**2,mgl,mst1))-dble(SU_B0(mt**2,mgl,mst2))))
c
        dyuk = yt**2/2*
     .(sa**2*(dble(SU_B1(mt**2,rmt,mh)) +dble(SU_B0(mt**2,rmt,mh))) +
     . ca**2*(dble(SU_B1(mt**2,rmt,ml)) +dble(SU_B0(mt**2,rmt,ml))) +
     . cbet**2*(dble(SU_B1(mt**2,rmt,ma)) -dble(SU_B0(mt**2,rmt,ma))) +
     . sbet**2*(dble(SU_B1(mt**2,rmt,mz)) -dble(SU_B0(mt**2,rmt,mz))) )
     . +.5*( (yb**2*sbet**2+yt**2*cbet**2)*dble(SU_B1(mt**2,mb,mhp)) +
     . (g2**2+yb**2*cbet**2+yt**2*sbet**2)*dble(SU_B1(mt**2,mb,mw)) ) +
     .  yb**2*cbet**2*(dble(SU_B0(mt**2,mb,mhp))
     . -dble(SU_B0(mt**2,mb,mw)) )
c
        dgauge= -e**2*(2.d0/3)**2*(5.d0+3*dlog(mz**2/mt**2)) +
     . g2**2/cw**2*( (gtl**2+gtr**2)*dble(SU_B1(mt**2,mt,mz)) +
     . 4*gtl*gtr*dble(SU_B0(mt**2,mt,mz)) )
c
       sq2=dsqrt(2.d0)
       ytr = -4.d0/3 
       ytl =  1.d0/3
       ybr =  2.d0/3
       ybl = 1.d0/3
c
       ap1tL = 0.d0
       bp1tL = e/cw/sq2*ytl
       ap1tR = e/cw/sq2*ytr
       bp1tR = 0.d0
       ap2tL = 0.d0
       bp2tL = sq2*g2*(.5d0)
       ap2tR = 0.d0
       bp2tR = 0.d0
       ap3tL = 0.d0
       ap3tR = 0.d0
       bp3tL = 0.d0
       bp3tR = 0.d0
       ap4tL = yt
       ap4tR = 0.d0
       bp4tL = 0.d0
       bp4tR = yt
c
       do i=1,4
       aNtR(i) = Z(i,1)*ap1tR +Z(i,2)*ap2tR +Z(i,3)*ap3tR +Z(i,4)*ap4tR	  
       bNtR(i) = Z(i,1)*bp1tR +Z(i,2)*bp2tR +Z(i,3)*bp3tR +Z(i,4)*bp4tR	  
       aNtL(i) = Z(i,1)*ap1tL +Z(i,2)*ap2tL +Z(i,3)*ap3tL +Z(i,4)*ap4tL	  
       bNtL(i) = Z(i,1)*bp1tL +Z(i,2)*bp2tL +Z(i,3)*bp3tL +Z(i,4)*bp4tL	  
       enddo
c
       do i=1,4
       aNt1(i) = ct*aNtL(i) +st*aNtR(i)
       bNt1(i) = ct*bNtL(i) +st*bNtR(i)
       aNt2(i) = -st*aNtL(i) +ct*aNtR(i)
       bNt2(i) = -st*bNtL(i) +ct*bNtR(i)
       enddo
c
       aX1tbL = 0.d0
       bX1tbL = g2
       aX1tbR = 0.d0
       bX1tbR = 0.d0
       aX2tbL = -yt
       bX2tbL = 0.d0
       aX2tbR = 0.d0
       bX2tbR = -yb
c
       do i=1,2
       aCtbL(i) = V(i,1)*aX1tbL +V(i,2)*aX2tbL
       bCtbL(i) = U(i,1)*bX1tbL +U(i,2)*bX2tbL
       aCtbR(i) = V(i,1)*aX1tbR +V(i,2)*aX2tbR
       bCtbR(i) = U(i,1)*bX1tbR +U(i,2)*bX2tbR
       enddo
c
       do i=1,2
       aCtb1(i) = cb*aCtbL(i) +sb*aCtbR(i)
       bCtb1(i) = cb*bCtbL(i) +sb*bCtbR(i)
       aCtb2(i) = -sb*aCtbL(i) +cb*aCtbR(i)
       bCtb2(i) = -sb*bCtbL(i) +cb*bCtbR(i)
       enddo
c
       dNino = 0.d0
       do i=1,4
       dNino = dNino +
     . .5*( (aNt1(i)**2+bNt1(i)**2)*dble(SU_B1(mt**2,dxmn(i),mst1)) +
     . 2*aNt1(i)*bNt1(i)*dxmn(i)/mt*dble(SU_B0(mt**2,dxmn(i),mst1)) ) + 
     . .5*( (aNt2(i)**2+bNt2(i)**2)*dble(SU_B1(mt**2,dxmn(i),mst2)) +
     . 2*aNt2(i)*bNt2(i)*dxmn(i)/mt*dble(SU_B0(mt**2,dxmn(i),mst2)) ) 
       enddo
c
       mc(1)=mc1
       mc(2)=mc2
       dCino =0.d0
       do i=1,2
       dCino = dCino +
     . .5*( (aCtb1(i)**2+bCtb1(i)**2)*dble(SU_B1(mt**2,mc(i),msb1)) +
     . 2*aCtb1(i)*bCtb1(i)*mc(i)/mt*dble(SU_B0(mt**2,mc(i),msb1)) ) + 
     . .5*( (aCtb2(i)**2+bCtb2(i)**2)*dble(SU_B1(mt**2,mc(i),msb2)) +
     . 2*aCtb2(i)*bCtb2(i)*mc(i)/mt*dble(SU_B0(mt**2,mc(i),msb2)) )   
       enddo
c pure QCD correction (including logs):
        mtlog = dlog((rmt/mz)**2)
        delmt = alphas/pi*(5.d0/3 -mtlog)
     .  +alphas**2*(0.538d0 -0.1815*mtlog +0.038*mtlog**2)

c SUSY Contributions added:
      delmtop= -delmt +(dqcd +dyuk +dgauge +dNino +dCino)/(16*pi**2)
       scale=scalsave
       end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SU_TAUMSCR(tgbeta,mu,m2,mnstau, delmtau)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the dominant SUSY radiative corrections to the tau  mass 
c  with the contribution of charginos/stau sneutrinos without re-summation. 
c  The input are respectively: tan(beta), the higgsino mass parameter mu, 
c  the SU(2) gaugino mass parameter and the 3d generation sneutrino mass. 
c  The output delmtau is the radiative correction to the tau lepton mass.
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B0
      COMMON/SU_param/gf,alph,mz,mw
      COMMON/SU_cpl/g12,g22,sw2   ! added to get correct sw2
      COMMON/SU_renscale/scale
c   fix ren. scale (used in B0, B1): should be M_Z in this case
      scalsave=scale
      scale = mz   
      mtau=1.7771d0
      pi=4*datan(1.d0)
c
       cinostau=  
     . g22/16/pi**2 *mu*m2*tgbeta/(mu**2-m2**2)*  
     . (dble(SU_B0(mtau**2,m2,mnstau))-dble(SU_B0(mtau**2,mu,mnstau)) ) 
       delmtau= cinostau
       scale=scalsave
       end
c 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c      The Passarino-Veltman one (A) and two points (B0,B1) functions   
c   needed for the evalution of the radiative corrections (and also V_loop).    
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real*8 FUNCTION SU_A(m)
      implicit real*8 (a-h,m,o-z)
      COMMON/SU_renscale/scale
      if(m.ne.0.d0) then
      SU_A = m*m * (1.d0-dlog(m*m/scale/scale)) 
      else
      SU_A = 0.d0
      endif
      end
c
c*****************************************************************      
* the main scalar two-point function B0
c from Looptools http://www.feynarts.de/looptools/
c
	double complex function su_B0(p, mm1, mm2)
	implicit none
	double precision p, m1, m2, mm1, mm2
	double precision mudim2, divergence, lambda2, scale
	double precision acc, eps
	double complex Ieps, onePeps, oneMeps
	common /su_cutoff/ mudim2, divergence, lambda2
        COMMON/SU_renscale/scale
	parameter (acc = 1D-12)
	parameter (eps = 1D-20)
	parameter (Ieps = (0,1)*eps)
	parameter (onePeps = 1 + Ieps)
	parameter (oneMeps = 1 - Ieps)

	double complex fpv, xlogx
	external fpv, xlogx

	double complex x1, x2, y1, y2, r
	double precision minacc
        m1 = mm1**2
        m2 = mm2**2
        divergence=0.d0
        lambda2=0.d0
        mudim2 = scale**2
	minacc = acc*(m1 + m2)
* general case
	if(abs(p) .gt. minacc) then
	  call roots(p, m1, m2, x1, x2, y1, y2, r)
	  if(abs(y1) .gt. .5D0 .and. abs(y2) .gt. .5D0) then
	    su_B0 = -log(m2/mudim2) - 
     +        fpv(1, x1, y1) - fpv(1, x2, y2)
	  else if(abs(x1) .lt. 10 .and. abs(x2) .lt. 10) then
	    su_B0 = 2 - log(p*oneMeps/mudim2) +
     +        xlogx(-x1) + xlogx(-x2) - xlogx(y1) - xlogx(y2)
	  else if(abs(x1) .gt. .5D0 .and. abs(x2) .gt. .5D0) then
	    su_B0 = -log(m1/mudim2) -
     +        fpv(1, y1, x1) - fpv(1, y2, x2)
	  else
c	    print *, "B0(", p, ",", m1, ",", m2, ") not defined"
	    su_B0 = 999D300
	  endif
* zero momentum
	else if(abs(m1 - m2) .gt. minacc) then
	  x2 = oneMeps*m1/(m1 - m2)
	  y2 = oneMeps*m2/(m2 - m1)
	  if(abs(y2) .gt. .5D0) then
	    su_B0 = -log(m2/mudim2) - fpv(1, x2, y2)
	  else
	    su_B0 = -log(m1/mudim2) - fpv(1, y2, x2)
	  endif
	else
	  su_B0 = -log(m2/mudim2)
	endif
        su_B0 = su_B0 + divergence
	end
c------------------------
c auxiliary functions used by the B0,B1 two-point functions
c from Looptools http://www.feynarts.de/looptools/
	subroutine roots(p, m1, m2, x1, x2, y1, y2, r)
	implicit none
	double precision p, m1, m2
	double complex x1, x2, y1, y2, r
	double precision mudim2, divergence, lambda2
	common /su_cutoff/ mudim2, divergence, lambda2

	double precision acc, eps
	double complex Ieps, onePeps, oneMeps
	parameter (acc = 1D-12)
	parameter (eps = 1D-20)
	parameter (Ieps = (0,1)*eps)
	parameter (onePeps = 1 + Ieps)
	parameter (oneMeps = 1 - Ieps)

	double precision q

	r = sqrt(dcmplx(p*(p - 2*(m1 + m2)) + (m1 - m2)**2))
	q = p + m1 - m2
	x1 = (q + r)/2D0/p
	x2 = (q - r)/2D0/p
	if(abs(x2) .gt. abs(x1)) then
	  x1 = m1/p/x2
	else if(abs(x1) .gt. abs(x2)) then
	  x2 = m1/p/x1
	endif
	x1 = x1 + abs(p*x1)/p*Ieps
	x2 = x2 - abs(p*x2)/p*Ieps

	q = p - m1 + m2
	y2 = (q + r)/2D0/p
	y1 = (q - r)/2D0/p
	if(abs(y2) .gt. abs(y1)) then
	  y1 = m2/p/y2
	else if(abs(y1) .gt. abs(y2)) then
	  y2 = m2/p/y1
	endif
	y1 = y1 - abs(p*y1)/p*Ieps
	y2 = y2 + abs(p*y2)/p*Ieps
	end
c************************************************************************
	double complex function fpv(n, x, y)
c from LoopTools http://www.feynarts.de/looptools/
	implicit none
	integer n
	double complex x, y
	double precision mudim2, divergence, lambda2
	common /su_cutoff/ mudim2, divergence, lambda2
	double precision acc, eps
	double complex Ieps, onePeps, oneMeps
	parameter (acc = 1D-12)
	parameter (eps = 1D-20)
	parameter (Ieps = (0,1)*eps)
	parameter (onePeps = 1 + Ieps)
	parameter (oneMeps = 1 - Ieps)

	integer m
	double complex xm

	if(abs(x) .lt. 10) then
	  if(n .eq. 0) then
	    fpv = -log(-y/x)
	  else if(abs(x) .lt. acc) then
	    fpv = -1D0/n
	  else
	    fpv = 0
	    xm = 1
	    do m = 0, n - 1
	      fpv = fpv - xm/(n - m)
	      xm = xm*x
	    enddo
	    fpv = fpv - xm*log(-y/x)
	  endif
	else
	  fpv = 0
	  xm = 1
	  do m = 1, 30
	    xm = xm/x
	    fpv = fpv + xm/(m + n)
	    if(abs(xm/fpv) .lt. acc**2) return
	  enddo
	endif
	end
c************************************************************************
	double complex function yfpv(n, x, y)
c from Looptools http://www.feynarts.de/looptools/
	implicit none
	integer n
	double complex x, y

	double complex fpv
	external fpv

	if(abs(y) .eq. 0) then
	  yfpv = 0
	else
	  yfpv = y*fpv(n, x, y)
	endif
	end
c************************************************************************
	double complex function xlogx(x)
	implicit none
	double complex x

	if(abs(x) .eq. 0) then
	  xlogx = 0
	else
	  xlogx = x*log(x)
	endif
	end
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
      COMPLEX*16 FUNCTION SU_B1(s,mi,mj)
      implicit real*8 (a-h,m-z)
      complex*16 SU_B0
      COMMON/SU_renscale/scale
c
      if(mi.eq.mj) then
      su_b1 = su_b0(s,mi,mj)/2
      else
c      if(qsq.eq.0d0) then
c        su_B1 = (1d0-dLog(mj**2/scale**2)+
c     . mi**4/(mi**2-mj**2)**2*dLog(mj**2/mi**2)
c     .        +(mi**2+mj**2)/(mi**2-mj**2)/2)/2
c      else
      SU_B1= ( SU_A(mj)/s-SU_A(mi)/s+(1.d0+mi**2/s-mj**2/s)*
     . SU_B0(s,mi,mj) )/2  
c      endif
      endif
      end
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      COMPLEX*16 FUNCTION SU_BT22(S,MI,MJ)
      IMPLICIT real*8 (A-H,M-Z)
      COMPLEX*16 su_B0
      common/SU_renscale/scale
c      su_BT22 = 1.D0/6*( su_A(MI)/2 +su_A(MJ)/2 +
c     .   (MI**2+MJ**2-S/2)*su_B0(S,MI,MJ) 
c     .  +(MJ**2-MI**2)/S/2*( su_A(MJ)-su_A(MI)-(MJ**2-MI**2)*
c     .  su_B0(S,MI,MJ) ) +MI**2+MJ**2-S/3 ) 
c     .  -su_A(MI)/4 -su_A(MJ)/4
      su_BT22 = S/6*( su_A(MI)/S/2 +su_A(MJ)/S/2 +
     .   (MI**2/S+MJ**2/S-1.d0/2)*su_B0(S,MI,MJ) 
     .  +(MJ**2/S-MI**2/S)/2*( su_A(MJ)/S-su_A(MI)/S-(MJ**2/S-MI**2/S)*
     .  su_B0(S,MI,MJ) ) +MI**2/S+MJ**2/S-1.d0/3  
     .  -3*su_A(MI)/S/2 -3*su_A(MJ)/S/2  )
      END
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      COMPLEX*16 FUNCTION su_BH(S,MI,MJ)
      IMPLICIT real*8 (A-H,M-Z)
      COMPLEX*16 su_B0
      common/su_renscale/scale
c      su_BH = 4.d0/6*( su_A(MI)/2 +su_A(MJ)/2 +
c     .  (MI**2+MJ**2-S/2)*su_B0(S,MI,MJ) 
c     .  +(MJ**2-MI**2)/S/2*( su_A(MJ)-su_A(MI)-(MJ**2-MI**2)*
c     .  su_B0(S,MI,MJ) ) +MI**2+MJ**2-S/3 )
c     .  +(S-MI**2-MJ**2)*su_B0(S,MI,MJ)-su_A(MI)-su_A(MJ)
      su_BH = 4*S/6*( su_A(MI)/S/2 +su_A(MJ)/S/2 +
     .  (MI**2/S+MJ**2/S-1.d0/2)*su_B0(S,MI,MJ) 
     .  +(MJ**2-MI**2)/S/2*( su_A(MJ)/S-su_A(MI)/S-(MJ**2-MI**2)/S*
     .  su_B0(S,MI,MJ) ) +MI**2/S +MJ**2/S-1.d0/3 )
     . +S*((1.d0-MI**2/S-MJ**2/S)*su_B0(S,MI,MJ)-su_A(MI)/S-su_A(MJ)/S)
      END
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      COMPLEX*16 FUNCTION su_BF(S,MI,MJ)
      IMPLICIT real*8 (A-H,M-Z)
      COMPLEX*16 su_B0
      common/su_renscale/scale
      su_BF = s*(su_A(mi)/s-2*su_A(mj)/s-(2.d0+2*mi**2/s-mj**2/s)*
     .  su_B0(S,mi,mj) )
       end
c------------------------------------------------------
      COMPLEX*16 FUNCTION su_BG(S,MI,MJ)
      IMPLICIT real*8 (A-H,M-Z)
      COMPLEX*16 su_B0
      common/su_renscale/scale
      su_BG = s*( -su_A(mi)/s-su_A(mj)/s +(1.d0-mi**2/s-mj**2/s)*
     . su_B0(S,mi,mj) )
      end
c------------------------------------------------------
c   ++++++++++++++++++++++ End of the routines for models  ++++++++++++++
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  The following routines are for the evaluation of the chargino/neutralino
c  and the squark, slepton masses including the radiative corrections a la   
c  Pierce, Bagger, Matchev, Zhang, hep-ph/9606211. For these corrections, one
c  needs the one- and two-loop Passarino-Veltman functions discussed above.
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
           SUBROUTINE SU_GAUGINO(mu,m1,m2,m3,b,a,mc,mn,xmn)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
c  Calculates the chargino and neutralino masses and mixing angles (with 
c  analytical expressions) including radiative corrections in the higgsino 
c  and gaugino limits). The input parameters at EWSB scale are: 
c  mu,m1.m2,m3: Higgs mass parameter and gaugino mass parameters, 
c  b,a: datan(tan(beta)) and the mixing angle alpha in the Higgs sector.
c  The output parameters are: 
c  mc: the two chargino masses,
c  mn: the four neutralino masses (absolute values),
c  mx: the four neutralino masses (including signs).  
c  The  mass values are ordered with increasing value. The diagonalizing 
c  (ordered) mass matrices U,V for charginos and Z for neutralinos are 
c  given in the common block SU_MATINO/u,v,z
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++c
c
      implicit real*8(a-h,k-z)
      complex*16 cxa,cxb,cxc,cxd,cx1,cx2,cx3
      logical su_isNaN
      dimension mc(2),mn(4),xmn(4),z(4,4),zx(4,4),u(2,2),v(2,2),
     .          iord(4),irem(2)
      dimension x(2,2)
      dimension dxmn(4)
      dimension ymn(4),yz(4,4),xmc(2),bu(2),bv(2)
      COMMON/SU_param/gf,alph,mz,mw
      COMMON/SU_break/msl,mtaur,mq3,mu3,md3,al,au,ad,
     .            mudum,m1dum,m2dum,m3dum
      COMMON/SU_break2/ml1,merdum,mq1,murdum,mdrdum
      COMMON/SU_mssmhpar/dum1,dum2,ma,dumu
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
      COMMON/SU_matino/u,v,z,dxmn

      COMMON/run_p/pizz       
      COMMON/SU_yukaewsb/ytau,yb,ytop,alsewsb,g2ewsb,g1ewsb
      common/su_msugra/m0,mhalf,A0
      common/su_nonpert/inonpert
c

      mzsave = mz               
      mwsave = mw
      cw = 1.d0/dsqrt(1.d0+(g1ewsb/g2ewsb)**2)
      sw = g1ewsb/g2ewsb*cw
      sw2=sw**2
      if(su_isNaN(pizz).or.mz**2+pizz.le.0.d0) then
c !!! protections added
c non-pert or NaN pb, uses tree-level values temporarily:
      pizz = 0.d0
      if(irge.eq.irgmax) inonpert=-1    
      endif
      rmz = dsqrt(mz**2+pizz)   
      rmw = rmz*cw
      mz = rmz
      mw = rmw 
c
      pi=4.d0*datan(1.d0)
      sb=dsin(b)
      cb=dcos(b)
      tw=sw/cw
c  
c      if(inorc.eq.1.and.irge.ge.2) then
      if(inorc.eq.1) then   !  only at very end of calculation
c  Adding R.C to M1, M2, MU:
      m1save=m1
      m2save=m2
      musave=mu
      CALL SU_RADCINO(ml1,mq1,mq3,mu3,md3,ma,ytop,yb,m1,m2,-mu,tan(b),
     .            rcm1,rcm2,rcmu)
      m1 = m1 +rcm1             
      m2 = m2 +rcm2             ! R.C to M1,M2,MU at very last iter
      mu = mu -rcmu

      endif
c      
c     ============== Neutralino masses and matrix elements ==========
c
c   adding protection for the special(unrealistic) problem M1=M2:
      m1eqm2=0.d0
      if(dabs(M1-M2).lt.1.d-4) then
      m1eqm2=1.d0
      m1sav2=m1
      m1= m1+1.d-3
      endif
      eps=-1.d-10
      xc2=(m1*m2-mz**2-mu**2)-3.d0/8.d0*(m1+m2)**2
      xc3=-1.d0/8.d0*(m1+m2)**3+1.d0/2.d0*(m1+m2)*(m1*m2-mz**2
     .    -mu**2)+(m1+m2)*mu**2+(m1*cw**2+m2*sw**2)*mz**2
     .    -mu*mz**2*dsin(2.d0*b)
      xc4=+(m1*cw**2+m2*sw**2)*mu*mz**2*dsin(2.d0*b)-m1*m2*mu**2
     .    +1.d0/4.d0*(m1+m2)*( (m1+m2)*mu**2+(m1*cw**2+m2*sw**2)
     .    *mz**2-mu*mz**2*dsin(2.d0*b) )+1.d0/16.d0*(m1+m2)**2*
     .    (m1*m2-mz**2-mu**2)-3.d0/256.d0*(m1+m2)**4
      xs=-xc3**2-2.d0/27.d0*xc2**3+8.d0/3.d0*xc2*xc4
      xu=-1.d0/3.d0*xc2**2-4.d0*xc4
      cxd=(-4*xu**3-27*xs**2)*dcmplx(1.d0,eps)
      cxc=1.d0/2.d0*(-xs+dcmplx(0.d0,1.d0)*cdsqrt(cxd/27.d0))
      cxa=dble(cxc**(1.d0/3.d0))*dcmplx(1.d0,-eps)
      cxb=8.d0*cxa-8.d0/3.d0*xc2*dcmplx(1.d0,-eps)
C     =============== Masses and couplings:
      x0=(m1+m2)/4.d0
      cx1= cxa/2.d0-xc2/6.d0*dcmplx(1.d0,-eps)
      cx2=-cxa/2.d0-xc2/3.d0*dcmplx(1.d0,-eps)
      cx3=xc3*dcmplx(1.d0,-eps)/cdsqrt(cxb)
      xmn(1)=x0-cdabs(cdsqrt(cx1))+cdabs(cdsqrt(cx2+cx3))
      xmn(2)=x0+cdabs(cdsqrt(cx1))-cdabs(cdsqrt(cx2-cx3))
      xmn(3)=x0-cdabs(cdsqrt(cx1))-cdabs(cdsqrt(cx2+cx3))
      xmn(4)=x0+cdabs(cdsqrt(cx1))+cdabs(cdsqrt(cx2-cx3))
      do 10 i=1,4
      mn(i)=dabs(xmn(i))
      ymn(i)=xmn(i)
      zx(i,2)=-cw/sw*(m1-xmn(i))/(m2-xmn(i)) 
      zx(i,3)=(mu*(m2-xmn(i))*(m1-xmn(i))-mz**2*sb*cb*((m1-m2)*cw**2
     .       +m2-xmn(i)))/mz/(m2-xmn(i))/sw/(mu*cb+xmn(i)*sb)
      zx(i,4)=(-xmn(i)*(m2-xmn(i))*(m1-xmn(i))-mz**2*cb*cb*((m1-m2)
     .       *cw**2+m2-xmn(i)))/mz/(m2-xmn(i))/sw/(mu*cb+xmn(i)*sb)
      zx(i,1)=1.d0/dsqrt(1.d0+zx(i,2)**2+zx(i,3)**2+zx(i,4)**2) 
      yz(i,1)=zx(i,1)
      yz(i,2)=zx(i,2)*zx(i,1)
      yz(i,3)=zx(i,3)*zx(i,1)
      yz(i,4)=zx(i,4)*zx(i,1)
 10   continue
c     ============ Ordering the disorder
      if(mn(3).eq.mn(4)) mn(4)=mn(4)*(1d0+1.d-8)  !protection
c  (such a degeneracy (within d.p. accuracy) may happen for very large MU) 
      xx0 = dmin1(mn(1),mn(2),mn(3),mn(4))
      xx1 = dmax1(mn(1),mn(2),mn(3),mn(4))
      idummy = 1
      do i = 1,4
       if(mn(i).eq.xx0)then
        iord(1) = i
       elseif(mn(i).eq.xx1)then
        iord(4) = i
       else
        irem(idummy) = i
        idummy = idummy+1
       endif
      enddo
      if(mn(irem(1)).le.mn(irem(2)))then
       iord(2) = irem(1)
       iord(3) = irem(2)
      else
       iord(2) = irem(2)
       iord(3) = irem(1)
      endif
c 
      do 98 j=1,4
      i=iord(j)
      xmn(j)=ymn(i)
c  
      dxmn(j) = xmn(j)
      mn(j) =dabs(ymn(i))
        do i1=1,4
        z(j,i1)=yz(i,i1)
        enddo
 98   continue
c
c     =================== Chargino masses and matrix elements =============
c
	delta=dabs(b-.25*pi)
	ddd=mu*dcos(b)+m2*dsin(b)
	ccc=mu*dsin(b)+m2*dcos(b)
	if(delta.lt.0.01d0) then
	phim=pi/4.d0-.5d0*datan((m2-mu)/(2.d0*mw))
	phip=phim
	else if	(dabs(ccc).lt.1.d-5) then
	phim=0.d0
	phip=datan(dsqrt(2.d0)*mw*dsin(b)/(m2+1.d-5))
	else if	(dabs(ddd).lt.1.d-5) then
	phip=0.d0
	phim=datan(dsqrt(2.d0)*mw*dcos(b)/(m2+1.d-5))
	else
	rad=dsqrt((m2**2-mu**2)**2+4.d0*mw**4*dcos(2.d0*b)**2
     +	+4.d0*mw**2*(m2**2+mu**2+2.d0*m2*mu*dsin(2.d0*b)))
	phip=datan((rad-(m2**2-mu**2+2.d0*mw**2*dcos(2.d0*b)))
     +	/(2.d0*dsqrt(2.d0)*mw*(mu*dcos(b)+m2*dsin(b))))
	phim=datan((rad-(m2**2-mu**2-2.d0*mw**2*dcos(2.d0*b)))
     +	/(2.d0*dsqrt(2.d0)*mw*(mu*dsin(b)+m2*dcos(b))))
	endif
	cp=dcos(phip)
	sp=dsin(phip)
	cm=dcos(phim)
	sm=dsin(phim)
c  my convention
	u(2,2)=cm
	u(2,1)=-sm
	u(1,2)=sm
	u(1,1)=cm
	v(1,1)=cp
	v(1,2)=sp
	v(2,1)=-sp
	v(2,2)=cp
        x(1,1)=m2
        x(1,2)=dsqrt(2.d0)*mw*dsin(b)
        x(2,1)=dsqrt(2.d0)*mw*dcos(b)
        x(2,2)=mu
 555    continue
       xmc(1)=(u(1,1)*x(1,1)+u(1,2)*x(2,1))*v(1,1)
     .       +(u(1,1)*x(1,2)+u(1,2)*x(2,2))*v(1,2)
       xmc(2)=(u(2,1)*x(1,1)+u(2,2)*x(2,1))*v(2,1)
     .       +(u(2,1)*x(1,2)+u(2,2)*x(2,2))*v(2,2)
        if(xmc(1).lt.0.d0) then
c some corrections to deal with case where BOTH M1,M2 <0:
c	v(1,1)=-cp
c	v(1,2)=-sp
        v(1,1)=-v(1,1)
        v(1,2)=-v(1,2)
c	v(2,1)=-sp
c	v(2,2)=cp
        goto 555
        endif
        if(xmc(2).lt.0.d0) then
c	v(1,1)=cp
c	v(1,2)=sp
c	v(2,1)=sp
c	v(2,2)=-cp
        v(2,1)=-v(2,1)
        v(2,2)=-v(2,2)
        goto 555
        endif
        if(xmc(1).gt.xmc(2)) then
        mtemp=xmc(1)
        xmc(1)=xmc(2)
        xmc(2)=mtemp
        do j=1,2
        bu(j)=u(1,j)
        u(1,j)=u(2,j)
        u(2,j)=bu(j)
        bv(j)=v(1,j)
        v(1,j)=v(2,j)
        v(2,j)=bv(j)
        enddo
        endif        
        mc(1)=dabs(xmc(1))
        mc(2)=dabs(xmc(2))
c    Some saving
      if(m1eqm2.eq.1.d0) m1=m1sav2   ! added for m1=m2 pbs (see above)
c      if(inorc.eq.1.and.irge.ge.2) then
      if(inorc.eq.1) then
      m1 = m1save
      m2 = m2save
      mu = musave
      endif

      mz = mzsave               
      mw = mwsave

      return
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_RADCINO(ml1,mq1,mq3,mu3,md3,ma,yt,yb,m1,m2,mu,tb,
     .            rcm1,rcm2,rcmu)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the radiative corrections to the gaugino and MU mass 
c  parameters. Input parameters at EWSB scale are: 
c  ml1,mq1,mq3,mu3,md3: sfermion mass parameters of 1st and 3d generations
c  ma,tb: pseudoscalar Higgs boson mass and tan(beta)
c  yt,yb: top and bottom Yukawa couplings,
c  m1,m2,mu: bare gaugino and Higgs mass parameters.
c  The outputs are the radiative corrections  rcm1,rcm2,rcmu to m1,m2, mu
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      COMMON/SU_param/gf,alph,mz,mw 
      COMMON/SU_renscale/scale
      COMMON/SU_treesfer/msbtr1,msbtr2,msttr1,msttr2
      COMMON/SU_yukaewsb/ytauewsb,ybewsb,ytewsb,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_stepwi/wistep,h1,kpole
      complex*16 SU_B1,SU_B0
c  fix re. scale (used in B0, B1):
      if(kpole.eq.1) scale = dsqrt(msttr1*msttr2)
      pi=4*datan(1.d0)
      cw = g2ewsb/dsqrt(g1ewsb**2+g2ewsb**2)
      sw=dsqrt(1.d0-cw**2)
      b=datan(tb)
      ep=1.d-5
      amu =dabs(mu)
c
      alphewsb = (g2ewsb*sw)**2/(4*pi)
      rm1=11.d0*dble(SU_B1(m1**2,ep,mq1))+9.d0*dble(SU_B1(m1**2,ep,ml1))
     .+mu/m1*dsin(2*b)*(dble(SU_B0(m1**2,amu,ma))
     . -dble(SU_B0(m1**2,amu,mz)))
     . +dble(SU_B1(m1**2,amu,ma))+dble(SU_B1(m1**2,amu,mz))
      rcm1=-alphewsb/(4.d0*pi*cw**2)*rm1   *m1
c
      rm2=9.d0*dble(SU_B1(m2**2,ep,mq1))+3.d0*dble(SU_B1(m2**2,ep,ml1))
     .+mu/m2*dsin(2*b)*(dble(SU_B0(m2**2,amu,ma))
     . -dble(SU_B0(m2**2,amu,mz)))
     . +dble(SU_B1(m2**2,amu,ma))+dble(SU_B1(m2**2,amu,mz))
     . -4.d0*(2.d0*dble(SU_B0(m2**2,m2,mw))-dble(SU_B1(m2**2,m2,mw)))
      rcm2=-alphewsb/(4.d0*pi*sw**2)*rm2    *m2
c
      rmu1=(ybewsb**2+ytewsb**2)*dble(SU_B1(amu**2,ep,mq3))
     .    + ytewsb**2*dble(SU_B1(amu**2,ep,mu3))    
     .    + ybewsb**2*dble(SU_B1(amu**2,ep,md3))    
      rmu2=dble(SU_B1(amu**2,m2,ma))+dble(SU_B1(amu**2,m2,mz))
     .    +2.d0*dble(SU_B1(amu**2,amu,mz))
     . -4.d0*dble(SU_B0(amu**2,amu,mz))
      rcmu=(-3.d0/(32*pi**2)*rmu1-3*alphewsb/(16*pi*sw**2)*rmu2 ) *mu
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_SFERMION(mql,mur,mdr,mel,mer,al,at,ab,mu, 
     .                        mst,msb,msl,msu,msd,mse,msn)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the sfermion masses including corrections, and the mixing angles 
c  for the 3d generation sfermions. The input parameters at EWSB scale are:
c  mql,mur,mdr,mel,mer,mql1,mur1,mdr1,mel1,mer1: sfermion mass terms,
c  al,at,ab,mu: 3d generation trilinear couplings and the parameter mu
c  The outputs are the sfermions masses: mst,msb,msl,msu,msd,mse,msn.
c  The masses are ordered such that the lightest is 1 and the heaviest is 2.
c  The mixing angles of 3 generation sfermion are in the common block
c  COMMON/SU_OUTMIX/thet,theb,thel (to be treated with care because of the
c  ordering of the sfermion masses, when compared to other calculations).
c  NB this routine also calculates sfermion masses and mixing
c  in a different (BPMZ) convention used in several other subroutines
c  (the latter are passed via common/su_bpew/..) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,k-z)
      logical su_isNaN
      dimension mst(2),msb(2),msl(2),msu(2),msd(2),mse(2),msn(4)
      COMMON/SU_runmasses/amtau,amb,amt
      COMMON/SU_param/gf,alph,mz,mw
c      COMMON/SU_cpl/g12,g22,sw2        
      COMMON/SU_hmix/b,a
      COMMON/SU_break/msldum,mtaurdum,msqdum,mtrdum,mbrdum,aldum,audum,
     . addum,mudum,m1dum,m2dum,m3   
      COMMON/SU_break2/meldum,merdum,muqdum,murdum,mdrdum
      COMMON/SU_outmix/thet,theb,thel
      COMMON/SU_bpew/msu1bp,msu2bp,msd1bp,msd2bp,
     . mse1bp,mse2bp,msn1bp,msntau,
     . msta1bp,msta2bp,msb1bp,msb2bp,mst1bp,mst2bp,
     . thetbp,thebbp,thelbp
      COMMON/SU_errsf/sterr,sberr,stauerr,stnuerr
      COMMON/SU_mssmhpar/dum1,dum2,ma,dumu
      COMMON/SU_outginos/dmc1,dmc2,dmn1,dmn2,dmn3,dmn4,mgluino
      COMMON/SU_yukaewsb/ytauewsb,ybewsb,ytewsb,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_yuka/ytau,yb,ytop
      COMMON/SU_treesfer/msb1,msb2,mst1,mst2
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
      common/su_mixflip/istflip,isbflip
      COMMON/SU_renscale/scale     
      integer kpole
      COMMON/SU_stepwi/wistep,h1,kpole
      common/su_nonpert/inonpert
      COMMON/run_p/pizz         
C
      sferr = 0.d0
      mql1=muqdum
      mur1=murdum
      mdr1=mdrdum
      mel1=meldum
      mer1=merdum
      pi = 4*datan(1.d0)
      tb=dtan(b)

c Redefining s^2_w at EWSB scale:  
       cw = 1.d0/dsqrt(1.d0+(g1ewsb/g2ewsb)**2)
       sw = g1ewsb/g2ewsb *cw
       sw2ew=sw**2
       sw2 =sw2ew      
      if(su_isNaN(pizz).or.mz**2+pizz.le.0.d0) then
c !!! protections added
c non-pert or NaN pb, uses tree-level values temporarily:
      pizz = 0.d0
      if(irge.eq.irgmax) inonpert=-1    
      endif
       rmz = dsqrt(mz**2+pizz)  
       rmw = cw*rmz
       mzsave = mz
       mwsave = mw
       mz = rmz
       mw = rmw

c first two generations:  no mixing included 
c up squarks: 
      mstl2=mql1**2+(0.5d0-2.d0/3.d0*sw2)*mz**2*dcos(2.d0*b)
      mstr2=mur1**2+2.d0/3.d0*sw2*mz**2*dcos(2.d0*b) 
      msu(1)=dsqrt(mstl2)
      msu(2)=dsqrt(mstr2)
      msu1bp = msu(1)
      msu2bp =msu(2)
c
c down squarks
      msbl2=mql1**2+(-0.5d0+1.d0/3.d0*sw2)*mz**2*dcos(2.d0*b)
      msbr2=mdr1**2-1.d0/3.d0*sw2*mz**2*dcos(2.d0*b) 
      msd(1)=dsqrt(msbl2)
      msd(2)=dsqrt(msbr2)
      msd1bp = msd(1)
      msd2bp =msd(2)
c sleptons
      msel2=mel1**2+(-0.5d0+sw2)*mz**2*dcos(2.d0*b)
      mser2=mer1**2- sw2*mz**2*dcos(2.d0*b) 
      msnl2=mel1**2+0.5d0*mz**2*dcos(2.d0*b)
      mse(1)=dsqrt(msel2)
      mse(2)=dsqrt(mser2)
      mse1bp = mse(1)
      mse2bp =mse(2)
      if(msnl2.lt.0.d0) then   
      msn(1)= 1.d0
      if(irge.eq.irgmax) stnuerr=-1.d0
      else
      msn(1)=dsqrt(msnl2)
      endif
      msn1bp = msn(1)
      msn(2)=1.d+15
c
c  add radiative corrections to first gen. squarks
      if(isfrc.eq.1.and.irge.ge.2) then
      CALL SU_SQCR(alsewsb,m3,msu(1),   dmsu1)  
      CALL SU_SQCR(alsewsb,m3,msu(2),   dmsu2)  
      CALL SU_SQCR(alsewsb,m3,msd(1),   dmsd1)
      CALL SU_SQCR(alsewsb,m3,msd(2),   dmsd2)
      msu(1)=msu(1)+dmsu1
      msu(2)=msu(2)+dmsu2
      msd(1)=msd(1)+dmsd1
      msd(2)=msd(2)+dmsd2
      endif
c
c now the third generation sfermions:
c stop masses/mixing
c     first some reinitializations:
       ifirst=0
       crll=0.d0
       crrr=0.d0
       crlr=0.d0
       istflip = 0
c Mb, Mt, Ml used in sfermion matrix elements should be running masses 
c at EWSB scale, including SUSY radiative corrections

       vd = dsqrt(2*rmz**2/(g1ewsb**2+g2ewsb**2)/(1.d0+tb**2))
       vu = vd*tb
       rmb= ybewsb*vd           
       rmt= ytewsb*vu
       rml= ytauewsb*vd
c
 1    mstl2=mql**2+(0.5d0-2.d0/3*sw2ew)*mz**2*dcos(2*b)  +crll
      mstr2=mur**2+2.d0/3*sw2ew*mz**2*dcos(2*b)          +crrr
      mlrt=at-mu/tb                                          +crlr/rmt
      delt=(mstl2-mstr2)**2+4*rmt**2*mlrt**2
      mst12=rmt**2+0.5d0*(mstl2+mstr2-dsqrt(delt))
      mst22=rmt**2+0.5d0*(mstl2+mstr2+dsqrt(delt))
        if(mst12.lt.0.d0)then   
c  tachyonic sfermion 1 mass
      mst(1)= 1.d0
      if(irge.eq.irgmax) sterr=-1.d0
      else
      mst(1)=dsqrt(mst12)
      endif
      mst(2)=dsqrt(mst22)
      thet= datan(2*rmt*mlrt / (mstl2-mstr2) )/2
      if(ifirst.eq.1) mst1true=mst(1)
      if(ifirst.eq.2) mst2true=mst(2)
      if(ifirst.eq.3) then
      mst(1)=mst1true
      mst(2)=mst2true
      endif
      ct=dcos(thet)
      st=dsin(thet)

c defining stop parameters at EWSB scale in bpmz conventions
         if(ifirst.eq.0) then
      mst1bp= dsqrt(ct**2*(rmt**2+mstl2)+st**2*(rmt**2+mstr2)
     .       +2*ct*st*rmt*mlrt)
      mst2bp= dsqrt(st**2*(rmt**2+mstl2)+ct**2*(rmt**2+mstr2)
     .       -2*ct*st*rmt*mlrt)
      if(su_isNaN(mst1bp)) mst1bp=1.d0   !!added protection 
      if(su_isNaN(mst2bp)) mst2bp=1.d0  
         thetbp = thet         
         endif
        if(mstl2.gt.mstr2) then
        thet = thet + pi/2
        istflip = 1
        endif       
c 
      if(ifirst.eq.0) then
c save tree-level values for other uses:
      mst1=mst(1)
      mst2=mst(2)
      thettree=thet
      endif
c
c  adding rad. corr.
      if(isfrc.eq.1.and.irge.ge.2.and.ifirst.lt.3) then
      ifirst= ifirst +1
c calculating stop rad. corr with 3 different momenta scales:
      if(ifirst.eq.1) pscale= mst1
      if(ifirst.eq.2) pscale=mst2
      if(ifirst.eq.3) pscale=dsqrt(mst1*mst2)
c 
      CALL SU_STOPCR(pscale,MU,at,ab,m3,  crLL,crLR,crRR)  
       goto 1
       endif

       ifirst = 0              
c
c sbottom masses/mixing
c
      isbflip = 0
      msbl2=mql**2+(-0.5d0+1.d0/3*sw2ew)*mz**2*dcos(2*b)
      msbr2=mdr**2-1.d0/3*sw2ew*mz**2*dcos(2*b) 
      mlrb=ab-mu*tb
      delb=(msbl2-msbr2)**2+4*rmb**2*mlrb**2
      msb12=rmb**2+0.5d0*(msbl2+msbr2-dsqrt(delb))
      msb22=rmb**2+0.5d0*(msbl2+msbr2+dsqrt(delb))
        if(msb12.lt.0.d0)then  
c   tachyonic sfermion mass
      msb(1)=1.d0
      if(irge.eq.irgmax) sberr=-1.d0
        else
      msb(1)=dsqrt(msb12)
      endif
      msb(2)=dsqrt(msb22)
      theb= datan(2*rmb*mlrb / (msbl2-msbr2) )/2
      cb=dcos(theb)
      sb=dsin(theb)
c defining sbottom parameters at EWSB scale in bpmz conventions
         if(ifirst.eq.0) then
      msb1bp= dsqrt(cb**2*(rmb**2+msbl2)+sb**2*(rmb**2+msbr2)
     .       +2*cb*sb*rmb*mlrb)
      msb2bp= dsqrt(sb**2*(rmb**2+msbl2)+cb**2*(rmb**2+msbr2)
     .       -2*cb*sb*rmb*mlrb)
      if(su_isNaN(msb1bp)) msb1bp=1.d0   !!added protection 
      if(su_isNaN(msb2bp)) msb2bp=1.d0   
         thebbp = theb
         endif
        if(msbl2.gt.msbr2) then
        theb = theb + pi/2
        isbflip =1
        endif
c       
       if(ifirst.eq.0) then
c save tree-level values for other uses:
      msb1=msb(1)
      msb2=msb(2)
      thebtree=theb
      endif
c
c  add radiative corrections to sbottom  quarks
      if(isfrc.eq.1.and.irge.ge.2) then
      CALL SU_SQCR(alsewsb,m3,msb(1),   dmsb1)  
      CALL SU_SQCR(alsewsb,m3,msb(2),   dmsb2)
      msb(1)=msb(1)+dmsb1
      msb(2)=msb(2)+dmsb2
      endif      
c
c  stau masses/mixing
c
      msel2=mel**2+(-0.5d0+sw2ew)*mz**2*dcos(2*b)
      mser2=mer**2- sw2ew*mz**2*dcos(2*b) 
      msntau2=mel**2+0.5d0*mz**2*dcos(2*b)
      mlre=al-mu*tb
      dele=(msel2-mser2)**2+4*rml**2*mlre**2
      mse12=rml**2+0.5d0*(msel2+mser2-dsqrt(dele))
      mse22=rml**2+0.5d0*(msel2+mser2+dsqrt(dele))
        if(mse12.lt.0.d0)then   
c   tachyonic sfermion mass
      msl(1)=1.d0
      if(irge.eq.irgmax) stauerr=-1.d0
      else
      msl(1)=dsqrt(mse12)
        endif
      msl(2)=dsqrt(mse22)
      thel= datan(2*rml*mlre / (msel2-mser2) )/2
      cl=dcos(thel)
      sl=dsin(thel)
        if(msntau2.lt.0.d0) then  
       stnuerr = -1.d0
       msn(3)=1.d0
       if(irge.eq.irgmax) goto 111
        else
c   tau sneutrino:
      msn(3)=dsqrt(msntau2) 
      endif
      msntau = msn(3)
      msn(4)=1.d+15
c defining stau parameters at EWSB scale in bpmz conventions
         if(ifirst.eq.0) then
      msta1bp= dsqrt(cl**2*(rml**2+msel2)+sl**2*(rml**2+mser2)
     .       +2*cl*sl*rml*mlre)
      msta2bp= dsqrt(sl**2*(rml**2+msel2)+cl**2*(rml**2+mser2)
     .       -2*cl*sl*rml*mlre)
      if(su_isNaN(msta1bp)) msta1bp=1.d0   !!added protection
      if(su_isNaN(msta2bp)) msta2bp=1.d0  
         thelbp = thel
         endif
      if(msel2.gt.mser2) thel = thel + pi/2
c        endif  
c nb: for convenience msn(1--4) contains: 
c msn_{e,mu}(1),msn_{e,mu}(2), msn_{tau}(1),msn_{tau}(2) 

      mz=mzsave                 
      mw=mwsave
      
 111  return 
      end 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_SFBPMZ(pizz,mql,mur,mdr,mel,mer,mql1,mur1,mdr1,
     . mel1,mer1,al,at,ab,mu,B_mz,tb,rmtau,rmb,rmt)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the sfermion masses and the mixing angles at scale=MZ 
c  for the 3d generation sfermions in the BPMZ conventions.
c
      implicit real*8(a-h,k-z)
      logical su_isNaN
      dimension mst(2),msb(2)
      COMMON/SU_param/gf,alph,mz,mw
      COMMON/SU_cpl/g12,g22,sw2
      COMMON/SU_bpmz/msu1,msu2,msd1,msd2,
     . mse1,mse2,msn1,msntau,
     . msta1,msta2,msb1,msb2,mst1,mst2,
     . thet,theb,thel
      COMMON/SU_higgsrunz/mlrunz,mhrunz,marunz,mchrunz,alpharunz
      COMMON/pietro/mApole,mCHpole 
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc   
      common/su_polemz/ipolemz
      common/su_errma/errma2z
      COMMON/SU_errsf/sterr,sberr,stauerr,stnuerr
C
      errma2z=0.d0
      b=datan(tb)
      cw=dsqrt(1.d0-sw2)
      rmz = dsqrt(mz**2+pizz) 
      rmw = cw*rmz
      mzsave = mz
      mwsave = mw
      mz = rmz
      mw = rmw
      mz2 = mz**2
      if(ipolemz.eq.0) then    ! running Higgs masses in loops at mZ
      ma2 = mu*B_mZ/sin(b)/cos(b) 
      if(ma2.lt.0.d0.and.irge.lt.irgmax) ma2=.1d0   ! temp protection 
      if(ma2.lt.0.d0) errma2z=-1.d0    !final ma2<0->error flag
      ma=dsqrt(dabs(ma2))
      else if(ipolemz.eq.1) then  ! Pole Higgs masses in loops at mZ
      ma2 =mapole**2 
      ma = dsqrt(dabs(ma2))  
      endif
cc       
      marunz = ma
      mchrunz = sqrt(abs(ma2+mw**2))
      mhht2=1.d0/2*(ma2+mz2+sqrt((ma2+mz2)**2-(2*ma*mz*cos(2d0*b))**2))
      mht2=1.d0/2*(ma2+mz2-sqrt((ma2+mz2)**2-(2*ma*mz*cos(2d0*b))**2))
      mhrunz=sqrt(abs(mhht2))
      mlrunz=sqrt(abs(mht2))
      pi = 4*datan(1.d0)
      s2alt = -(ma2+mz2)*dsin(2d0*b)
      c2alt = -(ma2-mz2)*dcos(2d0*b)
      t2alt = s2alt/c2alt
      if(c2alt.gt.0) then
         alpharunz = 0.5d0*atan(t2alt)
      elseif(c2alt.lt.0) then
         alpharunz = 0.5d0*atan(t2alt) -pi/2d0
      else
         alpharunz = -pi/4d0
      endif
c
c first two generations:  no mixing included 
c up squarks: 
      mstl2=mql1**2+(0.5d0-2.d0/3.d0*sw2)*mz**2*dcos(2.d0*b)
      mstr2=mur1**2+2.d0/3.d0*sw2*mz**2*dcos(2.d0*b) 
      msu1=dsqrt(mstl2)
      msu2=dsqrt(mstr2)
c down squarks
      msbl2=mql1**2+(-0.5d0+1.d0/3.d0*sw2)*mz**2*dcos(2.d0*b)
      msbr2=mdr1**2-1.d0/3.d0*sw2*mz**2*dcos(2.d0*b) 
      msd1=dsqrt(msbl2)
      msd2=dsqrt(msbr2)
c sleptons
      msel2=mel1**2+(-0.5d0+sw2)*mz**2*dcos(2.d0*b)
      mser2=mer1**2- sw2*mz**2*dcos(2.d0*b) 
      msnl2=mel1**2+0.5d0*mz**2*dcos(2.d0*b)
      mse1=dsqrt(msel2)
      mse2=dsqrt(mser2)
      msn1=dsqrt(msnl2)
c stop parameters
c
      mstl2=mql**2+(0.5d0-2.d0/3*sw2)*mz**2*dcos(2*b)  
      mstr2=mur**2+2.d0/3*sw2*mz**2*dcos(2*b)    
      mlrt=at-mu/tb                      
      thet= datan(2*rmt*mlrt / (mstl2-mstr2) )/2
      ct=dcos(thet)
      st=dsin(thet) 
      mst1= dsqrt(ct**2*(rmt**2+mstl2)+st**2*(rmt**2+mstr2)
     .       +2*ct*st*rmt*mlrt)
      mst2= dsqrt(st**2*(rmt**2+mstl2)+ct**2*(rmt**2+mstr2)
     .       -2*ct*st*rmt*mlrt)       
      if(su_isNaN(mst1)) then    !!added protection
      mst1=92.d0    
      if(irge.eq.irgmax) sterr=-1.d0    
      endif
      if(su_isNaN(mst2)) then    
      mst2=92.d0   
      if(irge.eq.irgmax) sterr=-1.d0    
      endif
c sbottom parameters:
c
      msbl2=mql**2+(-0.5d0+1.d0/3*sw2)*mz**2*dcos(2*b)
      msbr2=mdr**2-1.d0/3*sw2*mz**2*dcos(2*b) 
      mlrb=ab-mu*tb
      theb= datan(2*rmb*mlrb / (msbl2-msbr2) )/2
      cb=dcos(theb)
      sb=dsin(theb)
      msb1= dsqrt(cb**2*(rmb**2+msbl2)+sb**2*(rmb**2+msbr2)
     .       +2*cb*sb*rmb*mlrb)
      msb2= dsqrt(sb**2*(rmb**2+msbl2)+cb**2*(rmb**2+msbr2)
     .       -2*cb*sb*rmb*mlrb)
      if(su_isNaN(msb1)) then    !!added protection 
      msb1=92.d0   
      if(irge.eq.irgmax) sberr=-1.d0    
      endif
      if(su_isNaN(msb2)) then   
      msb2=92.d0  
      if(irge.eq.irgmax) sberr=-1.d0    
      endif
c
c  stau parameters 
c
      msel2=mel**2+(-0.5d0+sw2)*mz**2*dcos(2*b)
      mser2=mer**2- sw2*mz**2*dcos(2*b) 
      msntau2=mel**2+0.5d0*mz**2*dcos(2*b)
      mlre=al-mu*tb
      thel= datan(2*rmtau*mlre / (msel2-mser2) )/2
      cl=dcos(thel)
      sl=dsin(thel)
      msta1= dsqrt(cl**2*(rmtau**2+msel2)+sl**2*(rmtau**2+mser2)
     .       +2*cl*sl*rmtau*mlre)
      msta2= dsqrt(sl**2*(rmtau**2+msel2)+cl**2*(rmtau**2+mser2)
     .       -2*cl*sl*rmtau*mlre)
      if(su_isNaN(msta1)) then    !!added protection 
      msta1=92.d0 
      if(irge.eq.irgmax) stauerr=-1.d0    
      endif
      if(su_isNaN(msta2)) then   
      msta2=92.d0  
      if(irge.eq.irgmax) stauerr=-1.d0    
      endif
c   tau sneutrino:
      msntau=dsqrt(msntau2) 
c
c
      mz=mzsave                
      mw=mwsave     
      end 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_SQCR(alphas,mgluino,msquark,dmsquark)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the QCD (standard+SUSY) correction to squark (except stop)
c  masses. The input are: the strong coupling constant alphas, the gluino 
c  and tree-level squark masses and the output is the correction to the 
c  squark mass dmsquark). Squark mixing and Yukawa's are neglected. 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      COMMON/SU_renscale/scale
      COMMON/SU_tachyrc/tachsqrc
      COMMON/SU_param/gf,alpha,mz,mw
      COMMON/SU_treesfer/msbtr1,msbtr2,msttr1,msttr2
      COMMON/SU_stepwi/wistep,h1,kpole
c   fix ren. scale (used in B0, B1):
       if(kpole.eq.1) scale = dsqrt(msttr1*msttr2)    
c
       pi=4*datan(1.d0)
       x=(mgluino/msquark)**2
       corr2=2*alphas/pi/3*(1.d0+3*x+(x-1.d0)**2*dlog(dabs(x-1.d0))
     . -x**2*dlog(x)+4*x*dlog(scale/msquark) )
       if(corr2.gt.-1.d0) then
       corr = dsqrt(1.d0+corr2)-1.d0
       dmsquark = msquark*corr
       else
       dmsquark =0.d0
       tachsqrc= -1.d0
       endif
       end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SU_STOPCR(pscale,mu,at,ab,m3,crLL,crLR,crRR)  !added m3
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates radiative corrections to the two stop masses, including 
c  standard and SUSY-QCD corrections and the Yukawa corrections a la PBMZ. 
c  The input at the EWSB scale are, respectively: the strong coupling 
c  constant, the gluino mass, mu parameter, pseudoscalar Higgs boson mass,
c  and trilinear couplings of the top and the bottom quarks. The outputs are 
c  the radiative corrections to the LL,LR,RR entries of the stop mass matrix. 
c---------------------------------------------------------
c
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B0,SU_BG,SU_BF
      dimension gmc(2),gmn(4),dxmn(4),u(2,2),v(2,2),z(4,4),
     .  antR(4),antL(4),bntL(4),bntR(4),actL(2),actR(2),bctL(2),bctR(2),
     .  fttLL(4),gttLL(4),fttLR(4),gttLR(4),fttRR(4),gttRR(4), 
     .  fbtLL(2),gbtLL(2),fbtLR(2),gbtLR(2),fbtRR(2),gbtRR(2)   
      common/SU_runhiggsewsb/ma,ml,mh,mch,alfa
      COMMON/SU_outginos/mc1,mc2,mn1,mn2,mn3,mn4,mgluino
      COMMON/SU_matino/u,v,z,dxmn
      COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_fmasses/mtau,mb,mt
      COMMON/SU_renscale/scale     
      COMMON/SU_bpew/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     .      msta1,msta2,msb1,msb2,mst1,mst2,thet,theb,thel
      COMMON/SU_cte/xnf,cpi,mz,mw,tbeta
      COMMON/SU_stepwi/wistep,h1,kpole
      COMMON/run_p/pizz        

c   fix ren. scale (used in B0 functions):
      if(kpole.eq.1) scale = dsqrt(mst1*mst2)     
c
      if(mst1.eq.0.d0.or.msb1.eq.0.d0) goto 100
c   (NB protection: means mst1 or msb1 undefined yet)        
       sq2=dsqrt(2.d0)
       pi = 4*datan(1.d0)
       g=g2ewsb
       g1=g1ewsb
       alphas=alsewsb
       cw = 1.d0/dsqrt(1.d0+(g1ewsb/g2ewsb)**2)
       sw = g1ewsb/g2ewsb*cw
       sw2=sw**2
       cw2= cw**2
       cwm2 =1.d0/cw2
       e=g1*cw
c
       vd2 = 2*(mz**2+pizz)/(g1ewsb**2+g2ewsb**2)/(1.d0+tbeta**2)
       vu2 = vd2*tbeta**2
       vd= dsqrt(vd2)
       vu= dsqrt(vu2)
       rmb = yb*vd
       rmt = yt*vu
       rmtau = ytau*vd
c
       zero=1.d-2
       gmn(1)=dabs(dxmn(1))
       gmn(2)=dabs(dxmn(2))
       gmn(3)=dabs(dxmn(3))
       gmn(4)=dabs(dxmn(4))
       gmc(1)=mc1
       gmc(2)=mc2
c
       B=datan(tbeta)
       cbeta2=1.d0/(1.d0+tbeta**2)
       cbet= dsqrt(cbeta2)
       sbet=dsqrt(1.d0-cbeta2)
       c2b =2*cbeta2-1.d0
       sal=dsin(alfa)
       cal=dcos(alfa)
       s2a = 2*sal*cal 
       c2a=dcos(2*alfa)
c       
       ct=dcos(thet)
       st=dsin(thet)
       cb=dcos(theb)
       sb=dsin(theb)
       cta=dcos(thel)
       sta=dsin(thel)
c
c----------- Higgs couplings
c 
       s2tLtL = -g*mz/cw*(.5d0 -2*sw2/3)*sbet + sq2*yt*rmt
       s2tRtR = -g*mz/cw*(2*sw2/3)*sbet + sq2*yt*rmt
       s2tLtR = yt/sq2*At
       s1tLtL = g*mz/cw*(.5d0 -2*sw2/3)*cbet
       s1tRtR = g*mz/cw*(2*sw2/3)*cbet
       s1tLtR = -yt/sq2*mu
c
       s2tLt1=ct*s2tLtL+st*s2tLtR
       s2tLt2=-st*s2tLtL+ct*s2tLtR
       s2tRt1=ct*s2tLtR+st*s2tRtR
       s2tRt2=-st*s2tLtR+ct*s2tRtR
       s1tLt1=ct*s1tLtL+st*s1tLtR
       s1tLt2=-st*s1tLtL+ct*s1tLtR
       s1tRt1=ct*s1tLtR+st*s1tRtR
       s1tRt2=-st*s1tLtR+ct*s1tRtR
c
       ghtLt1= cal*s1tLt1+sal*s2tLt1
       gltLt1=-sal*s1tLt1+cal*s2tLt1
       ghtLt2= cal*s1tLt2+sal*s2tLt2
       gltLt2=-sal*s1tLt2+cal*s2tLt2
       ghtRt1= cal*s1tRt1+sal*s2tRt1
       gltRt1=-sal*s1tRt1+cal*s2tRt1
       ghtRt2= cal*s1tRt2+sal*s2tRt2
       gltRt2=-sal*s1tRt2+cal*s2tRt2
                    
c
       atLtR=-yt/sq2*(-mu*sbet-At*cbet)
       gtLtR=+yt/sq2*(-mu*cbet+At*sbet)
c
       gatLt1=st*atLtR
       gatLt2=ct*atLtR
       gatRt1=-ct*atLtR
       gatRt2=st*atLtR
c
       ggtLt1=st*gtLtR
       ggtLt2=ct*gtLtR
       ggtRt1=-ct*gtLtR
       ggtRt2=st*gtLtR
c
       gctLbL = g*mw/sq2*dsin(2*b)-yt*rmt*cbet-yb*rmb*sbet
       gctRbR = -yt*rmb*cbet-yb*rmt*sbet
       gctLbR = yb*(-mu*cbet-Ab*sbet)
       gctRbL = yt*(-mu*sbet-At*cbet)
c
       ggtLbL = -g*mw/sq2*dcos(2*b)-yt*rmt*sbet+yb*rmb*cbet
       ggtRbR = 0.d0
       ggtLbR = yb*(-mu*sbet+Ab*cbet)
       ggtRbL = -yt*(-mu*cbet+At*sbet)
c
       gctLb1=cb*gctLbL+sb*gctLbR ! Corrected by P. Slavich
       gctLb2=-sb*gctLbL+cb*gctLbR
       gctRb1=cb*gctRbL+sb*gctRbR
       gctRb2=-sb*gctRbL+cb*gctRbR
c
       ggtLb1=cb*ggtLbL+sb*ggtLbR ! Corrected by P. Slavich
       ggtLb2=-sb*ggtLbL+cb*ggtLbR
       ggtRb1=cb*ggtRbL+sb*ggtRbR
       ggtRb2=-sb*ggtRbL+cb*ggtRbR
c----------- neutralino/chargino couplings:
c
       ap1tL = 0.d0
       bp1tL = g1/dsqrt(2.d0)*(1.d0/3.d0)
       ap1tR = g1/dsqrt(2.d0)*(-4.d0/3.d0)
       bp1tR = 0.d0
       ap2tL = 0.d0
       bp2tL = dsqrt(2.d0)*g*(.5d0)
       ap2tR = 0.d0
       bp2tR = 0.d0
       ap3tL = 0.d0
       ap3tR = 0.d0
       bp3tL = 0.d0
       bp3tR = 0.d0
       ap4tL = yt
       ap4tR = 0.d0
       bp4tL = 0.d0
       bp4tR = yt
c
       aw1tL=g
       bw1tL=0.d0
       aw1tR=0.d0
       bw1tR=0.d0
       aw2tL=0.d0
       bw2tL=-yb
       aw2tR=-yt
       bw2tR=0.d0
      
       do i=1,4
       aNtR(i) = Z(i,1)*ap1tR +Z(i,2)*ap2tR +Z(i,3)*ap3tR +Z(i,4)*ap4tR	  
       bNtR(i) = Z(i,1)*bp1tR +Z(i,2)*bp2tR +Z(i,3)*bp3tR +Z(i,4)*bp4tR	  
       aNtL(i) = Z(i,1)*ap1tL +Z(i,2)*ap2tL +Z(i,3)*ap3tL +Z(i,4)*ap4tL	  
       bNtL(i) = Z(i,1)*bp1tL +Z(i,2)*bp2tL +Z(i,3)*bp3tL +Z(i,4)*bp4tL	  
       enddo
c
       do i=1,2
       aCtR(i) = V(i,1)*aw1tR +V(i,2)*aw2tR   
       bCtR(i) = U(i,1)*bw1tR +U(i,2)*bw2tR   
       aCtL(i) = V(i,1)*aw1tL +V(i,2)*aw2tL  
       bCtL(i) = U(i,1)*bw1tL +U(i,2)*bw2tL  
       enddo
c
       do i=1,4
       fttLL(i) = aNtL(i)*aNtL(i) + bNtL(i)*bNtL(i) 
       gttLL(i) = bNtL(i)*aNtL(i) + aNtL(i)*bNtL(i) 
       fttLR(i) = aNtL(i)*aNtR(i) + bNtL(i)*bNtR(i) 
       gttLR(i) = bNtL(i)*aNtR(i) + aNtL(i)*bNtR(i) 
       fttRR(i) = aNtR(i)*aNtR(i) + bNtR(i)*bNtR(i) 
       gttRR(i) = bNtR(i)*aNtR(i) + aNtR(i)*bNtR(i) 
       enddo      
c 
       do i=1,2
       fbtLL(i) = aCtL(i)*aCtL(i) + bCtL(i)*bCtL(i) 
       gbtLL(i) = bCtL(i)*aCtL(i) + aCtL(i)*bCtL(i) 
       fbtLR(i) = aCtL(i)*aCtR(i) + bCtL(i)*bCtR(i) 
       gbtLR(i) = bCtL(i)*aCtR(i) + aCtL(i)*bCtR(i) 
       fbtRR(i) = aCtR(i)*aCtR(i) + bCtR(i)*bCtR(i) 
       gbtRR(i) = bCtR(i)*aCtR(i) + aCtR(i)*bCtR(i) 
       enddo      
c       
c-------------------- LL contribution: 
       crLLqcd=16*pi*alphas/3*(2*dble(SU_BG(pscale**2,m3,rmt))
     . +ct**2*(dble(SU_BF(pscale**2,mst1,zero))+SU_A(mst1) )
     . +st**2*(dble(SU_BF(pscale**2,mst2,zero))+SU_A(mst2) ) )
c       
       crLLyuk=yt**2*(st**2*SU_A(mst1)+ct**2*SU_A(mst2) )
     .        +yb**2*(sb**2*SU_A(msb1)+cb**2*SU_A(msb2) )     
     . +0.5d0*(yt**2*sal**2-g**2*(.5d0-2*sw2/3)/2.d0/cw**2*
     .        (-dcos(2*alfa)) )*SU_A(mh) 
     . +0.5d0*(yt**2*cal**2-g**2*(.5d0-2*sw2/3)/2.d0/cw**2*
     .        (+dcos(2*alfa)) )*SU_A(ml) 
     . +0.5d0*(yt**2*sbet**2-g**2*(.5d0-2*sw2/3)/2.d0/cw**2*
     .        (-dcos(2*b))    )*SU_A(mz) 
     . +0.5d0*(yt**2*cbet**2-g**2*(.5d0-2*sw2/3)/2.d0/cw**2*
     .        (dcos(2*b))     )*SU_A(ma) 
     . +(yb**2*sbet**2+g**2*( (.5d0-2*sw2/3)/2.d0/cw**2-.5d0)*
     .        (-dcos(2*b)) )*SU_A(mch)
     . +(yb**2*cbet**2+g**2*( (.5d0-2*sw2/3)/2.d0/cw**2-.5d0)*
     .        (dcos(2*b)) )*SU_A(mw)
     . +ghtLt1**2*dble(SU_B0(pscale**2,mh,mst1))
     . +ghtLt2**2*dble(SU_B0(pscale**2,mh,mst2))
     . +gltLt1**2*dble(SU_B0(pscale**2,ml,mst1))
     . +gltLt2**2*dble(SU_B0(pscale**2,ml,mst2))
     . +ggtLt1**2*dble(SU_B0(pscale**2,mz,mst1))
     . +ggtLt2**2*dble(SU_B0(pscale**2,mz,mst2))
     . +gatLt1**2*dble(SU_B0(pscale**2,ma,mst1))
     . +gatLt2**2*dble(SU_B0(pscale**2,ma,mst2))          
     . +gctLb1**2*dble(SU_B0(pscale**2,mch,msb1))
     . +gctLb2**2*dble(SU_B0(pscale**2,mch,msb2))
     . +ggtLb1**2*dble(SU_B0(pscale**2,mw,msb1))
     . +ggtLb2**2*dble(SU_B0(pscale**2,mw,msb2))
c
      crLLgau=4*g**2/cw**2*(.5d0-2*sw2/3)**2*SU_A(mz)+2*g**2*SU_A(mw)
     .  +(2*g1*cw/3)**2*(ct**2*dble(SU_BF(pscale**2,mst1,zero))
     .                  +st**2*dble(SU_BF(pscale**2,mst2,zero)) )
     .  +g**2/cw**2*(.5d0-2*sw2/3)**2*(
     .                   ct**2*dble(SU_BF(pscale**2,mst1,mz))
     .                  +st**2*dble(SU_BF(pscale**2,mst2,mz)) )    
     .  +g**2*( cb**2*dble(SU_BF(pscale**2,msb1,mw))
     .         +sb**2*dble(SU_BF(pscale**2,msb2,mw)) )    
     .  +g**2/4*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) 
     .       +2*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) ) )       
c
      crLLhyp=g**2*0.5d0*(
     .        3.d0*(+.5d0)*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) )    
     .       +3.d0*(-.5d0)*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) )    
     .       +(-.5d0)*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) )    
     .       +6.d0*(+.5d0)*SU_A(msu1)+6.d0*(-.5d0)*SU_A(msd1)
     .       +2.d0*(-.5d0)*SU_A(mse1)
     .       +2.d0*(.5d0)*SU_A(msn1) +(.5d0)*SU_A(msntau) )
     .  +g1**2/4*(1.d0/3.d0)**2*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) ) 
     .  +g1**2/4*(1.d0/3.d0)*(
     .        3.d0*(1.d0/3.d0)*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) )    
     .       +3.d0*(1.d0/3.d0)*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) )    
     .       +(-1.d0)*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) )    
     .       +6.d0*(1.d0/3.d0)*SU_A(msu1)+6.d0*(1.d0/3.d0)*SU_A(msd1)
     .       +2.d0*(-1.d0)*SU_A(mse1)
     .       +2.d0*(-1.d0)*SU_A(msn1) + (-1.d0)*SU_A(msntau) )
     .  +g1**2/4*(1.d0/3.d0)*(
     .        3.d0*(-4.d0/3.d0)*(st**2*SU_A(mst1)+ct**2*SU_A(mst2) )    
     .       +3.d0*(2.d0/3.d0)*(sb**2*SU_A(msb1)+cb**2*SU_A(msb2) )    
     .       +(2.d0)*(sta**2*SU_A(msta1)+cta**2*SU_A(msta2) )    
     .       +6.d0*(-4.d0/3.d0)*SU_A(msu2)+6.d0*(2.d0/3.d0)*SU_A(msd2)
     .       +2.d0*(2.d0)*SU_A(mse2)  )
c     
      crLLnino=0.d0    
      do i=1,4
      crLLnino=crLLnino+ fttLL(i)*dble(SU_BG(pscale**2,gmn(i),rmt))
     .    -2.d0*rmt*dxmn(i)*gttLL(i)*dble(SU_B0(pscale**2,gmn(i),rmt))
      enddo
c
      crLLcino=0.d0
      do i=1,2
      crLLcino=crLLcino+ fbtLL(i)*dble(SU_BG(pscale**2,gmc(i),rmb))
     .    -2.d0*rmb*gmc(i)*gttLL(i)*dble(SU_B0(pscale**2,gmc(i),rmb))     
      enddo 
c
      crLL=-cpi*(crLLqcd+crLLyuk+crLLgau+crLLhyp+crLLnino+crLLcino)	    
c       
c-------------------- RR contribution: 
c
       crRRqcd=16*pi*alphas/3*(2*dble(SU_BG(pscale**2,m3,rmt))
     . +st**2*(dble(SU_BF(pscale**2,mst1,zero))+SU_A(mst1) )
     . +ct**2*(dble(SU_BF(pscale**2,mst2,zero))+SU_A(mst2) ) )
c       
       crRRyuk=yt**2*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) )
     .        +yt**2*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) )     
     . +0.5d0*(yt**2*sal**2-g**2*(2*sw2/3)/2.d0/cw**2*
     .        (-dcos(2*alfa)) )*SU_A(mh) 
     . +0.5d0*(yt**2*cal**2-g**2*(2*sw2/3)/2.d0/cw**2*
     .        (+dcos(2*alfa)) )*SU_A(ml) 
     . +0.5d0*(yt**2*sbet**2-g**2*(2*sw2/3)/2.d0/cw**2*
     .        (-dcos(2*b))    )*SU_A(mz) 
     . +0.5d0*(yt**2*cbet**2-g**2*(2*sw2/3)/2.d0/cw**2*
     .        (dcos(2*b))     )*SU_A(ma) 
     . +(yt**2*cbet**2+g**2*( (2*sw2/3)/2.d0/cw**2)*
     .        (-dcos(2*b)) )*SU_A(mch)
     . +(yt**2*sbet**2+g**2*( (2*sw2/3)/2.d0/cw**2)*
     .        (dcos(2*b)) )*SU_A(mw)
     . +ghtRt1**2*dble(SU_B0(pscale**2,mh,mst1))
     . +ghtRt2**2*dble(SU_B0(pscale**2,mh,mst2))
     . +gltRt1**2*dble(SU_B0(pscale**2,ml,mst1))
     . +gltRt2**2*dble(SU_B0(pscale**2,ml,mst2))
     . +ggtRt1**2*dble(SU_B0(pscale**2,mz,mst1))
     . +ggtRt2**2*dble(SU_B0(pscale**2,mz,mst2))
     . +gatRt1**2*dble(SU_B0(pscale**2,ma,mst1))
     . +gatRt2**2*dble(SU_B0(pscale**2,ma,mst2))          
     . +gctRb1**2*dble(SU_B0(pscale**2,mch,msb1))
     . +gctRb2**2*dble(SU_B0(pscale**2,mch,msb2))
     . +ggtRb1**2*dble(SU_B0(pscale**2,mw,msb1))
     . +ggtRb2**2*dble(SU_B0(pscale**2,mw,msb2))
c
      crRRgau=4*g**2/cw**2*(2*sw2/3)**2*SU_A(mz)
     .  +(2*g1*cw/3)**2*(st**2*dble(SU_BF(pscale**2,mst1,zero))
     .                  +ct**2*dble(SU_BF(pscale**2,mst2,zero)) )
     .  +g**2/cw**2*(2*sw2/3)**2*(
     .                   st**2*dble(SU_BF(pscale**2,mst1,mz))
     .                  +ct**2*dble(SU_BF(pscale**2,mst2,mz)) )    
c
      crRRhyp=
     .   g1**2/4*(-4.d0/3.d0)**2*(st**2*SU_A(mst1)+ct**2*SU_A(mst2) ) 
     .  +g1**2/4*(-4.d0/3.d0)*(
     .        3.d0*(1.d0/3.d0)*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) )    
     .       +3.d0*(1.d0/3.d0)*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) )    
     .       +(-1.d0)*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) )    
     .       +6.d0*(1.d0/3.d0)*SU_A(msu1)+6.d0*(1.d0/3.d0)*SU_A(msd1)
     .       +2.d0*(-1.d0)*SU_A(mse1)
     .       +2.d0*(-1.d0)*SU_A(msn1) +(-1.d0)*SU_A(msntau) )
     .  +g1**2/4*(1.d0/3.d0)*(
     .        3.d0*(-4.d0/3.d0)*(st**2*SU_A(mst1)+ct**2*SU_A(mst2) )    
     .       +3.d0*(2.d0/3.d0)*(sb**2*SU_A(msb1)+cb**2*SU_A(msb2) )    
     .       +(2.d0)*(sta**2*SU_A(msta1)+cta**2*SU_A(msta2) )    
     .       +6.d0*(-4.d0/3.d0)*SU_A(msu2)+6.d0*(2.d0/3.d0)*SU_A(msd2)
     .       +2.d0*(2.d0)*SU_A(mse2)  )
c     
      crRRnino=0.d0    
      do i=1,4
      crRRnino=crRRnino+ fttRR(i)*dble(SU_BG(pscale**2,gmn(i),rmt))
     .    -2.d0*rmt*dxmn(i)*gttRR(i)*dble(SU_B0(pscale**2,gmn(i),rmt))
      enddo
c
      crRRcino=0.d0
      do i=1,2
      crRRcino=crRRcino+ fbtRR(i)*dble(SU_BG(pscale**2,gmc(i),rmb))
     .    -2.d0*rmb*gmc(i)*gttRR(i)*dble(SU_B0(pscale**2,gmc(i),rmb))     
      enddo 
c
      crRR=-cpi*(crRRqcd+crRRyuk+crRRgau+crRRhyp+crRRnino+crRRcino)	    
c       
c-------------------- LR contribution: 
c
       crLRqcd=16*pi*alphas/3*(
     . 4*rmt*m3*dble(SU_B0(pscale**2,m3,rmt))
     . +ct*st*(dble(SU_BF(pscale**2,mst1,zero))-SU_A(mst1) 
     .        -dble(SU_BF(pscale**2,mst2,zero))+SU_A(mst2) ) )
c       
       crLRyuk=3*yt**2*st*ct*(SU_A(mst1)-SU_A(mst2) )
     . +ghtLt1*ghtRt1*dble(SU_B0(pscale**2,mh,mst1))
     . +ghtLt2*ghtRt2*dble(SU_B0(pscale**2,mh,mst2))
     . +gltLt1*gltRt1*dble(SU_B0(pscale**2,ml,mst1))
     . +gltLt2*gltRt2*dble(SU_B0(pscale**2,ml,mst2))
     . +ggtLt1*ggtRt1*dble(SU_B0(pscale**2,mz,mst1))
     . +ggtLt2*ggtRt2*dble(SU_B0(pscale**2,mz,mst2))
     . +gatLt1*gatRt1*dble(SU_B0(pscale**2,ma,mst1))
     . +gatLt2*gatRt2*dble(SU_B0(pscale**2,ma,mst2))          
     . +gctLb1*gctRb1*dble(SU_B0(pscale**2,mch,msb1))
     . +gctLb2*gctRb2*dble(SU_B0(pscale**2,mch,msb2))
     . +ggtLb1*ggtRb1*dble(SU_B0(pscale**2,mw,msb1))
     . +ggtLb2*ggtRb2*dble(SU_B0(pscale**2,mw,msb2))
c
      crLRgau=(2*g1*cw/3)**2*ct*st*(dble(SU_BF(pscale**2,mst1,zero))
     .                          -dble(SU_BF(pscale**2,mst2,zero)) )
     .  -g**2/cw**2*(.5d0-2*sw2/3)*(2*sw2/3)*st*ct*(
     .                           dble(SU_BF(pscale**2,mst1,mz))
     .                          -dble(SU_BF(pscale**2,mst2,mz)) )    
c
      crLRhyp=g1**2/4*(1.d0/3.d0)*(-4.d0/3.d0)*st*ct*(
     .                SU_A(mst1)-SU_A(mst2) ) 
c     
      crLRnino=0.d0    
      do i=1,4
      crLRnino=crLRnino+ fttLR(i)*dble(SU_BG(pscale**2,gmn(i),rmt))
     .    -2.d0*rmt*dxmn(i)*gttLR(i)*dble(SU_B0(pscale**2,gmn(i),rmt))
      enddo
c
      crLRcino=0.d0
      do i=1,2
      crLRcino=crLRcino+ fbtLR(i)*dble(SU_BG(pscale**2,gmc(i),rmb))
     .    -2.d0*rmb*gmc(i)*gttLR(i)*dble(SU_B0(pscale**2,gmc(i),rmb))     
      enddo 
c
      crLR=-cpi*(crLRqcd+crLRyuk+crLRgau+crLRhyp+crLRnino+crLRcino)	              
c
 100  continue
      end
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE SU_GINOCR(alphas,m3, mb,mt, delgino)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Calculates the radiative correction to the gluino mass, delgino.
c  The input parameters at EWSB scale are: 
c  alphas,m3: the strong coupling constant and the SU(3) gaugino mass,
c  msu1,msu2,msd1,msd2,msb1,msb2,mst1,mst2: the squark masses.       
c
      implicit real*8(a-h,m,o-z)
      complex*16 SU_B1,SU_B0
      COMMON/SU_renscale/scale
      COMMON/SU_param/gf,alpha,mz,mw
      COMMON/SU_bpew/msu1,msu2,msd1,msd2,
     . mse1,mse2,msn1,msntau,
     . msta1,msta2,msb1,msb2,mst1,mst2,
     . thet,theb,thel
      COMMON/SU_stepwi/wistep,h1,kpole
c
      if(kpole.eq.1) scale = dsqrt(mst1*mst2)
       pi=4*datan(1.d0)
       mu=.005d0
       md=.015d0
       ms=.19d0
       mc=1.40d0
       msc1=msu1
       msc2=msu2
       mss1=msd1
       mss2=msd2
c       
       sumB1= 
     .        dble(SU_B1(M3**2,mu,msu1))+ dble(SU_B1(M3**2,mu,msu2) )
     .       +dble(SU_B1(M3**2,md,msd1) )+dble(SU_B1(M3**2,md,msd2) )
     .       +dble(SU_B1(M3**2,ms,mss1) )+dble(SU_B1(M3**2,ms,mss2) )
     .       +dble(SU_B1(M3**2,mc,msc1) )+dble(SU_B1(M3**2,mc,msc2) )
     .       +dble(SU_B1(M3**2,mb,msb1) )+dble(SU_B1(M3**2,mb,msb2) )
     .       +dble(SU_B1(M3**2,mt,mst1) )+dble(SU_B1(M3**2,mt,mst2) )
c
       sumB0= mb*dsin(2*theb)*
     .       (dble(SU_B0(M3**2,mb,msb1))-dble(SU_B0(M3**2,mb,msb2)) )
     .       +mt*dsin(2*thet)*        
     .       (dble(SU_B0(M3**2,mt,mst1))-dble(SU_B0(M3**2,mt,mst2)) )  
c
        delgino =3*alphas/pi/4*M3*(5.d0-3*dlog(M3**2/scale**2))
     .          -alphas/pi/4*M3*sumB1
     .          -alphas/pi/4*sumB0
       end
c
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   +++++++ End of the routines for gaugino sfermion massses ++++++++++++
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  The following routine is for the evaluation of the Higgs boson masses:
cc  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
                 SUBROUTINE SU_SUSYCP(TGBET)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c  Calculates the MSSM Higgs bosons masses and the angle alpha including
c  radiative corrections for a given input value of the parameter tan(beta).       
c  The other input parameters (soft-SUSY breaking parameters, sparticle
c  masses and mixing angles, SM parameters, are called via common blocks. 
c  It returns the masses of the pole masses of the CP-odd (ama), lighter 
c  CP-even (aml), heavier CP-even (amh), charged Higgs boson (amch) as well as 
c  the running CP-odd (amar) Higgs masses, which are given in the block: 
c                common/su_HMASS/ama,aml,amh,amch,amar.
c  It gives also the couplings of the angle beta at the EWSB scale, the mixing
c  alpha and the Higgs boson couplings to standard particles in:
c        COMMON/COUP_hcoup/gat,gab,glt,glb,ght,ghb,glvv,ghvv,b,a
c  It returns also the couplings of the Higgs bosons to sfermions 
c        COMMON/SU_cplhsf/gcen,gctb,glee,gltt,glbb,ghee,ghtt,ghbb
c      .                   gatt,gabb,gaee
c  and the Higgs couplings to charginos and neutralinos:
c         COMMON/SU_cplhino/ac1,ac2,ac3,an1,an2,an3,acnl,acnr
c  For the radiative correction of Higgs masses, there is imodel=0
c  option where the calculation is made in an approximation based on the
c  of work Heinemeyer, Hollik, Weiglein (hep-ph/0002213), which is fast
c  BUT approximate, or 
c  ichoice(10)(=imodel)=1:  Full one-loop Higgs masses  a la PBMZ 
c                      =2:  Full one-loop PBMZ + two-loop BDSZ corrections
c
      implicit double precision (a-h,m,o-z)
      double precision la1,la2,la3,la4,la5,la6,la7,la3t
      complex*16 F0_TLH
      logical su_isNaN
      dimension mst(2),gltt(2,2),ghtt(2,2),
     .          msb(2),glbb(2,2),ghbb(2,2),
     .          msl(2),glee(2,2),ghee(2,2),
     .          gctb(2,2),gcen(2,2)
      dimension dxmn(4),z(4,4),uu(2,2),vv(2,2),
     .          qqn(4,4),ssn(4,4),ssc(2,2),qqc(2,2),ac1(2,2),ac2(2,2),
     .          ac3(2,2),an1(4,4),an2(4,4),an3(4,4),acnl(2,4),acnr(2,4)
      COMMON/SU_hflag/imodel
      COMMON/SU_fmasses/amtau,amb,amt
      COMMON/SU_hmass/ama,aml,amh,amch,amar
      COMMON/SU_break/amel,amer,amsq,amur,amdr,al,au,ad,amu,am1,am2,am3
      COMMON/SU_mssmhpar/mhu2,mhd2,madum,mudum
      COMMON/SU_param/gf,alph,amz,amw
      COMMON/SU_cplhsf/gcen,gctb,glee,gltt,glbb,ghee,ghtt,ghbb,
     .                 gatt,gabb,gaee
      COMMON/SU_matino/uu,vv,z,dxmn
      COMMON/SU_cplhino/ac1,ac2,ac3,an1,an2,an3,acnl,acnr  
      COMMON/SU_hcoup/tbdum,alfa,xgat,xgab,xglt,xglb,xght,xghb,
     .                xghvv,xglvv
c  Commons needed for the full one+two--loop calculation
      COMMON/SU_bpew/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     .      msta1,msta2,msb1,msb2,mst1,mst2,thett,thetb,thetl
      COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_tbewsb/vuewsb,vdewsb 
c      
c  Commons needed for interface with routines (from HDECAY3.0):
          COMMON/HSELF_HDEC/LA1,LA2,LA3,LA4,LA5,LA6,LA7
          COMMON/HMASSR_HDEC/AMLR,AMHR
          COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A     

      common/su_MAinput/piaa,tadba,DMA,kMAflag    
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
c  commons added
          COMMON/pietro/mApole,mCHpole  
          COMMON/SU_renscale/scale 
          COMMON/run_p/pizzp     
          COMMON/runhiggs/ama0,aml0,amh0,amch0
      COMMON/SU_runhiggsewsb/marunp,mlrunp,mhrunp,mchrunp,alfarunp
      common/su_savemar/madr2save
      common/su_runmavev/madr2,vev2
      common/su_runewsb/rmz,rmw,sw2,tbeta
c
c Some definitions:
      pi=4*datan(1.d0)
      v=1.d0/dsqrt(dsqrt(2.d0)*gf)
      tbeta=vuewsb/vdewsb
      bet=datan(tbeta)
      if(b.eq.0.d0) b=bet
      sb = dsin(bet)
      cb = dcos(bet)
      amar = ama
      marsave=ama
c nb at this stage ama is in fact running MA
      mt=runm(amt,6)
      mb=runm(amt,5)
      als=alphas(amt,2)
      sw2 = g1ewsb**2/(g1ewsb**2+g2ewsb**2)
c
C =================  Calculate the masses in an approximation =======
c (based on Heinemeyer, Hollik, Weiglein hep-ph/0002213 )
c====================================================================
      if(imodel.eq.0) then
      vev2=2.d0*(vuewsb**2+vdewsb**2)
      madr2= marsave**2
      endif
      amglu=am3
      ams2=dsqrt(amsq**2*amur**2+mt**2*(amsq**2+amur**2)+mt**4)
      xlam=1.d0/8.d0-sw2/3.d0+4*sw2**2/9.d0
      xt=au-amu/tbeta
      xr=mt**2/ams2
      xfac=gf*dsqrt(2.d0)/pi**2
      s11=xfac*amz**4*xlam*cb**2*dlog(xr)
      s12=-xfac*amz**2/tbeta*(-3*mt**2/8.d0+amz**2*xlam*sb**2)
     .     *dlog(xr)
      s22one=xfac*mt**4/8.d0/sb**2*(
     .     -2*amz**2/mt**2+11*amz**4/10.d0/mt**4
     .     +(12.d0-6*amz**2/mt**2*sb**2+8*amz**4/mt**4
     .     *xlam*sb**4)*dlog(xr)
     .  +xt**2/ams2*(-12.d0+4.d0*amz**2/mt**2+6.d0*xr)
     .  +xt**4/ams2**2*(1.d0-4*xr+3*xr**2)
     .  +xt**6/ams2**3*(3*xr/5.d0-12*xr**2/5.d0+2*xr**3)
     .  +xt**8/ams2**4*(3*xr**2/7.d0-12*xr**3/7.d0+3*xr**4/2.d0) )
      s22qcd=xfac*als/pi*mt**4/sb**2*(4.d0+3*dlog(xr)**2
     .  +2*dlog(xr)-6*xt/dsqrt(ams2)-xt**2/ams2*(3*dlog(xr)+8.d0)
     .  +17*xt**4/12.d0/ams2**2 )
      s22ew=-9*xfac**2/32.d0/sb**2*mt**6*(dlog(xr)**2
c !  next term changes sign (-2*xt**2/ams2->+2 xt..) wrt ref.: corrected: 
     .   +2*xt**2/ams2*dlog(xr)+xt**4/6.d0/ams2**2*dlog(xr))            
      s22=s22one+s22qcd+s22ew     
      
      xm11=ama**2*sb**2+amz**2*cb**2-s11
      xm12=-(ama**2+amz**2)*sb*cb-s12    
      xm22=ama**2*cb**2+amz**2*sb**2-s22       
      xml2=0.5d0*(xm11+xm22-dsqrt((xm11-xm22)**2+4*xm12**2))
      xmh2=0.5d0*(xm11+xm22+dsqrt((xm11-xm22)**2+4*xm12**2))    
      amlr=dsqrt(xml2)
      amhr=dsqrt(xmh2)
      ama=amar
      a=datan(xm12/(amz**2*cb**2+ama**2*sb**2-s11-aml**2))
c      
      sa=dsin(a)
      ca=dcos(a)
      if(ca.eq.0)then
       a = pi/2
      else
       a=datan(sa/ca)
      endif
      if(ca.lt.0d0)then
       if(sa.lt.0d0)then
        a = a-pi
       else
        a = a+pi
       endif
      endif
      sa=dsin(a)
      ca=dcos(a) 
c =====================================================================
C===== Now calculate the Higgs boson coupling to sfermions and gauginos:
C =====================================================================
      sbma = sb*ca-cb*sa
      cbma = cb*ca+sb*sa
      sbpa = sb*ca+cb*sa
      cbpa = cb*ca-sb*sa
       mstl2=amsq**2+(0.5d0-2.d0/3.d0*sw2)*amz**2*dcos(2.d0*b)
       mstr2=amur**2+2.d0/3.d0*sw2*amz**2*dcos(2.d0*b)
       mlrt=au-amu/tgbet
       delt=(mstl2-mstr2)**2+4*mt**2*mlrt**2
       mst12=mt**2+0.5d0*(mstl2+mstr2-dsqrt(delt))
       mst22=mt**2+0.5d0*(mstl2+mstr2+dsqrt(delt))
        if(mst12.lt.0.d0.and.imodel.eq.0)goto 111
       mst(1)=dsqrt(mst12)
       mst(2)=dsqrt(mst22)
       if(mstl2.eq.mstr2) then
        thet = pi/4
       else
        thet=0.5d0*datan(2.d0*mt*mlrt / (mstl2-mstr2) )
        if(mstl2.gt.mstr2) thet = thet + pi/2
       endif
       cst= dcos(thet)
       sst= dsin(thet)
c===== sbottom masses
       msbl2=amsq**2+(-0.5d0+1.d0/3.d0*sw2)*amz**2*dcos(2.d0*b)
       msbr2=amdr**2-1.d0/3.d0*sw2*amz**2*dcos(2.d0*b)
       mlrb=ad-amu*tgbet
       delb=(msbl2-msbr2)**2+4*mb**2*mlrb**2
       msb12=mb**2+0.5d0*(msbl2+msbr2-dsqrt(delb))
       msb22=mb**2+0.5d0*(msbl2+msbr2+dsqrt(delb))
        if(msb12.lt.0.d0.and.imodel.eq.0)goto 111
       msb(1)=dsqrt(msb12)
       msb(2)=dsqrt(msb22)
       if(msbl2.eq.msbr2) then
        theb = pi/4
       else
        theb=0.5d0*datan(2.d0*mb*mlrb / (msbl2-msbr2) )
        if(msbl2.gt.msbr2) theb = theb + pi/2
       endif
       csb= dcos(theb)
       ssb= dsin(theb)
c===== stau masses
      msel2=amel**2+(-0.5d0+sw2)*amz**2*dcos(2.d0*b)
      mser2=amer**2- sw2*amz**2*dcos(2.d0*b) 
      msn2=amel**2+0.5d0*amz**2*dcos(2.d0*b)
      mlre=al-amu*tgbet
      dele=(msel2-mser2)**2+4*amtau**2*mlre**2
      mse12=amtau**2+0.5d0*(msel2+mser2-dsqrt(dele))
      mse22=amtau**2+0.5d0*(msel2+mser2+dsqrt(dele))
        if(mse12.lt.0.d0.and.imodel.eq.0)goto 111 
      msl(1)=dsqrt(mse12)
      msl(2)=dsqrt(mse22)
      msn   =dsqrt(msn2)
       if(msel2.eq.mser2) then
        thel = pi/4
       else
      thel=0.5d0*datan(2.d0*amtau*mlre / (msel2-mser2) )
      if(msel2.gt.mser2) thel = thel + pi/2
        endif  
      csl= dcos(thel)
      ssl= dsin(thel) 
c===== light higgs couplings to sfermions 
       glt=ca/sb
       glb=-sa/cb
       gltt(1,1)=-sbpa*(0.5d0*cst**2-2.d0/3.d0*sw2*dcos(2*thet) )
     .     +mt**2/amz**2*glt + mt*sst*cst/amz**2*(au*glt+amu*ght)
       gltt(2,2)=-sbpa*(0.5d0*sst**2+2.d0/3.d0*sw2*dcos(2*thet) )
     .     +mt**2/amz**2*glt - mt*sst*cst/amz**2*(au*glt+amu*ght)
       gltt(1,2)=-2*sbpa*sst*cst*(2.d0/3.d0*sw2-0.25d0)
     .     + mt*dcos(2*thet)/2.d0/amz**2*(au*glt+amu*ght)
       gltt(2,1)=-2*sbpa*sst*cst*(2.d0/3.d0*sw2-0.25d0)
     .     + mt*dcos(2*thet)/2.d0/amz**2*(au*glt+amu*ght)
       glbb(1,1)=-sbpa*(-0.5d0*csb**2+1.d0/3.d0*sw2*dcos(2*theb))
     .     +mb**2/amz**2*glb + mb*ssb*csb/amz**2*(ad*glb-amu*ghb)
       glbb(2,2)=-sbpa*(-0.5d0*ssb**2-1.d0/3.d0*sw2*dcos(2*theb))
     .     +mb**2/amz**2*glb - mb*ssb*csb/amz**2*(ad*glb-amu*ghb)
       glbb(1,2)=-2*sbpa*ssb*csb*(-1.d0/3.d0*sw2+0.25d0)
     .    + mb*dcos(2*theb)/2.d0/amz**2*(ad*glb-amu*ghb)
       glbb(2,1)=-2*sbpa*ssb*csb*(-1.d0/3.d0*sw2+0.25d0)
     .     + mb*dcos(2*theb)/2.d0/amz**2*(ad*glb-amu*ghb)
       glee(1,1)=-sbpa*(-0.5d0*csl**2+sw2*dcos(2*thel))
     .     +amtau**2/amz**2*glb+amtau*ssl*csl/amz**2*(al*glb-amu*ghb)
       glee(2,2)=-sbpa*(-0.5d0*ssl**2-sw2*dcos(2*thel))
     .     +amtau**2/amz**2*glb-amtau*ssl*csl/amz**2*(al*glb-amu*ghb)
       glee(1,2)=-2*sbpa*ssl*csl*(-sw2+0.25d0)
     .    + amtau*dcos(2*thel)/2.d0/amz**2*(al*glb-amu*ghb)
       glee(2,1)=-2*sbpa*ssl*csl*(-sw2+0.25d0)
     .     + amtau*dcos(2*thel)/2.d0/amz**2*(al*glb-amu*ghb)

c===== heavy higgs couplings to sfermions 
       ght=sa/sb
       ghb=ca/cb
       ghtt(1,1)=cbpa*(0.5d0*cst**2-2.d0/3.d0*sw2*dcos(2*thet))
     .     +mt**2/amz**2*ght + mt*sst*cst/amz**2*(au*ght-amu*glt)
       ghtt(2,2)=cbpa*(0.5d0*sst**2+2.d0/3.d0*sw2*dcos(2*thet))
     .     +mt**2/amz**2*ght - mt*sst*cst/amz**2*(au*ght-amu*glt)
       ghtt(1,2)=2*cbpa*sst*cst*(2.d0/3.d0*sw2-0.25d0)
     .     +mt*dcos(2*thet)/2.d0/amz**2*(au*ght-amu*glt)
       ghtt(2,1)=2*cbpa*sst*cst*(2.d0/3.d0*sw2-0.25d0)
     .     + mt*dcos(2*thet)/2.d0/amz**2*(au*ght-amu*glt)
       ghbb(1,1)=cbpa*(-0.5d0*csb**2+1.d0/3.d0*sw2*dcos(2*theb))
     .     +mb**2/amz**2*ghb + mb*ssb*csb/amz**2*(ad*ghb+amu*glb)
       ghbb(2,2)=cbpa*(-0.5d0*ssb**2-1.d0/3.d0*sw2*dcos(2*theb))
     .     + mb**2/amz**2*ghb - mb*ssb*csb/amz**2*(ad*ghb+amu*glb)
       ghbb(1,2)=2*cbpa*ssb*csb*(-1.d0/3.d0*sw2+0.25d0)
     .     + mb*dcos(2*theb)/2.d0/amz**2*(ad*ghb+amu*glb)
       ghbb(2,1)=2*cbpa*ssb*csb*(-1.d0/3.d0*sw2+0.25d0)
     .     + mb*dcos(2*theb)/2.d0/amz**2*(ad*ghb+amu*glb)
       ghee(1,1)=cbpa*(-0.5d0*csl**2+sw2*dcos(2*thel))
     .     +amtau**2/amz**2*ghb+amtau*ssl*csl/amz**2*(al*ghb+amu*glb)
       ghee(2,2)=cbpa*(-0.5d0*ssb**2-sw2*dcos(2*thel))
     .     + amtau**2/amz**2*ghb-amtau*ssl*csl/amz**2*(al*ghb+amu*glb)
       ghee(1,2)=2*cbpa*ssl*csl*(-sw2+0.25d0)
     .     + amtau*dcos(2*thel)/2.d0/amz**2*(al*ghb+amu*glb)
       ghee(2,1)=2*cbpa*ssl*csl*(-sw2+0.25d0)
     .     + amtau*dcos(2*thel)/2.d0/amz**2*(al*ghb+amu*glb)
c===== pseudoscalar higgs couplings to sfermions
       gat=1.d0/tgbet
       gab=tgbet
       gatt=-mt/2.d0/amz**2*(amu+au*gat) 
       gabb=-mb/2.d0/amz**2*(amu+ad*gab)
       gaee=-amtau/2.d0/amz**2*(amu+al*gab) 
c===== charged higgs couplings sfermions 
      cll3=(amw**2*dsin(2*b)-mt**2*gat-mb**2*gab)/dsqrt(2.d0)/amw**2
      crr3=-mt*mb*(gat+gab)/dsqrt(2.d0)/amw**2
      clr3=-mb*(amu+ad*gab)/dsqrt(2.d0)/amw**2
      crl3=-mt*(amu+au*gat)/dsqrt(2.d0)/amw**2
      gctb(1,1)=+cst*csb*cll3+sst*ssb*crr3+cst*ssb*clr3+sst*csb*crl3
      gctb(1,2)=-cst*ssb*cll3+sst*csb*crr3+cst*csb*clr3-sst*ssb*crl3
      gctb(2,1)=-sst*csb*cll3+cst*ssb*crr3-sst*ssb*clr3+cst*csb*crl3
      gctb(2,2)=+sst*ssb*cll3+cst*csb*crr3-sst*csb*clr3-cst*ssb*crl3
      cll1=(amw**2*dsin(2*b)-amtau**2*gab)/dsqrt(2.d0)/amw**2
      clr1=-amtau*(amu+al*gab)/dsqrt(2.d0)/amw**2
      gcen(1,1)=csl*cll1+ssl*clr1
      gcen(1,2)=-ssl*cll1+csl*clr1
      gcen(2,1)=0.d0
      gcen(2,2)=0.d0 
c=====  neutral higgs couplings to neutralinos
        tanw=dsqrt(sw2)/dsqrt(1.d0-sw2) 
	do 11 i=1,4
	do 11 j=1,4
	qqn(i,j)=1.d0/2.d0*(z(i,3)*(z(j,2)-tanw*z(j,1))+z(j,3)*
     .		(z(i,2)-tanw*z(i,1)))
	ssn(i,j)=1.d0/2.d0*(z(i,4)*(z(j,2)-tanw*z(j,1))+z(j,4)*
     .		(z(i,2)-tanw*z(i,1)))
 11	continue
	do 21 i=1,4
	do 21 j=1,4
	an1(i,j)= qqn(i,j)*dcos(a)-ssn(i,j)*dsin(a)
	an2(i,j)=-qqn(i,j)*dsin(a)-ssn(i,j)*dcos(a)
	an3(i,j)= qqn(i,j)*dsin(bet)-ssn(i,j)*dcos(bet)
 21	continue
c=====  neutral higgs couplings to charginos 
	do 12 i=1,2
	do 12 j=1,2
	qqc(i,j)=dsqrt(1.d0/2.d0)*uu(j,2)*vv(i,1)
	ssc(i,j)=dsqrt(1.d0/2.d0)*uu(j,1)*vv(i,2)
 12	continue
	do 22 i=1,2
	do 22 j=1,2	
	ac1(i,j)= qqc(i,j)*dcos(a)+ssc(i,j)*dsin(a)
	ac2(i,j)=-qqc(i,j)*dsin(a)+ssc(i,j)*dcos(a)	
        ac3(i,j)= qqc(i,j)*dsin(bet)+ssc(i,j)*dcos(bet)
 22	continue
c=====  charged higgs couplings to charginos-neutralinos 
	do 13 i=1,2
	do 13 j=1,4
        acnl(i,j)=dcos(bet)*(z(j,4)*vv(i,1)+(z(j,2)+z(j,1)*tanw)
     .       *vv(i,2)/dsqrt(2.d0)) 
        acnr(i,j)=dsin(bet)*(z(j,3)*uu(i,1)-(z(j,2)+z(j,1)*tanw)
     .       *uu(i,2)/dsqrt(2.d0)) 
 13     continue 
c
      if(imodel.ge.1)then
C ============ gluino and heaviest chargino mass needed for subh ======
      amchi2=am2**2+amu**2+2.d0*amw**2+dsqrt((am2**2-amu**2)**2
     .      +4.d0*amw**4*dcos(2.d0*bet)**2+4.d0*amw**2*
     .      (am2**2+amu**2+2.d0*amu*am2*dsin(2.d0*bet) ) ) 
      amchi=dsqrt(0.5d0*amchi2)
      amglu=am3
c
c--use Carena et al. for some things not included in the Higgs routine:
       CALL SUBH_TLH(ama,tgbet,amsq,amur,amdr,amt,au,ad,amu,amchi,
     .            amlr,amhr,amchr,sa,ca,tanba,amglu)

c- Now call the routine for the full one-loop or 2-loop calculation:
c=======================================================================
       q2 = scale**2
c
      tbeta = vuewsb/vdewsb
      if(su_isNaN(pizzp).or.amz**2+pizzp.le.0.d0) then
c !!! protections added 
c non-pert or NaN pb, uses tree-level values temporarily:
      pizzp = 0.d0
      if(irge.eq.irgmax) inonpert=-1    
      endif
      rmz=dsqrt(amz**2+pizzp)  
      rmw= rmz*dsqrt(1d0-sw2) 
       ama = marsave
       if(kmaflag.eq.0.and.irge.ge.2) then  
       if(madr2save.gt.0.d0) ama =dsqrt(madr2save)
       amar=ama
       endif
c
       amdelta=(ama**2+rmz**2)**2-4.d0*ama**2*rmz**2*dcos(2.d0*bet)**2
       aml0=dsqrt(0.5d0*(ama**2+rmz**2-dsqrt(amdelta)))
       amh0=dsqrt(0.5d0*(ama**2+rmz**2+dsqrt(amdelta)))
       amch0=dsqrt(ama**2+rmw**2)
c defining running higgs masses here for common (added):
       marunp = ama
       mlrunp = aml0
       mhrunp = amh0
       mchrunp = amch0
       alfarunp =0.5*atan(tan(2d0*bet)*(ama**2+rmz**2)/(ama**2-rmz**2))
       if(cos(2d0*bet)*(ama**2-rmz**2).gt.0)
     $      alfarunp = alfarunp - pi/2d0
c
       if(aml.eq.0.d0) then
       aml=aml0
       mlpole = aml0
       endif      
       if(amh.eq.0.d0) then
       amh=amh0
       mhpole = amh0
       endif        
       if(amch.eq.0.d0) then
       amch=amch0
       mchpole = amch0
       endif
       amlight=aml
       amheavy=amh

       aml = aml0               
       amh = amh0               
       amch = amch0
       ama0 = ama

cccccccccccccccccccccccc
       CALL SU_HLOOP(q2,amlight,amu,Au,Ad,Al,
     .     pis1s1l,pis1s2l,pis2s2l,piaa,picc,pizz,piww,tad1,tad2) 
       CALL SU_HLOOP(q2,amheavy,amu,Au,Ad,Al,
     .     pis1s1h,pis1s2h,pis2s2h,piaa,picc,pizz,piww,tad1,tad2) 
c
      vd2 = 2*(amz**2+pizz)/(g1ewsb**2+g2ewsb**2)/(1.d0+tbeta**2)
      vu2 = vd2*tbeta**2      
      vev2=2.d0*(vu2+vd2)
       rmlDR = ytau*dsqrt(vd2)
       rmbDR = yb*dsqrt(vd2)
       rmtDR = yt*dsqrt(vu2)
       gstrong=dsqrt(4.d0*pi*alsewsb)
       sxt=dsin(thett)
       cxt=dcos(thett)
       sxb=dsin(thetb)
       cxb=dcos(thetb)      
       cxl=dcos(thetl)
       sxl=dsin(thetl)

       pizzp = pizz              
c 
      ihdr=0.d0
c%%%%%%%%%%%%%%%%%%%%%%%%%%%   two--loop alphas corrections (P. Slavich)
 1     if(imodel.ge.2) then
      call SU_DSZHiggs(rmtDR**2,am3,mst1**2,mst2**2,sxt,cxt,scale**2,
     .      -amu,tbeta,vev2,gstrong,0,S11s,S22s,S12s)
c
      call SU_DSZHiggs(rmbDR**2,am3,msb1**2,msb2**2,sxb,cxb,scale**2,
     .     -amu,1.d0/tbeta,vev2,gstrong,0,S22b,S11b,S12b)
           
      call SU_DSZodd(rmtDR**2,am3,mst1**2,mst2**2,sxt,cxt,scale**2,
     .     -amu,tbeta,vev2,gstrong,P2s)
c
      call SU_DSZodd(rmbDR**2,am3,msb1**2,msb2**2,sxb,cxb,scale**2,
     .     -amu,1.d0/tbeta,vev2,gstrong,P2b)
c
c%%%%%%%%  two-loop electroweak corrections (P. Slavich routines)
c       
       
      call SU_DDSHiggs(rmtDR**2,rmbDR**2,amar**2,mst1**2,mst2**2,
     .     msb1**2,msb2**2,sxt,cxt,sxb,cxb,scale**2,-amu,tbeta,vev2,
     .     S11w,S12w,S22w)
c
      call SU_DDSodd(rmtDR**2,rmbDR**2,amar**2,mst1**2,mst2**2,
     .    msb1**2,msb2**2,sxt,cxt,sxb,cxb,scale**2,-amu,tbeta,vev2,P2w)
        
c%%%%%%%%%%%%%%%%%%%% Now add the tau-lepton contributions. 
      call SU_taubot(rmlDR**2,rmbDR**2,msta1**2,msta2**2,msb1**2,
     .               msb2**2,sxl,cxl,sxb,cxb,scale**2,-amu,tbeta,vev2,
     .               S11bl,S12bl,S22bl)
      
      call SU_taubotodd(rmlDR**2,rmbDR**2,msta1**2,msta2**2,msb1**2,
     .          msb2**2,sxl,cxl,sxb,cxb,scale**2,-amu,tbeta,vev2, P2bl) 

      call SU_tausqHiggs(rmlDR**2,amar**2,msntau**2,msta1**2,msta2**2,
     .       sxl,cxl,scale**2,-amu,tbeta,vev2,0,
     $        S11l,S22l,S12l)

      call SU_tausqodd(rmlDR**2,amar**2,msntau**2,msta1**2,msta2**2,
     .       sxl,cxl,scale**2,-amu,tbeta,vev2,P2l)

c%%%%%%%%%%% 2-loop tadpole corrections (P. slavich routines)

      call SU_ewsb2loop(rmtDR**2,am3,mst1**2,mst2**2,sxt,cxt,
     .     scale**2,-amu,tbeta,vev2,gstrong,tad1st,tad2st)
     
      call SU_ewsb2loop(rmbDR**2,am3,msb1**2,msb2**2,sxb,cxb,
     .     scale**2,-amu,1.d0/tbeta,vev2,gstrong,tad2sb,tad1sb)

      call SU_DDStad(rmtDR**2,rmbDR**2,amar**2,mst1**2,mst2**2,
     .     msb1**2,msb2**2,sxt,cxt,sxb,cxb,scale**2,-amu,tbeta,vev2,
     .     tad1w,tad2w)

      call SU_taubottad(rmlDR**2,rmbDR**2,msta1**2,msta2**2,msb1**2,
     .      msb2**2,sxl,cxl,sxb,cxb,scale**2,-amu,tbeta,vev2,
     $      tad1bl,tad2bl) 

      call SU_tausqtad(rmlDR**2,amar**2,msntau**2,msta1**2,msta2**2,
     .      sxl,cxl,scale**2,-amu,tbeta,vev2,
     $      tad1l,tad2l)

      else
c full one-loop Higgs calculation but neglecting any two-loops
      S11s = 0.d0
      S12s= 0.d0
      S22s=0.d0
      S11b = 0.d0
      S12b= 0.d0
      S22b=0.d0
      P2s=0.d0
      P2b=0.d0
      S11w = 0.d0
      S12w= 0.d0
      S22w=0.d0
      P2w =0.d0
      S11bl = 0.d0
      S12bl= 0.d0
      S22bl=0.d0
      P2bl=0.d0
      S11l = 0.d0
      S12l= 0.d0
      S22l=0.d0
      P2l=0.d0
c
      endif
c     add two-loop tadpoles in running mA:
       if(imodel.ge.2) tad1loop= tad1st+tad1sb+tad1w+tad1l+tad1bl   	
      dVdvd2=-tad1
      if(imodel.ge.2) dVdvd2=dVdvd2+tad1loop

       if(imodel.ge.2) tad2loop=tad2st+tad2sb+tad2w+tad2l+tad2bl	
      dVdvu2=-tad2
      if(imodel.ge.2) dVdvu2=dVdvu2+tad2loop

      mz2dr=amz**2+pizz
      madr2=(mhu2+dVdvu2 -mhd2-dVdvd2)/dcos(2*bet)-mz2dr 

      DMA=P2s+P2w+P2b+P2l+P2bl
      if(kmaflag.eq.0) then    !! then ama is really MA pole input
      ama=marsave   
        madr2 =ama**2 +piaa -sb**2*tad1-cb**2*tad2 -DMA
      madr2save=madr2    
      endif
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      mL11loop = mz2dr*cb**2+madr2*sb**2-pis1s1l+tad1+
     .          S11s+S11w+S11b+S11l+S11bl+DMA*sb**2
      mL22loop = mz2dr*sb**2+madr2*cb**2-pis2s2l+tad2+
     .          S22s+S22w+S22b+S22l+S22bl+DMA*cb**2
      mL12loop = -(mz2dr+madr2)*sb*cb-pis1s2l+
     .           S12s+S12w+S12b+S12l+S12bl-DMA*sb*cb

      mH11loop = mz2dr*cb**2+madr2*sb**2-pis1s1h+tad1+
     .           S11s+S11w+S11b+S11l+S11bl+DMA*sb**2
      mH22loop = mz2dr*sb**2+madr2*cb**2-pis2s2h+tad2+
     .           S22s+S22w+S22b+S22l+S22bl+DMA*cb**2
      mH12loop = -(mz2dr+madr2)*sb*cb-pis1s2h+
     .           S12s+S12w+S12b+S12l+S12bl-DMA*sb*cb

      mLcr2=0.5d0*(mL11loop+mL22loop-dsqrt((mL11loop-mL22loop)**2
     .            +4*mL12loop**2) )
      mHcr2=0.5d0*(mH11loop+mH22loop+dsqrt((mH11loop-mH22loop)**2
     .            +4*mH12loop**2) )
c
      if(mLcr2.ge.0.d0) then
         mLpole=dsqrt(mLcr2)
      else
         mLpole=aml0
      if(irge.eq.irgmax.and.ifix.ge.3) mlpole=dsqrt(mlcr2)
      endif

      if(mHcr2.ge.0.d0) then
         mHpole=dsqrt(mHcr2)
      else
         mHpole=amh0
      if(irge.eq.irgmax.and.ifix.ge.3) mhpole=dsqrt(mhcr2)
      endif

      mH2dum=0.5d0*(mL11loop+mL22loop+dsqrt((mL11loop-mL22loop)**2
     .     +4*mL12loop**2) )    

      s2alfa=2.d0*mL12loop/(mH2dum-mLpole**2) 
      c2alfa= (mL11loop-mL22loop)/(mH2dum-mLpole**2)
       t2alfa=s2alfa/c2alfa   
c next is to have correct alpha angle convention:
       if(c2alfa.gt.0.d0) then
       a=0.5d0*datan(t2alfa)
       endif
       if(c2alfa.lt.0.d0) then
       if(s2alfa.lt.0.d0) then  
       a=0.5d0*datan(t2alfa)-pi/2
       else
       a=0.5d0*datan(t2alfa)+pi/2
       endif
       endif
c      
      tadba=sb**2*tad1+cb**2*tad2
      mAcr2 =madr2-piaa+tadba+DMA
      mCHcr2=macr2+amw**2+piaa-picc+piww      
      if(macr2.ge.0.d0) then
      mApole = dsqrt(mAcr2)
      else
      mApole = amA
      if(irge.eq.irgmax.and.ifix.ge.3) mapole=dsqrt(macr2)
      endif
      if(mCHcr2.ge.0.d0) then 
      mCHpole = dsqrt(mCHcr2)      
      else
      mCHpole = amch0
      if(irge.eq.irgmax.and.ifix.ge.3) mchpole=dsqrt(mchcr2)
      endif
c=========  end of the full calculation
      endif
c
      la3t=la3+la4+la5
      ama2=amar**2
      aml2=amlr**2
      amh2=amhr**2
      amp2=amchr**2
c ========== higgs couplings 
      sbma = sb*ca-cb*sa
      cbma = cb*ca+sb*sa
      sbpa = sb*ca+cb*sa
      cbpa = cb*ca-sb*sa
      s2a = 2*sa*ca
      c2a = ca**2-sa**2
      s2b = 2*sb*cb
      c2b = cb**2-sb**2
      glzz = 1.d0/v/2*aml2*sbma
      ghzz = 1.d0/v/2*amh2*cbma
      glww = 2*glzz
      ghww = 2*ghzz
      glaz = 1.d0/v*(aml2-ama2)*cbma
      ghaz = -1.d0/v*(amh2-ama2)*sbma
      glpw = -1.d0/v*(amp2-aml2)*cbma
      glmw = glpw
      ghpw = 1.d0/v*(amp2-amh2)*sbma
      ghmw = ghpw
      gapw = 1.d0/v*(amp2-ama2)
      gamw = -gapw
      ghhh = v/2*(la1*ca**3*cb + la2*sa**3*sb + la3t*sa*ca*sbpa
     .     + la6*ca**2*(3*sa*cb+ca*sb) + la7*sa**2*(3*ca*sb+sa*cb))
      glll = -v/2*(la1*sa**3*cb - la2*ca**3*sb + la3t*sa*ca*cbpa
     .     - la6*sa**2*(3*ca*cb-sa*sb) + la7*ca**2*(3*sa*sb-ca*cb))
      glhh = -3*v/2*(la1*ca**2*cb*sa - la2*sa**2*sb*ca
     .     + la3t*(sa**3*cb-ca**3*sb+2*sbma/3)
     .     - la6*ca*(cb*c2a-sa*sbpa) - la7*sa*(c2a*sb+ca*sbpa))
      ghll = 3*v/2*(la1*sa**2*cb*ca + la2*ca**2*sb*sa
     .     + la3t*(sa**3*sb+ca**3*cb-2*cbma/3)
     .     - la6*sa*(cb*c2a+ca*cbpa) + la7*ca*(c2a*sb+sa*cbpa))
      glaa = -v/2*(la1*sb**2*cb*sa - la2*cb**2*sb*ca
     .     - la3t*(sb**3*ca-cb**3*sa) + 2*la5*sbma
     .     - la6*sb*(cb*sbpa+sa*c2b) - la7*cb*(c2b*ca-sb*sbpa))
      ghaa = v/2*(la1*sb**2*cb*ca + la2*cb**2*sb*sa
     .     + la3t*(sb**3*sa+cb**3*ca) - 2*la5*cbma
     .     - la6*sb*(cb*cbpa+ca*c2b) + la7*cb*(sb*cbpa+sa*c2b))
      glpm = 2*glaa + v*(la5 - la4)*sbma
      ghpm = 2*ghaa + v*(la5 - la4)*cbma
      glzz = 2*glzz
      ghzz = 2*ghzz
      glll = 6*glll
      ghhh = 6*ghhh
      glhh = 2*glhh
      ghll = 2*ghll
      glaa = 2*glaa
      ghaa = 2*ghaa
      xnorm = amz**2/v
      glll = glll/xnorm
      ghll = ghll/xnorm
      glhh = glhh/xnorm
      ghhh = ghhh/xnorm
      ghaa = ghaa/xnorm
      glaa = glaa/xnorm
      glpm = glpm/xnorm
      ghpm = ghpm/xnorm
      gat=1.d0/tgbet
      gab=tgbet
      glt=ca/sb
      glb=-sa/cb
      ght=sa/sb
      ghb=ca/cb
      gzal=-cbma
      gzah=sbma
      glvv=sbma
      ghvv=cbma
      b=bet
c
C   Higgs couplings needed in SUSPECT:
      alfa = a
      xgat = gat
      xgab = gab
      xglt = glt
      xglb = glb
      xght = ght
      xghb = ghb
      xghvv= ghvv
      xglvv= glvv
C ===============================================================
C ========== Pole Higgs masses (but at one-loop)
c      if(imodel.eq.1.or.imodel.eq.0)then 
c ! affects only pure 1-loop Higgs mass choice:
      if(imodel.eq.0) then 
       xdlt=gf/(2.d0*dsqrt(2.d0)*pi**2)*glt**2*(-2.d0*mt**2+0.5d0*aml2)
     .     *dble(F0_TLH(mt,mt,aml2))
     .     *3*mt**2
       xdlb=gf/(2.d0*dsqrt(2.d0)*pi**2)*glb**2*(-2.d0*mb**2+0.5d0*aml2)
     .     *dble(F0_TLH(mb,mb,aml2))
     .     *3*mb**2
     .     +gf/(2.d0*dsqrt(2.d0)*pi**2)*glb**2*(0.5d0*aml2)
     .     *dlog(mb**2/mt**2)
     .     *3*mb**2
       xdht=gf/(2.d0*dsqrt(2.d0)*pi**2)*ght**2*(-2.d0*mt**2+0.5d0*amh2)
     .     *dble(F0_TLH(mt,mt,amh2))
     .     *3*mt**2
       xdhb=gf/(2.d0*dsqrt(2.d0)*pi**2)*ghb**2*(-2.d0*mb**2+0.5d0*amh2)
     .     *dble(F0_TLH(mb,mb,amh2))
     .     *3*mb**2
     .     +gf/(2.d0*dsqrt(2.d0)*pi**2)*ghb**2*(0.5d0*amh2)
     .     *dlog(mb**2/mt**2)
     .     *3*mb**2
       xdat=gf/(2.d0*dsqrt(2.d0)*pi**2)*gat**2*(-0.5d0*ama2)
     .     *dble(F0_TLH(mt,mt,ama2))
     .     *3*mt**2
       xdab=gf/(2.d0*dsqrt(2.d0)*pi**2)*gab**2*(-0.5d0*ama2)
     .     *dble(F0_TLH(mb,mb,ama2))
     .     *3*mb**2
     .     +gf/(2.d0*dsqrt(2.d0)*pi**2)*gab**2*(-0.5d0*ama2)
     .     *dlog(mb**2/mt**2)
     .     *3*mb**2
       xdlst=0.d0
       xdlsb=0.d0
       xdhst=0.d0
       xdhsb=0.d0
         do 311 i=1,2
         do 311 j=1,2
       xdlst=xdlst+gf/(2.d0*dsqrt(2.d0)*pi**2)*gltt(i,j)**2*
     .       dble(F0_TLH(mst(i),mst(j),aml2))
     .     *3*amz**4
       xdlsb=xdlsb+gf/(2.d0*dsqrt(2.d0)*pi**2)*glbb(i,j)**2*
     .       dble(F0_TLH(msb(i),msb(j),aml2))
     .    *3*amz**4
       xdhst=xdhst+gf/(2.d0*dsqrt(2.d0)*pi**2)*ghtt(i,j)**2*
     .       dble(F0_TLH(mst(i),mst(j),amh2))
     .     *3*amz**4
       xdhsb=xdhsb+gf/(2.d0*dsqrt(2.d0)*pi**2)*ghbb(i,j)**2*
     .       dble(F0_TLH(msb(i),msb(j),amh2))
     .     *3*amz**4
311    continue
       xdast=gf/(1.d0*dsqrt(2.d0)*pi**2)*gatt**2*
     .       dble(F0_TLH(mst(1),mst(2),ama2))
     .     *3*amz**4
       xdasb=gf/(1.d0*dsqrt(2.d0)*pi**2)*gabb**2*
     .       dble(F0_TLH(msb(1),msb(2),ama2))
     .     *3*amz**4
      
       aml=dsqrt(aml2+xdlt+xdlb+xdlst+xdlsb)
       amh=dsqrt(amh2+xdht+xdhb+xdhst+xdhsb)  
       ama=dsqrt(ama2+xdat+xdab+xdast+xdasb)  
       amch=dsqrt(ama**2+amw**2)     
      else
       aml=mlpole
       amh=mhpole     
       ama=mapole
       amch=mchpole
       alfa=a        
      endif 
      return
111   return
      end
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SU_HLOOP(q2,mhiggs,MU,AT,AB,AL,
     .           pis1s1,pis1s2,pis2s2,piaa,picc,pizz,piww,tad1,tad2) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
c  The main subroutine for the EWSB and calculates the tadpole corrections to
c  the Higgs mass terms squared. The inputs are:
c  q2: the scale at which EWSB is supposed to happen,
c  MU: the higgsino parameter mu ar EWSB scale
c  AT,AB,AL: the third generation trilinear couplings at EWSB scale
c  MQL,MUR,MDR,MEL,MER,MQL1,MUR1,MDR1,MEL1,MER1: the soft parameters at EWSB
c  Other important input parameters, such as the Higgs, chargino, neutralino 
c  masses and couplings as well as SM parameters are called via commons.
c  The output are dVdvd2, dVdvu2, which are (up to some appropriate overall 
c  constants) the derivatives of the full one-loop scalar potential including
c  the contributions of all SM and SUSY particles a la PBMZ (hep-ph/9606211).
c  The consistency of the EWSB mechanism is performed by the main program
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit real*8(a-h,m,o-z)
      real*8 nf
      complex*16 SU_B0,SU_BH,SU_BT22,SU_BG,SU_BF
      logical su_isNaN
      dimension u(2,2),v(2,2),z(4,4),dxmn(4),gmn(4),gmc(2)
      COMMON/SU_cte/nf,cpi,mz,mw,tbetdum
      COMMON/SU_hmass/ma,ml,mh,mch,mar
      COMMON/SU_outhiggs/dml,dmh,dmch,alfa 
      COMMON/SU_outginos/mc1,mc2,mn1,mn2,mn3,mn4,mgluino
      COMMON/SU_matino/u,v,z,dxmn
      COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_tbewsb/vuewsb,vdewsb 
      COMMON/SU_fmasses/mtau,mb,mt
      COMMON/SU_renscale/scale     
      COMMON/SU_bpew/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     .      msta1,msta2,msb1,msb2,mst1,mst2,thet,theb,thel
          COMMON/pietro/mApole,mCHpole
          COMMON/run_p/pizzp     
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
       common/su_nonpert/inonpert

c basic parameters and definitions used:
       sq2=dsqrt(2.d0)
       pi = 4*datan(1.d0)
       scale= dsqrt(q2)
       g=g2ewsb
       g1=g1ewsb
c defining s^2_w at EWSB scale:
       cw = 1.d0/dsqrt(1.d0+(g1ewsb/g2ewsb)**2)
       sw = g1ewsb/g2ewsb*cw
       sw2=sw**2
       cw2= cw**2
       cwm2 =1.d0/cw2
c defining other parameters for the Higgs mass calculation
       alph= (g*sw)**2/(4*pi)
       tbeta = vuewsb/vdewsb
       B=datan(tbeta)
       beta=B
       mup=1.d-2
       mdo=1.d-2
       me=0.5d-3
       mmu=0.106d0
       ms=0.190d0
       mcq=1.40d0
       eps=1.d-5
       eps0=eps**2
c   
       gmn(1)=dabs(dxmn(1))
       gmn(2)=dabs(dxmn(2))
       gmn(3)=dabs(dxmn(3))
       gmn(4)=dabs(dxmn(4))
       gmc(1)=mc1
       gmc(2)=mc2
c
       ct=dcos(thet)
       st=dsin(thet)
       cb=dcos(theb)
       sb=dsin(theb)
       cta=dcos(thel)
       sta=dsin(thel)
c
       cbeta2=1.d0/(1.d0+tbeta**2)
       cbet= dsqrt(cbeta2)
       sbet=dsqrt(1.d0-cbeta2)
       c2b =2*cbeta2-1.d0
c        
       if(su_isNaN(pizzp).or.mz**2+pizzp.le.0.d0) then  !!added protection
       pizzp=0.d0
       if(irge.eq.irgmax)  inonpert=-1
       endif
c
       alfasave = alfa          !  (alfa running)
       alfa =0.5*atan(tan(2d0*b)*(mar**2+mz**2+pizzp)
     .      /(mar**2-mz**2-pizzp))
c  ! take into account correct alpha sign convention:
       if(cos(2d0*b)*(mar**2-mz**2-pizzp).gt.0.d0) alfa = alfa-pi/2
       sal=dsin(alfa)
       cal=dcos(alfa)
       s2a = 2*sal*cal 
       s2al=s2a
c     compute running masses 
       vd2 = 2*(mz**2+pizzp)/(g1ewsb**2+g2ewsb**2)/(1.d0+tbeta**2)
       rmz= dsqrt(mz**2+pizzp)
       vu2 = vd2*tbeta**2
       vd= dsqrt(vd2)
       vu= dsqrt(vu2)
       rmt=yt*vu
       rmb=yb*vd
       rmtau=ytau*vd
       rmw=rmz*cw
c     use running masses everywhere
       mzsave = mz
       mwsave = mw
       mtsave = mt
       mbsave = mb 
       mtausave = mtau
       mz = rmz
       mw = rmw
       mt = rmt
       mb = rmb
       mtau = rmtau

c-----------------------------------------------------------------
c                 Z boson self-energy at q**2=mz**2 and q**2=0
c-----------------------------------------------------------------
      qsz=mzsave**2           
c
c      
      pizzf = 3*( (.5d0-2*sw2/3)**2+(2*sw2/3)**2)
     .*(dble(SU_BH(qsz,mt,mt))+dble(SU_BH(qsz,mcq,mcq))
     .       +dble(SU_BH(qsz,mup,mup)))
     .  + 3*((-.5d0+sw2/3)**2+(-sw2/3)**2)
     .*(dble(SU_BH(qsz,mb,mb))+dble(SU_BH(qsz,ms,ms))
     .       +dble(SU_BH(qsz,mdo,mdo)))
     .  + ((-.5d0+sw2)**2+(-sw2)**2)*(dble(SU_BH(qsz,me,me))
     .    +dble(SU_BH(qsz,mmu,mmu))+dble(SU_BH(qsz,mtau,mtau)))
     .  + .5d0**2*3*dble(SU_BH(qsz,eps0,eps0))
     .  -12*(.5d0-2*sw2/3)*(2*sw2/3)
     .  *(mt**2*dble(SU_B0(qsz,mt,mt))+mcq**2*dble(SU_B0(qsz,mcq,mcq))
     .  +mup**2*dble(SU_B0(qsz,mup,mup))) 
     .  -12*(-.5d0+sw2/3)*(-sw2/3)
     .  *(mb**2*dble(SU_B0(qsz,mb,mb))+ms**2*dble(SU_B0(qsz,ms,ms))
     .  +mdo**2*dble(SU_B0(qsz,mdo,mdo))) 
     .  -4*(-.5d0+sw2)*(-sw2)*(me**2*dble(SU_B0(qsz,me,me))+mmu**2
     .  *dble(SU_B0(qsz,mmu,mmu))+mtau**2*dble(SU_B0(qsz,mtau,mtau)))
c     
      pizzb = -2*cw**4*(2*qsz+mw**2-mz**2*sw**4/cw**2)
     . *dble(SU_B0(qsz,mw,mw))
     . -(8*cw**4+(cw2-sw2)**2)*dble(SU_BT22(qsz,mw,mw))
c
      pizzh0=-dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml))
      pizzhS = -dsin(beta-alfa)**2*(dble(SU_BT22(qsz,ma,mh))
     .  + dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsz,mz,mh))
     .  + dble(SU_BT22(qsz,ma,ml))-mz**2*dble(SU_B0(qsz,mz,mh)) )
     .  -(cw**2-sw**2)**2*dble(SU_BT22(qsz,mch,mch))
     .  -pizzh0
c
       pizzsu= -12*( (.5d0-2*sw2/3)*dcos(thet)**2
     .-(2*sw2/3)*dsin(thet)**2 )**2*dble(SU_BT22(qsz,mst1,mst1))
     .         -12*(-(.5d0-2*sw2/3)*dsin(thet)**2
     .+(2*sw2/3)*dcos(thet)**2 )**2*dble(SU_BT22(qsz,mst2,mst2))
     .      -24*( (.5d0)*dsin(thet)*dcos(thet) )**2
     .  *dble(SU_BT22(qsz,mst1,mst2))
     .    -24*(.5d0-2*sw2/3)**2*dble(SU_BT22(qsz,msu1,msu1))
     .    -24*(+2*sw2/3)**2*dble(SU_BT22(qsz,msu2,msu2))
c
      pizzsd= -12*( (-.5d0+sw2/3)*dcos(theb)**2
     .-(-sw2/3)*dsin(theb)**2)**2*dble(SU_BT22(qsz,msb1,msb1))
     .       -12*( -(-.5d0+sw2/3)*dsin(theb)**2
     .+(-sw2/3)*dcos(theb)**2)**2*dble(SU_BT22(qsz,msb2,msb2))
     .      -24*((-0.5d0)*dsin(theb)*dcos(theb))**2
     .  *dble(SU_BT22(qsz,msb1,msb2))
     .    -24*(-.5d0+sw2/3)**2*dble(SU_BT22(qsz,msd1,msd1))
     .    -24*(-sw2/3)**2*dble(SU_BT22(qsz,msd2,msd2))
c
      pizzsl=-4*( (-.5d0+sw2)*dcos(thel)**2
     .- (-sw2)*dsin(thel)**2 )**2*dble(SU_BT22(qsz,msta1,msta1))
     .       -4*( -(-.5d0+sw2)*dsin(thel)**2
     .  +(-sw2)*dcos(thel)**2 )**2*dble(SU_BT22(qsz,msta2,msta2))
     .      -8*((-.5d0)*dsin(thel)*dcos(thel))**2
     .  *dble(SU_BT22(qsz,msta1,msta2))
     .      -8*(-.5d0+sw2)**2*dble(SU_BT22(qsz,mse1,mse1))
     .       -8*(-sw2)**2*dble(SU_BT22(qsz,mse2,mse2))
     .       -8*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
     .       -4*(.5d0)**2*dble(SU_BT22(qsz,msntau,msntau))
c
      pizzs=pizzsl+pizzsd+pizzsu
c
      pizzn=0.d0
      do  i=1,4
      do  j=1,4
      pizzn = pizzn + 1.d0/4*(Z(i,3)*Z(j,3) -Z(i,4)*Z(j,4))**2*
     .       (dble(SU_BH(qsz,gmn(i),gmn(j)))
     .       -2*dxmn(i)*dxmn(j)*dble(SU_B0(qsz,gmn(i),gmn(j))) )
      enddo
      enddo
c
      pizzc=0.d0
      do i=1,2
      do j=1,2
      pizzc = pizzc +1.d0/4*( 
     .( ( 2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2) )**2+
     .  ( 2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2) )**2 )
     .            *dble(SU_BH(qsz,gmc(i),gmc(j))) 
     .     +4*(2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2))*
     .        (2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2))*
     .            gmc(i)*gmc(j)*dble(SU_B0(qsz,gmc(i),gmc(j))) )
      enddo
      enddo
c
c Sum of the susy contributions for pizz and final pizz(MZ**2) 
      pizzsm=alph/4.d0/pi/sw2/cw2*(pizzf+pizzb+pizzh0)
      pizzsusy=alph/4.d0/pi/sw2/cw2*
     .        (pizzhS+pizzs+pizzn+pizzc)
      pizz=pizzsm+pizzsusy
c----------------------------------------------------------------
      qsz=eps 
c
      pizzf0 = 3*( (.5d0-2*sw2/3)**2+(2*sw2/3)**2)
     .*(dble(SU_BH(qsz,mt,mt))+dble(SU_BH(qsz,mcq,mcq))
     .  +dble(SU_BH(qsz,mup,mup)))
     .  + 3*((-.5d0+sw2/3)**2+(-sw2/3)**2)
     .*(dble(SU_BH(qsz,mb,mb))+dble(SU_BH(qsz,ms,ms))
     .  +dble(SU_BH(qsz,mdo,mdo)))
     .  + ((-.5d0+sw2)**2+(-sw2)**2)*(dble(SU_BH(qsz,me,me))
     .    +dble(SU_BH(qsz,mmu,mmu))+dble(SU_BH(qsz,mtau,mtau)))
     .  + .5d0**2*3*dble(SU_BH(qsz,eps,eps))
     .  -12*(.5d0-2*sw2/3)*(2*sw2/3)
     .  *(mt**2*dble(SU_B0(qsz,mt,mt))+mcq**2*dble(SU_B0(qsz,mcq,mcq))
     .  +mup**2*dble(SU_B0(qsz,mup,mup))) 
     .  -12*(-.5d0+sw2/3)*(-sw2/3)
     .  *(mb**2*dble(SU_B0(qsz,mb,mb))+ms**2*dble(SU_B0(qsz,ms,ms))
     .  +mdo**2*dble(SU_B0(qsz,mdo,mdo))) 
     .  -4*(-.5d0+sw2)*(-sw2)*(me**2*dble(SU_B0(qsz,me,me))+mmu**2
     .  *dble(SU_B0(qsz,mmu,mmu))+mtau**2*dble(SU_B0(qsz,mtau,mtau)))
c
      pizzb0 = -2*cw**4*(2*qsz+mw**2-mz**2*sw**4/cw**2)
     . *dble(SU_B0(qsz,mw,mw))
     . -(8*cw**4+(cw2-sw2)**2)*dble(SU_BT22(qsz,mw,mw))
c
      pizzh00=-dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml))
      pizzhS0 = -dsin(beta-alfa)**2*(dble(SU_BT22(qsz,ma,mh))
     .  + dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsz,mz,mh))
     .  + dble(SU_BT22(qsz,ma,ml))-mz**2*dble(SU_B0(qsz,mz,mh)) )
     .  -(cw**2-sw**2)**2*dble(SU_BT22(qsz,mch,mch))
     .  -pizzh00
c
       pizzsu0= -12*( (.5d0-2*sw2/3)*dcos(thet)**2
     .-(2*sw2/3)*dsin(thet)**2 )**2*dble(SU_BT22(qsz,mst1,mst1))
     .         -12*(-(.5d0-2*sw2/3)*dsin(thet)**2
     .+(2*sw2/3)*dcos(thet)**2 )**2*dble(SU_BT22(qsz,mst2,mst2))
     .      -24*( (.5d0)*dsin(thet)*dcos(thet) )**2
     .  *dble(SU_BT22(qsz,mst1,mst2))
     .    -24*(.5d0-2*sw2/3)**2*dble(SU_BT22(qsz,msu1,msu1))
     .    -24*(+2*sw2/3)**2*dble(SU_BT22(qsz,msu2,msu2))
c
      pizzsd0= -12*( (-.5d0+sw2/3)*dcos(theb)**2
     .-(-sw2/3)*dsin(theb)**2)**2*dble(SU_BT22(qsz,msb1,msb1))
     .       -12*( -(-.5d0+sw2/3)*dsin(theb)**2
     .+(-sw2/3)*dcos(theb)**2)**2*dble(SU_BT22(qsz,msb2,msb2))
     .      -24*((-0.5d0)*dsin(theb)*dcos(theb))**2
     .  *dble(SU_BT22(qsz,msb1,msb2))
     .    -24*(-.5d0+sw2/3)**2*dble(SU_BT22(qsz,msd1,msd1))
     .    -24*(-sw2/3)**2*dble(SU_BT22(qsz,msd2,msd2))
c
      pizzsl0=-4*( (-.5d0+sw2)*dcos(thel)**2
     .- (-sw2)*dsin(thel)**2 )**2*dble(SU_BT22(qsz,msta1,msta1))
     .       -4*( -(-.5d0+sw2)*dsin(thel)**2
     .  +(-sw2)*dcos(thel)**2 )**2*dble(SU_BT22(qsz,msta2,msta2))
     .      -8*((-.5d0)*dsin(thel)*dcos(thel))**2
     .  *dble(SU_BT22(qsz,msta1,msta2))
     .      -8*(-.5d0+sw2)**2*dble(SU_BT22(qsz,mse1,mse1))
     .       -8*(-sw2)**2*dble(SU_BT22(qsz,mse2,mse2))
     .       -8*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
     .       -4*(.5d0)**2*dble(SU_BT22(qsz,msntau,msntau))
c
      pizzs0=pizzsl0+pizzsd0+pizzsu0
c
      pizzn0=0.d0
      do  i=1,4
      do  j=1,4
      pizzn0 = pizzn0 + 1.d0/4*(Z(i,3)*Z(j,3) -Z(i,4)*Z(j,4))**2*
     .       (dble(SU_BH(qsz,gmn(i),gmn(j)))
     .       -2*dxmn(i)*dxmn(j)*dble(SU_B0(qsz,gmn(i),gmn(j))) )
      enddo
      enddo
c
      pizzc0=0.d0
      do i=1,2
      do j=1,2
      pizzc0 = pizzc0 +1.d0/4*( 
     .( ( 2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2) )**2+
     .  ( 2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2) )**2 )
     .            *dble(SU_BH(qsz,gmc(i),gmc(j))) 
     .     +4*(2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2))*
     .        (2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2))*
     .            gmc(i)*gmc(j)*dble(SU_B0(qsz,gmc(i),gmc(j))) )
      enddo
      enddo
c
c Sum of the susy contributions for pizz and final pizz(MZ**2) 
      pizzsm0=alph/4.d0/pi/sw2/cw2*(pizzf0+pizzb0+pizzh00)
      pizzsusy0=alph/4.d0/pi/sw2/cw2*
     .        (pizzhS0+pizzs0+pizzn0+pizzc0)
      pizz0=pizzsm0+pizzsusy0
c
c-----------------------------------------------------------------
c                W boson self-energy at q**2=mw**2 and q**2=0
c-----------------------------------------------------------------
       qsw=mwsave**2 
c
      piwwf=3.d0/2*(dble(SU_BH(qsw,mt,mb))+dble(SU_BH(qsw,mcq,ms))
     . +dble(SU_BH(qsw,mup,mdo)))+0.5d0*(dble(SU_BH(qsw,me,eps))
     . +dble(SU_BH(qsw,mmu,eps))+dble(SU_BH(qsw,mtau,eps)))
c
      piwwb=-(1.d0+8*cw**2)*dble(SU_BT22(qsw,mz,mw))-sw**2*(
     . 8*dble(SU_BT22(qsw,mw,eps))+4*qsw*dble(SU_B0(qsw,mw,eps)))
     . -((4*qsw+mz**2+mw**2)*cw**2-mz**2*sw**4)
     . *dble(SU_B0(qsw,mz,mw))
c
      piwwh0=-   dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw))         
      piwwhS = -dsin(beta-alfa)**2*(dble(SU_BT22(qsw,mh,mch))
     .  + dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsw,ml,mch))
     .  + dble(SU_BT22(qsw,mh,mw))-mw**2*dble(SU_B0(qsw,mh,mw)) )
     .  -dble(SU_BT22(qsw,ma,mch)) -piwwh0
c
      piwws =-2*3*( 2*dble(SU_BT22(qsw,msu1,msd1)) 
     .+dcos(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst1,msb1))
     .+dcos(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst1,msb2))
     .+dsin(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst2,msb1))
     .+dsin(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst2,msb2)) )
     .       -2*(  2*dble(SU_BT22(qsw,msn1,mse1)) 
     . + dcos(thel)**2*dble(SU_BT22(qsw,msntau,msta1))
     . + dsin(thel)**2*dble(SU_BT22(qsw,msntau,msta2)) )
c 
      piwwnc=0.d0
       do i=1,4
       do j=1,2
       piwwnc= piwwnc +
     . ( (-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)**2+
     .   (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)**2 )*
     . dble(SU_BH(qsw,gmn(i),gmc(j))) 
     . + 4*(-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)*
     .        (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)*
     . dxmn(i)*gmc(j)*dble(SU_B0(qsw,gmn(i),gmc(j)))
       enddo
       enddo   
c
c Sum of the susy contributions for piww and final piww(Mw**2)  
      piwwsm=alph/4.d0/pi/sw2*(piwwf+piwwb+piwwh0)
      piwwsusy=alph/4.d0/pi/sw2*(piwwhS+piwws+piwwnc)
      piww=piwwsm+piwwsusy       
c
c-----------------------------------------------------------------
      qsw=eps
c
      piwwf0=3.d0/2*(dble(SU_BH(qsw,mt,mb))+dble(SU_BH(qsw,mcq,ms))
     . +dble(SU_BH(qsw,mup,mdo)))+0.5d0*(dble(SU_BH(qsw,me,eps))
     . +dble(SU_BH(qsw,mmu,eps))+dble(SU_BH(qsw,mtau,eps)))
c
      piwwb0=-(1.d0+8*cw**2)*dble(SU_BT22(qsw,mz,mw))-sw**2*(
     . 8*dble(SU_BT22(qsw,mw,eps))+4*qsw*dble(SU_B0(qsw,mw,eps)))
     . -((4*qsw+mz**2+mw**2)*cw**2-mz**2*sw**4)
     . *dble(SU_B0(qsw,mz,mw))
c
      piwwh00=-   dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw))         
      piwwhS0 = -dsin(beta-alfa)**2*(dble(SU_BT22(qsw,mh,mch))
     .  + dble(SU_BT22(qsw,ml,mw))-mw**2*dble(SU_B0(qsw,ml,mw)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsw,ml,mch))
     .  + dble(SU_BT22(qsw,mh,mw))-mw**2*dble(SU_B0(qsw,mh,mw)) )
     .  -dble(SU_BT22(qsw,ma,mch)) -piwwh00
c
      piwws0 =-2*3*( 2*dble(SU_BT22(qsw,msu1,msd1)) 
     .+dcos(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst1,msb1))
     .+dcos(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst1,msb2))
     .+dsin(thet)**2*dcos(theb)**2*dble(SU_BT22(qsw,mst2,msb1))
     .+dsin(thet)**2*dsin(theb)**2*dble(SU_BT22(qsw,mst2,msb2)) )
     .       -2*(  2*dble(SU_BT22(qsw,msn1,mse1)) 
     . + dcos(thel)**2*dble(SU_BT22(qsw,msntau,msta1))
     . + dsin(thel)**2*dble(SU_BT22(qsw,msntau,msta2)) )
c 
      piwwnc0=0.d0
       do i=1,4
       do j=1,2
       piwwnc0= piwwnc0 +
     . ( (-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)**2+
     .   (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)**2 )*
     . dble(SU_BH(qsw,gmn(i),gmc(j))) 
     . + 4*(-Z(i,2)*V(j,1)+Z(i,4)*V(j,2)/sq2)*
     .        (-Z(i,2)*U(j,1)-Z(i,3)*U(j,2)/sq2)*
     . dxmn(i)*gmc(j)*dble(SU_B0(qsw,gmn(i),gmc(j)))
       enddo
       enddo   
c
c Sum of the susy contributions for piww and final piww(0)  
      piwwsm0=alph/4.d0/pi/sw2*(piwwf0+piwwb0+piwwh00)
      piwwsusy0=alph/4.d0/pi/sw2*(piwwhS0+piwws0+piwwnc0)
      piww0=piwwsm0+piwwsusy0
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%% Now define vu,vd and running masses
c defining mtau,mb,mt running masses at ewsb scale:

      mz = mzsave             
      mw = mwsave              

       vd2 = 2*(mz**2+pizz)/(g1ewsb**2+g2ewsb**2)/(1.d0+tbeta**2)
       vu2 = vd2*tbeta**2
       vd= dsqrt(vd2)
       vu= dsqrt(vu2)
       rmt=yt*vu
       rmb=yb*vd
       rmtau=ytau*vd

       rmz = dsqrt(mz**2+pizz)  ! (USE RUNNING MW,MZ)
       rmw = rmz*cw
       mzsave = mz
       mwsave = mw
       mz = rmz
       mw = rmw

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       s1bLbL = g*mz/cw*(-.5d0 +sw2/3)*cbet +sq2*yb*rmb
       s1bRbR = g*mz/cw*(-sw2/3)*cbet +sq2*yb*rmb
       s1bLbR = yb/sq2*Ab
       s2bLbL = -g*mz/cw*(-.5d0 +sw2/3)*sbet
       s2bRbR = -g*mz/cw*(-sw2/3)*sbet
       s2bLbR = -yb/sq2*mu
c
       gs1d1d1 = g*mz/cw*(-.5d0 +sw2/3)*cbet 
       gs1d2d2 = g*mz/cw*(-sw2/3)*cbet 
       gs2d1d1 = -g*mz/cw*(-.5d0 +sw2/3)*sbet
       gs2d2d2 = -g*mz/cw*(-sw2/3)*sbet
c  
       s1tauLL = g*mz/cw*(-.5d0 +sw2)*cbet +sq2*ytau*rmtau
       s1tauRR = g*mz/cw*(-sw2)*cbet +sq2*ytau*rmtau
       s1tauLR = ytau/sq2*Al
       s2tauLL = -g*mz/cw*(-.5d0 +sw2)*sbet
       s2tauRR = -g*mz/cw*(-sw2)*sbet
       s2tauLR = -ytau/sq2*mu
c
       gs1e1e1 = g*mz/cw*(-.5d0 +sw2)*cbet 
       gs1e2e2 = g*mz/cw*(-sw2)*cbet 
       gs2e1e1 = -g*mz/cw*(-.5d0 +sw2)*sbet
       gs2e2e2 = -g*mz/cw*(-sw2)*sbet
c
       s2tLtL = -g*mz/cw*(.5d0 -2*sw2/3)*sbet +sq2*yt*rmt
       s2tRtR = -g*mz/cw*(2*sw2/3)*sbet +sq2*yt*rmt
       s2tLtR = yt/sq2*At
       s1tLtL = g*mz/cw*(.5d0 -2*sw2/3)*cbet
       s1tRtR = g*mz/cw*(2*sw2/3)*cbet
       s1tLtR = -yt/sq2*mu
c
       gs2u1u1 = -g*mz/cw*(.5d0 -2*sw2/3)*sbet 
       gs2u2u2 = -g*mz/cw*(2*sw2/3)*sbet
       gs1u1u1 = g*mz/cw*(.5d0 -2*sw2/3)*cbet
       gs1u2u2 = g*mz/cw*(2*sw2/3)*cbet
c
       gs2n1n1 = -g*mz/cw*(.5d0 )*sbet 
       gs1n1n1 = g*mz/cw*(.5d0 )*cbet
c
       gs1b1b1 = cb**2*s1bLbL +2*cb*sb*s1bLbR +sb**2*s1bRbR 
       gs1b2b2 = sb**2*s1bLbL -2*cb*sb*s1bLbR +cb**2*s1bRbR 
       gs1b1b2 = sb*cb*(s1bRbR-s1bLbL) +(cb**2-sb**2)*s1bLbR
c
       gs2b1b1 = cb**2*s2bLbL +2*cb*sb*s2bLbR +sb**2*s2bRbR 
       gs2b2b2 = sb**2*s2bLbL -2*cb*sb*s2bLbR +cb**2*s2bRbR 
       gs2b1b2 = sb*cb*(s2bRbR-s2bLbL) +(cb**2-sb**2)*s2bLbR
c
       gs1t1t1 = ct**2*s1tLtL +2*ct*st*s1tLtR +st**2*s1tRtR 
       gs1t2t2 = st**2*s1tLtL -2*ct*st*s1tLtR +ct**2*s1tRtR 
       gs1t1t2 = st*ct*(s1tRtR-s1tLtL) +(ct**2-st**2)*s1tLtR
	
c
       gs2t1t1 = ct**2*s2tLtL +2*ct*st*s2tLtR +st**2*s2tRtR 
       gs2t2t2 = st**2*s2tLtL -2*ct*st*s2tLtR +ct**2*s2tRtR 
       gs2t1t2 = st*ct*(s2tRtR-s2tLtL) +(ct**2-st**2)*s2tLtR
c
       gs1tau11 = cta**2*s1tauLL +2*cta*sta*s1tauLR +sta**2*s1tauRR 
       gs1tau22 = sta**2*s1tauLL -2*cta*sta*s1tauLR +cta**2*s1tauRR 
       gs1tau12 = sta*cta*(s1tauRR-s1tauLL) +(cta**2-sta**2)*s1tauLR
c
       gs2tau11 = cta**2*s2tauLL +2*cta*sta*s2tauLR +sta**2*s2tauRR 
       gs2tau22 = sta**2*s2tauLL -2*cta*sta*s2tauLR +cta**2*s2tauRR 
       gs2tau12 = sta*cta*(s2tauRR-s2tauLL) +(cta**2-sta**2)*s2tauLR
c
       gat1t2=-yt/sq2*(-mu*sbet-At*cbet)
       gab1b2=-yb/sq2*(-mu*cbet-Ab*sbet)
       gatau12=-ytau/sq2*(-mu*cbet-Al*sbet)
c       
       gctLbL = g*mw/sq2*dsin(2*b)-yt*rmt*cbet-yb*rmb*sbet
       gctRbR = -yt*rmb*cbet-yb*rmt*sbet
       gctLbR = yb*(-mu*cbet-Ab*sbet)
       gctRbL = yt*(-mu*sbet-At*cbet)
       gct1b1 = ct*(cb*gctLbL +sb*gctLbR) +st*(cb*gctRbL +sb*gctRbR)
       gct1b2 = ct*(-sb*gctLbL +cb*gctLbR) +st*(-sb*gctRbL+cb*gctRbR)
       gct2b1 = -st*(cb*gctLbL +sb*gctLbR) +ct*(cb*gctRbL +sb*gctRbR)
       gct2b2 = -st*(-sb*gctLbL +cb*gctLbR) +ct*(-sb*gctRbL+cb*gctRbR)
c
       gctauLL = g*mw/sq2*dsin(2*b)-ytau*rmtau*sbet
       gctauLR = ytau*(-mu*cbet-AL*sbet)
       gctau11 = cta*gctauLL +sta*gctauLR
       gctau12 = -sta*gctauLL +cta*gctauLR
       gceLL = g*mw/sq2*dsin(2*b)

c------------------------------------------------------------------
c                              pis1s1
c-------------------------------------------------------------------      
       qs1 = mhiggs**2
c
       pis1s1f=3*yb**2*((qs1-4*rmb**2)*dble(SU_B0(qs1,rmb,rmb))
     .        -2*SU_A(rmb))
     .+ytau**2*((qs1-4*rmtau**2)*dble(SU_B0(qs1,rmtau,rmtau))
     .        -2*SU_A(rmtau))
c
       pis1s1s=3*yb**2*(SU_A(msb1)+SU_A(msb2))+
     .          ytau**2*(SU_A(msta1)+SU_A(msta2)) 
     .           + g**2/(2*cw**2)*( 
     .   3*(.5d0-2*sw2/3)*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) 
     .                     +2*SU_A(msu1) )+
     .   3*(-.5d0 +sw2/3)*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) 
     .                     +2*SU_A(msd1) )+
     .   3*(2*sw2/3)*(ct**2*SU_A(mst2)+st**2*SU_A(mst1)
     .                +2*SU_A(msu2) )+
     .   3*( -sw2/3)*(cb**2*SU_A(msb2)+sb**2*SU_A(msb1) 
     .                +2*SU_A(msd2) ) +
     .   2*(.5d0)*SU_A(msn1) + (0.5d0)*SU_A(msntau) +
     .   (-.5d0 +sw2)*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) 
     .                 +2*SU_A(mse1) )+
     .   (-sw2)*(cta**2*SU_A(msta2)+sta**2*SU_A(msta1) 
     .           +2*SU_A(mse2) )  )  +
     .  3*(gs1t1t1**2*dble(SU_B0(qs1,mst1,mst1))+
     .     gs1t2t2**2*dble(SU_B0(qs1,mst2,mst2)) +
     .    2*gs1t1t2**2*dble(SU_B0(qs1,mst1,mst2)) ) +
     .  3*(gs1b1b1**2*dble(SU_B0(qs1,msb1,msb1))+
     .     gs1b2b2**2*dble(SU_B0(qs1,msb2,msb2)) +
     .    2*gs1b1b2**2*dble(SU_B0(qs1,msb1,msb2)) ) +
     .  3*2*(gs1u1u1**2*dble(SU_B0(qs1,msu1,msu1))+
     .     gs1u2u2**2*dble(SU_B0(qs1,msu2,msu2)) ) +
     .  3*2*(gs1d1d1**2*dble(SU_B0(qs1,msd1,msd1))+
     .     gs1d2d2**2*dble(SU_B0(qs1,msd2,msd2)) ) +
     .    gs1tau11**2*dble(SU_B0(qs1,msta1,msta1))+
     .     gs1tau22**2*dble(SU_B0(qs1,msta2,msta2)) +
     .    2*gs1tau12**2*dble(SU_B0(qs1,msta1,msta2)) +
     .  2*gs1n1n1**2*dble(SU_B0(qs1,msn1,msn1)) +
     .  1*gs1n1n1**2*dble(SU_B0(qs1,msntau,msntau)) +
     .    2*gs1e1e1**2*dble(SU_B0(qs1,mse1,mse1))+
     .    2*gs1e2e2**2*dble(SU_B0(qs1,mse2,mse2)) 
          ms=0.d0

      pis1s1v = g**2/4* (sbet**2*(2*dble(SU_BF(qs1,mch,mw))+
     . dble(SU_BF(qs1,mA,mz))/cw**2 ) +
     . cbet**2*(2*dble(SU_BF(qs1,mw,mw))+dble(SU_BF(qs1,mz,mz))/cw**2 )+
     . 7*cbet**2*(2*mw**2*dble(SU_B0(qs1,mw,mw))+mz**2/cw**2*
     . dble(SU_B0(qs1,mz,mz)) )  +4*(2*SU_A(mw) +SU_A(mz)/cw**2)   )
c
      pis1s1h3 = g**2*mz**2/(4*cw**2)/2*(
     . (cbet*(3*cal**2-sal**2)-sbet*s2al)**2*dble(SU_B0(qs1,MH,MH)) +
     . (-2*cbet*s2al-sbet*dcos(2*alfa))**2*dble(SU_B0(qs1,MH,ml)) +
     . (-2*cbet*s2al-sbet*dcos(2*alfa))**2*dble(SU_B0(qs1,ml,MH)) +
     . (cbet*(3*sal**2-cal**2)+sbet*s2al)**2*dble(SU_B0(qs1,ml,ml)) +
     .  (cbet*dcos(2*b))**2*dble(SU_B0(qs1,mz,mz)) +
     .  (cbet*dsin(2*b))**2*dble(SU_B0(qs1,mz,MA)) +
     .  (cbet*dsin(2*b))**2*dble(SU_B0(qs1,MA,mz)) +
     .  (cbet*dcos(2*b))**2*dble(SU_B0(qs1,MA,MA)) +   
     .  2*( (cbet*dcos(2*b))**2*dble(SU_B0(qs1,mw,mw)) +
     .     (-cbet*dsin(2*b)+cw**2*sbet)**2*dble(SU_B0(qs1,mw,mch))+
     .    (-cbet*dsin(2*b)+cw**2*sbet)**2*dble(SU_B0(qs1,mch,mw)) +
     .    (-cbet*dcos(2*b)+2*cw**2*cbet)**2*dble(SU_B0(qs1,mch,mch)) ) )
c
       pis1s1h4 = g**2/(4*cw**2)/2*( (3*cal**2-sal**2)*SU_A(MH)  + 
     . (3*sal**2-cal**2)*SU_A(ml)  + dcos(2*b)*SU_A(mw) 
     . -dcos(2*b)*SU_A(MA) + 2*( (cw**2 +sw2*dcos(2*b))*SU_A(mw)+
     .  (cw**2 -sw2*dcos(2*b))*SU_A(mch)))
c
       pis1s1b=pis1s1v+pis1s1h3+pis1s1h4
c      
      pis1s1Nino =0.d0
      do i=1,4
      do j=1,4      
    
      pis1s1Nino =  pis1s1Nino + .5d0/4.d0*2.d0*(
     . ( - ( Z(i,1)*Z(j,3) + Z(j,1)*Z(i,3) ) *g1 +
     .     ( Z(i,2)*Z(j,3) + Z(j,2)*Z(i,3) ) *g )**2 *
     .  (dble(SU_BG(qs1,gmn(i),gmn(j))) -      
     .   2*dxmn(i)*dxmn(j)*dble(SU_B0(qs1,gmn(i),gmn(j)))  ) )
     
      enddo
      enddo
c      
      pis1s1Cino =0.d0
      do i=1,2
      do j=1,2      
       pis1s1Cino =     pis1s1Cino +g**2/2 *(
     . ((V(i,1)*U(j,2))**2 +(U(i,2)*V(j,1))**2)*
     . dble(SU_BG(qs1,gmc(i),gmc(j)))
     . -4*V(i,1)*U(j,2)*U(i,2)*V(j,1)*gmc(i)*gmc(j)*
     . dble(SU_B0(qs1,gmc(i),gmc(j)))  )
      enddo
      enddo

c ------ Sum everything: 
      pis1s1= 1.d0/(16*pi**2)*
     . (pis1s1f+pis1s1s+pis1s1b+pis1s1Nino+pis1s1Cino)
c
c-----------------------------------------------------------------
c                           pis2s2 
c-----------------------------------------------------------------
c
      pis2s2f =3*yt**2*((qs1-4*rmt**2)*dble(SU_B0(qs1,rmt,rmt))
     .        -2*SU_A(rmt))

c      
      pis2s2s = 3*yt**2*( SU_A(mst1)+SU_A(mst2) ) - g**2/(2*cw**2)*( 
     .   3*(.5d0-2*sw2/3)*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) 
     .                     +2*SU_A(msu1) )+
     .   3*(-.5d0 +sw2/3)*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) 
     .                     +2*SU_A(msd1) )+
     .   3*(+2*sw2/3)*(ct**2*SU_A(mst2)+st**2*SU_A(mst1) 
     .                +2*SU_A(msu2) )+
     .   3*(-sw2/3)*(cb**2*SU_A(msb2)+sb**2*SU_A(msb1) 
     .               +2*SU_A(msd2) ) +
     .   2*(.5d0)*SU_A(msn1) + 1*(.5d0)*SU_A(msntau)+
     .   (-.5d0 +sw2)*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) 
     .                 +2*SU_A(mse1) )+
     .   (-sw2)*(cta**2*SU_A(msta2)+sta**2*SU_A(msta1) 
     .          +2*SU_A(mse2) )  )  +
     .  3*(gs2t1t1**2*dble(SU_B0(qs1,mst1,mst1))+
     .     gs2t2t2**2*dble(SU_B0(qs1,mst2,mst2)) +
     .    2*gs2t1t2**2*dble(SU_B0(qs1,mst1,mst2)) ) +
     .  3*(gs2b1b1**2*dble(SU_B0(qs1,msb1,msb1))+
     .     gs2b2b2**2*dble(SU_B0(qs1,msb2,msb2)) +
     .    2*gs2b1b2**2*dble(SU_B0(qs1,msb1,msb2)) ) +
     .  3*2*(gs2u1u1**2*dble(SU_B0(qs1,msu1,msu1))+
     .     gs2u2u2**2*dble(SU_B0(qs1,msu2,msu2)) ) +
     .  3*2*(gs2d1d1**2*dble(SU_B0(qs1,msd1,msd1))+
     .     gs2d2d2**2*dble(SU_B0(qs1,msd2,msd2)) ) +
     .    gs2tau11**2*dble(SU_B0(qs1,msta1,msta1))+
     .     gs2tau22**2*dble(SU_B0(qs1,msta2,msta2)) +
     .    2*gs2tau12**2*dble(SU_B0(qs1,msta1,msta2)) +
     .  2*gs2n1n1**2*dble(SU_B0(qs1,msn1,msn1)) +
     .  1*gs2n1n1**2*dble(SU_B0(qs1,msntau,msntau)) +
     .   2*gs2e1e1**2*dble(SU_B0(qs1,mse1,mse1))+
     .   2*gs2e2e2**2*dble(SU_B0(qs1,mse2,mse2)) 
c
      pis2s2v = g**2/4* (
     . cbet**2*(2*dble(SU_BF(qs1,mch,mw))+dble(SU_BF(qs1,mA,mz))/cw**2)+
     . sbet**2*(2*dble(SU_BF(qs1,mw,mw)) +dble(SU_BF(qs1,mz,mz))/cw**2)+
     . 7*sbet**2*( 2*mw**2*dble(SU_B0(qs1,mw,mw))+mz**2/cw**2*
     . dble(SU_B0(qs1,mz,mz)) )  + 4*(2*SU_A(mw) +SU_A(mz)/cw**2)   )
c
      pis2s2h3 = g**2*mz**2/(4*cw**2)/2*(
     . (sbet*(3*sal**2-cal**2) -cbet*s2al)**2*dble(SU_B0(qs1,MH,MH)) +
     . (2*sbet*s2al-cbet*dcos(2*alfa))**2*dble(SU_B0(qs1,MH,ml)) +
     . (2*sbet*s2al-cbet*dcos(2*alfa))**2*dble(SU_B0(qs1,ml,MH)) +
     . (sbet*(3*cal**2-sal**2)+cbet*s2al)**2*dble(SU_B0(qs1,ml,ml)) +
     .  (sbet*dcos(2*b))**2*dble(SU_B0(qs1,mz,mz)) +
     .  (sbet*dsin(2*b))**2*dble(SU_B0(qs1,mz,MA)) +
     .  (sbet*dsin(2*b))**2*dble(SU_B0(qs1,MA,mz)) +
     .  (sbet*dcos(2*b))**2*dble(SU_B0(qs1,MA,MA)) +   
     .  2*( (-sbet*dcos(2*b))**2*dble(SU_B0(qs1,mw,mw)) +
     .     (sbet*dsin(2*b)-cw**2*cbet)**2*dble(SU_B0(qs1,mw,mch))+
     .    (sbet*dsin(2*b)-cw**2*cbet)**2*dble(SU_B0(qs1,mch,mw)) +
     .    (sbet*dcos(2*b)+2*cw**2*sbet)**2*dble(SU_B0(qs1,mch,mch)) ) )
c
      pis2s2h4 = g**2/(4*cw**2)/2*(
     . (3*sal**2-cal**2)*SU_A(MH)  + 
     . (3*cal**2-sal**2)*SU_A(ml)  - 
     .  dcos(2*b)*SU_A(mz) +dcos(2*b)*SU_A(MA) +
     . 2*( (cw**2 -sw2*dcos(2*b))*SU_A(mw) +
     . (cw**2 +sw2*dcos(2*b))*SU_A(mch) ) )
c
      pis2s2b=pis2s2v+pis2s2h3+pis2s2h4
c      
      pis2s2Nino =0.d0
      do i=1,4
      do j=1,4      
      pis2s2Nino =  pis2s2Nino + .5d0/4*2*(
     . ( ( Z(i,1)*Z(j,4) + Z(j,1)*Z(i,4) )*g1 - 
     .   ( Z(i,2)*Z(j,4) + Z(j,2)*Z(i,4) ) *g)**2 *
     .  (dble(SU_BG(qs1,gmn(i),gmn(j))) -      
     .   2*dxmn(i)*dxmn(j)*dble(SU_B0(qs1,gmn(i),gmn(j)))  ) )      
      enddo
      enddo
c      
      pis2s2Cino =0.d0
      do i=1,2
      do j=1,2      
      pis2s2Cino =     pis2s2Cino +g**2/2 *(
     .((V(i,2)*U(j,1))**2+(U(i,1)*V(j,2))**2)*
     .  dble(SU_BG(qs1,gmc(i),gmc(j)))
     . -4*V(i,2)*U(j,1)*U(i,1)*V(j,2)*gmc(i)*gmc(j)*
     . dble(SU_B0(qs1,gmc(i),gmc(j)))  )
      enddo
      enddo
c
c ------ Sum everything:
      pis2s2=1.d0/(16*pi**2)*
     . (pis2s2f+pis2s2s+pis2s2b+pis2s2Nino+pis2s2Cino)
c     
c----------------------------------------------------------------
c                        pis1s2 
c------------------------------------------------------------------
c
       pis1s2s = 
     .  3*(gs1t1t1*gs2t1t1*dble(SU_B0(qs1,mst1,mst1))+
     .     gs1t2t2*gs2t2t2*dble(SU_B0(qs1,mst2,mst2)) +
     .    2*gs1t1t2*gs2t1t2*dble(SU_B0(qs1,mst1,mst2)) ) +
     .  3*(gs1b1b1*gs2b1b1*dble(SU_B0(qs1,msb1,msb1))+
     .     gs1b2b2*gs2b2b2*dble(SU_B0(qs1,msb2,msb2)) +
     .    2*gs1b1b2*gs2b1b2*dble(SU_B0(qs1,msb1,msb2)) ) +
     .  3*2*(gs1u1u1*gs2u1u1*dble(SU_B0(qs1,msu1,msu1))+
     .     gs1u2u2*gs2u2u2*dble(SU_B0(qs1,msu2,msu2)) ) +
     .  3*2*(gs1d1d1*gs2d1d1*dble(SU_B0(qs1,msd1,msd1))+
     .     gs1d2d2*gs2d2d2*dble(SU_B0(qs1,msd2,msd2)) ) +
     .    gs1tau11*gs2tau11*dble(SU_B0(qs1,msta1,msta1))+
     .     gs1tau22*gs2tau22*dble(SU_B0(qs1,msta2,msta2)) +
     .    2*gs1tau12*gs2tau12*dble(SU_B0(qs1,msta1,msta2)) +
     .  2*gs1n1n1*gs2n1n1*dble(SU_B0(qs1,msn1,msn1)) +
     .  1*gs1n1n1*gs2n1n1*dble(SU_B0(qs1,msntau,msntau)) +
     .   2*gs1e1e1*gs2e1e1*dble(SU_B0(qs1,mse1,mse1))+
     .   2*gs1e2e2*gs2e2e2*dble(SU_B0(qs1,mse2,mse2)) 
c
      pis1s2v = g**2/4* sbet*cbet*(2*dble(SU_BF(qs1,mw,mw)) -
     . 2*dble(SU_BF(qs1,mch,mw))  +
     . (dble(SU_BF(qs1,mz,mz)) -dble(SU_BF(qs1,MA,mz)))/cw**2  +
     . 7*(2*mw**2*dble(SU_B0(qs1,mw,mw))+mz**2/cw**2*
     . dble(SU_B0(qs1,mz,mz)) )  ) 
c
      pis1s2h3 = g**2*mz**2/(4*cw**2)/2*(
     . (cbet*(3*cal**2-sal**2)-sbet*s2al)*
     . (sbet*(3*sal**2-cal**2)-cbet*s2al)*dble(SU_B0(qs1,MH,MH)) +
     . (-2*cbet*s2al-sbet*dcos(2*alfa))*
     . (2*sbet*s2al-cbet*dcos(2*alfa))*dble(SU_B0(qs1,MH,ml)) +
     . (-2*cbet*s2al-sbet*dcos(2*alfa))*
     . (2*sbet*s2al-cbet*dcos(2*alfa))*dble(SU_B0(qs1,ml,MH)) +
     . (cbet*(3*sal**2-cal**2)+sbet*s2al)*
     . (sbet*(3*cal**2-sal**2)+cbet*s2al)*dble(SU_B0(qs1,ml,ml)) +
     .  (cbet*dcos(2*b))*(-sbet*dcos(2*b))*dble(SU_B0(qs1,mz,mz)) +
     .  (-cbet*dsin(2*b))*(sbet*dsin(2*b))*dble(SU_B0(qs1,mz,MA)) +
     .  (-cbet*dsin(2*b))*(sbet*dsin(2*b))*dble(SU_B0(qs1,MA,mz)) +
     .  (-cbet*dcos(2*b))*(sbet*dcos(2*b))*dble(SU_B0(qs1,MA,MA)) +   
     .  2*( (cbet*dcos(2*b))*(-sbet*dcos(2*b))*dble(SU_B0(qs1,mw,mw))+
     .     (-cbet*dsin(2*b)+cw**2*sbet)*
     .     (sbet*dsin(2*b)-cw**2*cbet)*dble(SU_B0(qs1,mw,mch))+
     .     (-cbet*dsin(2*b)+cw**2*sbet)*
     .     (sbet*dsin(2*b)-cw**2*cbet)*dble(SU_B0(qs1,mch,mw))+
     .    (-cbet*dcos(2*b)+2*cw**2*cbet)*
     .    (sbet*dcos(2*b)+2*cw**2*sbet)*dble(SU_B0(qs1,mch,mch)) ) )
c
      pis1s2h4 = g**2/(4*cw**2)/2*(-s2al*SU_A(MH)  + s2al*SU_A(ml)
     .      -2*cw**2*dsin(2*b)*SU_A(mw)+2*cw**2*dsin(2*b)*SU_A(Mch) )      
c
      pis1s2b=pis1s2v+pis1s2h3+pis1s2h4
c      
      pis1s2Nino =0.d0
      do i=1,4
      do j=1,4      
      pis1s2Nino =  pis1s2Nino + .5d0/4*2*(
     . ( - ( Z(i,1)*Z(j,3) + Z(j,1)*Z(i,3) )*g1 + 
     .     ( Z(i,2)*Z(j,3) + Z(j,2)*Z(i,3) )*g )*
     . (   ( Z(i,1)*Z(j,4) + Z(j,1)*Z(i,4) )*g1 - 
     .     ( Z(i,2)*Z(j,4) + Z(j,2)*Z(i,4) )*g )*
     .  (dble(SU_BG(qs1,gmn(i),gmn(j))) -      
     . 2*dxmn(i)*dxmn(j)*dble(SU_B0(qs1,gmn(i),gmn(j)))  ) )       
      enddo
      enddo
c
      pis1s2Cino =0.d0
      do i=1,2
      do j=1,2      
      pis1s2Cino =     pis1s2Cino +g**2/2 *(
     .( V(i,1)*U(j,2)*V(i,2)*U(j,1)+
     .  U(i,1)*V(j,2)*U(i,2)*V(j,1) )*dble(SU_BG(qs1,gmc(i),gmc(j)))  
     . -2*( V(i,1)*U(j,2)*U(i,1)*V(j,2)
     .    + V(i,2)*U(j,1)*U(i,2)*V(j,1) )*gmc(i)*gmc(j)*      
     . dble(SU_B0(qs1,gmc(i),gmc(j)))  ) 
      enddo
      enddo
c
c ------ Sum everything:
      pis1s2=1.d0/(16*pi**2)*
     . (pis1s2s+pis1s2b+pis1s2Nino+pis1s2Cino)
c      
c-----------------------------------------------------------------
c                           piAA 
c-----------------------------------------------------------------

      if(mApole.eq.0d0) then   
         qsa= ma**2
      else
         qsa=mApole**2
      endif
c       
       piaaf=cbet**2*3*yt**2*(qsa*dble(SU_B0(qsa,rmt,rmt))
     .                      -2*SU_A(rmt))+
     .         sbet**2*3*yb**2*(qsa*dble(SU_B0(qsa,rmb,rmb))
     .                       -2*SU_A(rmb))+
     .   sbet**2*ytau**2*(qsa*dble(SU_B0(qsa,rmtau,rmtau))
     .                       -2*SU_A(rmtau))
c
        piaastop = 3*2*gat1t2**2*dble(SU_B0(qsa,mst1,mst2))  +
     .  3*(yt**2*cbet**2-g**2/cw**2*0.5d0*(0.5d0-2*sw2/3)*dcos(2*b))*
     .   (ct**2*SU_A(mst1)+st**2*SU_A(mst2) ) +
     .  3*(yt**2*cbet**2-g**2/cw**2*0.5d0*(2*sw2/3)*dcos(2*b))*
     .   (st**2*SU_A(mst1)+ct**2*SU_A(mst2) ) 
c     
        piaasup = 
     .  6*(-g**2/cw**2*0.5d0*(0.5d0-2*sw2/3)*dcos(2*b))*SU_A(msu1)+
     .  6*(-g**2/cw**2*0.5d0*(2*sw2/3)*dcos(2*b))*SU_A(msu2) 
c
        piaasbot=  3*2*gab1b2**2*dble(SU_B0(qsa,msb1,msb2))  +
     .  3*(yb**2*sbet**2 -g**2/cw**2*0.5d0*(-0.5d0+sw2/3)*dcos(2*b))*
     .   (cb**2*SU_A(msb1)+sb**2*SU_A(msb2) ) +
     .  3*(yb**2*sbet**2 -g**2/cw**2*0.5d0*(-sw2/3)*dcos(2*b))*
     .   (sb**2*SU_A(msb1)+cb**2*SU_A(msb2) ) 
c     
        piaasdo=  
     .  6*(-g**2/cw**2*0.5d0*(-0.5d0+sw2/3)*dcos(2*b))*SU_A(msd1)+
     .  6*(-g**2/cw**2*0.5d0*(-sw2/3)*dcos(2*b))*SU_A(msd2) 
c     
        piaastau= 2*gatau12**2*dble(SU_B0(qsa,msta1,msta2))+
     .  (ytau**2*sbet**2 -g**2/cw**2*0.5d0*(-0.5d0+sw2)*dcos(2*b))*
     .   (cta**2*SU_A(msta1)+sta**2*SU_A(msta2) ) +
     .  (ytau**2*sbet**2 -g**2/cw**2*0.5d0*(-sw2)*dcos(2*b))*
     .   (sta**2*SU_A(msta1)+cta**2*SU_A(msta2) ) 
c
        piaaslep=
     .  2*(-g**2/cw**2*0.5d0*(-0.5d0+sw2)*dcos(2*b))*SU_A(mse1)+
     .  2*(-g**2/cw**2*0.5d0*(-sw2)*dcos(2*b))*SU_A(mse2)
c
        piaasneu=
     .  2*(-g**2/cw**2*0.5d0*(0.5d0)*dcos(2*b))*SU_A(msn1)+
     .  1*(-g**2/cw**2*0.5d0*(0.5d0)*dcos(2*b))*SU_A(msntau)
c 
       piaas=piaastop+piaasbot+piaastau+piaasup+piaasdo+piaaslep
     .      +piaasneu
c      
       piaav = g**2/4* (2*dble(SU_BF(qsa,mch,mw))+
     . dsin(alfa-beta)**2*dble(SU_BF(qsa,MH,mz))/cw**2+
     . dcos(alfa-beta)**2*dble(SU_BF(qsa,ml,mz))/cw**2 )
     .       +g**2*mw**2/2*dble(SU_B0(qsa,mw,mch))
     .       +g**2*(2*SU_A(mw)+SU_A(mz)/cw**2) 
c
      piaah=g**2/(4*cw**2)/2*(
     . (3*dsin(2*beta)**2-1.d0)*SU_A(Mz)+3*dcos(2*beta)**2*SU_A(ma)+ 
     .   dcos(2*b)*dcos(2*alfa)*(SU_A(ml) -SU_A(MH)) +
     . 2*((cw**2*(1.d0+dsin(2*beta)**2)-sw**2*dcos(2*beta)**2)*SU_A(mw)
     .    +dcos(2*beta)**2*SU_A(mch) ) )
     .    +g**2*mz**2/(4*cw**2)/2*(
     .(cal*(-dcos(2*b)*cbet)+sal*dcos(2*b)*sbet)**2*
     . dble(SU_B0(qsa,ma,MH))+
     .(cal*(-dcos(2*b)*cbet)+sal*dcos(2*b)*sbet)**2*
     . dble(SU_B0(qsa,MH,ma))+
     .(sal*(dcos(2*b)*cbet)+cal*dcos(2*b)*sbet)**2*
     . dble(SU_B0(qsa,ma,ml))+
     .(sal*(dcos(2*b)*cbet)+cal*dcos(2*b)*sbet)**2*
     . dble(SU_B0(qsa,ml,ma))+
     .(cal*(-dsin(2*b)*cbet)+sal*dsin(2*b)*sbet)**2*
     . dble(SU_B0(qsa,mz,mh))+
     .(cal*(-dsin(2*b)*cbet)+sal*dsin(2*b)*sbet)**2*
     . dble(SU_B0(qsa,mh,mz))+
     .(sal*(dsin(2*b)*cbet)+cal*dsin(2*b)*sbet)**2*
     . dble(SU_B0(qsa,mz,ml))+
     .(sal*(dsin(2*b)*cbet)+cal*dsin(2*b)*sbet)**2*
     . dble(SU_B0(qsa,ml,mz))  )
c
      piaab=piaav+piaah
c

      piaaNino =0.d0
      do i=1,4
      do j=1,4      
      piaaNino =  piaaNino + .5d0/4*2*(
     .    sbet*g1*( Z(i,1)*Z(j,3) +Z(i,3)*Z(j,1) )
     . -  sbet*g* ( Z(i,2)*Z(j,3) +Z(i,3)*Z(j,2) ) 
     . -  cbet*g1*( Z(i,1)*Z(j,4) +Z(i,4)*Z(j,1) )
     . +  cbet*g* ( Z(i,2)*Z(j,4) +Z(i,4)*Z(j,2) ) )**2*   
     .  ( dble(SU_BG(qsa,gmn(i),gmn(j))) +      
     .    2*dxmn(i)*dxmn(j)*dble(SU_B0(qsa,gmn(i),gmn(j)))  )          
      enddo
      enddo
c
      piaaCino =0.d0
      do i=1,2
      do j=1,2      
      piaaCino =     piaaCino +g**2/2 *(
     .  ( (-sbet*V(i,1)*U(j,2)-cbet*V(i,2)*U(j,1))**2+
     .    ( sbet*U(i,2)*V(j,1)+cbet*U(i,1)*V(j,2))**2)*
     . dble(SU_BG(qsa,gmc(i),gmc(j)))
     . - 4.d0*(-sbet*V(i,1)*U(j,2)-cbet*V(i,2)*U(j,1))*
     .        ( sbet*U(i,2)*V(j,1)+cbet*U(i,1)*V(j,2))*
     .     gmc(i)*gmc(j)*dble(SU_B0(qsa,gmc(i),gmc(j)))  )
      enddo
      enddo
c
c ------ Sum everything:
      piaa= 1.d0/(16*pi**2)*(piaaf+piaas+piaab+piaaNino+piaaCino)
c       
c------------------------------------------------------------------
c                               piH+H-
c------------------------------------------------------------------ 
      if(mCHpole.eq.0d0) then  
         qsc= mch**2
      else
         qsc=mCHpole**2
      endif
c       
       piccf = 3*(yt**2*cbet**2 +yb**2*sbet**2)*dble(SU_BG(qsc,rmt,rmb)) 
     . -2*yt*yb*rmt*rmb*dsin(2*b)*dble(SU_B0(qsc,rmt,rmb))  +
     .  ytau**2*sbet**2*dble(SU_BG(qsc,eps,rmtau))       
c
       piccs= 
     . 3*gct1b1**2*dble(SU_B0(qsc,mst1,msb1)) +
     . 3*gct1b2**2*dble(SU_B0(qsc,mst1,msb2)) +
     . 3*gct2b1**2*dble(SU_B0(qsc,mst2,msb1)) +
     . 3*gct2b2**2*dble(SU_B0(qsc,mst2,msb2)) +
c     . gctau11**2*dble(SU_B0(qsc,msn1,msta1)) +
c     . gctau12**2*dble(SU_B0(qsc,msn1,msta2)) +
c correction msn1 -> msntau:
     . gctau11**2*dble(SU_B0(qsc,msntau,msta1)) +
     . gctau12**2*dble(SU_B0(qsc,msntau,msta2)) +
     . gceLL**2*(6*dble(SU_B0(qsc,msu1,msd1)) 
     .          +2*dble(SU_B0(qsc,msn1,mse1)))+
     .  3*(yb**2*sbet**2-g**2/cw**2*0.5d0*(0.5d0-2*sw2/3)*dcos(2*b)+
     . g**2/2*dcos(2*b))*(ct**2*SU_A(mst1)+st**2*SU_A(mst2) ) +
     .  3*(yt**2*cbet**2-g**2/cw**2*0.5d0*(2*sw2/3)*dcos(2*b))*
     .   (st**2*SU_A(mst1)+ct**2*SU_A(mst2) ) +
     .  6*(-g**2/cw**2*0.5d0*(0.5d0-2*sw2/3)*dcos(2*b)+g**2/2*
     .  dcos(2*b))*SU_A(msu1)+
     .  6*(-g**2/cw**2*0.5d0*(2*sw2/3)*dcos(2*b))*SU_A(msu2) +
     .  (-g**2/cw**2*0.5d0*(0.5d0)*dcos(2*b)+g**2/2*dcos(2*b))*
     .  2*SU_A(msn1) 
     . + (ytau**2*sbet**2-g**2/cw**2*0.5d0*(0.5d0)*dcos(2*b)
     . + g**2/2*dcos(2*b))*SU_A(msntau) +
     .  3*(yt**2*cbet**2 +g**2/cw**2*(-0.5d0)*(-0.5d0+sw2/3)*dcos(2*b)-
     . g**2/2*dcos(2*b))*(cb**2*SU_A(msb1)+sb**2*SU_A(msb2) ) +
     .  3*(yb**2*sbet**2 +g**2/cw**2*(-0.5d0)*(-sw2/3)*dcos(2*b))*
     .   (sb**2*SU_A(msb1)+cb**2*SU_A(msb2) ) +
     .  6*(g**2/cw**2*(-0.5d0)*(-0.5d0+sw2/3)*dcos(2*b) -g**2/2*
     .  dcos(2*b))*SU_A(msd1)+
     .  6*( g**2/cw**2*(-0.5d0)*(-sw2/3)*dcos(2*b))*SU_A(msd2) +
     .  (ytau**2*cbet**2 +g**2/cw**2*(-0.5d0)*(-0.5d0+sw2)*dcos(2*b)-
     . g**2/2*dcos(2*b))*(cta**2*SU_A(msta1)+sta**2*SU_A(msta2) ) +
     .  (ytau**2*sbet**2 +g**2/cw**2*(-0.5d0)*(-sw2)*dcos(2*b))*
     .   (sta**2*SU_A(msta1)+cta**2*SU_A(msta2) ) +
     .  2*( g**2/cw**2*(-0.5d0)*(-0.5d0+sw2)*dcos(2*b) -g**2/2*
     .  dcos(2*b))*SU_A(mse1) +
     .  2*( g**2/cw**2*(-0.5d0)*(-sw2)*dcos(2*b))*SU_A(mse2)
c     
       piccv = g**2/4*(dsin(alfa-b)**2*dble(SU_BF(qsc,MH,mw)) +
     .   dcos(alfa-b)**2*dble(SU_BF(qsc,ml,mw))+dble(SU_BF(qsc,MA,mw))
     . +(cw**2-sw2)**2/cw**2 *dble(SU_BF(qsc,mch,mz)) )
     . +4*pi*alph* dble(SU_BF(qsc,mch,eps))+2*g**2*SU_A(mw)+
     . g**2*(cw**2-sw2)**2/cw**2 *SU_A(mz) +
     . g**2*mw**2/4* dble(SU_B0(qsc,mw,ma)) 
c
       picch4= g**2/(4*cw**2)/2*(
     . (cw**2*(1.d0+dsin(2*b)**2) -sw**2*dcos(2*b)**2)*SU_A(mz) +
     . dcos(2*b)**2*SU_A(mA) +
     . (cw**2*(1.d0+dsin(2*b)*dsin(2*alfa))
     . -sw**2*dcos(2*b)*dcos(2*alfa))*SU_A(MH) +  
     . (cw**2*(1.d0-dsin(2*b)*dsin(2*alfa))
     . +sw**2*dcos(2*b)*dcos(2*alfa))*SU_A(ml) ) +
     . g**2/(4*cw**2)* ( 
     . (2*dsin(2*b)**2-1.d0)*SU_A(mw) +2*dcos(2*b)**2*SU_A(mch) )  
c 
      picch3= g**2*mz**2/(4*cw**2)/2 *(
     .(cal*(-dsin(2*b)*cbet+cw**2*sbet)+
     . sal*(dsin(2*b)*sbet-cw**2*cbet))**2*dble(SU_B0(qsc,MH,mw)) +
     .(-sal*(-dsin(2*b)*cbet+cw**2*sbet)+
     . cal*(dsin(2*b)*sbet-cw**2*cbet))**2*dble(SU_B0(qsc,ml,mw)) +
     .(cal*(-dcos(2*b)*cbet+2*cw**2*cbet)+
     . sal*(dcos(2*b)*sbet+2*cw**2*sbet))**2*dble(SU_B0(qsc,MH,mch)) +
     .(-sal*(-dcos(2*b)*cbet+2*cw**2*cbet)+
     . cal*(dcos(2*b)*sbet+2*cw**2*sbet))**2*dble(SU_B0(qsc,ml,mch)) )
c --> add the  H+AG- term not present in PBMZ: 
     . + g**2*mw**2/4.d0*dble(SU_B0(qsc,ma,mw))
c
      piccb=piccv+picch3+picch4
c      
      piccino =0.d0
      do i=1,4
      do j=1,2      
      piccino =  piccino + g**2/2*( (
     .(-sbet*(Z(i,1)*U(j,2)*sw/cw+Z(i,2)*U(j,2)-sq2*Z(i,3)*U(j,1)))**2+ 
     .(cbet*(Z(i,1)*V(j,2)*sw/cw+Z(i,2)*V(j,2)+sq2*Z(i,4)*V(j,1)))**2)* 
     .  dble(SU_BG(qsc,gmc(j),gmn(i))) -      
     . 4*(-sbet*(Z(i,1)*U(j,2)*sw/cw+Z(i,2)*U(j,2)-sq2*Z(i,3)*U(j,1)))* 
     .  (cbet*(Z(i,1)*V(j,2)*sw/cw+Z(i,2)*V(j,2)+sq2*Z(i,4)*V(j,1)))* 
     .   gmc(j)*dxmn(i)*dble(SU_B0(qsc,gmc(j),gmn(i)) )   )        
      enddo
      enddo
c
      piccino =0.d0
      do i=1,4
      do j=1,2      

      fcoeff = g**2/2*       
     .((-sbet*(Z(i,1)*U(j,2)*sw/cw+Z(i,2)*U(j,2)-sq2*Z(i,3)*U(j,1)))**2 
     .+(cbet*(Z(i,1)*V(j,2)*sw/cw+Z(i,2)*V(j,2)+sq2*Z(i,4)*V(j,1)))**2)

      gcoeff = g**2*   
     . (-sbet*(Z(i,1)*U(j,2)*sw/cw+Z(i,2)*U(j,2)-sq2*Z(i,3)*U(j,1)))* 
     . (cbet*(Z(i,1)*V(j,2)*sw/cw+Z(i,2)*V(j,2)+sq2*Z(i,4)*V(j,1)))
      
      piccino =  piccino +  fcoeff* dble(SU_BG(qsc,gmc(j),gmn(i))) 
     . - 2*gcoeff*gmc(j)*dxmn(i)*dble(SU_B0(qsc,gmc(j),gmn(i)))

      enddo
      enddo  
c          
c ------ Sum everything:
      picc= 1.d0/(16*pi**2)*(piccf+piccs+piccb+piccino)
c       
c---------------------------------------------------------------------
c                     Tadpoles t1/v1 and t2/v2
c---------------------------------------------------------------------
c       
       dt1v1f = -6*yb**2*SU_A(rmb)
     .          -2*ytau**2*SU_A(rmtau)
c
       dt1v1s=     
     . g/2.d0/mw/cbet*(3*gs1t1t1*SU_A(mst1)+3*gs1t2t2*SU_A(mst2)
     .                +3*gs1b1b1*SU_A(msb1)+3*gs1b2b2*SU_A(msb2)
     .                +gs1tau11*SU_A(msta1)+gs1tau22*SU_A(msta2)
     .                +6*gs1u1u1*SU_A(msu1)+6*gs1u2u2*SU_A(msu2)
     .                +6*gs1d1d1*SU_A(msd1)+6*gs1d2d2*SU_A(msd2)
     .                +2*gs1e1e1*SU_A(mse1)+2*gs1e2e2*SU_A(mse2)
     .                +2*gs1n1n1*SU_A(msn1)+gs1n1n1*SU_A(msntau)  ) 
c
       dt1v1h=-g**2*dcos(2*b)/(8*cw**2)*(SU_A(ma)+2*SU_A(mch))
     . +g**2/2*SU_A(mch)
     . +g**2/(8*cw**2)*(3*sal**2-cal**2+s2al*tbeta)*SU_A(ml)
     . +g**2/(8*cw**2)*(3*cal**2-sal**2-s2al*tbeta)*SU_A(mh)

       dt1v1v=
     . +3*g**2/4.d0*(2*SU_A(mw)+SU_A(mz)/cw**2) 
     . +g**2*dcos(2*b)/(8*cw**2)*(2*SU_A(mw)+SU_A(mz))

       dt1v1b=dt1v1h+dt1v1v
c
       dt1v1nino = 0.d0
       do i=1,4
       dt1v1nino = dt1v1nino - g**2/mw/cbet*dxmn(i)*
     .           Z(i,3)*( Z(i,2)-Z(i,1)*sw/cw )*SU_A(gmn(i))
       enddo
c
       dt1v1cino = 0.d0
       do j=1,2
       dt1v1cino = dt1v1cino - g**2*sq2/mw/cbet*gmc(j)*
     .           V(j,1)*U(j,2)*SU_A(gmc(j))
       enddo
c--sum       
       tad1=(dt1v1f+dt1v1s+dt1v1b+dt1v1nino+dt1v1cino)/(16*pi**2)
       dVdvd2=-tad1

c --------------------------------------------------------------------
c
       dt2v2f = -6*yt**2*SU_A(rmt)
c       
       dt2v2s=
     . g/2.d0/mw/sbet*(3*gs2t1t1*SU_A(mst1)+3*gs2t2t2*SU_A(mst2)
     .                +3*gs2b1b1*SU_A(msb1)+3*gs2b2b2*SU_A(msb2)
     .                +gs2tau11*SU_A(msta1)+gs2tau22*SU_A(msta2)
     .                +6*gs2u1u1*SU_A(msu1)+6*gs2u2u2*SU_A(msu2)
     .                +6*gs2d1d1*SU_A(msd1)+6*gs2d2d2*SU_A(msd2)
     .                +2*gs2e1e1*SU_A(mse1)+2*gs2e2e2*SU_A(mse2)
     .                +2*gs2n1n1*SU_A(msn1)+ gs2n1n1*SU_A(msntau) ) 
c
       dt2v2b=g**2*dcos(2*b)/(8*cw**2)*(SU_A(ma)+2*SU_A(mch))
     . +g**2/2*SU_A(mch)
     . +g**2/(8*cw**2)*(3*cal**2-sal**2+s2al/tbeta)*SU_A(ml)
     . +g**2/(8*cw**2)*(3*sal**2-cal**2-s2al/tbeta)*SU_A(mh)
     . +3*g**2/4.d0*(2*SU_A(mw)+SU_A(mz)/cw**2) 
     . -g**2*dcos(2*b)/(8*cw**2)*(2*SU_A(mw)+SU_A(mz)) 
c     
       dt2v2nino = 0.d0
       do i=1,4
       dt2v2nino = dt2v2nino + g**2/mw/sbet*dxmn(i)*
     .           Z(i,4)*( Z(i,2)-Z(i,1)*sw/cw)*SU_A(gmn(i))
       enddo
c
       dt2v2cino = 0.d0
       do j=1,2
       dt2v2cino = dt2v2cino - g**2*sq2/mw/sbet*gmc(j)*
     .           V(j,2)*U(j,1)*SU_A(gmc(j))
       enddo
c--sum       
       tad2=(dt2v2f+dt2v2s+dt2v2b+dt2v2nino+dt2v2cino)/(16*pi**2)
       dVdvu2=-tad2
c
       mz = mzsave             
       mw = mwsave
       mt = mtsave
       mb = mbsave
       mtau = mtausave
       alfa = alfasave

        end
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cc  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   The following routine is for the one loop effective scalar potential 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE SU_VLOOP2(q2,MU,AT,AB,AL,dVdvd2,dVdvu2,pizz) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
c  The main subroutine for the EWSB and calculates the tadpole corrections to
c  the Higgs mass terms squared. The input are:
c  q2: the scale at which EWSB is supposed to happen,
c  MU: the higgsino parameter mu ar EWSB scale
c  AT,AB,AL: the third generation trilinear couplings at EWSB scale
c  Ytau, Yt, Yb: the Yukawa couplings (at EWSB scale)
c  msta1,msta2,msb1,msb2,mst1,mst2,..,thet,theb,thel: masses and mixing of
c  tau,b,top,.. etc sfermions at EWSB scale (input via common/su_bpew/..)
c  Other important input parameters, such as the Higgs, chargino, neutralino 
c  masses and couplings as well as SM parameters are called via commons.
c  The output are dVdvd2, dVdvu2, which are (up to some appropriate overall 
c  constants) the derivatives of the full one-loop scalar potential including
c  the contributions of all SM and SUSY particles a la PBMZ (hep-ph/9606211).
c Another output is pizz which allow to calculated the RC to MZ**2. 
c  The consistency of the EWSB mechanism is performed by the main program
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,m,o-z)
      real*8 nf
      complex*16 su_b0, su_bh, su_bt22
      logical su_isNaN
      dimension u(2,2),v(2,2),z(4,4),dxmn(4),gmn(4),gmc(2)
      COMMON/SU_cte/nf,cpi,zm,wm,tbeta
      COMMON/SU_outhiggs/dml,dmh,dmch,alfa 
      COMMON/SU_outginos/mc1,mc2,mn1,mn2,mn3,mn4,mgluino
      COMMON/SU_matino/u,v,z,dxmn
      COMMON/SU_hflag/imodel
      COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_bpew/msu1,msu2,msd1,msd2,mse1,mse2,msn1,msntau,
     . msta1,msta2,msb1,msb2,mst1,mst2,thet,theb,thel
      COMMON/SU_break/ameldum,amerdum,amsq,amurdum,amdrdum,
     .       aldum,audum,addum,amudum,am1dum,am2dum,am3
      COMMON/SU_fmasses/mtau,mbpole,mtpole
      COMMON/SU_renscale/scale     

      COMMON/runhiggs/ma,ml,mh,mch 
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
      common/su_nonpert/inonpert
c
       pi = 4*datan(1.d0)
       scale= dsqrt(q2)      
c basic parameters and definitions used:
       g=g2ewsb
       gstrong=dsqrt(4.d0*pi*alsewsb)
       mz=zm
       mw=wm
c defining s^2_w at EWSB scale:
       cw = 1.d0/dsqrt(1.d0+(g1ewsb/g2ewsb)**2)
       sw = g1ewsb/g2ewsb *cw
       sw2=sw**2
       cw2= cw**2
       cwm2 =1.d0/cw2
       alph = (g*sw)**2/(4*pi)
c defining mtau,mb,mt running masses at ewsb scale:
       if(su_isNaN(pizz).or.mz**2+pizz.le.0.d0) then  !!added protection
       pizz=0.d0
       if(irge.eq.irgmax)  inonpert=-1
       endif
       vd2 = 2*(mz**2+pizz)/(g1ewsb**2+g2ewsb**2)/(1.d0+tbeta**2)
       rmz= dsqrt(mz**2+pizz)   ! (use RUNNING MW,MZ)
       vu2 = vd2*tbeta**2
       vev2=2.d0*(vu2+vd2)
       vd= dsqrt(vd2)
       vu= dsqrt(vu2)
       rmtau = ytau*vd
       rmb = yb*vd
       rmt = yt*vu
       rmw = rmz*cw
       mzsave = mz
       mwsave = mw
       mz = rmz
       mw = rmw
c
       ct=dcos(thet)
       st=dsin(thet)
       cb=dcos(theb)
       sb=dsin(theb)
       cta=dcos(thel)
       sta=dsin(thel)
c
       beta= datan(tbeta)
       cbeta2=1.d0/(1.d0+tbeta**2)
       cbet= dsqrt(cbeta2)
       sbet=dsqrt(1.d0-cbeta2)
       c2b =2*cbeta2-1.d0
c 
       alfasave = alfa          !  alfa running
       alfa =0.5*atan(tan(2d0*beta)*(ma**2+mz**2)
     .      /(ma**2-mz**2))
       if(cos(2d0*beta)*(ma**2-mz**2).gt.0) alfa = alfa - pi/2d0
c
       sal=dsin(alfa)
       cal=dcos(alfa)
       s2a = 2*sal*cal 
       tm= rmt
       bm= rmb
       taum= rmtau
c ! yt, yb, ytau at ewsb scale are taken from COMMON/SU_yukaewsb/..
c relevant Sfermion couplings contributions:
       sq2=dsqrt(2.d0)
c
       s1bLbL = g*mz/cw*(-.5d0 +sw2/3)*cbet +sq2*yb*bm
       s1bRbR = g*mz/cw*(-sw2/3)*cbet +sq2*yb*bm
       s1bLbR = yb/sq2*Ab
       s2bLbL = -g*mz/cw*(-.5d0 +sw2/3)*sbet
       s2bRbR = -g*mz/cw*(-sw2/3)*sbet
       s2bLbR = -yb/sq2*mu
c
       gs1d1d1 = g*mz/cw*(-.5d0 +sw2/3)*cbet 
       gs1d2d2 = g*mz/cw*(-sw2/3)*cbet 
       gs2d1d1 = -g*mz/cw*(-.5d0 +sw2/3)*sbet
       gs2d2d2 = -g*mz/cw*(-sw2/3)*sbet
c  
       s1tauLL = g*mz/cw*(-.5d0 +sw2)*cbet +sq2*ytau*taum
       s1tauRR = g*mz/cw*(-sw2)*cbet +sq2*ytau*taum
       s1tauLR = ytau/sq2*AL
       s2tauLL = -g*mz/cw*(-.5d0 +sw2)*sbet
       s2tauRR = -g*mz/cw*(-sw2)*sbet
       s2tauLR = -ytau/sq2*mu
c
       gs1e1e1 = g*mz/cw*(-.5d0 +sw2)*cbet 
       gs1e2e2 = g*mz/cw*(-sw2)*cbet 
       gs2e1e1 = -g*mz/cw*(-.5d0 +sw2)*sbet
       gs2e2e2 = -g*mz/cw*(-sw2)*sbet
c
       s2tLtL = -g*mz/cw*(.5d0 -2*sw2/3)*sbet +sq2*yt*tm
       s2tRtR = -g*mz/cw*(2*sw2/3)*sbet +sq2*yt*tm
       s2tLtR = yt/sq2*AT
       s1tLtL = g*mz/cw*(.5d0 -2*sw2/3)*cbet
       s1tRtR = g*mz/cw*(2*sw2/3)*cbet
       s1tLtR = -yt/sq2*mu
c
       gs2u1u1 = -g*mz/cw*(.5d0 -2*sw2/3)*sbet 
       gs2u2u2 = -g*mz/cw*(2*sw2/3)*sbet
       gs1u1u1 = g*mz/cw*(.5d0 -2*sw2/3)*cbet
       gs1u2u2 = g*mz/cw*(2*sw2/3)*cbet
c
       gs2n1n1 = -g*mz/cw*(.5d0 )*sbet 
       gs1n1n1 = g*mz/cw*(.5d0 )*cbet
c
       gs1b1b1 = cb**2*s1bLbL +2*cb*sb*s1bLbR +sb**2*s1bRbR 
       gs1b2b2 = sb**2*s1bLbL -2*cb*sb*s1bLbR +cb**2*s1bRbR 
       gs2b1b1 = cb**2*s2bLbL +2*cb*sb*s2bLbR +sb**2*s2bRbR 
       gs2b2b2 = sb**2*s2bLbL -2*cb*sb*s2bLbR +cb**2*s2bRbR 
c
       gs1t1t1 = ct**2*s1tLtL +2*ct*st*s1tLtR +st**2*s1tRtR 
       gs1t2t2 = st**2*s1tLtL -2*ct*st*s1tLtR +ct**2*s1tRtR 
       gs2t1t1 = ct**2*s2tLtL +2*ct*st*s2tLtR +st**2*s2tRtR 
       gs2t2t2 = st**2*s2tLtL -2*ct*st*s2tLtR +ct**2*s2tRtR 
c
       gs1tau11 = cta**2*s1tauLL +2*cta*sta*s1tauLR +sta**2*s1tauRR 
       gs1tau22 = sta**2*s1tauLL -2*cta*sta*s1tauLR +cta**2*s1tauRR 
       gs2tau11 = cta**2*s2tauLL +2*cta*sta*s2tauLR +sta**2*s2tauRR 
       gs2tau22 = sta**2*s2tauLL -2*cta*sta*s2tauLR +cta**2*s2tauRR 
c
c  down fermion contributions: 
       confd= -2*ytau*rmtau/vd*SU_A(rmtau)
     .      -2*3*yb*rmb/vd*SU_A(rmb)
c  up fermion contributions:
       confu = -2*3*yt*rmt/vu*SU_A(rmt)

c    vd sfermion contributions: 
       dvdsf= 1d0/sq2/vd*(3*gs1t1t1*SU_A(mst1)+3*gs1t2t2*SU_A(mst2) +
     . 3*gs1b1b1*SU_A(msb1)+3*gs1b2b2*SU_A(msb2) +   
     . gs1tau11*SU_A(msta1)+ gs1tau22*SU_A(msta2) +
     .  2*( 3*gs1u1u1*SU_A(msu1)+3*gs1u2u2*SU_A(msu2) +
     . 3*gs1d1d1*SU_A(msd1)+3*gs1d2d2*SU_A(msd2) +   
     . gs1e1e1*SU_A(mse1)+ gs1e2e2*SU_A(mse2) ) + 
     .  2*gs1n1n1*SU_A(msn1) +gs1n1n1*SU_A(msntau) )
c    vu sfermion contributions:
       dvusf= 1d0/sq2/vu*(3*gs2t1t1*SU_A(mst1)+3*gs2t2t2*SU_A(mst2) +
     . 3*gs2b1b1*SU_A(msb1)+3*gs2b2b2*SU_A(msb2) +   
     . gs2tau11*SU_A(msta1)+ gs2tau22*SU_A(msta2) +
     .  2*( 3*gs2u1u1*SU_A(msu1)+3*gs2u2u2*SU_A(msu2) +
     . 3*gs2d1d1*SU_A(msd1)+3*gs2d2d2*SU_A(msd2) +   
     . gs2e1e1*SU_A(mse1)+ gs2e2e2*SU_A(mse2) ) + 
     .  2*gs2n1n1*SU_A(msn1) +gs2n1n1*SU_A(msntau) )
c   vd Higgs contributions:
       dvdH = -g**2*c2b/cw**2/8*(SU_A(ma)+2*SU_A(mch))+g**2/2*SU_A(mch)
     . + g**2/cw**2/8*(3*sal**2-cal**2+s2a*tbeta)*SU_A(ml)
     . + g**2/cw**2/8*(3*cal**2-sal**2-s2a*tbeta)*SU_A(mh)
c   vd gauge contributions:  
       dvdWZ = 3*g**2/4*(2*SU_A(mw)+SU_A(mz)/cw**2) 
     . +g**2*c2b/cw**2/8*(2*SU_A(mw) +SU_A(mz) )
c   vd gaugino contributions:
       dvdino = -g**2/mw/cbet*
     .  (dxmn(1)*Z(1,3)*(Z(1,2)-Z(1,1)*sw/cw)*SU_A(mn1)
     . + dxmn(2)*Z(2,3)*(Z(2,2)-Z(2,1)*sw/cw)*SU_A(mn2)
     . + dxmn(3)*Z(3,3)*(Z(3,2)-Z(3,1)*sw/cw)*SU_A(mn3)
     . + dxmn(4)*Z(4,3)*(Z(4,2)-Z(4,1)*sw/cw)*SU_A(mn4)  )
     . -dsqrt(2.d0)*g**2/mw/cbet*(mc1*V(1,1)*U(1,2)*SU_A(mc1) 
     .                           +mc2*V(2,1)*U(2,2)*SU_A(mc2)  )
c
c   vu Higgs contributions:
       dvuH = g**2*c2b/cw**2/8*(SU_A(ma)+2*SU_A(mch))+g**2/2*SU_A(mch)
     . + g**2/cw**2/8*(3*cal**2-sal**2+s2a/tbeta)*SU_A(ml)
     . + g**2/cw**2/8*(3*sal**2-cal**2-s2a/tbeta)*SU_A(mh)
c   vu gauge contributions:  
       dvuWZ = 3*g**2/4*(2*SU_A(mw)+SU_A(mz)/cw**2) 
     . -g**2*c2b/cw**2/8*(2*SU_A(mw) +SU_A(mz) )
c   vu gaugino contributions:
       dvuino = g**2/mw/sbet*
     . (dxmn(1)*Z(1,4)*(Z(1,2)-Z(1,1)*sw/cw)*SU_A(mn1)
     . + dxmn(2)*Z(2,4)*(Z(2,2)-Z(2,1)*sw/cw)*SU_A(mn2)
     . + dxmn(3)*Z(3,4)*(Z(3,2)-Z(3,1)*sw/cw)*SU_A(mn3)
     . + dxmn(4)*Z(4,4)*(Z(4,2)-Z(4,1)*sw/cw)*SU_A(mn4)  )
     . -dsqrt(2.d0)*g**2/mw/sbet*(mc1*V(1,2)*U(1,1)*SU_A(mc1) 
     .                           +mc2*V(2,2)*U(2,1)*SU_A(mc2)  )
c     
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c 2-loop tadpole parts: 
      if(imodel.ge.2) then
c  !  1-loop versus 2-loop possible choice:
      call SU_ewsb2loop(rmt**2,am3,mst1**2,mst2**2,st,ct,
     .     scale**2,-mu,tbeta,vev2,gstrong,tad1st,tad2st)
     
      call SU_ewsb2loop(rmb**2,am3,msb1**2,msb2**2,sb,cb,
     .     scale**2,-mu,1.d0/tbeta,vev2,gstrong,tad2sb,tad1sb)

      call SU_DDStad(rmt**2,rmb**2,ma**2,mst1**2,mst2**2,
     .     msb1**2,msb2**2,st,ct,sb,cb,scale**2,-mu,tbeta,vev2,
     .     tad1w,tad2w)

      call SU_taubottad(rmtau**2,rmb**2,msta1**2,msta2**2,msb1**2,
     .      msb2**2,sta,cta,sb,cb,scale**2,-mu,tbeta,vev2,
     $      tad1bl,tad2bl) 

      call SU_tausqtad(rmtau**2,ma**2,msntau**2,msta1**2,msta2**2,
     .      sta,cta,scale**2,-mu,tbeta,vev2,
     $      tad1l,tad2l)
      else
      tad1st=0.d0
      tad2st=0.d0
      tad1sb=0.d0
      tad2sb=0.d0
      tad1w=0.d0
      tad2w=0.d0
      tad1bl=0.d0
      tad2bl=0.d0
      tad1l=0.d0
      tad2l=0.d0
      endif
c
c final contributions including 2-loop corrections:
       
        tad1 = -cpi*(confd + dvdsf +dvdH +dvdWZ +dvdino) 
        tad1loop= tad1st+tad1sb+tad1w+tad1l+tad1bl   	
        dVdvd2=tad1+tad1loop

        tad2 = -cpi*(confu + dvusf +dvuH +dvuWZ +dvuino)   
     	tad2loop=tad2st+tad2sb+tad2w+tad2l+tad2bl	
        dVdvu2=tad2+tad2loop	
c    
c-----------------------------------------------------------------
c                 Z boson self-energy at q**2=mz**2
c-----------------------------------------------------------------
       qsz=mzsave**2 
       mup=1.d-2
       mdo=1.d-2
       me=0.5d-3
       mmu=0.106d0
       ms=0.190d0
       mcq=1.40d0
       mt =rmt                 
       mb =rmb

       eps=1.d-2
       eps0=eps**2
       gmn(1)=dabs(dxmn(1))
       gmn(2)=dabs(dxmn(2))
       gmn(3)=dabs(dxmn(3))
       gmn(4)=dabs(dxmn(4))
       gmc(1)=mc1
       gmc(2)=mc2
c
      pizzf = 3*( (.5d0-2*sw2/3)**2+(2*sw2/3)**2)
     .*(dble(SU_BH(qsz,mt,mt))+dble(SU_BH(qsz,mcq,mcq))
     .       +dble(SU_BH(qsz,mup,mup)))
     .  + 3*((-.5d0+sw2/3)**2+(-sw2/3)**2)
     .*(dble(SU_BH(qsz,mb,mb))+dble(SU_BH(qsz,ms,ms))
     .       +dble(SU_BH(qsz,mdo,mdo)))
     .  + ((-.5d0+sw2)**2+(-sw2)**2)*(dble(SU_BH(qsz,me,me))
     .    +dble(SU_BH(qsz,mmu,mmu))+dble(SU_BH(qsz,mtau,mtau)))
     .  + .5d0**2*3*dble(SU_BH(qsz,eps,eps))
     .  -12*(.5d0-2*sw2/3)*(2*sw2/3)
     .  *(mt**2*dble(SU_B0(qsz,mt,mt))+mcq**2*dble(SU_B0(qsz,mcq,mcq))
     .  +mup**2*dble(SU_B0(qsz,mup,mup))) 
     .  -12*(-.5d0+sw2/3)*(-sw2/3)
     .  *(mb**2*dble(SU_B0(qsz,mb,mb))+ms**2*dble(SU_B0(qsz,ms,ms))
     .  +mdo**2*dble(SU_B0(qsz,mdo,mdo))) 
     .  -4*(-.5d0+sw2)*(-sw2)*(me**2*dble(SU_B0(qsz,me,me))+mmu**2
     .  *dble(SU_B0(qsz,mmu,mmu))+mtau**2*dble(SU_B0(qsz,mtau,mtau)))
c     
      pizzb = -2*cw**4*(2*qsz+mw**2-mz**2*sw**4/cw**2)
     . *dble(SU_B0(qsz,mw,mw))
     . -(8*cw**4+(cw2-sw2)**2)*dble(SU_BT22(qsz,mw,mw))
c
      pizzh0=-dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml))
      pizzhS = -dsin(beta-alfa)**2*(dble(SU_BT22(qsz,ma,mh))
     .  + dble(SU_BT22(qsz,mz,ml))-mz**2*dble(SU_B0(qsz,mz,ml)) )
     .        -dcos(beta-alfa)**2*(dble(SU_BT22(qsz,mz,mh))
     .  + dble(SU_BT22(qsz,ma,ml))-mz**2*dble(SU_B0(qsz,mz,mh)) )
     .  -(cw**2-sw**2)**2*dble(SU_BT22(qsz,mch,mch))
     .  -pizzh0
c
       pizzsu= -12*( (.5d0-2*sw2/3)*dcos(thet)**2
     .-(2*sw2/3)*dsin(thet)**2 )**2*dble(SU_BT22(qsz,mst1,mst1))
     .         -12*(-(.5d0-2*sw2/3)*dsin(thet)**2
     .+(2*sw2/3)*dcos(thet)**2 )**2*dble(SU_BT22(qsz,mst2,mst2))
     .      -24*( (.5d0)*dsin(thet)*dcos(thet) )**2
     .  *dble(SU_BT22(qsz,mst1,mst2))
     .    -24*(.5d0-2*sw2/3)**2*dble(SU_BT22(qsz,msu1,msu1))
     .    -24*(+2*sw2/3)**2*dble(SU_BT22(qsz,msu2,msu2))
c
      pizzsd= -12*( (-.5d0+sw2/3)*dcos(theb)**2
     .-(-sw2/3)*dsin(theb)**2)**2*dble(SU_BT22(qsz,msb1,msb1))
     .       -12*( -(-.5d0+sw2/3)*dsin(theb)**2
     .+(-sw2/3)*dcos(theb)**2)**2*dble(SU_BT22(qsz,msb2,msb2))
     .      -24*((-0.5d0)*dsin(theb)*dcos(theb))**2
     .  *dble(SU_BT22(qsz,msb1,msb2))
     .    -24*(-.5d0+sw2/3)**2*dble(SU_BT22(qsz,msd1,msd1))
     .    -24*(-sw2/3)**2*dble(SU_BT22(qsz,msd2,msd2))
c
      pizzsl=-4*( (-.5d0+sw2)*dcos(thel)**2
     .- (-sw2)*dsin(thel)**2 )**2*dble(SU_BT22(qsz,msta1,msta1))
     .       -4*( -(-.5d0+sw2)*dsin(thel)**2
     .  +(-sw2)*dcos(thel)**2 )**2*dble(SU_BT22(qsz,msta2,msta2))
     .      -8*((-.5d0)*dsin(thel)*dcos(thel))**2
     .  *dble(SU_BT22(qsz,msta1,msta2))
     .      -8*(-.5d0+sw2)**2*dble(SU_BT22(qsz,mse1,mse1))
     .       -8*(-sw2)**2*dble(SU_BT22(qsz,mse2,mse2))
     .       -8*(.5d0)**2*dble(SU_BT22(qsz,msn1,msn1))
     .       -4*(.5d0)**2*dble(SU_BT22(qsz,msntau,msntau))
c
      pizzs=pizzsl+pizzsd+pizzsu
c
      pizzn=0.d0
      do  i=1,4
      do  j=1,4
      pizzn = pizzn + 1.d0/4*(Z(i,3)*Z(j,3) -Z(i,4)*Z(j,4))**2*
     .       (dble(SU_BH(qsz,gmn(i),gmn(j)))
     .       -2*dxmn(i)*dxmn(j)*dble(SU_B0(qsz,gmn(i),gmn(j))) )
      enddo
      enddo
c
      pizzc=0.d0
      do i=1,2
      do j=1,2
      pizzc = pizzc +1.d0/4*( 
     .( ( 2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2) )**2+
     .  ( 2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2) )**2 )
     .            *dble(SU_BH(qsz,gmc(i),gmc(j))) 
     .     +4*(2*cw2*V(i,1)*V(j,1)+(cw2-sw2)*V(i,2)*V(j,2))*
     .        (2*cw2*U(i,1)*U(j,1)+(cw2-sw2)*U(i,2)*U(j,2))*
     .            gmc(i)*gmc(j)*dble(SU_B0(qsz,gmc(i),gmc(j))) )
      enddo
      enddo
c
c Sum of the susy contributions for pizz and final pizz(MZ**2) 
      pizzsm=alph/4.d0/pi/sw2/cw2*(pizzf+pizzb+pizzh0)
      pizzsusy=alph/4.d0/pi/sw2/cw2*
     .        (pizzhS+pizzs+pizzn+pizzc)
      pizz=pizzsm+pizzsusy
c
       mz = mzsave            
       mw = mwsave
       alfa = alfasave
        end   
c   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   ++++++++++++++ End of the routines for the effective potential ++++++++
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  The following routines are for the RGE evolution of the parameters 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SUBROUTINE SU_ODEINT(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     .                      SU_DERIVS,SU_RKQC) 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c This is the main subroutine (from Numerical Recipes) integrating
c (coupled) Ordinary Differential Equation 
c    
      implicit real*8(a-h,o-z)
      parameter (maxstp=100000,nmax=31,two=2.d0,zero=0.d0,tiny=1.d-30)
      COMMON/ode_PATH/kmax,kount,dxsav,xp(200),yp(31,200)
      dimension ystart(nvar),yscal(nmax),y(nmax),dydx(nmax)
      COMMON/SU_good/iflop
      external SU_DERIVS, SU_RKQC 
      x=x1
      h=dsign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      xsav=x-dxsav*two
      do 16 nstp=1,maxstp
        CALL SU_DERIVS(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+tiny
12      continue
        if(kmax.gt.0)then
          if(dabs(x-xsav).gt.dabs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.zero) h=x2-x
c  new modif to speed up RGE integration in the safe zone far from GUT : 
          if(x.lt.dlog(1d14)) then          
        CALL SU_RKQC(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,su_derivs)
          else
        CALL SU_RKQC(y,dydx,nvar,x,h1,eps,yscal,hdid,hnext,su_derivs)
	  endif
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.zero)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
c        if(dabs(hnext).lt.hmin) pause 'stepsize smaller than minimum.'
        iflop=0
        if(dabs(hnext).lt.hmin) then
c        write(*,'(a)') 'stepsize smaller than minimum.'
        iflop=1
        endif
        h=hnext
16    continue
      iflop=1
      return
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
      SUBROUTINE SU_RKQC(y,dydx,n,x,htry,eps,yscal,hdid,hnext,SU_DERIVS)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Fourth order Runge--Kutta numerical algorithms solving differential 
c  equations by Numerical Recipes. Needed by the SU_ODEINT above. 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      implicit real*8(a-h,o-z)
      parameter (nmax=31,fcor=.066666666666666667d0,
     *    one=1.d0,safety=0.9d0,errcon=6.d-6)
      COMMON/SU_good/iflop
      external SU_DERIVS
      dimension y(n),dydx(n),yscal(n),ytemp(nmax),ysav(nmax),dysav(nmax)
      pgrow=-0.20d0
      pshrnk=-0.25d0
      xsav=x
      do 11 i=1,n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
11    continue
      h=htry
1     hh= h/2
      CALL SU_RK4(ysav,dysav,n,xsav,hh,ytemp,su_derivs)
      x=xsav+hh
      CALL SU_DERIVS(x,ytemp,dydx)
      CALL SU_RK4(ytemp,dydx,n,x,hh,y,su_derivs)
      x=xsav+h
      if(x.eq.xsav) then
c      pause 'stepsize not significant in rkqc.'
c      write(*,'(a)') 'stepsize not significant in rkqc.'
      iflop =1
      return
      endif
      CALL SU_RK4(ysav,dysav,n,xsav,h,ytemp,su_derivs)
      errmax=0.d0
      do 12 i=1,n
        ytemp(i)=y(i)-ytemp(i)
        errmax=max(errmax,abs(ytemp(i)/yscal(i)))
12    continue
      errmax=errmax/eps
      if(errmax.gt.one) then
        h=safety*h*(errmax**pshrnk)
        goto 1
      else
        hdid=h
        if(errmax.gt.errcon)then
          hnext=safety*h*(errmax**pgrow)
        else
          hnext=4*h
        endif
      endif
      do 13 i=1,n
        y(i)=y(i)+ytemp(i)*fcor
13    continue
      return
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
          SUBROUTINE SU_RK4(y,dydx,n,x,h,yout,SU_DERIVS)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Fourth order Runge--Kutta numerical algorithms solving differential 
c  equations by Numerical Recipes. Needed by the routines above for the RGEs.
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
      implicit real*8(a-h,o-z)
      parameter (nmax=31)
      dimension y(n),dydx(n),yout(n),yt(nmax),dyt(nmax),dym(nmax)
      hh=h/2
      h6=h/6
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      CALL SU_DERIVS(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      CALL SU_DERIVS(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      CALL SU_DERIVS(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2*dym(i))
14    continue
      return
      end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Following are the main routine for the RGE evolution of parameter between
c  low and high energy scales. It returns a set of n mass and coupling
c parameters "y" at a specified scale exp(x2) when given at an initial
c scale exp(x2). Its uses the beta functions in the subroutine SU_DERIVS
c and solves the coupled differential equations with the SU_RKQC
c Runge-Kutta subroutine. Thus y(n) is a  vector containing all the n RG
c evolving parameters at  various possible scales depending on evolution
c stages. The parameters are  
c  y(1) = g1^2   U(1) gauge coupling squared
c  y(2) = g2^2   SU(2)_L gauge coupling squared
c  y(3) = g3^2   SU(3) gauge coupling squared
c  y(4) = Y_tau  tau lepton Yukawa coupling 
c  y(5) = Y_b    bottom  quark Yukawa coupling
c  y(6) = Y_top  top quark Yukawa coupling
c  y(7) = Ln(vu) logarithm of the vev vu 
c  y(8) = Ln(vd) logarithm of the vev vd 
c  y(9) = A_tau  trilinear coupling for stau 
c  y(10)= A_b    trilinear coupling for sbottom
c  y(11)= A_top  trilinear coupling for stop
c  y(12)= (m_phi_u)^2  scalar phi_u mass term squared
c  y(13)= (m_phi_d)^2  scalar phi_d mass term squared  
c  y(14)= mtaur^2 right-handed stau mass term squared
c  y(15)= msl^2   left-handed stau mass term squared
c  y(16)= mbr^2   right-handed sbottom mass term squared
c  y(17)= mtr^2   right-handed stop mass term squared 
c  y(18)= msq^2   left-handed stop mass term squared
c  y(19)= B       the (dimensionful) bilinear parameter B 
c  y(20)= Ln(|M1|) logarithm of the bino mass term
c  y(21)= Ln(|M2|) logarithm of the wino mass term
c  y(22)= Ln(|M3|) logarithm of the gluino mass term
c  y(23)= Ln(|mu|) logarithm of the |mu| parameter
c  y(24)= mer^2 right-handed selectron (smuon) mass term squared
c  y(25)= mel^2 left-handed selectron (smuon) mass term squared
c  y(26)= mdr^2 right-handed sdown (sstrange) mass term squared
c  y(27)= mur^2 right-handed sup (scharm) mass term squared
c  y(28)= muq^2 left-handed sup (scharm) mass term squared
c  y(29)= A_l   trilinear coupling for selectron (smuon) 
c  y(30)= A_d   trilinear coupling for sdown (sstrange)
c  y(21)= A_u   trilinear coupling for sup (scharm). 
c  Note that the number of running parameters consist of the 22 parameters
c  of the phenomenological MSSM; + the 3 gauge and the  3 Yukawa couplings, 
c  +3 parameters (vu, vd, B) which are in fact linearly dependent of others.
cc  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           SUBROUTINE SU_DERIV1(x,y,dydx)		  
c           SUBROUTINE SU_DERIV2(x,y,dydx)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  These are the derivatives of the RG running parameters y(xN), i.e the beta 
c  functions beta(y)=d(y)/dln(Q). The analytic expressions of the functions 
c  are taken from (up to some sign conventions which have been changed):
c  Castano, Ramond, Piard in Phys. Rev. D49 (1994) 4882, 
c  Barger, Berger, Ohmann in Phys.Rev. D49 (1994) 4908.
c  DERIV1 : includes only 1-loop RGE with full MSSM threshold. 
c  DERIV2 : includes 2-loop RGE for gauge, Yukawa cpls, gaugino masses, m_Hu,d
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                 SUBROUTINE SU_DERIV1(x,y,dydx)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                 
       implicit real*8(a-h,m,o-z)
       parameter (n=31)
       dimension y(n),dydx(n),ygut(n),yewsb(n)
       COMMON/SU_cte/nf,cpi,zm,wm,tbeta
       COMMON/SU_sthresh/rmtop,susym,egut
       COMMON/SU_fmasses/mtaudum,mbdum,mtop
       COMMON/SU_stepwi/wistep,h1,kpole
       COMMON/SU_stegut/ifirst,jfirst,ygut
       COMMON/SU_gunif/kunif
       COMMON/SU_sgnm123/sgnm1,sgnm2,sgnm3
       COMMON/SU_mssmsqua/msq,mtr,mbr,muq,mur,mdr
       COMMON/SU_treesfer/msbtr1,msbtr2,msttr1,msttr2
       COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
       COMMON/SU_tbewsb/vuewsb,vdewsb
       common/su_allewsb/yewsb
       COMMON/SU_renscale/scale
       real*8 nu,nd,ne,nn,nsq,n1sq,nsu,n1su,nsd,n1sd,nsl,n1sl,nse,n1se,
     . nhino,nh,nbino,nwino,ngino,nf,nhl,nhlino,nhh,nhhino
       data nn/3.d0/,nd/3.d0/,ne/3.d0/
      pi=4*datan(1.d0)
      Q = dexp(x)
c g1, g2 gauge unif.:
      if(kunif.ne.0.and.h1.gt.0.d0.and.Q.gt.1.d15) then
      if(y(1).ge.y(2)) then
       if(ifirst.eq.0) then
      egut=Q
c freeze out the gauge +Yukawa+vu,vd couplings after that
       do ii=1,31 
       ygut(ii)=y(ii)
       enddo
       endif
       ifirst=1      
       endif
       endif
c simple (unique scale) threshold in beta functions:
      st1= 1.d0                !  (full MSSM RGEs)
      st2=1.d0
      nu=2.d0 +st1
      nsq=3*st2
      n1sq=st2
      nsu=3*st2
      n1su=st2
      nsd=3*st2
      n1sd=st2
      nsl=3*st2
      n1sl=st2
      nse=3*st2
      n1se=st2      
      nhino=2.d0*st2
      nbino=st2
      nwino=st2
      ngino=st2
      nh=1.d0+st2
      nhl=st2
      nhh=st2
      nhlino=st2
      nhhino=st2
c       
c  coefficient of the beta functions for gauge couplings
       b10 = 2.d0/5*(17*nu/12+5*nd/12+5*ne/4+nn/4)
     . +nsq/30
     . +4*nsu/15+nsd/15+nsl/10+nse/5
     . +nhino/5 +nh/10
       b20 = -22.d0/3+(nu+nd)/2+(ne+nn)/6
     . +nsq/2 + nsl/6 + nhino/3
     . +nh/6 +4*nwino/3
       b30 = 2*(nu+nd)/3 +nsq/3+nsu/6+nsd/6+2*ngino -11.d0
c
c - gauge coupling beta functions (nb the variables are g^2):
c -y(1)--y(3): g1^2,g2^2,g3^2.

       dydx(1) = 2*cpi*b10*y(1)*y(1) 
       dydx(2) = 2*cpi*b20*y(2)*y(2)
       dydx(3) = 2*cpi*b30*y(3)*y(3)

c - Yukawa coupling beta function (only Ytau,Yb,Ytop included):
c -y(4)--y(6): Ytau,Yb,Ytop
       ytau2 =y(4)*y(4)
       yb2 = y(5)*y(5)
       ytop2 = y(6)*y(6)

      tbet=dexp(y(7)-y(8))
      cb2=1.d0/(1.d0+tbet*tbet)
      sb2=1.d0-cb2
c  
      Ytaubeta = 3*y(4)*ytau2 
     .  +(3*yb2+ytau2)*y(4) 
     . +(-3.d0/5*y(1)*(15.d0/4+3.d0/4-
     . (1.d0/4+1.d0 +1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3.d0*(2)/4))*y(4)

      Ybbeta = 3*y(5)*yb2 +
     . y(5)*ytop2+(3*yb2+ytau2)*y(5) 
     . +(-3.d0/5*y(1)*(5.d0/12+3.d0/4-
     . (1.d0/36+1.d0/9+1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3*(2.d0)/4)
     . -y(3)*(8.d0-4.d0*(2)/3))*y(5)

      Ytbeta =3*y(6)*ytop2 +
     . y(6)*yb2+3*y(6)*ytop2
     .+(-3.d0/5*y(1)*(17.d0/12+3.d0/4-
     . (1.d0/36 +4.d0/9 +1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3.d0*(2)/4)
     . -y(3)*(8.d0-4.d0*(2)/3))*y(6)

c  
       dydx(4) = cpi*y(4)*(4*ytau2 +3*yb2-9*y(1)/5.d0-3*y(2))

       dydx(5) = cpi*y(5)*(6*yb2 +ytau2 +ytop2
     .            -7*y(1)/15. -3*y(2) -16*y(3)/3.) 

       dydx(6) = cpi*y(6)*(6*ytop2 +yb2
     .            -13*y(1)/15.d0 -3*y(2) -16*y(3)/3.d0)

c - Higgs vev beta functions:
c - y(7), y(8) = Ln(vu), Ln(vd)

       dydx(7) = cpi*(3.d0/4*(y(1)/5 +y(2)) -3*ytop2)
       dydx(8) = cpi*(3.d0/4*(y(1)/5 +y(2)) -3*yb2-ytau2)

c - soft susy-breaking terms beta functions:
c - y(9)--y(11) : Atau, Ab, Atop

       dydx(9) =cpi*(8*ytau2*y(9) +6*yb2*y(10)
     . +6*(3*y(1)*sgnM1*dexp(y(20))/5 +sgnM2*y(2)*dexp(y(21))))

       dydx(10) =cpi*(12*y(10)*yb2 +2*y(9)*ytau2+2*y(11)*ytop2
     . +14*y(1)/15*sgnM1*dexp(y(20))+6*y(2)*sgnM2*dexp(y(21))
     . +32*y(3)/3*sgnM3*dexp(y(22)))

       dydx(11) =cpi*(12*y(11)*ytop2 +2*y(10)*yb2
     . +26*y(1)/15*sgnM1*dexp(y(20))+6*y(2)*sgnM2*dexp(y(21))
     . +32*y(3)/3*sgnM3*dexp(y(22)))

c - y(12)--y(13) : m^2(phi_u), m^2(phi_d)
       trym2 = y(18)-2*y(17)+y(16)-y(15)+y(14)+y(12)-y(13)
     .     +2*(y(28)-2*y(27)+y(26)-y(25)+y(24))
       dydx(12) =2*cpi*(3*ytop2*(y(12)+y(18)+y(17)+y(11)*y(11))
     . +3.d0/10*y(1)*trym2 -3*y(1)/5*dexp(2*y(20))
     . -3*y(2)*dexp(2*y(21)))
    
       dydx(13) = 2*cpi*(ytau2*(y(13)+y(15)+y(14)+y(9)*y(9))
     . +3*yb2*(y(13)+y(18)+y(16)+y(10)*y(10))
     . -3.d0/10*y(1)*trym2 -3*y(1)/5*dexp(2*y(20))
     . -3*y(2)*dexp(2*y(21)))

c - (1-loop) y(14)--y(19) : m^2_tau, m^2_L, m^2_b, m^2_top, m^2_Q, B

       dydx(14) = 2*cpi*(2*ytau2*(y(13)+y(14)+y(15)+y(9)*y(9))
     . +3*y(1)/5*trym2 -12*y(1)/5*dexp(2*y(20)))

       dydx(15) = 2*cpi*(ytau2*(y(13)+y(15)+y(14)+y(9)*y(9))
     . -3*y(1)/10*trym2 -3*y(1)/5*dexp(2*y(20))
     . -3*y(2)*dexp(2*y(21)))

       dydx(16) = 2*cpi*(2*yb2*(y(13)+y(16)+y(18)+y(10)*y(10))
     . +y(1)/5*trym2-4*y(1)/15*dexp(2*y(20))
     . -16*y(3)/3*dexp(2*y(22)))

       dydx(17) =2*cpi*(2*ytop2*(y(12)+y(17)+y(18)+y(11)*y(11))
     . -2*y(1)/5*trym2 -16*y(1)/15*dexp(2*y(20))
     . -16*y(3)/3*dexp(2*y(22)))

       dydx(18) =2*cpi*(ytop2*(y(12)+y(17)+y(18)+y(11)*y(11))
     . +yb2*(y(13)+y(18)+y(16)+y(10)*y(10))
     . +y(1)/10*trym2 -y(1)/15*dexp(2*y(20))-3*y(2)*dexp(2*y(21))
     . -16*y(3)/3*dexp(2*y(22)))

       dydx(19) = 2*cpi*(3*y(11)*ytop2 +3*y(10)*yb2 +y(9)*ytau2
     . +3*y(1)/5*sgnM1*dexp(y(20)) +3*y(2)*sgnM2*dexp(y(21)))

c - Gauginos masses beta functions:
c - y(20)--y(22) : Ln (M1,M2,M3)

       dydx(20) = -2*cpi*(-3.d0/5-nf)*y(1)

       dydx(21) = -2*cpi*(5.d0-nf)*y(2)

       dydx(22) = -2*cpi*(9.d0-nf)*y(3)

c -  the mu parameter:
c - y(23) = Ln mu

       dydx(23) = cpi*(3*ytop2 +3*yb2+ytau2 -3*y(1)/5-3*y(2))

c - y(24)--y(28) : 1st and 2d gen. sfermion mass^2 terms:
c   m^2_er, m^2_eL, m^2_dr, m^2_ur, m^2_uL
       dydx(24) = 2*cpi*(3*y(1)/5*trym2 -12*y(1)/5*dexp(2*y(20)))

       dydx(25) = 2*cpi*(-3*y(1)/10*trym2 -3*y(1)/5*dexp(2*y(20))     
     . -3*y(2)*dexp(2*y(21)))

       dydx(26) = 2*cpi*(
     . y(1)/5*trym2-4*y(1)/15*dexp(2*y(20))-16*y(3)/3
     . *dexp(2*y(22)))

       dydx(27) =2*cpi*(
     . -2*y(1)/5*trym2 -16*y(1)/15*dexp(2*y(20))
     . -16*y(3)/3*dexp(2*y(22)))

       dydx(28) =2*cpi*(
     . y(1)/10*trym2 -y(1)/15*dexp(2*y(20))-3*y(2)*dexp(2*y(21))
     . -16*y(3)/3*dexp(2*y(22)))

c - y(29)--y(31) : Ae (Anu), Ad (As), Au (Ac)

       dydx(29) =cpi*(2*ytau2*y(9) +6*yb2*y(10)
     . +6*(3*y(1)*sgnM1*dexp(y(20))/5 +sgnM2*y(2)*dexp(y(21))))

       dydx(30) =cpi*(2*ytau2*y(9) +6*yb2*y(10)
     . +14*y(1)/15*sgnM1*dexp(y(20))+6*y(2)*sgnM2*dexp(y(21))
     . +32*y(3)/3*sgnM3*dexp(y(22)))

       dydx(31) =cpi*(6*ytop2*y(11)
     . +26*y(1)/15*sgnM1*dexp(y(20))+6*y(2)*sgnM2*dexp(y(21))
     . +32*y(3)/3*sgnM3*dexp(y(22)))

c
       end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            SUBROUTINE SU_DERIV2(x,y,dydx)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  -main reference for two-loop MSSM RGE:
c   S.P. Martin and M.T. Vaughn, hep-ph/9311340, Phys.Rev. D50 (1994) 2282
       implicit real*8(a-h,m,n,o-z)
       integer*2 n
       parameter (n=31)
       dimension y(n),dydx(n),ygut(n),yewsb(n)
       COMMON/SU_cte/nf,cpi,zm,wm,tbeta
       COMMON/SU_sthresh/rmtop,susym,egut
       COMMON/SU_fmasses/mtaudum,mbdum,mtop
       COMMON/SU_stepwi/wistep,h1,kpole
       COMMON/SU_stegut/ifirst,jfirst,ygut
       COMMON/SU_gunif/kunif
       COMMON/SU_sgnm123/sgnm1,sgnm2,sgnm3
       COMMON/SU_mssmsqua/msq,mtr,mbr,muq,mur,mdr
       COMMON/SU_treesfer/msbtr1,msbtr2,msttr1,msttr2
       COMMON/SU_yukaewsb/ytau,yb,yt,alsewsb,g2ewsb,g1ewsb
       COMMON/SU_tbewsb/vuewsb,vdewsb
       common/su_allewsb/yewsb
       COMMON/SU_renscale/scale
       data nn/3.d0/,nd/3.d0/,ne/3.d0/
c
      pi=4*datan(1.d0)
      Q = dexp(x)
c g1, g2 gauge unif.:
      if(kunif.ne.0.and.h1.gt.0.d0.and.Q.gt.1.d15) then
      if(y(1).ge.y(2)) then
       if(ifirst.eq.0) then
       egut=Q
c freeze out the gauge +Yukawa+vu,vd couplings at GUT scale:
       do ii=1,31 
       ygut(ii)=y(ii)
       enddo
       endif
      ifirst= 1      
      endif
      endif

      st1=1.d0                  ! (full MSSM RGEs)
      st2=1.d0                  
      nu=2.d0 +st1
      nsq=3*st2
      n1sq=st2
      nsu=3*st2
      n1su=st2
      nsd=3*st2
      n1sd=st2
      nsl=3*st2
      n1sl=st2
      nse=3*st2
      n1se=st2      
      nhino=2.d0*st2
      nbino=st2
      nwino=st2
      ngino=st2
      nh=1.d0+st2
      nhl=st2
      nhh=st2
      nhlino=st2
      nhhino=st2
c
       b10 = 2.d0/5*(17*nu/12+5*nd/12+5*ne/4+nn/4)
     . +nsq/30
     . +4*nsu/15+nsd/15+nsl/10+nse/5
     . +nhino/5 +nh/10
       b20 = -22.d0/3+(nu+nd)/2+(ne+nn)/6
     . +nsq/2 + nsl/6 + nhino/3
     . +nh/6 +4*nwino/3
       b30 = 2*(nu+nd)/3 +nsq/3+nsu/6+nsd/6+2*ngino -11.d0
c - 2-loop gauge coupling beta functions (nb the variables are g^2):
c -y(1)--y(3): g1^2,g2^2,g3^2.
       ytau2 =y(4)*y(4)
       yb2 = y(5)*y(5)
       ytop2 = y(6)*y(6)
       mm1 = sgnM1*dexp(y(20))
       mm2 = sgnM2*dexp(y(21))
       mm3 = sgnM3*dexp(y(22))
c
       dydx(1) = 2*cpi*b10*y(1)*y(1) 
     . +2*cpi*cpi*y(1)*y(1)*((19*nf/15+9.d0/25)*y(1)
     . +(3*nf/5+9.d0/5)
     . *y(2)+44*nf/15*y(3)-26*ytop2/5-14*yb2/5-18*ytau2/5)

       dydx(2) = 2*cpi*b20*y(2)*y(2)
     . +2*cpi*cpi*y(2)*y(2)*((nf/5+3.d0/5)*y(1)
     . +(7*nf-17.d0)*y(2)
     . +4*nf*y(3) -6*ytop2 -6*yb2-2*ytau2 )

       dydx(3) = 2*cpi*b30*y(3)*y(3)
     . +2*cpi*cpi*y(3)*y(3)*(11*nf/30*y(1)+3*nf/2*y(2)
     . +(34*nf/3-54.d0)
     . *y(3) -4*ytop2 -4*yb2 )
c - 2-loop Yukawa coupling beta function (only Ytau,Yb,Ytop included):
c -y(4)--y(6): Ytau,Yb,Ytop
c 
      tbet=dexp(y(7)-y(8))
      cb2=1.d0/(1.d0+tbet*tbet)
      sb2=1.d0-cb2
c  
      Ytaubeta = 3*y(4)*ytau2 
     .  +(3*yb2+ytau2)*y(4) 
     . +(-3.d0/5*y(1)*(15.d0/4+3.d0/4-
     . (1.d0/4+1.d0 +1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3.d0*(2)/4))*y(4)

      Ybbeta = 3*y(5)*yb2 +
     . y(5)*ytop2+(3*yb2+ytau2)*y(5) 
     . +(-3.d0/5*y(1)*(5.d0/12+3.d0/4-
     . (1.d0/36+1.d0/9+1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3*(2.d0)/4)
     . -y(3)*(8.d0-4.d0*(2)/3))*y(5)

      Ytbeta =3*y(6)*ytop2 +
     . y(6)*yb2+3*y(6)*ytop2
     .+(-3.d0/5*y(1)*(17.d0/12+3.d0/4-
     . (1.d0/36 +4.d0/9 +1.d0/4))
     . -y(2)*(9.d0/4+9.d0/4-3.d0*(2)/4)
     . -y(3)*(8.d0-4.d0*(2)/3))*y(6)

c
       dydx(4) =cpi*(Ytaubeta
     . +cpi*y(4)*(-10*ytau2*ytau2-9*yb2*yb2 -9*yb2*ytau2-3*yb2*ytop2
     . +(6*y(2)+6*y(1)/5)*ytau2 +(-2*y(1)/5+16*y(3))*yb2
     . +(9*nf/5+27.d0/10)*y(1)*y(1)+(3*nf-21.d0/2)*y(2)*y(2)
     . +9*y(1)*y(2)/5 )  ) 


        dydx(5) = cpi*(Ybbeta
     . +cpi*y(5)*(-22*yb2*yb2-5*ytop2*ytop2-5*yb2*ytop2-3*yb2*ytau2 
     . -3*ytau2*ytau2 
     . +4*y(1)/5*ytop2+(2*y(1)/5+6*y(2)+16*y(3))*yb2
     . +6*y(1)/5*ytau2 
     . +(7*nf/15+7.d0/18)*y(1)*y(1)+(3*nf-21.d0/2)*y(2)*y(2)
     . +(16*nf/3-304.d0/9)*y(3)*y(3)+y(1)*y(2)+8*y(1)*y(3)/9
     . +8*y(2)*y(3) )  )

        dydx(6) = cpi*(Ytbeta
     . +cpi*y(6)*(-22*ytop2*ytop2-5*yb2*yb2-5*yb2*ytop2-yb2*ytau2
     . +(6*y(1)/5+6*y(2)+16*y(3))*ytop2+2*y(1)/5*yb2 
     . +(13*nf/15+403.d0/450)*y(1)*y(1)+(3*nf-21.d0/2)
     . *y(2)*y(2)
     . +(16*nf/3-304.d0/9)*y(3)*y(3) +y(1)*y(2)+136.d0/45
     . *y(1)*y(3)
     . +8*y(2)*y(3) ) )

c - (2-loop) Higgs vev beta functions:
c - y(7), y(8) = Ln vu, Ln vd

       dydx(7) = cpi*(3.d0/4*(y(1)/5 +y(2)) -3*ytop2)
     . +cpi*cpi*(9*ytop2*ytop2/4+9*ytop2*yb2/4
     . -(19*y(1)/10+9*y(2)/2+20*y(3))*ytop2
     . -(279.d0/800+1803*nf/3200)*y(1)*y(1)
     . -(207.d0/32+357*nf/128)*y(2)*y(2)
     . -(27.d0/80+9*nf/160)*y(1)*y(2))

       dydx(8) = cpi*(3.d0/4*(y(1)/5 +y(2)) -3*yb2-ytau2)
     . +cpi*cpi*(9*yb2*yb2/4+9*yb2*ytop2/4+3*ytau2*ytau2/4
     . -(2*y(1)/5+9*y(2)/2+20*y(3))*yb2
     . -(9*y(1)/5+3*y(2)/2)*ytau2
     . -(279.d0/800+1803*nf/3200)*y(1)*y(1)
     . -(207.d0/32+357*nf/128)*y(2)*y(2)
     . -(27.d0/80+9*nf/160)*y(1)*y(2))

c - (2-loop) soft susy-breaking terms beta functions:
c - y(9)--y(11) : Atau, Ab, Atop
       
        dhtau1loop = y(9)*(3*yb2+6*ytau2-3*y(2) -9*y(1)/5)
     . +6*y(10)*yb2 +6*y(9)*ytau2+6*mM2*y(2)
     . +18*y(1)/5*mM1

       dhtau2loop = y(9)*(-9*yb2*yb2-3*yb2*ytop2-14*ytau2*ytau2
     . -15*yb2*ytau2 +(16*y(3)-2*y(1)/5)*yb2 +6*y(1)/5*ytau2
     . +(12*y(2)-6*y(1)/5)*ytau2 
     . +15*y(2)*y(2)/2 +9*y(2)*y(1)/5 +27*y(1)*y(1)/2 )
     . -6*(6*y(10)*yb2*yb2 +(y(10)+y(11))*yb2*ytop2)
     . -36*y(9)*ytau2*ytau2 -yb2*ytau2*(12*y(9)+18*y(10))
     . +(32*y(3)-4*y(1)/5)*yb2*y(10) +12*y(1)/5*ytau2*y(9)
     . +(6*y(2)+6*y(1)/5)*ytau2*y(9)
     . -(32*y(3)*mm3-4*y(1)/5*mm1)*yb2 -12*y(1)/5*mm1*ytau2
     . -12*y(2)*mm2*ytau2 
     . -30*y(2)*y(2)*mm2 -18*y(2)*y(1)/5*(mm1+mm2) -54*y(1)*y(1)*mm1

       dydx(9) =cpi*dhtau1loop +cpi*cpi*dhtau2loop -y(9)/y(4)*dydx(4)


       dhb1loop= y(10)*(8*yb2+ytau2+ytop2 -16*y(3)/3-3*y(2)-7*y(1)/15)
     . +10*y(10)*yb2 +2*y(9)*ytau2 +2*y(11)*ytop2
     . +14*y(1)/15*mM1+6*y(2)*mM2
     . +32*y(3)/3*mM3

       dhb2loop= y(10)*(-30*yb2*yb2 -5*ytop2*ytop2 -7*yb2*ytop2
     . -5*yb2*ytau2 -3*ytau2*ytau2
     .  +(16*y(3)-2*y(1)/5)*yb2 +6*y(1)/5*ytau2
     . +4*y(1)/5*ytop2 +(12*y(2)+6*y(1)/5)*yb2 
     . -16*y(3)*y(3)/9 +8*y(3)*y(2)/9 +15*y(2)*y(2)/2 +y(1)*y(2)
     . +287*y(1)*y(1)/90 )
     . -80*y(10)*yb2*yb2 -20*y(11)*ytop2*ytop2
     . -(8*y(10)+10*y(11))*yb2*ytop2
     . -12*y(9)*ytau2*ytau2 -yb2*ytau2*(6*y(9)+4*y(10))
     . +(32*y(3)-4*y(1)/5)*yb2*y(10) +12*y(1)/5*ytau2*y(9)
     . +8*y(1)/5*y(11)*ytop2 +(6*y(2)+6*y(1)/5)*yb2*y(10)
     . -(32*y(3)*mm3-4*y(1)/5*mm1)*yb2 -12*y(1)/5*mm1*ytau2
     . -(12*y(2)*mm2+8*y(1)/5*mm1)*yb2 -8*y(1)/5*mm1*ytop2 
     . +64*y(3)*y(3)/9*mm3 -16*y(3)*y(2)*(mm3+mm2)
     . -16*y(3)*y(1)/9*(mm3+mm1)  
     . -30*y(2)*y(2)*mm2 -2*y(2)*y(1)*(mm1+mm2) -574*y(1)*y(1)/45*mm1

       dydx(10) =cpi*dhb1loop +cpi*cpi*dhb2loop -y(10)/y(5)*dydx(5)



      dht1loop= y(11)*(8*ytop2+yb2 -16*y(3)/3-3*y(2)-13*y(1)/15)
     . +10*y(11)*ytop2 +2*y(10)*yb2
     . +26*y(1)/15*mM1+6*y(2)*mM2
     . +32*y(3)/3*mM3

      dht2loop= y(11)*(-5*yb2*yb2 -30*ytop2*ytop2 -7*yb2*ytop2
     .  -yb2*ytau2
     .  +(16*y(3)+4*y(1)/5)*ytop2 +12*y(2)*ytop2 +2*y(1)/5*yb2
     . -16*y(3)*y(3)/9 +8*y(3)*y(2) +15*y(2)*y(2)/2 +136*y(3)*y(1)/45
     . + y(2)*y(1) +2743*y(1)*y(1)/450 )
     . -80*y(11)*ytop2*ytop2 -20*y(10)*yb2*yb2
     . -(8*y(11)+10*y(10))*yb2*ytop2
     . -2*yb2*ytau2*(y(9)+y(10))
     . +(32*y(3)+8*y(1)/5)*ytop2*y(11) +4*y(1)/5*yb2*y(10)
     . +(6*y(2)+6*y(1)/5)*ytop2*y(11)
     . -(32*y(3)*mm3+8*y(1)/5*mm1)*ytop2 -4*y(1)/5*mm1*yb2
     . -(12*y(2)*mm2+4*y(1)/5*mm1)*ytop2  
     . +64*y(3)*y(3)/9*mm3 -16*y(3)*y(2)*(mm3+mm2)
     . -272*y(3)*y(1)/45*(mm3+mm1)  
     . -30*y(2)*y(2)*mm2 -2*y(2)*y(1)*(mm1+mm2) -5486*y(1)*y(1)/225*mm1

       dydx(11) =cpi*dht1loop +cpi*cpi*dht2loop -y(11)/y(6)*dydx(6)

c - (1-loop) y(12)--y(13) : m^2(phi_u), m^2(phi_d)
       trym2 = y(18)-2*y(17)+y(16)-y(15)+y(14)+y(12)-y(13)
     .     +2*(y(28)-2*y(27)+y(26)-y(25)+y(24))
       
       dydx(12) =2*cpi*(3*ytop2*(y(12)+y(18)+y(17)+y(11)*y(11))
     . +3.d0/10*y(1)*trym2 -3*y(1)/5*mm1**2
     . -3*y(2)*mm2**2)
c 
       dydx(13) = 2*cpi*(ytau2*(y(13)+y(15)+y(14)+y(9)*y(9))
     . +3*yb2*(y(13)+y(18)+y(16)+y(10)*y(10))
     . -3.d0/10*y(1)*trym2 -3*y(1)/5*mm1**2
     . -3*y(2)*mm2**2)

c     two-loop part 

       sumt = y(12)+y(18)+y(17)+y(11)*y(11)
       sumb = y(13)+y(18)+y(16)+y(10)*y(10)
       sumtau = y(13)+y(15)+y(14)+y(9)*y(9)

       trmQ = 2*y(28)+y(18)
       trmU = 2*y(27)+y(17)
       trmD = 2*y(26)+y(16)
       trmL = 2*y(25)+y(15)
       trmE = 2*y(24)+y(14)
       
       curlySp = -ytop2*(3*y(12)+y(18) -4*y(17))
     $      +yb2*(3*y(13)-y(18)-2*y(16))
     $      +ytau2*(y(13)+y(15)-2*y(14))
     $      +(1.5*y(2)+0.3*y(1))*(y(12)-y(13)-trmL)
     $      +(8.d0/3*y(3)+1.5*y(2)+1.d0/30*y(1))*trmQ
     $      -(16.d0/3*y(3)+16.d0/15*y(1))*trmU
     $      +(8.d0/3*y(3)+2.d0/15*y(1))*trmD + 1.2*y(1)*trmE
        
       sig1 = y(1)/5*(3*(y(12)+y(13))+trmQ+3*trmL+8*trmU
     $      + 2*trmD+6*trmE)
       sig2 = y(2)*(y(12)+y(13)+3*trmQ+trmL)
       sig3 = y(3)*(2*trmQ +trmU+trmD)
c
       dydx(12) = dydx(12)+cpi**2*
     $      (-6d0*(6*ytop2**2*(sumt+y(11)*y(11))
     $      +(sumt+sumb+2*y(11)*y(10))*ytop2*yb2)
     $      +32*y(3)*ytop2*(sumt+2d0*mm3**2-2*y(11)*mm3)
     $      +1.6*y(1)*ytop2*(sumt+2d0*mm1**2-2*y(11)*mm1)
     $      +1.2*y(1)*curlySp+33*y(2)**2*mm2**2
     $      +3.6*y(1)*y(2)*(mm2**2+mm1**2+mm2*mm1)
     $      +621.d0/25*y(1)**2*mm1**2+3*y(2)*sig2+0.6*y(1)*sig1)
       
       dydx(13) = dydx(13)+cpi**2*
     $      (-6d0*(6*yb2**2*(sumb+y(10)*y(10))
     $      +(sumt+sumb+2d0*y(11)*y(10))*ytop2*yb2
     $      +2*(sumtau+y(9)*y(9))*ytau2**2)
     $      +32*y(3)*yb2*(sumb+2*mm3**2-2*y(10)*mm3)
     $      -0.8*y(1)*yb2*(sumb+2*mm1**2-2*y(10)*mm1)
     $      +2.4*y(1)*ytau2*(sumtau+2d0*mm1**2-2*y(9)*mm1)
     $      -1.2*y(1)*curlySp+33*y(2)**2*mm2**2
     $      +3.6*y(1)*y(2)*(mm2**2+mm1**2+mm2*mm1)
     $      +621.d0/25*y(1)**2*mm1**2+3*y(2)*sig2+0.6*y(1)*sig1)

c - (2-loop) y(14)--y(19) : m^2_tau, m^2_L, m^2_b, m^2_top, m^2_Q, B
c 
       dydx(14) = 2*cpi*(2*ytau2*(y(13)+y(14)+y(15)+y(9)*y(9))
     . +3*y(1)/5*trym2 -12*y(1)/5*mm1**2)

       dydx(14) = dydx(14) +cpi**2*
     .      (-16*ytau2**2*(sumtau+y(9)*y(9))
     .      -ytau2*yb2*(12*(sumtau+sumb)+8*y(9)*y(10))
     .      +2*sumtau*ytau2*(6*y(2)-6*y(1)/5)
     .      +12*y(2)*ytau2*(2*mm2**2-2*y(9)*mm2)
     .      -12*y(1)/5*ytau2*(2*mm1**2-2*y(9)*mm1)
     .      +12*y(1)/5*curlySp +2808d0/25*y(1)**2*mm1**2 
     .      +12*y(1)/5*sig1)

       dydx(15) = 2*cpi*(ytau2*(y(13)+y(15)+y(14)+y(9)*y(9))
     . -3*y(1)/10*trym2 -3*y(1)/5*mm1**2
     . -3*y(2)*mm2**2)

       dydx(15) = dydx(15) +cpi**2*
     .      (-12*ytau2**2*(sumtau+y(9)*y(9))
     .      -6*ytau2*yb2*(sumtau+sumb+2*y(9)*y(10))
     .      + 6*y(1)/5*ytau2*(2*sumtau-4*mm1*y(9)+4*mm1**2)
     .      -6*y(1)/5*curlySp +33*y(2)**2*mm2**2
     .      +18*y(1)*y(2)/5*(mm2**2+mm1**2+mm1*mm2)
     .       +621d0/25*y(1)**2*mm1**2 
     .      +3*y(1)/5*sig1 +3*y(2)*sig2)

       dydx(16) = 2*cpi*(2*yb2*(y(13)+y(16)+y(18)+y(10)*y(10))
     . +y(1)/5*trym2-4*y(1)/15*mm1**2
     . -16*y(3)/3*mm3**2)

       dydx(16) = dydx(16) +cpi**2*
     .      (-32*yb2**2*(sumb+y(10)*y(10))
     .      -4*ytop2*yb2*(sumt+sumb+2*y(11)*y(10))
     .      -4*ytau2*yb2*(sumtau+sumb+2*y(9)*y(10))
     .      +2*(6*y(2)+2*y(1)/5)*yb2*sumb 
     .      +12*y(2)*yb2*(2*mm2**2-2*y(10)*mm2)
     .      +4*y(1)/5*yb2*(2*mm1**2-2*y(10)*mm1)
     .      +4*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +128*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +808d0/75*y(1)**2*mm1**2 
     .      +4*y(1)/15*sig1 +16*y(3)/3*sig3)

       dydx(17) =2*cpi*(2*ytop2*(y(12)+y(17)+y(18)+y(11)*y(11))
     . -2*y(1)/5*trym2 -16*y(1)/15*mm1**2
     . -16*y(3)/3*mm3**2)

       dydx(17) = dydx(17) +cpi**2*
     .      (-32*ytop2**2*(sumt+y(11)*y(11))
     .      -4*ytop2*yb2*(sumt+sumb+2*y(11)*y(10))
     .      +2*(6*y(2)-2*y(1)/5)*ytop2*sumt 
     .      +12*y(2)*ytop2*(2*mm2**2-2*y(11)*mm2)
     .      -4*y(1)/5*ytop2*(2*mm1**2-2*y(11)*mm1)
     .      -8*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +512*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +3424d0/75*y(1)**2*mm1**2 
     .      +16*y(1)/15*sig1 +16*y(3)/3*sig3)

       dydx(18) =2*cpi*(ytop2*(y(12)+y(17)+y(18)+y(11)*y(11))
     . +yb2*(y(13)+y(18)+y(16)+y(10)*y(10))
     . +y(1)/10*trym2 -y(1)/15*mm1**2 -3*y(2)*mm2**2
     . -16*y(3)/3*mm3**2 )

       dydx(18) = dydx(18) +cpi**2*
     .      (-20*ytop2**2*(sumt+y(11)*y(11))
     .      -20*yb2**2*(sumb+y(10)*y(10))
     .      -2*ytau2*yb2*(sumtau+sumb+y(9)*y(10))
     .      +2*y(1)/5*(4*ytop2*(sumt-2*mm1*y(11)+2*mm1**2)
     .                 +2*yb2*(sumb-2*mm1*y(10)+2*mm1**2) )
     .      +2*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +32*y(2)*y(3)*(mm3**2+mm2**2+mm2*mm3)
     .      +32*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +33*y(2)**2*mm2**2 
     .      +2*y(1)*y(2)/5*(mm2**2+mm1**2+mm1*mm2)
     .       +199*y(1)**2/75*mm1**2 
     .      +y(1)/15*sig1 +16*y(3)/3*sig3 +3*y(2)*sig2)

       dydx(19) = 2*cpi*(3*y(11)*ytop2 +3*y(10)*yb2 +y(9)*ytau2
     . +3*y(1)/5*mM1 +3*y(2)*mM2)
c
      betB2 = 
     . -12*(3*y(11)*ytop2**2+3*y(10)*yb2**2+y(11)*yb2*ytop2+
     . y(10)*ytop2*yb2 +y(9)*ytau2**2)
     . +(32*y(3)+8*y(1)/5)*y(11)*ytop2+(32*y(3)-4*y(1)/5)*y(10)*yb2
     . +12*y(1)/5*y(9)*ytau2
     . -(32*y(3)*mm3+8*y(1)*mm1/5)*ytop2 
     . -(32*y(3)*mm3-4*y(1)*mm1/5)*yb2 -12*y(1)/5*mm1*ytau2
     . -30*y(2)**2*mm2 -18*y(1)*y(2)/5*(mm1+mm2)-414*y(1)**2/25*mm1
c
      dydx(19)= dydx(19) +cpi**2*betB2 
c
c - Gauginos masses beta functions (includes two-loop):
c - y(20)--y(22) : Ln (M1,M2,M3)
       dydx(20) = -2*cpi*(-3.d0/5-nf)*y(1) 
     . +2*cpi*cpi*y(1)*((19*nf/15+9.d0/25)*y(1)*(1.d0+1.d0)
     . +(3*nf/5+9.d0/5)*y(2)*(1.d0+Mm2/Mm1) 
     . +44*nf/15*y(3)*(1.d0+Mm3/Mm1)
     . -26*ytop2*(1.d0-y(11)/mm1)/5
     . -14*yb2*(1.d0-y(10)/mm1)/5
     . -18*ytau2*(1.d0-y(9)/mm1)/5)

       dydx(21) = -2*cpi*(5.d0-nf)*y(2)
     . +2*cpi*cpi*y(2)*((nf/5+3.d0/5)*y(1)*(1.d0+mm1/mm2)
     . +(7*nf-17.d0)*y(2)*(1.d0+1.d0)
     . +4*nf*y(3)*(1.d0+mm3/mm2)
     . -6*ytop2*(1.d0-y(11)/mm2)
     . -6*yb2*(1.d0 -y(10)/mm2)
     . -2*ytau2*(1.d0 -y(9)/mm2) )

       dydx(22) = -2*cpi*(9.d0-nf)*y(3)
     . +2*cpi*cpi*y(3)*(11*nf/30*y(1)*(1.d0+mm1/mm3)
     . +3*nf/2*y(2)*(1.d0+mm2/mm3)
     . +(34*nf/3-54.d0)*y(3)*(1.d0+1.d0) 
     . -4*ytop2*(1.d0 -y(11)/mm3)
     . -4*yb2*(1.d0 -y(10)/mm3)  )

c - the mu parameter:
c - y(23) = Ln mu

       dydx(23) = cpi*(3*ytop2 +3*yb2+ytau2 -3*y(1)/5-3*y(2))

c     two-loop part

       dydx(23) = dydx(23)+ cpi**2*(
     $      -3*(3*ytop2**2+3*yb2**2+2*ytop2*yb2+ytau2**2)
     $      +(16*y(3)+4.d0/5.*y(1))*ytop2
     $      +(16*y(3)-2.d0/5.*y(1))*yb2+6.d0/5.*y(1)*ytau2
     $      +7.5*y(2)**2+1.8*y(1)*y(2)+207.d0/50.*y(1)**2)


c - (2-loop) y(24)--y(28) : 1st and 2d gen. sfermion mass^2 terms:
c   m^2_er, m^2_eL, m^2_dr, m^2_ur, m^2_uL

       dydx(24) = 2*cpi*(3*y(1)/5*trym2 -12*y(1)/5*mm1**2)

       dydx(24) = dydx(24) +cpi**2*
     .      (12*y(1)/5*curlySp +2808d0/25*y(1)**2*mm1**2 
     .      +12*y(1)/5*sig1)

       dydx(25) = 2*cpi*(-3*y(1)/10*trym2 -3*y(1)/5*mm1**2     
     . -3*y(2)*mm2**2)

       dydx(25) = dydx(25) +cpi**2*
     .      (-6*y(1)/5*curlySp +33*y(2)**2*mm2**2
     .      +18*y(1)*y(2)/5*(mm2**2+mm1**2+mm1*mm2)
     .       +621d0/25*y(1)**2*mm1**2 
     .      +3*y(1)/5*sig1 +3*y(2)*sig2)

       dydx(26) = 2*cpi*(
     . y(1)/5*trym2-4*y(1)/15*mm1**2 -16*y(3)/3*mm3**2 )

       dydx(26) = dydx(26) +cpi**2*
     .      (4*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +128*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +808d0/75*y(1)**2*mm1**2 
     .      +4*y(1)/15*sig1 +16*y(3)/3*sig3)

       dydx(27) =2*cpi*(
     . -2*y(1)/5*trym2 -16*y(1)/15*mm1**2 -16*y(3)/3*mm3**2)

       dydx(27) = dydx(27) +cpi**2*
     .      (-8*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +512*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +3424d0/75*y(1)**2*mm1**2 
     .      +16*y(1)/15*sig1 +16*y(3)/3*sig3)

       dydx(28) =2*cpi*(
     . y(1)/10*trym2 -y(1)/15*mm1**2 -3*y(2)*mm2**2
     . -16*y(3)/3*mm3**2 )

       dydx(28) = dydx(28) +cpi**2*
     .      (2*y(1)/5*curlySp -128*y(3)**2/3*mm3**2
     .      +32*y(2)*y(3)*(mm3**2+mm2**2+mm2*mm3)
     .      +32*y(1)*y(3)/45*(mm3**2+mm1**2+mm1*mm3)
     .       +33*y(2)**2*mm2**2 
     .      +2*y(1)*y(2)/5*(mm2**2+mm1**2+mm1*mm2)
     .       +199*y(1)**2/75*mm1**2 
     .      +y(1)/15*sig1 +16*y(3)/3*sig3 +3*y(2)*sig2)

c - (2-loop) y(29)--y(31) : Ae (Anu), Ad (As), Au (Ac)

        dhe1loop = y(29)*(ytau2 +3*yb2 -3*y(2) -9*y(1)/5)
     . +6*y(10)*yb2 +2*y(9)*ytau2 +6*mM2*y(2)
     . +18*y(1)/5*mM1

       dhe2loop = y(29)*(-9*yb2*yb2-3*yb2*ytop2-3*ytau2*ytau2
     .  +(16*y(3)-2*y(1)/5)*yb2 +6*y(1)/5*ytau2 
     . +15*y(2)*y(2)/2 +9*y(2)*y(1)/5 +27*y(1)*y(1)/2 )
     . -6*(6*y(10)*yb2*yb2 +(y(10)+y(11))*yb2*ytop2)
     . -12*y(9)*ytau2*ytau2
     . +(32*y(3)-4*y(1)/5)*yb2*y(10) +12*y(1)/5*ytau2*y(9)
     . -(32*y(3)*mm3-4*y(1)/5*mm1)*yb2 -12*y(1)/5*mm1*ytau2 
     . -30*y(2)*y(2)*mm2 -18*y(2)*y(1)/5*(mm1+mm2) -54*y(1)*y(1)*mm1

      dyovery4 = cpi*( ytau2 +3*yb2  -9*y(1)/5  -3*y(2))
     . +cpi**2*(
     . -3*ytau2*ytau2-9*yb2*yb2 -3*yb2*ytop2
     . +6*y(1)/5*ytau2 +(-2*y(1)/5+16*y(3))*yb2
     . +(9*nf/5+27.d0/10)*y(1)*y(1)+(3*nf-21.d0/2)*y(2)*y(2)
     . +9*y(1)*y(2)/5 )   

       dydx(29) =cpi*dhe1loop +cpi*cpi*dhe2loop -y(29)*dyovery4

       dhd1loop= y(30)*(3*yb2+ytau2-16*y(3)/3-3*y(2)-7*y(1)/15)
     . +6*y(10)*yb2 +2*y(9)*ytau2 
     . +14*y(1)/15*mM1+6*y(2)*mM2
     . +32*y(3)/3*mM3

       dhd2loop= y(30)*(-9*yb2*yb2  -3*yb2*ytop2
     . -3*ytau2*ytau2
     .  +(16*y(3)-2*y(1)/5)*yb2 +6*y(1)/5*ytau2 
     . -16*y(3)*y(3)/9 +8*y(3)*y(2)/9 +15*y(2)*y(2)/2 +y(1)*y(2)
     . +287*y(1)*y(1)/90 )
     . -36*y(10)*yb2*yb2 -6*(y(10)+y(11))*yb2*ytop2 
     . -12*y(9)*ytau2*ytau2
     . +(32*y(3)-4*y(1)/5)*yb2*y(10) +12*y(1)/5*ytau2*y(9)
     . -(32*y(3)*mm3-4*y(1)/5*mm1)*yb2 -12*y(1)/5*mm1*ytau2 
     . +64*y(3)*y(3)/9*mm3 -16*y(3)*y(2)*(mm3+mm2)
     . -16*y(3)*y(1)/9*(mm3+mm1)  
     . -30*y(2)*y(2)*mm2 -2*y(2)*y(1)*(mm1+mm2) -574*y(1)*y(1)/45*mm1

      dyovery5 = cpi*( 3*yb2+ytau2
     . -7*y(1)/15 -3*y(2) -16*y(3)/3)
     . +cpi**2*( -9*yb2*yb2 -3*yb2*ytop2 -3*ytau2*ytau2 
     . +(-2*y(1)/5+16*y(3))*yb2 +6*y(1)/5*ytau2 
     . +(7*nf/15+7.d0/18)*y(1)*y(1)+(3*nf-21.d0/2)*y(2)*y(2)
     . +(16*nf/3-304.d0/9)*y(3)*y(3)+y(1)*y(2)+8*y(1)*y(3)/9
     . +8*y(2)*y(3)  )

       dydx(30) =cpi*dhd1loop +cpi*cpi*dhd2loop -y(30)*dyovery5

      dhu1loop= y(31)*(3*ytop2 -16*y(3)/3-3*y(2)-13*y(1)/15)
     . +6*y(11)*ytop2
     . +26*y(1)/15*mM1 +6*y(2)*mM2 +32*y(3)/3*mM3

      dhu2loop= y(31)*( -9*ytop2*ytop2 -3*yb2*ytop2
     .  +(16*y(3)+4*y(1)/5)*ytop2 
     . -16*y(3)*y(3)/9 +8*y(3)*y(2) +15*y(2)*y(2)/2 +136*y(3)*y(1)/45
     . + y(2)*y(1) +2743*y(1)*y(1)/450 )
     . -36*y(11)*ytop2*ytop2 
     . -6*(y(11)+y(10))*yb2*ytop2
     . +(32*y(3)+8*y(1)/5)*ytop2*y(11) 
     . -(32*y(3)*mm3+8*y(1)/5*mm1)*ytop2  
     . +64*y(3)*y(3)/9*mm3 -16*y(3)*y(2)*(mm3+mm2)
     . -272*y(3)*y(1)/45*(mm3+mm1)  
     . -30*y(2)*y(2)*mm2 -2*y(2)*y(1)*(mm1+mm2) -5486*y(1)*y(1)/225*mm1

      dyovery6 = cpi*( 3*ytop2
     . -13*y(1)/15 -3*y(2) -16*y(3)/3)
     . +cpi**2*( -9*ytop2*ytop2-3*yb2*ytop2
     . +(4*y(1)/5 +16*y(3))*ytop2
     . + (13*nf/15+403.d0/450)*y(1)*y(1)+(3*nf-21.d0/2)
     . *y(2)*y(2)
     . +(16*nf/3-304.d0/9)*y(3)*y(3) +y(1)*y(2)+136.d0/45
     . *y(1)*y(3)
     . +8*y(2)*y(3) )

       dydx(31) =cpi*dhu1loop +cpi*cpi*dhu2loop -y(31)*dyovery6
c
        end
c
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   ++++++++++++++ End of the routines for the RG evolution ++++++++
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
c  
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  These routines are for the QCD running of quark masses and couplings. 
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           SUBROUTINE ALSINI(ACC)
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine for initialization in the evaluation of the strong coupling 
c  constant alpha_s. It needs the two iteration functions to determine the
c  (improved) values of QCD scale Lambda for a given number of quark flavor
c  and masses, loop order, etc..: 
c          DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)
c          DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)
c  There are also two important functions for the calculation of the 
c  running of the QCD coupling at scale Q and perturbative order N:        
c          DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
c  and the running of the quark masses at scale Q and with NF quark flavors:
c          DOUBLE PRECISION FUNCTION RUNM(Q,NF)
c  These routines are borrowed from the program HDECAY version 2.2 
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ************************* FUNCTION RUNM ***************************
      DOUBLE PRECISION FUNCTION RUNM(Q,NF)
      IMPLICIT DOUBLE PRECISION (A-H,m,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/SU_ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/SU_fmasses/AMTAU,AMB,AMT
      COMMON/SU_RUN/XMSB,XMHAT,XKFAC
      COMMON/SU_QCDFLAG/NNLO,IDRflag
      common/su_mbmb/mbmb,imbmb  
      SAVE ISTRANGE
      B0(NF)=(33.D0-2.D0*NF)/12D0
      B1(NF) = (102D0-38D0/3D0*NF)/16D0
      B2(NF) = (2857D0/2D0-5033D0/18D0*NF+325D0/54D0*NF**2)/64D0
      G0(NF) = 1D0
      G1(NF) = (202D0/3D0-20D0/9D0*NF)/16D0
      G2(NF) = (1249D0-(2216D0/27D0+160D0/3D0*ZETA3)*NF
     .       - 140D0/81D0*NF**2)/64D0
      C1(NF) = G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF) = ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .       + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .       - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2D0
      TRAN(X,XK)=1D0+coeff1*alphas(X,2)/PI+XK*(alphas(X,2)/PI)**2
     .           +coeff3*(alphas(x,2)/pi)**3
c 3-loop coeff3 in M_pole/M_running relation added 10/12/03
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      nnlo =1  
c (always use NNLO now)
C     Define the light quark masses
      AMC=1.40D0
      AMSB=0.19d0
c      AMC=AMCA
      AMS=AMSB
c      
      PI=4D0*DATAN(1D0)
      ACC = 1.D-8
      if(idrflag.ne.1) then
      coeff1 = 4D0/3D0
      else
      coeff1 = 5D0/3D0
      endif
      if(nnlo.eq.0) then
      coeff3=0.d0
      else
      coeff3=101.45424d0
      endif
      AM(1) = 0.d0
      AM(2) = 0.d0
C--------------------------------------------
      IMSBAR = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(alphas(AMS,2)/PI,3)
     .            *CQ(alphas(1.D0,2)/PI,3)/TRAN(AMS,0.D0)
        DD = (XMSB-AMSB)/AMSB
        IF(DABS(DD).GE.ACC)THEN
         IF(DD.LE.0.D0)THEN
          AMSD = AM(3)
         ELSE

          AMSU = AM(3)
         ENDIF
         GOTO 123
        ENDIF
        ISTRANGE=1
       ENDIF
       AM(3) = AMSB
      ELSE
       AMS=AMSB
       AM(3) = AMS
      ENDIF
C--------------------------------------------
c-!! modifs jlk: to determine (perturbatively, at an order consistent
c with the pert. level used in RUNM) Mb(pole) from mb(mb)_MSbar input:
c    mbmb= mb(mb)_MSbar ; MBpole determined iteratively to acc. d-8
      if(imbmb.eq.0) then
c   imbmb is just a counter because this calculation is only needed once 
      do i=1,20
      if(i.eq.1) then
      mbsave=0.d0
      MBpole=mbmb
      endif
      if(nnlo.eq.0) then
      xkb=0.d0
      else
      xkb= 16.11d0 -1.04d0*(4.d0-(amsb+amc)/MBpole)
      endif
      if(i.ge.3) then 
      amba=mbpole
      call alsini(1.d-8)
      endif
      mbMBpole=mbmb*CQ(alphas(MBpole,2)/pi,4)/CQ(alphas(mbmb,2)/pi,4)
c  mbMBpole is mb(MBpole)
      MBpole= mbMBpole*tran(MBpole,xkb)
c tran(Q,xk) is the usual pert. relation between Mpole and mrun(Mpole),
c see its def. above
      if(dabs(1.d0-mbsave/MBpole).lt.1.d-8) goto 2
      mbsave=MBpole
      enddo
 2    AMB=MBpole
      imbmb=1
      endif
c rest of calculation follows as before: 
c---
      AM(3) = AMSB
      AM(4) = AMC
      AM(5) = AMB
      AM(6) = AMT
      XK = 16.11D0
      DO 1 I=1,NF-1
       XK = XK - 1.04D0*(1.D0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB = AM(NF)/TRAN(AM(NF),0D0)
       XMHAT = XMSB/CQ(alphas(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0.d0
       XMHAT = 0.d0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(alphas(AM(4),2)/PI,3)/
     .                   CQ(alphas(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(alphas(AM(5),2)/PI,4)/
     .                   CQ(alphas(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(alphas(AM(6),2)/PI,5)/
     .                   CQ(alphas(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(alphas(AM(5),2)/PI,4)/
     .                   CQ(alphas(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(alphas(AM(6),2)/PI,5)/
     .                   CQ(alphas(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(alphas(AM(4),2)/PI,4)/
     .                   CQ(alphas(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(alphas(AM(6),2)/PI,5)/
     .                   CQ(alphas(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(alphas(AM(5),2)/PI,5)/
     .                   CQ(alphas(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(alphas(AM(4),2)/PI,4)/
     .                   CQ(alphas(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0=3
       Q0 = 1.D0
      ELSEIF(Q.LE.AMB)THEN
       N0=4
       Q0 = AMC
      ELSEIF(Q.LE.AMT)THEN
       N0=5
       Q0 = AMB
      ELSE
       N0=6
       Q0 = AMT
      ENDIF
      IF(NNLO.EQ.1.AND.NF.GT.3)THEN
       XKFAC = TRAN(AM(NF),0D0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC = 1.D0
      ENDIF
      runm = YMSB(N0)*CQ(alphas(Q,2)/PI,N0)/
     .               CQ(alphas(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END

C ************************* FUNCTION ALPHAS ***************************
      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/SU_ALSLAM/XLB1(6),XLB2(6)
      COMMON/SU_ALS/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
c
      PI=4.D0*DATAN(1.D0)
      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1      CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2      CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMB)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        alphas=ALS1(NF,Q)
      ELSE
        alphas=ALS2(NF,Q)
      ENDIF
      RETURN
      END

C ************************* FUNCTION ALSINI ***************************
      SUBROUTINE ALSINI(ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/SU_ALSLAM/XLB1(6),XLB2(6)
      COMMON/SU_ALS/XLAMBDA,AMC,AMB,AMT,N0

      PI=4.D0*DATAN(1.D0)
      XLB1(1)=0.D0
      XLB1(2)=0.D0
      XLB2(1)=0.D0
      XLB2(2)=0.D0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
      ENDIF
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE
      RETURN
      END
C ************************* FUNCTION XITER ***************************
      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2*(1.D0+DSQRT(1D0-4*B*DLOG(X)))
      PI=4*DATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER=XLB2
      RETURN
      END
C ************************* FUNCTION XITLA ***************************
      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)
C--ITERATION ROUTINE TO DETERMINE IMPROVED LAMBDA'S
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SU_PARAM/GF,alph,AMZ,AMW
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4*B*DLOG(X)))
      PI=4*DATAN(1.D0)
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2.D0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      XX=XIT(A,B,X)
      XLB=Q*DEXP(-XX/2.D0)
      Y1=ALP
      Y2=ALS2(NF,Q,XLB)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB
      RETURN
      END
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   ++++++++++++++ End of the routines for the QCD running  ++++++++++++++++
c  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine su_gminus2(mel,mer,Amu,mu,tb,u,v,z,mn,mc1,mc2,gmuon)
c------------------------------------------------------------------
c Calculates leading (chargino and neutralino loops) SUSY contributions
c to g_mu -2 
c  INPUT: MEL,MER, AL: relevant soft terms (i.e. 2d generation muon sector);
c  MU, tb=tan(beta); 
c  U,V,Z, mn, mc1,mc2: chargino and neutralino masses and mixing matrices.
c  OUPTUT:  gmuon, is a_mu = g_mu -2 in standard units
c------------------------------------------------------------------
       implicit real*8(a-h,m,o-z)
       dimension u(2,2),v(2,2),z(4,4),mc(2),mn(4),msl(2),sl(2),sr(2),
     .           xcl(2),xcr(2),xnl(4,2),xnr(4,2),anl(4,2),acl(2)
       COMMON/su_PARAM/GF,ALPH,MZ,MW
c
       fgm2a(x) =-(1.d0-6*x+3*x**2+2*x**3-6*x**2*dlog(x))
     .       /(1.d0-x)**4/6
c
       fgm2b(x) = (1.d0-x**2+2*x*dlog(x))/(1.d0-x)**3
c
       fgm2c(x) =(1.d0+1.5*x-3*x**2+0.5d0*x**3+3*x*dlog(x))
     .       /(1.d0-x)**4/3
c
       fgm2d(x) =-3*(1.d0-4.d0/3*x+x**2/3+2.d0/3*dlog(x))
     .       /(1.d0-x)**3
c
        ml=0.105658357d0
        mc(1)=mc1
        mc(2)=mc2
        b=datan(tb)
        cw=mw/mz
        sw=dsqrt(1.d0-cw**2)
        sw2=sw**2
        pi = 4*datan(1.d0)
C    
C     calculation of the slepton masses and mixing
      dt = dcos(2.d0*b)*mz*mz
      MSEL2 = MEL**2 + (-0.5D0+SW2)*dt
      MSER2 = MER**2 - SW2*dt
      MSNL2 = MEL**2 + 0.5D0*dt
      MLRE = Amu - MU*TB
      DELE = (MSEL2-MSER2)**2 + 4.d0*(ML*MLRE)**2
      MSE12 = ML**2 + 0.5D0*( MSEL2 + MSER2 + DSQRT(DELE) )
      MSE22 = ML**2 + 0.5D0*( MSEL2 + MSER2 - DSQRT(DELE) )
cc      MSL(1) = DSQRT(MSE12)
cc      MSL(2) = DSQRT(MSE22)
      MSN = DSQRT(MSNL2)      
      THEL = 0.5D0 * DATAN( 2.D0*ML*MLRE/(MSEL2-MSER2) )
      CCL = DCOS(THEL)
      SSL = DSIN(THEL) 
c def. of mass eingenvalues in terms of angle: 
      msl(1)= dsqrt(ccl**2*(ml**2+msel2)+ssl**2*(ml**2+mser2)
     .       +2*ccl*ssl*ml*mlre)
      msl(2)= dsqrt(ssl**2*(ml**2+msel2)+ccl**2*(ml**2+mser2)
     .       -2*ccl*ssl*ml*mlre)

C     Calculation of the chargino and neutralino couplings
      rt2 = dsqrt(2.d0)
      yuk = ml / ( rt2 * mw * sw * dcos(b) )
      xcl(1) = yuk*u(1,2)
      xcl(2) = yuk*u(2,2)
      xcr(1) = -v(1,1)/sw
      xcr(2) = -v(2,1)/sw
c 
c i.e. for small mixing angle, slepton_1 should be mostly slepton_L
        do ii = 1,4
      xnr(ii,1) = yuk*z(ii,3)*ccl + rt2/cw*z(ii,1)*ssl
      xnr(ii,2) = -yuk*z(ii,3)*ssl + rt2/cw*z(ii,1)*ccl
      gau = ( z(ii,1)/cw + z(ii,2)/sw ) / rt2
      xnl(ii,1) = -yuk*z(ii,3)*ssl + gau*ccl
      xnl(ii,2) = -yuk*z(ii,3)*ccl - gau*ssl
        enddo    

	do 13 i=1,4
	do 13 j=1,2
        anl(i,j)=ml/msl(j)**2*(xnl(i,j)*xnl(i,j)+xnr(i,j)*xnr(i,j))
     .            *fgm2a(mn(i)**2/msl(j)**2)
     .   +mn(i)/msl(j)**2*xnl(i,j)*xnr(i,j)*fgm2b(mn(i)**2/msl(j)**2)       
 13	continue
       anltot=anl(1,1)+anl(2,1)+anl(3,1)+anl(4,1)
     .       +anl(1,2)+anl(2,2)+anl(3,2)+anl(4,2)

        do 12 i=1,2
       acl(i)=ml/msn**2*(xcl(i)**2+xcr(i)**2)*fgm2c(mc(i)**2/msn**2)
     .        +mc(i)/msn**2*xcl(i)*xcr(i)*fgm2d(mc(i)**2/msn**2)
 12     continue
       acltot=acl(1)+acl(2)
       gmuon = ml/4.d0/pi/137.d0 * (anltot+acltot) 
C  *** Include leading-log 2-loop QED correction  ***
       MSUSY = .2D0*( MSL(1) + MSL(2) + MSN + mc(1) + mc(2) )
       GMUON = GMUON * ( 1.D0 - 4.D0/(PI*137.D0)*DLOG(MSUSY/ML) )
       return
       end
c--------------------------------------------------
       subroutine su_delrho(mt,gmst,gmsb,gmstau,msn,thetat,thetab,thel,
     . drho) 
c--------------------------------------------------
c   calculates leading one-loop SUSY delta_rho contributions of 3rd gen
c sfermions (plus leading two-loop QCD contributions) 
c  INPUT: MT, gmst(2), gmsb(2),gmstau(2),msn: top,stop,sbottom,
c  stau, stau neutrino masses and stop, sbottom, stau mixing angles
c  OUTPUT: drho = rho-1 
c--------------------------------------------------
       implicit real*8(a-h,m,o-z)
       dimension gmst(2),gmsb(2),gmstau(2)
       common/SU_PARAM/GF,alph,mz,mw
       su_fr(x,y) = x+y-2*x*y/(x-y)*dlog(x/y)
c
       pi = 4*datan(1.d0)
       ct=dcos(thetat)
       st=dsin(thetat)
       cb=dcos(thetab)
       sb=dsin(thetab)
       ctau =dcos(thel)
       stau =dsin(thel)
       cta2=ctau**2
       sta2=stau**2
       ct2=ct**2
       st2=st**2
       cb2=cb**2
       sb2=sb**2
       mt1=gmst(1)**2
       mt2=gmst(2)**2
       mb1=gmsb(1)**2
       mb2=gmsb(2)**2
       mta1=gmstau(1)**2
       mta2=gmstau(2)**2
c
       drhotb= (ct2*(cb2*su_fr(mt1,mb1)+sb2*su_fr(mt1,mb2)) +
     .       st2*(cb2*su_fr(mt2,mb1)+sb2*su_fr(mt2,mb2)) -
     .       ct2*st2*su_fr(mt1,mt2)-cb2*sb2*su_fr(mb1,mb2))
       drhotau= -cta2*sta2*su_fr(mta1,mta2)+cta2*su_fr(mta1,msn**2) +
     . sta2*su_fr(mta2,msn**2)
       drho = 3*drhotb*(1.d0 +2*0.12/3/pi*(1.d0+pi**2/3))+drhotau
       drho = GF/(8* pi**2* dsqrt(2.d0))*drho
       end
c-----------------------------------------------------
      subroutine su_finetune(mu,tb,mhd2,mhu2, 
     . czmu,czbmu,ctmu,ctbmu)
c--------------------------------------------------------------------------
c Calculates the degree of fine-tuning in a given model 
c  (at the moment with respect to MZ and Mtop only).
c input: mu,tbeta, mHd^2, mHu^2 (at the EWSB scale)
c output:   czmu,czbmu,ctmu,ctbmu are (dimensionless) measures of
c the degree of fine-tuning   on MU and B*MU with respect to MZ and Mtop,
c respectively. The larger those numbers (>>1), the more it is "fine-tuned"
c--------------------------------------------------------------------------
      implicit real*8(a-h,m,o-z)
       COMMON/su_PARAM/GF,ALPH,MZ,MW
c
      czmu = 2*mu**2/mz**2*(1.d0 + (tb**2+1.d0)/(tb**2-1.d0)**2*
     . 4*tb**2*(mhd2-mhu2)/((mhd2-mhu2)*(tb**2+1.d0)-
     . mz**2*(tb**2-1.d0)) )
c
       czbmu = 4*tb**2*(tb**2+1.d0)/(tb**2-1.d0)**3*(mhd2-mhu2)/mz**2
c
       ctmu = czmu/2 +2*mu**2/(mhd2+ mhu2+2*mu**2)/(tb**2-1.d0)
c
       ctbmu = czbmu/2 +1.d0/(1.d0-tb**2)
       end
c--------------------------------------------------------------------------
c -------------------------------------------------------------------- c
c ------ This and the following subroutines read in the spectrum ----- c
c ------ file given in the SUSY Les Houches Accord format        ----- c
c ------ hep-ph/0311123.                                         ----- c
c ------ Thanks to Tilman Plehn for the first version which has  ----- c
c ------ been expanded and changed here.                         ----- c
c ------Modified and adapted for SuSpect ver >=2.3 by J-L Kneur  ------c
c ------Updated for ver 2.4 (JLK) to correct some mismatch between 
c ------ SLHA and SuSpect for model choice input conventions
c -------------------------------------------------------------------- c

      subroutine SU_read_leshouches(input,ninlha,ichoice,imod)

      implicit double precision (a-h,m,o-z)
      double precision minval(1:20),extval(0:60),smval(1:20),
     .            massval(1:50),nmixval(4,4),umixval(2,2),vmixval(2,2),
     .                 stopmixval(2,2),sbotmixval(2,2),staumixval(2,2),
     .                 hmixval(1:10),gaugeval(1:3),msoftval(1:30),
     .                 auval(3,3),adval(3,3),aeval(3,3),yuval(3,3),
     .                 ydval(3,3),yeval(3,3)
       double precision nl,nq
      integer ichoice(1:11),check(1:24),check_final,imod(1:2)
      character line1*6,line2*100,
     .          spinfo1*100,spinfo2*100,modselval*100,mincom(1:20)*20,
     .          extcom(0:60)*20
      logical done,endfile

      COMMON/SU_leshouches1/spinfo1,spinfo2,modselval,mincom,extcom
      COMMON/SU_leshouches2/minval,extval,smval,massval,nmixval,umixval,
     .                      vmixval,stopmixval,sbotmixval,staumixval,
     .                      hmixval,gaugeval,msoftval,auval,adval,
     .                      aeval,yuval,ydval,yeval,alphaval,Qvalhmix,
     .                      Qvalgauge,Qvalmsoft,Qvalau,Qvalad,Qvalae,
     .                      Qvalyu,Qvalyd,Qvalye
       COMMON/SU_SMPAR/alfinv,sw2,alphas,mt,mb,mc,mtau
       COMMON/SU_param/gf,alpha,mz,mw
       COMMON/SU_RGSCAL/qewsb,ehigh,elow
       COMMON/SU_MSSMHPAR/mhu2,mhd2,ma,mu
       COMMON/SU_MSSMGPAR/m1,m2,m3 
       COMMON/SU_MSSMSLEP/msl,mtaur,mel,mer
       COMMON/SU_MSSMSQUA/msq,mtr,mbr,muq,mur,mdr
       COMMON/SU_ATRI3/atau,at,ab
       COMMON/SU_ATRI12/al,au,ad
       COMMON/SU_MSUGRA/m0,mhalf,a0
       COMMON/SU_RADEWSB/sgnmu0,tgbeta
       COMMON/SU_GMSB/mgmmess,mgmsusy,nl,nq
       COMMON/SU_AMSB/m32,am0,cq,cu,cd,cl,ce,chu,chd
       common/SU_slha_warn/smin_warn,extpar_warn,muma_warn,algo_warn
c
      if(input.eq.2) then
c = case where the slha format input file suspect2_lha.in is actually read:
      unlikely = -123456789d0    ! will be used for protection see later
      smin_warn=0d0
      extpar_warn=0d0
      muma_warn=0d0
c what follows is in case specific SuSpect block SU_ALGO undefined in file:
c then take defaut values of these algorithm control parameters
      algo_warn=-1d0
             ichoice(2)= 21     ! 2-loop RGE by defaut
	     ichoice(3)=1       ! GUT gauge coupling unif. request by defaut
	     ichoice(4)=2     ! RGE accuracy high by defaut
	     ichoice(5)=1     ! radiative EWSB by defaut
	     ichoice(6)=1     ! Mhu,Mhd input by defaut
	     ichoice(7)=2     ! rad. corr. for all sparticles by defaut
	     ichoice(8)=1    ! EWSB scale =sqrt(mst1*mst2) by defaut
	     ichoice(9)=2    ! final spectrum accuracy high by defaut
	     ichoice(10)=2   ! 2-loop rad. corr. to Higgs masses by defaut
	     ichoice(11)=0  ! higher order Higgs scheme: DRbar masses in loops
      do ism=1,7
      smval(ism)=unlikely    !protection against undefined input
      enddo
c -- start from the beginning of the file suspect_slha.in --
      rewind(ninlha)

c -- initialization of the check array --
      do i1=1,24,1
         check(i1) = 0
      end do

c ------------------------------------------------------------------- c
      do i=1,10000,1

c -- check if routine can be left --
         check_final = 1
         do i1=1,24,1
            check_final = check_final*check(i1)
         end do
         if(check_final.eq.1) then
            return
         endif

c -- read in new line --
         line1=' '
         read(ninlha,'(a6,a100)',end=9900) line1,line2
         
c -- rewrite line1(1:6) and line2(1:20) to upper case --
         do j1=1,6,1
            if(line1(j1:j1).ne.'#') then
               do j2=97,122,1
                  if(line1(j1:j1).eq.char(j2)) line1(j1:j1)=char(j2-32)
               end do
            endif
         end do

         do j1=1,20,1
            if(line2(j1:j1).ne.'#') then
               do j2=97,122,1
                  if(line2(j1:j1).eq.char(j2)) line2(j1:j1)=char(j2-32)
               end do
            endif
         end do

c -- looks for blocks and reads them in one after the other --
         if(line1(1:1).eq.'B') then

c -- look for Block MODSEL --
            if(line2(1:6).eq.'MODSEL') then
               call SU_READ_MODSEL(ninlha,modselval,imod,done,endfile)
c translate SLHA model choice imode values into SUSpect model ichoice values:
                if(imod(2).eq.0) ichoice(1)=0    ! general MSSM at low scale
                if(imod(2).eq.1) ichoice(1)=10   ! mSUGRA
                if(imod(2).eq.2) ichoice(1)=11   ! GMSB
                if(imod(2).eq.3) ichoice(1)=12   ! AMSB
                if(imod(2).eq.10) ichoice(1)=1   ! general MSSM at high scale
                if(imod(2).eq.-1) ichoice(1)=2   ! EWSB input and bottom-up RGE
               if (done) then
                  check(21) = 1
                  if(endfile) then
                     goto 9900
                  else
                     goto 1111
                  endif
               else
                  print*,'SU_read_leshouches: problem in MODSEL'
               endif

c -- look for Block SU_ALGO --(SuSpect algorithm control parameters)
            elseif(line2(1:7).eq.'SU_ALGO') then
               call SU_READ_SU_ALGO(ninlha,ichoice,done,endfile)
               if (done) then
                  check(22) = 1
		  algo_warn=0d0
                  if (endfile) then
                     goto 9900
                  else
                     goto 1111
                  endif
               else
            algo_warn=-1d0
c case where specific SuSpect bloc SU_ALGO undefined: take defaut values
               endif

c -- look for Block SMINPUTS --
            elseif(line2(1:8).eq.'SMINPUTS') then
               call SU_READ_SMINPUTS(ninlha,smval,done,endfile)
      u=unlikely
      if(smval(1).ne.u) alfinv = smval(1)
      if(smval(1).eq.u.or.smval(1).eq.0d0) alfinv= 127.934d0  
c  this and following similar def. are protection againt undefined values
      if(smval(2).ne.u) gf = smval(2)
      if(smval(2).eq.u.or.smval(2).eq.0d0) gf = 1.16639d-5
      if(smval(3).ne.u) alphas = smval(3)
      if(smval(3).eq.u.or.smval(3).eq.0d0) alphas = .1172d0
c
      if(smval(4).ne.u) mz = smval(4)
      if(smval(4).eq.u.or.smval(4).eq.0d0) mz = 91.187d0
      if(smval(5).ne.u) mb     = smval(5)    ! mb is mb(mb)_MSbar input
      if(smval(5).eq.u.or.smval(5).eq.0d0) mb     = 4.25d0   
      if(smval(6).ne.u) mt     = smval(6)
      if(smval(6).eq.u.or.smval(6).eq.0d0) mt     = 175d0
      if(smval(7).ne.u) mtau   = smval(7)
      if(smval(7).eq.u.or.smval(7).eq.0d0) mtau   = 1.777d0
      do ism=1,7
      if(smval(ism).eq.u) smin_warn=-1d0
      enddo
c
               
               if (done) then
                  check(1) = 1
                  if(endfile) then
                     goto 9900
                  else
                     goto 1111
                  endif
               else
      smin_warn=-1d0
      alfinv = 127.934d0  
      gf = 1.16639d-5
      alphas = .1172d0
      mz = 91.187d0
      mb = 4.25d0
      mt = 175d0
      mtau = 1.777d0
c --then create Block SMINPUTS --
       smval(1)= alfinv
       smval(2)= gf 
       smval(3)= alphas 
       smval(4)= mz
       smval(5)= mb  
       smval(6)= mt     
       smval(7)= mtau 
               endif

c -- look for Block MINPAR --
            elseif(line2(1:6).eq.'MINPAR') then
               call SU_READ_MINPAR(ninlha,minval,mincom,done,endfile)
      if(ichoice(1).eq.10) then
c minimal SUGRA models with full universality: if non-universality, values
c are supersed by block EXTPAR below
      m0 = minval(1)
c
        mhd2= m0**2
	mhu2= m0**2
        MSL      = m0
        MTAUR    = m0
        MSQ      = m0
        MTR      = m0
        MBR      = m0
c
        MEL      = m0
        MER      = m0
        MUQ      = m0
        MUR      = m0
        MDR      = m0
c
      mhalf = minval(2)
c
      m1=mhalf         
      m2=m1         
      m3=m1         
c
      A0 = minval(5)
c
      At=A0         
      Ab=A0         
      Atau=A0         
c
      Au=A0         
      Ad=A0         
      Al=A0         
c
      tgbeta = minval(3)
      sgnmu0 = minval(4)
      elseif(ichoice(1).eq.11) then
c GMSB models
      mgmmess = minval(2)
      mgmsusy = minval(1)
      tgbeta = minval(3)
      sgnmu0 = minval(4)
      if(minval(5).eq.unlikely) minval(5)=1d0
      if(minval(6).eq.unlikely) minval(6)=1d0
      nl     = minval(5)
      nq     = minval(6)
c
      elseif(ichoice(1).eq.12) then
c AMSB models
      m32 = minval(2)
      am0 = minval(1)
      tgbeta = minval(3)
      sgnmu0 = minval(4)
      do ii=5,11
      if(minval(ii).eq.unlikely) minval(ii)=1d0
      enddo
      cq     = minval(5)
      cu     = minval(6)
      cd     = minval(7)
      cl     = minval(8)
      ce     = minval(9)
      chu    = minval(10)
      chd    = minval(11)
      endif   
c added 21/05/08 (jlk): read tbeta(mZ) from MINPAR even if modsel-0 or -1:
      if(imod(2).le.0.and.minval(3).ne.u) tgbeta=minval(3)
               if (done) then
                  check(2) = 1
                  if (endfile) then
                     goto 9900
                  else
                     goto 1111
                  endif
               else
                  print*,'SU_read_leshouches: problem in MINPAR'
               endif

c -- look for Block EXTPAR --
            elseif(line2(1:6).eq.'EXTPAR') then
               call SU_READ_EXTPAR(ninlha,extval,extcom,done,endfile)
      if(ichoice(8).eq.0.and.extval(0).ne.unlikely) Qewsb = extval(0)
c! essai      if(ichoice(1).eq.2) ehigh = extval(10)
        if(extval(10).ne.unlikely) ehigh = extval(10) 
c special case (new) to adapt non-universal SUGRA to former SuSpect:
      u=unlikely
c what follows is to avoid soft terms values with silly input. I didn't find
c a more astute way.      
      if(extval(1).ne.u) m1=extval(1)         
      if(extval(2).ne.u) m2=extval(2)         
      if(extval(3).ne.u) m3=extval(3)         
c
      if(extval(11).ne.u) At=extval(11)         
      if(extval(12).ne.u) Ab=extval(12)         
      if(extval(13).ne.u) Atau=extval(13)         
c
      if(extval(14).ne.u) Au=extval(14)         
      if(extval(15).ne.u) Ad=extval(15)         
      if(extval(16).ne.u) Al=extval(16)         
c
      if(extval(21).ne.u) mhd2 = extval(21)
      if(extval(22).ne.u) mhu2 = extval(22)
      if(extval(23).ne.u) then
      mu = extval(23)
        if(mu.ne.0d0) then
	sgnmu0 = mu/dsqrt(mu)
        else
	sgnmu0 = minval(4)! modif to define sgn(mu) if MU not input in EXTPAR
	endif
      endif
      if(extval(24).ne.u) madr2 = extval(24)   ! running DRbar m^2_A
      if(extval(25).ne.u) tgbeta = extval(25)
      if(extval(26).ne.u) ma = extval(26)     ! pole A mass
      if(extval(23).ne.u.and.extval(26).ne.u) then 
      ichoice(6)=0 
      else
      ichoice(6)=1
      endif
c      i.e impose MA,MU input algorithm in this case if user didn't specify
c      and impose Mhu,Mhd input otherwise
      if(extval(23).ne.u.and.extval(22).ne.u) mamu_warn=-1d0
c
      if(extval(33).ne.u) MSL = extval(33)
      if(extval(36).ne.u) MTAUR = extval(36)
      if(extval(43).ne.u) MSQ = extval(43)
      if(extval(46).ne.u) MTR = extval(46)
      if(extval(49).ne.u) MBR = extval(49)
c
      if(extval(31).ne.u) MEL = extval(31)
      if(extval(34).ne.u) MER = extval(34)
      if(extval(41).ne.u) MUQ = extval(41)
      if(extval(44).ne.u) MUR = extval(44)
      if(extval(47).ne.u) MDR = extval(47)
c

      if(endfile) then
         goto 9900
      endif

               if (done) then
                  check(23) = 1
                  if(endfile) then
                     goto 9900
                  else
                     goto 1111
                  endif
               else
                  print*,'SU_read_leshouches: problem in EXTPAR'
               endif

c -- continue if the Block is not interesting --
            else
               goto 1111
            endif

c -- continue if it is not a Block statement --
         else
            goto 1111
         endif

c -- maximum number of lines exhausted --
 1111    continue
      end do

c 9900 print*,'SU_read_leshouches: end of file'
 9900  continue
       return
c 
      else 
c --Case when no SLHA input file is read, but still prepare output in 
c  SLHA format if needed:      
c --create Block MODSEL --
       imod(1) = 1
       if(ichoice(1).eq.0)  imod(2) = 0
       if(ichoice(1).eq.10) imod(2) = 1
       if(ichoice(1).eq.11) imod(2) = 2
       if(ichoice(1).eq.12) imod(2) = 3
       if(ichoice(1).eq.1)  imod(2) = 10
       if(ichoice(1).eq.2)  imod(2) = 11
c
       if(imod(2).eq.0) modselval = 'general MSSM low scale'
       if(imod(2).eq.1) modselval = 'SUGRA'
       if(imod(2).eq.2) modselval = 'GMSB'
       if(imod(2).eq.3) modselval = 'AMSB'
       if(imod(2).eq.10) modselval = 'general MSSM High scale'
       if(imod(2).eq.11) modselval = 'general MSSM low scale'
c --create Block SMINPUTS --
       smval(1)= alfinv 
       smval(2)= gf
       smval(3)= alphas 
       smval(4)= mz
       smval(5)= mb  
       smval(6)= mt     
       smval(7)= mtau 
c
c create Block MINPAR --
      if(ichoice(1).eq.10) then
c mSUGRA models
      minval(1) = m0 
      minval(2) = mhalf
      minval(5) = A0 
      minval(3) = tgbeta  
      minval(4) = sgnmu0  
c
      mincom(1) = 'm0' 
      mincom(2) = 'm1%2'
      mincom(5) = 'A0' 
      mincom(3) = 'tanbeta'  
      mincom(4) = 'sign(mu)'  

      elseif(ichoice(1).eq.11) then
c GMSB models
      minval(2) = mgmmess  
      minval(1) = mgmsusy  
      minval(3) = tgbeta  
      minval(4) = sgnmu0  
      minval(5) = nl     
      minval(6) = nq   
c
      mincom(2) = 'Lambda_mess'  
      mincom(1) = 'Lambda_susy'  
      mincom(3) = 'tanbeta'  
      mincom(4) = 'sign(mu)'  
      mincom(5) = 'Nl_mes'     
      mincom(6) = 'Nq_mes'  

      elseif(ichoice(1).eq.12) then
c AMSB models
      minval(2) = m32 
      minval(1) = am0  
      minval(3) = tgbeta  
      minval(4) = sgnmu0 
      minval(5) = cq     
      minval(6) = cu   
      minval(7) = cd
      minval(8) = cl
      minval(9) = ce
      minval(10) = chu
      minval(11) = chd
c
      mincom(2) = 'M_3%2' 
      mincom(1) = 'm0'  
      mincom(3) = 'tanbeta'  
      mincom(4) = 'sgn(mu)' 
      mincom(5) = 'cQ'     
      mincom(6) = 'cuR'   
      mincom(7) = 'cdR'
      mincom(8) = 'cL'
      mincom(9) = 'ceR'
      mincom(10) = 'cHu'
      mincom(11) = 'cHd'
      endif   
c -- create Block EXTPAR --
c a trick to jump over undefined parameters in subsequent writings:
      unlikely = -123456789D0
      do i=0,60,1
         extval(i) = unlikely
      end do
      if(ichoice(8).eq.0) then
      extval(0) = Qewsb
      endif
      extcom(0) = 'EWSB_scale'
c!essai      if(ichoice(1).eq.2) then
      if(ehigh.ne.0d0) then
      extval(10) = ehigh
      if(ichoice(1).eq.2) then
      extcom(10) = 'RGE final scale'
      else
      extcom(10) = 'RGE initial scale'
      endif
      endif
c
      extval(1) = m1        
      extval(2) = m2        
      extval(3) = m3        
c
      extval(11) = At        
      extval(12) = Ab        
      extval(13) = Atau        
c
      extval(14) = Au        
      extval(15) = Ad        
      extval(16) = Al        
c
      extval(21) = mhd2
      extval(22) = mhu2
      extval(23) = mu
      extval(24) = madr2 
      extval(25) = tgbeta
      extval(26) = ma
c
      extval(33) = MSL
      extval(36) = MTAUR
      extval(43) = MSQ
      extval(46) = MTR
      extval(49) = MBR
c
      extval(31) = MEL
      extval(34) = MER
      extval(41) = MUQ
      extval(44) = MUR
      extval(47) = MDR
c
      extcom(1) = 'M1'        
      extcom(2) = 'M2'        
      extcom(3) = 'M3'        
c
      extcom(11) = 'A_t'        
      extcom(12) = 'A_b'        
      extcom(13) = 'A_tau'        
c
      extcom(14) = 'A_u'        
      extcom(15) = 'A_d'        
      extcom(16) = 'A_e'        
c
      extcom(21) = 'M^2_Hd'
      extcom(22) = 'M^2_Hu'
      extcom(23) = 'MU(EWSB scale)'
      extcom(24) = 'M^2_A(run,EWSB)' 
      extcom(25) = 'tanbeta'
      extcom(26) = 'M_A(pole)'
c
      extcom(33) = 'M_tau_L'
      extcom(36) = 'M_tau_R'
      extcom(43) = 'M_Q_L'
      extcom(46) = 'M_t_R'
      extcom(49) = 'M_b_R'
c
      extcom(31) = 'M_e_L'
      extcom(34) = 'M_e_R'
      extcom(41) = 'M-qu_L'
      extcom(44) = 'M_u_R'
      extcom(47) = 'M_u_R'
c
      endif
      end

c -------------------------------------------------------------------- c

      subroutine SU_READ_MODSEL(ninlha,modselval,imod,done,endfile)

      implicit double precision (a-h,m,o-z)
      integer imod(1:2)
      character line1*1,line2*1,line3*100,modselval*100
      logical done,endfile

      done=.false.
      endfile=.false.

      modselval = ' '

      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
            read(ninlha,*) idum1,idum2   ! removed: ,line2,line3

            if(idum1.eq.1) then
               imod(1) = idum1
               imod(2) = idum2
c               modselval = line3
             if(imod(2).eq.0) modselval = 'general MSSM low scale'
             if(imod(2).eq.1) modselval = 'SUGRA'
             if(imod(2).eq.2) modselval = 'GMSB'
             if(imod(2).eq.3) modselval = 'AMSB'
             if(imod(2).eq.-1) modselval = 'bottom-up MSSM'
            endif

         elseif(line1.eq.'#') then
            go to 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done =.true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end

c -------------------------------------------------------------------- c
       subroutine SU_READ_SU_ALGO(ninlha,ichoice,done,endfile)

      implicit double precision (a-h,m,o-z)
      integer ichoice(1:11)
      character line1*1
      logical done,endfile

      done=.false.
      endfile=.false.

      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
            read(ninlha,*) idum,val

c -- The different suspect options ichoice(2)-ichoice(11)
            if(idum.eq.2) then
               ichoice(2) = val
c -- 
            elseif(idum.eq.3) then
               ichoice(3) = val
            elseif(idum.eq.4) then
               ichoice(4) = val
            elseif(idum.eq.6) then
               ichoice(6) = val
            elseif(idum.eq.7) then
               ichoice(7) = val
            elseif(idum.eq.8) then
               ichoice(8) = val
            elseif(idum.eq.9) then
               ichoice(9) = val
            elseif(idum.eq.10) then
               ichoice(10) = val
            elseif(idum.eq.11) then
               ichoice(11) = val
            endif
         elseif(line1.eq.'#') then
            go to 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done =.true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end
c---------------------------------------------------------------------
      subroutine SU_READ_SMINPUTS(ninlha,smval,done,endfile)

      implicit double precision (a-h,m,o-z)
      double precision smval(20)
      character line1*1
      logical done,endfile

      done=.false.
      endfile=.false.

      do i=1,20,1
         smval(i) = 0.D0
      end do

      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
            read(ninlha,*) idum,val

c -- inverse EM coupling at the Z pole in the MS_bar scheme (with --
c -- five active flavours) --
            if(idum.eq.1) then
               smval(1) = val
c -- G_F, Fermi constant (in units of GeV^-2)
            elseif(idum.eq.2) then
               smval(2) = val
c -- Strong coupling at the Z pole in the MS_bar scheme (with five --
c -- active flavours) --
            elseif(idum.eq.3) then
               smval(3) = val
c -- M_Z, pole mass --
            elseif(idum.eq.4) then
               smval(4) = val
c -- mb(mb)^MS_bar. b quark running mass in the MS_bar scheme --
            elseif(idum.eq.5) then
               smval(5) = val
c -- mt, pole mass --
            elseif(idum.eq.6) then
               smval(6) = val
c -- mtau, pole mass --
            elseif(idum.eq.7) then
               smval(7) = val
            endif
            
         elseif(line1.eq.'#') then
            go to 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done =.true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end

c -------------------------------------------------------------------- c

      subroutine SU_READ_MINPAR(ninlha,minval,mincom,done,endfile)

      implicit double precision (a-h,m,o-z)
      double precision minval(20)
      character line1*1,line2*1,line3*20,mincom(1:20)*20
      logical done,endfile

      done= .false.
      endfile= .false.
      unlikely=-123456789d0
      do i=1,20,1
         minval(i) = unlikely
      end do

      do i=1,20,1
         mincom(i) = ' '
      end do

      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
c!jlk            read(ninlha,*) idum,val,line2,line3
            read(ninlha,*) idum,val
               do ii=1,11,1
               if(idum.eq.ii) then
                  minval(ii) = val
c                  mincom(ii) = line3
               endif
            end do

c -- i=3: value for tanbeta(MZ) --
            
         elseif(line1.eq.'#') then
            goto 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done = .true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end

c -------------------------------------------------------------------- c
       subroutine SU_READ_EXTPAR(ninlha,extval,extcom,done,endfile)

      implicit double precision (a-h,m,o-z)
      dimension extval(0:60)
      character line1*1,line2*1,line3*20,extcom(0:60)*20
      logical done,endfile

      done=.false.
      endfile=.false.
c a trick to jump over undefined parameters:
      unlikely = -123456789D0
      do i=0,60,1
         extval(i) = unlikely
      end do

      do i=0,60,1
         extcom(i) = ' '
      end do
      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
c!jlk            read(ninlha,*) idum,val,line2,line3
           read(ninlha,*) idum,val
c -- The general MSSM model parameters according to SLHA nomenclature:
            do ii=0,60,1
            if(idum.eq.ii) then
            extval(ii) = val
c!jlk            extcom(ii) = line3
            endif
            enddo
c -- 
            
         elseif(line1.eq.'#') then
            go to 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done =.true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end
c -------------------------------------------------------------------- c



      subroutine SU_READ_SPINFO(ninlha,spinfo1,spinfo2,done,endfile)

      implicit double precision (a-h,m,o-z)
      character line1*1,line2*100,spinfo1*100,spinfo2*100
      logical done,endfile

      done= .false.
      endfile=.false.

      spinfo1 = ' '
      spinfo2 = ' '

      do i=1,200,1
         read(ninlha,'(a1)',end=9900) line1

c -- decide what it is and read the line if anything of interest --
         if (line1.eq.' ') then
            backspace ninlha
            read(ninlha,'(1x,i5,3x,a100)') idum,line2

c -- the name of the spectrum calculator --
            if(idum.eq.1) then
               spinfo1 = line2

c -- the version number of the spectrum calculator --
            elseif(idum.eq.2) then
               spinfo2 = line2
            endif

         elseif(line1.eq.'#') then
            goto 1111
         elseif(line1.eq.'b'.or.line1.eq.'B'.or.line1.eq.'d'.or.line1.eq
     ..'D') then
            backspace ninlha
            done = .true.
            return
         endif

 1111    continue
      end do

 9900 print*,'SU_read_leshouches: end of file'
      done = .true.
      endfile = .true.

      end
c--------------------------------------------------------------------------
c%% routine for writing SuSpect ver >= 2.3 output in SLHA form
c  released J-L Kneur 06/12/2004 
c--thanks to Margarete Muhlleitner for adapting simply from her writing --
c-----------------------------------------------------
      subroutine su_lhaout(nout,ichoice,errmess,imod)
      implicit real*8 (a-h,m,o-z)
      real*8 nl,nq
      double precision minval(1:20),extval(0:60),smval(1:20),
     .       massval(1:50),
     .       nmixval(4,4),umixval(2,2),vmixval(2,2),stopmixval(2,2),
     .       sbotmixval(2,2),staumixval(2,2),hmixval(1:10),
     .       gaugeval(1:3),msoftval(1:30),auval(3,3),adval(3,3),
     .       aeval(3,3),yuval(3,3),ydval(3,3),yeval(3,3)
      integer nx1t,ny1t,nnlo,imod(1:2)
      character spinfo1*100,spinfo2*100,modselval*100,mincom(1:20)*20,
     .          extcom(0:60)*20

      dimension ichoice(11),errmess(10)
      dimension amneut(4),xmneut(4),amchar(2)
      dimension uu(2,2),vv(2,2),zz(4,4),zp(4,4)
      COMMON/SU_strc/irge,irgmax,ifix,isfrc,inorc
      COMMON/SU_SMPAR/dalfinv,dsw2,dalphas,dmt,dmb,dmc,dmtau
      COMMON/SU_RGSCAL/qewsb,ehigh,elow
      COMMON/SU_MSSMHPAR/mhu2,mhd2,dma,dmu
      COMMON/SU_MSSMGPAR/dm1,dm2,dm3 
      COMMON/SU_MSSMSLEP/dmsl,dmtaur,dmel,dmer
      COMMON/SU_MSSMSQUA/dmsq,dmtr,dmbr,dmuq,dmur,dmdr
      COMMON/SU_ATRI3/dal,dau,dad
      COMMON/SU_ATRI12/dal1,dau1,dad1
      COMMON/SU_MSUGRA/m0,mhalf,a0
      COMMON/SU_RADEWSB/sgnmu0,tgbeta
      COMMON/SU_GMSB/mgmmess,mgmsusy,nl,nq
      COMMON/SU_AMSB/m32,am0,cq,cu,cd,cl,ce,chu,chd
      COMMON/SU_matino/uu,vv,zz,xmneut
      COMMON/SU_outhiggs/aml,amh,amch,alfa
c  light, heavy, charged Higgs masses, neutral (h,H) mix angle alpha 
      COMMON/SU_outginos/dmc1,dmc2,dmn1,dmn2,dmn3,dmn4,mgluino
c   charginos 1,2 masses, neutralinos 1-4 masses, gluino mass 
      COMMON/SU_outsqu/dmst1,dmst2,dmsu1,dmsu2
c  stop 1,2 and sup 1,2 = scharm 1,2 masses
      COMMON/SU_outsqd/dmsb1,dmsb2,dmsd1,dmsd2
c  sbottom 1,2 and sdown 1,2 = sstrange 1,2 masses
      COMMON/SU_outslep/dmsl1,dmsl2,dmse1,dmse2,dmsn1,dmsntau
c  stau 1,2 ; selectron (=smuon) 1,2; sneut_e,mu, sneut_tau masses
      COMMON/SU_outmix/thet,theb,thel
c  stop, sbottom, stau mixing angles
      COMMON/SU_param/gf,alpha,mz,mw
      COMMON/SU_fmasses/mtau,mbpole,mtpole
      COMMON/SU_yukaewsb/ytauewsb,ybewsb,ytewsb,alsewsb,g2ewsb,g1ewsb
      COMMON/SU_tbewsb/vuewsb,vdewsb 
      COMMON/SU_renscale/scale
      common/SU_ftune/czmu,czbmu,ctmu,ctbmu
c low-energy contrained parameter values: rho-1, g_mu-2, Br(b->s gamma):
      COMMON/SU_lowen/crho,gmuon,brsg
      common/su_runmavev/madr2,vev2
c -------------- common block given by SD_read_leshouches ------------ c
      COMMON/SU_leshouches1/spinfo1,spinfo2,modselval,mincom,extcom
      COMMON/SU_leshouches2/minval,extval,smval,massval,nmixval,umixval,
     .                      vmixval,stopmixval,sbotmixval,staumixval,
     .                      hmixval,gaugeval,msoftval,auval,adval,
     .                      aeval,yuval,ydval,yeval,alphaval,Qvalhmix,
     .                      Qvalgauge,Qvalmsoft,Qvalau,Qvalad,Qvalae,
     .                      Qvalyu,Qvalyd,Qvalye

      common/SU_slha_warn/smin_warn,extpar_warn,muma_warn,algo_warn
        pi=4*datan(1.d0)
c completing input/output in slha variables:
      smval(2)=gf
      smval(4)=mz
c PDG values:
      id =1
      idb=-1
      iu =2
      iub=-2
      is =3
      isb=-3
      ic =4
      icb=-4
      ib =5
      ibb=-5
      it =6
      itb=-6

      ie   =11
      ine  =12
      imu  =13
      inmu =14
      itau =15
      intau=16

      ihl=25
      ihh=35
      iha=36
      ihc=37
      igl=21
      iga=22
      iz =23
      iwc=24

      isdl=1000001
      isdr=2000001
      isul=1000002
      isur=2000002
      issl=1000003
      issr=2000003
      iscl=1000004
      iscr=2000004
      isb1=1000005
      isb2=2000005
      ist1=1000006
      ist2=2000006

      iglo=1000021
      in1 =1000022
      in2 =1000023
      in3 =1000025
      in4 =1000035
      ic1 =1000024
      ic2 =1000037

      intau1=1000016 
      intau2=2000016 
      inel  =1000012
      iner  =2000012
      inmul =1000014
      inmur =2000014
      
      isell =1000011
      iselr =2000011
      ismul =1000013
      ismur =2000013
      istau1=1000015
      istau2=2000015

      igrav =1000039
c hardcoding input names of block MINPAR and EXTPAR for some platforms  
c compatibility (and to eventually avoid input file typing confusion):
c added 02/06/2008 jlk
      mincom(3) = ' tanbeta(mZ)'
      mincom(4) = ' sign(mu)'
      if(imod(2).eq.1) then
c input for sugra models
      mincom(1) = ' m0'
      mincom(2) = ' m_1/2'
      mincom(5) = ' A0'
      else if(imod(2).eq.2) then
c    input for GMSB models:
      mincom(1) = ' Lambda_susy'
      mincom(2) = ' Lambda_mess'
      mincom(5) = ' Nl_mes'  
      mincom(6) = ' Nq_mes'  
      else if(imod(2).eq.3) then
c    input for AMSB models:
      mincom(1) = ' m0'
      mincom(2) = ' M_3/2 gravino'
      mincom(5) = ' cQ coeff m0 Q_L' 
      mincom(6) = ' cuR coeff m0 u_R'
      mincom(7) = ' cdR coeff m0 d_R' 
      mincom(8) = ' cL coeff m0 L' 
      mincom(9) = ' ceR coeff m0 e_R'
      mincom(10) = ' cHu coeff m0 Hu' 
      mincom(11) = ' cHd coeff m0 Hd' 
      endif
c def. of EXTPAR names:
        extcom(0) = ' EWSB scale'          
        extcom(10) = ' GUT scale'
	extcom(23) = ' mu(EWSB)'
	extcom(24) = ' m^2_A_run(EWSB)'
        extcom(25) = ' tanbeta(in)'
	extcom(26) = ' MA_pole'
        extcom(1) = ' M_1'
        extcom(2) = ' M_2'
        extcom(3) = ' M_3'
        extcom(21) = ' M^2_Hd'
        extcom(22) = ' M^2_Hu'
        extcom(31) = ' M_eL'
        extcom(32) = ' M_muL'
        extcom(33) = ' M_tauL'
        extcom(34) = ' M_eR'
        extcom(35) = ' M_muR'
        extcom(36) = ' M_tauR'
        extcom(41) = ' M_q1L'
        extcom(42) = ' M_q2L'
        extcom(43) = ' M_q3L'
        extcom(44) = ' M_uR'
        extcom(45) = ' M_cR'
        extcom(46) = ' M_tR'
        extcom(47) = ' M_dR'
        extcom(48) = ' M_sR'
        extcom(49) = ' M_bR'
        extcom(11) = ' A_t'
        extcom(12) = ' A_b'
        extcom(13) = ' A_tau'
        extcom(14) = ' A_u'
        extcom(15) = ' A_d'
        extcom(16) = ' A_e'

      write(nout,105)
      write(nout,50) "                              ====================
     .="
      write(nout,50) "                             | SuSpect 2.41 OUTPUT 
     . |"
      write(nout,50) "                              ====================
     .="
      write(nout,105)
      write(nout,105)

      write(nout,50)'             --------------------------------------
     .---------------'
      write(nout,50)'             |  SUSY Les Houches Accord - MSSM Spec
     .trum          |'
      write(nout,50)'             |                                     
     .              |'
      write(nout,50)'             |                     SuSpect 2.41      
     .              |'
      write(nout,50)'             |                                     
     .              |'
      write(nout,50)'             |  Authors: A.Djouadi, J.-L. Kneur and 
     . G. Moultaka  |'
      write(nout,50)'             |  Ref.:    hep-ph/0211331            
     .              |'
      write(nout,50)'             |                                     
     .              |'
      write(nout,50)'             --------------------------------------
     .---------------'
      write(nout,105)

c -------------------------------------------- c
c Information about the RGE + spectrum program c
c -------------------------------------------- c

      write(nout,105)
      write(nout,51) 'SPINFO','Spectrum Program information'
      write(nout,61) 1,'SuSpect     # RGE +Spectrum calculator'
      write(nout,61) 2,'2.41        # version number'
c
c The SuSpect warning/error flag section
c
       warnerr=0.d0
      do ii=1,10
      if(errmess(ii).eq.-1.d0) warnerr=1d0
      enddo
      if(irgmax.eq.50) warnerr=1d0
      if(warnerr.eq.0d0) then
      write(nout,'(a)')'# nothing to signal: output a priori reliable'
      else
      write(nout,'(a)')'# Caution: warning or error message follows '
      endif
      if(errmess(1).eq.-1.d0) then
      write(nout,61) 4,'Bad input: one m^2(3rd gen. sf) <0 from RGE '
      endif
      if(errmess(2).eq.-1.d0) then
      write(nout,61) 4,'Bad input: one m^2(1,2 gen. sf) <0 from RGE '
      endif
      if(errmess(3).eq.-1.d0) then
      write(nout,61) 3,'Warning:  MA^2(Q) <0 at a scale MZ<Q<EWSB ! '
      endif
      if(errmess(4).eq.-1.d0) then
      write(nout,61) 4,'STOP: one tachyonic m^2(3rd gen. sf) <0 '
      endif
      if(errmess(5).eq.-1.d0) then
      write(nout,61) 3,' Warning: MU unstable after many iter'
      endif
      if(errmess(6).eq.-1.d0) then
      write(nout,61) 3,'WARNING: EWSB unconvergent after 20 iter.' 
      endif
      if(errmess(7).eq.-1.d0) then
      write(nout,61) 4,'EWSB  maybe unconsistent/not realized '
      endif
      if(errmess(8).eq.-1.d0) then
      write(nout,61) 3, 'RG-improved V_eff has CCB or UFB problems '
      endif
      if(errmess(9).eq.-1.d0) then
      write(nout,61) 4, ' PROBLEM: some Higgs masses are NaN! '
      endif
      if(errmess(10).eq.-1.d0) then
      write(nout,61) 4,'STOP: non-pert. R.C., or Landau pole in RGE'
      endif
      if(irgmax.eq.50) then
      write(nout,61) 3,'Pb: non-convergent spectrum after 50 iter.!'
      endif
      if(smin_warn.eq.-1d0) then
      write(nout,61) 3,'warning: some SM input undefined: defaut taken'
      endif
      if(muma_warn.eq.-1d0) then
      write(nout,61) 3,'warning: redondent input MU,MA,Mhu,Mhd'
      endif
      if(algo_warn.eq.-1d0) then
      write(nout,61) 3,'warning: miss block SU_ALGO: defaut val. taken'
      write(nout,61) 3,'see original SuSpect input file for values!'
      endif

c ------------------------------------------------ c
c Information on the model which has been selected c
c ------------------------------------------------ c

      write(nout,105)
      write(nout,51) 'MODSEL','Model selection'
         write(nout,62) imod(1),imod(2),modselval(1:50)

c ------------------------------------------------------------------- c
c The input parameters for the different model choice given above     c
c ------------------------------------------------------------------- c
      unlikely=-123456789d0
      write(nout,105)
      write(nout,51) 'MINPAR','Input parameters'
         if(ichoice(1).eq.1) goto 1
         if(ichoice(1).eq.10) iimax=5 
         if(ichoice(1).eq.11) iimax=6 
         if(ichoice(1).eq.12) iimax=11 
         if(ichoice(1).eq.0.or.ichoice(1).eq.2) iimax=4 !modif to read tb(mz)
         do i=1,iimax,1                                 !still from MINPAR     
      if(minval(i).ne.unlikely) write(nout,52) i,minval(i),mincom(i)
         end do
	 if(ichoice(1).ge.10) then
      write(nout,105)
      write(nout,51) 'EXTPAR','Input parameters'
      if(ichoice(8).eq.1) extval(0) = scale
      extcom(0) ='EWSB scale'
      write(nout,72) 0,extval(0),extcom(0)
        if(extval(10).ne.unlikely.and.extval(10).ne.0d0) then
	extcom(10)='High boundary scale'
      write(nout,72) 10,extval(10),extcom(10)
        endif
         endif
c!         else
c!!          if(ichoice(1).eq.0.or.ichoice(1).eq.2) then
 1         if(ichoice(1).le.2) then
      write(nout,51) 'EXTPAR','Input parameters'
         unlikely = -123456789d0
c!        if(ichoice(1).ne.2) extval(10) = unlikely
      if(ichoice(8).eq.1) extval(0) = scale
      extcom(0) ='EWSB scale'
      if(ichoice(6).eq.0) then
c (case of MA_pole, MU input in general MSSM)
         extval(22) = unlikely
         extval(21) = unlikely
      else if(ichoice(6).eq.1) then
c (case of M^2_Hu, M^2_Hd input in general MSSM)
        extval(26) = unlikely
        extval(23) = unlikely
      endif
         do i=0,60,1
      if(extval(i).ne.unlikely) write(nout,72) i,extval(i),extcom(i)
         end do
         endif
c ----------------------- c
c The SM input parameters c
c ----------------------- c

      write(nout,105)
      write(nout,51) 'SMINPUTS','Standard Model inputs'
      write(nout,52) 1,dalfinv,'alpha_em^-1(M_Z)^MSbar'
      write(nout,52) 2,gf,'G_F [GeV^-2]'
      write(nout,52) 3,dalphas,'alpha_S(M_Z)^MSbar'
      write(nout,52) 4,mz,'M_Z pole mass'
      write(nout,52) 5,dmb,'mb(mb)^MSbar'
      write(nout,52) 6,dmt,'mt pole mass'
      write(nout,52) 7,dmtau,'mtau pole mass'
c ----------------- c
c The mass spectrum c
c ----------------- c

      write(nout,105)
      write(nout,51) 'MASS','Mass Spectrum'
      write(nout,50) 'PDG code           mass       particle'
      write(nout,52) iwc,mw,'W+'
      write(nout,52) ihl,aml,'h'
      write(nout,52) ihh,amh,'H'
      write(nout,52) iha,dma,'A'
      write(nout,52) ihc,amch,'H+'
      write(nout,52) ib,mbpole,'b pole mass calculated from mb(mb)_MSbar
     .'
      write(nout,52) isdl,dmsd1,'~d_L'
      write(nout,52) isdr,dmsd2,'~d_R'
      write(nout,52) isul,dmsu1,'~u_L'
      write(nout,52) isur,dmsu2,'~u_R'
      write(nout,52) issl,dmsd1,'~s_L'
      write(nout,52) issr,dmsd2,'~s_R'
      write(nout,52) iscl,dmsu1,'~c_L'
      write(nout,52) iscr,dmsu2,'~c_R'
      write(nout,52) isb1,dmsb1,'~b_1'
      write(nout,52) isb2,dmsb2,'~b_2'
      write(nout,52) ist1,dmst1,'~t_1'
      write(nout,52) ist2,dmst2,'~t_2'
      write(nout,52) isell,dmse1,'~e_L'
      write(nout,52) iselr,dmse2,'~e_R'
      write(nout,52) inel,dmsn1,'~nu_eL'
      write(nout,52) ismul,dmse1,'~mu_L'
      write(nout,52) ismur,dmse2,'~mu_R'
      write(nout,52) inmul,dmsn1,'~nu_muL'
      write(nout,52) istau1,dmsl1,'~tau_1'
      write(nout,52) istau2,dmsl2,'~tau_2'
      write(nout,52) intau1,dmsntau,'~nu_tauL'
      write(nout,52) iglo,mgluino,'~g'
      write(nout,52) in1,xmneut(1),'~chi_10'
      write(nout,52) in2,xmneut(2),'~chi_20'
      write(nout,52) in3,xmneut(3),'~chi_30'
      write(nout,52) in4,xmneut(4),'~chi_40'
      write(nout,52) ic1,dmc1,'~chi_1+'
      write(nout,52) ic2,dmc2,'~chi_2+'

c The constrained low-energy or LEP2 parameter value for info:
c
      write(nout,105)
      write(nout,51) 'SU_LOWPAR','Values constrained by exp data '
      write(nout,52) 1,crho,'Delta rho parameter'
      write(nout,52) 2,gmuon,'g_mu -2'
      write(nout,52) 3,brsg,'Br(b -> s gamma)'

c The main fine-tuning parameter values for info:
c
      write(nout,105)
      write(nout,51) 'SU_FINETUNE','Fine-tuning info: fine-tuned if >>1'
      write(nout,52) 1,czmu,'delta mZ^2/mZ^2 (mu^2)'
      write(nout,52) 2,czbmu,'delta mZ^2/mZ^2 (B.mu)'
      write(nout,52) 3,ctmu,'delta mt/mt (mu^2)'
      write(nout,52) 4,ctbmu,'delta mt/mt  (B.mu)' 

c ------------------------------------------------------------------- c
c The neutralino mixing matrix N and the chargino mixing matrices U,V c
c ------------------------------------------------------------------- c
      
      write(nout,105)
      write(nout,51) 'NMIX','Neutralino Mixing Matrix'
      write(nout,53) 1,1,zz(1,1),'N_11'
      write(nout,53) 1,2,zz(1,2),'N_12'
      write(nout,53) 1,3,zz(1,3),'N_13'
      write(nout,53) 1,4,zz(1,4),'N_14'
      write(nout,53) 2,1,zz(2,1),'N_21'
      write(nout,53) 2,2,zz(2,2),'N_22'
      write(nout,53) 2,3,zz(2,3),'N_23'
      write(nout,53) 2,4,zz(2,4),'N_24'
      write(nout,53) 3,1,zz(3,1),'N_31'
      write(nout,53) 3,2,zz(3,2),'N_32'
      write(nout,53) 3,3,zz(3,3),'N_33'
      write(nout,53) 3,4,zz(3,4),'N_34'
      write(nout,53) 4,1,zz(4,1),'N_41'
      write(nout,53) 4,2,zz(4,2),'N_42'
      write(nout,53) 4,3,zz(4,3),'N_43'
      write(nout,53) 4,4,zz(4,4),'N_44'

      write(nout,105)
      write(nout,51) 'UMIX','Chargino Mixing Matrix U'
      write(nout,53) 1,1,uu(1,1),'U_11'
      write(nout,53) 1,2,uu(1,2),'U_12'
      write(nout,53) 2,1,uu(2,1),'U_21'
      write(nout,53) 2,2,uu(2,2),'U_22'

      write(nout,105)
      write(nout,51) 'VMIX','Chargino Mixing Matrix V'
      write(nout,53) 1,1,vv(1,1),'V_11'
      write(nout,53) 1,2,vv(1,2),'V_12'
      write(nout,53) 2,1,vv(2,1),'V_21'
      write(nout,53) 2,2,vv(2,2),'V_22'

c ------------------------------------------ c
c The stop, sbottom and stau mixing matrices c
c ------------------------------------------ c

      write(nout,105)
      write(nout,51) 'STOPMIX','Stop Mixing Matrix'
      write(nout,53) 1,1,dcos(thet),'cos(theta_t)'
      write(nout,53) 1,2,dsin(thet),'sin(theta_t)'
      write(nout,53) 2,1,-dsin(thet),'-sin(theta_t)'
      write(nout,53) 2,2,dcos(thet),'cos(theta_t)'

      write(nout,105)
      write(nout,51) 'SBOTMIX','Sbottom Mixing Matrix'
      write(nout,53) 1,1,dcos(theb),'cos(theta_b)'
      write(nout,53) 1,2,dsin(theb),'sin(theta_b)'
      write(nout,53) 2,1,-dsin(theb),'-sin(theta_b)'
      write(nout,53) 2,2,dcos(theb),'cos(theta_b)'

      write(nout,105)
      write(nout,51) 'STAUMIX','Stau Mixing Matrix'
      write(nout,53) 1,1,dcos(thel),'cos(theta_tau)'
      write(nout,53) 1,2,dsin(thel),'sin(theta_tau)'
      write(nout,53) 2,1,-dsin(thel),'-sin(theta_tau)'
      write(nout,53) 2,2,dcos(thel),'cos(theta_tau)'

c ------------------------------------------------------------------- c
c The angle alpha in the Higgs sector and the Higgs mixing parameters c
c ------------------------------------------------------------------- c

      write(nout,105)
      write(nout,51) 'ALPHA','Higgs mixing'
      write(nout,60) alfa,'Mixing angle in the neutral Higgs boson secto
     .r'

      write(nout,105)
      write(nout,54) 'HMIX Q=',scale,'DRbar Higgs Parameters'
      write(nout,55) 1,dmu,'mu(Q)'
      write(nout,55) 2,vuewsb/vdewsb,'tanbeta(Q)'
      write(nout,55) 3,dsqrt(vev2),'vev(Q)'
      write(nout,55) 4,madr2,'MA^2(Q)'

c ------------------- c
c The gauge couplings c
c ------------------- c

      write(nout,105)      
      write(nout,54) 'GAUGE Q=',scale,'The gauge couplings'
      write(nout,55) 1,g1ewsb,'gprime(Q) DRbar'
      write(nout,55) 2,g2ewsb,'g(Q) DRbar'
      write(nout,55) 3,dsqrt(4*pi*alsewsb),'g_3(Q) DRbar'
c ------------------------------------- c
c The trilinear couplings Au, Ad and Ae c
c ------------------------------------- c
      scalesave=scale
      if(ichoice(1).eq.2) scale=ehigh
      write(nout,105)
      write(nout,54) 'Au Q=',scale,'The trilinear couplings'
      write(nout,53) 1,1,dau1, 'A_u(Q) DRbar'
      write(nout,53) 2,2,dau1, 'A_c(Q) DRbar'
      write(nout,53) 3,3,dau,'A_t(Q) DRbar'

      write(nout,105)
      write(nout,54) 'Ad Q=',scale,'The trilinear couplings'
      write(nout,53) 1,1,dad1,'A_d(Q) DRbar'
      write(nout,53) 2,2,dad1,'A_s(Q) DRbar'
      write(nout,53) 3,3,dad ,'A_b(Q) DRbar'

      write(nout,105)
      write(nout,54) 'Ae Q=',scale,'The trilinear couplings'
      write(nout,53) 1,1,dal1 ,'A_e(Q) DRbar'
      write(nout,53) 2,2,dal1 ,'A_mu(Q) DRbar'
      write(nout,53) 3,3,dal,'A_tau(Q) DRbar'

c ---------------------------------- c
c The Yukawa couplings Yu, Yd and Ye c
c ---------------------------------- c

      write(nout,105)
      write(nout,54) 'Yu Q=',scalesave,'The Yukawa couplings'
         write(nout,53) 3,3,ytewsb,'y_top(Q) DRbar'
c
      write(nout,105)
      write(nout,54) 'Yd Q=',scalesave,'The Yukawa couplings'
         write(nout,53) 3,3,ybewsb,'y_b(Q) DRbar'
c
      write(nout,105)
      write(nout,54) 'Ye Q=',scalesave,'The Yukawa couplings'
         write(nout,53) 3,3,ytauewsb,'y_tau(Q) DRbar'

c ----------------------------- c
c The soft SUSY breaking masses c
c ----------------------------- c

      write(nout,105)
      write(nout,54) 'MSOFT Q=',scale,'soft SUSY breaking masses at the 
     .scale Q'
      write(nout,52) 1,dm1,'M_1'
      write(nout,52) 2,dm2,'M_2'
      write(nout,52) 3,dm3,'M_3'
         write(nout,52) 21,mhd2,'M^2_Hd'
         write(nout,52) 22,mhu2,'M^2_Hu'
      write(nout,52) 31,dmel,'M_eL'
      write(nout,52) 32,dmel,'M_muL'
      write(nout,52) 33,dmsl,'M_tauL'
      write(nout,52) 34,dmer,'M_eR'
      write(nout,52) 35,dmer,'M_muR'
      write(nout,52) 36,dmtaur,'M_tauR'
      write(nout,52) 41,dmuq,'M_q1L'
      write(nout,52) 42,dmuq,'M_q2L'
      write(nout,52) 43,dmsq,'M_q3L'
      write(nout,52) 44,dmur,'M_uR'
      write(nout,52) 45,dmur,'M_cR'
      write(nout,52) 46,dmtr,'M_tR'
      write(nout,52) 47,dmdr,'M_dR'
      write(nout,52) 48,dmdr,'M_sR'
      write(nout,52) 49,dmbr,'M_bR'
c
      write(*,'(a)')' OUTPUT in SLHA format in slhaspectrum.in '

 50   format('#',1x,A)
 51   format('BLOCK',1x,A,2x,'#',1x,A)
 52   format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 53   format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
 54   format('BLOCK',1x,A,1P,E16.8,2x,'#',1x,A)
 55   format(1x,I5,3x,1P,E16.8,0P,3x,'#',1x,A)
 60   format(9x,1P,E16.8,0P,3x,'#',1x,A)
 61   format(1x,I5,3x,A)
 62   format(1x,I5,1x,I5,3x,'#',A)
 72   format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 105  format('#') 
 
c -------------------------------------------------------------------------
      end
c%%%%%%%%%%%%%%%%%%%%   END OF THE PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
