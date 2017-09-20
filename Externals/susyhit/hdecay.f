C ==================================================================
C ================= PROGRAM HDECAY: COMMENTS =======================
C ==================================================================
C The program version HDECAY 3.4 has been modified to be linked 
C to SDECAY in order to complete the program package
C SUSYHIT - SU(spect)-S(deca)Y-H(decay)-I(n)Terface
C
C
C         Last modification on December 28 2008 by M.Muehlleitner
C         Last modification on November 24th 2008 by M.S.
C ==================================================================
C ================= PROGRAM HDECAY: COMMENTS =======================
C ==================================================================
C
C                       *****************
C                       * VERSION 3.4   *
C                       *****************
C
C
C  This program calculates the total decay widths and the branching 
C  ratios of the C Standard Model Higgs boson (HSM) as well as those 
C  of the neutral (HL= the light CP-even, HH= the heavy CP-even, HA= 
C  the pseudoscalar) and the charged (HC) Higgs bosons of the Minimal
C  Supersymmetric extension of the Standard Model (MSSM). It includes:
C
C - All the decay channels which are kinematically allowed and which
C   have branching ratios larger than 10**(-4). 
C
C - All QCD corrections to the fermionic and gluonic decay modes.
C   Most of these corrections are mapped into running masses in a
C   consistent way with some freedom for including high order terms. 
C
C - Below--threshold three--body decays with off--shell top quarks
C   or ONE off-shell gauge boson, as well as some decays with one
C   off-shell Higgs boson in the MSSM. 
C
C - Double off-shell decays: HSM,HL,HH --> W*W*,Z*Z* -->4 fermions,
C   which could be important for Higgs masses close to MW or MZ.
C
C - In the MSSM, the radiative corrections with full squark mixing and 
C   uses the RG improved values of Higgs masses and couplings with the 
C   main NLO corrections implemented (based on M.Carena, M. Quiros and
C   C.E.M. Wagner, Nucl. Phys. B461 (1996) 407, hep-ph/9508343). 
C
C - In the MSSM, all the decays into CHARGINOS, NEUTRALINOS, SLEPTONS 
C   and SQUARKS (with mixing in the stop and sbottom sectors). 
C
C - Chargino, slepton and squark loops in the 2 photon decays and squark
C   loops in the gluonic decays (including QCD corrections). 
C
C  ===================================================================
C  This program has been written by A.Djouadi, J.Kalinowski, M.
C  Muehlleitner and M.Spira. For details on how to use the program see:
C  Comp. Phys. Commun. 108 (1998) 56, hep-ph/9704448. For any question,
C  comment, suggestion or complaint, please contact us at:
C          Abdelhak.Djouadi@th.u-psud.fr
C          kalino@fuw.edu.pl
C          muehl@lapp.in2p3.fr
C          Michael.Spira@psi.ch


C ================ IT USES AS INPUT PARAMETERS:
C
C   SLHAIN: =0: READ FROM hdecay.in
c change susyhit
C           =1: READ SUSY LES HOUCHES ACCORD INPUT (slhaspectrum.in)
c end change susyhit
C
C  SLHAOUT: =0: WRITE BR TABLES
c change susyhit
C           =1: WRITE SUSY LES HOUCHES ACCORD OUTPUT (hdecay_slha.out)
c end change susyhit
C
C   IHIGGS: =0: CALCULATE BRANCHING RATIOS OF SM HIGGS BOSON
C           =1: CALCULATE BRANCHING RATIOS OF MSSM h BOSON
C           =2: CALCULATE BRANCHING RATIOS OF MSSM H BOSON
C           =3: CALCULATE BRANCHING RATIOS OF MSSM A BOSON
C           =4: CALCULATE BRANCHING RATIOS OF MSSM H+ BOSON
C           =5: CALCULATE BRANCHING RATIOS OF ALL MSSM HIGGS BOSONS
C
C   IMODEL: USE SPECIFIC SUBROUTINE FOR MSSM HIGSS MASSES AND COUPLINGS
C           =1: CARENA ET AL., NUCL. PHYS. B461 (1996) 407 (SUBHPOLE)
C           =2: CARENA ET AL., PHYS. LETT. B355 (1995) 209 (SUBH)
C           =3: HABER ET AL.
C           =4: HEINEMEYER ET AL., HEP-PH/0002213 (FEYNHIGGSFAST1.2.2)
C
C TGBET:    TAN(BETA) FOR MSSM
C MABEG:    START VALUE OF M_A FOR MSSM AND M_H FOR SM
C MAEND:    END VALUE OF M_A FOR MSSM AND M_H FOR SM
C NMA:      NUMBER OF ITERATIONS FOR M_A
C ALS(MZ):  VALUE FOR ALPHA_S(M_Z)
C MSBAR(1): MSBAR MASS OF STRANGE QUARK AT SCALE Q=1 GEV
C MC:       CHARM POLE MASS
C MB:       BOTTOM POLE MASS
C MT:       TOP POLE MASS
C MTAU:     TAU MASS
C MMUON:    MUON MASS
C ALPH:     INVERSE QED COUPLING
C GF:       FERMI CONSTANT
C GAMW:     W WIDTH
C GAMZ:     Z WIDTH
C MZ:       Z MASS
C MW:       W MASS
C VUS:      CKM PARAMETER V_US
C VCB:      CKM PARAMETER V_CB
C VUB/VCB:  RATIO V_UB/V_CB
C 1ST AND 2ND GENERATION:
C MSL1:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED SLEPTONS 
C MER1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SLEPTONS 
C MQL1:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED SUPS
C MUR1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SUPS
C MDR1:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SDOWNS 
C 3RD GENERATION:
C MSL:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED STAUS 
C MER:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED STAUS 
C MSQ:      SUSY BREAKING MASS PARAMETERS OF LEFT HANDED STOPS
C MUR:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED STOPS
C MDR:      SUSY BREAKING MASS PARAMETERS OF RIGHT HANDED SBOTTOMS 
C AL:       STAU TRILINEAR SOFT BREAKING TERMS 
C AU:       STOP TRILINEAR SOFT BREAKING TERMS
C AD:       SBOTTOM TRILINEAR SOFT BREAKING TERMS
C MU:       SUSY HIGGS MASS PARAMETER
C M2:       GAUGINO MASS PARAMETER
C MGLUINO:  GLUINO MASS
C
C NNLO (M): =0: USE O(ALPHA_S) FORMULA FOR POLE MASS --> MSBAR MASS
C           =1: USE O(ALPHA_S**2) FORMULA FOR POLE MASS --> MSBAR MASS
C
C ON-SHELL: =0: INCLUDE OFF_SHELL DECAYS H,A --> T*T*, A --> Z*H,
C               H --> W*H+,Z*A, H+ --> W*A, W*H, T*B
C           =1: EXCLUDE THE OFF-SHELL DECAYS ABOVE
C
C ON-SH-WZ: =0: INCLUDE DOUBLE OFF-SHELL PAIR DECAYS PHI --> W*W*,Z*Z*
C           =1: INCLUDE ONLY SINGLE OFF-SHELL DECAYS PHI --> W*W,Z*Z
C
C IPOLE:    =0 COMPUTES RUNNING HIGGS MASSES (FASTER) 
C           =1 COMPUTES POLE HIGGS MASSES 
C
C OFF-SUSY: =0: INCLUDE DECAYS (AND LOOPS) INTO SUPERSYMMETRIC PARTICLES
C           =1: EXCLUDE DECAYS (AND LOOPS) INTO SUPERSYMMETRIC PARTICLES
C
C INDIDEC:  =0: PRINT OUT SUMS OF CHARGINO/NEUTRALINO/SFERMION DECAYS
C           =1: PRINT OUT INDIVIDUAL CHARGINO/NEUTRALINO/SFERMION DECAYS
C
C NF-GG:    NUMBER OF LIGHT FLAVORS INCLUDED IN THE GLUONIC DECAYS 
C            PHI --> GG* --> GQQ (3,4 OR 5)
C           
C IGOLD:    =0: EXCLUDE DECAYS INTO GRAVITINO + GAUGINO
C           =1: INCLUDE DECAYS INTO GRAVITINO + GAUGINO
C
C MPLANCK:  PLANCK MASS FOR DECAYS INTO GRAVITINO + GAUGINO
C MGOLD:    GRAVITINO MASS FOR DECAYS INTO GRAVITINO + GAUGINO
C
C =======================================================================
C ============== BEGINNING OF THE MAIN PROGRAM ==========================
C =======================================================================
C
c change susyhit
      subroutine HDECAY
c end change susyhit
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/FLAGS_HDEC/INDIDEC
      COMMON/SLHA_vals_HDEC/islhai,islhao

      CALL READ_HDEC(TGBET,AMABEG,AMAEND,NMA)
      if(islhao.ne.1) then
         CALL HEAD_HDEC(TGBET,AMABEG)
      endif

      DO 9999 II=1,NMA
       IF(NMA.NE.1)THEN
        AMAR = AMABEG + (AMAEND-AMABEG)/(NMA-1D0)*(II-1D0)
       ELSE
        AMAR = AMABEG
       ENDIF
       AMSM = AMAR
       AMA = AMAR
       CALL HDEC(TGBET)
c change susyhit
c       CALL WRITE_HDEC(TGBET)
c end change susyhit
 9999  CONTINUE

      CALL CLOSE_HDEC

c change susyhit
c      STOP
c end change susyhit
      END

      SUBROUTINE READ_HDEC(TGBET,AMABEG,AMAEND,NMA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60,
     .          ninlha=22)
      double precision minval(1:20),smval(1:20),massval(1:50),
     .                 nmixval(4,4),umixval(2,2),vmixval(2,2),
     .                 stopmixval(2,2),sbotmixval(2,2),staumixval(2,2),
     .                 hmixval(1:10),gaugeval(1:3),msoftval(1:100),
     .                 auval(3,3),adval(3,3),aeval(3,3),yuval(3,3),
     .                 ydval(3,3),yeval(3,3),qvalue(1:20),
     .                 extval(0:100),m_softval(1:100)
      double precision slhaneut(1:4),slhaxneut(1:4),slhachar(1:2),
     .                 slhaxchar(1:2),
     .                 slhau(2,2),slhav(2,2),slhaz(4,4),
     .                 slhast(2),slhasb(2),slhasu(2),slhasd(2),
     .                 slhase(2),slhasl(2),slhasn(2),slhasnl(2),
     .                 warning(1:10)
      integer   imod(1:2)
      integer check(1:22)
      double precision mbmsbar,mbl,mbu
      character spinfo1*100,spinfo2*100,modselval*100,mincom(1:20)*20,
     .          extcom(0:100)*20,softcom(1:100)*20,hmixcom(1:10)*20,
     .          m_softcom(1:100)*20
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
c -------------- common block given by read_leshouches ------------ c
      COMMON/SLHA_leshouches1_HDEC/spinfo1,spinfo2,modselval,mincom,
     .                             extcom,softcom,hmixcom
      COMMON/SLHA_leshouches2_HDEC/minval,extval,smval,massval,nmixval,
     .                      umixval,vmixval,stopmixval,sbotmixval,
     .                      staumixval,hmixval,gaugeval,msoftval,auval,
     .                      adval,aeval,yuval,ydval,yeval,alphaval,
     .                      qvalue,imod
c -------------- common blocks needed in HDECAY subroutines ---------- c
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/SLHA_m1_HDEC/am1
      COMMON/SLHA_gaug_HDEC/slhaneut,slhaxneut,slhachar,slhau,slhav,
     .                      slhaz,slhaxchar
      COMMON/SLHA_sfer_HDEC/slhast,slhasb,slhasu,slhasd,slhase,slhasl,
     .                 slhasn,slhasnl,slhacot,slhasit,slhacob,slhasib,
     .                 slhacol,slhasil
      COMMON/SLHA_hmass_HDEC/slhaml,slhamh,slhamc,slha_alpha
      COMMON/SLHAVAL_HDEC/g1ew,g2ew
      COMMON/SLHA_checkval_HDEC/check
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .              ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD
      COMMON/FLAGS_HDEC/INDIDEC
c change susyhit
      COMMON/SUSYHITIN/flagshsin,amsin,amcin,ammuonin,alphin,gamwin,
     .                 gamzin,vusin,vcbin,rvubin
c end change susyhit

      unlikely = -123456789D0

      PI = 4*DATAN(1D0)

c change susyhit
      if(flagshsin.gt.2.D0) then
         OPEN(NPAR,FILE='br.input')
      endif

c -- needed parameters, read in via sdecay from susyhit.in, which are not 
c -- given by the slhaspectrum.in file --

      ams    = amsin
      amc    = amcin
      ammuon = ammuonin
      alph   = alphin
      gamw   = gamwin
      gamz   = gamzin
      vus    = vusin
      vcb    = vcbin
      rvub   = rvubin

c -- parameters to be defined which are not given by the --
c -- slhaspectrum.in file --

      islhai  = 1
      islhao  = 1
      ihiggs  = 5
      imodel  = 1
      nnlo    = 1
      ionsh   = 0
      ionwz   = 0
      ipole   = 0
      iofsusy = 0
      indidec = 0
      nfgg    = 5
      igold   = 0
      axmpl   = 2.4d18
      axmgd   = 1.d-13

c -- read in the parameters from slhaspectrum.in --
c end change susyhit

c -- initialization of the check array --
      do i1=1,22,1
         check(i1) = 0
      end do

      if(islhai.eq.1) then
c change susyhit
         open(ninlha,file='slhaspectrum.in')
c end change susyhit
         call SLHA_read_leshouches_HDEC(ninlha)

c -- G_F --
         if(smval(2).ne.0.D0) then
            GF = smval(2)
         endif
c -- the strong coupling constant alphas_MSbar at the scale MZ --
         if(smval(3).ne.0.D0) then
            alsmz = smval(3)
         endif
         alphasmzms = alsmz
c -- Z pole mass --
         if(smval(4).ne.0.D0) then
            AMZ = smval(4)
         endif
c -- W pole mass --
         if(massval(1).ne.0.D0) then
            AMW = massval(1)
         endif
c -- the MSbar couplings g1,g2 at the scale Q --
         if(gaugeval(1).ne.0.D0) then
            g1ew  = gaugeval(1)
         endif
         if(gaugeval(2).ne.0.D0) then
            g2ew  = gaugeval(2)*(1-gaugeval(2)**2/96/pi**2*2)
         endif
         cw2calc = amw**2/amz**2
         sw2calc = 1-cw2calc
         cwcalc  = dsqrt(cw2calc)
         swcalc  = dsqrt(sw2calc)
c -- in case the gauge couplings are not given at the scale Q --
         if(gaugeval(1).eq.0.D0.or.gaugeval(2).eq.0.D0) then
c -- v at the scale Q --         
          if(smval(2).eq.0.D0.and.hmixval(3).ne.unlikely) then
           vewsb = hmixval(3)
           gf = 1/dsqrt(2.D0)/vewsb**2
          else
           vewsb = 1.D0/dsqrt(dsqrt(2.D0)*gf)
          endif
          g2ew = 2*amw/vewsb
          g1ew = g2ew*swcalc/cwcalc
         endif

c -- neutralino and chargino masses --      

         slhaneut(1) =dabs(massval(28))
         slhaneut(2) =dabs(massval(29))
         slhaneut(3) =dabs(massval(30))
         slhaneut(4) =dabs(massval(31))
         slhaxneut(1)=massval(28)
         slhaxneut(2)=massval(29)
         slhaxneut(3)=massval(30)
         slhaxneut(4)=massval(31)
         slhachar(1) =dabs(massval(32))
         slhachar(2) =dabs(massval(33))
         slhaxchar(1)=massval(32)
         slhaxchar(2)=massval(33)

c -- the chargino and neutralino mixing matrix elements --
      do i=1,2,1
         do j=1,2,1
            slhau(i,j)=umixval(i,j)
            slhav(i,j)=vmixval(i,j)
         end do
      end do
      do i=1,4,1
         do j=1,4,1
            slhaz(i,j)=nmixval(i,j)
         end do
      end do

c -- sfermion masses --

      slhast(1) = massval(16)
      slhast(2) = massval(17)
      slhasb(1) = massval(14)
      slhasb(2) = massval(15)

      slhasu(1) = massval(8)
      slhasu(2) = massval(9)
      slhasd(1) = massval(10)
      slhasd(2) = massval(11)

      slhase(1) = massval(18)
      slhase(2) = massval(19)
      slhasl(1) = massval(24)
      slhasl(2) = massval(25)

      slhasn(1) = massval(20)
      slhasn(2) = 1.D15
      slhasnl(1) = massval(26)
      slhasnl(2) = 1.D15

c -- the sfermion mixing angles --

      slhacot=stopmixval(1,1)
      slhasit=stopmixval(1,2)

      slhacob=sbotmixval(1,1)
      slhasib=sbotmixval(1,2)

      slhacol=staumixval(1,1)
      slhasil=staumixval(1,2)

c -- the gluino mass --

      AMGLUINO = massval(27)

c -- the Higgs masses --

      slhaml = massval(2)
      slhamh = massval(3)
      slhamc = massval(5)

      if(massval(4).ne.0.D0) then
         slhama = massval(4)
      elseif(extval(26).ne.unlikely) then
         slhama = extval(26)
      elseif(extval(24).ne.unlikely) then
         slhama = dsqrt(extval(24))
      endif

      amabeg = slhama
      amaend = slhama
      nma    = 1

c -- the MSSM mixing angle alpha in the Higgs sector --
c -- Attention: It might be that alphaval is not the DRbar value at
c -- the scale Q.
      slha_alpha = alphaval

c -- the fermion pole masses --

      if(smval(6).ne.0.D0) then
         AMT = smval(6)
      endif
      if(smval(7).ne.0.D0) then
         AMTAU = smval(7)
      endif

c -- the mass mb(mb)_MSbar --
      if(smval(5).ne.0.D0) then
         mbmsbar = smval(5)
      endif

      fmt = amt
      fmtau = amtau
      fms = ams
      fmc = amc
c -- calculation of the mb pole mass from mb(mb)_MSbar --
      if(smval(5).ne.0.D0) then
       del = 1.d-10
       mbl = mbmsbar
       mbu = 2*mbmsbar
       fmb = (mbl+mbu)/2
       amsb = ams
       amc0=amc
       amt0=amt
       acc=1.d-8
       nloop=2
11     amb=fmb
       amb0=amb
       xlambda=xitla_hdec(nloop,alsmz,acc)
       n0=5
       call alsini_hdec(acc)
c      xmb = runm_hdec(fmb,5)
       xmb = runm_hdec(mbmsbar,5)
       if(xmb.eq.mbmsbar)then
        mbl = fmb
        mbu = fmb
       elseif(xmb.gt.mbmsbar)then
        mbu = fmb
       else
        mbl = fmb
       endif
       fmb = (mbl+mbu)/2
       if(dabs(xmb/mbmsbar-1).gt.del) goto 11
      endif
      amb = fmb

c -- DRbar value of tanbeta at the scale Q --

      if(hmixval(2).ne.unlikely) then
         TGBET = hmixval(2)
      endif

c -- If no DRbar value at the scale Q has been given for tanbeta --

      if(hmixval(2).eq.unlikely) then
         if(extval(25).ne.0.D0.and.extval(25).ne.unlikely) then
            TGBET = extval(25)
         elseif(minval(3).ne.0.D0.and.minval(3).ne.unlikely) then
            TGBET = minval(3)
         endif
      endif

c -- The soft SUSY breaking parameters: DRbar values at the scale Q --

      do i=1,100,1
         m_softval(i) = unlikely
      end do

      do i=1,99,1
         if(msoftval(i).ne.unlikely) then
            m_softval(i) = msoftval(i)
            m_softcom(i) = softcom(i)
         elseif(extval(i).ne.unlikely) then
            m_softval(i) = extval(i)
            m_softcom(i) = extcom(i)
         endif
      end do

      AU=auval(3,3)
      AD=adval(3,3)
      AL=aeval(3,3)

c The mixing parameter mu in the MS_bar scheme
      amudrbar = hmixval(1)
      if(amudrbar.ne.unlikely)then
       AMU = amudrbar*(1.D0+g1ew**2/16.D0/pi**2*3.D0/5.D0+
     .                 g2ew**2/16.D0/pi**2*3.D0/4.D0)
      endif


c The soft SUSY breaking parameters M1, M2 in the MS_bar scheme 
      am1msbar = m_softval(1)*(1.D0+g1ew**2/16.D0/pi**2*0.D0)
      am2msbar = m_softval(2)*(1.D0+g2ew**2/16.D0/pi**2*2.D0)

      if(am1msbar.ne.0.d0)am1 = am1msbar
      if(am2msbar.ne.0.d0)AM2 = am2msbar

      if(m_softval(31).ne.unlikely.and.m_softval(32).ne.unlikely) then
         if(m_softval(31).ne.0.D0.and.m_softval(32).ne.0.D0) then
            AMEL1 = (m_softval(31)+m_softval(32))/2.D0
         elseif(m_softval(31).ne.0.D0) then
            AMEL1 = m_softval(31)
         elseif(m_softval(32).ne.0.D0) then
            AMEL1 = m_softval(32)
         else
            AMEL1 = 0.D0
         endif
      elseif(m_softval(31).ne.unlikely) then
         AMEL1 = m_softval(31)
      elseif(m_softval(32).ne.unlikely) then
         AMEL1 = m_softval(32)
      endif

      if(m_softval(34).ne.unlikely.and.m_softval(35).ne.unlikely) then
         if(m_softval(34).ne.0.D0.and.m_softval(35).ne.0.D0) then
            AMER1 = (m_softval(34)+m_softval(35))/2.D0
         elseif(m_softval(34).ne.0.D0) then
            AMER1 = m_softval(34)
         elseif(m_softval(35).ne.0.D0) then
            AMER1 = m_softval(35)
         else
            AMER1 = 0.D0
         endif
      elseif(m_softval(34).ne.unlikely) then
         AMER1 = m_softval(34)
      elseif(m_softval(35).ne.unlikely) then
         AMER1 = m_softval(35)
      endif

      if(m_softval(41).ne.unlikely.and.m_softval(42).ne.unlikely) then
         if(m_softval(41).ne.0.D0.and.m_softval(42).ne.0.D0) then
            AMQL1 = (m_softval(41)+m_softval(42))/2.D0
         elseif(m_softval(41).ne.0.D0) then
            AMQL1 = m_softval(41)
         elseif(m_softval(42).ne.0.D0) then
            AMQL1 = m_softval(42)
         else
            AMQL1 = 0.D0
         endif
      elseif(m_softval(41).ne.unlikely) then
         AMQL1 = m_softval(41)
      elseif(m_softval(42).ne.unlikely) then
         AMQL1 = m_softval(42)
      endif

      if(m_softval(44).ne.unlikely.and.m_softval(45).ne.unlikely) then
         if(m_softval(44).ne.0.D0.and.m_softval(45).ne.0.D0) then
            AMUR1 = (m_softval(44)+m_softval(45))/2.D0
         elseif(m_softval(44).ne.0.D0) then
            AMUR1 = m_softval(44)
         elseif(m_softval(45).ne.0.D0) then
            AMUR1 = m_softval(45)
         else
            AMUR1 = 0.D0
         endif
      elseif(m_softval(44).ne.unlikely) then
         AMUR1 = m_softval(44)
      elseif(m_softval(45).ne.unlikely) then
         AMUR1 = m_softval(45)
      endif

      if(m_softval(47).ne.unlikely.and.m_softval(48).ne.unlikely) then
         if(m_softval(47).ne.0.D0.and.m_softval(48).ne.0.D0) then
            AMDR1 = (m_softval(47)+m_softval(48))/2.D0
         elseif(m_softval(47).ne.0.D0) then
            AMDR1 = m_softval(47)
         elseif(m_softval(48).ne.0.D0) then
            AMDR1 = m_softval(48)
         else
            AMDR1 = 0.D0
         endif
      elseif(m_softval(47).ne.unlikely) then
         AMDR1 = m_softval(47)
      elseif(m_softval(48).ne.unlikely) then
         AMDR1 = m_softval(48)
      endif

      if(m_softval(33).ne.unlikely) then
         AMEL = m_softval(33)
      endif
      if(m_softval(36).ne.unlikely) then
         AMER = m_softval(36)
      endif
      if(m_softval(43).ne.unlikely) then
         AMSQ = m_softval(43)
      endif
      if(m_softval(46).ne.unlikely) then
         AMUR = m_softval(46)
      endif
      if(m_softval(49).ne.unlikely) then
         AMDR = m_softval(49)
      endif

      endif

      IF(IMODEL.EQ.3)THEN
       WRITE(6,*)'MU (UP TO THE SIGN) WILL BE IDENTIFIED WITH M_SQ...'
      ENDIF

      B = DATAN(TGBET)
      AMGLU = AMGLUINO

      VUB=RVUB*VCB
      ALPH=1.D0/ALPH
      AMSB = AMS

      AMC0=AMC
      AMB0=AMB
      AMT0=AMT
      ACC=1.D-8
      NLOOP=2
      XLAMBDA=XITLA_HDEC(NLOOP,ALSMZ,ACC)
      N0=5
      CALL ALSINI_HDEC(ACC)

C--INITIALIZE COEFFICIENTS FOR POLYLOGARITHMS
      NBER = 18
      CALL BERNINI_HDEC(NBER)

C--CHECK NFGG
      IF(NFGG.GT.5.OR.NFGG.LT.3)THEN
       WRITE(6,*)'NF-GG NOT VALID. TAKING THE DEFAULT NF-GG = 3....'
       NFGG = 3
      ENDIF

100   FORMAT(10X,G30.20)
101   FORMAT(10X,I30)

C--WRITE THE INPUT PARAMTERS TO A DATA-FILE

c change susyhit
      if(flagshsin.gt.2.D0) then
c end change susyhit
      WRITE(NPAR,8)'SLHAIN   = ',ISLHAI
      WRITE(NPAR,8)'SLHAOUT  = ',ISLHAO
      WRITE(NPAR,8)'HIGGS    = ',IHIGGS
      WRITE(NPAR,8)'MODEL    = ',IMODEL
      WRITE(NPAR,9)'TGBET    = ',TGBET
      WRITE(NPAR,9)'MABEG    = ',AMABEG
      WRITE(NPAR,9)'MAEND    = ',AMAEND
      WRITE(NPAR,7)'NMA      = ',NMA
      WRITE(NPAR,9)'ALS(MZ)  = ',ALSMZ
      WRITE(NPAR,9)'MSBAR(1) = ',AMS
      WRITE(NPAR,9)'MC       = ',AMC
      WRITE(NPAR,9)'MB       = ',AMB
      WRITE(NPAR,9)'MT       = ',AMT
      WRITE(NPAR,9)'MTAU     = ',AMTAU
      WRITE(NPAR,9)'MMUON    = ',AMMUON
      WRITE(NPAR,9)'ALPH     = ',1.D0/ALPH
      WRITE(NPAR,9)'GF       = ',GF
      WRITE(NPAR,9)'GAMW     = ',GAMW
      WRITE(NPAR,9)'GAMZ     = ',GAMZ
      WRITE(NPAR,9)'MZ       = ',AMZ
      WRITE(NPAR,9)'MW       = ',AMW
      WRITE(NPAR,9)'VUS      = ',VUS
      WRITE(NPAR,9)'VCB      = ',VCB
      WRITE(NPAR,9)'VUB/VCB  = ',RVUB
      WRITE(NPAR,9)'MU       = ',AMU
      WRITE(NPAR,9)'M2       = ',AM2
      WRITE(NPAR,9)'MGLUINO  = ',AMGLUINO
      WRITE(NPAR,9)'MSL1     = ',AMEL1
      WRITE(NPAR,9)'MER1     = ',AMER1
      WRITE(NPAR,9)'MQL1     = ',AMQL1
      WRITE(NPAR,9)'MUR1     = ',AMUR1
      WRITE(NPAR,9)'MDR1     = ',AMDR1
      WRITE(NPAR,9)'MSL      = ',AMEL
      WRITE(NPAR,9)'MER      = ',AMER
      WRITE(NPAR,9)'MSQ      = ',AMSQ
      WRITE(NPAR,9)'MUR      = ',AMUR
      WRITE(NPAR,9)'MDR      = ',AMDR
      WRITE(NPAR,9)'AL       = ',AL
      WRITE(NPAR,9)'AU       = ',AU
      WRITE(NPAR,9)'AD       = ',AD
      WRITE(NPAR,8)'NNLO (M) = ',NNLO
      WRITE(NPAR,8)'ON-SHELL = ',IONSH
      WRITE(NPAR,8)'ON-SH-WZ = ',IONWZ
      WRITE(NPAR,8)'IPOLE    = ',IPOLE 
      WRITE(NPAR,8)'OFF-SUSY = ',IOFSUSY
      WRITE(NPAR,8)'INDIDEC  = ',INDIDEC
      WRITE(NPAR,8)'NF-GG    = ',NFGG
      WRITE(NPAR,8)'IGOLD    = ',IGOLD
      WRITE(NPAR,9)'MPLANCK  = ',AXMPL
      WRITE(NPAR,9)'MGOLD    = ',AXMGD
C     WRITE(NPAR,9)'LAMBDA_5 = ',XLAMBDA

      CLOSE(NPAR)
c change susyhit
      endif
c end change susyhit

7     FORMAT(A11,I7)
8     FORMAT(A11,I4)
9     FORMAT(A11,G15.8)

c change susyhit
c      CLOSE(NI)
c end change susyhit

      RETURN
      END

      SUBROUTINE HEAD_HDEC(TGBET,AMABEG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60)
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .              ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD
      COMMON/FLAGS_HDEC/INDIDEC

      PI = 4*DATAN(1D0)

      IF(IHIGGS.EQ.0) THEN
       OPEN(NSA,FILE='br.sm1')
       OPEN(NSB,FILE='br.sm2')
      ENDIF
      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
       OPEN(NLA,FILE='br.l1')
       OPEN(NLB,FILE='br.l2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYL,FILE='br.ls')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYLA,FILE='br.ls1')
        OPEN(NSUSYLB,FILE='br.ls2')
        OPEN(NSUSYLC,FILE='br.ls3')
        OPEN(NSUSYLD,FILE='br.ls4')
        OPEN(NSUSYLE,FILE='br.ls5')
        OPEN(NSUSYLF,FILE='br.ls6')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       OPEN(NHA,FILE='br.h1')
       OPEN(NHB,FILE='br.h2')
       OPEN(NHC,FILE='br.h3')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYH,FILE='br.hs')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYHA,FILE='br.hs1')
        OPEN(NSUSYHB,FILE='br.hs2')
        OPEN(NSUSYHC,FILE='br.hs3')
        OPEN(NSUSYHD,FILE='br.hs4')
        OPEN(NSUSYHE,FILE='br.hs5')
        OPEN(NSUSYHF,FILE='br.hs6')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
       OPEN(NAA,FILE='br.a1')
       OPEN(NAB,FILE='br.a2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYA,FILE='br.as')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYAA,FILE='br.as1')
        OPEN(NSUSYAB,FILE='br.as2')
        OPEN(NSUSYAC,FILE='br.as3')
        OPEN(NSUSYAD,FILE='br.as4')
       ENDIF
      ENDIF
      ENDIF
      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
       OPEN(NCA,FILE='br.c1')
       OPEN(NCB,FILE='br.c2')
      IF(IOFSUSY.EQ.0)THEN 
       OPEN(NSUSYC,FILE='br.cs')
       IF(INDIDEC.NE.0)THEN 
        OPEN(NSUSYCA,FILE='br.cs1')
        OPEN(NSUSYCB,FILE='br.cs2')
        OPEN(NSUSYCC,FILE='br.cs3')
        OPEN(NSUSYCD,FILE='br.cs4')
       ENDIF
      ENDIF
      ENDIF

C--SETUP THE HEADS OF THE TABLES IN THE DATA-FILES

      IF(IHIGGS.EQ.0) THEN
      WRITE(NSA,70)'MHSM  ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NSA,69)
      WRITE(NSA,*)
      WRITE(NSB,70)'MHSM  ','GG ','GAM GAM','Z GAM ','WW ','ZZ ','WIDTH'
      WRITE(NSB,69)
      WRITE(NSB,*)
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
      WRITE(NLA,70)'MHL   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NLA,69)
      WRITE(NLA,*)
      WRITE(NLB,70)'MHL   ','GG ','GAM GAM','Z GAM ','WW ','ZZ ','WIDTH'
      WRITE(NLB,69)
      WRITE(NLB,*)
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
      WRITE(NHA,70)'MHH   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NHA,69)
      WRITE(NHA,*)
      WRITE(NHB,72)'MHH   ','GG ','GAM GAM','Z GAM ','WW ','ZZ '
      WRITE(NHB,69)
      WRITE(NHB,*)
      WRITE(NHC,72)'MHH   ','hh ','AA ','Z A ','W+- H-+','WIDTH '
      WRITE(NHC,69)
      WRITE(NHC,*)
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
      WRITE(NAA,70)'MHA   ','BB   ','TAU TAU','MU MU ','SS ','CC ','TT '
      WRITE(NAA,69)
      WRITE(NAA,*)
      WRITE(NAB,72)'MHA   ','GG ','GAM GAM','Z GAM ','Z HL ','WIDTH '
      WRITE(NAB,69)
      WRITE(NAB,*)
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
      WRITE(NCA,70)'MHC   ','BC   ','TAU NU ','MU NU ','SU ','CS ','TB '
      WRITE(NCA,69)
      WRITE(NCA,*)
      WRITE(NCB,70)'MHC   ','BU   ','hW ','AW ','WIDTH '
      WRITE(NCB,69)
      WRITE(NCB,*)
      ENDIF

69    FORMAT(79('_'))
70    FORMAT(A9,6(1X,A10))
71    FORMAT(A9,4(1X,A10))
72    FORMAT(A9,5(1X,A10))
73    FORMAT(A9,3(1X,A10))

      AMAR = AMABEG
      AMSM = AMAR
      AMA = AMAR

      IF(IHIGGS.NE.0)THEN 
C *******************************  SUSY OUTPUT 

       CALL GAUGINO_HDEC(AMU,AM2,B,A,GMC,GMN,XMN,AC1,AC2,AC3,
     .              AN1,AN2,AN3,ACNL,ACNR,AGDL,AGDA,AGDH,AGDC)
      TSC = (AMSQ+AMUR+AMDR)/3
      BSC = (AMSQ+AMUR+AMDR)/3
      AMT00 = AMT0
      AMT0 = 3.D8
      CALL SFERMION_HDEC(TSC,BSC,AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .               GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN, 
     .               GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .               GAEE,GATT,GABB,GCEN,GCTB)
      AMT0 = AMT00
      CALL SUSYCP_HDEC(TGBET)
c     write(6,*)'MZ,MW,SW2,alpha: ',AMZ,AMW,1-AMW**2/AMZ**2,A
c     write(6,*)'tan(beta),Ab,mu: ',TGBET,AD,AMU
c     write(6,*)'M_A, M_h, M_H, M_H+: ',AMA,AML,AMH,AMCH
c     write(6,*)'Lambda_hhh/Lambda_SM: ',GLLL/AML**2*AMZ**2/3
c     write(6,*)
c     write(96,*)ama,aml,amh,amch
c     write(97,*)glb,glt,glvv,ghb,ght,ghvv

      IF(IOFSUSY.EQ.0)THEN
C--WRITE THE GAUGINO MASSES/ TB, MU AND M2 IN THE SUSY DATA-FILE
C--WRITE THE SFERMION MASSES/ SUSY MASSES AND COUPLINGS IN SUSY DATA-FILE
C 
       IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYL,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYL,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYL,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYL,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYL,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYL,*)
       WRITE(NSUSYL,*)'   MHL        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS  GRAVITINO+GAUGINO'
       WRITE(NSUSYL,69)
       WRITE(NSUSYL,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYLA,73)'MHL   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYLA,69)
         WRITE(NSUSYLA,*)
         WRITE(NSUSYLB,71)'MHL   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYLB,69)
         WRITE(NSUSYLB,*)
         WRITE(NSUSYLC,70)'MHL   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYLC,69)
         WRITE(NSUSYLC,*)
         WRITE(NSUSYLD,*)'   MHL        SNL SNL    SEL SEL    '//
     .   'SER SER    STA1 STA1  STA1 STA2  STA2 STA2' 
         WRITE(NSUSYLD,69)
         WRITE(NSUSYLD,*)
         WRITE(NSUSYLE,*)'   MHL        SUL SUL    SUR SUR    '//
     .   'SDL SDL    SDR SDR'
         WRITE(NSUSYLE,69)
         WRITE(NSUSYLE,*)
         WRITE(NSUSYLF,*)'   MHL        SB1 SB1    SB1 SB2    '//
     .   'SB2 SB2    ST1 ST1    ST1 ST2    ST2 ST2'
         WRITE(NSUSYLF,69)
         WRITE(NSUSYLF,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYH,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYH,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYH,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYH,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYH,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYH,*)
       WRITE(NSUSYH,*)'   MHH        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS  GRAVITINO+GAUGINO'
       WRITE(NSUSYH,69)
       WRITE(NSUSYH,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYHA,73)'MHH   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYHA,69)
         WRITE(NSUSYHA,*)
         WRITE(NSUSYHB,71)'MHH   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYHB,69)
         WRITE(NSUSYHB,*)
         WRITE(NSUSYHC,70)'MHH   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYHC,69)
         WRITE(NSUSYHC,*)
         WRITE(NSUSYHD,*)'   MHH        SNL SNL    SEL SEL    '//
     .   'SER SER    STA1 STA1  STA1 STA2  STA2 STA2' 
         WRITE(NSUSYHD,69)
         WRITE(NSUSYHD,*)
         WRITE(NSUSYHE,*)'   MHH        SUL SUL    SUR SUR    '//
     .   'SDL SDL    SDR SDR'
         WRITE(NSUSYHE,69)
         WRITE(NSUSYHE,*)
         WRITE(NSUSYHF,*)'   MHH        SB1 SB1    SB1 SB2    '//
     .   'SB2 SB2    ST1 ST1    ST1 ST2    ST2 ST2'
         WRITE(NSUSYHF,69)
         WRITE(NSUSYHF,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYA,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYA,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYA,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYA,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYA,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYA,*)
       WRITE(NSUSYA,*)'   MHA        CHARGINOS  NEUTRALS   '//
     . 'SLEPTONS   SQUARKS  GRAVITINO+GAUGINO'
       WRITE(NSUSYA,69)
       WRITE(NSUSYA,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYAA,73)'MHA   ','C1 C1 ','C2 C2 ','C1 C2 '
         WRITE(NSUSYAA,69)
         WRITE(NSUSYAA,*)
         WRITE(NSUSYAB,71)'MHA   ','N1 N1 ','N2 N2 ','N3 N3 ','N4 N4 '
         WRITE(NSUSYAB,69)
         WRITE(NSUSYAB,*)
         WRITE(NSUSYAC,70)'MHA   ','N1 N2 ','N1 N3 ','N1 N4 ','N2 N3 ',
     .                    'N2 N4 ','N3 N4 '
         WRITE(NSUSYAC,69)
         WRITE(NSUSYAC,*)
         WRITE(NSUSYAD,*)
         WRITE(NSUSYAD,*)'   MHA        STA1 STA2  SB1 SB2    ST1 ST2'
         WRITE(NSUSYAD,69)
         WRITE(NSUSYAD,*)
        ENDIF
       ENDIF

       IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
       WRITE(NSUSYC,347) TGBET,AM2,AMU,AMSQ
       WRITE(NSUSYC,348) GMC(1),GMC(2),GMN(1),GMN(2),GMN(3),GMN(4)
       WRITE(NSUSYC,349) GMST(1),GMST(2),GMSU(1),GMSU(2)
       WRITE(NSUSYC,350) GMSB(1),GMSB(2),GMSD(1),GMSD(2)
       WRITE(NSUSYC,351) GMSL(1),GMSL(2),GMSE(1),GMSE(2),GMSN(1)
       WRITE(NSUSYC,*)
       WRITE(NSUSYC,*)'   MHC        CHARG/NEU  SLEPTONS   SQUARKS',
     .                '  GRAVITINO+GAUGINO'
       WRITE(NSUSYC,69)
       WRITE(NSUSYC,*)
        IF(INDIDEC.NE.0)THEN
         WRITE(NSUSYCA,70)'MHC   ','C1 N1 ','C1 N2 ','C1 N3 ','C1 N4 '
         WRITE(NSUSYCA,69)
         WRITE(NSUSYCA,*)
         WRITE(NSUSYCB,70)'MHC   ','C2 N1 ','C2 N2 ','C2 N3 ','C2 N4 '
         WRITE(NSUSYCB,69)
         WRITE(NSUSYCB,*)
         WRITE(NSUSYCC,*)'   MHC        SEL SNL    STAU1 SNL  STAU2 SNL'
         WRITE(NSUSYCC,69)
         WRITE(NSUSYCC,*)
         WRITE(NSUSYCD,*)'   MHC        SUL SDL    ST1 SB1    '//
     .   'ST1 SB2    ST2 SB1    ST2 SB2'
         WRITE(NSUSYCD,69)
         WRITE(NSUSYCD,*)
        ENDIF
       ENDIF

347    FORMAT('TB=',G12.6,1X,'M2=',G12.6,1X,'MU=',G12.6,1X,
     .        'MSQ=',G12.6)
348    FORMAT('C1=',F7.3,1X,'C2=',F8.3,1X,'N1=',F7.3,1X,'N2=',F7.3,1X,
     .        'N3=',F8.3,1X,'N4=',F8.3)
349    FORMAT('MST1=',G12.6,1X,'MST2=',G12.6,1X,
     .        'MSUL=',G12.6,1X,'MSUR=',G12.6) 
350    FORMAT('MSB1=',G12.6,1X,'MSB2=',G12.6,1X,
     .        'MSDL=',G12.6,1X,'MSDR=',G12.6) 
351    FORMAT('TAU1=',F8.3,1X,'TAU2=',F8.3,1X,'EL=',F8.3,1X,
     .        'ER=',F8.3,1X,'NL=',F8.3)
C
C
C **************************************************************
      ENDIF
      ENDIF

      RETURN
      END

      SUBROUTINE WRITE_HDEC(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60)
      parameter (nout=44)
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      dimension hlbrsn1(4,4),hhbrsn1(4,4),habrsn1(4,4)
      double precision minval(1:20),smval(1:20),massval(1:50),
     .                 nmixval(4,4),umixval(2,2),vmixval(2,2),
     .                 stopmixval(2,2),sbotmixval(2,2),staumixval(2,2),
     .                 hmixval(1:10),gaugeval(1:3),msoftval(1:100),
     .                 auval(3,3),adval(3,3),aeval(3,3),yuval(3,3),
     .                 ydval(3,3),yeval(3,3),qvalue(1:20),extval(0:100),
     .                 m_softval(1:100)
      double precision slhaneut(1:4),slhaxneut(1:4),slhachar(1:2),
     .                 slhau(2,2),slhav(2,2),slhaz(4,4),slhaxchar(1:2),
     .                 slhast(2),slhasb(2),slhasu(2),slhasd(2),
     .                 slhase(2),slhasl(2),slhasn(2),slhasnl(2),
     .                 warning(1:10)
      integer imod(1:2)
      integer check(1:22)
      double precision mbmsbar,mbl,mbu
      character spinfo1*100,spinfo2*100,modselval*100,mincom(1:20)*20,
     .          extcom(0:100)*20,softcom(1:100)*20,hmixcom(1:10)*20,
     .          m_softcom(1:100)*20
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .              ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD
      COMMON/FLAGS_HDEC/INDIDEC
c -------------- common block given by read_leshouches ------------ c
      COMMON/SLHA_leshouches1_HDEC/spinfo1,spinfo2,modselval,mincom,
     .                             extcom,softcom,hmixcom
      COMMON/SLHA_leshouches2_HDEC/minval,extval,smval,massval,nmixval,
     .                      umixval,vmixval,stopmixval,sbotmixval,
     .                      staumixval,hmixval,gaugeval,msoftval,auval,
     .                      adval,aeval,yuval,ydval,yeval,alphaval,
     .                      qvalue,imod
c -------------- common blocks needed in HDECAY subroutines ---------- c
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/SLHA_m1_HDEC/am1
      COMMON/SLHA_gaug_HDEC/slhaneut,slhaxneut,slhachar,slhau,slhav,
     .                      slhaz,slhaxchar
      COMMON/SLHA_sfer_HDEC/slhast,slhasb,slhasu,slhasd,slhase,slhasl,
     .                 slhasn,slhasnl,slhacot,slhasit,slhacob,slhasib,
     .                 slhacol,slhasil
      COMMON/SLHA_hmass_HDEC/slhaml,slhamh,slhamc,slha_alpha
      COMMON/GAUGINOMIX_HDEC/ZZ(4,4),UU(2,2),VV(2,2)
      COMMON/TAUMIX_HDEC/CL,SL
      COMMON/SLHAVAL_HDEC/g1ew,g2ew
      COMMON/SLHA_checkval_HDEC/check

      PI = 4*DATAN(1D0)

      if(islhao.eq.1) then
c change susyhit
         open(nout,file='hdecay_slha.out')
c end change susyhit

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

c ----------------------------------- c
c Information about the decay program c
c ----------------------------------- c

      write(nout,105)
      write(nout,51) 'DCINFO','Decay Program information'
      write(nout,61) 1,'HDECAY      # decay calculator'
c maggie changed 28/12/08
      write(nout,61) 2,'3.400       # version number'
c end maggie changed 28/12/08

c ----------------------------------------------------------------- c
c The program information: Which spectrum calculator has been used. c
c ----------------------------------------------------------------- c

      if(check(22).eq.1) then
         write(nout,105)
         write(nout,51) 'SPINFO','Spectrum calculator information'
         write(nout,61) 1,spinfo1(1:50)
         write(nout,61) 2,spinfo2(1:50)
      endif

c ------------------------------------------------ c
c Information on the model which has been selected c
c ------------------------------------------------ c

      write(nout,105)
      write(nout,51) 'MODSEL','Model selection'
c change susyhit
      write(nout,62) imod(1),imod(2),modselval(1:50)
c end change susyhit

c ----------------------- c
c The SM input parameters c
c ----------------------- c

c     if(smval(1).ne.0.d0)then
c      salpha_MS = 1/smval(1)
c     else
c      salpha_MS = 1/127.934D0
c     endif
c -- calculation of mb(mb)_MSbar from the mb pole mass --
      del = 1.d-8
      rmb0 = amb
444   rmb = rmb0
      rmb0 = runm_hdec(rmb,5)
      if(dabs(rmb0/rmb-1).gt.del)goto 444
      rmb = rmb0
      alsmz = alphas_hdec(amz,2)
      write(nout,105)
      write(nout,51) 'SMINPUTS','Standard Model inputs'
c     write(nout,52) 1,1.D0/salpha_MS,'alpha_em^-1(M_Z)^MSbar'
      write(nout,52) 2,gf,'G_F [GeV^-2]'
      write(nout,52) 3,alsmz,'alpha_S(M_Z)^MSbar'
      write(nout,52) 4,amz,'M_Z pole mass'
      write(nout,52) 5,rmb,'mb(mb)^MSbar'
      write(nout,52) 6,amt,'mt pole mass'
      write(nout,52) 7,amtau,'mtau pole mass'

c ------------------------------------------------ c
c Input parameters for minimal/default SUSY models c
c ------------------------------------------------ c

      if(check(3).eq.1) then
         write(nout,105)
         write(nout,51) 'MINPAR','Input parameters - minimal models'
         unlikely = -123456789D0
         do ii=1,20,1
            if(minval(ii).ne.unlikely) then
               write(nout,52) ii,minval(ii),mincom(ii)
            endif
         end do
      endif

c ------------------------------------------------------------------- c
c Optional input parameters for non-minimal/non-universal SUSY models c
c ------------------------------------------------------------------- c

      if(check(4).eq.1) then
         write(nout,105)
         write(nout,51) 'EXTPAR','Input parameters - non-minimal models'
         unlikely = -123456789D0
         do ii=1,100,1
            if(extval(ii-1).ne.unlikely) then
               write(nout,52) ii-1,extval(ii-1),extcom(ii-1)
            endif
         end do
      endif

c ----------------- c
c The mass spectrum c
c ----------------- c

      write(nout,105)
      write(nout,51) 'MASS','Mass Spectrum'
      write(nout,50) 'PDG code           mass       particle'
      write(nout,52) iwc,amw,'W+'
      write(nout,52) ihl,aml,'h'
      write(nout,52) ihh,amh,'H'
      write(nout,52) iha,ama,'A'
      write(nout,52) ihc,amch,'H+'
      write(nout,52) ib,amb,
     .'b-quark pole mass calculated from mb(mb)_Msbar'
      write(nout,52) isdl,gmsd(1),'~d_L'
      write(nout,52) isdr,gmsd(2),'~d_R'
      write(nout,52) isul,gmsu(1),'~u_L'
      write(nout,52) isur,gmsu(2),'~u_R'
      write(nout,52) issl,gmsd(1),'~s_L'
      write(nout,52) issr,gmsd(2),'~s_R'
      write(nout,52) iscl,gmsu(1),'~c_L'
      write(nout,52) iscr,gmsu(2),'~c_R'
      write(nout,52) isb1,GMSB(1),'~b_1'
      write(nout,52) isb2,GMSB(2),'~b_2'
      write(nout,52) ist1,GMST(1),'~t_1'
      write(nout,52) ist2,GMST(2),'~t_2'
      write(nout,52) isell,GMSE(1),'~e_L'
      write(nout,52) iselr,GMSE(2),'~e_R'
      write(nout,52) inel,GMSN(1),'~nu_eL'
      write(nout,52) ismul,GMSE(1),'~mu_L'
      write(nout,52) ismur,GMSE(2),'~mu_R'
      write(nout,52) inmul,GMSN(1),'~nu_muL'
      write(nout,52) istau1,GMSL(1),'~tau_1'
      write(nout,52) istau2,GMSL(2),'~tau_2'
      write(nout,52) intau1,GMSN(1),'~nu_tauL'
      write(nout,52) iglo,amgluino,'~g'
      write(nout,52) in1,xmn(1),'~chi_10'
      write(nout,52) in2,xmn(2),'~chi_20'
      write(nout,52) in3,xmn(3),'~chi_30'
      write(nout,52) in4,xmn(4),'~chi_40'
      write(nout,52) ic1,gmc(1),'~chi_1+'
      write(nout,52) ic2,gmc(2),'~chi_2+'

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
      write(nout,53) 1,1,ctht,'cos(theta_t)'
      write(nout,53) 1,2,stht,'sin(theta_t)'
      write(nout,53) 2,1,-stht,'-sin(theta_t)'
      write(nout,53) 2,2,ctht,'cos(theta_t)'

      write(nout,105)
      write(nout,51) 'SBOTMIX','Sbottom Mixing Matrix'
      write(nout,53) 1,1,cthb,'cos(theta_b)'
      write(nout,53) 1,2,sthb,'sin(theta_b)'
      write(nout,53) 2,1,-sthb,'-sin(theta_b)'
      write(nout,53) 2,2,cthb,'cos(theta_b)'

      write(nout,105)
      write(nout,51) 'STAUMIX','Stau Mixing Matrix'
      write(nout,53) 1,1,cl,'cos(theta_tau)'
      write(nout,53) 1,2,sl,'sin(theta_tau)'
      write(nout,53) 2,1,-sl,'-sin(theta_tau)'
      write(nout,53) 2,2,cl,'cos(theta_tau)'

c ------------------------------------------------------------------- c
c The angle alpha in the Higgs sector and the Higgs mixing parameters c
c ------------------------------------------------------------------- c

      alphaval = A
      write(nout,105)
      write(nout,51) 'ALPHA','Higgs mixing'
      write(nout,60) alphaval,
     .'Mixing angle in the neutral Higgs boson sector'

      amudrbar = AMU/(1.D0+g1ew**2/16.D0/pi**2*3.D0/5.D0+
     .                 g2ew**2/16.D0/pi**2*3.D0/4.D0)
      if(qvalue(1).ne.0.d0)then
       qq = qvalue(1)
      else
       qq = amt
      endif
      write(nout,105)
      write(nout,54) 'HMIX Q=',qq,'DRbar Higgs Parameters'
      write(nout,52) 1,amudrbar,'mu(Q)'
      write(nout,52) 2,tgbet,'tanbeta(Q)'

c ------------------- c
c The gauge couplings c
c ------------------- c
 
c     del = 1.d-8
c     g2ew0 = g2ew
c80   g2test  = g2ew/(1-g2ew0**2/96/pi**2*2)
c     g2ew1  = g2test*(1-g2test**2/96/pi**2*2)
c     write(6,*)g2ew,g2ew1,g2test
c     g2ew0 = g2test
c     if(dabs(g2ew1/g2ew-1).gt.del)goto 80
c     g2drbar = g2test

      if(qvalue(2).ne.0.d0)then
       write(nout,105)
       write(nout,54) 'GAUGE Q=',qvalue(2),'The gauge couplings'
       if(gaugeval(1).ne.0.D0) then
          write(nout,55) 1,gaugeval(1),'gprime(Q) DRbar'
       endif
       if(gaugeval(2).ne.0.D0) then
          write(nout,55) 2,gaugeval(2),'g(Q) DRbar'
       endif
      endif

c ------------------------------------- c
c The trilinear couplings Au, Ad and Ae c
c ------------------------------------- c

      qq = amt
      if(qvalue(4).ne.0.d0)qq = qvalue(4)
      write(nout,105)
      write(nout,54) 'AU Q=',qq,'The trilinear couplings'
c change susyhit
      write(nout,53) 1,1,auval(1,1),'A_u(Q) DRbar'
      write(nout,53) 2,2,auval(2,2),'A_c(Q) DRbar'
      write(nout,53) 3,3,auval(3,3),'A_t(Q) DRbar'
c end change susyhit

      qq = amt
      if(qvalue(5).ne.0.d0)qq = qvalue(5)
      write(nout,105)
      write(nout,54) 'AD Q=',qq,'The trilinear couplings'
c change susyhit
      write(nout,53) 1,1,adval(1,1),'A_d(Q) DRbar'
      write(nout,53) 2,2,adval(2,2),'A_s(Q) DRbar'
      write(nout,53) 3,3,adval(3,3),'A_b(Q) DRbar'
c end change susyhit

      qq = amt
      if(qvalue(6).ne.0.d0)qq = qvalue(6)
      write(nout,105)
      write(nout,54) 'AE Q=',qq,'The trilinear couplings'
c change susyhit
      write(nout,53) 1,1,aeval(1,1),'A_e(Q) DRbar'
      write(nout,53) 2,2,aeval(2,2),'A_mu(Q) DRbar'
      write(nout,53) 3,3,aeval(3,3),'A_tau(Q) DRbar'
c end change susyhit

c ----------------------------- c
c The soft SUSY breaking masses c
c ----------------------------- c

      if(check(15).eq.1) then
         write(nout,105)
         write(nout,54) 'MSOFT Q=',scaleofewsb,'The soft SUSY breaking m
     .asses at the scale Q'
         unlikely = -123456789D0
         do ii=1,99,1
            if(msoftval(ii).ne.unlikely) then
               if(ii.ne.11.and.ii.ne.12.and.ii.ne.13.and.ii.ne.23.and.
     .            ii.ne.24.and.ii.ne.25.and.ii.ne.26) then
                  write(nout,52) ii,msoftval(ii),softcom(ii)
               endif
            endif
         end do
      else
         write(nout,105)
         write(nout,54) 'MSOFT Q=',scaleofewsb,'The soft SUSY breaking m
     .asses at the scale Q'
         cw=amw/amz
         sw=dsqrt(1-cw**2)
         tw=sw/cw
         am1=5.D0/3.D0*tw**2*am2
         am2=am2/(1.D0+g2ew**2/16.D0/pi**2*2.D0)
         write(nout,52) 1,am1,'M_1(Q)'
         write(nout,52) 2,am2,'M_2(Q)'
         write(nout,52) 31,amel1,'AMEL1'
         write(nout,52) 33,amel,'AMEL'
         write(nout,52) 34,amer1,'AMER1'
         write(nout,52) 36,amer,'AMER'
         write(nout,52) 41,amql1,'AMQL1'
         write(nout,52) 43,amsq,'AMSQ'
         write(nout,52) 44,amur1,'AMUR1'
         write(nout,52) 46,amur,'AMUR'
         write(nout,52) 47,amdr1,'AMDR1'
         write(nout,52) 49,amdr,'AMDR'
      endif

         if(ihiggs.eq.0) then
            if(smwdth.ne.0.D0) then
               write(nout,99)
               write(nout,100) 25,smwdth,'SM Higgs decays'

               write(nout,101)
      
               if(smbrb.ne.0.D0) then
      write(nout,102) smbrb,2,ib,ibb        ,'BR(H -> b       bb     )'
               endif
               if(smbrl.ne.0.D0) then
      write(nout,102) smbrl,2,-itau,itau    ,'BR(H -> tau+    tau-   )'
               endif
               if(smbrm.ne.0.D0) then
      write(nout,102) smbrm,2,-imu,imu      ,'BR(H -> mu+     mu-    )'
               endif
               if(smbrs.ne.0.D0) then
      write(nout,102) smbrs,2,is,isb        ,'BR(H -> s       sb     )'
               endif
               if(smbrc.ne.0.D0) then
      write(nout,102) smbrc,2,ic,icb        ,'BR(H -> c       cb     )'
               endif
               if(smbrt.ne.0.D0) then
      write(nout,102) smbrt,2,it,itb        ,'BR(H -> t       tb     )' 
               endif
               if(smbrg.ne.0.D0) then
      write(nout,102) smbrg,2,igl,igl       ,'BR(H -> g       g      )' 
               endif
               if(smbrga.ne.0.D0) then
      write(nout,102) smbrga,2,iga,iga      ,'BR(H -> gam     gam    )'  
               endif
               if(smbrzga.ne.0.D0) then
      write(nout,102) smbrzga,2,iga,iz      ,'BR(H -> Z       gam    )' 
               endif
               if(smbrw.ne.0.D0) then
      write(nout,102) smbrw,2,iwc,-iwc      ,'BR(H -> W+      W-     )' 
               endif
               if(smbrz.ne.0.D0) then
      write(nout,102) smbrz,2,iz,iz         ,'BR(H -> Z       Z      )' 
               endif

            elseif(smwdth.eq.0.D0) then
               write(nout,99)
               write(nout,100) 25,0.000000000E+00,'SM Higgs decays'
               
            endif
         endif

         if(ihiggs.eq.1.or.ihiggs.eq.5) then
            write(nout,105)

      if(hlwdth.ne.0.D0) then

      write(nout,99)
      write(nout,100) 26,hlwdth,'h decays'

      write(nout,101)
      if(hlbrb.ne.0.D0) then
      write(nout,102) hlbrb,2,ib,ibb        ,'BR(h -> b       bb     )'
      endif
      if(hlbrl.ne.0.D0) then
      write(nout,102) hlbrl,2,-itau,itau    ,'BR(h -> tau+    tau-   )'
      endif
      if(hlbrm.ne.0.D0) then
      write(nout,102) hlbrm,2,-imu,imu      ,'BR(h -> mu+     mu-    )'
      endif
      if(hlbrs.ne.0.D0) then
      write(nout,102) hlbrs,2,is,isb        ,'BR(h -> s       sb     )'
      endif
      if(hlbrc.ne.0.D0) then
      write(nout,102) hlbrc,2,ic,icb        ,'BR(h -> c       cb     )'
      endif
      if(hlbrt.ne.0.D0) then
      write(nout,102) hlbrt,2,it,itb        ,'BR(h -> t       tb     )' 
      endif
      if(hlbrg.ne.0.D0) then
      write(nout,102) hlbrg,2,igl,igl       ,'BR(h -> g       g      )' 
      endif
      if(hlbrga.ne.0.D0) then
      write(nout,102) hlbrga,2,iga,iga      ,'BR(h -> gam     gam    )' 
      endif
      if(hlbrzga.ne.0.D0) then
      write(nout,102) hlbrzga,2,iga,iz      ,'BR(h -> Z       gam    )' 
      endif
      if(hlbrw.ne.0.D0) then
      write(nout,102) hlbrw,2,iwc,-iwc      ,'BR(h -> W+      W-     )' 
      endif
      if(hlbrz.ne.0.D0) then
      write(nout,102) hlbrz,2,iz,iz         ,'BR(h -> Z       Z      )' 
      endif
      if(hlbrsc(1,1).ne.0.D0) then
      write(nout,102) hlbrsc(1,1),2,ic1,-ic1,'BR(h -> ~chi_1+ ~chi_1-)' 
      endif
      if(hlbrsc(2,2).ne.0.D0) then
      write(nout,102) hlbrsc(2,2),2,ic2,-ic2,'BR(h -> ~chi_2+ ~chi_2-)' 
      endif
      if(hlbrsc(1,2).ne.0.D0) then
      write(nout,102) hlbrsc(1,2),2,ic1,-ic2,'BR(h -> ~chi_1+ ~chi_2-)' 
      endif
      if(hlbrsc(2,1).ne.0.D0) then
      write(nout,102) hlbrsc(2,1),2,ic2,-ic1,'BR(h -> ~chi_2+ ~chi_1-)' 
      endif
      hlbrsn1(1,2) = 2.D0*hlbrsn(1,2) 
      hlbrsn1(1,3) = 2.D0*hlbrsn(1,3) 
      hlbrsn1(1,4) = 2.D0*hlbrsn(1,4)
      hlbrsn1(2,3) = 2.D0*hlbrsn(2,3) 
      hlbrsn1(2,4) = 2.D0*hlbrsn(2,4)  
      hlbrsn1(3,4) = 2.D0*hlbrsn(3,4)  
      if(hlbrsn(1,1).ne.0.D0) then
      write(nout,102) hlbrsn(1,1),2,in1,in1 ,'BR(h -> ~chi_10 ~chi_10)' 
      endif
      if(hlbrsn(2,2).ne.0.D0) then
      write(nout,102) hlbrsn(2,2),2,in2,in2 ,'BR(h -> ~chi_20 ~chi_20)' 
      endif
      if(hlbrsn(3,3).ne.0.D0) then
      write(nout,102) hlbrsn(3,3),2,in3,in3 ,'BR(h -> ~chi_30 ~chi_30)' 
      endif
      if(hlbrsn(4,4).ne.0.D0) then
      write(nout,102) hlbrsn(4,4),2,in4,in4 ,'BR(h -> ~chi_40 ~chi_40)' 
      endif
      if(hlbrsn(1,2).ne.0.D0) then
      write(nout,102) hlbrsn1(1,2),2,in1,in2,'BR(h -> ~chi_10 ~chi_20)' 
      endif
      if(hlbrsn(1,3).ne.0.D0) then
      write(nout,102) hlbrsn1(1,3),2,in1,in3,'BR(h -> ~chi_10 ~chi_30)' 
      endif
      if(hlbrsn(1,4).ne.0.D0) then
      write(nout,102) hlbrsn1(1,4),2,in1,in4,'BR(h -> ~chi_10 ~chi_40)' 
      endif
      if(hlbrsn(2,3).ne.0.D0) then
      write(nout,102) hlbrsn1(2,3),2,in2,in3,'BR(h -> ~chi_20 ~chi_30)' 
      endif
      if(hlbrsn(2,4).ne.0.D0) then
      write(nout,102) hlbrsn1(2,4),2,in2,in4,'BR(h -> ~chi_20 ~chi_40)' 
      endif
      if(hlbrsn(3,4).ne.0.D0) then
      write(nout,102) hlbrsn1(3,4),2,in3,in4,'BR(h -> ~chi_30 ~chi_40)' 
      endif
      bhlslnl1 = bhlslnl/3.D0
      bhlslel1 = bhlslel/2.D0
      bhlsler1 = bhlsler/2.D0
      bhlsqul1 = bhlsqul/2.d0
      bhlsqur1 = bhlsqur/2.d0
      bhlsqdl1 = bhlsqdl/2.d0
      bhlsqdr1 = bhlsqdr/2.d0
      if(bhlsqul1.ne.0.D0) then
      write(nout,102) bhlsqul1,2,isul,-isul  ,'BR(h -> ~u_L    ~u_L*  )'
      endif
      if(bhlsqur1.ne.0.D0) then
      write(nout,102) bhlsqur1,2,isur,-isur  ,'BR(h -> ~u_R    ~u_R*  )'
      endif
      if(bhlsqul1.ne.0.D0) then
      write(nout,102) bhlsqul1,2,iscl,-iscl  ,'BR(h -> ~c_L    ~c_L*  )'
      endif
      if(bhlsqur1.ne.0.D0) then
      write(nout,102) bhlsqur1,2,iscr,-iscr  ,'BR(h -> ~c_R    ~c_R*  )'
      endif
      if(bhlst(1,1).ne.0.D0) then
      write(nout,102) bhlst(1,1),2,ist1,-ist1,'BR(h -> ~t_1    ~t_1*  )'
      endif
      if(bhlst(2,2).ne.0.D0) then
      write(nout,102) bhlst(2,2),2,ist2,-ist2,'BR(h -> ~t_2    ~t_2*  )'
      endif
      if(bhlst(1,2).ne.0.D0) then
      write(nout,102) bhlst(1,2),2,ist1,-ist2,'BR(h -> ~t_1    ~t_2*  )'
      endif
      if(bhlst(2,1).ne.0.D0) then
      write(nout,102) bhlst(2,1),2,ist2,-ist1,'BR(h -> ~t_2    ~t_1*  )'
      endif
      if(bhlsqdl1.ne.0.D0) then
      write(nout,102) bhlsqdl1,2,isdl,-isdl  ,'BR(h -> ~d_L    ~d_L*  )'
      endif
      if(bhlsqdr1.ne.0.D0) then
      write(nout,102) bhlsqdr1,2,isdr,-isdr  ,'BR(h -> ~d_R    ~d_R*  )'
      endif
      if(bhlsqdl1.ne.0.D0) then
      write(nout,102) bhlsqdl1,2,issl,-issl  ,'BR(h -> ~s_L    ~s_L*  )'
      endif
      if(bhlsqdr1.ne.0.D0) then
      write(nout,102) bhlsqdr1,2,issr,-issr  ,'BR(h -> ~s_R    ~s_R*  )'
      endif
      if(bhlsb(1,1).ne.0.D0) then
      write(nout,102) bhlsb(1,1),2,isb1,-isb1,'BR(h -> ~b_1    ~b_1*  )'
      endif
      if(bhlsb(2,2).ne.0.D0) then
      write(nout,102) bhlsb(2,2),2,isb2,-isb2,'BR(h -> ~b_2    ~b_2*  )'
      endif
      if(bhlsb(1,2).ne.0.D0) then
      write(nout,102) bhlsb(1,2),2,isb1,-isb2,'BR(h -> ~b_1    ~b_2*  )'
      endif
      if(bhlsb(2,1).ne.0.D0) then
      write(nout,102) bhlsb(2,1),2,isb2,-isb1,'BR(h -> ~b_2    ~b_1*  )'
      endif
      if(bhlslel1.ne.0.D0) then
      write(nout,102) bhlslel1,2,isell,-isell,'BR(h -> ~e_L-   ~e_L+  )'
      endif
      if(bhlsler1.ne.0.D0) then
      write(nout,102) bhlsler1,2,iselr,-iselr,'BR(h -> ~e_R-   ~e_R+  )'
      endif
      if(bhlslel1.ne.0.D0) then
      write(nout,102) bhlslel1,2,ismul,-ismul,'BR(h -> ~mu_L-  ~mu_L+ )'
      endif
      if(bhlsler1.ne.0.D0) then
      write(nout,102) bhlsler1,2,ismur,-ismur,'BR(h -> ~mu_R-  ~mu_R+ )'
      endif
      if(bhlstau(1,1).ne.0.D0) then
      write(nout,102) bhlstau(1,1),2,istau1,-istau1,'BR(h -> ~tau_1- ~ta
     .u_1+)'
      endif
      if(bhlstau(2,2).ne.0.D0) then
      write(nout,102) bhlstau(2,2),2,istau2,-istau2,'BR(h -> ~tau_2- ~ta
     .u_2+)'
      endif
      if(bhlstau(1,2).ne.0.D0) then
      write(nout,102) bhlstau(1,2),2,istau1,-istau2,'BR(h -> ~tau_1- ~ta
     .u_2+)'
      endif
      if(bhlstau(2,1).ne.0.D0) then
      write(nout,102) bhlstau(2,1),2,istau2,-istau1,'BR(h -> ~tau_2- ~ta
     .u_1+)'
      endif
      if(bhlslnl1.ne.0.D0) then
      write(nout,102) bhlslnl1,2,inel,-inel  ,'BR(h -> ~nu_eL  ~nu_eL*  
     . )'
      write(nout,102) bhlslnl1,2,inmul,-inmul,'BR(h -> ~nu_muL ~nu_muL* 
     . )'
      write(nout,102) bhlslnl1,2,intau1,-intau1,'BR(h -> ~nu_tauL ~nu_ta
     .uL*)'
      endif

      elseif(hlwdth.eq.0.D0) then
      write(nout,99)
      write(nout,100) 26,0.000000000E+00,'h decays'

      endif
      endif

         if(ihiggs.eq.2.or.ihiggs.eq.5) then
            write(nout,105)

      if(hhwdth.ne.0.D0) then
      write(nout,99)
      write(nout,100) 35,hhwdth,'H decays'

      write(nout,101)
      if(hhbrb.ne.0.D0) then
      write(nout,102) hhbrb,2,ib,ibb        ,'BR(H -> b       bb     )'
      endif
      if(hhbrl.ne.0.D0) then
      write(nout,102) hhbrl,2,-itau,itau    ,'BR(H -> tau+    tau-   )'
      endif
      if(hhbrm.ne.0.D0) then
      write(nout,102) hhbrm,2,-imu,imu      ,'BR(H -> mu+     mu-    )'
      endif
      if(hhbrs.ne.0.D0) then
      write(nout,102) hhbrs,2,is,isb        ,'BR(H -> s       sb     )'
      endif
      if(hhbrc.ne.0.D0) then
      write(nout,102) hhbrc,2,ic,icb        ,'BR(H -> c       cb     )'
      endif
      if(hhbrt.ne.0.D0) then
      write(nout,102) hhbrt,2,it,itb        ,'BR(H -> t       tb     )' 
      endif
      if(hhbrg.ne.0.D0) then
      write(nout,102) hhbrg,2,igl,igl       ,'BR(H -> g       g      )' 
      endif
      if(hhbrga.ne.0.D0) then
      write(nout,102) hhbrga,2,iga,iga      ,'BR(H -> gam     gam    )' 
      endif
      if(hhbrzga.ne.0.D0) then
      write(nout,102) hhbrzga,2,iz,iga      ,'BR(H -> Z       gam    )' 
      endif
      if(hhbrw.ne.0.D0) then
      write(nout,102) hhbrw,2,iwc,-iwc      ,'BR(H -> W+      W-     )' 
      endif
      if(hhbrz.ne.0.D0) then
      write(nout,102) hhbrz,2,iz,iz         ,'BR(H -> Z       Z      )' 
      endif
      if(hhbrh.ne.0.D0) then
      write(nout,102) hhbrh,2,ihl,ihl       ,'BR(H -> h       h      )' 
      endif
      if(hhbra.ne.0.D0) then
      write(nout,102) hhbra,2,iha,iha       ,'BR(H -> A       A      )' 
      endif
      if(hhbraz.ne.0.D0) then
      write(nout,102) hhbraz,2,iz,iha       ,'BR(H -> Z       A      )' 
      endif
      if(hhbrhw.ne.0.D0) then
      write(nout,102) hhbrhw/2.D0,2,iwc,-ihc,'BR(H -> W+      H-     )'
      write(nout,102) hhbrhw/2.D0,2,-iwc,ihc,'BR(H -> W-      H+     )'
      endif
      if(hhbrsc(1,1).ne.0.D0) then
      write(nout,102) hhbrsc(1,1),2,ic1,-ic1,'BR(H -> ~chi_1+ ~chi_1-)' 
      endif
      if(hhbrsc(2,2).ne.0.D0) then
      write(nout,102) hhbrsc(2,2),2,ic2,-ic2,'BR(H -> ~chi_2+ ~chi_2-)' 
      endif
      if(hhbrsc(1,2).ne.0.D0) then
      write(nout,102) hhbrsc(1,2),2,ic1,-ic2,'BR(H -> ~chi_1+ ~chi_2-)' 
      endif
      if(hhbrsc(2,1).ne.0.D0) then
      write(nout,102) hhbrsc(2,1),2,ic2,-ic1,'BR(H -> ~chi_2+ ~chi_1-)' 
      endif
      if(hhbrsn(1,1).ne.0.D0) then
      write(nout,102) hhbrsn(1,1),2,in1,in1 ,'BR(H -> ~chi_10 ~chi_10)' 
      endif
      if(hhbrsn(2,2).ne.0.D0) then
      write(nout,102) hhbrsn(2,2),2,in2,in2 ,'BR(H -> ~chi_20 ~chi_20)' 
      endif
      if(hhbrsn(3,3).ne.0.D0) then
      write(nout,102) hhbrsn(3,3),2,in3,in3 ,'BR(H -> ~chi_30 ~chi_30)' 
      endif
      if(hhbrsn(4,4).ne.0.D0) then
      write(nout,102) hhbrsn(4,4),2,in4,in4 ,'BR(H -> ~chi_40 ~chi_40)' 
      endif
      hhbrsn1(1,2) = 2.D0*hhbrsn(1,2) 
      hhbrsn1(1,3) = 2.D0*hhbrsn(1,3) 
      hhbrsn1(1,4) = 2.D0*hhbrsn(1,4)
      hhbrsn1(2,3) = 2.D0*hhbrsn(2,3) 
      hhbrsn1(2,4) = 2.D0*hhbrsn(2,4)  
      hhbrsn1(3,4) = 2.D0*hhbrsn(3,4)  
      if(hhbrsn1(1,2).ne.0.D0) then
      write(nout,102) hhbrsn1(1,2),2,in1,in2,'BR(H -> ~chi_10 ~chi_20)' 
      endif
      if(hhbrsn1(1,3).ne.0.D0) then
      write(nout,102) hhbrsn1(1,3),2,in1,in3,'BR(H -> ~chi_10 ~chi_30)' 
      endif
      if(hhbrsn1(1,4).ne.0.D0) then
      write(nout,102) hhbrsn1(1,4),2,in1,in4,'BR(H -> ~chi_10 ~chi_40)' 
      endif
      if(hhbrsn1(2,3).ne.0.D0) then
      write(nout,102) hhbrsn1(2,3),2,in2,in3,'BR(H -> ~chi_20 ~chi_30)' 
      endif
      if(hhbrsn1(2,4).ne.0.D0) then
      write(nout,102) hhbrsn1(2,4),2,in2,in4,'BR(H -> ~chi_20 ~chi_40)' 
      endif
      if(hhbrsn1(3,4).ne.0.D0) then
      write(nout,102) hhbrsn1(3,4),2,in3,in4,'BR(H -> ~chi_30 ~chi_40)' 
      endif
      bhhslnl1 = bhhslnl/3.D0
      bhhslel1 = bhhslel/2.D0
      bhhsler1 = bhhsler/2.D0
      bhhsqul1 = bhhsqul/2.d0
      bhhsqur1 = bhhsqur/2.d0
      bhhsqdl1 = bhhsqdl/2.d0
      bhhsqdr1 = bhhsqdr/2.d0
      if(bhhsqul1.ne.0.D0) then
      write(nout,102) bhhsqul1,2,isul,-isul  ,'BR(H -> ~u_L    ~u_L*  )'
      endif
      if(bhhsqur1.ne.0.D0) then
      write(nout,102) bhhsqur1,2,isur,-isur  ,'BR(H -> ~u_R    ~u_R*  )'
      endif
      if(bhhsqul1.ne.0.D0) then
      write(nout,102) bhhsqul1,2,iscl,-iscl  ,'BR(H -> ~c_L    ~c_L*  )'
      endif
      if(bhhsqur1.ne.0.D0) then
      write(nout,102) bhhsqur1,2,iscr,-iscr  ,'BR(H -> ~c_R    ~c_R*  )'
      endif
      if(bhhst(1,1).ne.0.D0) then
      write(nout,102) bhhst(1,1),2,ist1,-ist1,'BR(H -> ~t_1    ~t_1*  )'
      endif
      if(bhhst(2,2).ne.0.D0) then
      write(nout,102) bhhst(2,2),2,ist2,-ist2,'BR(H -> ~t_2    ~t_2*  )'
      endif
      if(bhhst(1,2).ne.0.D0) then
      write(nout,102) bhhst(1,2),2,ist1,-ist2,'BR(H -> ~t_1    ~t_2*  )'
      endif
      if(bhhst(2,1).ne.0.D0) then
      write(nout,102) bhhst(2,1),2,ist2,-ist1,'BR(H -> ~t_2    ~t_1*  )'
      endif
      if(bhhsqdl1.ne.0.D0) then
      write(nout,102) bhhsqdl1,2,isdl,-isdl  ,'BR(H -> ~d_L    ~d_L*  )'
      endif
      if(bhhsqdr1.ne.0.D0) then
      write(nout,102) bhhsqdr1,2,isdr,-isdr  ,'BR(H -> ~d_R    ~d_R*  )'
      endif
      if(bhhsqdl1.ne.0.D0) then
      write(nout,102) bhhsqdl1,2,issl,-issl  ,'BR(H -> ~s_L    ~s_L*  )'
      endif
      if(bhhsqdr1.ne.0.D0) then
      write(nout,102) bhhsqdr1,2,issr,-issr  ,'BR(H -> ~s_R    ~s_R*  )'
      endif
      if(bhhsb(1,1).ne.0.D0) then
      write(nout,102) bhhsb(1,1),2,isb1,-isb1,'BR(H -> ~b_1    ~b_1*  )'
      endif
      if(bhhsb(2,2).ne.0.D0) then
      write(nout,102) bhhsb(2,2),2,isb2,-isb2,'BR(H -> ~b_2    ~b_2*  )'
      endif
      if(bhhsb(1,2).ne.0.D0) then
      write(nout,102) bhhsb(1,2),2,isb1,-isb2,'BR(H -> ~b_1    ~b_2*  )'
      endif
      if(bhhsb(2,1).ne.0.D0) then
      write(nout,102) bhhsb(2,1),2,isb2,-isb1,'BR(H -> ~b_2    ~b_1*  )'
      endif
      if(bhhslel1.ne.0.D0) then
      write(nout,102) bhhslel1,2,isell,-isell,'BR(H -> ~e_L-   ~e_L+  )'
      endif
      if(bhhsler1.ne.0.D0) then
      write(nout,102) bhhsler1,2,iselr,-iselr,'BR(H -> ~e_R-   ~e_R+  )'
      endif
      if(bhhslel1.ne.0.D0) then
      write(nout,102) bhhslel1,2,ismul,-ismul,'BR(H -> ~mu_L-  ~mu_L+ )'
      endif
      if(bhhsler1.ne.0.D0) then
      write(nout,102) bhhsler1,2,ismur,-ismur,'BR(H -> ~mu_R-  ~mu_R+ )'
      endif
      if(bhhstau(1,1).ne.0.D0) then
      write(nout,102) bhhstau(1,1),2,istau1,-istau1,'BR(H -> ~tau_1- ~ta
     .u_1+)'
      endif
      if(bhhstau(2,2).ne.0.D0) then
      write(nout,102) bhhstau(2,2),2,istau2,-istau2,'BR(H -> ~tau_2- ~ta
     .u_2+)'
      endif
      if(bhhstau(1,2).ne.0.D0) then
      write(nout,102) bhhstau(1,2),2,istau1,-istau2,'BR(H -> ~tau_1- ~ta
     .u_2+)'
      endif
      if(bhhstau(2,1).ne.0.D0) then
      write(nout,102) bhhstau(2,1),2,istau2,-istau1,'BR(H -> ~tau_2- ~ta
     .u_1+)'
      endif
      if(bhhslnl1.ne.0.D0) then
      write(nout,102) bhhslnl1,2,inel,-inel  ,'BR(H -> ~nu_eL  ~nu_eL*  
     . )'
      write(nout,102) bhhslnl1,2,inmul,-inmul,'BR(H -> ~nu_muL ~nu_muL* 
     . )'
      write(nout,102) bhhslnl1,2,intau1,-intau1,'BR(H -> ~nu_tauL ~nu_ta
     .uL*)'
      endif

      elseif(hhwdth.eq.0.D0) then
      write(nout,99)
      write(nout,100) 35,0.000000000E+00,'H decays'

      endif
      endif

      if(ihiggs.eq.3.or.ihiggs.eq.5) then
            write(nout,105)

      if(awdth.ne.0.D0) then
      write(nout,99)
      write(nout,100) 36,awdth,'A decays'

      write(nout,101)
      if(abrb.ne.0.D0) then
      write(nout,102) abrb,2,ib,ibb         ,'BR(A -> b       bb     )'
      endif
      if(abrl.ne.0.D0) then
      write(nout,102) abrl,2,-itau,itau     ,'BR(A -> tau+    tau-   )'
      endif
      if(abrm.ne.0.D0) then
      write(nout,102) abrm,2,-imu,imu       ,'BR(A -> mu+     mu-    )'
      endif
      if(abrs.ne.0.D0) then
      write(nout,102) abrs,2,is,isb         ,'BR(A -> s       sb     )'
      endif
      if(abrc.ne.0.D0) then
      write(nout,102) abrc,2,ic,icb         ,'BR(A -> c       cb     )'
      endif
      if(abrt.ne.0.D0) then
      write(nout,102) abrt,2,it,itb         ,'BR(A -> t       tb     )' 
      endif
      if(abrg.ne.0.D0) then
      write(nout,102) abrg,2,igl,igl        ,'BR(A -> g       g      )' 
      endif
      if(abrga.ne.0.D0) then
      write(nout,102) abrga,2,iga,iga       ,'BR(A -> gam     gam    )' 
      endif
      if(abrzga.ne.0.D0) then
      write(nout,102) abrzga,2,iz,iga       ,'BR(A -> Z       gam    )' 
      endif
      if(abrz.ne.0.D0) then
      write(nout,102) abrz,2,iz,ihl         ,'BR(A -> Z       h      )' 
      endif
      if(habrsc(1,1).ne.0.D0) then
      write(nout,102) habrsc(1,1),2,ic1,-ic1,'BR(A -> ~chi_1+ ~chi_1-)' 
      endif
      if(habrsc(2,2).ne.0.D0) then
      write(nout,102) habrsc(2,2),2,ic2,-ic2,'BR(A -> ~chi_2+ ~chi_2-)' 
      endif
      if(habrsc(1,2).ne.0.D0) then
      write(nout,102) habrsc(1,2),2,ic1,-ic2,'BR(A -> ~chi_1+ ~chi_2-)' 
      endif
      if(habrsc(2,1).ne.0.D0) then
      write(nout,102) habrsc(2,1),2,ic2,-ic1,'BR(A -> ~chi_2+ ~chi_1-)' 
      endif
      habrsn1(1,2) = 2.D0*habrsn(1,2) 
      habrsn1(1,3) = 2.D0*habrsn(1,3) 
      habrsn1(1,4) = 2.D0*habrsn(1,4)
      habrsn1(2,3) = 2.D0*habrsn(2,3) 
      habrsn1(2,4) = 2.D0*habrsn(2,4)  
      habrsn1(3,4) = 2.D0*habrsn(3,4)  
      if(habrsn(1,1).ne.0.D0) then
      write(nout,102) habrsn(1,1),2,in1,in1 ,'BR(A -> ~chi_10 ~chi_10)' 
      endif
      if(habrsn(2,2).ne.0.D0) then
      write(nout,102) habrsn(2,2),2,in2,in2 ,'BR(A -> ~chi_20 ~chi_20)' 
      endif
      if(habrsn(3,3).ne.0.D0) then
      write(nout,102) habrsn(3,3),2,in3,in3 ,'BR(A -> ~chi_30 ~chi_30)' 
      endif
      if(habrsn(4,4).ne.0.D0) then
      write(nout,102) habrsn(4,4),2,in4,in4 ,'BR(A -> ~chi_40 ~chi_40)' 
      endif
      if(habrsn1(1,2).ne.0.D0) then
      write(nout,102) habrsn1(1,2),2,in1,in2,'BR(A -> ~chi_10 ~chi_20)' 
      endif
      if(habrsn1(1,3).ne.0.D0) then
      write(nout,102) habrsn1(1,3),2,in1,in3,'BR(A -> ~chi_10 ~chi_30)' 
      endif
      if(habrsn1(1,4).ne.0.D0) then
      write(nout,102) habrsn1(1,4),2,in1,in4,'BR(A -> ~chi_10 ~chi_40)' 
      endif
      if(habrsn1(2,3).ne.0.D0) then
      write(nout,102) habrsn1(2,3),2,in2,in3,'BR(A -> ~chi_20 ~chi_30)' 
      endif
      if(habrsn1(2,4).ne.0.D0) then
      write(nout,102) habrsn1(2,4),2,in2,in4,'BR(A -> ~chi_20 ~chi_40)' 
      endif
      if(habrsn1(3,4).ne.0.D0) then
      write(nout,102) habrsn1(3,4),2,in3,in4,'BR(A -> ~chi_30 ~chi_40)' 
      endif
      if(habrst.ne.0.D0) then
      write(nout,102) habrst/2.D0,2,ist1,-ist2,'BR(A -> ~t_1    ~t_2*  )
     .'
      write(nout,102) habrst/2.D0,2,-ist1,ist2,'BR(A -> ~t_1*   ~t_2   )
     .'
      endif
      if(habrsb.ne.0.D0) then
      write(nout,102) habrsb/2.D0,2,isb1,-isb2,'BR(A -> ~b_1    ~b_2*  )
     .'
      write(nout,102) habrsb/2.D0,2,-isb1,isb2,'BR(A -> ~b_1*   ~b_2   )
     .'
      endif
      if(habrsl.ne.0.D0) then
      write(nout,102) habrsl/2.D0,2,istau1,-istau2,'BR(A -> ~tau_1- ~tau
     ._2+)'
      write(nout,102) habrsl/2.D0,2,-istau1,istau2,'BR(A -> ~tau_1+ ~tau
     ._2-)'
      endif

      elseif(awdth.eq.0.D0) then
      write(nout,99)
      write(nout,100) 36,0.000000000E+00,'A decays'

      endif
      endif

      if(ihiggs.eq.4.or.ihiggs.eq.5) then
            write(nout,105)

      if(hcwdth.ne.0.D0) then
      write(nout,99)
      write(nout,100) 37,hcwdth,'H+ decays'

      write(nout,101)
      if(hcbrb.ne.0.D0) then
      write(nout,102) hcbrb,2,ic,ibb        ,'BR(H+ -> c       bb     )'
      endif
      if(hcbrl.ne.0.D0) then
      write(nout,102) hcbrl,2,-itau,intau   ,'BR(H+ -> tau+    nu_tau )'
      endif
      if(hcbrm.ne.0.D0) then
      write(nout,102) hcbrm,2,-imu,inmu     ,'BR(H+ -> mu+     nu_mu  )'
      endif
      if(hcbrbu.ne.0.D0) then
      write(nout,102) hcbrbu,2,iu,ibb       ,'BR(H+ -> u       bb     )'
      endif
      if(hcbrs.ne.0.D0) then
      write(nout,102) hcbrs,2,iu,isb        ,'BR(H+ -> u       sb     )'
      endif
      if(hcbrc.ne.0.D0) then
      write(nout,102) hcbrc,2,ic,isb        ,'BR(H+ -> c       sb     )'
      endif
      if(hcbrt.ne.0.D0) then
      write(nout,102) hcbrt,2,it,ibb        ,'BR(H+ -> t       bb     )'
      endif
      if(hcbrw.ne.0.D0) then
      write(nout,102) hcbrw,2,iwc,ihl       ,'BR(H+ -> W+      h      )'
      endif
      if(hcbra.ne.0.D0) then
      write(nout,102) hcbra,2,iwc,iha       ,'BR(H+ -> W+      A      )'
      endif
      if(hcbrsu(1,1).ne.0.D0) then
      write(nout,102) hcbrsu(1,1),2,ic1,in1 ,'BR(H+ -> ~chi_1+ ~chi_10)'
      endif
      if(hcbrsu(1,2).ne.0.D0) then
      write(nout,102) hcbrsu(1,2),2,ic1,in2 ,'BR(H+ -> ~chi_1+ ~chi_20)'
      endif
      if(hcbrsu(1,3).ne.0.D0) then
      write(nout,102) hcbrsu(1,3),2,ic1,in3 ,'BR(H+ -> ~chi_1+ ~chi_30)'
      endif
      if(hcbrsu(1,4).ne.0.D0) then
      write(nout,102) hcbrsu(1,4),2,ic1,in4 ,'BR(H+ -> ~chi_1+ ~chi_40)'
      endif
      if(hcbrsu(2,1).ne.0.D0) then
      write(nout,102) hcbrsu(2,1),2,ic2,in1 ,'BR(H+ -> ~chi_2+ ~chi_10)'
      endif
      if(hcbrsu(2,2).ne.0.D0) then
      write(nout,102) hcbrsu(2,2),2,ic2,in2 ,'BR(H+ -> ~chi_2+ ~chi_20)'
      endif
      if(hcbrsu(2,3).ne.0.D0) then
      write(nout,102) hcbrsu(2,3),2,ic2,in3 ,'BR(H+ -> ~chi_2+ ~chi_30)'
      endif
      if(hcbrsu(2,4).ne.0.D0) then
      write(nout,102) hcbrsu(2,4),2,ic2,in4 ,'BR(H+ -> ~chi_2+ ~chi_40)'
      endif
      bhcsl02=bhcsl00/2.D0
      if(bhcsl02.ne.0.D0) then
      write(nout,102) bhcsl02,2,-isell,inel ,'BR(H+ -> ~e_L+   ~nu_eL )'
      write(nout,102) bhcsl02,2,-ismul,inmul,'BR(H+ -> ~mu_L+  ~nu_muL)'
      endif
      if(bhcsl11.ne.0.D0) then
      write(nout,102) bhcsl11,2,-istau1,intau1,'BR(H+ -> ~tau_1+ ~nu_tau
     .L)'
      endif
      if(bhcsl21.ne.0.D0) then
      write(nout,102) bhcsl21,2,-istau2,intau1,'BR(H+ -> ~tau_2+ ~nu_tau
     .L)'
      endif
      hcbrsq1=hcbrsq/2.D0
      if(hcbrsq1.ne.0.D0) then
      write(nout,102) hcbrsq1,2,isul,-isdl  ,'BR(H+ -> ~u_L    ~d_L*  )'
      write(nout,102) hcbrsq1,2,iscl,-issl  ,'BR(H+ -> ~c_L    ~s_L*  )'
      endif
      if(hcbrstb(1,1).ne.0.D0) then
      write(nout,102) hcbrstb(1,1),2,ist1,-isb1,'BR(H+ -> ~t_1    ~b_1* 
     . )'
      endif
      if(hcbrstb(2,2).ne.0.D0) then
      write(nout,102) hcbrstb(2,2),2,ist2,-isb2,'BR(H+ -> ~t_2    ~b_2* 
     . )'
      endif
      if(hcbrstb(1,2).ne.0.D0) then
      write(nout,102) hcbrstb(1,2),2,ist1,-isb2,'BR(H+ -> ~t_1    ~b_2* 
     . )'
      endif
      if(hcbrstb(2,1).ne.0.D0) then
      write(nout,102) hcbrstb(2,1),2,ist2,-isb1,'BR(H+ -> ~t_2    ~b_1* 
     . )'
      endif

      elseif(hcwdth.eq.0.D0) then
      write(nout,99)
      write(nout,100) 37,0.000000000E+00,'H+ decays'

      endif
      endif

 49   format('#',1x,A,E16.8)
 50   format('#',1x,A)
 51   format('BLOCK',1x,A,2x,'#',1x,A)
 551  format(1x,A,2x,'#',1x,A)
 52   format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 552  format(2x,E16.8,0P,3x,A)
 53   format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
 54   format('BLOCK',1x,A,1P,E16.8,2x,'#',1x,A)
 554  format(2x,A,1P,E16.8,2x,1x,A)
 55   format(1x,I5,3x,1P,E16.8,0P,3x,'#',1x,A)
 56   format(1x,I4,3x,'#',1x,A,E16.8)
 57   format(1x,I5,3x,1P,E16.8,0P,3x,'#',1x,A,E16.8)
 58   format(1x,I2,1x,I2,3x,'#',1x,A)
 59   format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A,E16.8)
 60   format(9x,1P,E16.8,0P,3x,'#',1x,A)
 61   format(1x,I5,3x,A)
 661  format(2x,A)
 62   format(1x,I5,1x,I5,3x,A)
 662  format(2x,A)
 63   format(1x,I5,3x,A,1x,'#',1x,A)

 99   format('#',9x,'PDG',12x,'Width')
 100  format('DECAY',1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 101  format('#',10x,'BR',9x,'NDA',6x,'ID1',7x,'ID2')
 102  format(3x,1P,E16.8,0P,3x,I2,3x,(I9,1x),(I9,1x),2x,'#',1x,A)
 103  format('#',11x,'BR',9x,'NDA',6x,'ID1',7x,'ID2',7x,'ID3')
 107  format('#',11x,'BR',9x,'NDA',6x,'ID1',7x,'ID2',7x,'ID3',7x,'ID4')
 104  format(3x,1P,E16.8,0P,3x,I2,3x,(I9,1x),(I9,1x),(I9,1x),2x,'#',
     .1x,A)
 106  format(3x,1P,E16.8,0P,3x,I2,3x,(I9,1x),(I9,1x),(I9,1x),(I9,1x),
     .2x,'#',1x,A)
 105  format('#') 

       close(nout)

      else

      IF(IHIGGS.EQ.0)THEN
      WRITE(NSA,20)AMSM,SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT
      WRITE(NSB,20)AMSM,SMBRG,SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5)THEN
      WRITE(NLA,20)AML,HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT
      WRITE(NLB,20)AML,HLBRG,HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLWDTH
      IF(IOFSUSY.EQ.0)THEN 
       WRITE(NSUSYL,22)AML,HLBRCHT,HLBRNET,HLBRSL,HLBRSQT,HLBRGD
       IF(INDIDEC.NE.0)THEN
        WRITE(NSUSYLA,23)AML,HLBRSC(1,1),HLBRSC(2,2),
     .                   HLBRSC(1,2)+HLBRSC(2,1)
        WRITE(NSUSYLB,21)AML,HLBRSN(1,1),HLBRSN(2,2),HLBRSN(3,3),
     .                   HLBRSN(4,4)
        WRITE(NSUSYLC,20)AML,HLBRSN(1,2)+HLBRSN(2,1),
     .                   HLBRSN(1,3)+HLBRSN(3,1),
     .                   HLBRSN(1,4)+HLBRSN(4,1),
     .                   HLBRSN(2,3)+HLBRSN(3,2),
     .                   HLBRSN(2,4)+HLBRSN(4,2),
     .                   HLBRSN(3,4)+HLBRSN(4,3)
        WRITE(NSUSYLD,20)AML,BHLSLNL,BHLSLEL,BHLSLER,BHLSTAU(1,1),
     .                   BHLSTAU(1,2)+BHLSTAU(2,1),BHLSTAU(2,2)
        WRITE(NSUSYLE,21)AML,BHLSQUL,BHLSQUR,BHLSQDL,BHLSQDR
      WRITE(NSUSYLF,20)AML,BHLSB(1,1),BHLSB(1,2)+BHLSB(2,1),BHLSB(2,2),
     .                   BHLST(1,1),BHLST(1,2)+BHLST(2,1),BHLST(2,2)
       ENDIF 
      ENDIF
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5)THEN
      WRITE(NHA,20)AMH,HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT
      WRITE(NHB,20)AMH,HHBRG,HHBRGA,HHBRZGA,HHBRW,HHBRZ
      WRITE(NHC,20)AMH,HHBRH,HHBRA,HHBRAZ,HHBRHW,HHWDTH
      IF(IOFSUSY.EQ.0)THEN 
       WRITE(NSUSYH,22)AMH,HHBRCHT,HHBRNET,HHBRSL,HHBRSQT,HHBRGD
       IF(INDIDEC.NE.0)THEN
        WRITE(NSUSYHA,23)AMH,HHBRSC(1,1),HHBRSC(2,2),
     .                  HHBRSC(1,2)+HHBRSC(2,1)
        WRITE(NSUSYHB,21)AMH,HHBRSN(1,1),HHBRSN(2,2),HHBRSN(3,3),
     .                   HHBRSN(4,4)
        WRITE(NSUSYHC,20)AMH,HHBRSN(1,2)+HHBRSN(2,1),
     .                   HHBRSN(1,3)+HHBRSN(3,1),
     .                   HHBRSN(1,4)+HHBRSN(4,1),
     .                   HHBRSN(2,3)+HHBRSN(3,2),
     .                   HHBRSN(2,4)+HHBRSN(4,2),
     .                   HHBRSN(3,4)+HHBRSN(4,3)
        WRITE(NSUSYHD,20)AMH,BHHSLNL,BHHSLEL,BHHSLER,BHHSTAU(1,1),
     .                   BHHSTAU(1,2)+BHHSTAU(2,1),BHHSTAU(2,2)
        WRITE(NSUSYHE,21)AMH,BHHSQUL,BHHSQUR,BHHSQDL,BHHSQDR
      WRITE(NSUSYHF,20)AMH,BHHSB(1,1),BHHSB(1,2)+BHHSB(2,1),BHHSB(2,2),
     .                   BHHST(1,1),BHHST(1,2)+BHHST(2,1),BHHST(2,2)
       ENDIF
      ENDIF
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5)THEN
      WRITE(NAA,20)AMA,ABRB,ABRL,ABRM,ABRS,ABRC,ABRT
      WRITE(NAB,22)AMA,ABRG,ABRGA,ABRZGA,ABRZ,AWDTH
      IF(IOFSUSY.EQ.0)THEN 
       WRITE(NSUSYA,22)AMA,HABRCHT,HABRNET,HABRSL,HABRST+HABRSB,HABRGD
       IF(INDIDEC.NE.0)THEN
        WRITE(NSUSYAA,23)AMA,HABRSC(1,1),HABRSC(2,2),
     .                   HABRSC(1,2)+HABRSC(2,1)
        WRITE(NSUSYAB,21)AMA,HABRSN(1,1),HABRSN(2,2),HABRSN(3,3),
     .                   HABRSN(4,4)
        WRITE(NSUSYAC,20)AMA,HABRSN(1,2)+HABRSN(2,1),
     .                   HABRSN(1,3)+HABRSN(3,1),
     .                   HABRSN(1,4)+HABRSN(4,1),
     .                   HABRSN(2,3)+HABRSN(3,2),
     .                   HABRSN(2,4)+HABRSN(4,2),
     .                   HABRSN(3,4)+HABRSN(4,3)
        WRITE(NSUSYAD,23)AMA,BHASTAU,BHASB,BHAST
       ENDIF
      ENDIF
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5)THEN
      WRITE(NCA,20)AMCH,HCBRB,HCBRL,HCBRM,HCBRS,HCBRC,HCBRT
      WRITE(NCB,22)AMCH,HCBRBU,HCBRW,HCBRA,HCWDTH
      IF(IOFSUSY.EQ.0)THEN 
       WRITE(NSUSYC,21)AMCH,HCBRCNT,HCBRSL,HCBRSQT,HCBRGD
       IF(INDIDEC.NE.0)THEN
        WRITE(NSUSYCA,21)AMCH,HCBRSU(1,1),HCBRSU(1,2),
     .                   HCBRSU(1,3),HCBRSU(1,4)
        WRITE(NSUSYCB,21)AMCH,HCBRSU(2,1),HCBRSU(2,2),
     .                   HCBRSU(2,3),HCBRSU(2,4)
        WRITE(NSUSYCC,23)AMCH,BHCSL00,BHCSL11,BHCSL21
        WRITE(NSUSYCD,22)AMCH,BHCSQ,BHCSTB(1,1),BHCSTB(1,2),
     .                   BHCSTB(2,1),BHCSTB(2,2)
       ENDIF
      ENDIF
      ENDIF

20    FORMAT(G12.6,6(1X,G10.4))
21    FORMAT(G12.6,4(1X,G10.4))
22    FORMAT(G12.6,5(1X,G10.4))
23    FORMAT(G12.6,3(1X,G10.4))
      endif

      RETURN
      END

      SUBROUTINE CLOSE_HDEC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(K=6,NI=87,NSA=85,NSB=86,NLA=88,NLB=89,NHA=90,NHB=91,
     .          NHC=92,NAA=93,NAB=94,NCA=95,NCB=96,NRA=97,NRB=98,
     .          NSUSYL=81,NSUSYA=82,NSUSYH=83,NSUSYC=84,NPAR=80,
     .          NSUSYLA=79,NSUSYLB=78,NSUSYLC=77,NSUSYLD=76,NSUSYLE=75,
     .          NSUSYLF=59,NSUSYHF=58,
     .          NSUSYHA=74,NSUSYHB=73,NSUSYHC=72,NSUSYHD=71,NSUSYHE=70,
     .          NSUSYAA=69,NSUSYAB=68,NSUSYAC=67,NSUSYAD=66,NSUSYAE=65,
     .          NSUSYCA=64,NSUSYCB=63,NSUSYCC=62,NSUSYCD=61,NSUSYCE=60)
      DIMENSION GMN(4),XMN(4),GMC(2),GMST(2),GMSB(2),GMSL(2),
     .          GMSU(2),GMSD(2),GMSE(2),GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),HHBRSN(4,4),
     .          HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4)
      DIMENSION GLTT(2,2),GLBB(2,2),GHTT(2,2),GHBB(2,2),GCTB(2,2),
     .          GLEE(2,2),GHEE(2,2),GCEN(2,2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/SFER1ST_HDEC/AMQL1,AMUR1,AMDR1,AMEL1,AMER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .               ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/GMN,XMN,GMC,GMST,GMSB,GMSL,GMSU,GMSD,GMSE,GMSN 
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD

      IF(IHIGGS.EQ.0) THEN
       CLOSE(NSA)
       CLOSE(NSB)
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5) THEN
       CLOSE(NLA)
       CLOSE(NLB) 
       CLOSE(NSUSYL)
      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5) THEN
       CLOSE(NHA)
       CLOSE(NHB) 
       CLOSE(NHC)
       CLOSE(NSUSYH)
      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5) THEN
       CLOSE(NAA)
       CLOSE(NAB) 
       CLOSE(NSUSYA)
      ENDIF

      IF(IHIGGS.EQ.4.OR.IHIGGS.EQ.5) THEN
       CLOSE(NCA)
       CLOSE(NCB) 
       CLOSE(NSUSYC)
      ENDIF

      RETURN
      END

C =====================================================================
C =========== BEGINNING OF THE SUBROUTINE FOR THE DECAYS ==============
C !!!!!!!!!!!!!! Any change below this line is at your own risk!!!!!!!!
C =====================================================================

      SUBROUTINE HDEC(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB_HDEC
      COMPLEX*16 CFACQ_HDEC,CFACSQ_HDEC
      DIMENSION XX(4),YY(4)
      DIMENSION AMCHAR(2),AMNEUT(4),XMNEUT(4),
     .          AC1(2,2),AC2(2,2),AC3(2,2),
     .          AN1(4,4),AN2(4,4),AN3(4,4),
     .          ACNL(2,4),ACNR(2,4),
     .          AMST(2),AMSB(2),AMSL(2),
     .          AMSU(2),AMSD(2),AMSE(2),AMSN(2),
     .          GLTT(2,2),GLBB(2,2),GLEE(2,2),
     .          GHTT(2,2),GHBB(2,2),GHEE(2,2),
     .          GCTB(2,2),GCEN(2,2)
      DIMENSION GMST(2),GMSB(2),GMSL(2),GMSU(2),GMSD(2),GMSE(2),
     .          GMSN(2)
      DIMENSION HLBRSC(2,2),HLBRSN(4,4),HHBRSC(2,2),
     .          HHBRSN(4,4),HABRSC(2,2),HABRSN(4,4),HCBRSU(2,4),
     .          HHBRST(2,2),HHBRSB(2,2),HCBRSTB(2,2) 
      DIMENSION WHLCH(2,2),WHLNE(4,4),WHHCH(2,2),WHHNE(4,4),
     .          WHACH(2,2),WHANE(4,4),WHCCN(2,4),
     .          WHHST(2,2),WHHSB(2,2),WHHSTAU(2,2),WHCSTB(2,2), 
     .          WHLST(2,2),WHLSB(2,2),WHLSTAU(2,2)
      DIMENSION WHLGD(4),WHCGD(2),WHHGD(4),WHAGD(4)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      DIMENSION slhaneut(4),slhaxneut(4),slhachar(2),slhau(2,2),
     .          slhav(2,2),slhaz(4,4),xmchar(2)
      COMPLEX*16 CF,CG,CI1,CI2,CA,CB,CTT,CTB,CTC,CTW,CLT,CLB,CLW,
     .           CAT,CAB,CAC,CAW,CAH,CTH,CLH,CX1,CX2,CAX1,CAX2,CTL,CAL,
     .           CSL,CSQ,CSB1,CSB2,CST1,CST2,CSL1,CSL2,
     .           CXL,CXQ,CXB1,CXB2,CXT1,CXT2,CXL1,CXL2
      COMPLEX*16 CSEL,CSER,CSUL,CSUR,CSDL,CSDR,
     .           CXEL,CXER,CXUL,CXUR,CXDL,CXDR
      COMPLEX*16 CAT0,CAB0,CAC0,CXUL0,CXUR0,CXDL0,CXDR0,CXB10,CXB20,
     .           CXT10,CXT20
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/HMASSR_HDEC/AMLR,AMHR
      COMMON/CHIMASS_HDEC/AMCHI
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/ALS_HDEC/XLAMBDA,AMC0,AMB0,AMT0,N0
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/CKMPAR_HDEC/VUS,VCB,VUB
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      COMMON/ONSHELL_HDEC/IONSH,IONWZ,IOFSUSY
      COMMON/OLDFASH_HDEC/NFGG
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/WIDTHSM_HDEC/SMBRB,SMBRL,SMBRM,SMBRS,SMBRC,SMBRT,SMBRG,
     .               SMBRGA,SMBRZGA,SMBRW,SMBRZ,SMWDTH
      COMMON/WIDTHA_HDEC/ABRB,ABRL,ABRM,ABRS,ABRC,ABRT,ABRG,ABRGA,
     .              ABRZGA,ABRZ,AWDTH
      COMMON/WIDTHHL_HDEC/HLBRB,HLBRL,HLBRM,HLBRS,HLBRC,HLBRT,HLBRG,
     .               HLBRGA,HLBRZGA,HLBRW,HLBRZ,HLBRA,HLBRAZ,HLBRHW,
     .               HLWDTH
      COMMON/WIDTHHH_HDEC/HHBRB,HHBRL,HHBRM,HHBRS,HHBRC,HHBRT,HHBRG,
     .               HHBRGA,HHBRZGA,HHBRW,HHBRZ,HHBRH,HHBRA,HHBRAZ,
     .               HHBRHW,HHWDTH
      COMMON/WIDTHHC_HDEC/HCBRB,HCBRL,HCBRM,HCBRBU,HCBRS,HCBRC,HCBRT,
     .               HCBRW,HCBRA,HCWDTH
      COMMON/WISUSY_HDEC/HLBRSC,HLBRSN,HHBRSC,HHBRSN,HABRSC,HABRSN,
     .              HCBRSU,HLBRCHT,HHBRCHT,HABRCHT,HLBRNET,HHBRNET,
     .              HABRNET,HCBRCNT,HLBRSL,HHBRSL,HCBRSL,HABRSL,HABRST,
     .              HABRSB,HHBRSQ,HHBRST,HHBRSB,HHBRSQT,HCBRSQ,HCBRSTB,
     .              HCBRSQT,HLBRSQ,HLBRSQT
      COMMON/WISFER_HDEC/BHLSLNL,BHLSLEL,BHLSLER,BHLSQUL,BHLSQUR,
     .              BHLSQDL,BHLSQDR,BHLST(2,2),BHLSB(2,2),BHLSTAU(2,2),
     .              BHHSLNL,BHHSLEL,BHHSLER,BHHSQUL,BHHSQUR,BHHSQDL,
     .              BHHSQDR,BHHST(2,2),BHHSB(2,2),BHHSTAU(2,2),
     .              BHASTAU,BHASB,BHAST,
     .              BHCSL00,BHCSL11,BHCSL21,BHCSQ,BHCSTB(2,2)
      COMMON/SMASS_HDEC/AMNEUT,XMNEUT,AMCHAR,AMST,AMSB,AMSL,
     .              AMSU,AMSD,AMSE,AMSN 
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/GOLDST_HDEC/AXMPL,AXMGD,IGOLD
      COMMON/WIGOLD_HDEC/HLBRGD,HABRGD,HHBRGD,HCBRGD
      COMMON/SLHA_gaug_HDEC/slhaneut,slhaxneut,slhachar,slhau,slhav,
     .                      slhaz,xmchar
      COMMON/SLHA_vals_HDEC/islhai,islhao
      HVV(X,Y)= GF/(4.D0*PI*DSQRT(2.D0))*X**3/2.D0*BETA_HDEC(Y)
     .            *(1.D0-4.D0*Y+12.D0*Y**2)
      AFF(X,Y)= GF/(4*PI*DSQRT(2.D0))*X**3*Y*(BETA_HDEC(Y))
      HFF(X,Y)= GF/(4*PI*DSQRT(2.D0))*X**3*Y*(BETA_HDEC(Y))**3
      CFF(Z,TB,X,Y)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB_HDEC(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2+Y/TB**2)-4.D0*X*Y)
      HV(V)=3.D0*(1.D0-8.D0*V+20.D0*V**2)/DSQRT((4.D0*V-1.D0))
     .      *DACOS((3.D0*V-1.D0)/2.D0/DSQRT(V**3))
     .      -(1.D0-V)*(47.D0/2.D0*V-13.D0/2.D0+1.D0/V)
     .      -3.D0/2.D0*(1.D0-6.D0*V+4.D0*V**2)*DLOG(V)
      HVH(X,Y)=0.25D0*( (1-X)*(-2+4*X-2*X**2+9*Y+9*X*Y-6*Y**2)
     .        /(3*Y)-2*(1-X-X**2+X**3-3*Y-2*X*Y-3*X**2*Y+3*Y**2
     .        +3*X*Y**2-Y**3)*(-PI/2- DATAN((1-2*X+X**2-Y-X*Y)/
     .         ((1-X)*DSQRT(-1.D0+2*X+2*Y-(X-Y)**2))))/DSQRT(-1.D0
     .         +2*X-(X-Y)**2+2*Y)-(1+X**2-2*Y-2*X*Y+Y**2)*DLOG(X))
      QCD0(X) = (1+X**2)*(4*SP_HDEC((1-X)/(1+X))+2*SP_HDEC((X-1)/(X+1))
     .        - 3*DLOG((1+X)/(1-X))*DLOG(2/(1+X))
     .        - 2*DLOG((1+X)/(1-X))*DLOG(X))
     .        - 3*X*DLOG(4/(1-X**2)) - 4*X*DLOG(X)
      HQCDM(X)=QCD0(X)/X+(3+34*X**2-13*X**4)/16/X**3*DLOG((1+X)/(1-X))
     .        + 3.D0/8/X**2*(7*X**2-1)
      AQCDM(X)=QCD0(X)/X + (19+2*X**2+3*X**4)/16/X*DLOG((1+X)/(1-X))
     .        + 3.D0/8*(7-X**2)
      HQCD(X)=(4.D0/3*HQCDM(BETA_HDEC(X))
     .        +2*(4.D0/3-DLOG(X))*(1-10*X)/(1-4*X))*ASH/PI
     .       + (29.14671D0 + RATCOUP*(1.570D0 - 2*DLOG(HIGTOP)/3
     .                                     + DLOG(X)**2/9))*(ASH/PI)**2
     .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
      AQCD(X)=(4.D0/3*AQCDM(BETA_HDEC(X))
     .        +2*(4.D0/3-DLOG(X))*(1-6*X)/(1-4*X))*ASH/PI
     .       + (29.14671D0 + RATCOUP*(23/6.D0 - DLOG(HIGTOP)
     .                                     + DLOG(X)**2/6))*(ASH/PI)**2
     .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
c     HQCD(X)=5.67d0*ASH/PI
c    .       + (29.14D0 + RATCOUP*(1.570D0 - 2*DLOG(HIGTOP)/3
c    .                                     + DLOG(X)**2/9))*(ASH/PI)**2
c    .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
c     AQCD(X)=5.67d0*ASH/PI
c    .       + (29.14D0 + RATCOUP*(3.83D0 - DLOG(HIGTOP)
c    .                                     + DLOG(X)**2/6))*(ASH/PI)**2
c    .       + (164.14D0 - 25.77D0*5 + 0.259D0*5**2)*(ASH/PI)**3
      QCDH(X)=1.D0+HQCD(X)
      QCDA(X)=1.D0+AQCD(X)
      TQCDH(X)=1.D0+4.D0/3*HQCDM(BETA_HDEC(X))*ASH/PI
      TQCDA(X)=1.D0+4.D0/3*AQCDM(BETA_HDEC(X))*ASH/PI
      QCDC(X,Y)=1.D0+4/3.D0*ASH/PI*(9/4.D0 + (3-2*X+2*Y)/4*DLOG(X/Y)
     .         +((1.5D0-X-Y)*LAMB_HDEC(X,Y)**2+5*X*Y)/2/LAMB_HDEC(X,Y)
     .         /(1-X-Y)*DLOG(XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .         + BIJ_HDEC(X,Y))
     .         + ASH/PI*(2*(4/3.D0-DLOG(X))
     .         - (X*2*(4/3.D0-DLOG(X)) + Y*2*(4/3.D0-DLOG(Y)))/(1-X-Y)
     .         - (X*2*(4/3.D0-DLOG(X))*(1-X+Y)
     .           +Y*2*(4/3.D0-DLOG(Y))*(1+X-Y))/LAMB_HDEC(X,Y)**2)
      QCDCI(X,Y)=1.D0+4/3.D0*ASH/PI*(3 + (Y-X)/2*DLOG(X/Y)
     .         +(2*(1-X-Y)+LAMB_HDEC(X,Y)**2)/2/LAMB_HDEC(X,Y)
     .         *DLOG(XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .         + BIJ_HDEC(X,Y))
     .         + ASH/PI*(2*(4/3.D0-DLOG(X)) + 2*(4/3.D0-DLOG(Y))
     .         - (X*2*(4/3.D0-DLOG(X))*(1-X+Y)
     .           +Y*2*(4/3.D0-DLOG(Y))*(1+X-Y))/LAMB_HDEC(X,Y)**2)
      QCDCM(X,Y)=1.D0+4/3.D0*ASH/PI*(9/4.D0 + (3-2*X+2*Y)/4*DLOG(X/Y)
     .         +((1.5D0-X-Y)*LAMB_HDEC(X,Y)**2+5*X*Y)/2/LAMB_HDEC(X,Y)
     .         /(1-X-Y)*DLOG(4*X*Y/(1-X-Y+LAMB_HDEC(X,Y))**2)
     .         + BIJ_HDEC(X,Y))
      QCDCMI(X,Y)=1.D0+4/3.D0*ASH/PI*(3 + (Y-X)/2*DLOG(X/Y)
     .         +(2*(1-X-Y)+LAMB_HDEC(X,Y)**2)/2/LAMB_HDEC(X,Y)
     .         *DLOG(4*X*Y/(1-X-Y+LAMB_HDEC(X,Y))**2)
     .         + BIJ_HDEC(X,Y))
      CQCD(Z,TB,X,Y,R)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB_HDEC(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2*R**2*QCDC(X,Y)
     .                           +Y/TB**2*QCDC(Y,X))
     .               -4.D0*X*Y*QCDCI(X,Y))
      CQCDM(Z,TB,X,Y,R)= GF/(4*PI*DSQRT(2.D0))*Z**3*LAMB_HDEC(X,Y)
     .              *((1.D0-X-Y)*(X*TB**2*R**2*QCDCM(X,Y)
     .                           +Y/TB**2*QCDCM(Y,X))
     .               -4.D0*X*Y*QCDCMI(X,Y))
      ELW(AMH,AMF,QF,ACF)=ALPH/PI*3.D0/2*QF**2
     .                              *(3.D0/2-DLOG(AMH**2/AMF**2))
     .      +GF/8/DSQRT(2.D0)/PI**2*(ACF*AMT**2
     .        +AMW**2*(3*DLOG(CS)/SS-5)+AMZ**2*(0.5D0
     .          -3*(1-4*SS*DABS(QF))**2))
      CF(CA) = -CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))**2/4
      CG(CA) = CDSQRT(1-CA)/2*CDLOG(-(1+CDSQRT(1-CA))/(1-CDSQRT(1-CA)))
      CI1(CA,CB) = CA*CB/2/(CA-CB)
     .           + CA**2*CB**2/2/(CA-CB)**2*(CF(CA)-CF(CB))
     .           + CA**2*CB/(CA-CB)**2*(CG(CA)-CG(CB))
      CI2(CA,CB) = -CA*CB/2/(CA-CB)*(CF(CA)-CF(CB))
      HGGQCD(ASG,NF)=1.D0+ASG/PI*(95.D0/4.D0-NF*7.D0/6.D0)
      SGGQCD(ASG)=ASG/PI*17.D0/6.D0
      AGGQCD(ASG,NF)=1.D0+ASG/PI*(97.D0/4.D0-NF*7.D0/6.D0)
      HFFSELF(AMH)=1.D0+GF*AMH**2/16.D0/PI**2/DSQRT(2.D0)*2.117203D0
     .            -(GF*AMH**2/16.D0/PI**2/DSQRT(2.D0))**2*32.6567D0
      HVVSELF(AMH)=1.D0+GF*AMH**2/16.D0/PI**2/DSQRT(2.D0)*2.800952D0
     .            +(GF*AMH**2/16.D0/PI**2/DSQRT(2.D0))**2*62.0308D0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     ELW(AMH,AMF,QF,ACF)=0
c     HFFSELF(AMH)=1.D0
c     HVVSELF(AMH)=1.D0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      PI=4D0*DATAN(1D0)
      SS=1.D0-(AMW/AMZ)**2
      CS=1.D0-SS

      IF(IHIGGS.NE.0)THEN
       CALL SUSYCP_HDEC(TGBET)
      ENDIF

C--DECOUPLING THE TOP QUARK FROM ALPHAS
      AMT0=3.D8

C--TOP QUARK DECAY WIDTH
      GAMT0 = GF*AMT**3/8/DSQRT(2D0)/PI*(1-AMW**2/AMT**2)**2
     .                                 *(1+2*AMW**2/AMT**2)
      IF(IHIGGS.NE.0.AND.AMT.GT.AMCH+AMB)THEN
       GAMT1 = GF*AMT**3/8/DSQRT(2D0)/PI*(1-AMCH**2/AMT**2)**2
     .        *((AMB/AMT)**2*TGBET**2 + 1/TGBET**2)
      ELSE
       GAMT1 = 0
      ENDIF
      GAMT1 = GAMT0+GAMT1

      IF(IHIGGS.EQ.0)THEN

C        =========================================================
C                              SM HIGGS DECAYS
C        =========================================================
      AMXX=AMH
      AMH=AMSM
C     =============  RUNNING MASSES 
      RMS = RUNM_HDEC(AMH,3)
      RMC = RUNM_HDEC(AMH,4)
      RMB = RUNM_HDEC(AMH,5)
      RMT = RUNM_HDEC(AMH,6)
      RATCOUP = 1
      HIGTOP = AMH**2/AMT**2

      ASH=ALPHAS_HDEC(AMH,2)
      AMC0=1.D8
      AMB0=2.D8
      AS3=ALPHAS_HDEC(AMH,2)
      AMC0=AMC
      AS4=ALPHAS_HDEC(AMH,2)
      AMB0=AMB
C     AMT0=AMT
C     =============== PARTIAL WIDTHS 
C  H ---> G G
C
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       HGG=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       DCC=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG
C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB)**2*FQCD
       DBB=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF
c      write(6,*)'XFAC(',AMH,') = ',xfac
c      write(6,*)'BR(H -> gg) = ',HGG

C  H ---> MU MU
      IF(AMH.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AMH,(AMMUON/AMH)**2)
     .    *(1+ELW(AMH,AMMUON,-1.D0,7.D0))
     .    *HFFSELF(AMH)
      ENDIF
C  H ---> TAU TAU
      IF(AMH.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AMH,(AMTAU/AMH)**2)
     .    *(1+ELW(AMH,AMTAU,-1.D0,7.D0))
     .    *HFFSELF(AMH)
      ENDIF
C  H --> SS
      IF(AMH.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS2=3.D0*HFF(AMH,(RMS/AMH)**2)
     .    *QCDH(RMS**2/AMH**2)
     .    *(1+ELW(AMH,RMS,-1.D0/3.D0,7.D0))
     .    *HFFSELF(AMH)
       IF(HS2.LT.0.D0) HS2 = 0
       HS1=3.D0*HFF(AMH,(AMS/AMH)**2)
     .    *TQCDH(AMS**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMS/AMH
       HSS = QQINT_HDEC(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      IF(AMH.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC2=3.D0*HFF(AMH,(RMC/AMH)**2)
     .    *QCDH(RMC**2/AMH**2)
     .    *(1+ELW(AMH,RMC,2.D0/3.D0,7.D0))
     .    *HFFSELF(AMH)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       HC1=3.D0*HFF(AMH,(AMC/AMH)**2)
     .    *TQCDH(AMC**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMC/AMH
       HCC = QQINT_HDEC(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      IF(AMH.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB2=3.D0*HFF(AMH,(RMB/AMH)**2)
     .    *QCDH(RMB**2/AMH**2)
     .    *(1+ELW(AMH,RMB,-1.D0/3.D0,1.D0))
     .    *HFFSELF(AMH)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       HB1=3.D0*HFF(AMH,(AMB/AMH)**2)
     .    *TQCDH(AMB**2/AMH**2)
     .    *HFFSELF(AMH)
       RAT = 2*AMB/AMH
       HBB = QQINT_HDEC(RAT,HB1,HB2)
      ENDIF

c     HB1X=3.D0*HFF(AMH,(AMB/AMH)**2)
c    .    *TQCDH(AMB**2/AMH**2)
c    .    /(BETA_HDEC(AMB**2/AMH**2))**3
c     HB2X=3.D0*HFF(AMH,(RMB/AMH)**2)
c    .    *QCDH(RMB**2/AMH**2)
c    .    /(BETA_HDEC(RMB**2/AMH**2))**3
c     RATCOUP = 0
c     deltaqcd = QCDH(RMB**2/AMH**2)
c     RATCOUP = 1
c     deltat = QCDH(RMB**2/AMH**2) - deltaqcd
c     write(6,*)'SM: MH     = ',AMH
c     write(6,*)'alphas(MZ) = ',ALPHAS_HDEC(91.d0,2)
c     write(6,*)'alphas(MH) = ',ALPHAS_HDEC(AMH,2)
c     write(6,*)'alphas(mb) = ',ALPHAS_HDEC(AMB,2)
c     write(6,*)'mb,mb(mb)  = ',AMB,RUNM_HDEC(AMB,5)
c     write(6,*)'mb(MH,100) = ',RMB,RUNM_HDEC(100.D0,5)
c     write(6,*)'deltaqcd,t = ',deltaqcd,deltat
c     write(6,*)'Gamma(0)   = ',HB2X,HB1X
c     write(6,*)'Gamma(mb)  = ',HB2,HB1
C  H ---> TT
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMH.LE.AMT+AMW+AMB) THEN
       HTT=0.D0
       ELSEIF (AMH.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMH**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS_HDEC(AMH,AMT,AMB,AMW,HTTS)
        HTT=FACTT*HTTS
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS_HDEC(XX(1),AMT,AMB,AMW,HTTS)
        YY(1)=FACTT*HTTS
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTTS_HDEC(XX(2),AMT,AMB,AMW,HTTS)
        YY(2)=FACTT*HTTS
        XMT = RUNM_HDEC(XX(3),6)
        XY2=3.D0*HFF(XX(3),(XMT/XX(3))**2)
     .    *QCDH(XMT**2/XX(3)**2)
     .    *HFFSELF(XX(3))
        IF(XY2.LT.0.D0) XY2 = 0
        XY1=3.D0*HFF(XX(3),(AMT/XX(3))**2)
     .    *TQCDH(AMT**2/XX(3)**2)
     .    *HFFSELF(XX(3))
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT_HDEC(RAT,XY1,XY2)
        XMT = RUNM_HDEC(XX(4),6)
        XY2=3.D0*HFF(XX(4),(XMT/XX(4))**2)
     .    *QCDH(XMT**2/XX(4)**2)
     .    *HFFSELF(XX(4))
        IF(XY2.LT.0.D0) XY2 = 0
        XY1=3.D0*HFF(XX(4),(AMT/XX(4))**2)
     .    *TQCDH(AMT**2/XX(4)**2)
     .    *HFFSELF(XX(4))
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT_HDEC(RAT,XY1,XY2)
        HTT = FINT_HDEC(AMH,XX,YY)
       ELSE
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)
     .    *QCDH(RMT**2/AMH**2)
     .    *HFFSELF(AMH)
        IF(HT2.LT.0.D0) HT2 = 0
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)
     .    *TQCDH(AMT**2/AMH**2)
     .    *HFFSELF(AMH)
        RAT = 2*AMT/AMH
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMH.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)
     .    *QCDH(RMT**2/AMH**2)
     .    *HFFSELF(AMH)
        IF(HT2.LT.0.D0) HT2 = 0
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)
     .    *TQCDH(AMT**2/AMH**2)
     .    *HFFSELF(AMH)
        RAT = 2*AMT/AMH
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       XRMC = RUNM_HDEC(AMH/2,4)*AMC/RUNM_HDEC(AMC,4)
       XRMB = RUNM_HDEC(AMH/2,5)*AMB/RUNM_HDEC(AMB,5)
       XRMT = RUNM_HDEC(AMH/2,6)*AMT/RUNM_HDEC(AMT,6)
       CTT = 4*XRMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*XRMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTC = 4*XRMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))
     .     * CFACQ_HDEC(0,AMH,XRMT)
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))
     .     * CFACQ_HDEC(0,AMH,XRMB)
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))
     .     * CFACQ_HDEC(0,AMH,XRMC)
       CAL =         2*CTL*(1+(1-CTL)*CF(CTL))
       XFAC = CDABS(CAT+CAB+CAC+CAL+CAW)**2
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      CAT0 = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))
c      CAB0 = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))
c      CAC0 = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))
c      XFACLO = CDABS(CAT0+CAB0+CAC0+CAL+CAW)**2
c      write(54,('4(1X,E12.6)'))AMH,HGA,HGA*XFACLO/XFAC,XFAC/XFACLO-1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       EPS=1.D-8
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))
       XFAC = CDABS(CAT+CAB+CAW)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV_HDEC(AMH,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AMH**3*HTWW
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV_HDEC(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
     .       *HVVSELF(XX(3))
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
     .       *HVVSELF(XX(4))
        HWW = FINT_HDEC(AMH,XX,YY)
       ELSE
        HWW=HVV(AMH,AMW**2/AMH**2)
     .     *HVVSELF(AMH)
       ENDIF
      ELSE
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
      IF (AMH.LE.AMW) THEN
       HWW=0
      ELSE IF (AMH.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AMH**2)*CWW*AMH
      ELSE IF (AMH.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
     .      *HVVSELF(XX(3))
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
     .      *HVVSELF(XX(4))
       HWW = FINT_HDEC(AMH,XX,YY)
      ELSE
       HWW=HVV(AMH,AMW**2/AMH**2)
     .     *HVVSELF(AMH)
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV_HDEC(AMH,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AMH**3*HTZZ
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV_HDEC(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
     .       *HVVSELF(XX(3))
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
     .       *HVVSELF(XX(4))
        HZZ = FINT_HDEC(AMH,XX,YY)
       ELSE
        HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0
     .      *HVVSELF(AMH)
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AMH.LE.AMZ) THEN
      HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH
      ELSE IF (AMH.LT.XM2) THEN
      CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
      XX(1) = XM1-1D0
      XX(2) = XM1
      XX(3) = XM2
      XX(4) = XM2+1D0
      YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
      YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
      YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
     .      *HVVSELF(XX(3))
      YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
     .      *HVVSELF(XX(4))
      HZZ = FINT_HDEC(AMH,XX,YY)
      ELSE
      HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0
     .   *HVVSELF(AMH)
      ENDIF
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
C
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ
      SMBRT=HTT/WTOT
      SMBRB=HBB/WTOT
      SMBRL=HLL/WTOT
      SMBRM=HMM/WTOT
      SMBRC=HCC/WTOT
      SMBRS=HSS/WTOT
      SMBRG=HGG/WTOT
      SMBRGA=HGA/WTOT
      SMBRZGA=HZGA/WTOT
      SMBRW=HWW/WTOT
      SMBRZ=HZZ/WTOT
      SMWDTH=WTOT

      AMH=AMXX

      endif

      IF(IHIGGS.GT.0)THEN

C +++++++++++++++++++++++  SUSY HIGGSSES +++++++++++++++++++++++
C
      CALL GAUGINO_HDEC(AMU,AM2,B,A,AMCHAR,AMNEUT,XMNEUT,AC1,AC2,AC3,
     .             AN1,AN2,AN3,ACNL,ACNR,AGDL,AGDA,AGDH,AGDC)
C
      TSC = (AMSQ+AMUR+AMDR)/3
      BSC = (AMSQ+AMUR+AMDR)/3
      CALL SFERMION_HDEC(TSC,BSC,AMSQ,AMUR,AMDR,AMEL,AMER,AL,AU,AD,AMU,
     .               AMST,AMSB,AMSL,AMSU,AMSD,AMSE,AMSN, 
     .               GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .               GAEE,GATT,GABB,GCEN,GCTB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     QSUSY = 1.D0/3
c     QSUSY = 1
c     QSUSY = 3
c     LOOP = 1
c     LOOP = 2
c     write(6,*)'Loop, Factor = ?'
c     read(5,*)LOOP,QSUSY
c     QSUSY = DMIN1(AMSB(1),AMSB(2),AMGLU)*QSUSY
c     QSUSY = (AMSB(1)+AMSB(2)+AMGLU)/3*QSUSY
c     QSUSY = +0.8204315362167340D3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      QSUSY = 1
      LOOP = 2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      FACTOR = 1
c     write(6,*)'Factor?'
c     read(5,*)FACTOR
      QSQ = FACTOR*AMST(1)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      ENDIF

      IF(IHIGGS.EQ.1.OR.IHIGGS.EQ.5)THEN
C        =========================================================
C                           LIGHT CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM_HDEC(AML,3)
      RMC = RUNM_HDEC(AML,4)
      RMB = RUNM_HDEC(AML,5)
      RMT = RUNM_HDEC(AML,6)
      RATCOUP = GLT/GLB
      HIGTOP = AML**2/AMT**2

      ASH=ALPHAS_HDEC(AML,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS_HDEC(AML,2)
      AMC0=AMC
      AS4=ALPHAS_HDEC(AML,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS 
C  H ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))*GLT
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))*GLB
       CTC = 4*AMC**2/AML**2*DCMPLX(1D0,-EPS)
       CAC = 2*CTC*(1+(1-CTC)*CF(CTC))*GLT
C
       IF(IOFSUSY.EQ.0) THEN 
       CSB1= 4*AMSB(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXB1=-AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GLBB(1,1)
       CXB2=-AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GLBB(2,2)
       CXT1=-AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GLTT(1,1)
       CXT2=-AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GLTT(2,2)

       CSUL = 4*AMSU(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AML**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AML**2*DCMPLX(1D0,-EPS)
       CXUL=2*(1.D0/2.D0-2.D0/3.D0*SS)*AMZ**2/AMSU(1)**2*DSIN(A+B)
     .      *CSUL*(1-CSUL*CF(CSUL))
       CXUR=2*(2.D0/3.D0*SS)*AMZ**2/AMSU(2)**2*DSIN(A+B)
     .      *CSUR*(1-CSUR*CF(CSUR))
       CXDL=2*(-1.D0/2.D0+1.D0/3.D0*SS)*AMZ**2/AMSD(1)**2*DSIN(A+B)
     .      *CSDL*(1-CSDL*CF(CSDL))
       CXDR=2*(-1.D0/3.D0*SS)*AMZ**2/AMSD(2)**2*DSIN(A+B)
     .      *CSDR*(1-CSDR*CF(CSDR))

       ELSE
       CXB1=0.D0 
       CXB2=0.D0 
       CXT1=0.D0 
       CXT2=0.D0 

       CXUL=0.D0 
       CXUR=0.D0 
       CXDL=0.D0 
       CXDR=0.D0 
       ENDIF

       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       HGG=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8

c      write(6,*)'glb, glt: ',glb,glt

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       DCC=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8 - HGG

C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       DBB=HVV(AML,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  H ---> MU MU
      IF(AML.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AML,(AMMUON/AML)**2)*GLB**2
      ENDIF
C  H ---> TAU TAU
      IF(AML.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AML,(AMTAU/AML)**2)*GLB**2
      ENDIF
C  H --> SS
      IF(AML.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*HFF(AML,(AMS/AML)**2)
     .    *GLB**2
     .    *TQCDH(AMS**2/AML**2)
       HS2=3.D0*HFF(AML,(RMS/AML)**2)*GLB**2
     .    *QCDH(RMS**2/AML**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AML
       HSS = QQINT_HDEC(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      RATCOUP = 1
      IF(AML.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*HFF(AML,(AMC/AML)**2)
     .    *GLT**2
     .    *TQCDH(AMC**2/AML**2)
       HC2=3.D0*HFF(AML,(RMC/AML)**2)*GLT**2
     .    *QCDH(RMC**2/AML**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AML
       HCC = QQINT_HDEC(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      QQ = AMB
      SUSY = 0
      XGLB = GLB
c     SSUSY = AML
      SSUSY = (AMSB(1)+AMSB(2)+AMGLU)/3*QSUSY
      AS0 = ALPHAS_HDEC(SSUSY,2)
      IF(IOFSUSY.EQ.0) THEN
       I0 = 1
       CALL DMBAPP_HDEC(I0,DGLB,DGHB,DGAB,QSUSY,LOOP)
       I0 = 1
       BSC = (AMSQ+AMUR+AMDR)/3
c      XMB = RUNM_HDEC(BSC,5)
       XMB = AMB
       SUSY = COFSUSY_HDEC(I0,AMB,XMB,QQ)*AS0/PI - 2*DGLB
       CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,LOOP)
      ENDIF
      RATCOUP = GLT/XGLB
      IF(AML.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*HFF(AML,(AMB/AML)**2)
     .    *(XGLB**2+XGLB*GLB*SUSY)
     .    *TQCDH(AMB**2/AML**2)
       HB2=3.D0*HFF(AML,(RMB/AML)**2)
     .    *(XGLB**2+XGLB*GLB*SUSY)
     .    *QCDH(RMB**2/AML**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AML
       HBB = QQINT_HDEC(RAT,HB1,HB2)

c      write(6,*)HB1,HB2,HBB,XGLB**2+XGLB*GLB*SUSY,XGLB**2,
c    .           XGLB*GLB*SUSY,SUSY+2*DGLB

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      XB0=3.D0*HFF(AML,(AMB/AML)**2)
c    .    *GLB**2
c      XB1=3.D0*HFF(AML,(RMB/AML)**2)
c    .    *GLB**2
c    .    *QCDH(RMB**2/AML**2)
c    .   + DBB
c      XB2=3.D0*HFF(AML,(RMB/AML)**2)
c    .    *(XGLB**2+XGLB*GLB*SUSY)
c    .    *QCDH(RMB**2/AML**2)
c    .   + DBB
c      write(51,('5(1X,G15.8)'))AMA,AML,XB0,XB1,XB2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      write(51,('4(1X,G15.8)'))AMA,AML,SUSY+2*DGLB,SUSY/(SUSY+2*DGLB)
c      write(51,('4(1X,G15.8)'))AMA,AML,HBB,2*DGLB,XGLB,SUSY-1+2*DLGB,
c    .                          DSIN(A),DCOS(A)
c      write(51,*)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c      X1 = (QCDH(RMB**2/AML**2)*HFF(AML,(RMB/AML)**2)/
c    .       HFF(AML,(AMB/AML)**2)-1)
c      X2 = (SUSY-1)

c     RATCOUP = GLT/XGLB
c      HB1X=3.D0*HFF(AML,(AMB/AML)**2)
c    .    *XGLB**2
c    .    *TQCDH(AMB**2/AML**2)
c    .    /(BETA_HDEC(AMB**2/AML**2))**3
c    .    *SUSY
c      HB2X=3.D0*HFF(AML,(RMB/AML)**2)*XGLB**2
c    .    *QCDH(RMB**2/AML**2)
c    .    /(BETA_HDEC(RMB**2/AML**2))**3
c    .    *SUSY
c      HB1X=3.D0*HFF(AML,(RMB/AML)**2)*GLB**2
c    .    *QCDH(RMB**2/AML**2)
c    .    /(BETA_HDEC(RMB**2/AML**2))**3
c    .    *(SUSY+2*DGLB)

c     RATCOUP = 0
c     deltaqcd = QCDH(RMB**2/AML**2)
c     RATCOUP = GLT/XGLB
c     deltat = QCDH(RMB**2/AML**2) - deltaqcd

c      write(6,*)
c      write(6,*)'h:'
c      write(6,*)'MB,RUNMB,alpha_s: ',AMB,RMB,ASH
c      write(6,*)'Mh =              ',AML
c      write(6,*)'MA =              ',AMA
c      write(6,*)'Delta(mb) = ',-DGAB
c      write(6,*)'QCD           SUSY          APPROX',
c    .           '        APPROX/FULL  Gbh(QCD)    Gbh(SQCD):'
c      write(6,*)X1,X2+2*DGLB,2*DGLB,2*DGLB/(X2+2*DGLB),GLB,XGLB
c      write(6,*)'Resummation: ',(XGLB/GLB)**2-1
c      write(6,*)'Rest:        ',SUSY-1
c      write(6,*)'Rest:        ',SUSY-1,dtan(a),tgbet
c      write(6,*)AMSQ,AMUR,AMDR,(SUSY-1)/(X2+2*DGLB)
c      write(6,*)'Total SUSY:  ',(XGLB/GLB)**2*SUSY-1
c      write(6,*)'deltaqcd,t = ',deltaqcd,deltat
c      write(6,*)'Gamma(0)   = ',AMA,HB2X,HB1X
c      write(6,*)'Gamma(mb)  = ',HB2,HB1
c      write(6,*)
c      write(9,*)AMA,AML,HB2X,HB2X/SUSY,GLB,XGLB
c      write(6,*)'Rest: h      ',AMA,AML,(SUSY-1)/(X2+2*DGLB)
c      write(51,*)AMA,AML,(SUSY-1)/(X2+2*DGLB)
      ENDIF
C  H ---> TT
      RATCOUP = 0
      IF (AML.LE.2*AMT) THEN
       HTT=0.D0
      ELSE
       HT1=3.D0*HFF(AML,(AMT/AML)**2)*GLT**2
     .    *TQCDH(AMT**2/AML**2)
       HT2=3.D0*HFF(AML,(RMT/AML)**2)*GLT**2
     .    *QCDH(RMT**2/AML**2)
       IF(HT2.LT.0.D0) HT2 = 0
       RAT = 2*AMT/AML
       HTT = QQINT_HDEC(RAT,HT1,HT2)
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       XRMC = RUNM_HDEC(AML/2,4)*AMC/RUNM_HDEC(AMC,4)
       XRMB = RUNM_HDEC(AML/2,5)*AMB/RUNM_HDEC(AMB,5)
       XRMT = RUNM_HDEC(AML/2,6)*AMT/RUNM_HDEC(AMT,6)
       CTT = 4*XRMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*XRMB**2/AML**2*DCMPLX(1D0,-EPS)
       CTC = 4*XRMC**2/AML**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AML**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AML**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AML**2*DCMPLX(1D0,-EPS)
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GLT
     .     * CFACQ_HDEC(0,AML,XRMT)
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GLB
     .     * CFACQ_HDEC(0,AML,XRMB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,LOOP)
c      CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*XGLB
c    .     * CFACQ_HDEC(0,AML,XRMB)
c      write(6,*)CTB,XGLB,CFACQ_HDEC(0,AML,XRMB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GLT
     .     * CFACQ_HDEC(0,AML,XRMC)
       CAL = 1.D0  * 2*CTL*(1+(1-CTL)*CF(CTL))*GLB
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))*GLVV
       CAH = -AMZ**2/2/AMCH**2*CTH*(1-CTH*CF(CTH))*GLPM
       IF(IOFSUSY.EQ.0) THEN 
        RMSU1 = RUNMS_HDEC(AML/2,AMSU(1))
        RMSU2 = RUNMS_HDEC(AML/2,AMSU(2))
        RMSD1 = RUNMS_HDEC(AML/2,AMSD(1))
        RMSD2 = RUNMS_HDEC(AML/2,AMSD(2))
        RMSB1 = RUNMS_HDEC(AML/2,AMSB(1))
        RMSB2 = RUNMS_HDEC(AML/2,AMSB(2))
        RMST1 = RUNMS_HDEC(AML/2,AMST(1))
        RMST2 = RUNMS_HDEC(AML/2,AMST(2))
        CX1 = 4*AMCHAR(1)**2/AML**2*DCMPLX(1D0,-EPS)
        CX2 = 4*AMCHAR(2)**2/AML**2*DCMPLX(1D0,-EPS)
        CSB1= 4*RMSB1**2/AML**2*DCMPLX(1D0,-EPS)
        CSB2= 4*RMSB2**2/AML**2*DCMPLX(1D0,-EPS)
        CST1= 4*RMST1**2/AML**2*DCMPLX(1D0,-EPS)
        CST2= 4*RMST2**2/AML**2*DCMPLX(1D0,-EPS)
        CSL1= 4*AMSL(1)**2/AML**2*DCMPLX(1D0,-EPS)
        CSL2= 4*AMSL(2)**2/AML**2*DCMPLX(1D0,-EPS)
        CAX1= AMW/XMCHAR(1) * 2*CX1*(1+(1-CX1)*CF(CX1))*2*AC2(1,1) 
        CAX2= AMW/XMCHAR(2) * 2*CX2*(1+(1-CX2)*CF(CX2))*2*AC2(2,2) 

        CSEL = 4*AMSE(1)**2/AML**2*DCMPLX(1D0,-EPS)
        CSER = 4*AMSE(2)**2/AML**2*DCMPLX(1D0,-EPS)
        CSUL = 4*RMSU1**2/AML**2*DCMPLX(1D0,-EPS)
        CSUR = 4*RMSU2**2/AML**2*DCMPLX(1D0,-EPS)
        CSDL = 4*RMSD1**2/AML**2*DCMPLX(1D0,-EPS)
        CSDR = 4*RMSD2**2/AML**2*DCMPLX(1D0,-EPS)
        CXEL=2*(-1/2D0+SS)*AMZ**2/AMSE(1)**2*DSIN(A+B)
     .       *CSEL*(1-CSEL*CF(CSEL))
        CXER=-2*(SS)*AMZ**2/AMSE(2)**2*DSIN(A+B)
     .       *CSER*(1-CSER*CF(CSER))
        CXUL=2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .       *AMZ**2/AMSU(1)**2*DSIN(A+B)*CSUL*(1-CSUL*CF(CSUL))
     .      * CFACSQ_HDEC(AML,RMSU1)
        CXUR=2*4.D0/3.D0*(2.D0/3.D0*SS)
     .       *AMZ**2/AMSU(2)**2*DSIN(A+B)*CSUR*(1-CSUR*CF(CSUR))
     .      * CFACSQ_HDEC(AML,RMSU2)
        CXDL=2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .       *AMZ**2/AMSD(1)**2*DSIN(A+B)*CSDL*(1-CSDL*CF(CSDL))
     .      * CFACSQ_HDEC(AML,RMSD1)
        CXDR=2/3.D0*(-1.D0/3.D0*SS)
     .       *AMZ**2/AMSD(2)**2*DSIN(A+B)*CSDR*(1-CSDR*CF(CSDR))
     .      * CFACSQ_HDEC(AML,RMSD2)

        CXB1=-1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GLBB(1,1)
     .      * CFACSQ_HDEC(AML,RMSB1)
        CXB2=-1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GLBB(2,2)
     .      * CFACSQ_HDEC(AML,RMSB2)
        CXT1=-4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GLTT(1,1)
     .      * CFACSQ_HDEC(AML,RMST1)
        CXT2=-4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GLTT(2,2)
     .      * CFACSQ_HDEC(AML,RMST2)
        CSL1= 4*AMSL(1)**2/AML**2*DCMPLX(1D0,-EPS)
        CSL2= 4*AMSL(2)**2/AML**2*DCMPLX(1D0,-EPS)
        CXL1=      -AMZ**2/AMSL(1)**2*CSL1*(1-CSL1*CF(CSL1))*GLEE(1,1)
        CXL2=      -AMZ**2/AMSL(2)**2*CSL2*(1-CSL2*CF(CSL2))*GLEE(2,2)
        XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXEL+CXER+CXUL+CXUR+CXDL+CXDR
     .      +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       ELSE 
        XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       ENDIF
       HGA=HVV(AML,0.D0)*(ALPH/PI)**2/16.D0*XFAC
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       XFACQ = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       XFACS = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXL1+CXL2)**2
       XFACSQ = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       HGA0 = HGA*XFACSQ/XFAC
       CAC0 = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GLT
       CAT0 = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GLT
       CAB0 = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GLB
       CXB10= -1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GLBB(1,1)
       CXB20= -1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GLBB(2,2)
       CXT10= -4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GLTT(1,1)
       CXT20= -4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GLTT(2,2)
       XFACLOQ = CDABS(CAT0+CAB0+CAC0+CAL+CAW+CAH)**2
       CXUL0=2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(1)**2*DSIN(A+B)*CSUL*(1-CSUL*CF(CSUL))
       CXUR0=2*4.D0/3.D0*(2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(2)**2*DSIN(A+B)*CSUR*(1-CSUR*CF(CSUR))
       CXDL0=2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(1)**2*DSIN(A+B)*CSDL*(1-CSDL*CF(CSDL))
       CXDR0=2/3.D0*(-1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(2)**2*DSIN(A+B)*CSDR*(1-CSDR*CF(CSDR))
       XFACLO = CDABS(CAT0+CAB0+CAC0+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXEL+CXER+CXUL0+CXUR0+CXDL0+CXDR0
     .      +  CXB10+CXB20+CXT10+CXT20+CXL1+CXL2)**2
       CSQ = 1+3*ALPHAS_HDEC(AML,2)
       XFACSQL = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXEL+CXER+(CXUL0+CXUR0+CXDL0+CXDR0
     .      +  CXB10+CXB20+CXT10+CXT20)*CSQ+CXL1+CXL2)**2
c      write(54,('6(1X,E12.6)'))AML,HGA0,HGA0*XFACQ/XFACSQ,
c    .       XFACSQ/XFACLO-1,XFACQ/XFACLOQ-1,(XFACSQ-XFACSQL)/XFACSQL
c      write(54,('7(1X,E12.6)'))AML,HGA,HGA*XFACQ/XFAC,HGA*XFACS/XFAC,
c    .       XFAC/XFACLO-1,XFACQ/XFACLOQ-1,(XFAC-XFACSQL)/XFACSQL
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  H ---> Z GAMMA
      IF(AML.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GLT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GLB
       EPS=1.D-8
       CTT = 4*AMT**2/AML**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AML**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AML**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AML**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GLVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMZ**2/2/AMCH**2*CI1(CTH,CLH)*GLPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AML**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AML**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AML.LE.XM1) THEN
        CALL HTOVV_HDEC(AML,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AML**3*HTWW*GLVV**2
       ELSEIF (AML.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV_HDEC(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
        HWW = FINT_HDEC(AML,XX,YY)*GLVV**2
       ELSE
        HWW=HVV(AML,AMW**2/AML**2)*GLVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMW-DLD
      XM2 = 2D0*AMW+DLU
      IF (AML.LE.AMW) THEN
       HWW=0
      ELSE IF (AML.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AML**2)*CWW*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
       HWW = FINT_HDEC(AML,XX,YY)*GLVV**2
      ELSE
       HWW=HVV(AML,AMW**2/AML**2)*GLVV**2
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AML.LE.XM1) THEN
        CALL HTOVV_HDEC(AML,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AML**3*HTZZ*GLVV**2
       ELSEIF (AML.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV_HDEC(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
        HZZ = FINT_HDEC(AML,XX,YY)*GLVV**2
       ELSE
        HZZ=HVV(AML,AMZ**2/AML**2)/2.D0*GLVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AML.LE.AMZ) THEN
       HZZ=0
      ELSE IF (AML.LE.XM1) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       HZZ=HV(AMZ**2/AML**2)*CZZ*AML*GLVV**2
      ELSE IF (AML.LT.XM2) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
       YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
       YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2D0
       YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2D0
       HZZ = FINT_HDEC(AML,XX,YY)*GLVV**2
      ELSE
       HZZ=HVV(AML,AMZ**2/AML**2)/2.D0*GLVV**2
      ENDIF
      ENDIF
C  H ---> A A
      IF (AML.LE.2.D0*AMA) THEN
      HAA=0
      ELSE
      HAA=GF/16.D0/DSQRT(2D0)/PI*AMZ**4/AML
     .   *BETA_HDEC(AMA**2/AML**2)*GLAA**2
      ENDIF
C  H ---> A Z
      IF (AML.LE.AMZ+AMA) THEN
      HAZ=0
      ELSE
      CAZ=LAMB_HDEC(AMA**2/AML**2,AMZ**2/AML**2)
     .   *LAMB_HDEC(AML**2/AMZ**2,AMA**2/AMZ**2)**2
      HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AML*CAZ*GZAL**2
      ENDIF
C  H ---> H+ W+
      IF (AML.LE.AMW+AMCH) THEN
      HHW=0
      ELSE
      CHW=LAMB_HDEC(AMCH**2/AML**2,AMW**2/AML**2)
     .   *LAMB_HDEC(AML**2/AMW**2,AMCH**2/AMW**2)**2
      HHW=GF/8.D0/DSQRT(2D0)/PI*AMZ**2*AMW**2/AML*CHW*GHVV**2
      ENDIF

C  ============================ SUSY DECAYS 
      IF(IOFSUSY.EQ.0) THEN
C
C  HL ----> CHARGINOS
C
      DO 711 I=1,2
      DO 711 J=1,2
      IF (AML.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHLCH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AML 
     .     *LAMB_HDEC(AMCHAR(I)**2/AML**2,AMCHAR(J)**2/AML**2)
     .     *( (AC2(I,J)**2+AC2(J,I)**2)*(AML**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)-4.D0*AC2(I,J)*AC2(J,I)* 
     .         XMCHAR(I)*XMCHAR(J) ) 
      ELSE
      WHLCH(I,J)=0.D0
      ENDIF
      WHLCHT=WHLCH(1,1)+WHLCH(1,2)+WHLCH(2,1)+WHLCH(2,2)
 711  CONTINUE
C
C  HL ----> NEUTRALINOS 
C
      DO 712 I=1,4
      DO 712 J=1,4
      IF (AML.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHLNE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AML 
     .         *AN2(I,J)**2*(AML**2-(XMNEUT(I)+XMNEUT(J))**2)
     .         *LAMB_HDEC(AMNEUT(I)**2/AML**2,AMNEUT(J)**2/AML**2)
      ELSE 
      WHLNE(I,J)=0.D0
      ENDIF
 712  CONTINUE
      WHLNET= WHLNE(1,1)+WHLNE(1,2)+WHLNE(1,3)+WHLNE(1,4)
     .       +WHLNE(2,1)+WHLNE(2,2)+WHLNE(2,3)+WHLNE(2,4)
     .       +WHLNE(3,1)+WHLNE(3,2)+WHLNE(3,3)+WHLNE(3,4)
     .       +WHLNE(4,1)+WHLNE(4,2)+WHLNE(4,3)+WHLNE(4,4)
CCC
C  HL ----> SLEPTONS 
C
      IF (AML.GT.2.D0*AMSE(1)) THEN
      WHLSLEL=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSE(1)**2/AML**2)*(-0.5D0+SS)**2
      ELSE
      WHLSLEL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSE(2)) THEN
      WHLSLER=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSE(2)**2/AML**2)*SS**2
      ELSE
      WHLSLER=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSN(1)) THEN
      WHLSLNL=3*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSN(1)**2/AML**2)*0.5D0**2
      ELSE
      WHLSLNL=0.D0
      ENDIF

      DO 718 I=1,2
      DO 718 J=1,2
      IF(AML.GT.AMSL(I)+AMSL(J)) THEN
      WHLSTAU(I,J)=GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLEE(I,J)**2*
     .      LAMB_HDEC(AMSL(I)**2/AML**2,AMSL(J)**2/AML**2)/AML
      ELSE
      WHLSTAU(I,J)=0.D0
      ENDIF
 718  CONTINUE

      WHLSLT=WHLSTAU(1,1)+WHLSTAU(2,1)+WHLSTAU(1,2)+WHLSTAU(2,2) 
     .       +WHLSLEL+WHLSLER+WHLSLNL
C
C  HL ----> SQUARKS 
C
      IF (AML.GT.2.D0*AMSU(1)) THEN
      WHLSQUL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSU(1)**2/AML**2)*(0.5D0-2.D0/3.D0*SS)**2
      ELSE
      WHLSQUL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSU(2)) THEN
      WHLSQUR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSU(2)**2/AML**2)*(-2.D0/3.D0*SS)**2
      ELSE
      WHLSQUR=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSD(1)) THEN
      WHLSQDL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSD(1)**2/AML**2)*(-0.5D0+1.D0/3.D0*SS)**2
      ELSE
      WHLSQDL=0.D0
      ENDIF

      IF (AML.GT.2.D0*AMSD(2)) THEN
      WHLSQDR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AML*DSIN(B+A)**2
     .      *BETA_HDEC(AMSD(2)**2/AML**2)*(+1.D0/3.D0*SS)**2
      ELSE
      WHLSQDR=0.D0
      ENDIF

      WHLSQ=WHLSQUL+WHLSQUR+WHLSQDL+WHLSQDR
      
C
C  HL ----> STOPS 
      SUSY = 1
      DO 713 I=1,2
      DO 713 J=1,2
      IF(AML.GT.AMST(I)+AMST(J)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(1,1,I,J,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHLST(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLTT(I,J)**2*
     .      LAMB_HDEC(AMST(I)**2/AML**2,AMST(J)**2/AML**2)/AML
     .          *SUSY
c      write(6,*)'h -> stop: ',I,J,AML,AMST(I),AMST(J),SUSY-1,
c    .           WHLST(I,J)/SUSY,WHLST(I,J)
c      write(6,*)'h -> stop: ',I,J,AML,AMST(I),AMST(J),SUSY-1
      ELSE
      WHLST(I,J)=0.D0
      ENDIF
 713  CONTINUE
C
C  HL ----> SBOTTOMS 
      SUSY = 1
      DO 714 I=1,2
      DO 714 J=1,2
      IF(AML.GT.AMSB(I)+AMSB(J)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(1,2,I,J,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHLSB(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GLBB(I,J)**2*
     .       LAMB_HDEC(AMSB(I)**2/AML**2,AMSB(J)**2/AML**2)/AML
     .      *SUSY
c      write(6,*)'h -> sbot: ',I,J,AML,AMSB(I),AMSB(J),SUSY-1,
c    .           WHLSB(I,J)/SUSY,WHLSB(I,J)
c      write(6,*)'h -> sbot: ',I,J,AML,AMSB(I),AMSB(J),SUSY-1
      ELSE
      WHLSB(I,J)=0.D0
      ENDIF
 714  CONTINUE
C
      WHLSTT=WHLST(1,1)+WHLST(1,2)+WHLST(2,1)+WHLST(2,2) 
      WHLSBB=WHLSB(1,1)+WHLSB(1,2)+WHLSB(2,1)+WHLSB(2,2) 
      WHLSQT=WHLSTT+WHLSBB+WHLSQ

      ELSE 
      WHLCHT=0.D0
      WHLNET=0.D0
      WHLSLT=0.D0
      WHLSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHLCH(I,J)=0.D0
        WHLST(I,J)=0.D0
        WHLSB(I,J)=0.D0
        WHLSTAU(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHLNE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF

      IF(IGOLD.NE.0)THEN
C   HL ---> GOLDSTINOS
       DO 710 I=1,4
       IF (AML.GT.AMNEUT(I)) THEN
        WHLGD(I)=AML**5/AXMPL**2/AXMGD**2/48.D0/PI*
     .           (1.D0-AMNEUT(I)**2/AML**2)**4*AGDL(I)**2
       ELSE
        WHLGD(I)=0.D0
       ENDIF
 710   CONTINUE
       WHLGDT=WHLGD(1)+WHLGD(2)+WHLGD(3)+WHLGD(4)
      ELSE
       WHLGDT=0
      ENDIF

C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HAA+HAZ+HHW
     .    +WHLCHT+WHLNET+WHLSLT+WHLSQT + WHLGDT
      HLBRT=HTT/WTOT
      HLBRB=HBB/WTOT
      HLBRL=HLL/WTOT
      HLBRM=HMM/WTOT
      HLBRS=HSS/WTOT
      HLBRC=HCC/WTOT
      HLBRG=HGG/WTOT
      HLBRGA=HGA/WTOT
      HLBRZGA=HZGA/WTOT
      HLBRW=HWW/WTOT
      HLBRZ=HZZ/WTOT
      HLBRA=HAA/WTOT
      HLBRAZ=HAZ/WTOT
      HLBRHW=HHW/WTOT
      DO 811 I=1,2
      DO 811 J=1,2
      HLBRSC(I,J)=WHLCH(I,J)/WTOT
811   CONTINUE
      DO 812 I=1,4
      DO 812 J=1,4
      HLBRSN(I,J)=WHLNE(I,J)/WTOT
812   CONTINUE
      HLBRCHT=WHLCHT/WTOT 
      HLBRNET=WHLNET/WTOT 
      HLBRSL=WHLSLT/WTOT 
      HLBRSQ=WHLSQ/WTOT 
      HLBRSQT=WHLSQT/WTOT 
      HLBRGD =WHLGDT/WTOT
      HLWDTH=WTOT

      BHLSLNL = WHLSLNL/WTOT
      BHLSLEL = WHLSLEL/WTOT
      BHLSLER = WHLSLER/WTOT
      BHLSQUL = WHLSQUL/WTOT
      BHLSQUR = WHLSQUR/WTOT
      BHLSQDL = WHLSQDL/WTOT
      BHLSQDR = WHLSQDR/WTOT
      DO I = 1,2
       DO J = 1,2
        BHLST(I,J) = WHLST(I,J)/WTOT
        BHLSB(I,J) = WHLSB(I,J)/WTOT
        BHLSTAU(I,J) = WHLSTAU( I,J)/WTOT
       ENDDO
      ENDDO

      ENDIF

      IF(IHIGGS.GT.1)THEN
      

C        =========================================================
C                       CHARGED HIGGS DECAYS
C        =========================================================
      TB=TGBET
C     =============  RUNNING MASSES 
      RMS = RUNM_HDEC(AMCH,3)
      RMC = RUNM_HDEC(AMCH,4)
      RMB = RUNM_HDEC(AMCH,5)
      RMT = RUNM_HDEC(AMCH,6)
      ASH=ALPHAS_HDEC(AMCH,2)
C     =============== PARTIAL WIDTHS 
C  H+ ---> MU NMU
      IF(AMCH.LE.AMMUON) THEN
       HMN = 0
      ELSE
      HMN=CFF(AMCH,TB,(AMMUON/AMCH)**2,0.D0)
      ENDIF
C  H+ ---> TAU NTAU
      IF(AMCH.LE.AMTAU) THEN
       HLN = 0
      ELSE
      HLN=CFF(AMCH,TB,(AMTAU/AMCH)**2,0.D0)
      ENDIF
C  H+ --> SU
      EPS = 1.D-12
      RATX = 1
      IF(AMCH.LE.AMS+EPS) THEN
       HSU = 0
      ELSE
       HSU1=3.D0*VUS**2*CQCDM(AMCH,TB,(AMS/AMCH)**2,EPS,RATX)
       HSU2=3.D0*VUS**2*CQCD(AMCH,TB,(RMS/AMCH)**2,EPS,RATX)
       IF(HSU2.LT.0.D0) HSU2 = 0
       RAT = AMS/AMCH
       HSU = QQINT_HDEC(RAT,HSU1,HSU2)
      ENDIF
C  H+ --> CS
      IF(AMCH.LE.AMS+AMC) THEN
       HSC = 0
      ELSE
       HSC1=3.D0*CQCDM(AMCH,TB,(AMS/AMCH)**2,(AMC/AMCH)**2,RATX)
       HSC2=3.D0*CQCD(AMCH,TB,(RMS/AMCH)**2,(RMC/AMCH)**2,RATX)
       IF(HSC2.LT.0.D0) HSC2 = 0
       RAT = (AMS+AMC)/AMCH
       HSC = QQINT_HDEC(RAT,HSC1,HSC2)
      ENDIF
C  H+ --> CB
      QQ = AMB
      SUSY = 0
      XGAB = GAB
c     SSUSY = AMCH
      SSUSY = (AMSB(1)+AMSB(2)+AMGLU)/3*QSUSY
      AS0 = ALPHAS_HDEC(SSUSY,2)
      IF(IOFSUSY.EQ.0) THEN
       I0 = 1
       CALL DMBAPP_HDEC(I0,DGLB,DGHB,DGAB,QSUSY,LOOP)
       I0 = 1
       BSC = (AMSQ+AMUR+AMDR)/3
c      XMB = RUNM_HDEC(BSC,5)
       XMB = AMB
c      SUSY = COFSUSY_HDEC(I0,AMB,XMB,QQ)*AS0/PI - 2*DGLB
       CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,LOOP)
      ENDIF
      RATX = XGAB/GAB
c     write(6,*)'ratio = ',ratx
      IF(AMCH.LE.AMB+AMC) THEN
       HBC = 0
      ELSE
       HBC1=3.D0*VCB**2*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMC/AMCH)**2,RATX)
       HBC2=3.D0*VCB**2*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMC/AMCH)**2,RATX)
       IF(HBC2.LT.0.D0) HBC2 = 0
       RAT = (AMB+AMC)/AMCH
       HBC = QQINT_HDEC(RAT,HBC1,HBC2)
      ENDIF
C  H+ --> BU
      EPS = 1.D-12
      IF(AMCH.LE.AMB+EPS) THEN
       HBU = 0
      ELSE
       HBU1=3.D0*VUB**2*CQCDM(AMCH,TB,(AMB/AMCH)**2,EPS,RATX)
       HBU2=3.D0*VUB**2*CQCD(AMCH,TB,(RMB/AMCH)**2,EPS,RATX)
       IF(HBU2.LT.0.D0) HBU2 = 0
       RAT = AMB/AMCH
       HBU = QQINT_HDEC(RAT,HBU1,HBU2)
      ENDIF
C  H+ --> TB :
      IF(IONSH.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = AMT+AMB-DLD
       XM2 = AMT+AMB+DLU
       IF (AMCH.LE.AMW+2*AMB) THEN
        HBT=0.D0
       ELSEIF (AMCH.LE.XM1) THEN
        FACTB=3.D0*GF**2*AMCH*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT_HDEC(AMCH,AMT,AMB,AMW,CTT0)
        HBT=FACTB*CTT0
       ELSEIF (AMCH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTB=3.D0*GF**2*XX(1)*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT_HDEC(XX(1),AMT,AMB,AMW,CTT0)
        YY(1)=FACTB*CTT0
        FACTB=3.D0*GF**2*XX(2)*AMT**4/32.D0/PI**3/TB**2
        CALL CTOTT_HDEC(XX(2),AMT,AMB,AMW,CTT0)
        YY(2)=FACTB*CTT0
        XMB = RUNM_HDEC(XX(3),5)
        XMT = RUNM_HDEC(XX(3),6)
        XYZ2 = 3.D0*CQCD(XX(3),TB,(XMB/XX(3))**2,(XMT/XX(3))**2,RATX)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        XYZ1 = 3.D0*CQCDM(XX(3),TB,(AMB/XX(3))**2,(AMT/XX(3))**2,RATX)
        RAT = (AMB+AMT)/XX(3)
        YY(3) = QQINT_HDEC(RAT,XYZ1,XYZ2)
        XMB = RUNM_HDEC(XX(4),5)
        XMT = RUNM_HDEC(XX(4),6)
        XYZ2 = 3.D0*CQCD(XX(4),TB,(XMB/XX(4))**2,(XMT/XX(4))**2,RATX)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        XYZ1 = 3.D0*CQCDM(XX(4),TB,(AMB/XX(4))**2,(AMT/XX(4))**2,RATX)
        RAT = (AMB+AMT)/XX(4)
        YY(4) = QQINT_HDEC(RAT,XYZ1,XYZ2)
        HBT = FINT_HDEC(AMCH,XX,YY)
       ELSE
        HBT2=3.D0*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMT/AMCH)**2,RATX)
        IF(HBT2.LT.0.D0) HBT2 = 0
        HBT1=3.D0*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMT/AMCH)**2,RATX)
        RAT = (AMB+AMT)/AMCH
        HBT = QQINT_HDEC(RAT,HBT1,HBT2)
       ENDIF
      ELSE
       IF (AMCH.LE.AMT+AMB) THEN
        HBT=0.D0
       ELSE
        HBT2=3.D0*CQCD(AMCH,TB,(RMB/AMCH)**2,(RMT/AMCH)**2,RATX)
        IF(HBT2.LT.0.D0) HBT2 = 0
        HBT1=3.D0*CQCDM(AMCH,TB,(AMB/AMCH)**2,(AMT/AMCH)**2,RATX)
        RAT = (AMB+AMT)/AMCH
        HBT = QQINT_HDEC(RAT,HBT1,HBT2)
       ENDIF
      ENDIF
C  H+ ---> W H
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = AMW+AML-DLD
       XM2 = AMW+AML+DLU
       IF (AMCH.LT.AML) THEN
        HWH=0
       ELSEIF (AMCH.LE.XM1) THEN
        IF(AMCH.LE.DABS(AMW-AML))THEN
         HWH=0
        ELSE
         HWH=9.D0*GF**2/16.D0/PI**3*AMW**4*AMCH*GHVV**2	
     .      *HVH((AML/AMCH)**2,(AMW/AMCH)**2)
        ENDIF
       ELSEIF (AMCH.LT.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1) = 9.D0*GF**2/16.D0/PI**3*AMW**4*XX(1)
     .         *HVH((AML/XX(1))**2,(AMW/XX(1))**2)
        YY(2) = 9.D0*GF**2/16.D0/PI**3*AMW**4*XX(2)
     .         *HVH((AML/XX(2))**2,(AMW/XX(2))**2)
        CWH=LAMB_HDEC(AML**2/XX(3)**2,AMW**2/XX(3)**2)
     .     *LAMB_HDEC(XX(3)**2/AMW**2,AML**2/AMW**2)**2
        YY(3)=GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(3)*CWH
        CWH=LAMB_HDEC(AML**2/XX(4)**2,AMW**2/XX(4)**2)
     .     *LAMB_HDEC(XX(4)**2/AMW**2,AML**2/AMW**2)**2
        YY(4)=GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(4)*CWH
        HWH = FINT_HDEC(AMCH,XX,YY)*GHVV**2
       ELSE
        CWH=LAMB_HDEC(AML**2/AMCH**2,AMW**2/AMCH**2)
     .     *LAMB_HDEC(AMCH**2/AMW**2,AML**2/AMW**2)**2
        HWH=GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GHVV**2*CWH
       ENDIF
      ELSE
       IF (AMCH.LT.AMW+AML) THEN
        HWH=0
       ELSE
        CWH=LAMB_HDEC(AML**2/AMCH**2,AMW**2/AMCH**2)
     .     *LAMB_HDEC(AMCH**2/AMW**2,AML**2/AMW**2)**2
        HWH=GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GHVV**2*CWH
       ENDIF
      ENDIF
C  H+ ---> W A
      IF(IONSH.EQ.0)THEN
       IF (AMCH.LT.AMA) THEN
        HWA=0
       ELSEIF (AMCH.LT.AMW+AMA) THEN
        IF(AMCH.LE.DABS(AMW-AMA))THEN
         HWA=0
        ELSE
         HWA=9.D0*GF**2/16.D0/PI**3*AMW**4*AMCH	
     .      *HVH((AMA/AMCH)**2,(AMW/AMCH)**2)
        ENDIF
       ELSE
        HWA=0.D0
       ENDIF
      ELSE
       IF (AMCH.LT.AMW+AMA) THEN
        HWA=0
       ELSE
        HWA=0.D0
       ENDIF
      ENDIF

C  ======================= SUSY DECAYS 
      IF(IOFSUSY.EQ.0) THEN
C
C  H+ ----> CHARGINOS+NEUTRALINOS
C
      DO 751 I=1,2
      DO 751 J=1,4
      IF (AMCH.GT.AMCHAR(I)+AMNEUT(J)) THEN
      WHCCN(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMCH
     .   *LAMB_HDEC(AMCHAR(I)**2/AMCH**2,AMNEUT(J)**2/AMCH**2)*(
     .   (ACNL(I,J)**2+ACNR(I,J)**2)*(AMCH**2-AMCHAR(I)**2-XMNEUT(J)
     .   **2)-4.D0*ACNL(I,J)*ACNR(I,J)*XMCHAR(I)*XMNEUT(J) )
      ELSE
      WHCCN(I,J)=0.D0
      ENDIF
 751  CONTINUE

      WHCCNT=WHCCN(1,1)+WHCCN(1,2)+WHCCN(1,3)+WHCCN(1,4)
     .      +WHCCN(2,1)+WHCCN(2,2)+WHCCN(2,3)+WHCCN(2,4)
C
C  H+ ----> SLEPTONS 
C
      IF (AMCH.GT.AMSE(1)+AMSN(1)) THEN
      WHCSL00=2*GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*DSIN(2.D0*B)**2
     .     *LAMB_HDEC(AMSE(1)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL00=0.D0
      ENDIF

      IF (AMCH.GT.AMSL(1)+AMSN(1)) THEN
      WHCSL11=GF/2.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GCEN(1,1)**2
     .     *LAMB_HDEC(AMSL(1)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL11=0.D0
      ENDIF

      IF (AMCH.GT.AMSL(2)+AMSN(1)) THEN
      WHCSL21=GF/2.D0/DSQRT(2D0)/PI*AMW**4/AMCH*GCEN(1,2)**2
     .     *LAMB_HDEC(AMSL(2)**2/AMCH**2,AMSN(1)**2/AMCH**2)
      ELSE 
      WHCSL21=0.D0
      ENDIF

      WHCSLT=WHCSL00+WHCSL11+WHCSL21

C
C  H+ ----> SQUARKS 
C
      IF (AMCH.GT.AMSU(1)+AMSD(1)) THEN
      WHCSQ=6*GF/4.D0/DSQRT(2D0)/PI*AMW**4/AMCH*DSIN(2.D0*B)**2
     .     *LAMB_HDEC(AMSU(1)**2/AMCH**2,AMSD(1)**2/AMCH**2)
      ELSE 
      WHCSQ=0.D0
      ENDIF
C
      DO 753 I=1,2
      DO 753 J=1,2
      IF(AMCH.GT.AMST(I)+AMSB(J)) THEN
      WHCSTB(I,J)=3*GF*AMW**4/2.D0/DSQRT(2.D0)/PI*GCTB(I,J)**2
     .      *LAMB_HDEC(AMST(I)**2/AMCH**2,AMSB(J)**2/AMCH**2)/AMCH
      ELSE
      WHCSTB(I,J)=0.D0
      ENDIF

 753  CONTINUE
C
      WHCSQT=WHCSQ+WHCSTB(1,1)+WHCSTB(1,2)+WHCSTB(2,1)+WHCSTB(2,2) 

      ELSE 
      WHCCNT=0.D0
      WHCSLT=0.D0
      WHCSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHCSTB(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,2
       DO J=1,4
        WHCCN(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF

      IF(IGOLD.NE.0)THEN
C   HC ---> GOLDSTINOS
       DO 750 I=1,2
       IF (AMCH.GT.AMCHAR(I)) THEN
        WHCGD(I)=AMCH**5/AXMPL**2/AXMGD**2/48.D0/PI*
     .           (1.D0-AMCHAR(I)**2/AMCH**2)**4*AGDC(I)**2
       ELSE
        WHCGD(I)=0.D0
       ENDIF
 750   CONTINUE
       WHCGDT=WHCGD(1)+WHCGD(2)
      ELSE
       WHCGDT=0
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
C
      WTOT=HLN+HMN+HSU+HBU+HSC+HBC+HBT+HWH+HWA+WHCCNT+WHCSLT+WHCSQT
     .    +WHCGDT

      HCBRL=HLN/WTOT
      HCBRM=HMN/WTOT
      HCBRS=HSU/WTOT
      HCBRBU=HBU/WTOT
      HCBRC=HSC/WTOT
      HCBRB=HBC/WTOT
      HCBRT=HBT/WTOT
      HCBRW=HWH/WTOT
      HCBRA=HWA/WTOT
      DO 851 I=1,2
      DO 851 J=1,4
      HCBRSU(I,J)=WHCCN(I,J)/WTOT
851   CONTINUE
      HCBRCNT=WHCCNT/WTOT
      HCBRSL=WHCSLT/WTOT 
      HCBRSQ=WHCSQ/WTOT 
      HCBRSQT=WHCSQT/WTOT 
      DO 853 I=1,2
      DO 853 J=1,2
      HCBRSTB(I,J)=WHCSTB(I,J)/WTOT
853   CONTINUE
      HCBRGD=WHCGDT/WTOT
      HCWDTH=WTOT

      BHCSL00 = WHCSL00/WTOT
      BHCSL11 = WHCSL11/WTOT
      BHCSL21 = WHCSL21/WTOT
      BHCSQ = WHCSQ/WTOT
      DO I = 1,2
       DO J = 1,2
        BHCSTB(I,J) = WHCSTB(I,J)/WTOT
       ENDDO
      ENDDO

      GAMC0 = WTOT

      ENDIF

      IF(IHIGGS.EQ.2.OR.IHIGGS.EQ.5)THEN
     
C        =========================================================
C                       HEAVY CP EVEN HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM_HDEC(AMH,3)
      RMC = RUNM_HDEC(AMH,4)
      RMB = RUNM_HDEC(AMH,5)
      RMT = RUNM_HDEC(AMH,6)
      RATCOUP = GHT/GHB
      HIGTOP = AMH**2/AMT**2

      ASH=ALPHAS_HDEC(AMH,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS_HDEC(AMH,2)
      AMC0=AMC
      AS4=ALPHAS_HDEC(AMH,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS 
C  H ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CAT = 2*CTT*(1+(1-CTT)*CF(CTT))*GHT
       CAB = 2*CTB*(1+(1-CTB)*CF(CTB))*GHB
       CTC = 4*AMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CAC = 2*CTC*(1+(1-CTC)*CF(CTC))*GHT
C
       IF(IOFSUSY.EQ.0) THEN 
       CSB1= 4*AMSB(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSB2= 4*AMSB(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST1= 4*AMST(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CST2= 4*AMST(2)**2/AMH**2*DCMPLX(1D0,-EPS)
C
       CXB1=-AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GHBB(1,1)
       CXB2=-AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GHBB(2,2)
       CXT1=-AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GHTT(1,1)
       CXT2=-AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GHTT(2,2)
C
       CSUL = 4*AMSU(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSUR = 4*AMSU(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDL = 4*AMSD(1)**2/AMH**2*DCMPLX(1D0,-EPS)
       CSDR = 4*AMSD(2)**2/AMH**2*DCMPLX(1D0,-EPS)
       CXUL=-2*(1.D0/2.D0-2.D0/3.D0*SS)*AMZ**2/AMSU(1)**2*DCOS(A+B)
     .      *CSUL*(1-CSUL*CF(CSUL))
       CXUR=-2*(2.D0/3.D0*SS)*AMZ**2/AMSU(2)**2*DCOS(A+B)
     .      *CSUR*(1-CSUR*CF(CSUR))
       CXDL=-2*(-1.D0/2.D0+1.D0/3.D0*SS)*AMZ**2/AMSD(1)**2*DCOS(A+B)
     .      *CSDL*(1-CSDL*CF(CSDL))
       CXDR=-2*(-1.D0/3.D0*SS)*AMZ**2/AMSD(2)**2*DCOS(A+B)
     .      *CSDR*(1-CSDR*CF(CSDR))
       ELSE
       CXB1=0.D0 
       CXB2=0.D0 
       CXT1=0.D0 
       CXT2=0.D0 
       CXUL=0.D0
       CXUR=0.D0
       CXDL=0.D0
       CXDR=0.D0
       ENDIF

       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       HGG=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8

c      write(6,*)'ghb, ght: ',ghb,ght

C  H ---> G G* ---> G CC   TO BE ADDED TO H ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       DCC=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG

C  H ---> G G* ---> G BB   TO BE ADDED TO H ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=HGGQCD(ASG,NFEXT)
       SQCD=SGGQCD(ASG)
       XFAC = CDABS(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .             +CXUL+CXUR+CXDL+CXDR)**2*FQCD
     .      + DREAL(DCONJG(CAT+CAB+CAC+CXB1+CXB2+CXT1+CXT2
     .                    +CXUL+CXUR+CXDL+CXDR)
     .             *(CXB1+CXB2+CXT1+CXT2+CXUL+CXUR+CXDL+CXDR))*SQCD
       DBB=HVV(AMH,0.D0)*(ASG/PI)**2*XFAC/8 - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

C  H ---> MU MU
      IF(AMH.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=HFF(AMH,(AMMUON/AMH)**2)*GHB**2
      ENDIF
C  H ---> LL
      IF(AMH.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=HFF(AMH,(AMTAU/AMH)**2)*GHB**2
      ENDIF
C  H --> SS
      IF(AMH.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*HFF(AMH,(AMS/AMH)**2)
     .    *GHB**2
     .    *TQCDH(AMS**2/AMH**2)
       HS2=3.D0*HFF(AMH,(RMS/AMH)**2)*GHB**2
     .    *QCDH(RMS**2/AMH**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AMH
       HSS = QQINT_HDEC(RAT,HS1,HS2)
      ENDIF
C  H --> CC
      RATCOUP = 1
      IF(AMH.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*HFF(AMH,(AMC/AMH)**2)
     .    *GHT**2
     .    *TQCDH(AMC**2/AMH**2)
       HC2=3.D0*HFF(AMH,(RMC/AMH)**2)*GHT**2
     .    *QCDH(RMC**2/AMH**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AMH
       HCC = QQINT_HDEC(RAT,HC1,HC2)
      ENDIF
C  H --> BB :
      QQ = AMB
      SUSY = 0
      XGHB = GHB
c     SSUSY = AMH
      SSUSY = (AMSB(1)+AMSB(2)+AMGLU)/3*QSUSY
      AS0 = ALPHAS_HDEC(SSUSY,2)
      IF(IOFSUSY.EQ.0) THEN
       I0 = 1
       CALL DMBAPP_HDEC(I0,DGLB,DGHB,DGAB,QSUSY,LOOP)
       I0 = 2
       BSC = (AMSQ+AMUR+AMDR)/3
c      XMB = RUNM_HDEC(BSC,5)
       XMB = AMB
       SUSY = COFSUSY_HDEC(I0,AMB,XMB,QQ)*AS0/PI - 2*DGHB
       CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,LOOP)
      ENDIF
      RATCOUP = GHT/XGHB
      IF(AMH.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*HFF(AMH,(AMB/AMH)**2)
     .    *(XGHB**2+XGHB*GHB*SUSY)
     .    *TQCDH(AMB**2/AMH**2)
       HB2=3.D0*HFF(AMH,(RMB/AMH)**2)
     .    *(XGHB**2+XGHB*GHB*SUSY)
     .    *QCDH(RMB**2/AMH**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AMH
       HBB = QQINT_HDEC(RAT,HB1,HB2)

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      XB0=3.D0*HFF(AMH,(AMB/AMH)**2)
c    .    *GHB**2
c      XB1=3.D0*HFF(AMH,(RMB/AMH)**2)
c    .    *GHB**2
c    .    *QCDH(RMB**2/AMH**2)
c    .   + DBB
c      XB2=3.D0*HFF(AMH,(RMB/AMH)**2)
c    .    *(XGHB**2+XGHB*GHB*SUSY)
c    .    *QCDH(RMB**2/AMH**2)
c    .   + DBB
c      write(52,('5(1X,G15.8)'))AMA,AMH,XB0,XB1,XB2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      write(52,('4(1X,G15.8)'))AMA,AMH,HBB,SUSY/(SUSY+2*DGHB)
c      write(52,('4(1X,G15.8)'))AMA,AMH,SUSY+2*DGHB,SUSY/(SUSY+2*DGHB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c      X1 = (QCDH(RMB**2/AMH**2)*HFF(AMH,(RMB/AMH)**2)/
c    .       HFF(AMH,(AMB/AMH)**2)-1)
c      X2 = (SUSY-1)

c     RATCOUP = GHT/XGHB
c      HB1X=3.D0*HFF(AMH,(AMB/AMH)**2)
c    .    *XGHB**2
c    .    *TQCDH(AMB**2/AMH**2)
c    .    /(BETA_HDEC(AMB**2/AMH**2))**3
c    .    *SUSY
c      HB2X=3.D0*HFF(AMH,(RMB/AMH)**2)*XGHB**2
c    .    *QCDH(RMB**2/AMH**2)
c    .    /(BETA_HDEC(RMB**2/AMH**2))**3
c    .    *SUSY

c     RATCOUP = 0
c     deltaqcd = QCDH(RMB**2/AMH**2)
c     RATCOUP = GHT/XGHB
c     deltat = QCDH(RMB**2/AMH**2) - deltaqcd

c      write(6,*)
c      write(6,*)'H:'
c      write(6,*)'MB,RUNMB,alpha_s: ',AMB,RMB,ASH
c      write(6,*)'MH =              ',AMH
c      write(6,*)'QCD           SUSY        APPROX      APPROX/FULL',
c    .           '  GbH(QCD)    GbH(SQCD):'
c      write(6,*)X1,X2+2*DGHB,2*DGHB,2*DGHB/(X2+2*DGHB),GHB,XGHB
c      write(6,*)'Resummation: ',(XGHB/GHB)**2-1
c      write(6,*)'Rest:        ',SUSY-1
c      write(6,*)'Total SUSY:  ',(XGHB/GHB)**2*SUSY-1
c      write(6,*)'deltaqcd,t = ',deltaqcd,deltat
c      write(6,*)'Gamma(0)   = ',HB2X,HB1X
c      write(6,*)'Gamma(mb)  = ',HB2,HB1
c      write(6,*)'Rest: H      ',AMA,AMH,(SUSY-1)/(X2+2*DGHB)
c      write(52,*)AMA,AMH,(SUSY-1)/(X2+2*DGHB)
      ENDIF
C  H ---> TT
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMH.LE.AMT+AMW+AMB) THEN
        HTT=0.D0
       ELSEIF (AMH.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMH**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT_HDEC(AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        HTT=FACTT*HTT0
       ELSEIF (AMH.LE.XM2) THEN
        ZZMA=AMAR
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL AMHAMA_HDEC(2,XX(1),TGBET)
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT_HDEC(XX(1),AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        YY(1)=FACTT*HTT0
        CALL AMHAMA_HDEC(2,XX(2),TGBET)
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL HTOTT_HDEC(XX(2),AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
        YY(2)=FACTT*HTT0
        CALL AMHAMA_HDEC(2,XX(3),TGBET)
        XMT = RUNM_HDEC(XX(3),6)
        HT1=3.D0*HFF(XX(3),(AMT/XX(3))**2)*GHT**2
     .    *TQCDH(AMT**2/XX(3)**2)
        HT2=3.D0*HFF(XX(3),(XMT/XX(3))**2)*GHT**2
     .    *QCDH(XMT**2/XX(3)**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT_HDEC(RAT,HT1,HT2)
        CALL AMHAMA_HDEC(2,XX(4),TGBET)
        XMT = RUNM_HDEC(XX(4),6)
        HT1=3.D0*HFF(XX(4),(AMT/XX(4))**2)*GHT**2
     .    *TQCDH(AMT**2/XX(4)**2)
        HT2=3.D0*HFF(XX(4),(XMT/XX(4))**2)*GHT**2
     .    *QCDH(XMT**2/XX(4)**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT_HDEC(RAT,HT1,HT2)
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
        HTT=FINT_HDEC(AMH,XX,YY)
       ELSE
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)*GHT**2
     .    *TQCDH(AMT**2/AMH**2)
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)*GHT**2
     .    *QCDH(RMT**2/AMH**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMH
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMH.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT1=3.D0*HFF(AMH,(AMT/AMH)**2)*GHT**2
     .    *TQCDH(AMT**2/AMH**2)
        HT2=3.D0*HFF(AMH,(RMT/AMH)**2)*GHT**2
     .    *QCDH(RMT**2/AMH**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMH
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  H ---> GAMMA GAMMA
       EPS=1.D-8
       XRMC = RUNM_HDEC(AMH/2,4)*AMC/RUNM_HDEC(AMC,4)
       XRMB = RUNM_HDEC(AMH/2,5)*AMB/RUNM_HDEC(AMB,5)
       XRMT = RUNM_HDEC(AMH/2,6)*AMT/RUNM_HDEC(AMT,6)
       CTT = 4*XRMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*XRMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTL = 4*AMTAU**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AMH**2*DCMPLX(1D0,-EPS)
       CTC = 4*XRMC**2/AMH**2*DCMPLX(1D0,-EPS)
       CAC = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GHT
     .     * CFACQ_HDEC(0,AMH,XRMC)
       CAT = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GHT
     .     * CFACQ_HDEC(0,AMH,XRMT)
       CAB = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GHB
     .     * CFACQ_HDEC(0,AMH,XRMB)
       CAL = 1.D0  * 2*CTL*(1+(1-CTL)*CF(CTL))*GHB
       CAW = -(2+3*CTW+3*CTW*(2-CTW)*CF(CTW))*GHVV
       CAH = -AMZ**2/2/AMCH**2*CTH*(1-CTH*CF(CTH))*GHPM
       IF(IOFSUSY.EQ.0) THEN 
        RMSU1 = RUNMS_HDEC(AMH/2,AMSU(1))
        RMSU2 = RUNMS_HDEC(AMH/2,AMSU(2))
        RMSD1 = RUNMS_HDEC(AMH/2,AMSD(1))
        RMSD2 = RUNMS_HDEC(AMH/2,AMSD(2))
        RMSB1 = RUNMS_HDEC(AMH/2,AMSB(1))
        RMSB2 = RUNMS_HDEC(AMH/2,AMSB(2))
        RMST1 = RUNMS_HDEC(AMH/2,AMST(1))
        RMST2 = RUNMS_HDEC(AMH/2,AMST(2))
        CX1 = 4*AMCHAR(1)**2/AMH**2*DCMPLX(1D0,-EPS)
        CX2 = 4*AMCHAR(2)**2/AMH**2*DCMPLX(1D0,-EPS)
        CAX1= AMW/XMCHAR(1) * 2*CX1*(1+(1-CX1)*CF(CX1))*2*AC1(1,1) 
        CAX2= AMW/XMCHAR(2) * 2*CX2*(1+(1-CX2)*CF(CX2))*2*AC1(2,2) 
        CSL1= 4*AMSL(1)**2/AMH**2*DCMPLX(1D0,-EPS)
        CSL2= 4*AMSL(2)**2/AMH**2*DCMPLX(1D0,-EPS)
        CSB1= 4*RMSB1**2/AMH**2*DCMPLX(1D0,-EPS)
        CSB2= 4*RMSB2**2/AMH**2*DCMPLX(1D0,-EPS)
        CST1= 4*RMST1**2/AMH**2*DCMPLX(1D0,-EPS)
        CST2= 4*RMST2**2/AMH**2*DCMPLX(1D0,-EPS)

        CSEL = 4*AMSE(1)**2/AMH**2*DCMPLX(1D0,-EPS)
        CSER = 4*AMSE(2)**2/AMH**2*DCMPLX(1D0,-EPS)
        CSUL = 4*RMSU1**2/AMH**2*DCMPLX(1D0,-EPS)
        CSUR = 4*RMSU2**2/AMH**2*DCMPLX(1D0,-EPS)
        CSDL = 4*RMSD1**2/AMH**2*DCMPLX(1D0,-EPS)
        CSDR = 4*RMSD2**2/AMH**2*DCMPLX(1D0,-EPS)
        CXEL=-2*(-1/2D0+SS)*AMZ**2/AMSE(1)**2*DCOS(A+B)
     .       *CSEL*(1-CSEL*CF(CSEL))
        CXER=2*(SS)*AMZ**2/AMSE(2)**2*DCOS(A+B)
     .       *CSER*(1-CSER*CF(CSER))
        CXUL=-2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .       *AMZ**2/AMSU(1)**2*DCOS(A+B)*CSUL*(1-CSUL*CF(CSUL))
     .      * CFACSQ_HDEC(AMH,RMSU1)
        CXUR=-2*4.D0/3.D0*(2.D0/3.D0*SS)
     .       *AMZ**2/AMSU(2)**2*DCOS(A+B)*CSUR*(1-CSUR*CF(CSUR))
     .      * CFACSQ_HDEC(AMH,RMSU2)
        CXDL=-2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .       *AMZ**2/AMSD(1)**2*DCOS(A+B)*CSDL*(1-CSDL*CF(CSDL))
     .      * CFACSQ_HDEC(AMH,RMSD1)
        CXDR=-2/3.D0*(-1.D0/3.D0*SS)
     .       *AMZ**2/AMSD(2)**2*DCOS(A+B)*CSDR*(1-CSDR*CF(CSDR))
     .      * CFACSQ_HDEC(AMH,RMSD2)

        CXB1= -1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GHBB(1,1)
     .      * CFACSQ_HDEC(AMH,RMSB1)
        CXB2= -1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GHBB(2,2)
     .      * CFACSQ_HDEC(AMH,RMSB2)
        CXT1= -4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GHTT(1,1)
     .      * CFACSQ_HDEC(AMH,RMST1)
        CXT2= -4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GHTT(2,2)
     .      * CFACSQ_HDEC(AMH,RMST2)
        CXL1=       -AMZ**2/AMSL(1)**2*CSL1*(1-CSL1*CF(CSL1))*GHEE(1,1)
        CXL2=       -AMZ**2/AMSL(2)**2*CSL2*(1-CSL2*CF(CSL2))*GHEE(2,2)
        XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .       +  CXEL+CXER+CXUL+CXUR+CXDL+CXDR
     .       +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       ELSE 
        XFAC = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       ENDIF
       HGA=HVV(AMH,0.D0)*(ALPH/PI)**2/16.D0*XFAC
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       XFACQ = CDABS(CAT+CAB+CAC+CAL+CAW+CAH)**2
       XFACS = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXL1+CXL2)**2
       XFACSQ = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXB1+CXB2+CXT1+CXT2+CXL1+CXL2)**2
       HGA0 = HGA*XFACSQ/XFAC
       CAC0 = 4/3D0 * 2*CTC*(1+(1-CTC)*CF(CTC))*GHT
       CAT0 = 4/3D0 * 2*CTT*(1+(1-CTT)*CF(CTT))*GHT
       CAB0 = 1/3D0 * 2*CTB*(1+(1-CTB)*CF(CTB))*GHB
       CXB10= -1/3D0*AMZ**2/AMSB(1)**2*CSB1*(1-CSB1*CF(CSB1))*GHBB(1,1)
       CXB20= -1/3D0*AMZ**2/AMSB(2)**2*CSB2*(1-CSB2*CF(CSB2))*GHBB(2,2)
       CXT10= -4/3D0*AMZ**2/AMST(1)**2*CST1*(1-CST1*CF(CST1))*GHTT(1,1)
       CXT20= -4/3D0*AMZ**2/AMST(2)**2*CST2*(1-CST2*CF(CST2))*GHTT(2,2)
       XFACLOQ = CDABS(CAT0+CAB0+CAC0+CAL+CAW+CAH)**2
       CXUL0=-2*4.D0/3.D0*(1.D0/2.D0-2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(1)**2*DCOS(A+B)*CSUL*(1-CSUL*CF(CSUL))
       CXUR0=-2*4.D0/3.D0*(2.D0/3.D0*SS)
     .      *AMZ**2/AMSU(2)**2*DCOS(A+B)*CSUR*(1-CSUR*CF(CSUR))
       CXDL0=-2/3.D0*(-1.D0/2.D0+1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(1)**2*DCOS(A+B)*CSDL*(1-CSDL*CF(CSDL))
       CXDR0=-2/3.D0*(-1.D0/3.D0*SS)
     .      *AMZ**2/AMSD(2)**2*DCOS(A+B)*CSDR*(1-CSDR*CF(CSDR))
       XFACLO = CDABS(CAT0+CAB0+CAC0+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXEL+CXER+CXUL0+CXUR0+CXDL0+CXDR0
     .      +  CXB10+CXB20+CXT10+CXT20+CXL1+CXL2)**2
       CSQ = 1+3*ALPHAS_HDEC(AMH,2)
       XFACSQL = CDABS(CAT+CAB+CAC+CAL+CAW+CAH+CAX1+CAX2
     .      +  CXEL+CXER+(CXUL0+CXUR0+CXDL0+CXDR0
     .      +  CXB10+CXB20+CXT10+CXT20)*CSQ+CXL1+CXL2)**2
c      write(55,('6(1X,E12.6)'))AMH,HGA0,HGA0*XFACQ/XFACSQ,
c    .       XFACSQ/XFACLO-1,XFACQ/XFACLOQ-1,(XFACSQ-XFACSQL)/XFACSQL
c      write(55,('7(1X,E12.6)'))AMH,HGA,HGA*XFACQ/XFAC,HGA*XFACS/XFAC,
c    .       XFAC/XFACLO-1,XFACQ/XFACLOQ-1,(XFAC-XFACSQL)/XFACSQL
c      write(6,*)AMCH,AMST,AMSB,AMSL,AMCHAR
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  H ---> Z GAMMA
      IF(AMH.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GHT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GHB
       EPS=1.D-8
       CTT = 4*AMT**2/AMH**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMH**2*DCMPLX(1D0,-EPS)
       CTW = 4*AMW**2/AMH**2*DCMPLX(1D0,-EPS)
       CTH = 4*AMCH**2/AMH**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLW = 4*AMW**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLH = 4*AMCH**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(CI1(CTT,CLT) - CI2(CTT,CLT))
       CAB = FB*(CI1(CTB,CLB) - CI2(CTB,CLB))
       CAW = -1/DSQRT(TS)*(4*(3-TS)*CI2(CTW,CLW)
     .     + ((1+2/CTW)*TS - (5+2/CTW))*CI1(CTW,CLW))*GHVV
       CAH = (1-2*SS)/DSQRT(SS*CS)*AMZ**2/2/AMCH**2*CI1(CTH,CLH)*GHPM
       XFAC = CDABS(CAT+CAB+CAW+CAH)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMH**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMH**2)**3
      ENDIF
C  H ---> W W
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMW-DLD
       XM2 = 2D0*AMW+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV_HDEC(AMH,AMW,GAMW,HTWW)
        HWW = 3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/AMH**3*HTWW*GHVV**2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMW,GAMW,HTWW)
        YY(1)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(1)**3*HTWW
        CALL HTOVV_HDEC(XX(2),AMW,GAMW,HTWW)
        YY(2)=3D0/2D0*GF*AMW**4/DSQRT(2D0)/PI/XX(2)**3*HTWW
        YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
        YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
        HWW = FINT_HDEC(AMH,XX,YY)*GHVV**2
       ELSE
        HWW=HVV(AMH,AMW**2/AMH**2)*GHVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMW-DLD
      XM2 = 2D0*AMW+DLU
      IF (AMH.LE.AMW) THEN
       HWW=0
      ELSE IF (AMH.LE.XM1) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       HWW=HV(AMW**2/AMH**2)*CWW*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
       CWW=3.D0*GF**2*AMW**4/16.D0/PI**3
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMW**2/XX(1)**2)*CWW*XX(1)
       YY(2)=HV(AMW**2/XX(2)**2)*CWW*XX(2)
       YY(3)=HVV(XX(3),AMW**2/XX(3)**2)
       YY(4)=HVV(XX(4),AMW**2/XX(4)**2)
       HWW = FINT_HDEC(AMH,XX,YY)*GHVV**2
      ELSE
       HWW=HVV(AMH,AMW**2/AMH**2)*GHVV**2
      ENDIF
      ENDIF
C  H ---> Z Z
      IF(IONWZ.EQ.0)THEN
       DLD=2D0
       DLU=2D0
       XM1 = 2D0*AMZ-DLD
       XM2 = 2D0*AMZ+DLU
       IF (AMH.LE.XM1) THEN
        CALL HTOVV_HDEC(AMH,AMZ,GAMZ,HTZZ)
        HZZ = 3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/AMH**3*HTZZ*GHVV**2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        CALL HTOVV_HDEC(XX(1),AMZ,GAMZ,HTZZ)
        YY(1)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(1)**3*HTZZ
        CALL HTOVV_HDEC(XX(2),AMZ,GAMZ,HTZZ)
        YY(2)=3D0/4D0*GF*AMZ**4/DSQRT(2D0)/PI/XX(2)**3*HTZZ
        YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2
        YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2
        HZZ = FINT_HDEC(AMH,XX,YY)*GHVV**2
       ELSE
        HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0*GHVV**2
       ENDIF
      ELSE
      DLD=2D0
      DLU=2D0
      XM1 = 2D0*AMZ-DLD
      XM2 = 2D0*AMZ+DLU
      IF (AMH.LE.AMZ) THEN
       HZZ=0
      ELSE IF (AMH.LE.XM1) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       HZZ=HV(AMZ**2/AMH**2)*CZZ*AMH*GHVV**2
      ELSE IF (AMH.LT.XM2) THEN
       CZZ=3.D0*GF**2*AMZ**4/192.D0/PI**3*(7-40/3.D0*SS+160/9.D0*SS**2)
       XX(1) = XM1-1D0
       XX(2) = XM1
       XX(3) = XM2
       XX(4) = XM2+1D0
       YY(1)=HV(AMZ**2/XX(1)**2)*CZZ*XX(1)
       YY(2)=HV(AMZ**2/XX(2)**2)*CZZ*XX(2)
       YY(3)=HVV(XX(3),AMZ**2/XX(3)**2)/2D0
       YY(4)=HVV(XX(4),AMZ**2/XX(4)**2)/2D0
       HZZ = FINT_HDEC(AMH,XX,YY)*GHVV**2
      ELSE
       HZZ=HVV(AMH,AMZ**2/AMH**2)/2.D0*GHVV**2
      ENDIF
      ENDIF
C  H ---> h h
      IF(IONSH.EQ.0)THEN
      if(islhai.eq.0)then
       ZZMA = AMAR
       AMREAL = AMH
       AMA = 1.D0
       AMLOW = AMH
12345  CALL SUSYCP_HDEC(TGBET)
       IF(AMLR.LT.0.D0)THEN
        AMA = AMAR + 1
        GOTO 12345
       ENDIF
       AMLOW = AMH
       AMDEL = AMREAL - AMLOW
       DLD = 0.3D0*(TGBET-1.3D0)
       DLD = DMAX1(0.1D0,DLD)
       DLU=DLD
       AMA = ZZMA
       CALL SUSYCP_HDEC(TGBET)
       XM1 = 2*AML-DLD
       XM2 = 2*AML+DLU
       IF (AMH.LE.AML) THEN
        HHH=0
       ELSEIF (AMH.LT.XM1) THEN
        XH=AML**2/AMH**2
        XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .    *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
        XH2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GHLL**2*GLB**2*AMB**2
        HHH=XH1*XH2
       ELSEIF (AMH.LT.XM2) THEN
        IFLON0 = 0
        IFLON1 = 0
        ZZMA=AMAR
        AMACRIT = AMAR
        AMA0 = AMAR
        AMA1 = AMAR
510     AMA0 = AMA0 - 1
        AMA1 = AMA1 + 1
        AMA = AMA0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.2*AML) THEN
         IFLON0 = -1
        ELSE
         IFLON0 = 1
        ENDIF
        AMA = AMA1
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.2*AML) THEN
         IFLON1 = -1
        ELSE
         IFLON1 = 1
        ENDIF
        IF(IFLON0*IFLON1.NE.-1) GOTO 510
501     AMA = (AMA0+AMA1)/2
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.2*AML) THEN
         IF(IFLON0.EQ.-1) THEN
          AMA0 = AMAR
         ELSE
          AMA1 = AMAR
         ENDIF
        ELSE
         IF(IFLON0.EQ.-1) THEN
          AMA1 = AMAR
         ELSE
          AMA0 = AMAR
         ENDIF
        ENDIF
        AMACRIT = (AMA0+AMA1)/2
        DEL = 1.D-8
        AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
        IF(AMDEL.GT.DEL) GOTO 501
       AMA = AMACRIT
       CALL SUSYCP_HDEC(TGBET)
       YM1 = AMACRIT
       YM2 = AMACRIT
       AMA0 = AMACRIT
       AMA1 = AMACRIT
       DELSTEP = 1.D0
511    AMA0 = AMA0 - DELSTEP
       AMA1 = AMA1 + DELSTEP
       AMA = AMACRIT
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLONC = -1
       ELSE
        IFLONC = 1
       ENDIF
       AMA = AMA0
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLON0 = -1
       ELSE
        IFLON0 = 1
       ENDIF
       AMA = AMA1
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IFLON1 = -1
       ELSE
        IFLON1 = 1
       ENDIF
       IF(IFLON0*IFLONC.NE.-1.AND.IFLONC*IFLON1.NE.-1) GOTO 511
       IF(IFLON0*IFLONC.EQ.-1) THEN
         AMA1 = AMACRIT
         IFLON1 = IFLONC
       ELSE
         AMA0 = AMACRIT
         IFLON0 = IFLONC
       ENDIF
512    AMA = (AMA0+AMA1)/2
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML-DLD) THEN
        IF(IFLON0.EQ.-1) THEN
         AMA0 = AMAR
        ELSE
         AMA1 = AMAR
        ENDIF
       ELSE
        IF(IFLON0.EQ.-1) THEN
         AMA1 = AMAR
        ELSE
         AMA0 = AMAR
        ENDIF
       ENDIF
       YM1 = (AMA0+AMA1)/2
       DEL = 1.D-8
       AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
       IF(AMDEL.GT.DEL) GOTO 512
       AMA = YM1
       CALL SUSYCP_HDEC(TGBET)
       AMA0 = AMACRIT
       AMA1 = AMACRIT
       DELSTEP = 1.D0
513    AMA0 = AMA0 - DELSTEP
       AMA1 = AMA1 + DELSTEP
       AMA = AMACRIT
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLONC = -1
       ELSE
        IFLONC = 1
       ENDIF
       AMA = AMA0
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLON0 = -1
       ELSE
        IFLON0 = 1
       ENDIF
       AMA = AMA1
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IFLON1 = -1
       ELSE
        IFLON1 = 1
       ENDIF
       IF(IFLON0*IFLONC.NE.-1.AND.IFLONC*IFLON1.NE.-1) GOTO 513
       IF(IFLON0*IFLONC.EQ.-1) THEN
         AMA1 = AMACRIT
         IFLON1 = IFLONC
       ELSE
         AMA0 = AMACRIT
         IFLON0 = IFLONC
       ENDIF
514    AMA = (AMA0+AMA1)/2
       CALL SUSYCP_HDEC(TGBET)
       IF(AMH.LT.2*AML+DLU) THEN
        IF(IFLON0.EQ.-1) THEN
         AMA0 = AMAR
        ELSE
         AMA1 = AMAR
        ENDIF
       ELSE
        IF(IFLON0.EQ.-1) THEN
         AMA1 = AMAR
        ELSE
         AMA0 = AMAR
        ENDIF
       ENDIF
       YM2 = (AMA0+AMA1)/2
       DEL = 1.D-8
       AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
       IF(AMDEL.GT.DEL) GOTO 514
       AMA = YM2
       CALL SUSYCP_HDEC(TGBET)
       DEL = 1.D-4
        XX(1) = YM1 - DEL
        XX(2) = YM1
        XX(3) = YM2
        XX(4) = YM2 + DEL
        AMAR = ZZMA
        DO J=1,4
         AMA = XX(J)
         CALL SUSYCP_HDEC(TGBET)
         XX(J) = AMH
         IF(AMH.GE.2*AML)THEN
          YY(J)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(J)
     .          *BETA_HDEC(AML**2/XX(J)**2)
         ELSEIF(AMH.LE.AML)THEN
          YY(J) = 0
         ELSE
          XH=AML**2/XX(J)**2
          XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .    *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
          XH2=3*GF**2/16.D0/PI**3*AMZ**4/XX(J)*GLB**2*AMB**2
          YY(J)=XH1*XH2
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
        HHH = FINT_HDEC(AMH,XX,YY)*GHLL**2
       ELSE
        HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AML**2/AMH**2)
     .      *GHLL**2
       ENDIF
      else
       DLD=0.1D0
       DLU=0.1D0
       XM1 = 2D0*AML-DLD
       XM2 = 2D0*AML+DLU
       IF (AMH.LE.AML) THEN
        HHH = 0D0
       ELSEIF (AMH.LE.XM1) THEN
        XH=AML**2/AMH**2
        XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .  *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .   -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
        XH2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GLB**2*AMB**2
        HHH=XH1*XH2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        XH=AML**2/XX(1)**2
        XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .  *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .   -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
        XH2=3*GF**2/16.D0/PI**3*AMZ**4/XX(1)*GLB**2*AMB**2
        YY(1)=XH1*XH2
        XH=AML**2/XX(1)**2
        XH1=(XH-1.D0)*(2.D0-.5D0*DLOG(XH))+(1.D0-5.D0*XH)
     .  *(DATAN((2.D0*XH-1.D0)/DSQRT(4.D0*XH-1.D0))
     .   -DATAN(1.D0/DSQRT(4.D0*XH-1.D0)))/DSQRT(4.D0*XH-1.D0)
        XH2=3*GF**2/16.D0/PI**3*AMZ**4/XX(2)*GLB**2*AMB**2
        YY(2)=XH1*XH2
        YY(3)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(3)
     .        *BETA_HDEC(AML**2/XX(3)**2)
        YY(4)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(4)
     .        *BETA_HDEC(AML**2/XX(4)**2)
        HHH = FINT_HDEC(AMH,XX,YY)*GHLL**2
       ELSE
        HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AML**2/AMH**2)
     .      *GHLL**2
       ENDIF
      endif
      ELSE
       IF (AMH.LE.2*AML) THEN
        HHH=0
       ELSE
        HHH=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AML**2/AMH**2)
     .      *GHLL**2
       ENDIF
      ENDIF
C  H ---> A A
      IF(IONSH.EQ.0)THEN
      if(islhai.eq.0)then
       DLD = 0.3D0*(TGBET-1.3D0)
       DLD = DMAX1(0.1D0,DLD)
       DLU=DLD
       ALD = DLD/2
       ALU = DLU/2
       XM1 = 2*AMA-DLD
       XM2 = 2*AMA+DLU
       IF (AMH.LE.AMA) THEN
        HAA=0
       ELSEIF (AMH.LT.XM1) THEN
        XA=AMA**2/AMH**2
        XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .    *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
        XA2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GHAA**2*GAB**2*AMB**2
        HAA=XA1*XA2
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
        AMACRIT = AMAR
        AMA0 = 10.D0
        AMA1 = AMAR + 50.D0
        AMA = AMA0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.2*AMA) THEN
         IFLON0 = -1
        ELSEIF(AMH.EQ.2*AMA) THEN
         IFLON0 = 0
         AMACRIT = AMAR
        ELSE
         IFLON0 = 1
        ENDIF
        AMA = AMA1
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.2*AMA) THEN
         IFLON1 = -1
        ELSEIF(AMH.EQ.2*AMA) THEN
         IFLON1 = 0
         AMACRIT = AMAR
        ELSE
         IFLON1 = 1
        ENDIF
        IF(IFLON0*IFLON1.EQ.0)THEN
         IFLON0 = 0
         IFLON1 = 0
        ENDIF
        IF(IFLON0.NE.IFLON1)THEN
502      AMA = (AMA0+AMA1)/2
         CALL SUSYCP_HDEC(TGBET)
         IF(AMH.LT.2*AMA) THEN
          IF(IFLON0.EQ.-1) THEN
           AMA0 = AMAR
          ELSE
           AMA1 = AMAR
          ENDIF
         ELSEIF(AMH.EQ.2*AMA) THEN
          IFLON0 = 0
          IFLON1 = 0
          AMACRIT = AMAR
         ELSE
          IF(IFLON0.EQ.-1) THEN
           AMA1 = AMAR
          ELSE
           AMA0 = AMAR
          ENDIF
         ENDIF
         IF(IFLON0.NE.0)THEN
          AMACRIT = (AMA0+AMA1)/2
          DEL = 1.D-8
          AMDEL = 2*DABS(AMA1-AMA0)/(AMA1+AMA0)
          IF(AMDEL.GT.DEL) GOTO 502
         ENDIF
        ENDIF
        DEL = 1.D-4
        XX(1) = AMACRIT - ALD - DEL
        XX(2) = AMACRIT - ALD
        XX(3) = AMACRIT + ALU
        XX(4) = AMACRIT + ALU + DEL
        DO J=1,4
         AMA = XX(J)
         CALL SUSYCP_HDEC(TGBET)
         XX(J) = AMH
         IF(AMH.GE.2*AMA)THEN
          YY(J)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(J)
     .          *BETA_HDEC(AMA**2/XX(J)**2)
         ELSEIF(AMH.LE.AMA)THEN
          YY(J) = 0
         ELSE
          XA=AMA**2/XX(J)**2
          XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .    *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
          XA2=3*GF**2/16.D0/PI**3*AMZ**4/XX(J)*GAB**2*AMB**2
          YY(J)=XA1*XA2
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
        HAA = FINT_HDEC(AMH,XX,YY)*GHAA**2
       ELSE
        HAA=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AMA**2/AMH**2)
     .       *GHAA**2
       ENDIF
      else
       DLD=0.1D0
       DLU=0.1D0
       XM1 = 2D0*AMA-DLD
       XM2 = 2D0*AMA+DLU
       IF (AMH.LE.AMA) THEN
        HAA = 0D0
       ELSEIF (AMH.LE.XM1) THEN
        XA=AMA**2/AMH**2
        XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .    *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .     -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
        XA2=3*GF**2/16.D0/PI**3*AMZ**4/AMH*GHAA**2*GAB**2*AMB**2
        HAA=XA1*XA2
       ELSEIF (AMH.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        XA=AMA**2/XX(1)**2
        XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .  *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .   -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
        XA2=3*GF**2/16.D0/PI**3*AMZ**4/XX(1)*GAB**2*AMB**2
        YY(1)=XA1*XA2
        XA=AMA**2/XX(2)**2
        XA1=(XA-1.D0)*(2.D0-.5D0*DLOG(XA))+(1.D0-5.D0*XA)
     .  *(DATAN((2.D0*XA-1.D0)/DSQRT(4.D0*XA-1.D0))
     .   -DATAN(1.D0/DSQRT(4.D0*XA-1.D0)))/DSQRT(4.D0*XA-1.D0)
        XA2=3*GF**2/16.D0/PI**3*AMZ**4/XX(2)*GAB**2*AMB**2
        YY(2)=XA1*XA2
        YY(3)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(3)
     .        *BETA_HDEC(AMA**2/XX(3)**2)
        YY(4)=GF/16D0/DSQRT(2D0)/PI*AMZ**4/XX(4)
     .        *BETA_HDEC(AMA**2/XX(4)**2)
        HAA = FINT_HDEC(AMH,XX,YY)*GHAA**2
       ELSE
        HAA=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AMA**2/AMH**2)
     .       *GHAA**2
       ENDIF
      endif
      ELSE
       IF (AMH.LE.2*AMA) THEN
        HAA=0
       ELSE
        HAA=GF/16D0/DSQRT(2D0)/PI*AMZ**4/AMH*BETA_HDEC(AMA**2/AMH**2)
     .       *GHAA**2
       ENDIF
      ENDIF
C  H ---> A Z
      IF(IONSH.EQ.0)THEN
      if(islhai.eq.0)then
       DLD=1D0
       DLU=8D0
       XM1 = AMA+AMZ-DLD
       XM2 = AMA+AMZ+DLU
       IF (AMH.LT.AMA) THEN
        HAZ=0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMZ-AMA))THEN
         HAZ=0
        ELSE
        HAZ=9.D0*GF**2/8.D0/PI**3*AMZ**4*AMH*GZAH**2*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AMA/AMH)**2,(AMZ/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
165     AMA = AMAR - 1.D0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.AMA+AMZ+DLU.AND.AMH.GT.AMA+AMZ-DLD) GOTO 165
        XX(1) = AMAR-1D0
        XX(2) = AMAR
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
166     AMA = AMAR + 1.D0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.AMA+AMZ+DLU.AND.AMH.GT.AMA+AMZ-DLD) GOTO 166
        XX(3) = AMAR
        XX(4) = AMAR+1D0
        DO IJ=1,4
         AMA = XX(IJ)
         CALL SUSYCP_HDEC(TGBET)
         XX(IJ) = AMH
         IF(AMH.LE.AMA+AMZ) THEN
          YY(IJ)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(IJ)*
     .          (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .          *HVH((AMA/XX(IJ))**2,(AMZ/XX(IJ))**2)
         ELSE
          CAZ=LAMB_HDEC(AMA**2/XX(IJ)**2,AMZ**2/XX(IJ)**2)
     .       *LAMB_HDEC(XX(IJ)**2/AMZ**2,AMA**2/AMZ**2)**2
          YY(IJ)=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/XX(IJ)*CAZ
         ENDIF
        ENDDO
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
        HAZ = FINT_HDEC(AMH,XX,YY)*GZAH**2
       ELSE
        CAZ=LAMB_HDEC(AMA**2/AMH**2,AMZ**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
        HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
       ENDIF
      else
       DLD=1D0
       DLU=8D0
       XM1 = AMA+AMZ-DLD
       XM2 = AMA+AMZ+DLU
       IF (AMH.LT.AMA) THEN
        HAZ=0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMZ-AMA))THEN
         HAZ=0
        ELSE
        HAZ=9.D0*GF**2/8.D0/PI**3*AMZ**4*AMH*GZAH**2*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AMA/AMH)**2,(AMZ/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(1)*
     .        (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .        *HVH((AMA/XX(1))**2,(AMZ/XX(1))**2)
        YY(2)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(2)*
     .        (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .        *HVH((AMA/XX(2))**2,(AMZ/XX(2))**2)
        CAZ=LAMB_HDEC(AMA**2/XX(3)**2,AMZ**2/XX(3)**2)
     .     *LAMB_HDEC(XX(3)**2/AMZ**2,AMA**2/AMZ**2)**2
        YY(3)=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/XX(3)*CAZ
        CAZ=LAMB_HDEC(AMA**2/XX(4)**2,AMZ**2/XX(4)**2)
     .     *LAMB_HDEC(XX(4)**2/AMZ**2,AMA**2/AMZ**2)**2
        YY(4)=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/XX(4)*CAZ
        HAZ = FINT_HDEC(AMH,XX,YY)*GZAH**2
       ELSE
        CAZ=LAMB_HDEC(AMA**2/AMH**2,AMZ**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
        HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
       ENDIF
      endif
      ELSE
       IF (AMH.LT.AMZ+AMA) THEN
        HAZ=0
       ELSE
        CAZ=LAMB_HDEC(AMA**2/AMH**2,AMZ**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMZ**2,AMA**2/AMZ**2)**2
        HAZ=GF/8.D0/DSQRT(2D0)/PI*AMZ**4/AMH*CAZ*GZAH**2
       ENDIF
      ENDIF
C  H ---> H+ W+
      IF(IONSH.EQ.0)THEN
      if(islhai.eq.0)then
       DLD=3D0
       DLU=9D0
       XM1 = AMCH+AMW-DLD
       XM2 = AMCH+AMW+DLU
       IF (AMH.LT.AMCH) THEN
        HHW=0.D0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMW-AMCH))THEN
         HHW=0
        ELSE
        HHW=9.D0*GF**2/8.D0/PI**3*AMW**4*AMH*GLVV**2*2
     .      *HVH((AMCH/AMH)**2,(AMW/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        ZZMA=AMAR
167     AMA = AMAR - 1.D0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.LT.AMCH+AMW+DLU) GOTO 167
        XX(1) = AMAR-1D0
        XX(2) = AMAR
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
168     AMA = AMAR + 1.D0
        CALL SUSYCP_HDEC(TGBET)
        IF(AMH.GT.AMCH+AMW-DLD) GOTO 168
        XX(3) = AMAR
        XX(4) = AMAR+1D0
        AMA = XX(1)
        CALL SUSYCP_HDEC(TGBET)
        XX(1) = AMH
        CHW=LAMB_HDEC(AMCH**2/XX(1)**2,AMW**2/XX(1)**2)
     .     *LAMB_HDEC(XX(1)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(1)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(1)*CHW
        AMA = XX(2)
        CALL SUSYCP_HDEC(TGBET)
        XX(2) = AMH
        CHW=LAMB_HDEC(AMCH**2/XX(2)**2,AMW**2/XX(2)**2)
     .     *LAMB_HDEC(XX(2)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(2)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(2)*CHW
        AMA = XX(3)
        CALL SUSYCP_HDEC(TGBET)
        XX(3) = AMH
        YY(3)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(3)*2
     .       *HVH((AMCH/XX(3))**2,(AMW/XX(3))**2)
        AMA = XX(4)
        CALL SUSYCP_HDEC(TGBET)
        XX(4) = AMH
        YY(4)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(4)*2
     .       *HVH((AMCH/XX(4))**2,(AMW/XX(4))**2)
        AMA = ZZMA
        CALL SUSYCP_HDEC(TGBET)
        HHW=FINT_HDEC(AMH,XX,YY)*GLVV**2
       ELSE
        CHW=LAMB_HDEC(AMCH**2/AMH**2,AMW**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMW**2,AMCH**2/AMW**2)**2
        HHW=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMH*CHW*GLVV**2
       ENDIF
      else
       DLD=3D0
       DLU=9D0
       XM1 = AMCH+AMW-DLD
       XM2 = AMCH+AMW+DLU
       IF (AMH.LT.AMCH) THEN
        HHW=0.D0
       ELSEIF (AMH.LT.XM1) THEN
        IF(AMH.LE.DABS(AMW-AMCH))THEN
         HHW=0
        ELSE
        HHW=9.D0*GF**2/8.D0/PI**3*AMW**4*AMH*GLVV**2*2
     .      *HVH((AMCH/AMH)**2,(AMW/AMH)**2)
        ENDIF
       ELSEIF (AMH.LT.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(1)*2
     .       *HVH((AMCH/XX(1))**2,(AMW/XX(1))**2)
        YY(2)=9.D0*GF**2/8.D0/PI**3*AMW**4*XX(2)*2
     .       *HVH((AMCH/XX(2))**2,(AMW/XX(2))**2)
        CHW=LAMB_HDEC(AMCH**2/XX(3)**2,AMW**2/XX(3)**2)
     .     *LAMB_HDEC(XX(3)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(3)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(3)*CHW
        CHW=LAMB_HDEC(AMCH**2/XX(4)**2,AMW**2/XX(4)**2)
     .     *LAMB_HDEC(XX(4)**2/AMW**2,AMCH**2/AMW**2)**2
        YY(4)=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/XX(4)*CHW
        HHW=FINT_HDEC(AMH,XX,YY)*GLVV**2
       ELSE
        CHW=LAMB_HDEC(AMCH**2/AMH**2,AMW**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMW**2,AMCH**2/AMW**2)**2
        HHW=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMH*CHW*GLVV**2
       ENDIF
      endif
      ELSE
       IF (AMH.LT.AMW+AMCH) THEN
        HHW=0.D0
       ELSE
        CHW=LAMB_HDEC(AMCH**2/AMH**2,AMW**2/AMH**2)
     .     *LAMB_HDEC(AMH**2/AMW**2,AMCH**2/AMW**2)**2
        HHW=2*GF/8.D0/DSQRT(2D0)/PI*AMW**4/AMH*CHW*GLVV**2
       ENDIF
      ENDIF

C ========================== SUSY DECAYS 
C
      IF(IOFSUSY.EQ.0) THEN
C  HH ----> CHARGINOS
      DO 741 I=1,2
      DO 741 J=1,2
      IF (AMH.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHHCH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMH 
     .     *LAMB_HDEC(AMCHAR(I)**2/AMH**2,AMCHAR(J)**2/AMH**2)
     .     *( (AC1(I,J)**2+AC1(J,I)**2)*(AMH**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)-4.D0*AC1(I,J)*AC1(J,I)* 
     .         XMCHAR(I)*XMCHAR(J) ) 
      ELSE
      WHHCH(I,J)=0.D0
      ENDIF
 741  CONTINUE
      WHHCHT=WHHCH(1,1)+WHHCH(1,2)+WHHCH(2,1)+WHHCH(2,2)
C
C  HH ----> NEUTRALINOS 
      DO 742 I=1,4
      DO 742 J=1,4
      IF (AMH.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHHNE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMH
     .         *AN1(I,J)**2*(AMH**2-(XMNEUT(I)+XMNEUT(J))**2)
     .         *LAMB_HDEC(AMNEUT(I)**2/AMH**2,AMNEUT(J)**2/AMH**2)
      ELSE 
      WHHNE(I,J)=0.D0
      ENDIF
 742  CONTINUE
      WHHNET= WHHNE(1,1)+WHHNE(1,2)+WHHNE(1,3)+WHHNE(1,4)
     .       +WHHNE(2,1)+WHHNE(2,2)+WHHNE(2,3)+WHHNE(2,4)
     .       +WHHNE(3,1)+WHHNE(3,2)+WHHNE(3,3)+WHHNE(3,4)
     .       +WHHNE(4,1)+WHHNE(4,2)+WHHNE(4,3)+WHHNE(4,4)
C
C  HH ----> SLEPTONS 
C
      IF (AMH.GT.2.D0*AMSE(1)) THEN
      WHHSLEL=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSE(1)**2/AMH**2)*(-0.5D0+SS)**2
      ELSE
      WHHSLEL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSE(2)) THEN
      WHHSLER=2*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSE(2)**2/AMH**2)*SS**2
      ELSE
      WHHSLER=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSN(1)) THEN
      WHHSLNL=3*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSN(1)**2/AMH**2)*0.5D0**2
      ELSE
      WHHSLNL=0.D0
      ENDIF

      DO 748 I=1,2
      DO 748 J=1,2
      IF(AMH.GT.AMSL(I)+AMSL(J)) THEN
      WHHSTAU(I,J)=GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHEE(I,J)**2*
     .      LAMB_HDEC(AMSL(I)**2/AMH**2,AMSL(J)**2/AMH**2)/AMH
      ELSE
      WHHSTAU(I,J)=0.D0
      ENDIF
 748  CONTINUE

      WHHSLT=WHHSTAU(1,1)+WHHSTAU(1,2)+WHHSTAU(2,1)+WHHSTAU(2,2) 
     .       +WHHSLEL+WHHSLER+WHHSLNL
C
C  HH ----> SQUARKS 
C
      IF (AMH.GT.2.D0*AMSU(1)) THEN
      WHHSQUL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSU(1)**2/AMH**2)*(0.5D0-2.D0/3.D0*SS)**2
      ELSE
      WHHSQUL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSU(2)) THEN
      WHHSQUR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSU(2)**2/AMH**2)*(-2.D0/3.D0*SS)**2
      ELSE
      WHHSQUR=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSD(1)) THEN
      WHHSQDL=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSD(1)**2/AMH**2)*(-0.5D0+1.D0/3.D0*SS)**2
      ELSE
      WHHSQDL=0.D0
      ENDIF

      IF (AMH.GT.2.D0*AMSD(2)) THEN
      WHHSQDR=6*GF/2.D0/DSQRT(2D0)/PI*AMZ**4/AMH*DCOS(B+A)**2
     .      *BETA_HDEC(AMSD(2)**2/AMH**2)*(+1.D0/3.D0*SS)**2
      ELSE
      WHHSQDR=0.D0
      ENDIF

      WHHSQ=WHHSQUL+WHHSQUR+WHHSQDL+WHHSQDR
C
C  HH ----> STOPS 
      SUSY = 1
      DO 743 I=1,2
      DO 743 J=1,2
      IF(AMH.GT.AMST(I)+AMST(J)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(2,1,I,J,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHHST(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHTT(I,J)**2*
     .       LAMB_HDEC(AMST(I)**2/AMH**2,AMST(J)**2/AMH**2)/AMH
     .      *SUSY
c      write(6,*)'H -> stop: ',I,J,AMH,AMST(I),AMST(J),SUSY-1,
c    .           WHHST(I,J)/SUSY,WHHST(I,J)
c      write(6,*)'H -> stop: ',I,J,AMH,AMST(I),AMST(J),SUSY-1
      ELSE
      WHHST(I,J)=0.D0
      ENDIF
 743  CONTINUE
C
C  HH ----> SBOTTOMS 
      SUSY = 1
      DO 744 I=1,2
      DO 744 J=1,2
      IF(AMH.GT.AMSB(I)+AMSB(J)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(2,2,I,J,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHHSB(I,J)=3*GF*AMZ**4/2.D0/DSQRT(2.D0)/PI*GHBB(I,J)**2*
     .       LAMB_HDEC(AMSB(I)**2/AMH**2,AMSB(J)**2/AMH**2)/AMH
     .      *SUSY
       write(6,*)'H -> sbot: ',I,J,AMH,AMSB(I),AMSB(J),SUSY-1,
     .           WHHSB(I,J)/SUSY,WHHSB(I,J)
c      write(6,*)'H -> sbot: ',I,J,AMH,AMSB(I),AMSB(J),SUSY-1
      ELSE
      WHHSB(I,J)=0.D0
      ENDIF
 744  CONTINUE
C
      WHHSTT=WHHST(1,1)+WHHST(1,2)+WHHST(2,1)+WHHST(2,2) 
      WHHSBB=WHHSB(1,1)+WHHSB(1,2)+WHHSB(2,1)+WHHSB(2,2) 
      WHHSQT=WHHSTT+WHHSBB+WHHSQ
C
      ELSE 
      WHHCHT=0.D0
      WHHNET=0.D0
      WHHSLT=0.D0
      WHHSQT=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHHCH(I,J)=0.D0
        WHHST(I,J)=0.D0
        WHHSB(I,J)=0.D0
        WHHSTAU(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHHNE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF

      IF(IGOLD.NE.0)THEN
C   HH ---> GOLDSTINOS
       DO 740 I=1,4
       IF (AMH.GT.AMNEUT(I)) THEN
        WHHGD(I)=AMH**5/AXMPL**2/AXMGD**2/48.D0/PI*
     .           (1.D0-AMNEUT(I)**2/AMH**2)**4*AGDH(I)**2
       ELSE
        WHHGD(I)=0.D0
       ENDIF
 740   CONTINUE
       WHHGDT=WHHGD(1)+WHHGD(2)+WHHGD(3)+WHHGD(4)
      ELSE
       WHHGDT=0
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HTT+HGG+HGA+HZGA+HWW+HZZ+HHH+HAA+HAZ
     .    +HHW+WHHCHT+WHHNET+WHHSLT+WHHSQT + WHHGDT
      HHBRT=HTT/WTOT
      HHBRB=HBB/WTOT
      HHBRL=HLL/WTOT
      HHBRM=HMM/WTOT
      HHBRS=HSS/WTOT
      HHBRC=HCC/WTOT
      HHBRG=HGG/WTOT
      HHBRGA=HGA/WTOT
      HHBRZGA=HZGA/WTOT
      HHBRW=HWW/WTOT
      HHBRZ=HZZ/WTOT
      HHBRH=HHH/WTOT
      HHBRA=HAA/WTOT
      HHBRAZ=HAZ/WTOT
      HHBRHW=HHW/WTOT
      DO 841 I=1,2
      DO 841 J=1,2
      HHBRSC(I,J)=WHHCH(I,J)/WTOT
841   CONTINUE
      DO 842 I=1,4
      DO 842 J=1,4
      HHBRSN(I,J)=WHHNE(I,J)/WTOT
842   CONTINUE
      HHBRCHT=WHHCHT/WTOT 
      HHBRNET=WHHNET/WTOT 
      HHBRSL=WHHSLT/WTOT
      HHBRSQ=WHHSQ/WTOT
      HHBRSQT=WHHSQT/WTOT
      DO 843 I=1,2
      DO 843 J=1,2
      HHBRST(I,J)=WHHST(I,J)/WTOT
843   CONTINUE
      DO 844 I=1,2
      DO 844 J=1,2
      HHBRSB(I,J)=WHHSB(I,J)/WTOT
844   CONTINUE
      HHBRGD =WHHGDT/WTOT
      HHWDTH=WTOT

      BHHSLNL = WHHSLNL/WTOT
      BHHSLEL = WHHSLEL/WTOT
      BHHSLER = WHHSLER/WTOT
      BHHSQUL = WHHSQUL/WTOT
      BHHSQUR = WHHSQUR/WTOT
      BHHSQDL = WHHSQDL/WTOT
      BHHSQDR = WHHSQDR/WTOT
      DO I = 1,2
       DO J = 1,2
        BHHST(I,J) = WHHST(I,J)/WTOT
        BHHSB(I,J) = WHHSB(I,J)/WTOT
        BHHSTAU(I,J) = WHHSTAU( I,J)/WTOT
       ENDDO
      ENDDO

      ENDIF

      IF(IHIGGS.EQ.3.OR.IHIGGS.EQ.5)THEN 
C
C        =========================================================
C                       CP ODD  HIGGS DECAYS
C        =========================================================
C     =============  RUNNING MASSES 
      RMS = RUNM_HDEC(AMA,3)
      RMC = RUNM_HDEC(AMA,4)
      RMB = RUNM_HDEC(AMA,5)
      RMT = RUNM_HDEC(AMA,6)
      RATCOUP = GAT/GAB
      HIGTOP = AMA**2/AMT**2

      ASH=ALPHAS_HDEC(AMA,2)
      AMC0=1.D8
      AMB0=2.D8
C     AMT0=3.D8
      AS3=ALPHAS_HDEC(AMA,2)
      AMC0=AMC
      AS4=ALPHAS_HDEC(AMA,2)
      AMB0=AMB
C     AMT0=AMT

C     =============== PARTIAL WIDTHS
C  A ---> G G
       EPS=1.D-8
       NFEXT = 3
       ASG = AS3
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CAT = CTT*CF(CTT)*GAT
       CAB = CTB*CF(CTB)*GAB
       CTC = 4*AMC**2/AMA**2*DCMPLX(1D0,-EPS)
       CAC = CTC*CF(CTC)*GAT
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       HGG=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC

C  A ---> G G* ---> G CC   TO BE ADDED TO A ---> CC
       NFEXT = 4
       ASG = AS4
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       DCC=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC
     .     - HGG

C  A ---> G G* ---> G BB   TO BE ADDED TO A ---> BB
       NFEXT = 5
       ASG = ASH
       FQCD=AGGQCD(ASG,NFEXT)
       XFAC = CDABS(CAT+CAB+CAC)**2*FQCD
       DBB=GF/(16.D0*PI*DSQRT(2.D0))*AMA**3*(ASG/PI)**2*XFAC
     .     - HGG - DCC

      IF(NFGG.EQ.5)THEN
       HGG = HGG + DBB + DCC
       DBB = 0
       DCC = 0
      ELSEIF(NFGG.EQ.4)THEN
       HGG = HGG + DCC
       DCC = 0
      ENDIF

c      XFAC0 = CDABS(CAT+CAB)**2*FQCD
c      write(6,*)'XFAC (',AMA,') = ',xfac
c      write(6,*)'XFAC0(',AMA,') = ',xfac0,xfac0/xfac
c      write(6,*)'gab, gat: ',gab,gat
c      write(6,*)'BR(A -> gg) = ',HGG

C  A ---> MU MU
      IF(AMA.LE.2*AMMUON) THEN
       HMM = 0
      ELSE
      HMM=AFF(AMA,(AMMUON/AMA)**2)*GAB**2
      ENDIF
C  A ---> LL
      IF(AMA.LE.2*AMTAU) THEN
       HLL = 0
      ELSE
      HLL=AFF(AMA,(AMTAU/AMA)**2)*GAB**2
      ENDIF
C  A --> SS
      IF(AMA.LE.2*AMS) THEN
       HSS = 0
      ELSE
       HS1=3.D0*AFF(AMA,(AMS/AMA)**2)
     .    *GAB**2
     .    *TQCDA(AMS**2/AMA**2)
       HS2=3.D0*AFF(AMA,(RMS/AMA)**2)
     .    *GAB**2
     .    *QCDA(RMS**2/AMA**2)
       IF(HS2.LT.0.D0) HS2 = 0
       RAT = 2*AMS/AMA
       HSS = QQINT_HDEC(RAT,HS1,HS2)
      ENDIF
C  A --> CC
      RATCOUP = 1
      IF(AMA.LE.2*AMC) THEN
       HCC = 0
      ELSE
       HC1=3.D0*AFF(AMA,(AMC/AMA)**2)
     .    *GAT**2
     .    *TQCDA(AMC**2/AMA**2)
       HC2=3.D0*AFF(AMA,(RMC/AMA)**2)
     .    *GAT**2
     .    *QCDA(RMC**2/AMA**2)
     .   + DCC
       IF(HC2.LT.0.D0) HC2 = 0
       RAT = 2*AMC/AMA
       HCC = QQINT_HDEC(RAT,HC1,HC2)
      ENDIF
C  A --> BB :
      QQ = AMB
      SUSY = 0
      XGAB = GAB
c     SSUSY = AMA
      SSUSY = (AMSB(1)+AMSB(2)+AMGLU)/3*QSUSY
      AS0 = ALPHAS_HDEC(SSUSY,2)
      IF(IOFSUSY.EQ.0) THEN
       I0 = 1
       CALL DMBAPP_HDEC(I0,DGLB,DGHB,DGAB,QSUSY,LOOP)
       I0 = 3
       BSC = (AMSQ+AMUR+AMDR)/3
c      XMB = RUNM_HDEC(BSC,5)
       XMB = AMB
       SUSY = COFSUSY_HDEC(I0,AMB,XMB,QQ)*AS0/PI - 2*DGAB
       CALL BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,QSUSY,LOOP)
      ENDIF
      RATCOUP = GAT/XGAB
      IF(AMA.LE.2*AMB) THEN
       HBB = 0
      ELSE
       HB1=3.D0*AFF(AMA,(AMB/AMA)**2)
     .    *(XGAB**2+XGAB*GAB*SUSY)
     .    *TQCDA(AMB**2/AMA**2)
       HB2=3.D0*AFF(AMA,(RMB/AMA)**2)
     .    *(XGAB**2+XGAB*GAB*SUSY)
     .    *QCDA(RMB**2/AMA**2)
     .   + DBB
       IF(HB2.LT.0.D0) HB2 = 0
       RAT = 2*AMB/AMA
       HBB = QQINT_HDEC(RAT,HB1,HB2)

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      XB0=3.D0*AFF(AMA,(AMB/AMA)**2)
c    .    *GAB**2
c      XB1=3.D0*AFF(AMA,(RMB/AMA)**2)
c    .    *GAB**2
c    .    *QCDH(RMB**2/AMA**2)
c    .   + DBB
c      XB2=3.D0*AFF(AMA,(RMB/AMA)**2)
c    .    *(XGAB**2+XGAB*GAB*SUSY)
c    .    *QCDH(RMB**2/AMA**2)
c    .   + DBB
c      write(53,('5(1X,G15.8)'))AMA,AMA,XB0,XB1,XB2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      write(53,('4(1X,G15.8)'))AMA,AMA,HBB,SUSY/(SUSY+2*DGAB)
c      write(53,('4(1X,G15.8)'))AMA,AMA,SUSY+2*DGAB,SUSY/(SUSY+2*DGAB)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c      X1 = (QCDA(RMB**2/AMA**2)*AFF(AMA,(RMB/AMA)**2)/
c    .       AFF(AMA,(AMB/AMA)**2)-1)
c      X2 = (SUSY-1)

c     RATCOUP = GAT/XGAB
c      HB1X=3.D0*AFF(AMA,(AMB/AMA)**2)
c    .    *XGAB**2
c    .    *TQCDA(AMB**2/AMA**2)
c    .    /(BETA_HDEC(AMB**2/AMA**2))
c    .    *SUSY
c      HB2X=3.D0*AFF(AMA,(RMB/AMA)**2)*XGAB**2
c    .    *QCDA(RMB**2/AMA**2)
c    .    /(BETA_HDEC(RMB**2/AMA**2))
c    .    *SUSY

c     RATCOUP = 0
c     deltaqcd = QCDA(RMB**2/AMA**2)
c     RATCOUP = GAT/XGAB
c     deltat = QCDA(RMB**2/AMA**2) - deltaqcd

c      write(6,*)
c      write(6,*)'A:'
c      write(6,*)'MB,RUNMB,alpha_s: ',AMB,RMB,ASH
c      write(6,*)'MA =              ',AMA
c      write(6,*)'QCD           SUSY        APPROX     APPROX/FULL',
c    .           ' GbA(QCD) GbA(SQCD):'
c      write(6,*)X1,X2+2*DGAB,2*DGAB,2*DGAB/(X2+2*DGAB),GAB,XGAB
c      write(6,*)'Resummation: ',(XGAB/GAB)**2-1
c      write(6,*)'Rest:        ',SUSY-1
c      write(6,*)'Total SUSY:  ',(XGAB/GAB)**2*SUSY-1
c      write(6,*)'deltaqcd,t = ',deltaqcd,deltat
c      write(6,*)'Gamma(0)   = ',HB2X,HB1X
c      write(6,*)'Gamma(mb)  = ',HB2,HB1
c      write(6,*)'Rest: A      ',AMA,AMA,(SUSY-1)/(X2+2*DGAB)
c      write(53,*)AMA,AMA,(SUSY-1)/(X2+2*DGAB)
c      write(6,*)
      ENDIF
C  A --> TT :
      RATCOUP = 0
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=4D0
       XM1 = 2D0*AMT-DLD
       XM2 = 2D0*AMT+DLU
       IF (AMA.LE.AMT+AMW+AMB) THEN
        HTT=0.D0
       ELSEIF (AMA.LE.XM1) THEN
        FACTT=6.D0*GF**2*AMA**3*AMT**2/2.D0/128.D0/PI**3*GAT**2
        CALL ATOTT_HDEC(AMA,AMT,AMB,AMW,AMCH,ATT0)
        HTT=FACTT*ATT0
       ELSEIF (AMA.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        FACTT=6.D0*GF**2*XX(1)**3*AMT**2/2.D0/128.D0/PI**3
        CALL ATOTT_HDEC(XX(1),AMT,AMB,AMW,AMCH,ATT0)
        YY(1)=FACTT*ATT0
        FACTT=6.D0*GF**2*XX(2)**3*AMT**2/2.D0/128.D0/PI**3
        CALL ATOTT_HDEC(XX(2),AMT,AMB,AMW,AMCH,ATT0)
        YY(2)=FACTT*ATT0
        XMT = RUNM_HDEC(XX(3),6)
        XYZ1 =3.D0*AFF(XX(3),(AMT/XX(3))**2)
     .    *TQCDA(AMT**2/XX(3)**2)
        XYZ2 =3.D0*AFF(XX(3),(XMT/XX(3))**2)
     .    *QCDA(XMT**2/XX(3)**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        RAT = 2*AMT/XX(3)
        YY(3) = QQINT_HDEC(RAT,XYZ1,XYZ2)
        XMT = RUNM_HDEC(XX(4),6)
        XYZ1 =3.D0*AFF(XX(4),(AMT/XX(4))**2)
     .    *TQCDA(AMT**2/XX(4)**2)
        XYZ2 =3.D0*AFF(XX(4),(XMT/XX(4))**2)
     .    *QCDA(XMT**2/XX(4)**2)
        IF(XYZ2.LT.0.D0) XYZ2 = 0
        RAT = 2*AMT/XX(4)
        YY(4) = QQINT_HDEC(RAT,XYZ1,XYZ2)
        HTT = FINT_HDEC(AMA,XX,YY)*GAT**2
       ELSE
        HT1=3.D0*AFF(AMA,(AMT/AMA)**2)*GAT**2
     .    *TQCDA(AMT**2/AMA**2)
        HT2=3.D0*AFF(AMA,(RMT/AMA)**2)*GAT**2
     .    *QCDA(RMT**2/AMA**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMA
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ELSE
       IF (AMA.LE.2.D0*AMT) THEN
        HTT=0.D0
       ELSE
        HT1=3.D0*AFF(AMA,(AMT/AMA)**2)*GAT**2
     .    *TQCDA(AMT**2/AMA**2)
        HT2=3.D0*AFF(AMA,(RMT/AMA)**2)*GAT**2
     .    *QCDA(RMT**2/AMA**2)
        IF(HT2.LT.0.D0) HT2 = 0
        RAT = 2*AMT/AMA
        HTT = QQINT_HDEC(RAT,HT1,HT2)
       ENDIF
      ENDIF
C  A ---> GAMMA GAMMA
       EPS=1.D-8
       XRMC = RUNM_HDEC(AMA/2,4)*AMC/RUNM_HDEC(AMC,4)
       XRMB = RUNM_HDEC(AMA/2,5)*AMB/RUNM_HDEC(AMB,5)
       XRMT = RUNM_HDEC(AMA/2,6)*AMT/RUNM_HDEC(AMT,6)
       CTT = 4*XRMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*XRMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CAT = 4/3D0 * CTT*CF(CTT)*GAT
     .     * CFACQ_HDEC(1,AMA,XRMT)
       CAB = 1/3D0 * CTB*CF(CTB)*GAB
     .     * CFACQ_HDEC(1,AMA,XRMB)
       CTC = 4*XRMC**2/AMA**2*DCMPLX(1D0,-EPS)
       CAC = 4/3D0 * CTC*CF(CTC)*GAT
     .     * CFACQ_HDEC(1,AMA,XRMC)
       CTL = 4*AMTAU**2/AMA**2*DCMPLX(1D0,-EPS)
       CAL = 1.D0  * CTL*CF(CTL)*GAB
       IF(IOFSUSY.EQ.0) THEN 
        CX1 = 4*AMCHAR(1)**2/AMA**2*DCMPLX(1D0,-EPS)
        CX2 = 4*AMCHAR(2)**2/AMA**2*DCMPLX(1D0,-EPS)
        CAX1= AMW/XMCHAR(1) * CX1*CF(CX1) * 2*AC3(1,1) 
        CAX2= AMW/XMCHAR(2) * CX2*CF(CX2) * 2*AC3(2,2) 
        XFAC = CDABS(CAT+CAB+CAC+CAL+CAX1+CAX2)**2
       ELSE 
        XFAC = CDABS(CAT+CAB+CAC+CAL)**2
       ENDIF
       HGA=GF/(32.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)**2*XFAC
C  A ---> Z GAMMA
      IF(AMA.LE.AMZ)THEN
       HZGA=0
      ELSE
       TS = SS/CS
       FT = -3*2D0/3*(1-4*2D0/3*SS)/DSQRT(SS*CS)*GAT
       FB = 3*1D0/3*(-1+4*1D0/3*SS)/DSQRT(SS*CS)*GAB
       EPS=1.D-8
       CTT = 4*AMT**2/AMA**2*DCMPLX(1D0,-EPS)
       CTB = 4*AMB**2/AMA**2*DCMPLX(1D0,-EPS)
       CLT = 4*AMT**2/AMZ**2*DCMPLX(1D0,-EPS)
       CLB = 4*AMB**2/AMZ**2*DCMPLX(1D0,-EPS)
       CAT = FT*(- CI2(CTT,CLT))
       CAB = FB*(- CI2(CTB,CLB))
       XFAC = CDABS(CAT+CAB)**2
       ACOUP = DSQRT(2D0)*GF*AMZ**2*SS*CS/PI**2
       HZGA = GF/(4.D0*PI*DSQRT(2.D0))*AMA**3*(ALPH/PI)*ACOUP/16.D0
     .        *XFAC*(1-AMZ**2/AMA**2)**3
      ENDIF
C  A ---> H Z* ---> HFF
      IF(IONSH.EQ.0)THEN
       DLD=3D0
       DLU=5D0
       XM1 = AML+AMZ-DLD
       XM2 = AML+AMZ+DLU
       IF (AMA.LE.AML) THEN
        HAZ=0
       ELSEIF (AMA.LE.XM1) THEN
        IF (AMA.LE.DABS(AMZ-AML)) THEN
         HAZ=0
        ELSE
         HAZ=9.D0*GF**2/8.D0/PI**3*AMZ**4*AMA*GZAL**2*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/AMA)**2,(AMZ/AMA)**2)
        ENDIF
       ELSEIF (AMA.LE.XM2) THEN
        XX(1) = XM1-1D0
        XX(2) = XM1
        XX(3) = XM2
        XX(4) = XM2+1D0
        YY(1)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(1)*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/XX(1))**2,(AMZ/XX(1))**2)
        YY(2)=9.D0*GF**2/8.D0/PI**3*AMZ**4*XX(2)*
     .      (7.D0/12.D0-10.D0/9.D0*SS+40.D0/27.D0*SS**2)
     .      *HVH((AML/XX(2))**2,(AMZ/XX(2))**2)
        CAZ=LAMB_HDEC(AML**2/XX(3)**2,AMZ**2/XX(3)**2)
     .     *LAMB_HDEC(XX(3)**2/AMZ**2,AML**2/AMZ**2)**2
        YY(3)=GF/8D0/DSQRT(2D0)/PI*AMZ**4/XX(3)*CAZ
        CAZ=LAMB_HDEC(AML**2/XX(4)**2,AMZ**2/XX(4)**2)
     .     *LAMB_HDEC(XX(4)**2/AMZ**2,AML**2/AMZ**2)**2
        YY(4)=GF/8D0/DSQRT(2D0)/PI*AMZ**4/XX(4)*CAZ
        HAZ = FINT_HDEC(AMA,XX,YY)*GZAL**2
       ELSE
        CAZ=LAMB_HDEC(AML**2/AMA**2,AMZ**2/AMA**2)
     .     *LAMB_HDEC(AMA**2/AMZ**2,AML**2/AMZ**2)**2
        HAZ=GF/8D0/DSQRT(2D0)/PI*AMZ**4/AMA*GZAL**2*CAZ
       ENDIF
      ELSE
       IF (AMA.LE.AMZ+AML) THEN
        HAZ=0
       ELSE
        CAZ=LAMB_HDEC(AML**2/AMA**2,AMZ**2/AMA**2)
     .     *LAMB_HDEC(AMA**2/AMZ**2,AML**2/AMZ**2)**2
        HAZ=GF/8D0/DSQRT(2D0)/PI*AMZ**4/AMA*GZAL**2*CAZ
       ENDIF
      ENDIF
C
C ========================== SUSY DECAYS  
C
      IF(IOFSUSY.EQ.0) THEN 
C  A ----> CHARGINOS
      DO 731 I=1,2
      DO 731 J=1,2
      IF (AMA.GT.AMCHAR(I)+AMCHAR(J)) THEN
      WHACH(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMA
     .     *LAMB_HDEC(AMCHAR(I)**2/AMA**2,AMCHAR(J)**2/AMA**2)
     .     *( (AC3(I,J)**2+AC3(J,I)**2)*(AMA**2-AMCHAR(I)
     .         **2-AMCHAR(J)**2)+4.D0*AC3(I,J)*AC3(J,I)* 
     .         XMCHAR(I)*XMCHAR(J) ) 
      ELSE 
      WHACH(I,J)=0.D0
      ENDIF
 731  CONTINUE
      WHACHT=WHACH(1,1)+WHACH(1,2)+WHACH(2,1)+WHACH(2,2)
C  A ----> NEUTRALINOS 
      DO 732 I=1,4
      DO 732 J=1,4
      IF (AMA.GT.AMNEUT(I)+AMNEUT(J)) THEN
      WHANE(I,J)=GF*AMW**2/(2*PI*DSQRT(2.D0))/AMA
     .         *AN3(I,J)**2*(AMA**2-(XMNEUT(I)-XMNEUT(J))**2)
     .         *LAMB_HDEC(AMNEUT(I)**2/AMA**2,AMNEUT(J)**2/AMA**2)
      ELSE 
      WHANE(I,J)=0.D0
      ENDIF
 732  CONTINUE
      WHANET= WHANE(1,1)+WHANE(1,2)+WHANE(1,3)+WHANE(1,4)
     .       +WHANE(2,1)+WHANE(2,2)+WHANE(2,3)+WHANE(2,4)
     .       +WHANE(3,1)+WHANE(3,2)+WHANE(3,3)+WHANE(3,4)
     .       +WHANE(4,1)+WHANE(4,2)+WHANE(4,3)+WHANE(4,4)

C  A ----> STAU'S 
C
      IF(AMA.GT.AMSL(1)+AMSL(2)) THEN
      WHASL=GF*AMZ**4/DSQRT(2.D0)/PI*GAEE**2*
     .      LAMB_HDEC(AMSL(1)**2/AMA**2,AMSL(2)**2/AMA**2)/AMA
      ELSE
      WHASL=0.D0
      ENDIF
C
C  A ----> STOPS 
C
      SUSY = 1
      IF(AMA.GT.AMST(1)+AMST(2)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(3,1,1,2,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHAST=3*GF*AMZ**4/DSQRT(2.D0)/PI*GATT**2*
     .       LAMB_HDEC(AMST(1)**2/AMA**2,AMST(2)**2/AMA**2)/AMA
     .      *SUSY
c      write(6,*)'A -> stop: ',AMA,AMST(1),AMST(2),SUSY-1,
c    .           WHAST/2/SUSY,WHAST/2
c      write(6,*)'A -> stop: ',AMA,AMST(1),AMST(2),SUSY-1
      ELSE
      WHAST=0.D0
      ENDIF
C
C  A ----> SBOTTOMS 
C
      SUSY = 1
      IF(AMA.GT.AMSB(1)+AMSB(2)) THEN
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      SUSY = 1+SQSUSY_HDEC(3,2,1,2,QSQ)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       WHASB=3*GF*AMZ**4/DSQRT(2.D0)/PI*GABB**2*
     .       LAMB_HDEC(AMSB(1)**2/AMA**2,AMSB(2)**2/AMA**2)/AMA
     .      *SUSY
c      write(6,*)'A -> sbot: ',AMA,AMSB(1),AMSB(2),SUSY-1,
c    .           WHASB/2/SUSY,WHASB/2
c      write(6,*)'A -> sbot: ',AMA,AMSB(1),AMSB(2),SUSY-1
      ELSE
      WHASB=0.D0
      ENDIF
C
      ELSE 
      WHACHT=0.D0
      WHANET=0.D0
      WHASL=0.D0
      WHAST=0.D0
      WHASB=0.D0
C--Change thanks to Elzbieta Richter-Was
      DO I=1,2
       DO J=1,2
        WHACH(I,J)=0.D0
       ENDDO
      ENDDO
      DO I=1,4
       DO J=1,4
        WHANE(I,J)=0.D0
       ENDDO
      ENDDO
      ENDIF

      IF(IGOLD.NE.0)THEN
C   HA ---> GOLDSTINOS
       DO 730 I=1,4
       IF (AMA.GT.AMNEUT(I)) THEN
        WHAGD(I)=AMA**5/AXMPL**2/AXMGD**2/48.D0/PI*
     .           (1.D0-AMNEUT(I)**2/AMA**2)**4*AGDA(I)**2
       ELSE
        WHAGD(I)=0.D0
       ENDIF
 730   CONTINUE
       WHAGDT=WHAGD(1)+WHAGD(2)+WHAGD(3)+WHAGD(4)
      ELSE
       WHAGDT=0
      ENDIF
C
C    ==========  TOTAL WIDTH AND BRANCHING RATIOS 
      WTOT=HLL+HMM+HSS+HCC+HBB+HGG+HGA+HZGA+HAZ+HTT
     .    +WHACHT+WHANET+WHASL+WHAST+WHASB + WHAGDT
      ABRT=HTT/WTOT
      ABRB=HBB/WTOT
      ABRL=HLL/WTOT
      ABRM=HMM/WTOT
      ABRS=HSS/WTOT
      ABRC=HCC/WTOT
      ABRG=HGG/WTOT
      ABRGA=HGA/WTOT
      ABRZGA=HZGA/WTOT
      ABRZ=HAZ/WTOT
      DO 831 I=1,2
      DO 831 J=1,2
      HABRSC(I,J)=WHACH(I,J)/WTOT
831   CONTINUE
      DO 832 I=1,4
      DO 832 J=1,4
      HABRSN(I,J)=WHANE(I,J)/WTOT
832   CONTINUE
      HABRCHT=WHACHT/WTOT      
      HABRNET=WHANET/WTOT      
      HABRSL=WHASL/WTOT 
      HABRST=WHAST/WTOT 
      HABRSB=WHASB/WTOT 
      HABRGD=WHAGDT/WTOT
C 
      AWDTH=WTOT

      BHASTAU = WHASL/WTOT
      BHASB = WHASB/WTOT
      BHAST = WHAST/WTOT

C    ==============================================================
      ENDIF

      RETURN
      END
 
      DOUBLE PRECISION FUNCTION BIJ_HDEC(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB_HDEC
      BIJ_HDEC = (1-X-Y)/LAMB_HDEC(X,Y)*(
     .          4*SP_HDEC(XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .        - 2*SP_HDEC(-XI_HDEC(X,Y)) - 2*SP_HDEC(-XI_HDEC(Y,X))
     .        + 2*DLOG(XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .           *DLOG(1-XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .        - DLOG(XI_HDEC(X,Y))*DLOG(1+XI_HDEC(X,Y))
     .        - DLOG(XI_HDEC(Y,X))*DLOG(1+XI_HDEC(Y,X))
     .          )
     .        -4*(DLOG(1-XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .        +XI_HDEC(X,Y)*XI_HDEC(Y,X)/(1-XI_HDEC(X,Y)*XI_HDEC(Y,X))
     .          *DLOG(XI_HDEC(X,Y)*XI_HDEC(Y,X)))
     .        +(LAMB_HDEC(X,Y)+X-Y)/LAMB_HDEC(X,Y)*(DLOG(1+XI_HDEC(X,Y))
     .              - XI_HDEC(X,Y)/(1+XI_HDEC(X,Y))*DLOG(XI_HDEC(X,Y)))
     .        +(LAMB_HDEC(X,Y)-X+Y)/LAMB_HDEC(X,Y)*(DLOG(1+XI_HDEC(Y,X))
     .              - XI_HDEC(Y,X)/(1+XI_HDEC(Y,X))*DLOG(XI_HDEC(Y,X)))
      RETURN
      END

      DOUBLE PRECISION FUNCTION BETA_HDEC(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      BETA_HDEC=DSQRT(1.D0-4.D0*X)
      RETURN
      END

      DOUBLE PRECISION FUNCTION LAMB_HDEC(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LAMB_HDEC=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
      RETURN
      END

      DOUBLE PRECISION FUNCTION XI_HDEC(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB_HDEC
      XI_HDEC = 2*X/(1-X-Y+LAMB_HDEC(X,Y))
      RETURN
      END

C *****************************************************************
C ************* SUBROUTINE FOR THE SUSY COUPLINGS *****************
C *****************************************************************
      SUBROUTINE SUSYCP_HDEC(TGBET)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION LA1,LA2,LA3,LA4,LA5,LA6,LA7,LA3T
      COMPLEX*16 F0_HDEC
      DIMENSION MST(2),GLTT(2,2),GHTT(2,2),
     .          MSB(2),GLBB(2,2),GHBB(2,2)
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
      COMMON/MODEL_HDEC/IMODEL
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      COMMON/HMASSR_HDEC/AMLR,AMHR
      COMMON/CHIMASS_HDEC/AMCHI
      COMMON/HSELF_HDEC/LA1,LA2,LA3,LA4,LA5,LA6,LA7
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/BREAKGLU_HDEC/AMGLU
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/GLUINO_HDEC/AMGLUINO,AMSB1,AMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/SLHA_hmass_HDEC/slhaml,slhamh,slhamc,slha_alpha
      COMMON/SLHA_gaug_HDEC/slhaneut(4),slhaxneut(4),slhachar(2),
     .          slhau(2,2),slhav(2,2),slhaz(4,4),slhaxchar(2)

      PI=4*DATAN(1D0)
      V=1.D0/DSQRT(DSQRT(2.D0)*GF)
      BET=DATAN(TGBET)
      SB = DSIN(BET)
      CB = DCOS(BET)
      AMAR = AMA
C  ============ HEAVIEST CHARGINO MASS NEEDED FOR SUBH ========== 
      if(islhai.eq.0) then
         AMCHI2=AM2**2+AMU**2+2.D0*AMW**2+DSQRT((AM2**2-AMU**2)**2
     .        +4.D0*AMW**4*DCOS(2.D0*BET)**2+4.D0*AMW**2*
     .        (AM2**2+AMU**2+2.D0*AMU*AM2*DSIN(2.D0*BET) ) ) 
         AMCHI=DSQRT(0.5D0*AMCHI2)
      else
         amchi = slhachar(2)
      endif
C ===============================================================
C ========== RUNNING MASSES
      if(islhai.eq.0) then
      IF(IMODEL.EQ.1)THEN
       CALL SUBH1_HDEC(AMA,TGBET,AMSQ,AMUR,AMDR,AMT,AU,AD,AMU,AMCHI,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA,AMGLU)
      ELSEIF(IMODEL.EQ.2)THEN
       CALL SUBH2_HDEC(AMA,TGBET,AMSQ,AMUR,AMT,AU,AD,AMU,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA)
c change susyhit
c      ELSEIF(IMODEL.EQ.3)THEN
c       CALL HABER(TGBET,SA,CA)
c       AMLR = AML
c       AMHR = AMH
c end change susyhit
      ELSEIF(IMODEL.EQ.4)THEN
C--Use Carena et al. for everything not included in FeynHiggs....
c      CALL SUBH2_HDEC(AMA,TGBET,AMSQ,AMUR,AMT,AU,AD,AMU,
c    .            AMLR,AMHR,AMCH,SA,CA,TANBA)
       CALL SUBH1_HDEC(AMA,TGBET,AMSQ,AMUR,AMDR,AMT,AU,AD,AMU,AMCHI,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA,AMGLU)
       IF(CTHT.GE.0.D0)THEN
        XMST1 = AMST1
        XMST2 = AMST2
        STT = STHT
       ELSE
        XMST1 = AMST1
        XMST2 = AMST2
        STT = CTHT
       ENDIF
       IF(CTHB.GE.0.D0)THEN
        XMSB1 = AMSB1
        XMSB2 = AMSB2
        STB = STHB
       ELSE
        XMSB1 = AMSB1
        XMSB2 = AMSB2
        STB = CTHB
       ENDIF
c change susyhit
c       CALL FEYNHIGGS(AMA,TGBET,AMT,XMST1,XMST2,STT,XMSB1,
c     .                XMSB2,STB,AMU,AMGLU,AM2,AMLR,AMHR,SA,CA)
c end change susyhit
      ENDIF
      else
       CALL SUBH1_HDEC(AMA,TGBET,AMSQ,AMUR,AMDR,AMT,AU,AD,AMU,AMCHI,
     .            AMLR,AMHR,AMCH,SA,CA,TANBA,AMGLU)
       amlr = slhaml
       amhr = slhamh
       aml  = slhaml
       amh  = slhamh
       amch = slhamc
       sa   = dsin(slha_alpha)
       ca   = dcos(slha_alpha)
      endif
      LA3T=LA3+LA4+LA5
      AMA2=AMAR**2
      AML2=AMLR**2
      AMH2=AMHR**2
      AMP2=AMCH**2
C ========== HIGGS COUPLINGS 
      SBMA = SB*CA-CB*SA
      CBMA = CB*CA+SB*SA
      SBPA = SB*CA+CB*SA
      CBPA = CB*CA-SB*SA
      S2A = 2*SA*CA
      C2A = CA**2-SA**2
      S2B = 2*SB*CB
      C2B = CB**2-SB**2
      GLZZ = 1/V/2*AML2*SBMA
      GHZZ = 1/V/2*AMH2*CBMA
      GLWW = 2*GLZZ
      GHWW = 2*GHZZ
      GLAZ = 1/V*(AML2-AMA2)*CBMA
      GHAZ = -1/V*(AMH2-AMA2)*SBMA
      GLPW = -1/V*(AMP2-AML2)*CBMA
      GLMW = GLPW
      GHPW = 1/V*(AMP2-AMH2)*SBMA
      GHMW = GHPW
      GAPW = 1/V*(AMP2-AMA2)
      GAMW = -GAPW
      GHHH = V/2*(LA1*CA**3*CB + LA2*SA**3*SB + LA3T*SA*CA*SBPA
     .     + LA6*CA**2*(3*SA*CB+CA*SB) + LA7*SA**2*(3*CA*SB+SA*CB))
      GLLL = -V/2*(LA1*SA**3*CB - LA2*CA**3*SB + LA3T*SA*CA*CBPA
     .     - LA6*SA**2*(3*CA*CB-SA*SB) + LA7*CA**2*(3*SA*SB-CA*CB))
      GLHH = -3*V/2*(LA1*CA**2*CB*SA - LA2*SA**2*SB*CA
     .     + LA3T*(SA**3*CB-CA**3*SB+2*SBMA/3)
     .     - LA6*CA*(CB*C2A-SA*SBPA) - LA7*SA*(C2A*SB+CA*SBPA))
      GHLL = 3*V/2*(LA1*SA**2*CB*CA + LA2*CA**2*SB*SA
     .     + LA3T*(SA**3*SB+CA**3*CB-2*CBMA/3)
     .     - LA6*SA*(CB*C2A+CA*CBPA) + LA7*CA*(C2A*SB+SA*CBPA))
      GLAA = -V/2*(LA1*SB**2*CB*SA - LA2*CB**2*SB*CA
     .     - LA3T*(SB**3*CA-CB**3*SA) + 2*LA5*SBMA
     .     - LA6*SB*(CB*SBPA+SA*C2B) - LA7*CB*(C2B*CA-SB*SBPA))
      GHAA = V/2*(LA1*SB**2*CB*CA + LA2*CB**2*SB*SA
     .     + LA3T*(SB**3*SA+CB**3*CA) - 2*LA5*CBMA
     .     - LA6*SB*(CB*CBPA+CA*C2B) + LA7*CB*(SB*CBPA+SA*C2B))
      GLPM = 2*GLAA + V*(LA5 - LA4)*SBMA
      GHPM = 2*GHAA + V*(LA5 - LA4)*CBMA
      GLZZ = 2*GLZZ
      GHZZ = 2*GHZZ
      GLLL = 6*GLLL
      GHHH = 6*GHHH
      GLHH = 2*GLHH
      GHLL = 2*GHLL
      GLAA = 2*GLAA
      GHAA = 2*GHAA
      XNORM = AMZ**2/V
      GLLL = GLLL/XNORM
      GHLL = GHLL/XNORM
      GLHH = GLHH/XNORM
      GHHH = GHHH/XNORM
      GHAA = GHAA/XNORM
      GLAA = GLAA/XNORM
      GLPM = GLPM/XNORM
      GHPM = GHPM/XNORM
      GAT=1.D0/TGBET
      GAB=TGBET
      GLT=CA/SB
      GLB=-SA/CB
      GHT=SA/SB
      GHB=CA/CB
      GZAL=-CBMA
      GZAH=SBMA
      GLVV=SBMA
      GHVV=CBMA
      B=BET
      IF(CA.EQ.0)THEN
       A = PI/2
      ELSE
       A=DATAN(SA/CA)
      ENDIF
      IF(CA.LT.0D0)THEN
       IF(SA.LT.0D0)THEN
        A = A-PI
       ELSE
        A = A+PI
       ENDIF
      ENDIF
C ===============================================================
C ========== POLE MASSES 
      if(islhai.eq.0) then
      IF(IMODEL.EQ.1)THEN
      IF(IPOLE.EQ.1) THEN 
       MT=RUNM_HDEC(AMT,6)
       MB=RUNM_HDEC(AMT,5)
       SW2=1.D0-AMW**2/AMZ**2
C===== STOP MASSES
       MSTL2=AMSQ**2+(0.5D0-2.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
       MSTR2=AMUR**2+2.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
       MLRT=AU-AMU/TGBET
       DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
       MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
       MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)GOTO 111
       MST(1)=DSQRT(MST12)
       MST(2)=DSQRT(MST22)
       IF(MSTL2.EQ.MSTR2) THEN
        THET = PI/4
       ELSE
        THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
        IF(MSTL2.GT.MSTR2) THET = THET + PI/2
       ENDIF
       CST= DCOS(THET)
       SST= DSIN(THET)
C===== SBOTTOM MASSES
       MSBL2=AMSQ**2+(-0.5D0+1.D0/3.D0*SW2)*AMZ**2*DCOS(2.D0*B)
       MSBR2=AMDR**2-1.D0/3.D0*SW2*AMZ**2*DCOS(2.D0*B)
       MLRB=AD-AMU*TGBET
       DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
       MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
       MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)GOTO 111
       MSB(1)=DSQRT(MSB12)
       MSB(2)=DSQRT(MSB22)
       IF(MSBL2.EQ.MSBR2) THEN
        THEB = PI/4
       ELSE
        THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
        IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
       ENDIF
       CSB= DCOS(THEB)
       SSB= DSIN(THEB)
C===== LIGHT HIGGS COUPLINGS 
       GLTT(1,1)=-SBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET) )
     .     +MT**2/AMZ**2*GLT + MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(2,2)=-SBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET) )
     .     +MT**2/AMZ**2*GLT - MT*SST*CST/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(1,2)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
       GLTT(2,1)=-2*SBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GLT+AMU*GHT)
       GLBB(1,1)=-SBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GLB + MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(2,2)=-SBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GLB - MB*SSB*CSB/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(1,2)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
       GLBB(2,1)=-2*SBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GLB-AMU*GHB)
C===== HEAVY HIGGS COUPLINGS 
       GHTT(1,1)=CBPA*(0.5D0*CST**2-2.D0/3.D0*SW2*DCOS(2*THET))
     .     +MT**2/AMZ**2*GHT + MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(2,2)=CBPA*(0.5D0*SST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .     +MT**2/AMZ**2*GHT - MT*SST*CST/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(1,2)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     +MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
       GHTT(2,1)=2*CBPA*SST*CST*(2.D0/3.D0*SW2-0.25D0)
     .     + MT*DCOS(2*THET)/2.D0/AMZ**2*(AU*GHT-AMU*GLT)
       GHBB(1,1)=CBPA*(-0.5D0*CSB**2+1.D0/3.D0*SW2*DCOS(2*THEB))
     .     +MB**2/AMZ**2*GHB + MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(2,2)=CBPA*(-0.5D0*SSB**2-1.D0/3.D0*SW2*DCOS(2*THEB))
     .     + MB**2/AMZ**2*GHB - MB*SSB*CSB/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(1,2)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
       GHBB(2,1)=2*CBPA*SSB*CSB*(-1.D0/3.D0*SW2+0.25D0)
     .     + MB*DCOS(2*THEB)/2.D0/AMZ**2*(AD*GHB+AMU*GLB)
C===== PSEUDOSCALAR HIGGS COUPLINGS 
       GATT=MT/2.D0/AMZ**2*(AMU+AU*GAT) 
       GABB=MB/2.D0/AMZ**2*(AMU+AD*GAB) 
C======= LOOP CORRECTIONS  
       XDLT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLT**2*(-2.D0*MT**2+0.5D0*AML2)
     .     *DREAL(F0_HDEC(MT,MT,AML2))
     .     *3*MT**2
       XDLB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(-2.D0*MB**2+0.5D0*AML2)
     .     *DREAL(F0_HDEC(MB,MB,AML2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GLB**2*(0.5D0*AML2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDHT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHT**2*(-2.D0*MT**2+0.5D0*AMH2)
     .     *DREAL(F0_HDEC(MT,MT,AMH2))
     .     *3*MT**2
       XDHB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(-2.D0*MB**2+0.5D0*AMH2)
     .     *DREAL(F0_HDEC(MB,MB,AMH2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GHB**2*(0.5D0*AMH2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDAT=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAT**2*(-0.5D0*AMA2)
     .     *DREAL(F0_HDEC(MT,MT,AMA2))
     .     *3*MT**2
       XDAB=GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .     *DREAL(F0_HDEC(MB,MB,AMA2))
     .     *3*MB**2
C--BUG IN CARENA ET AL. FIXED
     .     +GF/(2.D0*DSQRT(2.D0)*PI**2)*GAB**2*(-0.5D0*AMA2)
     .     *DLOG(MB**2/MT**2)
     .     *3*MB**2
       XDLST=0.D0
       XDLSB=0.D0
       XDHST=0.D0
       XDHSB=0.D0
         DO 311 I=1,2
         DO 311 J=1,2
       XDLST=XDLST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLTT(I,J)**2*
     .       DREAL(F0_HDEC(MST(I),MST(J),AML2))
     .     *3*AMZ**4
       XDLSB=XDLSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GLBB(I,J)**2*
     .       DREAL(F0_HDEC(MSB(I),MSB(J),AML2))
     .    *3*AMZ**4
       XDHST=XDHST+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHTT(I,J)**2*
     .       DREAL(F0_HDEC(MST(I),MST(J),AMH2))
     .     *3*AMZ**4
       XDHSB=XDHSB+GF/(2.D0*DSQRT(2.D0)*PI**2)*GHBB(I,J)**2*
     .       DREAL(F0_HDEC(MSB(I),MSB(J),AMH2))
     .     *3*AMZ**4
311    CONTINUE
       XDAST=GF/(1.D0*DSQRT(2.D0)*PI**2)*GATT**2*
     .       DREAL(F0_HDEC(MST(1),MST(2),AMA2))
     .     *3*AMZ**4
       XDASB=GF/(1.D0*DSQRT(2.D0)*PI**2)*GABB**2*
     .       DREAL(F0_HDEC(MSB(1),MSB(2),AMA2))
     .     *3*AMZ**4
      
       AML=DSQRT(AML2+XDLT+XDLB+XDLST+XDLSB)
       AMH=DSQRT(AMH2+XDHT+XDHB+XDHST+XDHSB)  
       AMA=DSQRT(AMA2+XDAT+XDAB+XDAST+XDASB)  
      ELSE
       AML=AMLR
       AMH=AMHR     
       AMA=AMAR     
      ENDIF 
      ELSE
       AML=AMLR
       AMH=AMHR
       AMA=AMAR
      ENDIF
      endif
      RETURN
111   STOP
      END

C ===================== THE FUNCTION F0 ===============
      COMPLEX*16 FUNCTION F0_HDEC(M1,M2,QSQ)
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
       F0_HDEC = 0.D0
      ELSE
       IF(M1.EQ.M2) THEN
        F0_HDEC = -2.D0 + CR*CDLOG(-(1+CR)/(1-CR))
       ELSE
        CBET = CDSQRT(1-4*M1*M2/(CQ2 - (M1-M2)**2))
        CXX = (CBET-1)/(CBET+1)
        F0_HDEC = -1 + ((QSQ+M2SQ-M1SQ)/2/QSQ - M2SQ/(M2SQ-M1SQ))
     .                                           *DLOG(M2SQ/M1SQ)
     .     - (QSQ-(M1-M2)**2)/QSQ*CBET*CDLOG(CXX)
       ENDIF
      ENDIF
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR HSM ---> V*V* ---> 4F
C     ************************************************************
      SUBROUTINE HTOVV_HDEC(AMH,AMV,GAMV,HTVV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/VVOFF_HDEC/AMH1,AMV1,GAMV1
      COMMON/PREC_HDEC/IP
      EXTERNAL FTOVV1_HDEC
      IP=20
      AMH1=AMH
      AMV1=AMV
      GAMV1=GAMV
      DLT=1D0/IP
      SUM=0D0
      DO 1 I=1,IP
       UU=DLT*I
       DD=UU-DLT
       CALL QGAUS1_HDEC(FTOVV1_HDEC,DD,UU,RES)
       SUM=SUM+RES
1     CONTINUE
      HTVV=SUM
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV1_HDEC(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FIRST_HDEC/X1
      COMMON/PREC_HDEC/IP
      EXTERNAL FTOVV2_HDEC
      X1=XX
      DLT=1D0/IP
      SUM=0D0
      DO 1 I=1,IP
       UU=DLT*I
       DD=UU-DLT
       CALL QGAUS2_HDEC(FTOVV2_HDEC,DD,UU,RES)
       SUM=SUM+RES
1     CONTINUE
      FTOVV1_HDEC=SUM
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV2_HDEC(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION YY(2)
      COMMON/FIRST_HDEC/X1
      YY(1)=X1
      YY(2)=XX
      FTOVV2_HDEC=FTOVV_HDEC(YY)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FTOVV_HDEC(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LAMB
      DIMENSION XX(2)
      COMMON/VVOFF_HDEC/AMH,AMW,GAMW
      LAMB(X,Y)=DSQRT((1.D0-X-Y)**2-4.D0*X*Y)
      PI=4D0*DATAN(1D0)
      ICASE = 1
      IF(ICASE.EQ.0)THEN
       YY = AMH**2
       Y1 = DATAN((YY-AMW**2)/AMW/GAMW)
       Y2 = -DATAN((AMW**2)/AMW/GAMW)
       DJAC = Y1-Y2
       T1 = TAN(Y1*XX(1)+Y2*(1.D0-XX(1)))
       SP = AMW**2 + AMW*GAMW*T1
       YY = (AMH-DSQRT(SP))**2
       Y1 = DATAN((YY-AMW**2)/AMW/GAMW)
       Y2 = -DATAN((AMW**2)/AMW/GAMW)
       DJAC = DJAC*(Y1-Y2)
       T2 = TAN(Y1*XX(2)+Y2*(1.D0-XX(2)))
       SM = AMW**2 + AMW*GAMW*T2
       AM2=AMH**2
       GAM = AM2*LAMB(SP/AM2,SM/AM2)*(1+LAMB(SP/AM2,SM/AM2)**2*AMH**4
     .                               /SP/SM/12)
       PRO1 = SP/AMW**2
       PRO2 = SM/AMW**2
       FTOVV_HDEC = PRO1*PRO2*GAM*DJAC/PI**2
      ELSE
       SP = AMH**2*XX(1)
       SM = (AMH-DSQRT(SP))**2*XX(2)
       DJAC = AMH**2*(AMH-DSQRT(SP))**2/PI**2
       AM2=AMH**2
       GAM = AM2*LAMB(SP/AM2,SM/AM2)*(1+LAMB(SP/AM2,SM/AM2)**2*AMH**4
     .                               /SP/SM/12)
       PRO1 = SP*GAMW/AMW/((SP-AMW**2)**2+AMW**2*GAMW**2)
       PRO2 = SM*GAMW/AMW/((SM-AMW**2)**2+AMW**2*GAMW**2)
       FTOVV_HDEC = PRO1*PRO2*GAM*DJAC
      ENDIF
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR HSM ---> TT* ---> TBW
C     ************************************************************
      SUBROUTINE HTOTTS_HDEC(AMH,AMT,AMB,AMW,HTTS)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNSTT1_HDEC
      COMMON/IKSY0_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP0_HDEC/AMH0,AMT0,AMB0,AMW0
      AMH0=AMH
      AMT0=AMT
      AMB0=AMB
      AMW0=AMW
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMH) GOTO 12
      ECM=AMH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1_HDEC(FUNSTT1_HDEC,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      HTTS=XSEC
12    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNSTT1_HDEC(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY0_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNSTT2_HDEC
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNSTT1_HDEC=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2_HDEC(FUNSTT2_HDEC,D,U,SS)
      FUNSTT1_HDEC=FUNSTT1_HDEC+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNSTT2_HDEC(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY0_HDEC/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMSTT_HDEC(SS)
      FUNSTT2_HDEC=SS
      RETURN
      END

      SUBROUTINE ELEMSTT_HDEC(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY0_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP0_HDEC/AMH,AMT,AMB,AMW
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW0,GAMZ0
      GAMT=GAMT0**2*AMT**2/AMH**4
      GAMW=GAMW0**2*AMW**2/AMH**4
      W=AMW**2/AMH**2
      T=AMT**2/AMH**2
      Y1=1-X2
      Y2=1-X1
      X0=2.D0-X1-X2
      W1=(1.D0-X2)
      W3=(1.-X1-X2)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W33=1.D0/(W3**2+GAMW**2)
      W13=W1*W3*W11*W33

      R11=4*T*W-16.*T*W*Y1-4.*T*Y2*Y1+8.*T*Y1+32.*T*W**2-20
     . .*T*Y1**2+8.*W*Y2*Y1+4.*W*Y1**2-4.*Y2*Y1**2-16.*T**2*W-
     .  32.*T**2*Y1+4.*T**2-16.*T**3-8.*W**2+4.*Y1**2-4.*Y1**3
      R33=-4.*T*W+4.*T*W*Y2-2.*T*W*Y2*Y1+4.*T*W*Y1+T*W*Y2**2-
     .  3.*T*W*Y1**2+2.*T*Y2*Y1-3.*T*Y2*Y1**2+4.*T*W**2-4.*T*W**3
     .  +T*Y2**2-3.*T*Y2**2*Y1-T*Y2**3+T*Y1**2-T*Y1**3+4.*T**2
     .  *W-4.*T**2*W*Y2-4.*T**2*W*Y1-2.*T**2*Y2*Y1-4.*T**2*W**2-
     .  T**2*Y2**2-T**2*Y1**2+4.*W**2*Y2*Y1-8.*W**3*Y2-8.*W**3*Y1
     .  +4.*W**3+8.*W**4
      R13=8.*W-24.*T*W+16.*T*W*Y1 -4.*T*Y2+16.*T*Y2*Y1-4.*T*
     .  Y1+16.*T*W**2+4.*T*Y2**2+12.*T*Y1**2-8.*W*Y2-12.*W*Y2*Y1
     .  -8.*W*Y1+4.*W*Y1**2-4.*Y2*Y1+8.*Y2*Y1**2+16.*T**2*W+8.
     .  *T**2*Y2+8.*T**2*Y1+16.*W**2*Y2+24.*W**2*Y1+4.*Y2**2*Y1-
     .  32.*W**3-4.*Y1**2+4.*Y1**3
      RES=R11*W11+4.D0*R33*W33/T-2.D0*R13*W13
      RETURN
      END

C     **************************************************
C     SUBROUTINE FOR A -> TT* -> TBW
C     **************************************************

      SUBROUTINE ATOTT_HDEC(AMA,AMT,AMB,AMW,AMCH,ATT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNATT1_HDEC
      COMMON/IKSY1_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP1_HDEC/AMA1,AMT1,AMB1,AMW1,AMCH1
      AMA1=AMA
      AMT1=AMT
      AMB1=AMB
      AMW1=AMW
      AMCH1=AMCH
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C        FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMA) GOTO 12
      ECM=AMA
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1_HDEC(FUNATT1_HDEC,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      ATT0=XSEC
 12   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNATT1_HDEC(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY1_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNATT2_HDEC
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNATT1_HDEC=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2_HDEC(FUNATT2_HDEC,D,U,SS)
      FUNATT1_HDEC=FUNATT1_HDEC+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNATT2_HDEC(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY1_HDEC/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMATT_HDEC(SS)
      FUNATT2_HDEC=SS
      RETURN
      END

      SUBROUTINE ELEMATT_HDEC(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY1_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP1_HDEC/AMA,AMT,AMB,AMW,AMCH
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      GAMT=GAMT1**2*AMT**2/AMA**4
      GAMC=GAMC0**2*AMCH**2/AMA**4
      CH=AMCH**2/AMA**2
      W=AMW**2/AMA**2
      T=AMT**2/AMA**2
      Y1=1-X1
      Y2=1-X2
      X0=2.D0-X1-X2
      W1=(1.D0-x2)
      W2=(1.D0-X0+W-CH)
      W22=1.D0/ ((1.D0-X0+W-CH)**2+GAMC)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W12=W1*W2*W11*W22
      R11=4.D0*T*W-4.D0*T*Y1*Y2+8.D0*T*Y2-4.D0*T*Y2**2+8.D0*W*Y1*Y2+4.D0
     .  *W*Y2**2-4.D0*Y1*Y2**2+4.D0*T**2-8.D0*W**2+4.D0*Y2**2-4.D0*Y2**3
      R22=-16.D0*W+16.D0*T*W-8.D0*T*Y1*Y2-4.D0*T*Y1**2-4.D0*T*Y2**2+16.
     .D0*W*Y1+8.D0*W*Y1*Y2+16.D0*W*Y2+4.D0*W*Y1**2+4.D0*W*Y2**2+8.D0*Y1*
     . Y2-12.D0*Y1*Y2**2-12.D0*Y1**2*Y2-16.D0*W**2+4.D0*Y1**2-4.D0*Y1**3
     . +4.D0*Y2**2-4.D0*Y2**3
      R12=16.D0*W-16.D0*T*W-8.D0*T*Y1+16.D0*T*Y1*Y2-8.D0*T*Y2+8.D0*T*Y1
     . **2+8.D0*T*Y2**2-16.D0*W*Y1-8.D0*W*Y1*Y2-16.D0*W*Y2-8.D0*W*Y2**2-
     . 8.D0*Y1*Y2+16.D0*Y1*Y2**2+8.D0*Y1**2*Y2+16.D0*W**2-8.D0*Y2**2
     . +8.D0*Y2**3
      RES=R11*W11+R22*W22+R12*W12
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR H ---> TT* ---> TBW
C     ************************************************************
      SUBROUTINE HTOTT_HDEC(AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV,HTT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNHTT1_HDEC
      COMMON/IKSY2_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP2_HDEC/AMH2,AMT2,AMB2,AMW2,AMCH2,TB2,GHT2,GAT2,GHVV2
      AMH2=AMH
      AMT2=AMT
      AMB2=AMB
      AMW2=AMW
      AMCH2=AMCH
      TB2=TB
      GHT2=GHT
      GAT2=GAT
      GHVV2=GHVV
      IP=5
      M1=AMB
      M2=AMT
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMH) GOTO 12
      ECM=AMH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1_HDEC(FUNHTT1_HDEC,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      HTT0=XSEC
 12   CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNHTT1_HDEC(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY2_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNHTT2_HDEC
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNHTT1_HDEC=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2_HDEC(FUNHTT2_HDEC,D,U,SS)
      FUNHTT1_HDEC=FUNHTT1_HDEC+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNHTT2_HDEC(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY2_HDEC/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMHTT_HDEC(SS)
      FUNHTT2_HDEC=SS
      RETURN
      END

      SUBROUTINE ELEMHTT_HDEC(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY2_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP2_HDEC/AMH,AMT,AMB,AMW,AMCH,TB,GHT,GAT,GHVV
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW0,GAMZ0
      GAMT=GAMT1**2*AMT**2/AMH**4
      GAMC=GAMC0**2*AMCH**2/AMH**4
      GAMW=GAMW0**2*AMW**2/AMH**4
      CH=AMCH**2/AMH**2
      W=AMW**2/AMH**2
      T=AMT**2/AMH**2
      Y1=1-X2
      Y2=1-X1
      X0=2.D0-X1-X2
      W1=(1.D0-X2)
      W2=(1.D0-X0+W-CH)
      W3=-(1.-X1-X2)
      W22=1.D0/ ((1.D0-X0+W-CH)**2+GAMC)
      W11=1.D0/((1.D0-X2)**2+GAMT)
      W33=1.D0/(W3**2+GAMW**2)
      W12=W1*W2*W11*W22
      W13=W1*W3*W11*W33
      W23=W2*W3*W22*W33

      R11=4*T*W-16.*T*W*Y1-4.*T*Y2*Y1+8.*T*Y1+32.*T*W**2-20
     . .*T*Y1**2+8.*W*Y2*Y1+4.*W*Y1**2-4.*Y2*Y1**2-16.*T**2*W-
     .  32.*T**2*Y1+4.*T**2-16.*T**3-8.*W**2+4.*Y1**2-4.*Y1**3
      R22=-16.*W+16.*T*W-8.*T*Y2*Y1-4.*T*Y2**2-4.*T*Y1**2+16
     .  .*W*Y2 + 8.*W*Y2*Y1 + 16.*W*Y1 + 4.*W*Y2**2 + 4.*W*Y1**2+8.*Y2*
     .  Y1-12.*Y2*Y1**2-12.*Y2**2*Y1-16.*W**2+4.*Y2**2-4.*Y2**3
     .  +4.*Y1**2-4.*Y1**3
      R33=-4.*T*W+4.*T*W*Y2-2.*T*W*Y2*Y1+4.*T*W*Y1+T*W*Y2**2-
     .  3.*T*W*Y1**2+2.*T*Y2*Y1-3.*T*Y2*Y1**2+4.*T*W**2-4.*T*W**3
     .  +T*Y2**2-3.*T*Y2**2*Y1-T*Y2**3+T*Y1**2-T*Y1**3+4.*T**2
     .  *W-4.*T**2*W*Y2-4.*T**2*W*Y1-2.*T**2*Y2*Y1-4.*T**2*W**2-
     .  T**2*Y2**2-T**2*Y1**2+4.*W**2*Y2*Y1-8.*W**3*Y2-8.*W**3*Y1
     .  +4.*W**3+8.*W**4
      R12=-16.*W+48.*T*W-16.*T*W*Y2+16.*T*W*Y1+8.*T*Y2-32.*T
     .  *Y2*Y1+8.*T*Y1-8.*T*Y2**2 - 24.*T*Y1**2+16.*W*Y2+8.*W*Y2*
     .  Y1+16.*W*Y1+8.*W*Y1**2+8.*Y2*Y1-16.*Y2*Y1**2-16.*T**2*Y2
     .  -16.*T**2*Y1-8.*Y2**2*Y1-16.*W**2+8.*Y1**2-8.*Y1**3
      R13=8.*W-24.*T*W+16.*T*W*Y1 -4.*T*Y2+16.*T*Y2*Y1-4.*T*
     .  Y1+16.*T*W**2+4.*T*Y2**2+12.*T*Y1**2-8.*W*Y2-12.*W*Y2*Y1
     .  -8.*W*Y1+4.*W*Y1**2-4.*Y2*Y1+8.*Y2*Y1**2+16.*T**2*W+8.
     .  *T**2*Y2+8.*T**2*Y1+16.*W**2*Y2+24.*W**2*Y1+4.*Y2**2*Y1-
     .  32.*W**3-4.*Y1**2+4.*Y1**3
      R23=16.*W-16.*T*W+8.*T*W*Y2+8.*T*W*Y1+8.*T*Y2*Y1+4.*T*
     .  Y2**2+4.*T*Y1**2-16.*W*Y2-16.*W*Y1-4.*W*Y2**2+4.*W*Y1**2
     .  -8.*Y2*Y1+12.*Y2*Y1**2+8.*W**2*Y2-8.*W**2*Y1+12.*Y2**2*
     .  Y1-4.*Y2**2+4.*Y2**3-4.*Y1**2+4.*Y1**3
      GLVV=DSQRT(1.D0-GHVV**2)
      RES=GHT**2*R11*W11+GLVV**2*GAT**2*R22*W22+
     .    4.D0*GHVV**2*R33*W33/T+2.D0*GHT*GLVV*GAT*R12*W12+
     .    2.D0*GHT*GHVV*R13*W13+2.D0*GHVV*GLVV*GAT*R23*W23
      RETURN
      END

C     ************************************************************
C     SUBROUTINE FOR H+ ---> BT* ---> BBW
C     ************************************************************
      SUBROUTINE CTOTT_HDEC(AMCH,AMT,AMB,AMW,CTT0)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,K
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNCTT1_HDEC
      COMMON/IKSY3_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP3_HDEC/AMH3,AMT3,AMB3,AMW3
      AMH3=AMCH
      AMT3=AMT
      AMB3=AMB
      AMW3=AMW
      IP=5
      M1=AMB
      M2=AMB
      M3=AMW
C     FIRST INTEGRATE OVER X2, i.e. (1+3) SYSTEM
C        CHECK WHETHER ENOUGH PHASE SPACE
      MASTOT=M1+M2+M3
      IF(MASTOT.GE.AMCH) GOTO 12
      ECM=AMCH
      S=ECM**2
      U1=(ECM-M2)**2
      D1=(M1+M3)**2
      U=(S-D1+M2**2)/s
      D=(S-U1+M2**2)/s
      DEL=(U-D)/IP
      U=D+DEL
      XSEC=0.D0
      DO K=1,IP
      CALL QGAUS1_HDEC(FUNCTT1_HDEC,D,U,SS)
      D=U
      U=D+DEL
      XSEC=XSEC+SS
      ENDDO
      CTT0=XSEC
12    CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNCTT1_HDEC(XL)
      IMPLICIT REAL*8(A-Z)
      INTEGER IP,I
      COMMON/IKSY3_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/PREC1_HDEC/IP
      EXTERNAL FUNCTT2_HDEC
      X2=XL
      S13=S-S*X2+M2**2
      TEM=2.D0*DSQRT(S13)
      E2S=(S-S13-M2**2)/TEM
      E3S=(S13+M3**2-M1**2)/TEM
C     SECOND INTEGRAL OVER X1, i.e. (2+3) SYSTEM
      U1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)-DSQRT(E3S**2-M3**2))**2
      D1=(E2S+E3S)**2-(DSQRT(E2S**2-M2**2)+DSQRT(E3S**2-M3**2))**2
      U=(S-D1+M1**2)/s
      D=(S-U1+M1**2)/s
      DEL=(U-D)/IP
      FUNCTT1_HDEC=0.d0
      U=D+DEL
      DO I=1,IP
      CALL QGAUS2_HDEC(FUNCTT2_HDEC,D,U,SS)
      FUNCTT1_HDEC=FUNCTT1_HDEC+SS
      D=U
      U=D+DEL
      ENDDO
      RETURN
      END

      DOUBLE PRECISION FUNCTION FUNCTT2_HDEC(XK)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY3_HDEC/X1,X2,M1,M2,M3,ECM,S
      X1=XK
      CALL ELEMCTT_HDEC(SS)
      FUNCTT2_HDEC=SS
      RETURN
      END

      SUBROUTINE ELEMCTT_HDEC(RES)
      IMPLICIT REAL*8(A-Z)
      COMMON/IKSY3_HDEC/X1,X2,M1,M2,M3,ECM,S
      COMMON/TOP3_HDEC/AMCH,AMT,AMB,AMW
      COMMON/WZWDTH_HDEC/GAMC0,GAMT0,GAMT1,GAMW,GAMZ
      GAMT=GAMT1**2*AMT**2/AMCH**4
      W=AMW**2/AMCH**2
      T=AMT**2/AMCH**2
      B=AMB**2/AMCH**2
      RES=((1.D0-X1-W)*(1.D0-X2-W)+W*(X1+X2-1.D0+W))/
     .   ((1.D0-X2+B-T)**2+GAMT)
      RETURN
      END

C   *****************  INTEGRATION ROUTINE ***********************
C    Returns SS as integral of FUNC from A to B, by 10-point Gauss-
C    Legendre integration
      SUBROUTINE QGAUS1_HDEC(FUNC,A,B,SS)
      IMPLICIT REAL*8(A-Z)
      INTEGER J
      DIMENSION X(5),W(5)
      EXTERNAL FUNC
      DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     .  ,.8650633666D0,.9739065285D0/
      DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     .  ,.1494513491D0,.0666713443D0/
      XM=0.5D0*(B+A)
      XR=0.5D0*(B-A)
      SS=0.D0
      DO 11 J=1,5
        DX=XR*X(J)
        SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
11    CONTINUE
      SS=XR*SS
      RETURN
      END

C     Returns SS as integral of FUNC from A to B, by 10-point Gauss-
C      Legendre integration
      SUBROUTINE QGAUS2_HDEC(FUNC,A,B,SS)
      IMPLICIT REAL*8(A-Z)
      INTEGER J
      DIMENSION X(5),W(5)
      EXTERNAL FUNC
      DATA X/.1488743389D0,.4333953941D0,.6794095682D0
     .  ,.8650633666D0,.9739065285D0/
      DATA W/.2955242247D0,.2692667193D0,.2190863625D0
     .  ,.1494513491D0,.0666713443D0/
      XM=0.5D0*(B+A)
      XR=0.5D0*(B-A)
      SS=0.D0
      DO 11 J=1,5
        DX=XR*X(J)
        SS=SS+W(J)*(FUNC(XM+DX)+FUNC(XM-DX))
11    CONTINUE
      SS=XR*SS
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AMHAMA_HDEC(ICASE,MH,TANB)
C--CALCULATION OF PSEUDOSCALAR HIGGS MASS FROM HIGGS MASS MH
C--ICASE=0: MH=PSEUDOSCALAR MASS
C--ICASE=1: MH=LIGHT SCALAR MASS
C--ICASE=2: MH=HEAVY SCALAR MASS
C--ICASE=3: MH=CHARGED HIGGS MASS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/HMASS_HDEC/AMSM,AMA,AML,AMH,AMCH,AMAR
      IF(ICASE.EQ.0)THEN
       MA = MH
      ELSE
       DEL0 = 1.D-4
       MA0 = 1.D0
       MA1 = 1.D4
1      MA = (MA0+MA1)/2
C      CALL SUBH1_HDEC(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MA
       CALL SUSYCP_HDEC(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       DEL = DABS(MA1 - MA0)/MA
       IF(DEL.GT.DEL0) THEN
        IF(MX.GT.MH) MA1 = MA
        IF(MX.LT.MH) MA0 = MA
        GOTO 1
       ENDIF
       FAC = 1
       MAX = DINT(FAC*MA+0.5D0)/FAC
C      CALL SUBH1_HDEC(MAX,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                 MHP,HMP,MCH,SA,CA,TANBA)
       AMA = MAX
       CALL SUSYCP_HDEC(TANB)
       IF(ICASE.EQ.1)THEN
        MX = AML
       ELSEIF(ICASE.EQ.2)THEN
        MX = AMH
       ELSEIF(ICASE.EQ.3)THEN
        MX = AMCH
       ENDIF
       IF(MX.EQ.MH)THEN
        MA = MAX
       ELSE
        DEL0 = 1.D-8
2       MA = (MA0+MA1)/2
C       CALL SUBH1_HDEC(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
C    *                  MHP,HMP,MCH,SA,CA,TANBA)
        AMA = MA
        CALL SUSYCP_HDEC(TANB)
        IF(ICASE.EQ.1)THEN
         MX = AML
        ELSEIF(ICASE.EQ.2)THEN
         MX = AMH
        ELSEIF(ICASE.EQ.3)THEN
         MX = AMCH
        ENDIF
        DEL = DABS(MA1 - MA0)/MA
        IF(DEL.GT.DEL0) THEN
         IF(MX.GT.MH) MA1 = MA
         IF(MX.LT.MH) MA0 = MA
         GOTO 2
        ENDIF
       ENDIF
      ENDIF
      AMA = MA
      CALL SUSYCP_HDEC(TANB)
      RETURN
      END
C
C     ****************************************************************
C	  CHARGINO AND NEUTRALINO MASS MATRICES AND COUPLINGS
C     ****************************************************************
      SUBROUTINE GAUGINO_HDEC(MU,M2,B,A,MC,MN,XMN,AC1,AC2,AC3,AN1,AN2
     .                 ,AN3,ACNL,ACNR,AGDL,AGDA,AGDH,AGDC)            
      IMPLICIT REAL*8(A-H,K-Z)
      COMPLEX*16 CXA,CXB,CXC,CXD,CX1,CX2,CX3
      DIMENSION MC(2),MN(4),XMN(4),Z(4,4),ZX(4,4),U(2,2),V(2,2),
     .          QQ(4,4),SS(4,4),S(2,2),Q(2,2),AC1(2,2),AC2(2,2),
     .          AC3(2,2),AN1(4,4),AN2(4,4),AN3(4,4),ACNL(2,4),
     .          ACNR(2,4),IORD(4),IREM(2)
      DIMENSION X(2,2)
      DIMENSION YMN(4),YZ(4,4),XMC(2),BU(2),BV(2)
      DIMENSION AGDL(4),AGDA(4),AGDH(4),AGDC(2)
      DIMENSION slhaneut(4),slhaxneut(4),slhachar(2),slhau(2,2),
     .          slhav(2,2),slhaz(4,4),slhaxchar(2)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,MZ,MW
      COMMON/GAUGINOMIX_HDEC/ZX,U,V
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/SLHA_m1_HDEC/am1
      COMMON/SLHA_gaug_HDEC/slhaneut,slhaxneut,slhachar,slhau,slhav,
     .                      slhaz,slhaxchar
      CW=MW/MZ
      SW=DSQRT(1-CW**2)
      PI=4.D0*DATAN(1.D0)
      SB=DSIN(B)
      CB=DCOS(B)
      TW=SW/CW
      if(islhai.eq.0) then
         M1=5.D0/3.D0*TW**2*M2
      else
         M1 = am1
      endif
C     ************  NEUTRALINO MASSES AND MATRIX ELEMENTS ***********
      EPS=-1.D-10
      XC2=(M1*M2-MZ**2-MU**2)-3.D0/8.D0*(M1+M2)**2
      XC3=-1.D0/8.D0*(M1+M2)**3+1.D0/2.D0*(M1+M2)*(M1*M2-MZ**2
     .    -MU**2)+(M1+M2)*MU**2+(M1*CW**2+M2*SW**2)*MZ**2
     .    -MU*MZ**2*DSIN(2.D0*B)
      XC4=+(M1*CW**2+M2*SW**2)*MU*MZ**2*DSIN(2.D0*B)-M1*M2*MU**2
     .    +1.D0/4.D0*(M1+M2)*( (M1+M2)*MU**2+(M1*CW**2+M2*SW**2)
     .    *MZ**2-MU*MZ**2*DSIN(2.D0*B) )+1.D0/16.D0*(M1+M2)**2*
     .    (M1*M2-MZ**2-MU**2)-3.D0/256.D0*(M1+M2)**4
      XS=-XC3**2-2.D0/27.D0*XC2**3+8.D0/3.D0*XC2*XC4
      XU=-1.D0/3.D0*XC2**2-4.D0*XC4
      CXD=(-4*XU**3-27*XS**2)*DCMPLX(1.D0,EPS)
      CXC=1.D0/2.D0*(-XS+DCMPLX(0.D0,1.D0)*CDSQRT(CXD/27.D0))
      CXA=DREAL(CXC**(1.D0/3.D0))*DCMPLX(1.D0,-EPS)
      CXB=8.D0*CXA-8.D0/3.D0*XC2*DCMPLX(1.D0,-EPS)
C     *********** MASSES AND COUPLINGS:
      if(islhai.eq.0) then
         X0=(M1+M2)/4.D0
         CX1= CXA/2.D0-XC2/6.D0*DCMPLX(1.D0,-EPS)
         CX2=-CXA/2.D0-XC2/3.D0*DCMPLX(1.D0,-EPS)
         CX3=XC3*DCMPLX(1.D0,-EPS)/CDSQRT(CXB)
         XMN(1)=X0-CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2+CX3))
         XMN(2)=X0+CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2-CX3))
         XMN(3)=X0-CDABS(CDSQRT(CX1))-CDABS(CDSQRT(CX2+CX3))
         XMN(4)=X0+CDABS(CDSQRT(CX1))+CDABS(CDSQRT(CX2-CX3))
         DO 10 I=1,4
            MN(I)=DABS(XMN(I))
            YMN(I)=XMN(I)
            ZX(I,2)=-CW/SW*(M1-XMN(I))/(M2-XMN(I))
            ZX(I,3)=(MU*(M2-XMN(I))*
     .           (M1-XMN(I))-MZ**2*SB*CB*((M1-M2)*CW**2
     .           +M2-XMN(I)))/MZ/(M2-XMN(I))/SW/(MU*CB+XMN(I)*SB)
            ZX(I,4)=(-XMN(I)*(M2-XMN(I))*(M1-XMN(I))-MZ**2*CB*CB*
     .           ((M1-M2)*CW**2+M2-XMN(I)))/MZ/(M2-XMN(I))
     .           /SW/(MU*CB+XMN(I)*SB)
            ZX(I,1)=1.D0/DSQRT(1.D0+ZX(I,2)**2+ZX(I,3)**2+ZX(I,4)**2) 
            YZ(I,1)=ZX(I,1)
            YZ(I,2)=ZX(I,2)*ZX(I,1)
            YZ(I,3)=ZX(I,3)*ZX(I,1)
            YZ(I,4)=ZX(I,4)*ZX(I,1)
 10      CONTINUE
      else
         do i=1,4,1
            xmn(i) = slhaxneut(i)
            mn(i)  = dabs(xmn(i))
            ymn(i) = xmn(i)
            do j=1,4,1
               zx(i,j) = slhaz(i,j)
            end do
            yz(i,1)=zx(i,1)
            yz(i,2)=zx(i,2)
            yz(i,3)=zx(i,3)
            yz(i,4)=zx(i,4)
         end do
      endif
C     *************  ORDERING THE DISORDER ******************
      XX0 = DMIN1(MN(1),MN(2),MN(3),MN(4))
      XX1 = DMAX1(MN(1),MN(2),MN(3),MN(4))
      IDUMMY = 1
      DO I = 1,4
       IF(MN(I).EQ.XX0)THEN
        IORD(1) = I
       ELSEIF(MN(I).EQ.XX1)THEN
        IORD(4) = I
       ELSE
        IREM(IDUMMY) = I
        IDUMMY = IDUMMY+1
       ENDIF
      ENDDO
      IF(MN(IREM(1)).LE.MN(IREM(2)))THEN
       IORD(2) = IREM(1)
       IORD(3) = IREM(2)
      ELSE
       IORD(2) = IREM(2)
       IORD(3) = IREM(1)
      ENDIF
C 
      DO 98 J=1,4
      I=IORD(J)
      XMN(J)=YMN(I)
      MN(J) =DABS(YMN(I))
        DO I1=1,4
        Z(J,I1)=YZ(I,I1)
        ENDDO
 98   CONTINUE
C     ************  NEUTRALINO COUPLINGS TO HIGGS BOSONS ***********
	DO 11 I=1,4
	DO 11 J=1,4
	QQ(I,J)=1.D0/2.D0*(Z(I,3)*(Z(J,2)-TW*Z(J,1))+Z(J,3)*
     .		(Z(I,2)-TW*Z(I,1)))
	SS(I,J)=1.D0/2.D0*(Z(I,4)*(Z(J,2)-TW*Z(J,1))+Z(J,4)*
     .		(Z(I,2)-TW*Z(I,1)))
 11	CONTINUE
	DO 21 I=1,4
	DO 21 J=1,4
	AN1(I,J)= QQ(I,J)*DCOS(A)-SS(I,J)*DSIN(A)
	AN2(I,J)=-QQ(I,J)*DSIN(A)-SS(I,J)*DCOS(A)
	AN3(I,J)= QQ(I,J)*DSIN(B)-SS(I,J)*DCOS(B)
 21	CONTINUE

C       ************* CHARGINO MASSES AND MATRIX ELEMENTS ***********
        if(islhai.eq.0) then
           DELTA=DABS(B-.25*PI)
           DDD=MU*DCOS(B)+M2*DSIN(B)
           CCC=MU*DSIN(B)+M2*DCOS(B)
           IF(DELTA.LT.0.01D0) THEN
              PHIM=PI/4.D0-.5D0*DATAN((M2-MU)/(2.D0*MW))
              PHIP=PHIM
           ELSE IF	(DABS(CCC).LT.1.D-5) THEN
              PHIM=0.D0
              PHIP=DATAN(DSQRT(2.D0)*MW*DSIN(B)/(M2+1.D-5))
           ELSE IF	(DABS(DDD).LT.1.D-5) THEN
              PHIP=0.D0
              PHIM=DATAN(DSQRT(2.D0)*MW*DCOS(B)/(M2+1.D-5))
           ELSE
              RAD=DSQRT((M2**2-MU**2)**2+4.D0*MW**4*DCOS(2.D0*B)**2
     +             +4.D0*MW**2*(M2**2+MU**2+2.D0*M2*MU*DSIN(2.D0*B)))
              PHIP=DATAN((RAD-(M2**2-MU**2+2.D0*MW**2*DCOS(2.D0*B)))
     +             /(2.D0*DSQRT(2.D0)*MW*(MU*DCOS(B)+M2*DSIN(B))))
              PHIM=DATAN((RAD-(M2**2-MU**2-2.D0*MW**2*DCOS(2.D0*B)))
     +             /(2.D0*DSQRT(2.D0)*MW*(MU*DSIN(B)+M2*DCOS(B))))
           ENDIF
           CP=DCOS(PHIP)
           SP=DSIN(PHIP)
           CM=DCOS(PHIM)
           SM=DSIN(PHIM)
C MY CONVENTION
           U(2,2)=CM
           U(2,1)=-SM
           U(1,2)=SM
           U(1,1)=CM
           V(1,1)=CP
           V(1,2)=SP
           V(2,1)=-SP
           V(2,2)=CP
           X(1,1)=M2
           X(1,2)=DSQRT(2.D0)*MW*DSIN(B)
           X(2,1)=DSQRT(2.D0)*MW*DCOS(B)
           X(2,2)=MU
 555       CONTINUE
           XMC(1)=(U(1,1)*X(1,1)+U(1,2)*X(2,1))*V(1,1)
     .          +(U(1,1)*X(1,2)+U(1,2)*X(2,2))*V(1,2)
           XMC(2)=(U(2,1)*X(1,1)+U(2,2)*X(2,1))*V(2,1)
     .          +(U(2,1)*X(1,2)+U(2,2)*X(2,2))*V(2,2)
           IF(XMC(1).LT.0.D0) THEN
              V(1,1)=-CP
              V(1,2)=-SP
              V(2,1)=-SP
              V(2,2)=CP
              GOTO 555
           ENDIF
           IF(XMC(2).LT.0.D0) THEN
              V(1,1)=CP
              V(1,2)=SP
              V(2,1)=SP
              V(2,2)=-CP
              GOTO 555
           ENDIF
           IF(XMC(1).GT.XMC(2)) THEN
              MTEMP=XMC(1)
              XMC(1)=XMC(2)
              XMC(2)=MTEMP
              DO J=1,2
                 BU(J)=U(1,J)
                 U(1,J)=U(2,J)
                 U(2,J)=BU(J)
                 BV(J)=V(1,J)
                 V(1,J)=V(2,J)
                 V(2,J)=BV(J)
              ENDDO
           ENDIF        
           MC(1)=DABS(XMC(1))
           MC(2)=DABS(XMC(2))
           slhaxchar(1) = mc(1)
           slhaxchar(2) = mc(2)
        else
           mc(1) = slhachar(1)
           mc(2) = slhachar(2)
           do i=1,2,1
              do j=1,2,1
                 u(i,j) = slhau(i,j)
                 v(i,j) = slhav(i,j)
              end do
           end do
        endif

C     ************  CHARGINO COUPLINGS TO HIGGS BOSONS ***********
	DO 12 I=1,2
	DO 12 J=1,2
	Q(I,J)=DSQRT(1.D0/2.D0)*U(J,2)*V(I,1)
	S(I,J)=DSQRT(1.D0/2.D0)*U(J,1)*V(I,2)
 12	CONTINUE
	DO 22 I=1,2
	DO 22 J=1,2	
	AC1(I,J)= Q(I,J)*DCOS(A)+S(I,J)*DSIN(A)
	AC2(I,J)=-Q(I,J)*DSIN(A)+S(I,J)*DCOS(A)
	AC3(I,J)= Q(I,J)*DSIN(B)+S(I,J)*DCOS(B)
 22	CONTINUE
C     **** CHARGINO-NEUTRALINO COUPLINGS TO CHARGED HIGGS BOSONS 
	DO 13 I=1,2
	DO 13 J=1,4
        ACNL(I,J)=DCOS(B)*(Z(J,4)*V(I,1)+(Z(J,2)+Z(J,1)*TW)
     .       *V(I,2)/DSQRT(2.D0)) 
        ACNR(I,J)=DSIN(B)*(Z(J,3)*U(I,1)-(Z(J,2)+Z(J,1)*TW)
     .       *U(I,2)/DSQRT(2.D0)) 
 13     CONTINUE

C   ************* HIGGS--NEUTRALINO--GOLDSTINO COUPLINGS
      DO 51 I=1,4
      AGDL(I)=Z(I,3)*DSIN(A)-Z(I,4)*DCOS(A)
      AGDH(I)=Z(I,3)*DCOS(A)+Z(I,4)*DSIN(A)
      AGDA(I)=Z(I,3)*DSIN(B)+Z(I,4)*DCOS(B)
 51   CONTINUE
C
C   ************* CHARGED HIGGS--CHARGINO--GOLDSTINO COUPLINGS
      AGDC(1)=DSQRT( V(1,2)**2*DCOS(B)**2+ U(1,2)**2*DSIN(B)**2 )
      AGDC(2)=DSQRT( V(2,2)**2*DCOS(B)**2+ U(2,2)**2*DSIN(B)**2 )

       RETURN
       END

C   ****************************************************************
C     SUBROUTINE FOR SFERMION MASSES, MIXING AND COUPLINGS 
C   ****************************************************************

      SUBROUTINE SFERMION_HDEC(TSC,BSC,MQL,MUR,MDR,MEL,MER,AL,AT,AB,MU,
     .                    MST,MSB,MSL,MSU,MSD,MSE,MSN, 
     .                    GLEE,GLTT,GLBB,GHEE,GHTT,GHBB,
     .                    GAEE,GATT,GABB,GCEN,GCTB)

      IMPLICIT REAL*8(A-H,K-Z)
      DIMENSION MST(2),MSB(2),MSL(2),MSU(2),MSD(2),MSE(2),MSN(2),
     .          GCEN(2,2),GCTB(2,2),GLEE(2,2),GLTT(2,2),GLBB(2,2),
     .          GHEE(2,2),GHTT(2,2),GHBB(2,2)
      DIMENSION slhast(2),slhasb(2),slhasu(2),slhasd(2),slhase(2),
     .          slhasl(2),slhasn(2),slhasnl(2)
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,MZ,MW
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,GHHH,GLLL,GHLL,
     .            GLHH,GHAA,GLAA,GLVV,GHVV,GLPM,GHPM,B,A
      COMMON/SFER1ST_HDEC/MQL1,MUR1,MDR1,MEL1,MER1
      COMMON/GLUINO_HDEC/AMGLUINO,XMSB1,XMSB2,STHB,CTHB,
     .              XLBB(2,2),XHBB(2,2),XABB(2,2),
     .              XMST1,XMST2,STHT,CTHT,
     .              XLTT(2,2),XHTT(2,2),XATT(2,2)
      COMMON/TAUMIX_HDEC/CL,SL
      COMMON/SLHA_vals_HDEC/islhai,islhao
      COMMON/SLHA_sfer_HDEC/slhast,slhasb,slhasu,slhasd,slhase,slhasl,
     .                 slhasn,slhasnl,slhacot,slhasit,slhacob,slhasib,
     .                 slhacol,slhasil
C
      PI = 4*DATAN(1.D0)
      SW2=1.D0-MW**2/MZ**2
      TB=DTAN(B)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MT = AMT
      MB = AMB
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     MT = RUNM_HDEC(TSC,6)
c     MB = RUNM_HDEC(BSC,5)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     MT = RUNM_HDEC(AMT,6)
c     MB = RUNM_HDEC(AMT,5)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ML = AMTAU
C FIRST TWO GENERATIONS:  NO MIXING INCLUDED 
      if(islhai.eq.0) then
C UP SQUARKS: 
         MSTL2=MQL1**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
         MSTR2=MUR1**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
         MSU(1)=DSQRT(MSTL2)
         MSU(2)=DSQRT(MSTR2)
C DOWN SQUARKS
         MSBL2=MQL1**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
         MSBR2=MDR1**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
         MSD(1)=DSQRT(MSBL2)
         MSD(2)=DSQRT(MSBR2)
C SLEPTONS
         MSEL2=MEL1**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
         MSER2=MER1**2- SW2*MZ**2*DCOS(2.D0*B) 
         MSNL2=MEL1**2+0.5D0*MZ**2*DCOS(2.D0*B)
         MSE(1)=DSQRT(MSEL2)
         MSE(2)=DSQRT(MSER2)
         MSN(1)=DSQRT(MSNL2)
         MSN(2)=1.D+15

C NOW THE THIRD GENERATION
C
C STOP MASSES/MIXING
C
      MSTL2=MQL**2+(0.5D0-2.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSTR2=MUR**2+2.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRT=AT-MU/TB
      DELT=(MSTL2-MSTR2)**2+4*MT**2*MLRT**2
      MST12=MT**2+0.5D0*(MSTL2+MSTR2-DSQRT(DELT))
      MST22=MT**2+0.5D0*(MSTL2+MSTR2+DSQRT(DELT))
        IF(MST12.LT.0.D0)THEN 
      PRINT *, 'MSTOP**2 is negative!!!!'
      GOTO 111 
      ELSE 
      MST(1)=DSQRT(MST12)
      MST(2)=DSQRT(MST22)
      IF(MSTL2.EQ.MSTR2) THEN
       THET = PI/4
      ELSE
       THET=0.5D0*DATAN(2.D0*MT*MLRT / (MSTL2-MSTR2) )
       IF(MSTL2.GT.MSTR2) THET = THET + PI/2
      ENDIF
        ENDIF 
      CT= DCOS(THET)
      ST= DSIN(THET) 
C
C SBOTTOM MASSES/MIXING
C
      MSBL2=MQL**2+(-0.5D0+1.D0/3.D0*SW2)*MZ**2*DCOS(2.D0*B)
      MSBR2=MDR**2-1.D0/3.D0*SW2*MZ**2*DCOS(2.D0*B) 
      MLRB=AB-MU*TB
      DELB=(MSBL2-MSBR2)**2+4*MB**2*MLRB**2
      MSB12=MB**2+0.5D0*(MSBL2+MSBR2-DSQRT(DELB))
      MSB22=MB**2+0.5D0*(MSBL2+MSBR2+DSQRT(DELB))
        IF(MSB12.LT.0.D0)THEN
      PRINT *, 'MSBOT**2 is negative!!!!'
      GOTO 111
        ELSE
      MSB(1)=DSQRT(MSB12)
      MSB(2)=DSQRT(MSB22)
      IF(MSBL2.EQ.MSBR2) THEN
       THEB = PI/4
      ELSE
       THEB=0.5D0*DATAN(2.D0*MB*MLRB / (MSBL2-MSBR2) )
       IF(MSBL2.GT.MSBR2) THEB = THEB + PI/2
      ENDIF
        ENDIF  
      CB= DCOS(THEB)
      SB= DSIN(THEB) 

C
C  STAU MASSES/MIXING
C
      MSEL2=MEL**2+(-0.5D0+SW2)*MZ**2*DCOS(2.D0*B)
      MSER2=MER**2- SW2*MZ**2*DCOS(2.D0*B) 
      MSNL2=MEL**2+0.5D0*MZ**2*DCOS(2.D0*B)
      MLRE=AL-MU*TB
      DELE=(MSEL2-MSER2)**2+4*ML**2*MLRE**2
      MSE12=ML**2+0.5D0*(MSEL2+MSER2-DSQRT(DELE))
      MSE22=ML**2+0.5D0*(MSEL2+MSER2+DSQRT(DELE))
        IF(MSE12.LT.0.D0)THEN
      PRINT *, 'MSTAU**2 is negative!!!!'
      GOTO 111
        ELSE
      MSL(1)=DSQRT(MSE12)
      MSL(2)=DSQRT(MSE22)
      IF(MSEL2.EQ.MSER2) THEN
       THEL = PI/4
      ELSE
       THEL=0.5D0*DATAN(2.D0*ML*MLRE / (MSEL2-MSER2) )
       IF(MSEL2.GT.MSER2) THEL = THEL + PI/2
      ENDIF
        ENDIF  
      CL= DCOS(THEL)
      SL= DSIN(THEL) 

      else
         do i=1,2,1
            msu(i) = slhasu(i)
            msd(i) = slhasd(i)
            mse(i) = slhase(i)
            msn(i) = slhasn(i)
            mst(i) = slhast(i)
            msb(i) = slhasb(i)
            msl(i) = slhasl(i)
         end do
         ct = slhacot
         st = slhasit
         cb = slhacob
         sb = slhasib
         cl = slhacol
         sl = slhasil

c maggie changed 28/12/2008
         thet=dacos(ct)
         theb=dacos(cb)
         thel=dacos(cl)
         
         if(slhasit.le.0.D0) then
            if(thet.ge.0.D0) then
               thet = -1.D0*thet
            endif
         endif

         if(slhasib.le.0.D0) then
            if(theb.ge.0.D0) then
               theb = -1.D0*theb
            endif
         endif
         
         if(slhasil.le.0.D0) then
            if(thel.ge.0.D0) then
               thel = -1.D0*thel
            endif
         endif
c end maggie changed 28/12/2008
      endif
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STOPS
C 
      GLTT(1,1)=-DSIN(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GLT + MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(2,2)=-DSIN(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GLT - MT*ST*CT/MZ**2*(AT*GLT+MU*GHT)
      GLTT(1,2)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
      GLTT(2,1)=-2*DSIN(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GLT+MU*GHT) 
C
C LIGHT CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GLBB(1,1)=-DSIN(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB + MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(2,2)=-DSIN(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GLB - MB*SB*CB/MZ**2*(AB*GLB-MU*GHB)
      GLBB(1,2)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 
      GLBB(2,1)=-2*DSIN(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GLB-MU*GHB) 

C
C LIGHT CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GLEE(1,1)=-DSIN(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB + ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(2,2)=-DSIN(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GLB - ML*SL*CL/MZ**2*(AL*GLB-MU*GHB)
      GLEE(1,2)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
      GLEE(2,1)=-2*DSIN(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GLB-MU*GHB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STOPS
C
      GHTT(1,1)=DCOS(B+A)*(0.5D0*CT**2-2.D0/3.D0*SW2*DCOS(2*THET)) 
     .    + MT**2/MZ**2*GHT + MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(2,2)= DCOS(B+A)*(0.5D0*ST**2+2.D0/3.D0*SW2*DCOS(2*THET))
     .    + MT**2/MZ**2*GHT - MT*ST*CT/MZ**2*(AT*GHT-MU*GLT)
      GHTT(1,2)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
      GHTT(2,1)=2*DCOS(B+A)*ST*CT*(2.D0/3.D0*SW2-0.25D0)
     .    + MT*DCOS(2*THET)/2.D0/MZ**2*(AT*GHT-MU*GLT) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO SBOTTOMS
C
      GHBB(1,1)= DCOS(B+A)*(-0.5D0*CB**2+1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB + MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(2,2)= DCOS(B+A)*(-0.5D0*SB**2-1.D0/3.D0*SW2*DCOS(2*THEB)) 
     .    + MB**2/MZ**2*GHB - MB*SB*CB/MZ**2*(AB*GHB+MU*GLB)
      GHBB(1,2)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
      GHBB(2,1)=2*DCOS(B+A)*SB*CB*(-1.D0/3.D0*SW2+0.25D0)
     .    + MB*DCOS(2*THEB)/2.D0/MZ**2*(AB*GHB+MU*GLB) 
C
C HEAVY CP--EVEN HIGGS COUPLINGS TO STAU'S 
C
      GHEE(1,1)= DCOS(B+A)*(-0.5D0*CL**2+SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB + ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(2,2)= DCOS(B+A)*(-0.5D0*SL**2-SW2*DCOS(2*THEL)) 
     .    + ML**2/MZ**2*GHB - ML*SL*CL/MZ**2*(AL*GHB+MU*GLB)
      GHEE(1,2)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 
      GHEE(2,1)=2*DCOS(B+A)*SL*CL*(-SW2+0.25D0)
     .    + ML*DCOS(2*THEL)/2.D0/MZ**2*(AL*GHB+MU*GLB) 

C
C PSEUDOSCALAR COUPLINGS 
C
      GATT=MT/2.D0/MZ**2*(MU+AT*GAT) 
      GABB=MB/2.D0/MZ**2*(MU+AB*GAB) 
      GAEE=ML/2.D0/MZ**2*(MU+AL*GAB) 
C
C CHARGED HIGGS COUPLINGS STOPS/SBOTTOMS 
C
      CLL=(MW**2*DSIN(2*B)-MT**2*GAT-MB**2*GAB)/DSQRT(2.D0)/MW**2
      CRR=-MT*MB*(GAT+GAB)/DSQRT(2.D0)/MW**2
      CLR=-MB*(MU+AB*GAB)/DSQRT(2.D0)/MW**2
      CRL=-MT*(MU+AT*GAT)/DSQRT(2.D0)/MW**2
      GCTB(1,1)=+CT*CB*CLL+ST*SB*CRR+CT*SB*CLR+ST*CB*CRL
      GCTB(1,2)=-CT*SB*CLL+ST*CB*CRR+CT*CB*CLR-ST*SB*CRL
      GCTB(2,1)=-ST*CB*CLL+CT*SB*CRR-ST*SB*CLR+CT*CB*CRL
      GCTB(2,2)=+ST*SB*CLL+CT*CB*CRR-ST*CB*CLR-CT*SB*CRL

C
C CHARGED HIGGS COUPLINGS TAU'S AND NEUTRINOS 
C
      CLL=(MW**2*DSIN(2*B)-ML**2*GAB)/DSQRT(2.D0)/MW**2
      CLR=-ML*(MU+AL*GAB)/DSQRT(2.D0)/MW**2
      GCEN(1,1)=CL*CLL+SL*CLR
      GCEN(1,2)=-SL*CLL+CL*CLR
      GCEN(2,1)=0.D0
      GCEN(2,2)=0.D0 

C--FILL COMMON BLOCK GLUINO_HDEC FOR SUSY-QCD CORRECTIONS TO
C  HIGGS -> BB, SQUARKS
      XMST1 = MST(1)
      XMST2 = MST(2)
      XMSB1 = MSB(1)
      XMSB2 = MSB(2)
      STHT = ST
      CTHT = CT
      STHB = SB
      CTHB = CB
       DO I=1,2
        DO J=1,2
         XLBB(I,J) = GLBB(I,J)
         XHBB(I,J) = GHBB(I,J)
         XABB(I,J) = 0
         XLTT(I,J) = GLTT(I,J)
         XHTT(I,J) = GHTT(I,J)
         XATT(I,J) = 0
        ENDDO
       ENDDO
       XABB(1,2) = GABB
       XABB(2,1) = -GABB
       XATT(1,2) = GATT
       XATT(2,1) = -GATT

      RETURN 
111   STOP
      END 

C ******************************************************************

C      DOUBLE PRECISION FUNCTION RUNP_HDEC(Q,NF)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON/RUN_HDEC/XMSB,XMHAT,XKFAC
C      RUNP_HDEC = RUNM_HDEC(Q,NF)
C      RUNP_HDEC = RUNM_HDEC(Q/2.D0,NF)*XKFAC
C      RETURN
C      END

      DOUBLE PRECISION FUNCTION RUNM_HDEC(Q,NF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594D0)
      DIMENSION AM(NN),YMSB(NN)
      COMMON/ALS_HDEC/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/STRANGE_HDEC/AMSB
      COMMON/RUN_HDEC/XMSB,XMHAT,XKFAC
      COMMON/FLAG_HDEC/IHIGGS,NNLO,IPOLE
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
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     C1(NF) = 1.175d0
c     C2(NF) = 1.501d0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      TRAN(X,XK)=1D0+4D0/3D0*ALPHAS_HDEC(X,2)/PI
     .              +XK*(ALPHAS_HDEC(X,2)/PI)**2
      CQ(X,NF)=(2D0*B0(NF)*X)**(G0(NF)/B0(NF))
     .            *(1D0+C1(NF)*X+C2(NF)*X**2)
      DATA ISTRANGE/0/
      PI=4D0*DATAN(1D0)
      ACC = 1.D-8
      AM(1) = 0
      AM(2) = 0
C--------------------------------------------
      IMSBAR = 0
      IF(IMSBAR.EQ.1)THEN
       IF(ISTRANGE.EQ.0)THEN
C--STRANGE POLE MASS FROM MSBAR-MASS AT 1 GEV
        AMSD = XLAMBDA
        AMSU = 1.D8
123     AMS  = (AMSU+AMSD)/2
        AM(3) = AMS
        XMSB = AMS/CQ(ALPHAS_HDEC(AMS,2)/PI,3)
     .            *CQ(ALPHAS_HDEC(1.D0,2)/PI,3)/TRAN(AMS,0D0)
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
       XMHAT = XMSB/CQ(ALPHAS_HDEC(AM(NF),2)/PI,NF)
      ELSE
       XMSB = 0
       XMHAT = 0
      ENDIF
      YMSB(3) = AMSB
      IF(NF.EQ.3)THEN
       YMSB(4) = YMSB(3)*CQ(ALPHAS_HDEC(AM(4),2)/PI,3)/
     .                   CQ(ALPHAS_HDEC(1.D0,2)/PI,3)
       YMSB(5) = YMSB(4)*CQ(ALPHAS_HDEC(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS_HDEC(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS_HDEC(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS_HDEC(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4) = XMSB
       YMSB(5) = YMSB(4)*CQ(ALPHAS_HDEC(AM(5),2)/PI,4)/
     .                   CQ(ALPHAS_HDEC(AM(4),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS_HDEC(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS_HDEC(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5) = XMSB
       YMSB(4) = YMSB(5)*CQ(ALPHAS_HDEC(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS_HDEC(AM(5),2)/PI,4)
       YMSB(6) = YMSB(5)*CQ(ALPHAS_HDEC(AM(6),2)/PI,5)/
     .                   CQ(ALPHAS_HDEC(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6) = XMSB
       YMSB(5) = YMSB(6)*CQ(ALPHAS_HDEC(AM(5),2)/PI,5)/
     .                   CQ(ALPHAS_HDEC(AM(6),2)/PI,5)
       YMSB(4) = YMSB(5)*CQ(ALPHAS_HDEC(AM(4),2)/PI,4)/
     .                   CQ(ALPHAS_HDEC(AM(5),2)/PI,4)
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
       XKFAC = 1D0
      ENDIF
      RUNM_HDEC = YMSB(N0)*CQ(ALPHAS_HDEC(Q,2)/PI,N0)/
     .               CQ(ALPHAS_HDEC(Q0,2)/PI,N0)
     .       * XKFAC
      RETURN
      END

      DOUBLE PRECISION FUNCTION ALPHAS_HDEC(Q,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM_HDEC/XLB1(6),XLB2(6)
      COMMON/ALS_HDEC/XLAMBDA,AMC,AMB,AMT,N0
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS1(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2))
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
        ALPHAS_HDEC=ALS1(NF,Q)
      ELSE
        ALPHAS_HDEC=ALS2(NF,Q)
      ENDIF
      RETURN
      END

      SUBROUTINE ALSINI_HDEC(ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XLB(6)
      COMMON/ALSLAM_HDEC/XLB1(6),XLB2(6)
      COMMON/ALS_HDEC/XLAMBDA,AMC,AMB,AMT,N0
      PI=4.D0*DATAN(1.D0)
      XLB1(1)=0D0
      XLB1(2)=0D0
      XLB2(1)=0D0
      XLB2(2)=0D0
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
      DO 1 I=3,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2.D0/25.D0)
     .             *(2.D0*DLOG(AMC/XLB(3)))**(-107.D0/1875.D0)
       XLB(4)=XITER_HDEC(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER_HDEC(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_HDEC(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMB)**(2.D0/23.D0)
     .             *(2.D0*DLOG(AMB/XLB(4)))**(-963.D0/13225.D0)
       XLB(5)=XITER_HDEC(AMB,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_HDEC(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_HDEC(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER_HDEC(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_HDEC(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2.D0/21.D0)
     .            *(2.D0*DLOG(AMT/XLB(5)))**(-321.D0/3381.D0)
       XLB(6)=XITER_HDEC(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2.D0/23.D0)
     .            *(2.D0*DLOG(AMT/XLB(6)))**(321.D0/3703.D0)
       XLB(5)=XITER_HDEC(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMB)**(-2.D0/25.D0)
     .             *(2.D0*DLOG(AMB/XLB(5)))**(963.D0/14375.D0)
       XLB(4)=XITER_HDEC(AMB,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2.D0/27.D0)
     .             *(2.D0*DLOG(AMC/XLB(4)))**(107.D0/2025.D0)
       XLB(3)=XITER_HDEC(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=3,6
       XLB2(I)=XLB(I)
2     CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITER_HDEC(Q,XLB1,NF1,XLB,NF2,ACC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2.D0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER_HDEC=XLB2
      RETURN
      END

      DOUBLE PRECISION FUNCTION FINT_HDEC(Z,XX,YY)
C--ONE-DIMENSIONAL CUBIC INTERPOLATION
C--Z  = WANTED POINT
C--XX = ARRAY OF 4 DISCRETE X-VALUES AROUND Z
C--YY = ARRAY OF 4 DISCRETE FUNCTION-VALUES AROUND Z
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(4),YY(4)
      X = DLOG(Z)
      X0=DLOG(XX(1))
      X1=DLOG(XX(2))
      X2=DLOG(XX(3))
      X3=DLOG(XX(4))
      Y0=DLOG(YY(1))
      Y1=DLOG(YY(2))
      Y2=DLOG(YY(3))
      Y3=DLOG(YY(4))
      A0=(X-X1)*(X-X2)*(X-X3)/(X0-X1)/(X0-X2)/(X0-X3)
      A1=(X-X0)*(X-X2)*(X-X3)/(X1-X0)/(X1-X2)/(X1-X3)
      A2=(X-X0)*(X-X1)*(X-X3)/(X2-X0)/(X2-X1)/(X2-X3)
      A3=(X-X0)*(X-X1)*(X-X2)/(X3-X0)/(X3-X1)/(X3-X2)
      FINT_HDEC=DEXP(A0*Y0+A1*Y1+A2*Y2+A3*Y3)
      RETURN
      END

      DOUBLE PRECISION FUNCTION SP_HDEC(X)
C--REAL DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 CX,LI2_HDEC
      CX = DCMPLX(X,0.D0)
      SP_HDEC = DREAL(LI2_HDEC(CX))
      RETURN
      END
 
      COMPLEX*16 FUNCTION LI2_HDEC(X)
C--COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Y,CLI2_HDEC
      COMMON/CONST_HDEC/ZETA2,ZETA3
      ZERO=1.D-16
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      LI2_HDEC=0
      IF(R2.LE.ZERO)THEN
        LI2_HDEC=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1.D0.AND.XI.EQ.0.D0)THEN
        IF(XR.EQ.1.D0)THEN
          LI2_HDEC=DCMPLX(ZETA2)
        ELSE
          LI2_HDEC=-DCMPLX(ZETA2/2.D0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.GT.0.5D0)THEN
        Y=(X-1.D0)/X
        LI2_HDEC=CLI2_HDEC(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
     .          +0.5D0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1.D0.AND.RR.LE.0.5D0)THEN
        Y=1.D0/X
        LI2_HDEC=-CLI2_HDEC(Y)-ZETA2-0.5D0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.GT.0.5D0)THEN
        Y=1.D0-X
        LI2_HDEC=-CLI2_HDEC(Y)+ZETA2-CDLOG(X)*CDLOG(1.D0-X)
       RETURN
      ELSEIF(R2.LE.1.D0.AND.XR.LE.0.5D0)THEN
        Y=X
        LI2_HDEC=CLI2_HDEC(Y)
        RETURN
      ENDIF
      END
 
      COMPLEX*16 FUNCTION CLI2_HDEC(X)
C--TAYLOR-EXPANSION FOR COMPLEX DILOGARITHM (SPENCE-FUNCTION)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 X,Z
      COMMON/BERNOULLI_HDEC/B2(18),B12(18),B3(18)
      COMMON/POLY_HDEC/NBER
      N=NBER-1
      Z=-CDLOG(1.D0-X)
      CLI2_HDEC=B2(NBER)
      DO 111 I=N,1,-1
        CLI2_HDEC=Z*CLI2_HDEC+B2(I)
111   CONTINUE
      CLI2_HDEC=Z**2*CLI2_HDEC+Z
      RETURN
      END
 
      DOUBLE PRECISION FUNCTION FACTRL_HDEC(N)
C--DOUBLE PRECISION VERSION OF FACTORIAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FACTRL_HDEC=1.D0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        FACTRL_HDEC=FACTRL_HDEC*DFLOAT(I)
999   CONTINUE
      RETURN
      END
 
      SUBROUTINE BERNINI_HDEC(N)
C--INITIALIZATION OF COEFFICIENTS FOR POLYLOGARITHMS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(18),PB(19)
      COMMON/BERNOULLI_HDEC/B2(18),B12(18),B3(18)
      COMMON/CONST_HDEC/ZETA2,ZETA3
      COMMON/POLY_HDEC/NBER
 
      NBER=N
      PI=4.D0*DATAN(1.D0)
 
      B(1)=-1.D0/2.D0
      B(2)=1.D0/6.D0
      B(3)=0.D0
      B(4)=-1.D0/30.D0
      B(5)=0.D0
      B(6)=1.D0/42.D0
      B(7)=0.D0
      B(8)=-1.D0/30.D0
      B(9)=0.D0
      B(10)=5.D0/66.D0
      B(11)=0.D0
      B(12)=-691.D0/2730.D0
      B(13)=0.D0
      B(14)=7.D0/6.D0
      B(15)=0.D0
      B(16)=-3617.D0/510.D0
      B(17)=0.D0
      B(18)=43867.D0/798.D0
      ZETA2=PI**2/6.D0
      ZETA3=1.202056903159594D0
 
      DO 995 I=1,18
        B2(I)=B(I)/FACTRL_HDEC(I+1)
        B12(I)=DFLOAT(I+1)/FACTRL_HDEC(I+2)*B(I)/2.D0
        PB(I+1)=B(I)
        B3(I)=0.D0
995   CONTINUE
      PB(1)=1.D0
      DO 996 I=1,18
      DO 996 J=0,I
        B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACTRL_HDEC(I-J)/FACTRL_HDEC(J+1)
     .                                            /DFLOAT(I+1)
996   CONTINUE
 
      RETURN
      END

      DOUBLE PRECISION FUNCTION QQINT_HDEC(RAT,H1,H2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      N = 2
      QQINT_HDEC = RAT**N * H1 + (1-RAT**N) * H2
      RETURN
      END

      DOUBLE PRECISION FUNCTION XITLA_HDEC(NO,ALP,ACC)
C--ITERATION ROUTINE TO DETERMINE IMPROVED LAMBDAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      B0(NF)=33.D0-2.D0*NF
      B1(NF)=6.D0*(153.D0-19.D0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12.D0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .              *(1.D0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .              /DLOG(X**2/XLB**2))
      AA(NF)=12D0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2.D0*(1D0+DSQRT(1D0-4D0*B*DLOG(X)))
      PI=4.D0*DATAN(1.D0)
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
111   XITLA_HDEC=XLB
      RETURN
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DOUBLE PRECISION FUNCTION COFSUSY_HDEC(IHIGGS,AMB,RMB,QQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMPLEX*16 C03_HDEC
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HMASS_HDEC/AMSM,AMA,AMHL,AMHH,AMCH,AMAR
      COMMON/GLUINO_HDEC/AMG,AMSB1,AMSB2,STH,CTH,
     .              GLBB(2,2),GHBB(2,2),GABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              GLTT(2,2),GHTT(2,2),GATT(2,2)
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      FC1(VI,VJ,AI,AJ,BI,BJ,BIJ,CIJ,AMI,AMJ) = -1.D0/4*
     .  ((VI*VJ+AI*AJ)/(AMH**2-4*AMB**2)*(BI+BJ-2*BIJ
     .                  +(AMI**2+AMJ**2-2*AMG**2-2*AMB**2)*CIJ)
     .  + AMG/RMB*(VI*VJ-AI*AJ)*CIJ)
      FCA(VI,VJ,AI,AJ,BI,BJ,BIJ,CIJ,AMI,AMJ) = -1.D0/4*(
     .   (VI*AJ+AI*VJ)/AMH**2*(BJ-BI+(AMI**2-AMJ**2)*CIJ)
     .  + AMG/RMB*(AI*VJ-VI*AJ)*CIJ
     .  )
      FC2(VI,AI,A0I,A0G,BI,BPI,AMI) = -1.D0/8*(
     .   (VI**2+AI**2)/2/AMB**2*(A0G-A0I+(AMB**2-AMG**2+AMI**2)*BI
     .                          +2*(AMB**2+AMG**2-AMI**2)*AMB**2*BPI)
     .  + 2*AMG*AMB*(VI**2-AI**2)*BPI
     .  )
      FC3(VI,AI,A0I,A0G,BI,AMI,Q2) = 1.D0/8*(
     .   (VI**2+AI**2)/2/Q2*(A0I-A0G+(Q2+AMG**2-AMI**2)*BI)
     .   +AMG/RMB*(VI**2-AI**2)*BI
     .  )
      FC4(VI,AI,A0I,A0G,BI,BPI,AMI) = 1.D0/8*(
     .   (VI**2+AI**2)/2*((A0G-A0I)/(AMG**2-AMI**2)
     .   +(AMG**2-AMI**2)*BPI)
     .   +AMG/RMB*(VI**2-AI**2)*BI
     .  )
      CF = 4.D0/3.D0
      PI = 4*DATAN(1.D0)
c     write(6,*)
      IF(IHIGGS.EQ.1)THEN
       AMH = AMHL
       GLO = GLB/AMZ**2
       G11 = GLBB(1,1)/GLO
       G12 = GLBB(1,2)/GLO
       G21 = GLBB(2,1)/GLO
       G22 = GLBB(2,2)/GLO
c      write(6,*)'h:'
      ELSEIF(IHIGGS.EQ.2)THEN
       AMH = AMHH
       GLO = GHB/AMZ**2
       G11 = GHBB(1,1)/GLO
       G12 = GHBB(1,2)/GLO
       G21 = GHBB(2,1)/GLO
       G22 = GHBB(2,2)/GLO
c      write(6,*)'H:'
      ELSEIF(IHIGGS.EQ.3)THEN
       AMH = AMA
       GLO = GAB/AMZ**2
       G11 = GABB(1,1)/GLO
       G12 = GABB(1,2)/GLO
       G21 = GABB(2,1)/GLO
       G22 = GABB(2,2)/GLO
c      write(6,*)'A:'
      ENDIF
c     write(6,*)'=='
c     write(6,*)'MB, MSB1, MSB2, MG: ',AMB,AMSB1,AMSB2,AMG
c     write(6,*)'LO,11,12,21,22: ',
c    .          GLO*AMZ**2,G11,G12,G21,G22
c     write(6,*)'SIN(THETA), COS(THETA): ',STH,CTH
      XMU = AMH
      V1 = CTH-STH
      V2 = -CTH-STH
      A1 = CTH+STH
      A2 = CTH-STH
      CC11  = DREAL(C03_HDEC(AMB**2,AMB**2,AMH**2,AMSB1,AMG,AMSB1))
      CC12  = DREAL(C03_HDEC(AMB**2,AMB**2,AMH**2,AMSB1,AMG,AMSB2))
      CC21  = DREAL(C03_HDEC(AMB**2,AMB**2,AMH**2,AMSB2,AMG,AMSB1))
      CC22  = DREAL(C03_HDEC(AMB**2,AMB**2,AMH**2,AMSB2,AMG,AMSB2))
      BB1 = B02_HDEC(AMB**2,AMG,AMSB1,XMU**2)
      BB2 = B02_HDEC(AMB**2,AMG,AMSB2,XMU**2)
      BB11 = B02_HDEC(AMH**2,AMSB1,AMSB1,XMU**2)
      BB12 = B02_HDEC(AMH**2,AMSB1,AMSB2,XMU**2)
      BB21 = B02_HDEC(AMH**2,AMSB2,AMSB1,XMU**2)
      BB22 = B02_HDEC(AMH**2,AMSB2,AMSB2,XMU**2)
      BP1 = BP02_HDEC(AMB**2,AMG,AMSB1,XMU**2)
      BP2 = BP02_HDEC(AMB**2,AMG,AMSB2,XMU**2)
      AA1 = AMSB1**2*(1+DLOG(XMU**2/AMSB1**2))
      AA2 = AMSB2**2*(1+DLOG(XMU**2/AMSB2**2))
      AAG = AMG**2*(1+DLOG(XMU**2/AMG**2))
      BCT1 = B02_HDEC(QQ**2,AMG,AMSB1,XMU**2)
      BCT2 = B02_HDEC(QQ**2,AMG,AMSB2,XMU**2)
      BPCT1 = BP02_HDEC(QQ**2,AMG,AMSB1,XMU**2)
      BPCT2 = BP02_HDEC(QQ**2,AMG,AMSB2,XMU**2)
c     write(6,*)'A0: m1, m2, mg: ',AA1,AA2,AAG
c     write(6,*)'B0: g1, g2, 11, 12, 21, 22: ',BB1,BB2,BB11,BB12,BB21,BB22
c     write(6,*)'B''0: g1, g2: ',BP1,BP2
c     write(6,*)'B''0: g1, g2: ',BPCT1,BPCT2
c     write(6,*)'C0: 11, 12, 21, 22: ',CC11,CC12,CC21,CC22
      IF(IHIGGS.EQ.3)THEN
       COF1 = G11*FCA(V1,V1,A1,A1,BB1,BB1,BB11,CC11,AMSB1,AMSB1)
     .      + G12*FCA(V1,V2,A1,A2,BB1,BB2,BB12,CC12,AMSB1,AMSB2)
     .      + G21*FCA(V2,V1,A2,A1,BB2,BB1,BB21,CC21,AMSB2,AMSB1)
     .      + G22*FCA(V2,V2,A2,A2,BB2,BB2,BB22,CC22,AMSB2,AMSB2)
      ELSE
       COF1 = G11*FC1(V1,V1,A1,A1,BB1,BB1,BB11,CC11,AMSB1,AMSB1)
     .      + G12*FC1(V1,V2,A1,A2,BB1,BB2,BB12,CC12,AMSB1,AMSB2)
     .      + G21*FC1(V2,V1,A2,A1,BB2,BB1,BB21,CC21,AMSB2,AMSB1)
     .      + G22*FC1(V2,V2,A2,A2,BB2,BB2,BB22,CC22,AMSB2,AMSB2)
      ENDIF
      COF2 = FC2(V1,A1,AA1,AAG,BB1,BP1,AMSB1)
     .     + FC2(V2,A2,AA2,AAG,BB2,BP2,AMSB2)
      IF(QQ.EQ.0.D0)THEN
       COF3 = FC4(V1,A1,AA1,AAG,BCT1,BPCT1,AMSB1)
     .      + FC4(V2,A2,AA2,AAG,BCT2,BPCT2,AMSB2)
      ELSE
       COF3 = FC3(V1,A1,AA1,AAG,BCT1,AMSB1,QQ**2)
     .      + FC3(V2,A2,AA2,AAG,BCT2,AMSB2,QQ**2)
      ENDIF
      COF1 = 2*CF*COF1
      COF2 = 2*CF*COF2
      COF3 = 2*CF*COF3
      COFSUSY_HDEC = COF1 + COF2 + COF3
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     COFSUSY_HDEC = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     write(6,*)cof1,cof2+cof3
      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE BOTSUSY_HDEC(GLB,GHB,GAB,XGLB,XGHB,XGAB,SCALE,IL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ICASE = 0
      CALL DMBAPP_HDEC(ICASE,DGLB,DGHB,DGAB,SCALE,IL)
      XGLB = GLB*(1+DGLB)
      XGHB = GHB*(1+DGHB)
      XGAB = GAB*(1+DGAB)
      RETURN
      END
 
      SUBROUTINE DMBAPP_HDEC(ICASE,DGLB,DGHB,DGAB,SCALE,IL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AMCHAR(2),AMNEUT(4),XMNEUT(4),
     .          XMST(2),XMSB(2),AMSL(2),
     .          AMSU(2),AMSD(2),AMSE(2),AMSN(2)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/MASSES_HDEC/AMS,AMC,AMB,AMT
      COMMON/HMASS_HDEC/AMSM,AMA,AMHL,AMHH,AMCH,AMAR
      COMMON/GLUINO_HDEC/AMG,AMSB1,AMSB2,STH,CTH,
     .              GLBB(2,2),GHBB(2,2),GABB(2,2),
     .              AMST1,AMST2,STHT,CTHT,
     .              GLTT(2,2),GHTT(2,2),GATT(2,2)
      COMMON/COUP_HDEC/GAT,GAB,GLT,GLB,GHT,GHB,GZAH,GZAL,
     .            GHHH,GLLL,GHLL,GLHH,GHAA,GLAA,GLVV,GHVV,
     .            GLPM,GHPM,B,A
      COMMON/BREAK_HDEC/AMEL,AMER,AMSQ,AMUR,AMDR,AL,AU,AD,AMU,AM2
      COMMON/SMASS_HDEC/AMNEUT,XMNEUT,AMCHAR,XMST,XMSB,AMSL,
     .              AMSU,AMSD,AMSE,AMSN 
      PI = 4*DATAN(1.D0)
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      TANB = DTAN(B)
      TANA = DTAN(A)
      SB = TANB/DSQRT(1+TANB**2)
      AT = AU
      AB = AD
      SCALELW = SCALE*(AMST1+AMST2+AMU)/3
      SCALQCD = SCALE*(AMSB1+AMSB2+AMG)/3
      RMTOP   = RUNM_HDEC(SCALELW,6)
      HT = RMTOP/V/SB
      STOP1 = AMST1
      STOP2 = AMST2
      SBOT1 = AMSB1
      SBOT2 = AMSB2

      FELW = 1
      FQCD = 1
      IF(IL.EQ.2)THEN
       ASH = ALPHAS_HDEC(SCALELW,2)
       CELW = FELW_HDEC(SCALELW,AMU,AMG,SBOT1,SBOT2,STOP1,STOP2,AMT)
       FELW = 1+ASH/PI*CELW
       ASH = ALPHAS_HDEC(SCALQCD,2)
       CQCD = FQCD_HDEC(SCALQCD,AMT,AMG,SBOT1,SBOT2,STOP1,STOP2,
     .                  AMSU(1),AMSU(2),AMSD(1),AMSD(2))
       FQCD = 1+ASH/PI*CQCD
      ENDIF
 
      ASH = ALPHAS_HDEC(SCALQCD,2)
      IF(ICASE.EQ.0)THEN
       DELTAMB = 2*ASH/3/PI*AMG*AMU*TANB*T_HDEC(SBOT1,SBOT2,AMG)*FQCD
     *         /(1-2*ASH/3/PI*AMG*AB*T_HDEC(SBOT1,SBOT2,AMG))
     * + HT**2/(4*PI)**2*AT*AMU*TANB*T_HDEC(STOP1,STOP2,AMU)*FELW
       DGLB = -DELTAMB/(1+DELTAMB)*(1+1/TANA/TANB)
       DGHB = -DELTAMB/(1+DELTAMB)*(1-TANA/TANB)
       DGAB = -DELTAMB/(1+DELTAMB)*(1+1/TANB**2)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      write(6,*)deltamb,fqcd,felw
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      sc1 = 0.1d0
c      sc2 = 10
c      nsc = 101
c      open(77,file='fort.77')
c      open(78,file='fort.78')
c      do i = 1,nsc
c       scfac = sc1*(sc2/sc1)**((i-1)/(nsc-1.d0))
c       scalelw = scfac*(amst1+amst2+amu)/3
c       scalqcd = scfac*(amsb1+amsb2+amg)/3

c       ash = alphas_hdec(scalqcd,2)
c       dmbqcd = 2*ash/3/pi*amg*amu*tanb*t_hdec(sbot1,sbot2,amg)
c       cqcd = fqcd_hdec(scalqcd,amt,amg,sbot1,sbot2,stop1,stop2,
c    .                   amsu(1),amsu(2),amsd(1),amsd(2))
c       fqcd = 1+ash/pi*cqcd
c       dmbqcd1 = dmbqcd*fqcd
c       rmtop   = runm_hdec(scalelw,6)
c       ht = rmtop/v/sb
c       dmbelw = ht**2/(4*pi)**2*at*amu*tanb*t_hdec(stop1,stop2,amu)
c       ash = alphas_hdec(scalelw,2)
c       celw = felw_hdec(scalelw,amu,amg,sbot1,sbot2,stop1,stop2,amt)
c       felw = 1+ash/pi*celw
c       dmbelw1 = dmbelw*felw
c       write(77,*)scfac,dmbqcd,dmbqcd1
c       write(78,*)scfac,dmbelw,dmbelw1
c      enddo
c      close(77)
c      close(78)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ELSE
       DELTAMB = 2*ASH/3/PI*AMG*AMU*TANB*T_HDEC(SBOT1,SBOT2,AMG)
       DGLB = -DELTAMB*(1+1/TANA/TANB)
       DGHB = -DELTAMB*(1-TANA/TANB)
       DGAB = -DELTAMB*(1+1/TANB**2)
      ENDIF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     DGLB = 0
c     DGHB = 0
c     DGAB = 0
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     write(6,*)'delta_b: ',deltamb, ash

      RETURN
      END
 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

***********************************************************************
        FUNCTION ETA_HDEC(C1,C2)
***********************************************************************
*       COMPLEX ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       8.06.90    ANSGAR DENNER                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA_HDEC,C1,C2
        REAL*8     PI,IM1,IM2,IM12                                     
                                                                       
        PI     = 4D0*DATAN(1D0)                                        
        IM1    = DIMAG(C1)                                             
        IM2    = DIMAG(C2)                                             
        IM12   = DIMAG(C1*C2)                                          
                                                                       
        IF(IM1.LT.0D0.AND.IM2.LT.0D0.AND.IM12.GT.0D0) THEN             
            ETA_HDEC = DCMPLX(0D0,2D0*PI)
        ELSE IF (IM1.GT.0D0.AND.IM2.GT.0D0.AND.IM12.LT.0D0) THEN       
            ETA_HDEC = DCMPLX(0D0,-2D0*PI)
        ELSE                                                           
            ETA_HDEC = DCMPLX(0D0)
        END IF                                                         
        END                                                            

***********************************************************************
        FUNCTION ETAS_HDEC(Y,R,RS)
***********************************************************************
*       MODIFIED ETA-FUNKTION                                           
*---------------------------------------------------------------------*
*       18.1.94   SD                                       
***********************************************************************
        IMPLICIT   LOGICAL(A-Z)                                        
        COMPLEX*16 ETA_HDEC,ETAS_HDEC,Y,R,RS
        REAL*8     PI,IMY,IMRS
                                                                       
        PI     = 4D0*DATAN(1D0)                                        

	IF( DIMAG(R).NE.0D0 ) THEN
	    ETAS_HDEC = ETA_HDEC(Y,R)
	ELSE	    
	    IF( DREAL(R).GT.0D0 ) THEN
		ETAS_HDEC = DCMPLX(0D0,0D0)
	    ELSE
	 	IMY  = DIMAG(Y)
		IMRS = DIMAG(RS)
		ETAS_HDEC = 2D0*DCMPLX(0D0,PI)*(
     *			(1D0+SIGN(1D0,-IMY))*(1D0+SIGN(1D0,-IMRS))-
     *			(1D0+SIGN(1D0, IMY))*(1D0+SIGN(1D0, IMRS))
     *					  )/4D0
	    ENDIF
	ENDIF
        END                                                            

***********************************************************************
        FUNCTION SQE_HDEC(A,B,C)
***********************************************************************
*       SOLUTION OF QUADRATIC EQUATION				      *
*---------------------------------------------------------------------*
*       13.1.92  SD						      *
***********************************************************************
        IMPLICIT REAL*8 (A-Z)                                        
        COMPLEX*16 A,B,C,SQE_HDEC,X1,X2

	X1=(-B+SQRT(B**2-4D0*A*C))/2D0/A
	X2=(-B-SQRT(B**2-4D0*A*C))/2D0/A

	IF (ABS(X1).GT.ABS(X2)) THEN
	   SQE_HDEC=X1
	ELSE
	   SQE_HDEC=X2
	ENDIF

        END                                                            

************************************************************************
        FUNCTION D04_HDEC(P1,P2,P3,P4,P12,P23,M1,M2,M3,M4)
************************************************************************
*  SCALAR 4-POINT FUNCTION WITH AT LEAST ONE MASS ZERO                 *
*  P1,P2,P3,P4 = SQUARED EXTERNAL MOMENTA			       *
*  P12 = (p1+p2)**2,  P23 = (p2+p3)**2				       *
*----------------------------------------------------------------------*
*  2.1.92  SD	         					       *
************************************************************************
        IMPLICIT REAL*8 (A-Z)
	REAL*8 M(4),P(4,4),K(4,4)
	COMPLEX*16 A1,A2,A3,A4,SWAP
	COMPLEX*16 SS(4), XX(2), X(2,4),RS(4,4)
	COMPLEX*16 S0(4),XX0(2),X0(2,4), R(4,4),G(2)
        COMPLEX*16 C04,D04_HDEC,CSPEN_HDEC,ETA_HDEC,SQE_HDEC,ETAS_HDEC
	COMPLEX*16 AA,BB,CC,DD,IEPS,H,HH,L1,L2,L3,L4
        COMPLEX*16 Z2,B,SC,TC,WP,WM,BS,XS
	INTEGER GEN,I,J

        PI = 4*DATAN(1.D0)

        MM1=M1
        MM2=M2
        MM3=M3
        MM4=M4
        M12=M1*M1
        M22=M2*M2
        M32=M3*M3
        M42=M4*M4
        Q1=P1
        Q2=P2
        Q3=P3
	Q4=P4
        Q12=P12
        Q23=P23

C	IS AT LEAST ONE MASS ZERO ???
	IF (MM1*MM2*MM3*MM4.NE.0D0) GOTO 130

C	PERMUTATE UNTIL MM3=0D0
	GOTO 20
10	CONTINUE
	MM0=MM1
	MM1=MM2
	MM2=MM3
	MM3=MM4
	MM4=MM0
	M02=M12
	M12=M22
	M22=M32
	M32=M42
	M42=M02
	Q00=Q12
	Q12=Q23
	Q23=Q00
	Q0=Q1
	Q1=Q2
	Q2=Q3
	Q3=Q4
	Q4=Q0
20	IF (MM3.NE.0D0) GOTO 10
C	ONLY MM3 IS ZERO
	IF (MM1*MM2*MM4.NE.0D0) GOTO 30
C	ONLY MM3 AND MM4 ARE ZERO ==> 3->2, 4->3...
	IF ((MM1*MM2.NE.0D0).AND.(MM4.EQ.0D0)) GOTO 10
C	ONLY MM2 AND MM3 ARE ZERO
	IF ((MM1*MM4.NE.0D0).AND.(MM2.EQ.0D0)) GOTO 40
	WRITE(*,*)'CASE OF THIS SPECIAL D0-FUNCTION NOT IMPLEMENTED!'
	STOP

C	****** NO MASS EQUAL TO ZERO ******
130	CONTINUE
	EPS=1D-18
	IEPS=DCMPLX(0D0,EPS)

	IF( ABS((MM1**2+MM3**2-Q12)/MM1/MM3).LT.2D0 ) THEN
C	R13 WOULD BE NOT REAL. -> PERMUTATION! -> R(2,4) IS NOT REAL.
	   M(1)=MM2
	   M(2)=MM3
	   M(3)=MM4
	   M(4)=MM1
	   P(1,2)=Q2
	   P(1,3)=Q23
	   P(1,4)=Q1
	   P(2,3)=Q3
	   P(2,4)=Q12
	   P(3,4)=Q4
	ELSE
C	R(1,3) IS REAL.
	   M(1)=MM1
	   M(2)=MM2
	   M(3)=MM3
	   M(4)=MM4
	   P(1,2)=Q1
	   P(1,3)=Q12
	   P(1,4)=Q4
	   P(2,3)=Q2
	   P(2,4)=Q23
	   P(3,4)=Q3
	ENDIF

	DO 11 J=2,4
	DO 11 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	R(I,J) =SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
11	CONTINUE

	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	S0(1)=R(1,2)
	S0(2)=R(2,3)
	S0(3)=R(3,4)
	S0(4)=R(1,4)
	AA=K(3,4)/R(2,4)+R(1,3)*K(1,2)-K(1,4)*R(1,3)/R(2,4)-K(2,3)
	BB=(R(2,4)-1D0/R(2,4))*(R(1,3)-1D0/R(1,3))
     *		+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)/R(1,3)+R(2,4)*K(3,4)-K(1,4)*R(2,4)/R(1,3)-K(2,3)
	DD=K(2,3)-R(1,3)*K(1,2)-R(2,4)*K(3,4)+R(1,3)*R(2,4)*K(1,4)
	XX(1)=SQE_HDEC(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	XX0(1)=SQE_HDEC(AA,BB,CC)
	XX0(2)=CC/AA/XX0(1)
c	IF (ABS(DREAL(XX0(1)-XX(2))).LT.ABS(DREAL(XX0(1)-XX(1)))) THEN
	IF (ABS(XX0(1)-XX(2)).LT.ABS(XX0(1)-XX(1))) THEN
	  SWAP  =XX0(1)
	  XX0(1)=XX0(2)
	  XX0(2)=SWAP
	ENDIF

	DO 12 I=1,2
	G(I)  =SIGN( 1D0,DREAL(AA*(XX(I)-XX(3-I))) )
	 X(I,1)= XX(I)/R(2,4)
	X0(I,1)=XX0(I)/R(2,4)
	 X(I,2)= XX(I)/R(2,4)*R(1,3)
	X0(I,2)=XX0(I)/R(2,4)*R(1,3)
	 X(I,3)= XX(I)*R(1,3)
	X0(I,3)=XX0(I)*R(1,3)
	 X(I,4)= XX(I)
	X0(I,4)=XX0(I)
12	CONTINUE

	D04_HDEC = DCMPLX(0D0,0D0)
	DO 13 I=1,2
	DO 13 J=1,4
	A1 = 1D0+X0(I,J)*S0(J) + ABS(1D0+X0(I,J)*S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)*SS(J)))
	A2 = 1D0+X0(I,J)/S0(J) + ABS(1D0+X0(I,J)/S0(J))*IEPS*
     *				  SIGN(1D0,DIMAG(X(I,J)/SS(J)))
	D04_HDEC = D04_HDEC + (-1D0)**(I+J)*(
     *		CSPEN_HDEC(A1)+ETA_HDEC(-X(I,J),SS(J))*LOG(A1)
     *	       +CSPEN_HDEC(A2)+ETA_HDEC(-X(I,J),1D0/SS(J))*LOG(A2))
13	CONTINUE

	IF( DIMAG(R(1,3)).EQ.0D0 ) THEN
	DO 14 I=1,2
	   A1 = (K(1,3)-2D0*R(1,3))/XX0(I)
     *		      -R(1,3)*K(1,4)+K(3,4)
     	   A2 = ((K(2,4)-2D0*R(2,4))*R(1,3)*XX0(I)
     *		      -R(2,4)*K(3,4)+K(2,3))/DD
	   A3 = (K(1,3)-2D0*R(1,3))*R(2,4)/XX0(I)
     *		      -R(1,3)*K(1,2)+K(2,3)
	   A4 = ((K(2,4)-2D0*R(2,4))*XX0(I)
     *		      -R(2,4)*K(1,4)+K(1,2))/DD
	   L1 = LOG( A1-ABS(A1)*IEPS )
     	   L2 = LOG( A2+ABS(A2)*IEPS*G(I)*SIGN(1D0,DREAL(R(1,3))
     *				        	  *DIMAG(RS(2,4))) )
	   L3 = LOG( A3-ABS(A3)*IEPS )
	   L4 = LOG( A4+ABS(A4)*IEPS*G(I)*SIGN(1D0,DIMAG(RS(2,4))) )

	   D04_HDEC = D04_HDEC + (3D0-2D0*I)*(
     *		 ETAS_HDEC(-XX(I),R(1,3),RS(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L1 + L2 )
     *		+ETAS_HDEC(-XX(I),1D0/R(2,4),1D0/RS(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L3 + L4 )
     *		-( ETAS_HDEC(-XX(I),R(1,3)/R(2,4),RS(1,3)/RS(2,4))
     *		  +ETA_HDEC(RS(1,3),1D0/RS(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 + L2 )
     *	  	+ETA_HDEC(RS(1,3),1D0/RS(2,4))
     *		  *ETAS_HDEC(-XX(I),-R(1,3)/R(2,4),-RS(1,3)/RS(2,4))   )
14	CONTINUE
	ELSE
	DO 15 I=1,2
	   L1 = LOG( R(2,4)/XX0(I)+XX0(I)/R(2,4)+K(1,2)
     *		     -XX0(I)/R(2,4)*EPS*BB*G(I) )
	   L2 = LOG( R(1,3)*XX0(I)+1D0/XX0(I)/R(1,3)+K(3,4)
     *		     -XX0(I)*R(1,3)*EPS*BB*G(I) )
	   L3 = LOG( R(1,3)/R(2,4)*XX0(I)+R(2,4)/XX0(I)/R(1,3)+K(2,3)
     *		     -XX0(I)*R(1,3)/R(2,4)*EPS*BB*G(I) )

	   D04_HDEC = D04_HDEC + (3D0-2D0*I)*(
     *		+ETA_HDEC(-XX(I),1D0/R(2,4))
     *		   *( LOG(XX(I)/R(2,4)) + L1 )
     *		+ETA_HDEC(-XX(I),R(1,3))
     *		   *( LOG(R(1,3)*XX(I)) + L2 )
     *		-( ETA_HDEC(-XX(I),R(1,3)/R(2,4))
     *		  +ETA_HDEC(R(1,3),1D0/R(2,4)) )
     *		   *( LOG(XX(I)*R(1,3)/R(2,4)) + L3 )
     *	  	+ETA_HDEC(R(1,3),1D0/R(2,4))
     *		   *ETA_HDEC(-XX(I),-R(1,3)/R(2,4))
     *		   *(1D0-G(I)*SIGN(1D0,DREAL(BB)))	    )
15	CONTINUE
	ENDIF

	D04_HDEC = D04_HDEC/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN


C--->	***************** SPEZIELL ( --> T.SACK-PROMOTION )
C	D1=Q12-M12
C	D2=Q2 -M22
C	D3=Q3 -M42
C	IF ((D1*D2.LE.0D0).OR.(D2*D3.LE.0D0)) THEN
C	   WRITE(*,*) 'THE CASE OF DIFFERENT SIGNS OF THE D1,D2,D3'
C	   WRITE(*,*) 'IN D04(...) IS NOT IMPLEMENTED !!!'
C	   STOP
C	ENDIF
C	NM1=ABS(MM1/D1)
C	NM2=ABS(MM2/D2)
C	NM3=ABS(MM4/D3)
C	NP1=Q2/D2**2+Q12/D1**2+(Q1-Q2-Q12)/D1/D2
C	NP2=Q2/D2**2+ Q3/D3**2+(Q23-Q2-Q3)/D2/D3
C	NP3=Q3/D3**2+Q12/D1**2+(Q4-Q3-Q12)/D1/D3
C	D04_HDEC=C04(NP1,NP2,NP3,NM1,NM2,NM3)/D1/D2/D3

C	*************** ALLGEMEIN


C	****** ONLY MM3 IS ZERO ******
30	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)
	M(1)=MM1
	M(2)=MM2
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 1 J=2,4
	DO 1 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
1	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(3,4)/R(2,4)-K(2,3)
	BB=K(1,3)*(1D0/R(2,4)-R(2,4))+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(1,3)*K(1,4)*R(2,4)+R(2,4)*K(3,4)-K(2,3)
	DD=K(2,3)-R(2,4)*K(3,4)
	XX(1)=SQE_HDEC(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 2 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
2	CONTINUE
	D04_HDEC = DCMPLX(0D0,0D0)
	DO 3 I=1,2
	D04_HDEC = D04_HDEC + (2D0*I-3D0)*(
     *		CSPEN_HDEC(1D0+SS(4)*X(I,4))
     *	       -CSPEN_HDEC(1D0+SS(1)*X(I,1))
     *	       +CSPEN_HDEC(1D0+X(I,4)/SS(4))
     *	       -CSPEN_HDEC(1D0+X(I,1)/SS(1))
     *	       +ETA_HDEC(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       -ETA_HDEC(-X(I,1),SS(1))*LOG(1D0+SS(1)*X(I,1))
     *	       +ETA_HDEC(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -ETA_HDEC(-X(I,1),1D0/SS(1))*LOG(1D0+X(I,1)/SS(1))
     *	       -CSPEN_HDEC(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +CSPEN_HDEC(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA_HDEC(-X(I,4),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,4)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       +ETA_HDEC(-X(I,1),(K(2,3)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+X(I,1)*(K(2,3)-IEPS)/(K(1,3)-IEPS))   )
	IF (DIMAG(R(2,4)).NE.0D0) THEN
	   H=ETA_HDEC(-1D0/XX(I),R(2,4))
	ELSE
	   H=DCMPLX(0D0,0D0)
	   IF (DREAL(R(2,4)).LT.0D0) THEN
	      HH=-1D0/XX(I)
	      IM1=DIMAG(HH)
	      IM2=DIMAG(RS(2,4))
	      IF ((IM1.GT.0D0).AND.(IM2.GT.0D0)) THEN
	         H=-DCMPLX(0D0,2D0*PI)
	      ENDIF
	      IF ((IM1.LT.0D0).AND.(IM2.LT.0D0)) THEN
	         H=+DCMPLX(0D0,2D0*PI)
	      ENDIF
	   ENDIF
	ENDIF
	D04_HDEC = D04_HDEC + (2D0*I-3D0)*
     *	          H*( LOG( (K(1,2)-R(2,4)*K(1,4)
     *			  +XX(I)*(1D0/R(2,4)-R(2,4)))/DD )
     *		     +LOG(K(1,3)-IEPS) )
3	CONTINUE
	D04_HDEC = D04_HDEC/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))
	RETURN

C	****** ONLY MM2 AND MM3 ARE ZERO ******
40	CONTINUE
	EPS=1D-17
	IEPS=DCMPLX(0D0,EPS)

	M(1)=MM1
	M(2)=10D0
	M(3)=10D0
	M(4)=MM4
	P(1,2)=Q1
	P(1,3)=Q12
	P(1,4)=Q4
	P(2,3)=Q2
	P(2,4)=Q23
	P(3,4)=Q3
	DO 4 J=2,4
	DO 4 I=1,J-1
	K(I,J)=(M(I)**2+M(J)**2-P(I,J))/M(I)/M(J)
	IF (I.EQ.2) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.2) K(I,J)=K(I,J)-M(J)/M(I)
	IF (I.EQ.3) K(I,J)=K(I,J)-M(I)/M(J)
	IF (J.EQ.3) K(I,J)=K(I,J)-M(J)/M(I)
	R(I,J) =SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),0D0),
     *	            DCMPLX(1D0,0D0))
	IF( DIMAG(R(I,J)).EQ.0D0 ) THEN
	   RS(I,J)=SQE_HDEC(DCMPLX(1D0,0D0),DCMPLX(-K(I,J),EPS),
     *	               DCMPLX(1D0,0D0))
	ELSE
	   RS(I,J)=R(I,J)
	ENDIF
4	CONTINUE
	SS(1)=RS(1,2)
	SS(2)=RS(2,3)
	SS(3)=RS(3,4)
	SS(4)=RS(1,4)
	AA=K(2,4)*K(3,4)-K(2,3)
	BB=K(1,3)*K(2,4)+K(1,2)*K(3,4)-K(1,4)*K(2,3)
	CC=K(1,2)*K(1,3)-K(2,3)
	DD=K(2,3)
	XX(1)=SQE_HDEC(AA,BB,CC+IEPS*DD)
	XX(2)=(CC+IEPS*DD)/AA/XX(1)
	DO 5 I=1,2
	X(I,1)=XX(I)/R(2,4)
	X(I,2)=XX(I)/R(2,4)*R(1,3)
	X(I,3)=XX(I)*R(1,3)
	X(I,4)=XX(I)
5	CONTINUE
	D04_HDEC = DCMPLX(0D0,0D0)
	DO 6 I=1,2
	D04_HDEC = D04_HDEC + (2D0*I-3D0)*(
     *		CSPEN_HDEC(1D0+SS(4)*X(I,4))
     *	       +CSPEN_HDEC(1D0+X(I,4)/SS(4))
     *	       +ETA_HDEC(-X(I,4),SS(4))*LOG(1D0+SS(4)*X(I,4))
     *	       +ETA_HDEC(-X(I,4),1D0/SS(4))*LOG(1D0+X(I,4)/SS(4))
     *	       -CSPEN_HDEC(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -CSPEN_HDEC(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       -ETA_HDEC(-XX(I),(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(3,4)-IEPS)/(K(1,3)-IEPS))
     *	       -ETA_HDEC(-XX(I),(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	           *LOG(1D0+XX(I)*(K(2,4)-IEPS)/(K(1,2)-IEPS))
     *	       +LOG(-XX(I))*( LOG(K(1,2)-IEPS)
     *			     +LOG(K(1,3)-IEPS)-LOG(K(2,3)-IEPS) ) )
6	CONTINUE
	D04_HDEC = D04_HDEC/M(1)/M(2)/M(3)/M(4)/AA/(XX(1)-XX(2))

	RETURN

	END

************************************************************************
        FUNCTION C03_HDEC(P1,P2,P3,M1,M2,M3)
************************************************************************
*  SCALAR 3-POINT FUNCTION                                             *
*  P1,P2,P3 = SQUARED EXTERNAL MOMENTA  			       *
*----------------------------------------------------------------------*
*  5.12.96  M. SPIRA    					       *
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3
      REAL*8 R(0:2)
      COMPLEX*16 C03_HDEC,CSPEN_HDEC,ETA_HDEC,IEPS,IM
      COMPLEX*16 ALP(0:2),X(0:2,2),Y0(0:2),Y(0:2,2)
      COMPLEX*16 CDUM
C     REAL*8 KAPPA
      COMPLEX*16 KAPPA
C     KAPPA(A,B,C) = DSQRT(A**2+B**2+C**2-2*(A*B+A*C+B*C))
C     KAPPA(A,B,C) = DSQRT(DABS(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
c     KAPPA(A,B,C) = CDSQRT(DCMPLX(A**2+B**2+C**2-2*(A*B+A*C+B*C)))
      KAPPA(A,B,C,D) = CDSQRT((A**2+B**2+C**2-2*(A*B+A*C+B*C))
     .               * (1+IEPS*D))
      EPS = 1.D-8*(P1+P2+P3)
      IM = DCMPLX(0.D0,1.D0)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IEPS = DCMPLX(0.D0,1.D-17)
c     IEPS = DCMPLX(0.D0,1.D-20)
      PI = 4*DATAN(1.D0)
      XX = 0.D0
C     IF(P1.LT.0.D0.OR.P2.LT.0.D0.OR.P3.LT.0.D0) XX=1.D0
      IF(P1.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q10 = P1
      ELSE
       Q10 = EPS
      ENDIF
      IF(P3.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q20 = P3
      ELSE
       Q20 = EPS
      ENDIF
      IF(P2.NE.0.D0.OR.XX.NE.0.D0)THEN
       Q21 = P2
      ELSE
       Q21 = EPS
      ENDIF
      R(0) = P2
      R(1) = P3
      R(2) = P1
      SM0 = M1**2
      SM1 = M2**2
      SM2 = M3**2
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     ALPHA  = KAPPA(Q10,Q21,Q20)
c     ALP(0) = KAPPA(Q21,SM1,SM2)*(1+IEPS*Q21)
c     ALP(1) = KAPPA(Q20,SM2,SM0)*(1+IEPS*Q20)
c     ALP(2) = KAPPA(Q10,SM0,SM1)*(1+IEPS*Q10)
      ALPHA  = KAPPA(Q10,Q21,Q20,1.D0)
      ALP(0) = KAPPA(Q21,SM1,SM2,DSIGN(1.D0,Q21))
      ALP(1) = KAPPA(Q20,SM2,SM0,DSIGN(1.D0,Q20))
      ALP(2) = KAPPA(Q10,SM0,SM1,DSIGN(1.D0,Q10))
      X(0,1) = (Q21 - SM1 + SM2 + ALP(0))/2/Q21
      X(0,2) = (Q21 - SM1 + SM2 - ALP(0))/2/Q21
      X(1,1) = (Q20 - SM2 + SM0 + ALP(1))/2/Q20
      X(1,2) = (Q20 - SM2 + SM0 - ALP(1))/2/Q20
      X(2,1) = (Q10 - SM0 + SM1 + ALP(2))/2/Q10
      X(2,2) = (Q10 - SM0 + SM1 - ALP(2))/2/Q10
      Y0(0) = (Q21*(Q21-Q20-Q10+2*SM0-SM1-SM2) - (Q20-Q10)*(SM1-SM2)
     .      + ALPHA*(Q21-SM1+SM2))/2/ALPHA/Q21
      Y0(1) = (Q20*(Q20-Q10-Q21+2*SM1-SM2-SM0) - (Q10-Q21)*(SM2-SM0)
     .      + ALPHA*(Q20-SM2+SM0))/2/ALPHA/Q20
      Y0(2) = (Q10*(Q10-Q21-Q20+2*SM2-SM0-SM1) - (Q21-Q20)*(SM0-SM1)
     .      + ALPHA*(Q10-SM0+SM1))/2/ALPHA/Q10
      Y(0,1) = Y0(0) - X(0,1)
      Y(0,2) = Y0(0) - X(0,2)
      Y(1,1) = Y0(1) - X(1,1)
      Y(1,2) = Y0(1) - X(1,2)
      Y(2,1) = Y0(2) - X(2,1)
      Y(2,2) = Y0(2) - X(2,2)
      CDUM=0.D0
      DO I=0,2
       DO J=1,2
        CDUM = CDUM + CSPEN_HDEC((Y0(I)-1)/Y(I,J))
     .              - CSPEN_HDEC(Y0(I)/Y(I,J))
        CX = ETA_HDEC(1-X(I,J),1/Y(I,J))
        IF(CX.NE.0.D0)THEN
         CDUM = CDUM + CX*CDLOG((Y0(I)-1)/Y(I,J))
        ENDIF
        CY = ETA_HDEC(-X(I,J),1/Y(I,J))
        IF(CY.NE.0.D0)THEN
         CDUM = CDUM - CY*CDLOG(Y0(I)/Y(I,J))
        ENDIF
       ENDDO
       CX = ETA_HDEC(-X(I,1),-X(I,2))
       IF(CX.NE.0.D0)THEN
        CDUM = CDUM - CX*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       CY = ETA_HDEC(Y(I,1),Y(I,2))
       IF(CY.NE.0.D0)THEN
        CDUM = CDUM + CY*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
       A = -R(I)
       B = -DIMAG(Y(I,1)*Y(I,2))
       IF(A.GT.0.D0.AND.B.GT.0.D0) THEN
        CDUM = CDUM + 2*PI*IM*CDLOG((1-Y0(I))/(-Y0(I)))
       ENDIF
      ENDDO
      C03_HDEC = CDUM/ALPHA
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C        SUBROUTINE CALCULATING THE FINITE REAL PART OF THE           C
C          GENERAL MASSIVE TWO POINT FUNCTION                         C
C                                                                     C
C           B02(P.P,M1,M2,MU**2)                                      C
C           BP02(P.P,M1,M2,MU**2)                                     C
C                                                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c ---------------------------------------------------------------------
      real*8 function B02_HDEC(s,m1,m2,mu2)

      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2 

      m12 = m1**2 
      m22 = m2**2 

      zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     &                     -2.D0*(s*m12+s*m22+m12*m22)))

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            B02_HDEC=-dlog(m12/mu2) 
         else
            B02_HDEC=1.D0 - m12/(m12-m22)*dlog(m12/mu2)
     &                 + m22/(m12-m22)*dlog(m22/mu2) 
         endif
      else 
         if ((m12.eq.0.D0).and.(m22.eq.0.D0)) then 
            B02_HDEC=2.D0 - dlog(s/mu2)
         elseif ((m12.eq.s).and.(m22.eq.0.D0)) then 
            B02_HDEC=2.D0 - dlog(m12/mu2)
         elseif ((m22.eq.s).and.(m12.eq.0.D0)) then 
            B02_HDEC=2.D0 - dlog(m22/mu2)
         elseif (m12.eq.0.D0) then
            B02_HDEC=2.D0 - (s-m22)/s*dlog( dabs(m22-s)/m22 )
     &                 - dlog(m22/mu2)
         elseif (m22.eq.0.D0) then
            B02_HDEC=2.D0 - (s-m12)/s*dlog( dabs(m12-s)/m12 ) 
     &                 - dlog(m12/mu2)
         else
            x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
            x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
            B02_HDEC=dreal( 2.D0+ dlog(mu2/m22) 
     &                       + x1*cdlog(1.D0-1.D0/x1) 
     &                       + x2*cdlog(1.D0-1.D0/x2))
         endif
      endif 

      return
      end



c ---------------------------------------------------------------------
      real*8 function BP02_HDEC(s,m1,m2,mu2)
      
      implicit none 

      real*8     s,m1,m2,mu2,m12,m22 
      complex*16 zkappa,x1,x2
      
      m12 = m1**2
      m22 = m2**2 

      zkappa=cdsqrt(dcmplx(s**2+m12**2+m22**2
     &                    -2.D0*(s*m12+s*m22+m12*m22)))

      if (s.eq.0.D0) then
         if (m12.eq.m22) then
            BP02_HDEC=1.D0/(6.D0*m12)
         else
            BP02_HDEC=( (m12+m22)/2.D0 
     &        - m12*m22/(m12-m22)*dlog(m12/m22) )/(m12-m22)**2 
         endif
      elseif ((s.eq.m12).and.(m22.eq.0.D0)) then 
         BP02_HDEC=( -1.D0 + dlog(m12/mu2)/2.D0 )/m12
      elseif ((s.eq.m22).and.(m12.eq.0.D0)) then 
         BP02_HDEC=( -1.D0 + dlog(m22/mu2)/2.D0 )/m22
      else 
         x1=dcmplx( (s-m22+m12+zkappa)/(2.D0*s) )
         x2=dcmplx( (s-m22+m12-zkappa)/(2.D0*s) )
         BP02_HDEC=dreal( -1.D0 + ( x1*(1.D0-x1)*cdlog(1.D0-1.D0/x1)
     &                     - x2*(1.D0-x2)*cdlog(1.D0-1.D0/x2) )  
     &                                                  /(x1-x2) )/s
      endif 

      return
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION CSPEN_HDEC(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C       20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX*16 CSPEN_HDEC,W,SUM,Z,U
        REAL*8 RZ,AZ,A1
        REAL*8 B(9)/
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
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
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
c      write(*,*) 'z:',z
      Z =Z*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c      write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
        CSPEN_HDEC=-CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen_HDEC
        RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
        CSPEN_HDEC=1.64493406684822643D0
c        write(*,*) 'cspen:', cspen_HDEC
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
 2    CSPEN_HDEC=SUM
c        write(*,*) 'cspen:', cspen_HDEC
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
12    CSPEN_HDEC=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
c        write(*,*) 'cspen:', cspen_HDEC
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
22    CSPEN_HDEC=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
c        write(*,*) 'cspen:', cspen_HDEC
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
32    CSPEN_HDEC=SUM+3.28986813369645287D0
     *               +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
50    CONTINUE
c        write(*,*) 'cspen:', cspen_HDEC
      END

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

      SUBROUTINE SUBH1_HDEC(MA,TANB,MQ,MUR,MD,MTOP,AU,AD,MU,MCHI0,
     *                 MHP,HMP,MCH,SA,CA,TANBA,MGLU)

      IMPLICIT REAL*8(A-H,L,M,O-Z)
      DIMENSION VH(2,2),M2(2,2),M2P(2,2)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HSELF_HDEC/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7

      MCHI = MCHI0
      TANBA = TANB
      TANBT = TANB
      
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS_HDEC(AMZ,2)
      ALPHA3  = ALPHAS_HDEC(MTOP,2)
      MB      = RUNM_HDEC(MTOP,5)
      RMTOP   = RUNM_HDEC(MTOP,6)

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
      CALL GFUN_HDEC(MA,TANBA,MQ,MUR,MD,MTOP,AU,AD,MU,MGLU,
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
        SUBROUTINE GFUN_HDEC(MA,TANB,MQ,MUR,MD,MTOP,AT,AB,MU,MGLU,
     *                     DLAM1,DLAM2,DLAM3,DLAM4,DLAM5,DLAM6,DLAM7)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        DIMENSION VH(2,2),VH1(2,2),VH2(2,2),
     *            VH3T(2,2),VH3B(2,2),AL(2,2)
        COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
        G(X,Y) = 2.D0 - (X+Y)/(X-Y)*DLOG(X/Y)

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      MB = RUNM_HDEC(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS_HDEC(AMZ,2)
      ALPHA3  = ALPHAS_HDEC(MTOP,2)

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

        RMTOP   = RUNM_HDEC(MTOP,6)

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

      CALL DELMB_HDEC(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
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

      FUNCTION T_HDEC(X,Y,Z)
      implicit real*8(a-h,l,m,o-z)
      if(x.eq.y) x = x - 0.00001
      if(x.eq.z) x = x - 0.00002
      if(y.eq.z) y = y - 0.00003
c       write(*,*) 'xyz',x,y,z
      T_HDEC = (X**2*Y**2*log(X**2/Y**2) + X**2*Z**2*log(Z**2/X**2)
     * + Y**2*Z**2*log(Y**2/Z**2))/((X**2-Y**2)*(Y**2-Z**2)*(X**2-Z**2))
      return
      end

      SUBROUTINE DELMB_HDEC(MA,TANB,MQ,MUR,MD,AT,AB,MU,MGLU,
     .           MTOP,DELTAMT,DELTAMB,STOP12,STOP22,SBOT12,SBOT22)
        IMPLICIT REAL*8 (A-H,L,M,O-Z)
        COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW

        IF(DABS(MU).LT.0.000001) MU = 0.000001
        MQ2   = MQ**2
        MUR2  = MUR**2
        MD2   = MD**2
        TANBA = TANB
        SINBA = TANBA/DSQRT(TANBA**2+1.D0)
        COSBA = SINBA/TANBA        
        SINB = TANB/DSQRT(TANB**2+1.D0)
        COSB = SINB/TANB

      RMTOP = RUNM_HDEC(MTOP,6)
      MB = RUNM_HDEC(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      MW = AMW
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS_HDEC(AMZ,2)
      ALPHA3  = ALPHAS_HDEC(MTOP,2)

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
     *  T_HDEC(sbot1,sbot2,mglu)
     *  + ht**2/(4.*pi)**2*(at-mu/tanb)*mu*tanb*
     *  T_HDEC(stop1,stop2,mu)


        deltamt = -2.*alpha3/3./pi*(at-mu/tanb)*mglu*
     *  T_HDEC(stop1,stop2,mglu)

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

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     THIS PROGRAM COMPUTES THE RENORMALIZATION GROUP IMPROVED
C     VALUES OF HIGGS MASSES AND COUPLINGS IN THE MSSM.
C
C     INPUT: MA,TANB = TAN(BETA),MQ,MUR,MTOP,AU,AD,MU.
C
C     ALL MASSES IN GEV UNITS. MA IS THE CP-ODD HIGGS MASS,
C     MTOP IS THE PHYSICAL TOP MASS, MQ AND MUR ARE THE SOFT
C     SUPERSYMMETRY BREAKING MASS PARAMETERS OF LEFT HANDED
C     AND RIGHT HANDED STOPS RESPECTIVELY, AU AND AD ARE THE
C     STOP AND SBOTTOM TRILINEAR SOFT BREAKING TERMS,
C     RESPECTIVELY,  AND MU IS THE SUPERSYMMETRIC
C     HIGGS MASS PARAMETER. WE USE THE  CONVENTIONS FROM
C     THE PHYSICS REPORT OF HABER AND KANE: LEFT RIGHT
C     STOP MIXING TERM PROPORTIONAL TO (AU - MU/TANB).
C
C     WE USE AS INPUT TANB DEFINED AT THE SCALE MTOP.
C
C     OUTPUT: MH,HM,MHCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C
C     WHERE MH AND HM ARE THE LIGHTEST AND HEAVIEST CP-EVEN
C     HIGGS MASSES, MHCH IS THE CHARGED HIGGS MASS AND
C     ALPHA IS THE HIGGS MIXING ANGLE.
C
C     TANBA IS THE ANGLE TANB AT THE CP-ODD HIGGS MASS SCALE.
C
C     RANGE OF VALIDITY:
C
C    (STOP1**2 - STOP2**2)/(STOP2**2 + STOP1**2) < 0.5
C    (SBOT1**2 - SBOT2**2)/(SBOT2**2 + SBOT2**2) < 0.5
C
C     WHERE STOP1, STOP2, SBOT1 AND SBOT2 ARE THE STOP AND
C     ARE THE SBOTTOM  MASS EIGENVALUES, RESPECTIVELY. THIS
C     RANGE AUTOMATICALLY EXCLUDES THE EXISTENCE OF TACHYONS.
C
C
C     FOR THE CHARGED HIGGS MASS COMPUTATION, THE METHOD IS
C     VALID IF
C
C     2 * |MB * AD* TANB|  < M_SUSY**2,  2 * |MTOP * AU| < M_SUSY**2
C
C     2 * |MB * MU * TANB| < M_SUSY**2,  2 * |MTOP * MU| < M_SUSY**2
C
C     WHERE M_SUSY**2 IS THE AVERAGE OF THE SQUARED STOP MASS
C     EIGENVALUES, M_SUSY**2 = (STOP1**2 + STOP2**2)/2. THE SBOTTOM
C     MASSES HAVE BEEN ASSUMED TO BE OF ORDER OF THE STOP ONES.
C
C     M_SUSY**2 = (MQ**2 + MUR**2)*0.5 + MTOP**2
C
C     PROGRAM BASED ON THE WORK BY M. CARENA, J.R. ESPINOSA,
C     M. QUIROS AND C.E.M. WAGNER, PHYS. LETT. B355 (1995) 209
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SUBH2_HDEC(MA,TANB,MQ,MUR,MTOP,AU,AD,MU,MH,HM,
     * MHCH,SA,CA,TANBA)
      IMPLICIT REAL*8(A-H,L,M,O-Z)
      COMMON/PARAM_HDEC/GF,ALPH,AMTAU,AMMUON,AMZ,AMW
      COMMON/HSELF_HDEC/LAMBDA1,LAMBDA2,LAMBDA3,LAMBDA4,LAMBDA5,
     .             LAMBDA6,LAMBDA7
C     MZ = 91.18
C     ALPHA1 = 0.0101
C     ALPHA2 = 0.0337
C     ALPHA3Z = 0.12
C     V = 174.1
C     PI = 3.14159
      TANBA = TANB
      TANBT = TANB

C     MBOTTOM(MTOP) = 3. GEV
C     MB = 3.
C     ALPHA3 = ALPHA3Z/(1. +(11. - 10./3.)/4./PI*ALPHA3Z*
C    *LOG(MTOP**2/MZ**2))

C     RMTOP= RUNNING TOP QUARK MASS
C     RMTOP = MTOP/(1.+4.*ALPHA3/3./PI)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MB = RUNM_HDEC(MTOP,5)
      PI = 4*DATAN(1D0)
      MZ = AMZ
      V  = 1/DSQRT(2*DSQRT(2D0)*GF)
      CW = AMW**2/AMZ**2
      SW = 1-CW
      ALPHA2  = (2*AMW/V/DSQRT(2D0))**2/4/PI
      ALPHA1  = ALPHA2*SW/CW
      ALPHA3Z = ALPHAS_HDEC(AMZ,2)
      ALPHA3  = ALPHAS_HDEC(MTOP,2)
      RMTOP   = RUNM_HDEC(MTOP,6)
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C      RMTOP=MTOP
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MS = ((MQ**2 + MUR**2)/2. + MTOP**2)**.5
      T = LOG(MS**2/MTOP**2)
      SINB = TANB/((1. + TANB**2)**.5)
      COSB = SINB/TANB
C      IF(MA.LE.MTOP) TANBA = TANBT
      IF(MA.GT.MTOP)
     *TANBA = TANBT*(1.-3./32./PI**2*
     *(RMTOP**2/V**2/SINB**2-MB**2/V**2/COSB**2)*
     *LOG(MA**2/MTOP**2))

      SINBT = TANBT/((1. + TANBT**2)**.5)
      COSBT = 1./((1. + TANBT**2)**.5)
      COS2BT = (TANBT**2 - 1.)/(TANBT**2 + 1.)
      G1 = (ALPHA1*4.*PI)**.5
      G2 = (ALPHA2*4.*PI)**.5
      G3 = (ALPHA3*4.*PI)**.5
      HU = RMTOP/V/SINBT
      HD =  MB/V/COSBT

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C      G3=0
C      HU=0
C      HD=0
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      XAU = (2.*AU**2/MS**2)*(1. - AU**2/12./MS**2)
      XAD = (2.*AD**2/MS**2)*(1. - AD**2/12./MS**2)
      AUD = (-6.*MU**2/MS**2 - ( MU**2- AD*AU)**2/MS**4.
     *+ 3.*(AU + AD)**2/MS**2)/6.
      LAMBDA1 = ((G1**2 + G2**2)/4.)*(1.-3.*HD**2*T/8./PI**2)
     *+(3.*HD**4/8./PI**2) * (T + XAD/2. + (3.*HD**2/2. + HU**2/2.
     *- 8.*G3**2) * (XAD*T + T**2)/16./PI**2)
     *-(3.*HU**4* MU**4/96./PI**2/MS**4) * (1+ (9.*HU**2 -5.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA2 = ((G1**2 + G2**2)/4.)*(1.-3.*HU**2*T/8./PI**2)
     *+(3.*HU**4/8./PI**2) * (T + XAU/2. + (3.*HU**2/2. + HD**2/2.
     *- 8.*G3**2) * (XAU*T + T**2)/16./PI**2)
     *-(3.*HD**4* MU**4/96./PI**2/MS**4) * (1+ (9.*HD**2 -5.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA3 = ((G2**2 - G1**2)/4.)*(1.-3.*
     *(HU**2 + HD**2)*T/16./PI**2)
     *+(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2)
     *+(3.*HU**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4)* (1.+ (6.*HU**2 -2.* HD**2/2.
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4)*(1.+ (6.*HD**2 -2.* HU**2
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA4 = (- G2**2/2.)*(1.-3.*(HU**2 + HD**2)*T/16./PI**2)
     *-(6.*HU**2*HD**2/16./PI**2) * (T + AUD/2. + (HU**2 + HD**2
     *- 8.*G3**2) * (AUD*T + T**2)/16./PI**2)
     *+(3.*HU**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AU**2/
     *MS**4)*
     *(1+ (6.*HU**2 -2.* HD**2
     *-  16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4/96./PI**2) * (3.*MU**2/MS**2 - MU**2*AD**2/
     *MS**4)*
     *(1+ (6.*HD**2 -2.* HU**2/2.
     *-  16.*G3**2) *T/16./PI**2)
      LAMBDA5 = -(3.*HU**4* MU**2*AU**2/96./PI**2/MS**4) *
     * (1- (2.*HD**2 -6.* HU**2 + 16.*G3**2) *T/16./PI**2)
     *-(3.*HD**4* MU**2*AD**2/96./PI**2/MS**4) *
     * (1- (2.*HU**2 -6.* HD**2 + 16.*G3**2) *T/16./PI**2)
      LAMBDA6 = (3.*HU**4* MU**3*AU/96./PI**2/MS**4) *
     * (1- (7.*HD**2/2. -15.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HD**4* MU *(AD**3/MS**3 - 6.*AD/MS )/96./PI**2/MS) *
     * (1- (HU**2/2. -9.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
      LAMBDA7 = (3.*HD**4* MU**3*AD/96./PI**2/MS**4) *
     * (1- (7.*HU**2/2. -15.* HD**2/2. + 16.*G3**2) *T/16./PI**2)
     *+(3.*HU**4* MU *(AU**3/MS**3 - 6.*AU/MS )/96./PI**2/MS) *
     * (1- (HD**2/2. -9.* HU**2/2. + 16.*G3**2) *T/16./PI**2)
      TRM2 = MA**2 + 2.*V**2* (LAMBDA1* COSBT**2 +
     *2.* LAMBDA6*SINBT*COSBT
     *+ LAMBDA5*SINBT**2 + LAMBDA2* SINBT**2 + 2.* LAMBDA7*SINBT*COSBT
     *+ LAMBDA5*COSBT**2)
      DETM2 = 4.*V**4*(-(SINBT*COSBT*(LAMBDA3 + LAMBDA4) +
     *LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2)**2 + (LAMBDA1* COSBT**2 +
     *2.* LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2)*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2)) + MA**2*2.*V**2 *
     *((LAMBDA1* COSBT**2 +2.*
     *LAMBDA6* COSBT*SINBT + LAMBDA5*SINBT**2)*COSBT**2 +
     *(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT + LAMBDA5*COSBT**2)
     **SINBT**2
     * +2.*SINBT*COSBT* (SINBT*COSBT*(LAMBDA3
     * + LAMBDA4) + LAMBDA6*COSBT**2
     *+ LAMBDA7* SINBT**2))

      MH2 = (TRM2 - (TRM2**2 - 4.* DETM2)**.5)/2.
      HM2 = (TRM2 + (TRM2**2 - 4.* DETM2)**.5)/2.
      HM = HM2**.5
      MH = MH2**.5
      MHCH2 = MA**2 + (LAMBDA5 - LAMBDA4)* V**2
      MHCH = MHCH2**.5
      MHCH = MHCH2**.5

      SINALPHA = SQRT(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.*
     *LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2)))/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5))/2.**.5

      COSALPHA = (2.*(2.*V**2*(SINBT*COSBT*(LAMBDA3 + LAMBDA4) +
     *LAMBDA6*COSBT**2 + LAMBDA7* SINBT**2) -
     *MA**2*SINBT*COSBT))/2.**.5/
     *SQRT(((TRM2**2 - 4.* DETM2)**.5)*
     *(((TRM2**2 - 4.* DETM2)**.5) -
     * ((2.*V**2*(LAMBDA1* COSBT**2 + 2.*
     *LAMBDA6* COSBT*SINBT
     *+ LAMBDA5*SINBT**2) + MA**2*SINBT**2)
     *- (2.*V**2*(LAMBDA2* SINBT**2 +2.* LAMBDA7* COSBT*SINBT
     *+ LAMBDA5*COSBT**2) + MA**2*COSBT**2))))

      SA = -SINALPHA
      CA = -COSALPHA

 2242 RETURN
      END

c change susyhit
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Comment by Maggie: For use of hdecay.f within SUSYHIT, I put the file 
c dmb.f, which is separate in the program package HDECAY, into the file
c hdecay.f in order to reduce the number of separate files within the 
c program package SUSYHIT.
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function t134p_hdec(am1,am2,mu)
      implicit double precision (a-h,o-z)
      double precision m1,m2,mu,ll1,ll2
      complex*16 li2_hdec
      sp(a) = dreal(li2_hdec(dcmplx(a,0.d0)))
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6
      if(am1.lt.am2)then
       m1 = am1
       m2 = am2
      else
       m1 = am2
       m2 = am1
      endif
      ll1 = dlog(mu**2/m1**2)
      ll2 = dlog(mu**2/m2**2)
      if(m1.eq.m2)then
       dummy = 7*(m1**2+m2**2)/2
     .       + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .       - m1**2/2*dlog(m1**2/m2**2)**2
      else
       dummy = 7*(m1**2+m2**2)/2
     .       + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .       + (m1**2-m2**2)*(dlog(m1**2/m2**2)*dlog(1-m1**2/m2**2)
     .                       + sp(m1**2/m2**2)-zeta2)
     .       - m1**2/2*dlog(m1**2/m2**2)**2
      endif
      t134p_hdec = dummy/mu**2
      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function t134_hdec(am1,am2,am3,mu)
      implicit double precision (a-h,o-z)
      double precision m1,m2,m3,mu,ll1,ll2,ll3
      complex*16 sp,li2_hdec,lam,ll,cx,cy,phi
      sp(cx) = li2_hdec(cx)
      lam(cx,cy) = (1-cx-cy)**2-4*cx*cy
      phi(cx,cy,ll) = (2*cdlog((1+cx-cy-ll)/2)*cdlog((1-cx+cy-ll)/2)
     .              -cdlog(cx)*cdlog(cy)+2*zeta2
     .              -2*sp((1+cx-cy-ll)/2)-2*sp((1-cx+cy-ll)/2))/ll
      eps = 1.d-15
      rim = dcmplx(1.d0,eps)
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6
      m1 = dmin1(am1,am2,am3)
      m3 = dmax1(am1,am2,am3)
      if(m1.eq.am2.and.m3.eq.am3.or.m1.eq.am3.and.m3.eq.am2) m2 = am1
      if(m1.eq.am1.and.m3.eq.am3.or.m1.eq.am3.and.m3.eq.am1) m2 = am2
      if(m1.eq.am1.and.m3.eq.am2.or.m1.eq.am2.and.m3.eq.am1) m2 = am3
      cx = m1**2/m3**2*rim
      cy = m2**2/m3**2*rim
      ll1 = dlog(mu**2/m1**2)
      ll2 = dlog(mu**2/m2**2)
      ll3 = dlog(mu**2/m3**2)
      ll = cdsqrt(lam(cx,cy))
      dummy = 7*(m1**2+m2**2+m3**2)/2
     .      + m1**2*(ll1**2+3*ll1) + m2**2*(ll2**2+3*ll2)
     .      + m3**2*(ll3**2+3*ll3)
     .      +  (m1**2-m2**2-m3**2)/4*dlog(m2**2/m3**2)**2
     .      + (-m1**2+m2**2-m3**2)/4*dlog(m1**2/m3**2)**2
     .      + (-m1**2-m2**2+m3**2)/4*dlog(m1**2/m2**2)**2
     .      + m3**2/2*lam(cx,cy)*phi(cx,cy,ll)
      t134_hdec = dummy/mu**2
      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function
     .       felw_hdec(scale,amu,amg,amsb1,amsb2,amst1,amst2,amt)
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amg,amsb1,amsb2,amst1,amst2,amu,amt
      double precision mg,mb1,mb2,mt1,mt2,mu,mt
      double precision anomalous
      double precision a
      complex*16 sp,li2_hdec,xgl,xt1,xt2,xt,xb1,xb2
      sp(r) = li2_hdec(r)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      fnorm = 1/amu**2

      cf = 4/3.d0

      rim = dcmplx(1.d0,eps)

      mt  = amt
      mu  = dabs(amu)
      mg  = amg
      mt1 = amst1
      mt2 = amst2
      mb1 = amsb1
      mb2 = amsb2

      xt  = amt**2/amu**2   * rim
      xgl = amg**2/amu**2   * rim
      xt1 = amst1**2/amu**2 * rim
      xt2 = amst2**2/amu**2 * rim
      xb1 = amsb1**2/amu**2 * rim
      xb2 = amsb2**2/amu**2 * rim

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      ans3=t134(mt2,mg,mt,mg)*xgl*(-xb1*xgl*xt1+xb1*xgl-xb1*xt*xt1+
     . xb1*xt+xb1*xt1*xt2-xb1*xt2-xb2*xgl*xt1+xb2*xgl-xb2*xt*xt1+xb2*
     . xt+xb2*xt1*xt2-xb2*xt2+2*xgl**2*xt1-2*xgl**2+2*xgl*xt*xt1-2*
     . xgl*xt-2*xgl*xt1*xt2+2*xgl*xt2)+t134(mu,mb1,mt,mg)*xgl*(-xb1*
     . xb2*xt1+xb1*xb2*xt2+xb1*xgl*xt1-xb1*xgl*xt2-xb2*xt*xt1+xb2*xt*
     . xt2+xb2*xt1-xb2*xt2+xgl*xt*xt1-xgl*xt*xt2-xgl*xt1+xgl*xt2)+
     . t134(mu,mb2,mt,mg)*xgl*(-xb1*xb2*xt1+xb1*xb2*xt2-xb1*xt*xt1+
     . xb1*xt*xt2+xb1*xt1-xb1*xt2+xb2*xgl*xt1-xb2*xgl*xt2+xgl*xt*xt1-
     . xgl*xt*xt2-xgl*xt1+xgl*xt2)+t134(mu,mg,mt,mg)*xgl*(xb1*xgl*xt1
     . -xb1*xgl*xt2+xb1*xt*xt1-xb1*xt*xt2-xb1*xt1+xb1*xt2+xb2*xgl*xt1
     . -xb2*xgl*xt2+xb2*xt*xt1-xb2*xt*xt2-xb2*xt1+xb2*xt2-2*xgl**2*
     . xt1+2*xgl**2*xt2-2*xgl*xt*xt1+2*xgl*xt*xt2+2*xgl*xt1-2*xgl*xt2
     . )
      ans2=t134(mt1,mb1,mt,mg)*xgl*(-xb1*xb2*xt2+xb1*xb2+xb1*xgl*xt2-
     . xb1*xgl-xb2*xt*xt2+xb2*xt+xb2*xt1*xt2-xb2*xt1+xgl*xt*xt2-xgl*
     . xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mb2,mt,mg)*xgl*(-xb1*xb2*xt2+
     . xb1*xb2-xb1*xt*xt2+xb1*xt+xb1*xt1*xt2-xb1*xt1+xb2*xgl*xt2-xb2*
     . xgl+xgl*xt*xt2-xgl*xt-xgl*xt1*xt2+xgl*xt1)+t134(mt1,mg,mt,mg)*
     . xgl*(xb1*xgl*xt2-xb1*xgl+xb1*xt*xt2-xb1*xt-xb1*xt1*xt2+xb1*xt1
     . +xb2*xgl*xt2-xb2*xgl+xb2*xt*xt2-xb2*xt-xb2*xt1*xt2+xb2*xt1-2*
     . xgl**2*xt2+2*xgl**2-2*xgl*xt*xt2+2*xgl*xt+2*xgl*xt1*xt2-2*xgl*
     . xt1)+t134(mt2,mb1,mt,mg)*xgl*(xb1*xb2*xt1-xb1*xb2-xb1*xgl*xt1+
     . xb1*xgl+xb2*xt*xt1-xb2*xt-xb2*xt1*xt2+xb2*xt2-xgl*xt*xt1+xgl*
     . xt+xgl*xt1*xt2-xgl*xt2)+t134(mt2,mb2,mt,mg)*xgl*(xb1*xb2*xt1-
     . xb1*xb2+xb1*xt*xt1-xb1*xt-xb1*xt1*xt2+xb1*xt2-xb2*xgl*xt1+xb2*
     . xgl-xgl*xt*xt1+xgl*xt+xgl*xt1*xt2-xgl*xt2)+ans3
      ans1=log(xb1)*log(xt1)*xb1*xt1*(-xb2*xt2+xb2+xgl*xt2-xgl)+log(
     . xb1)*log(xt2)*xb1*xt2*(xb2*xt1-xb2-xgl*xt1+xgl)+log(xb2)*log(
     . xt1)*xb2*xt1*(-xb1*xt2+xb1+xgl*xt2-xgl)+log(xb2)*log(xt2)*xb2*
     . xt2*(xb1*xt1-xb1-xgl*xt1+xgl)+log(xgl)*log(xt1)*xt1*(4*xb1*xb2
     . *xt2-4*xb1*xb2-3*xb1*xgl*xt2+3*xb1*xgl-3*xb2*xgl*xt2+3*xb2*xgl
     . +2*xgl**2*xt2-2*xgl**2)+log(xgl)*log(xt2)*xt2*(-4*xb1*xb2*xt1+
     . 4*xb1*xb2+3*xb1*xgl*xt1-3*xb1*xgl+3*xb2*xgl*xt1-3*xb2*xgl-2*
     . xgl**2*xt1+2*xgl**2)+log(xt1)**2*xt1*(-xb1*xb2*xt2+xb1*xb2+xb1
     . *xgl*xt2-xb1*xgl+xb2*xgl*xt2-xb2*xgl-xgl**2*xt2+xgl**2)+4*log(
     . xt1)*xt1*(xb1*xb2*xt2-xb1*xb2-xb1*xgl*xt2+xb1*xgl-xb2*xgl*xt2+
     . xb2*xgl+xgl**2*xt2-xgl**2)+log(xt2)**2*xt2*(xb1*xb2*xt1-xb1*
     . xb2-xb1*xgl*xt1+xb1*xgl-xb2*xgl*xt1+xb2*xgl+xgl**2*xt1-xgl**2)
     . +4*log(xt2)*xt2*(-xb1*xb2*xt1+xb1*xb2+xb1*xgl*xt1-xb1*xgl+xb2*
     . xgl*xt1-xb2*xgl-xgl**2*xt1+xgl**2)+ans2
      r52=ans1/(4*(xb1*xb2*xt1**2*xt2-xb1*xb2*xt1**2-xb1*xb2*xt1*xt2
     . **2+xb1*xb2*xt1+xb1*xb2*xt2**2-xb1*xb2*xt2-xb1*xgl*xt1**2*xt2+
     . xb1*xgl*xt1**2+xb1*xgl*xt1*xt2**2-xb1*xgl*xt1-xb1*xgl*xt2**2+
     . xb1*xgl*xt2-xb2*xgl*xt1**2*xt2+xb2*xgl*xt1**2+xb2*xgl*xt1*xt2
     . **2-xb2*xgl*xt1-xb2*xgl*xt2**2+xb2*xgl*xt2+xgl**2*xt1**2*xt2-
     . xgl**2*xt1**2-xgl**2*xt1*xt2**2+xgl**2*xt1+xgl**2*xt2**2-xgl**
     . 2*xt2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      ans5=2*((2*((14*xt2-15+14*xt1+5*xt)*xt-((xt2-3)*xt2+xt1**2+(4*
     . xt2-3)*xt1))*xgl**3+(4*xt2-3+4*xt1-8*xt)*xgl**4-((2*xt+3)*(xt-
     . xt1)**2*(xt-xt2)**2+2*xgl**5)-(2*(2*((3*xt2-2)*xt2+3*xt1**2)+(
     . 5*xt2-4)*xt1)*xt-((4*(xt2-3)*xt1-3*xt2)*xt2+(4*xt2-3)*xt1**2)-
     . 2*(24*xt2-23+24*xt1)*xt**2-10*xt**3)*xgl**2-2*((2*((3*xt2-2)*
     . xt2+3*xt1**2)+(5*xt2-4)*xt1)*xt**2+(3*((xt2-1)*xt1**2-xt2**2)+
     . (3*xt2+5)*xt1*xt2)*xt-(14*xt2-15+14*xt1)*xt**3+((xt2-3)*xt1-3*
     . xt2)*xt1*xt2+4*xt**4)*xgl)*(xt1-1)*(xt2-1)-(2*(xt+xt1)*xgl-(xt
     . -xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt-1+xgl
     . )*(xt2-2+xt1)*t134(mu,mt,mg,mg)*xgl)*(xt1-xt2)
      ans4=-2*((2*((xt-xt1)**2*(xt-xt2)**2+xgl**4-(6*xt**2+10*xt*xt1+
     . 10*xt*xt2-12*xt-xt1**2-4*xt1*xt2-xt2**2)*xgl**2-2*(xt1+xt2-xt)
     . *xgl**3+2*(((2*xt2-3)*xt2+2*xt1**2+(7*xt2-3)*xt1)*xt-(5*xt2-6+
     . 5*xt1)*xt**2-(xt1+xt2)*xt1*xt2+xt**3)*xgl)*(xt1-1)*(xt2-1)+(2*
     . (xt+xt1)*xgl-(xt-xt1)**2-xgl**2)*(2*(xt+xt2)*xgl-(xt-xt2)**2-
     . xgl**2)*(xt2-2+xt1)*log(xgl))*(xt1-xt2)+(2*(xt+xt1)*xgl+(xt-
     . xt1)**2*(xt1-2)-xgl**2*xt1)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2
     . )*(log(xgl)-log(xt1))*(xt2-1)**2-(2*(xt+xt1)*xgl-(xt-xt1)**2-
     . xgl**2)*(2*(xt+xt2)*xgl+(xt-xt2)**2*(xt2-2)-xgl**2*xt2)*(log(
     . xgl)-log(xt2))*(xt1-1)**2)*(log(xgl)-log(xt))*xt+ans5
      ans3=-(4*((xt-xt2**2+xt2)*(xt-xt2)**2+xgl**3-((xt2+1)*xt2+xt)*
     . xgl**2-((xt+4*xt2**2)*xt-(2*xt2-1)*xt2**2)*xgl)+(((xt+xt2)*(xt
     . -xt2)*(xt2-2)-4*xt*xt2)*xgl+(xt+2*xt2)*(xt2-2)*xgl**2-(xt-xt2)
     . **2*(xt2-2)*xt-(xt2-2)*xgl**3)*(log(xgl)-log(xt2)))*(2*(xt+xt1
     . )*xgl-(xt-xt1)**2-xgl**2)*(log(xgl)-log(xt2))*(xt1-1)**2+(4*((
     . xt-xt1**2+xt1)*(xt-xt1)**2+xgl**3-((xt1+1)*xt1+xt)*xgl**2-((xt
     . +4*xt1**2)*xt-(2*xt1-1)*xt1**2)*xgl)+(((xt+xt1)*(xt-xt1)*(xt1-
     . 2)-4*xt*xt1)*xgl+(xt+2*xt1)*(xt1-2)*xgl**2-(xt-xt1)**2*(xt1-2)
     . *xt-(xt1-2)*xgl**3)*(log(xgl)-log(xt1)))*(2*(xt+xt2)*xgl-(xt-
     . xt2)**2-xgl**2)*(log(xgl)-log(xt1))*(xt2-1)**2+ans4
      ans2=(2*(((xt+2*xt2)*xt+(xt2-4)*xt2)*xgl+(xt2+2+xt)*xgl**2-(xt+
     . xt2-2)*(xt-xt2)**2-xgl**3)*(xt1-1)**2*t134(mt2,mt,mg,mg)*xgl-(
     . 2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(xt2-2+xt1)*(log(xgl)+4)*(
     . xgl+xt)*(xt1-xt2)*log(xgl))*(2*(xt+xt1)*xgl-(xt-xt1)**2-xgl**2
     . )-2*(((xt+2*xt1)*xt+(xt1-4)*xt1)*xgl+(xt1+2+xt)*xgl**2-(xt+xt1
     . -2)*(xt-xt1)**2-xgl**3)*(2*(xt+xt2)*xgl-(xt-xt2)**2-xgl**2)*(
     . xt2-1)**2*t134(mt1,mt,mg,mg)*xgl+((2*(xt1+xt2)-xgl)*xgl**3-(xt
     . -xt1)**2*(xt-xt2)**2+(2*xt**2+6*xt*xt1+6*xt*xt2-8*xt-xt1**2-4*
     . xt1*xt2-xt2**2)*xgl**2-2*(((xt2-2)*xt2+xt1**2+2*(3*xt2-1)*xt1)
     . *xt-(3*xt2-4+3*xt1)*xt**2-(xt1+xt2)*xt1*xt2)*xgl)*(log(xgl)-
     . log(xt))**2*(xt1-xt2)*(xt1-1)*(xt2-1)*xt+ans3
      ans1=-ans2
      r82=ans1/(4*((xt-xt1)**2+xgl**2-2*(xt+xt1)*xgl)*((xt-xt2)**2+
     . xgl**2-2*(xt+xt2)*xgl)*(xt1-xt2)*(xt1-1)**2*(xt2-1)**2)

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      ans6=-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(((xt1-xt2)*xt1+
     . 4*xt**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+4*log(xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2-4*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-
     . xt1)**2+xgl**2)))*xt1+(xt1-xt2)*log(-(xt-xt1+xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1)*(xt-xt1+xgl)*xt2
      ans5=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(((xt1-xt2)*xt2-
     . 4*xt**2)*(log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . )-4*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*
     . xt**2+4*log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . ))*xt**2-(xt1-xt2)*log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2)))*xt2+(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt2)*(xt-xt2+xgl)*xt1+ans6
      ans4=4*(xgl+xt-xt1)*(xgl-xt-xt1)*log(xt-xt1-xgl+sqrt(-2*(xt+xt1
     . )*xgl+(xt-xt1)**2+xgl**2))*xt**2*xt2+4*(xgl+xt-xt1)*(xgl-xt-
     . xt1)*log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt**2*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl
     . +sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-
     . xt2)*(xgl-xt-xt2)*(xt1-xt2)*log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*
     . xgl+(xt-xt2)**2+xgl**2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(
     . xt1-xt2)*log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl
     . **2)))*xt1*xt2+(xgl+xt-xt1)*(xgl-xt-xt1)*(xt1-xt2)*log(-(xt-
     . xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*xt1*xt2+
     . ans5
      ans3=(2*((xt1+xt2-2*xt-2*xgl)*(xt1-xt2)+(xt1-xt2-8*xt)*log(xt2)
     . *xt2+(xt1-xt2+8*xt)*log(xt1)*xt1-6*(xt1-xt2)*log(xt)*xt-(xt1+
     . xt2+2*xt)*(xt1-xt2)*log(xgl))*xt*xt2+(xgl+xt-xt2)*(xgl-xt-xt2)
     . *(4*xt**2-xt1*xt2+xt2**2)*log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+
     . (xt-xt2)**2+xgl**2))+(xgl+xt-xt2)*(xgl-xt-xt2)*(4*xt**2-xt1*
     . xt2+xt2**2)*log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2)))*xt1-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*
     . xt2)*log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*
     . xt2-(xgl+xt-xt1)*(xgl-xt-xt1)*(4*xt**2+xt1**2-xt1*xt2)*log(xt+
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))*xt2-4*(xgl+
     . xt-xt2)*(xgl-xt-xt2)*log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-
     . xt2)**2+xgl**2))*xt**2*xt1-4*(xgl+xt-xt2)*(xgl-xt-xt2)*log(xt-
     . xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))*xt**2*xt1+
     . ans4
      ans7=((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)
      ans2=ans3*ans7
      ans1=-ans2
      rat2=ans1/(16*(xt1-xt2)**2*(xt1-1)*(xt2-1)*xt**2*xt1*xt2)

      ans5=((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2
     . )))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-xt2)*log(xgl)
      ans4=sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(-(xt-xt2+xgl+sqrt(-2*
     . (xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))+log(-(xt-xt2+xgl-sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))))*(xt-xt2+xgl)*((xt1-1)*log(
     . xt2)*xt2-(xt2-1)*log(xt1)*xt1)+2*((log(xgl)-log(xt2))**2*(2*xt
     . -xt2)*(xt1-1)*xt2-2*(xt1-xt2)*log(xgl)**2*xt-(log(xgl)-log(xt1
     . ))**2*(2*xt-xt1)*(xt2-1)*xt1)*xt+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)
     . **2+xgl**2)*(log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+
     . xgl**2))-log(xt+xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**
     . 2))-log(-(xt-xt1+xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))
     . )+log(-(xt-xt1+xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))))
     . *(xt-xt1+xgl)*((xt1-1)*log(xt2)*xt2-(xt2-1)*log(xt1)*xt1)+ans5
      ans3=-ans4
      ans2=-(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl
     . **2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)
     . )+log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt2-1)-2*((log(xgl)-log(xt2))*xt2-log(xgl))*(xt1-xt2)*xt)*(log
     . (xgl)-log(xt1))*xt1+ans3
      ans1=(((log(-(xt-xt2+xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2)))-log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+
     . log(-(xt-xt2+xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)-2*(xt1+xt2-2*xgl)*xt+(log(-(xt-xt1+xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl+sqrt
     . (-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(-(xt-xt1+xgl+sqrt(-2
     . *(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))-log(xt+xt1-xgl-sqrt(-2*(xt
     . +xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl))*(
     . xt1-1)-2*(xt1-xt2)*log(xgl)*xt)*(log(xgl)-log(xt2))*xt2+ans2
      rlt2=ans1/(8*(xt1-xt2)*(xt1-1)*(xt2-1)*xt**2)

      ans3=(2*(((log(xgl)-log(xt2))*(xt1-1)**2*xt2-(xt1-xt2)**2*log(
     . xgl)+(xt2+1-7*xt)*(xt2-1)*xt1)*xt1-(7*xt1**2-3*xt2-2*(xt2+1)*
     . xt1)*(xt2-1)*xgl)*xt1+(2*((2*(xt2+1)*xt1+3*xt2)*xt-2*xt1**3)*
     . xt1-(log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)*(xt1**2-xt2))*(xt2-1))*(log(xgl)-log(xt1))
      ans2=sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)*(log(xt+xt1-xgl-
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+sqrt(-2*(
     . xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt-xt1-xgl-sqrt(-2*(xt+
     . xt1)*xgl+(xt-xt1)**2+xgl**2)))*((xt1**2-xt2)*(xt2-1)*log(xt1)-
     . (xt1-1)**2*log(xt2)*xt2-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt1+xgl
     . )+((log(xt-xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-
     . log(xt+xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(
     . xt-xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+
     . xt1-xgl+sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl
     . )*(xt-xt1+xgl)+2*(5*xt+xt1+5*xgl)*xt1)*((log(xgl)-log(xt2))*(
     . xt1-1)**2*xt2-(xt1-xt2)**2*log(xgl))+ans3
      ans1=2*(((log(xgl)-log(xt2))**2*(xt1-1)**2*xt2-(xt1-xt2)**2*log
     . (xgl)**2)*(xgl+xt)+(xt2+1-5*xt)*(xt2-1)*xt1**2)*xt1-(log(xt-
     . xt1-xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-
     . xgl-sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))+log(xt-xt1-xgl+
     . sqrt(-2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2))-log(xt+xt1-xgl+sqrt(
     . -2*(xt+xt1)*xgl+(xt-xt1)**2+xgl**2)))*(xt+xt1-xgl)*(xt-xt1+xgl
     . )*(xt1-xt2)*(xt1-1)*(xt2-1)-2*(xt+xt1+xgl)*(log(xgl)-log(xt1))
     . **2*(xt1**2-xt2)*(xt2-1)*xt1-2*(5*xt*xt2+xt1**3+5*(xt1-xt2)*(
     . xt1-1)*xgl-(5*(xt2+1)*xt-xt2)*xt1)*(xt2-1)*xt1-2*((log(xgl)-
     . log(xt1))*(xt1**2-xt2)*(xt2-1)-(log(xgl)-log(xt2))*(xt1-1)**2*
     . xt2+((xt1-xt2)*log(xgl)+(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt1+ans2
      rmt12=ans1/(4*(xt1-xt2)**2*(xt1-1)**2*(xt2-1)*xt1)

      ans3=((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))
     . -log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(log(xgl)-log(xt2))*xt2**2+2*(5*xt+xt2+5*xgl)
     . *xt2)*(log(xgl)-log(xt1))*(xt2-1)**2*xt1+(2*(((2*xt2+3)*xt+xt2
     . **2)*xt1**2-(xt1-xt2)**2*log(xgl)*xt2-((7*xt2**2+3)*xt+2*xt2**
     . 3)*xt1)*xt2+(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+
     . xgl**2))-log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**
     . 2))+log(xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+
     . xt2-xgl)*(xt-xt2+xgl)*(xt1-xt2**2)*(xt1-1)+2*(((7*xt-1)*xt2-2*
     . (xt-xt2**2))*xt2-((7*xt2-2)*xt2-(2*xt2+3)*xt1)*(xt1-1)*xgl)*
     . xt2)*(log(xgl)-log(xt2))
      ans2=-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)*(log(xt+xt2-xgl-
     . sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+xt2-xgl+sqrt(
     . -2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(xt-xt2-xgl+sqrt(-2*(
     . xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt-xt2-xgl-sqrt(-2*(xt+
     . xt2)*xgl+(xt-xt2)**2+xgl**2)))*((xt1-xt2**2)*(xt1-1)*log(xt2)+
     . (xt2-1)**2*log(xt1)*xt1-(xt1-xt2)*(xt1-1)*(xt2-1))*(xt-xt2+xgl
     . )-((log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)+2*(5*xt+xt2+5*xgl)*xt2)*(xt1-xt2)**2*log(xgl)+
     . ans3
      ans1=(log(xt-xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-
     . log(xt+xt2-xgl-sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))+log(
     . xt-xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2))-log(xt+
     . xt2-xgl+sqrt(-2*(xt+xt2)*xgl+(xt-xt2)**2+xgl**2)))*(xt+xt2-xgl
     . )*(xt-xt2+xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)+2*((log(xgl)-log(xt1)
     . )**2*(xt2-1)**2*xt1-(xt1-xt2)**2*log(xgl)**2)*(xgl+xt)*xt2+2*(
     . xt+xt2+xgl)*(log(xgl)-log(xt2))**2*(xt1-xt2**2)*(xt1-1)*xt2+2*
     . (5*xt+xt2+5*xgl)*(xt1-xt2)*(xt1-1)*(xt2-1)*xt2+2*((log(xgl)-
     . log(xt1))*(xt2-1)**2*xt1+(log(xgl)-log(xt2))*(xt1-xt2**2)*(xt1
     . -1)-((xt1-xt2)*log(xgl)-(xt1-1)*(xt2-1))*(xt1-xt2))*(log(xgl)-
     . log(xt))*xt*xt2+ans2
      rmt22=ans1/(4*(xt1-xt2)**2*(xt1-1)*(xt2-1)**2*xt2)

      relw = r22+r42+r52+r72+r82+r92+rat2+rlt2+rmt12+rmt22
      bo  = fi(amst1**2,amst2**2,amu**2)
      relw = cf*relw/bo

      anomalous =-cf
      finscale = 2*dlog(scale**2/amg**2)
      felw_hdec = dreal(relw)*fnorm + anomalous + finscale

      return
      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      double precision function fqcd_hdec(scale,amt,amg,amsb1,amsb2,
     .                           amst1,amst2,amsu1,amsu2,amsd1,amsd2)
      implicit double precision (b-h,o-q,s-z), complex*16 (a,r)
      double precision amt,amg,amsb1,amsb2,amst1,amst2,
     .                 amsu1,amsu2,amsd1,amsd2
      double precision mt,mg,mb1,mb2,mt1,mt2,ms1,ms2,mu
      complex*16 sp,li2_hdec,xt,xgl,xt1,xt2,xb1,xb2,xs1,xs2,xq
      double precision m,mq
      double precision anomalous
      double precision a
      complex*16 xp,xm
      sp(r) = li2_hdec(r)
      xp(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . +cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      xm(m) = (amg**2+amt**2-m**2)/2/amg**2/rim
     . -cdsqrt(((amg**2+amt**2-m**2)/2/amg**2/rim)**2-amt**2/amg**2/rim)
      fi(a,b,c) = (a*b*log(a/b)+b*c*log(b/c)+c*a*log(c/a))
     .          / (a-b)/(b-c)/(a-c)
      t134p(a,b,c)  = t134p_hdec(a,b,c)
      t134(a,b,c,d) = t134_hdec(a,b,c,d)

      eps = 1.d-15
      pi = 4*datan(1.d0)
      zeta2 = pi**2/6

      ca = 3
      cf = 4/3.d0
      tr = 1/2.d0
      nu = 2
      nd = 2
      nf = nu+nd+1

      fnorm = 4/amg**2

      rim = dcmplx(1.d0,eps)

      mq  = amt
      mt  = amt
      mg  = amg
      mu  = amg
      mt1 = amsb1
      mt2 = amsb2
      mb1 = amsb1
      mb2 = amsb2

      xq  = amt**2/amg**2   * rim
      xt  = amt**2/amg**2   * rim
      xgl = amg**2/amg**2   * rim
      xt1 = amsb1**2/amg**2 * rim
      xt2 = amsb2**2/amg**2 * rim
      xb1 = amsb1**2/amg**2 * rim
      xb2 = amsb2**2/amg**2 * rim

      r12=(log(xb1)**2*xb1*(-2*xb1*xb2**2+4*xb1*xb2-2*xb1+xb2**2-2*
     . xb2+1)+4*log(xb1)*xb1*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1
     . )+log(xb2)**2*xb2*(2*xb1**2*xb2-xb1**2-4*xb1*xb2+2*xb1+2*xb2-1
     . )+4*log(xb2)*xb2*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1)+4*
     . t134p(mb1,mb1,mg)*xb1*(xb2**2-2*xb2+1)+4*t134p(mb1,mg,mg)*(-
     . xb1*xb2**2+2*xb1*xb2-xb1-xb2**2+2*xb2-1)+4*t134p(mb2,mb2,mg)*
     . xb2*(-xb1**2+2*xb1-1)+4*t134p(mb2,mg,mg)*(xb1**2*xb2+xb1**2-2*
     . xb1*xb2-2*xb1+xb2+1)+4*t134p(mg,mg,mg)*(-xb1**2+2*xb1+xb2**2-2
     . *xb2))/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*
     . xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2
     . **2-xb2))

      r22=(4*log(xgl)**2*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+4*log(xgl)*
     . log(xt1)*xt1*(-xt2**2+2*xt2-1)+4*log(xgl)*log(xt2)*xt2*(xt1**2
     . -2*xt1+1)+12*log(xgl)*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2)+log(xt1
     . )**2*xt1*(xt2**2-2*xt2+1)+2*log(xt1)*xt1*(-xt1*xt2**2+2*xt1*
     . xt2-xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(-xt1**2+2*xt1-1)+2*
     . log(xt2)*xt2*(xt1**2*xt2+2*xt1**2-2*xt1*xt2-4*xt1+xt2+2)+t134p
     . (mu,mt1,mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+
     . t134p(mu,mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)
     . +14*(-xt1**2*xt2+xt1*xt2**2+xt1-xt2))/(2*(xt1**3*xt2**2-2*xt1
     . **3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1**2+2*xt1*xt2**
     . 3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r32=(log(xb1)**2*xb1*(-xb2**2+2*xb2-1)+2*log(xb1)*xb1*(2*xb1*
     . xb2**2-4*xb1*xb2+2*xb1+xb2**2-2*xb2+1)+log(xb2)**2*xb2*(xb1**2
     . -2*xb1+1)+2*log(xb2)*xb2*(-2*xb1**2*xb2-xb1**2+4*xb1*xb2+2*xb1
     . -2*xb2-1)+2*t134p(mb1,mg,mg)*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2
     . **2+2*xb2-1)+2*t134p(mb2,mg,mg)*(-2*xb1**2*xb2+xb1**2+4*xb1*
     . xb2-2*xb1-2*xb2+1)+3*t134p(mg,mg,mg)*(xb1**2*xb2-xb1**2-xb1*
     . xb2**2+xb1+xb2**2-xb2)+14*(xb1**2*xb2-xb1*xb2**2-xb1+xb2))/(8*
     . (xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-
     . 2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      r42=(8*log(xgl)**2*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+8*log(xgl)*
     . log(xt1)*xt1*(xt2**2-2*xt2+1)+8*log(xgl)*log(xt2)*xt2*(-xt1**2
     . +2*xt1-1)+24*log(xgl)*(xt1**2*xt2-xt1*xt2**2-xt1+xt2)+2*log(
     . xt1)**2*xt1*(-xt2**2+2*xt2-1)+log(xt1)*xt1*(5*xt1*xt2**2-10*
     . xt1*xt2+5*xt1+7*xt2**2-14*xt2+7)+2*log(xt2)**2*xt2*(xt1**2-2*
     . xt1+1)+log(xt2)*xt2*(-5*xt1**2*xt2-7*xt1**2+10*xt1*xt2+14*xt1-
     . 5*xt2-7)+2*t134p(mt1,mu,mg)*xgl*(xt1*xt2**2-2*xt1*xt2+xt1+xt2
     . **2-2*xt2+1)+2*t134p(mt2,mu,mg)*xgl*(-xt1**2*xt2-xt1**2+2*xt1*
     . xt2+2*xt1-xt2-1)+28*(xt1**2*xt2-xt1*xt2**2-xt1+xt2))/(2*(xt1**
     . 3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*xt1**2*xt2-2*xt1
     . **2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2**2-xt2))

      r52=(log(xb1)**2*xb1*(2*xb1*xb2**2-4*xb1*xb2+2*xb1-xb2**2+2*xb2
     . -1)+4*log(xb1)*xb1*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+
     . log(xb2)**2*xb2*(-2*xb1**2*xb2+xb1**2+4*xb1*xb2-2*xb1-2*xb2+1)
     . +4*log(xb2)*xb2*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+2*
     . t134p(mb1,mb2,mg)*(xb1**2*xb2-xb1**2-xb1*xb2**2+xb1+xb2**2-xb2
     . )+2*t134p(mb1,mg,mg)*(-xb1**2*xb2+xb1**2+2*xb1*xb2-2*xb1-xb2+1
     . )+2*t134p(mb2,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1-xb2**2+2*xb2-1)
     . )/(16*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3*xb1**
     . 2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*xb2**2-
     . xb2))

      r62=(log(xb1)**2*xb1*(xb2**2-2*xb2+1)+2*log(xb1)*xb1*(-xb2**2+2
     . *xb2-1)+log(xb2)**2*xb2*(-xb1**2+2*xb1-1)+2*log(xb2)*xb2*(xb1
     . **2-2*xb1+1)+2*t134p(mb1,mg,mg)*(xb1*xb2**2-2*xb1*xb2+xb1+xb2
     . **2-2*xb2+1)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2-xb1**2+2*xb1*xb2+
     . 2*xb1-xb2-1)+2*t134p(mg,mg,mg)*(xb1**2*xb2+xb1**2-xb1*xb2**2-3
     . *xb1-xb2**2+3*xb2))/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1
     . **2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1
     . -xb2**3+2*xb2**2-xb2))

      ans1=8*log(xgl)**2*(xt1**2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2
     . )+4*log(xgl)*log(xt1)*xt1*(4*xt1*xt2**2-8*xt1*xt2+4*xt1-3*xt2
     . **2+6*xt2-3)+4*log(xgl)*log(xt2)*xt2*(-4*xt1**2*xt2+3*xt1**2+8
     . *xt1*xt2-6*xt1-4*xt2+3)+12*log(xgl)*(xt1**2*xt2-xt1**2-xt1*xt2
     . **2+xt1+xt2**2-xt2)+log(xt1)**2*xt1*(-8*xt1*xt2**2+16*xt1*xt2-
     . 8*xt1+5*xt2**2-10*xt2+5)+4*log(xt1)*xt1*(3*xt1*xt2**2-6*xt1*
     . xt2+3*xt1-2*xt2**2+4*xt2-2)+log(xt2)**2*xt2*(8*xt1**2*xt2-5*
     . xt1**2-16*xt1*xt2+10*xt1+8*xt2-5)+4*log(xt2)*xt2*(-3*xt1**2*
     . xt2+2*xt1**2+6*xt1*xt2-4*xt1-3*xt2+2)+4*t134p(mt1,mt1,mg)*xgl*
     . (xt1*xt2**2-2*xt1*xt2+xt1+xt2**2-2*xt2+1)+4*t134p(mt2,mt2,mg)*
     . xgl*(-xt1**2*xt2-xt1**2+2*xt1*xt2+2*xt1-xt2-1)+4*t134p(mu,mt1,
     . mg)*xgl*(-xt1*xt2**2+2*xt1*xt2-xt1-xt2**2+2*xt2-1)+4*t134p(mu,
     . mt2,mg)*xgl*(xt1**2*xt2+xt1**2-2*xt1*xt2-2*xt1+xt2+1)+8*(xt1**
     . 2*xt2-xt1**2-xt1*xt2**2+xt1+xt2**2-xt2)
      r72=ans1/(8*(xt1**3*xt2**2-2*xt1**3*xt2+xt1**3-xt1**2*xt2**3+3*
     . xt1**2*xt2-2*xt1**2+2*xt1*xt2**3-3*xt1*xt2**2+xt1-xt2**3+2*xt2
     . **2-xt2))

      r82=(log(xb1)**2*(xb1*xb2**2-2*xb1*xb2+xb1-2*xb2**2+4*xb2-2)+4*
     . log(xb1)*(-xb1**2*xb2**2+2*xb1**2*xb2-xb1**2+xb1*xb2**2-2*xb1*
     . xb2+xb1+xb2**2-2*xb2+1)+log(xb2)**2*(-xb1**2*xb2+2*xb1**2+2*
     . xb1*xb2-4*xb1-xb2+2)+4*log(xb2)*(xb1**2*xb2**2-xb1**2*xb2-xb1
     . **2-2*xb1*xb2**2+2*xb1*xb2+2*xb1+xb2**2-xb2-1)+2*t134p(mb1,mg,
     . mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2-2*xb2+1)+2*t134p(mb2,mg,
     . mg)*(xb1**2*xb2-xb1**2-2*xb1*xb2+2*xb1+xb2-1)+10*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2))/(8*(xb1**3*xb2**2-2*xb1**3*
     . xb2+xb1**3-xb1**2*xb2**3+3*xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*
     . xb1*xb2**2+xb1-xb2**3+2*xb2**2-xb2))

      ans2=4*log(xt2)*xt2*(-3*xt1**3*xt2+2*xt1**3+3*xt1**2*xt2**2+4*
     . xt1**2*xt2-4*xt1**2-6*xt1*xt2**2+xt1*xt2+2*xt1+3*xt2**2-2*xt2)
     . +4*(xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*
     . xt2**3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)
      ans1=4*log(xgl)*log(xt1)*xt1*(2*xt1**2*xt2**2-4*xt1**2*xt2+2*
     . xt1**2-2*xt1*xt2**3+3*xt1*xt2**2-xt1+xt2**3-2*xt2**2+xt2)+4*
     . log(xgl)*log(xt2)*xt2*(-2*xt1**3*xt2+xt1**3+2*xt1**2*xt2**2+3*
     . xt1**2*xt2-2*xt1**2-4*xt1*xt2**2+xt1+2*xt2**2-xt2)+4*log(xgl)*
     . (xt1**3*xt2-xt1**3-2*xt1**2*xt2**2+xt1**2*xt2+xt1**2+xt1*xt2**
     . 3+xt1*xt2**2-2*xt1*xt2-xt2**3+xt2**2)+log(xt1)**2*xt1*(-6*xt1
     . **2*xt2**2+12*xt1**2*xt2-6*xt1**2+2*xt1*xt2**3-xt1*xt2**2-4*
     . xt1*xt2+3*xt1+xt2**3-2*xt2**2+xt2)+4*log(xt1)*log(xt2)*xt1*xt2
     . *(xt1**2*xt2-xt1**2+xt1*xt2**2-4*xt1*xt2+3*xt1-xt2**2+3*xt2-2)
     . +4*log(xt1)*xt1*(3*xt1**2*xt2**2-6*xt1**2*xt2+3*xt1**2-3*xt1*
     . xt2**3+4*xt1*xt2**2+xt1*xt2-2*xt1+2*xt2**3-4*xt2**2+2*xt2)+log
     . (xt2)**2*xt2*(2*xt1**3*xt2+xt1**3-6*xt1**2*xt2**2-xt1**2*xt2-2
     . *xt1**2+12*xt1*xt2**2-4*xt1*xt2+xt1-6*xt2**2+3*xt2)+ans2
      r92=ans1/(8*(xt1**4*xt2**2-2*xt1**4*xt2+xt1**4-2*xt1**3*xt2**3+
     . 2*xt1**3*xt2**2+2*xt1**3*xt2-2*xt1**3+xt1**2*xt2**4+2*xt1**2*
     . xt2**3-6*xt1**2*xt2**2+2*xt1**2*xt2+xt1**2-2*xt1*xt2**4+2*xt1*
     . xt2**3+2*xt1*xt2**2-2*xt1*xt2+xt2**4-2*xt2**3+xt2**2))

      mt1 = amst1
      mt2 = amst2
      xt1 = amst1**2/amg**2 * rim
      xt2 = amst2**2/amg**2 * rim

      ralsca2=(-((9*log(xb1)-10)*(xb2-1)*log(xb1)*xb1-(9*log(xb2)-10)
     . *(xb1-1)*log(xb2)*xb2))/(48*(xb1-xb2)*(xb1-1)*(xb2-1))

      ralscf2=(-((xb1-1)*log(xb2)*xb2-(xb2-1)*log(xb1)*xb1))/(8*(xb1-
     . xb2)*(xb1-1)*(xb2-1))

      rmb12=(-((((2*xb1**3-3*xb2-(xb2-6)*xb1**2-2*(xb2+1)*xb1)*xb1+(
     . xb1**2-xb2)*(xb1-1)**2*log(-(xb1-1)))*(xb2-1)+(xb1-1)**2*log(
     . xb2)*xb1**2*xb2)*log(xb1)-(((xb1**2-xb2)*(xb1+1)*(xb2-1)*log(
     . xb1)**2-(xb1-1)**2*log(xb2)**2*xb2)*xb1+((xb1-xb2)*(xb2-1)+(
     . xb1-1)*log(xb2)*xb2)*((xb1+5)*xb1+(xb1-1)**2*log(-(xb1-1)))*(
     . xb1-1))))/(4*(xb1-xb2)**2*(xb1-1)**2*(xb2-1)*xb1)

      rmb22=(-(((xb1-xb2**2)*(xb1-1)*(xb2+1)*log(xb2)**2+(xb2-1)**2*
     . log(xb1)**2*xb1)*xb2+((xb2+5)*xb2+(xb2-1)**2*log(-(xb2-1)))*(
     . xb1-xb2)*(xb1-1)*(xb2-1)-((xb2-1)**2*log(-(xb2-1))-log(xb2)*
     . xb2**2+(xb2+5)*xb2)*(xb2-1)**2*log(xb1)*xb1+((2*(xb2**2+3*xb2-
     . 1)*xb2-(xb2**2+2*xb2+3)*xb1)*xb2-(xb1-xb2**2)*(xb2-1)**2*log(-
     . (xb2-1)))*(xb1-1)*log(xb2)))/(4*(xb1-xb2)**2*(xb1-1)*(xb2-1)**
     . 2*xb2)

      rmgca2=(-((3*log(xb1)-14)*(xb1+1)*(xb2-1)**2*log(xb1)*xb1-(3*
     . log(xb2)-14)*(xb1-1)**2*(xb2+1)*log(xb2)*xb2-28*(xb1-xb2)*(xb1
     . -1)*(xb2-1)))/(16*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsu1
      ms2 = amsu2
      xs1 = amsu1**2/amg**2 * rim
      xs2 = amsu2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      ru102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rualstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rumgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsd1
      ms2 = amsd2
      xs1 = amsd1**2/amg**2 * rim
      xs2 = amsd2**2/amg**2 * rim

      ans2=2*t134p(mb2,ms1,mg)*(xb1**2*xb2-xb1**2*xs1-2*xb1*xb2+2*xb1
     . *xs1+xb2-xs1)+2*t134p(mb2,ms2,mg)*(xb1**2*xb2-xb1**2*xs2-2*xb1
     . *xb2+2*xb1*xs2+xb2-xs2)+2*t134p(mg,ms1,mg)*(-2*xb1**2*xb2+xb1
     . **2*xs1+xb1**2+2*xb1*xb2**2-2*xb1*xs1-xb2**2*xs1-xb2**2+2*xb2*
     . xs1)+2*t134p(mg,ms2,mg)*(-2*xb1**2*xb2+xb1**2*xs2+xb1**2+2*xb1
     . *xb2**2-2*xb1*xs2-xb2**2*xs2-xb2**2+2*xb2*xs2)+12*(xb1**2*xb2*
     . xs1+xb1**2*xb2*xs2+xb1**2*xb2-xb1**2*xs1-xb1**2*xs2-xb1**2-xb1
     . *xb2**2*xs1-xb1*xb2**2*xs2-xb1*xb2**2+xb1*xs1+xb1*xs2+xb1+xb2
     . **2*xs1+xb2**2*xs2+xb2**2-xb2*xs1-xb2*xs2-xb2)
      ans1=log(xb1)**2*xb1*(-xb2**2*xs1-xb2**2*xs2+2*xb2*xs1+2*xb2*
     . xs2-xs1-xs2)+2*log(xb1)*log(xs1)*xb1*xs1*(-xb2**2+2*xb2-1)+2*
     . log(xb1)*log(xs2)*xb1*xs2*(-xb2**2+2*xb2-1)+4*log(xb1)*xb1*(
     . xb2**2*xs1+xb2**2*xs2-2*xb2*xs1-2*xb2*xs2+xs1+xs2)+log(xb2)**2
     . *xb2*(xb1**2*xs1+xb1**2*xs2-2*xb1*xs1-2*xb1*xs2+xs1+xs2)+2*log
     . (xb2)*log(xs1)*xb2*xs1*(xb1**2-2*xb1+1)+2*log(xb2)*log(xs2)*
     . xb2*xs2*(xb1**2-2*xb1+1)+4*log(xb2)*xb2*(-xb1**2*xs1-xb1**2*
     . xs2+2*xb1*xs1+2*xb1*xs2-xs1-xs2)+log(xs1)**2*xs1*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs1)*xs1*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+log(xs2)**2*xs2*(xb1**2*xb2-
     . xb1**2-xb1*xb2**2+xb1+xb2**2-xb2)+8*log(xs2)*xs2*(-xb1**2*xb2+
     . xb1**2+xb1*xb2**2-xb1-xb2**2+xb2)+2*t134p(mb1,ms1,mg)*(-xb1*
     . xb2**2+2*xb1*xb2-xb1+xb2**2*xs1-2*xb2*xs1+xs1)+2*t134p(mb1,ms2
     . ,mg)*(-xb1*xb2**2+2*xb1*xb2-xb1+xb2**2*xs2-2*xb2*xs2+xs2)+ans2
      rd102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rdalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rdmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amsb1
      ms2 = amsb2
      xs1 = amsb1**2/amg**2 * rim
      xs2 = amsb2**2/amg**2 * rim

      ans1=log(xb1)**2*xb1*(xb1**2*xb2-xb1**2-4*xb1*xb2**2+6*xb1*xb2-
     . 2*xb1-xb2**3+3*xb2**2-2*xb2)+2*log(xb1)*log(xb2)*xb1*xb2*(xb1
     . **2-2*xb1-xb2**2+2*xb2)+4*log(xb1)*xb1*(-2*xb1**2*xb2+2*xb1**2
     . +3*xb1*xb2**2-2*xb1*xb2-xb1+xb2**3-4*xb2**2+3*xb2)+log(xb2)**2
     . *xb2*(xb1**3+4*xb1**2*xb2-3*xb1**2-xb1*xb2**2-6*xb1*xb2+2*xb1+
     . xb2**2+2*xb2)+4*log(xb2)*xb2*(-xb1**3-3*xb1**2*xb2+4*xb1**2+2*
     . xb1*xb2**2+2*xb1*xb2-3*xb1-2*xb2**2+xb2)+2*t134p(mb1,mb2,mg)*(
     . -xb1**3+xb1**2*xb2+2*xb1**2-xb1*xb2**2-2*xb1+xb2**3-2*xb2**2+2
     . *xb2)+2*t134p(mb1,mg,mg)*(xb1**3-2*xb1**2*xb2-xb1**2+xb1*xb2**
     . 2+2*xb1*xb2-xb2**2)+2*t134p(mb2,mg,mg)*(-xb1**2*xb2+xb1**2+2*
     . xb1*xb2**2-2*xb1*xb2-xb2**3+xb2**2)+12*(xb1**3*xb2-xb1**3-xb1*
     . xb2**3+xb1+xb2**3-xb2)
      rb102=ans1/(8*(xb1**3*xb2**2-2*xb1**3*xb2+xb1**3-xb1**2*xb2**3+3
     . *xb1**2*xb2-2*xb1**2+2*xb1*xb2**3-3*xb1*xb2**2+xb1-xb2**3+2*
     . xb2**2-xb2))

      rbalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rbmgtr2=((xs2-6+xs1+log(xs2)+log(xs1)+(log(xs2-1)-log(xs2))*(xs2
     . -1)**2+(log(xs1-1)-log(xs1))*(xs1-1)**2)*((xb1+1)*(xb2-1)**2*
     . log(xb1)*xb1-(xb1-1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs2-1)-
     . log(xs2))*(xs2-1)**2+log(xs1)+log(xs2)+(log(xs1-1)-log(xs1))*(
     . xs1-1)**2)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*(xs2-6+xs1)*(
     . xb1-xb2)*(xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      ms1 = amst1
      ms2 = amst2
      xs1 = amst1**2/amg**2 * rim
      xs2 = amst2**2/amg**2 * rim

      ans14=-8*log(xs2)*xb1*xb2**2*xs2-4*log(xs2)*xb1*xq**2*xs2-4*log
     . (xs2)*xb1*xq*xs2**2-4*log(xs2)*xb1*xq*xs2+8*log(xs2)*xb1*xs2**
     . 3-16*log(xs2)*xb1*xs2**2+8*log(xs2)*xb1*xs2-4*log(xs2)*xb2**2*
     . xq**2*xs2-4*log(xs2)*xb2**2*xq*xs2**2-4*log(xs2)*xb2**2*xq*xs2
     . +8*log(xs2)*xb2**2*xs2**3-16*log(xs2)*xb2**2*xs2**2+8*log(xs2)
     . *xb2**2*xs2+4*log(xs2)*xb2*xq**2*xs2+4*log(xs2)*xb2*xq*xs2**2+
     . 4*log(xs2)*xb2*xq*xs2-8*log(xs2)*xb2*xs2**3+16*log(xs2)*xb2*
     . xs2**2-8*log(xs2)*xb2*xs2
      ans13=-log(xs2)**2*xb1*xb2**2*xq**2*xs2+log(xs2)**2*xb1*xb2**2*
     . xs2**3-2*log(xs2)**2*xb1*xb2**2*xs2**2+log(xs2)**2*xb1*xb2**2*
     . xs2+log(xs2)**2*xb1*xq**2*xs2-log(xs2)**2*xb1*xs2**3+2*log(xs2
     . )**2*xb1*xs2**2-log(xs2)**2*xb1*xs2+log(xs2)**2*xb2**2*xq**2*
     . xs2-log(xs2)**2*xb2**2*xs2**3+2*log(xs2)**2*xb2**2*xs2**2-log(
     . xs2)**2*xb2**2*xs2-log(xs2)**2*xb2*xq**2*xs2+log(xs2)**2*xb2*
     . xs2**3-2*log(xs2)**2*xb2*xs2**2+log(xs2)**2*xb2*xs2-4*log(xs2)
     . *xb1**2*xb2*xq**2*xs2-4*log(xs2)*xb1**2*xb2*xq*xs2**2-4*log(
     . xs2)*xb1**2*xb2*xq*xs2+8*log(xs2)*xb1**2*xb2*xs2**3-16*log(xs2
     . )*xb1**2*xb2*xs2**2+8*log(xs2)*xb1**2*xb2*xs2+4*log(xs2)*xb1**
     . 2*xq**2*xs2+4*log(xs2)*xb1**2*xq*xs2**2+4*log(xs2)*xb1**2*xq*
     . xs2-8*log(xs2)*xb1**2*xs2**3+16*log(xs2)*xb1**2*xs2**2-8*log(
     . xs2)*xb1**2*xs2+4*log(xs2)*xb1*xb2**2*xq**2*xs2+4*log(xs2)*xb1
     . *xb2**2*xq*xs2**2+4*log(xs2)*xb1*xb2**2*xq*xs2-8*log(xs2)*xb1*
     . xb2**2*xs2**3+16*log(xs2)*xb1*xb2**2*xs2**2+ans14
      ans12=-4*log(xb2)**2*xb1*xb2*xq*xs1-8*log(xb2)**2*xb1*xb2*xq*
     . xs2**2+4*log(xb2)**2*xb1*xb2*xq*xs2-4*log(xb2)**2*xb1*xb2*xq+2
     . *log(xb2)**2*xb1*xb2*xs1*xs2**2-4*log(xb2)**2*xb1*xb2*xs1*xs2+
     . 2*log(xb2)**2*xb1*xb2*xs1+2*log(xb2)**2*xb1*xb2*xs2**3-4*log(
     . xb2)**2*xb1*xb2*xs2**2+2*log(xb2)**2*xb1*xb2*xs2+2*log(xb2)**2
     . *xb2*xq**3-log(xb2)**2*xb2*xq**2*xs1-5*log(xb2)**2*xb2*xq**2*
     . xs2-4*log(xb2)**2*xb2*xq**2+2*log(xb2)**2*xb2*xq*xs1*xs2+2*log
     . (xb2)**2*xb2*xq*xs1+4*log(xb2)**2*xb2*xq*xs2**2-2*log(xb2)**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb2*xq-log(xb2)**2*xb2*xs1*xs2**2+2*
     . log(xb2)**2*xb2*xs1*xs2-log(xb2)**2*xb2*xs1-log(xb2)**2*xb2*
     . xs2**3+2*log(xb2)**2*xb2*xs2**2-log(xb2)**2*xb2*xs2+log(xs2)**
     . 2*xb1**2*xb2*xq**2*xs2-log(xs2)**2*xb1**2*xb2*xs2**3+2*log(xs2
     . )**2*xb1**2*xb2*xs2**2-log(xs2)**2*xb1**2*xb2*xs2-log(xs2)**2*
     . xb1**2*xq**2*xs2+log(xs2)**2*xb1**2*xs2**3-2*log(xs2)**2*xb1**
     . 2*xs2**2+log(xs2)**2*xb1**2*xs2+ans13
      ans11=4*log(xb1)**2*xb1*xq**2-2*log(xb1)**2*xb1*xq*xs1*xs2-2*
     . log(xb1)**2*xb1*xq*xs1-4*log(xb1)**2*xb1*xq*xs2**2+2*log(xb1)
     . **2*xb1*xq*xs2-2*log(xb1)**2*xb1*xq+log(xb1)**2*xb1*xs1*xs2**2
     . -2*log(xb1)**2*xb1*xs1*xs2+log(xb1)**2*xb1*xs1+log(xb1)**2*xb1
     . *xs2**3-2*log(xb1)**2*xb1*xs2**2+log(xb1)**2*xb1*xs2+2*log(xb2
     . )**2*xb1**2*xb2*xq**3-log(xb2)**2*xb1**2*xb2*xq**2*xs1-5*log(
     . xb2)**2*xb1**2*xb2*xq**2*xs2-4*log(xb2)**2*xb1**2*xb2*xq**2+2*
     . log(xb2)**2*xb1**2*xb2*xq*xs1*xs2+2*log(xb2)**2*xb1**2*xb2*xq*
     . xs1+4*log(xb2)**2*xb1**2*xb2*xq*xs2**2-2*log(xb2)**2*xb1**2*
     . xb2*xq*xs2+2*log(xb2)**2*xb1**2*xb2*xq-log(xb2)**2*xb1**2*xb2*
     . xs1*xs2**2+2*log(xb2)**2*xb1**2*xb2*xs1*xs2-log(xb2)**2*xb1**2
     . *xb2*xs1-log(xb2)**2*xb1**2*xb2*xs2**3+2*log(xb2)**2*xb1**2*
     . xb2*xs2**2-log(xb2)**2*xb1**2*xb2*xs2-4*log(xb2)**2*xb1*xb2*xq
     . **3+2*log(xb2)**2*xb1*xb2*xq**2*xs1+10*log(xb2)**2*xb1*xb2*xq
     . **2*xs2+8*log(xb2)**2*xb1*xb2*xq**2-4*log(xb2)**2*xb1*xb2*xq*
     . xs1*xs2+ans12
      ans10=5*log(xb1)**2*xb1*xb2**2*xq**2*xs2+4*log(xb1)**2*xb1*xb2
     . **2*xq**2-2*log(xb1)**2*xb1*xb2**2*xq*xs1*xs2-2*log(xb1)**2*
     . xb1*xb2**2*xq*xs1-4*log(xb1)**2*xb1*xb2**2*xq*xs2**2+2*log(xb1
     . )**2*xb1*xb2**2*xq*xs2-2*log(xb1)**2*xb1*xb2**2*xq+log(xb1)**2
     . *xb1*xb2**2*xs1*xs2**2-2*log(xb1)**2*xb1*xb2**2*xs1*xs2+log(
     . xb1)**2*xb1*xb2**2*xs1+log(xb1)**2*xb1*xb2**2*xs2**3-2*log(xb1
     . )**2*xb1*xb2**2*xs2**2+log(xb1)**2*xb1*xb2**2*xs2+4*log(xb1)**
     . 2*xb1*xb2*xq**3-2*log(xb1)**2*xb1*xb2*xq**2*xs1-10*log(xb1)**2
     . *xb1*xb2*xq**2*xs2-8*log(xb1)**2*xb1*xb2*xq**2+4*log(xb1)**2*
     . xb1*xb2*xq*xs1*xs2+4*log(xb1)**2*xb1*xb2*xq*xs1+8*log(xb1)**2*
     . xb1*xb2*xq*xs2**2-4*log(xb1)**2*xb1*xb2*xq*xs2+4*log(xb1)**2*
     . xb1*xb2*xq-2*log(xb1)**2*xb1*xb2*xs1*xs2**2+4*log(xb1)**2*xb1*
     . xb2*xs1*xs2-2*log(xb1)**2*xb1*xb2*xs1-2*log(xb1)**2*xb1*xb2*
     . xs2**3+4*log(xb1)**2*xb1*xb2*xs2**2-2*log(xb1)**2*xb1*xb2*xs2-
     . 2*log(xb1)**2*xb1*xq**3+log(xb1)**2*xb1*xq**2*xs1+5*log(xb1)**
     . 2*xb1*xq**2*xs2+ans11
      ans9=-4*t134(mb1,mq,ms2,mg)*xb1*xb2+2*t134(mb1,mq,ms2,mg)*xb1*
     . xq**2-4*t134(mb1,mq,ms2,mg)*xb1*xq*xs2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xq+2*t134(mb1,mq,ms2,mg)*xb1*xs2**2-4*t134(mb1,mq,ms2,mg)*
     . xb1*xs2+2*t134(mb1,mq,ms2,mg)*xb1+2*t134(mb1,mq,ms2,mg)*xb2**2
     . *xq**3-6*t134(mb1,mq,ms2,mg)*xb2**2*xq**2*xs2-4*t134(mb1,mq,
     . ms2,mg)*xb2**2*xq**2+6*t134(mb1,mq,ms2,mg)*xb2**2*xq*xs2**2+2*
     . t134(mb1,mq,ms2,mg)*xb2**2*xq-2*t134(mb1,mq,ms2,mg)*xb2**2*xs2
     . **3+4*t134(mb1,mq,ms2,mg)*xb2**2*xs2**2-2*t134(mb1,mq,ms2,mg)*
     . xb2**2*xs2-4*t134(mb1,mq,ms2,mg)*xb2*xq**3+12*t134(mb1,mq,ms2,
     . mg)*xb2*xq**2*xs2+8*t134(mb1,mq,ms2,mg)*xb2*xq**2-12*t134(mb1,
     . mq,ms2,mg)*xb2*xq*xs2**2-4*t134(mb1,mq,ms2,mg)*xb2*xq+4*t134(
     . mb1,mq,ms2,mg)*xb2*xs2**3-8*t134(mb1,mq,ms2,mg)*xb2*xs2**2+4*
     . t134(mb1,mq,ms2,mg)*xb2*xs2+2*t134(mb1,mq,ms2,mg)*xq**3-6*t134
     . (mb1,mq,ms2,mg)*xq**2*xs2-4*t134(mb1,mq,ms2,mg)*xq**2+6*t134(
     . mb1,mq,ms2,mg)*xq*xs2**2+2*t134(mb1,mq,ms2,mg)*xq-2*t134(mb1,
     . mq,ms2,mg)*xs2**3+4*t134(mb1,mq,ms2,mg)*xs2**2-2*t134(mb1,mq,
     . ms2,mg)*xs2-2*log(xb1)**2*xb1*xb2**2*xq**3+log(xb1)**2*xb1*xb2
     . **2*xq**2*xs1+ans10
      ans8=4*t134(mb1,mq,ms1,mg)*xb2*xq**2*xs1+8*t134(mb1,mq,ms1,mg)*
     . xb2*xq**2*xs2+8*t134(mb1,mq,ms1,mg)*xb2*xq**2-8*t134(mb1,mq,
     . ms1,mg)*xb2*xq*xs1*xs2-8*t134(mb1,mq,ms1,mg)*xb2*xq*xs1-4*t134
     . (mb1,mq,ms1,mg)*xb2*xq*xs2**2+8*t134(mb1,mq,ms1,mg)*xb2*xq*xs2
     . -4*t134(mb1,mq,ms1,mg)*xb2*xq+4*t134(mb1,mq,ms1,mg)*xb2*xs1*
     . xs2**2-8*t134(mb1,mq,ms1,mg)*xb2*xs1*xs2+4*t134(mb1,mq,ms1,mg)
     . *xb2*xs1+2*t134(mb1,mq,ms1,mg)*xq**3-2*t134(mb1,mq,ms1,mg)*xq
     . **2*xs1-4*t134(mb1,mq,ms1,mg)*xq**2*xs2-4*t134(mb1,mq,ms1,mg)*
     . xq**2+4*t134(mb1,mq,ms1,mg)*xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*
     . xq*xs1+2*t134(mb1,mq,ms1,mg)*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xq*xs2+2*t134(mb1,mq,ms1,mg)*xq-2*t134(mb1,mq,ms1,mg)*xs1*xs2
     . **2+4*t134(mb1,mq,ms1,mg)*xs1*xs2-2*t134(mb1,mq,ms1,mg)*xs1+2*
     . t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2**2*xq*xs2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*xq+2*t134(mb1,
     . mq,ms2,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms2,mg)*xb1*xb2**2*
     . xs2+2*t134(mb1,mq,ms2,mg)*xb1*xb2**2-4*t134(mb1,mq,ms2,mg)*xb1
     . *xb2*xq**2+8*t134(mb1,mq,ms2,mg)*xb1*xb2*xq*xs2+8*t134(mb1,mq,
     . ms2,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms2,mg)*xb1*xb2*xs2**2+8*t134
     . (mb1,mq,ms2,mg)*xb1*xb2*xs2+ans9
      ans7=2*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq**2-4*t134(mb1,mq,ms1,
     . mg)*xb1*xb2**2*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xb2**2*xq+2*
     . t134(mb1,mq,ms1,mg)*xb1*xb2**2*xs2**2-4*t134(mb1,mq,ms1,mg)*
     . xb1*xb2**2*xs2+2*t134(mb1,mq,ms1,mg)*xb1*xb2**2-4*t134(mb1,mq,
     . ms1,mg)*xb1*xb2*xq**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xq*xs2+8*
     . t134(mb1,mq,ms1,mg)*xb1*xb2*xq-4*t134(mb1,mq,ms1,mg)*xb1*xb2*
     . xs2**2+8*t134(mb1,mq,ms1,mg)*xb1*xb2*xs2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xb2+2*t134(mb1,mq,ms1,mg)*xb1*xq**2-4*t134(mb1,mq,ms1,mg)
     . *xb1*xq*xs2-4*t134(mb1,mq,ms1,mg)*xb1*xq+2*t134(mb1,mq,ms1,mg)
     . *xb1*xs2**2-4*t134(mb1,mq,ms1,mg)*xb1*xs2+2*t134(mb1,mq,ms1,mg
     . )*xb1+2*t134(mb1,mq,ms1,mg)*xb2**2*xq**3-2*t134(mb1,mq,ms1,mg)
     . *xb2**2*xq**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2**2*xq**2*xs2-4*
     . t134(mb1,mq,ms1,mg)*xb2**2*xq**2+4*t134(mb1,mq,ms1,mg)*xb2**2*
     . xq*xs1*xs2+4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs1+2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xq*xs2**2-4*t134(mb1,mq,ms1,mg)*xb2**2*xq*xs2+2
     . *t134(mb1,mq,ms1,mg)*xb2**2*xq-2*t134(mb1,mq,ms1,mg)*xb2**2*
     . xs1*xs2**2+4*t134(mb1,mq,ms1,mg)*xb2**2*xs1*xs2-2*t134(mb1,mq,
     . ms1,mg)*xb2**2*xs1-4*t134(mb1,mq,ms1,mg)*xb2*xq**3+ans8
      ans15=(xq**2-2*xq*xs1-2*xq+xs1**2-2*xs1+1)
      ans6=ans7*ans15
      ans5=-ans6
      ans4=(4*(2*(xs1-1)**2-xq**2-(xs1+1)*xq)-(xs1-1+xq)*(xs1-1-xq)*
     . log(xs1))*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*(xb1-1)*(
     . xb2-1)*log(xs1)*xs1+ans5
      ans3=-ans4
      ans17=4*(3*((xs1+1)*xs1+xs2**2+xs2-1)+2*(3*xq**2+2)*xq+((2*xq**
     . 2-5*xq+8)*xq-3*(xs1+xs2)**2)*(xs1+xs2)-(4*xq-3*xs1*xs2)*(2*xs1
     . **2+3*xs1*xs2+2*xs2**2)+2*(3*xs1**2-8*xs1*xs2+3*xs2**2)*(xs1+
     . xs2)*xq-(13*xq**2+4*xs1*xs2)*xq**2+((3*xs1**2-2*xs1*xs2+3*xs2
     . **2)*xq-2*(xs1+xs2)*xq**2-3*(xs1+xs2)*xs1*xs2+6*xq**3)*(xq-xs1
     . )*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)+4*(2*((xs2-1)**2+xs1**2+
     . 2*(2*xs2-1)*xs1+(5*xq**2-1)*xq)-(xs1**2-4*xs1*xs2+xs2**2)*xq-2
     . *(xq**2+2*xs1*xs2)*(xs1+xs2)+(xs1-1-xq)*(2*(xs2+1)*xq-(xs2-1)
     . **2-xq**2)*log(xs1)*xs1+(2*(xs1+1)*xq-(xs1-1)**2-xq**2)*(xs2-1
     . -xq)*log(xs2)*xs2+(xs1+xs2-6*xq)*xq+((xs1+xs2)*xq-4*xq**2+2*
     . xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq)*
     . xq
      ans16=-2*(2*((xs1-1)**2*xs1-xq**3-(3*xs1-1)*xq*xs1+(3*xs1+1)*xq
     . **2)+(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1)*xb2-(2*((xs1+1
     . )*xq-(xs1-1)**2)*xb2-(xq**2-2*xq*xs1+xs1**2-2*xs1+1)*(xq-xs1-1
     . ))*xb1)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)*(xb1-xb2)*t134(mg,mq,
     . ms1,mg)-2*((xs2-1)**2+xs1**2+2*(2*xs2-1)*xs1+2*(xq+1)*(xq-1)*
     . xq+(xq-2*xs1*xs2)*(xs1+xs2)-(xs1**2-4*xs1*xs2+xs2**2)*xq-(xq**
     . 2-xs1*xs2)*(xq-xs1)*(xq-xs2))*(xb1-xb2)*(xb1-1)*(xb2-1)*log(xq
     . )**2*xq-2*((2*((xs2-1)**2*xs2-xq**3-(3*xs2-1)*xq*xs2+(3*xs2+1)
     . *xq**2)+(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-xs2-1)*xb2-(2*((
     . xs2+1)*xq-(xs2-1)**2)*xb2-(xq**2-2*xq*xs2+xs2**2-2*xs2+1)*(xq-
     . xs2-1))*xb1)*(xb1-xb2)*t134(mg,mq,ms2,mg)+((xq-xs1+xb2)*t134(
     . mb2,mq,ms1,mg)+(xq-xs2+xb2)*t134(mb2,mq,ms2,mg))*(2*(xs2+1)*xq
     . -(xs2-1)**2-xq**2)*(xb1-1)**2)*(2*(xs1+1)*xq-(xs1-1)**2-xq**2)
     . +ans17
      ans2=2*(2*(xs1+xs2-2*xq)-log(xs2)*xs2-log(xs1)*xs1+2*log(xq)*xq
     . )*((xb1-1)**2*log(xb2)*xb2-(xb2-1)**2*log(xb1)*xb1)*(2*(xs1+1)
     . *xq-(xs1-1)**2-xq**2)*(2*(xs2+1)*xq-(xs2-1)**2-xq**2)+ans3+
     . ans16
      ans1=-ans2
      rt102=ans1/(8*((xs1-1)**2+xq**2-2*(xs1+1)*xq)*((xs2-1)**2+xq**2-
     . 2*(xs2+1)*xq)*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      rtalstr2=(-((log(xs2)-6+log(xs1))*((xb1-1)*log(xb2)*xb2-(xb2-1)*
     . log(xb1)*xb1)+3*((xb1-1)*log(xb2)**2*xb2-(xb2-1)*log(xb1)**2*
     . xb1)))/(24*(xb1-xb2)*(xb1-1)*(xb2-1))

      rtmgtr2=((xs2-6+xs1-2*xq+(xq+1)*log(xs2)+(xq+1)*log(xs1)-2*log(
     . xq)*xq-(xs2-1-xq)*log(xp(ms2))*xp(ms2)-(xs1-1-xq)*log(xp(ms1))
     . *xp(ms1)-(xs2-1-xq)*log(xm(ms2))*xm(ms2)-(xs1-1-xq)*log(xm(ms1
     . ))*xm(ms1)+(xs2-1-xq)*log(xp(ms2)-1)*xp(ms2)+(xs1-1-xq)*log(xp
     . (ms1)-1)*xp(ms1)+(xs2-1-xq)*log(xm(ms2)-1)*xm(ms2)+(xs1-1-xq)*
     . log(xm(ms1)-1)*xm(ms1))*((xb1+1)*(xb2-1)**2*log(xb1)*xb1-(xb1-
     . 1)**2*(xb2+1)*log(xb2)*xb2)+2*((log(xs1)+log(xs2))*(xq+1)-2*
     . log(xq)*xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+(xb1+1)*(xb2-1)**2*log(
     . xb1)**2*xb1-(xb1-1)**2*(xb2+1)*log(xb2)**2*xb2+2*((log(xm(ms2)
     . -1)-log(xm(ms2)))*xm(ms2)+(log(xp(ms2)-1)-log(xp(ms2)))*xp(ms2
     . ))*(xs2-1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*((log(xm(ms1)-1)-log
     . (xm(ms1)))*xm(ms1)+(log(xp(ms1)-1)-log(xp(ms1)))*xp(ms1))*(xs1
     . -1-xq)*(xb1-xb2)*(xb1-1)*(xb2-1)+2*(xs2-6+xs1-2*xq)*(xb1-xb2)*
     . (xb1-1)*(xb2-1))/(8*(xb1-xb2)*(xb1-1)**2*(xb2-1)**2)

      rctca = ralsca2 + rmgca2
      rctcf = ralscf2 + rmb12 + rmb22
      ructtr = rualstr2 + rumgtr2
      rdcttr = rdalstr2 + rdmgtr2
      rbcttr = rbalstr2 + rbmgtr2
      rtcttr = rtalstr2 + rtmgtr2

      bo  =-2*fi(amsb1**2,amsb2**2,amg**2)

      rca = r12 + r22/4 + r32 + r52 + r62 + rctca
      rcf = -r22/2 - r42/2 - 2*r52 - r72/2 + r82 - r92/2 + rctcf
      rtr = nu*(ru102 + ructtr) + nd*(rd102 + rdcttr)
      rtrb= rb102 + rbcttr
      rtrt= rt102 + rtcttr + bo/3*dlog(amg**2/amt**2)/fnorm

      rqcd = ca*rca + cf*rcf + tr*(rtr+rtrb+rtrt)
      rqcd = rqcd/bo
      anomalous = - cf/4
      finscale = (11*ca-4*tr*nf)/12*dlog(scale**2/amg**2)
      fqcd_hdec = dreal(rqcd)*fnorm + anomalous + finscale

      return
      end
 
c end change susyhit
