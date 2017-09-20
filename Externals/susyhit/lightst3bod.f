c the following routine was written by Ramona Gröber on 9/10/14
c it computes the light stop 3-body decay into Wbchi_0 allowing for FV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                     light stop 3-body decays			      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine SD_lightstop3bod(width_tot3bod,width_bchiW, 
     .                            width_jchiw)
      implicit none 
      integer i, j, k, s, ifavvio
      integer nxt1, nyt1
      double precision amsupq(6), amsdownq(6), 
     .amslepton(6), amsneutrino(3)
      double precision sdgf,amz,amw,pi,g2
      double precision width_bchiW, width_jchiw, width_tot3bod
      double precision amt, amb, amtau, mfd(3), mfu(3), mfe(3)
      double precision DcharL(3,6,2), EcharR(3,6,2), 
     . cchaneutL(4,2), cchaneutR(4,2), ccharl(3,6,2), fcharr(3,6,2)
      double precision gneutul(3,4,6), 
     . gneutur(3,4,6), gneutdl(3,4,6), gneutdr(3,4,6)
      double precision gsqsqW(6)
      double precision amneut(4),xmneut(4),amchar(2),xmchar(2)
      integer nx1t, ny1t
      double precision xmu1, xmu2, xmu3
      double precision sum1, sum2
      double precision SD_ay,SD_by,SD_ax,SD_bx
      external SD_ay,SD_by,SD_ax,SD_bx

      double precision SD_3bodbchiW, SD_3bodjchiW
c----- intialize parameters for couplings
      double precision uu(2,2), vv(2,2), zz(4,4), zp(4,4)
      double precision vckm(3,3), msq2(3,3), msu2(3,3), 
     .msd2(3,3), td(3,3), tu(3,3),
     . usqmix(6,6), dsqmix(6,6)
      double precision tanbeta, alp_mssm
      double precision alsew,g2ew,g1ew
      double precision amwp, amzp
c---- ramona changed 2/2/15
      double precision topwidth
c---- end ramona changed

      external SD_3bodbchiw, SD_3bodjchiw

      COMMON/msfermion/ amsupq, amsdownq, amslepton, amsneutrino
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      COMMON/SD_fermion/amt,amb,amtau
      COMMON/SD_massgino/amneut,xmneut,amchar,xmchar
      COMMON/SD_nx1/nx1t,ny1t
      common/charcoup/DcharL, EcharR, FcharR, CCharL
      common/neutcoup/Gneutul, Gneutur, Gneutdr, Gneutdl
      common/Wcharneutcoup/CchaneutL, CchaneutR
      common/squarkWcoup/gsqsqW
      COMMON/SD_mwmzpole/amwp,amzp

c---- common block for parameters entering the couplings
      COMMON/flavviolation/vckm, msq2, msd2, msu2, td, 
     .tu, usqmix, ifavvio, dsqmix
      COMMON/SD_mixmat/uu,vv,zz,zp
      COMMON/SD_mixang/alp_mssm,tanbeta
      COMMON/SD_coupewsb/alsew,g2ew,g1ew
c---- ramona cahnged 2/2/15
      COMMON/widthtop/topwidth
c----- end ramona changed

      
c---- set fermion masses
      mfd(1)=0d0
      mfd(2)=0d0
      mfd(3)=amb
      mfu(1)=0d0
      mfu(2)=0d0
      mfu(3)=amt
      mfe(1)=0d0
      mfe(2)=0d0
      mfe(3)=amtau

      topwidth=0d0
      do I=1,3
        topwidth=topwidth+(VCKM(3,i)**2)
     . *g2ew**2/(64d0*pi)*amt**3/amw**2*
     . dsqrt(1d0-2d0*(mfd(i)**2+amw**2)/amt**2
     . +(mfd(i)**2-amw**2)**2/amt**4)*(amw**2*(amt**2+mfd(i)**2)
     . +(amt**2-mfd(i)**2)**2-2d0*amw**4)/amt**4
      enddo



c---- couplings, notation see SPARTICLES book
c---- chargino-squark quark coupling
      do i=1,3
      do k=1,2
      do s=1,6
! Attention in sparticles book on which one is with left handed projector and which one is with right handed
! Attention! Definition changed compared to 4-body decay. Indices changed!
      DcharL(i,k,s)=g2ew*vv(k,2)/Sqrt(2d0)/amw/dsin(datan(tanbeta))
     .		    *mfu(3)*Vckm(3,i)*USQMix(s,6)-g2ew*vv(k,1)*
     .		    (Usqmix(s,1)*Vckm(1,i)+Usqmix(s,2)*Vckm(2,i)
     .		    +Usqmix(s,3)*Vckm(3,i))
      EcharR(i,k,s)=g2ew*uu(k,2)/dsqrt(2d0)/amw
     .		    /dcos(datan(tanbeta))*mfd(i)*usqmix(s,3)
     .		    *Vckm(3,i)


      end do
      end do
      end do

 
      
      
      do i=1,3
      do k=1,4
      do s=1,6
      Gneutul(i,k,s)=-Sqrt(2d0)*(1d0/2d0*ZZ(k,2)*g2ew
     .               +1d0/6d0*g1ew*ZZ(k,1))*Usqmix(s,i)-g2ew*mfu(i)
     .	             /(dsqrt(2d0)*amw*Sin(datan(tanbeta)))
     .	             *ZZ(k,4)*Usqmix(s,i+3)
      Gneutur(i,k,s)=2d0*dsqrt(2d0)/3d0*g1ew*ZZ(k,1)
     .		     *Usqmix(s,i+3)-g2ew*mfu(i)
     . 	      	     /(dsqrt(2d0)*amw*Sin(datan(tanbeta)))
     .		     *ZZ(k,4)*Usqmix(s,i)
      
      Gneutdl(i,k,s)=Sqrt(2d0)*(1d0/2d0*ZZ(k,2)*g2ew
     .      	      -1d0/6d0*g1ew*ZZ(k,1))*dsqmix(s,i)-g2ew*mfd(i)
     .	              /(dsqrt(2d0)*amw*Cos(datan(tanbeta)))
     .		      *ZZ(k,3)*dsqmix(s,i+3)
      Gneutdr(i,k,s)=-dsqrt(2d0)/3d0*g1ew*ZZ(k,1)
     .		     *dsqmix(s,i+3)-g2ew*mfd(i)
     . 	      	      /(dsqrt(2d0)*amw*Cos(datan(tanbeta)))
     .		      *ZZ(k,3)*dsqmix(s,i)
      end do
      end do
      end do






c--- neutralino-chargino-W-coupling
      do k=1,2
      do i=1,4
      CchaneutL(i,k)=g2ew*(ZZ(i,2)*VV(k,1)
     .		     -1d0/dsqrt(2d0)*ZZ(i,4)*VV(k,2))
      CChaneutR(i,k)=g2ew*(ZZ(i,2)*UU(k,1)
     .		     +1d0/dsqrt(2d0)*ZZ(i,3)*UU(k,2))
      end do
      end do


c---- W-squark-squark coupling
      do k=1,6
      gsqsqW(k)=-g2ew/dsqrt(2d0)*((USQmix(1,1)*VCKM(1,1)
     . +USQMIX(1,2)*VCKM(2,1)+USQMIX(1,3)*VCKM(3,1))*Dsqmix(k,1)
     . +(USQmix(1,1)*VCKM(1,2)
     . +USQMIX(1,2)*VCKM(2,2)+USQMIX(1,3)*VCKM(3,2))*Dsqmix(k,2)
     . +(USQmix(1,1)*VCKM(1,3)
     . +USQMIX(1,2)*VCKM(2,3)+USQMIX(1,3)*VCKM(3,3))*Dsqmix(k,3))

 
      enddo


      
      if(amw+amb+amneut(1).gt.amsupq(1))then
      print*, "warning in 3 body decay. No onshell W possible"
      width_bchiW=0d0
      width_jchiW=0d0
      width_tot3bod=0d0
      return
      endif

!       if(amsupq(1)-amneut(1).gt.amt)then
!       width_bchiW=0d0
!       width_jchiW=0d0
!       width_tot3bod=0d0
!       print*,"mstop>mt+mneut, decay stop-> top neutralino not implemente
!      .d for FV"
!       stop
!       endif 

  
c------ check whether stop is NLSP
	if(amsupq(1).gt.amsdownq(1)) then
	print*, "Stop not NLSP"
        stop
        endif
        if(amsupq(1).gt.amslepton(1))then
        print*, "Stop not NLSP"
        stop
        endif
       if(amsupq(1).gt.amsneutrino(1))then
        print*, "Stop not NLSP"
        stop
        endif
        if(amsupq(1).gt.amchar(1))then
        print*, "Stop not NLSP"
        stop
        endif
	if(amsupq(1).gt.amneut(2))then
        print*, "Stop not NLSP"
        stop
        endif

       xmu1=amb**2/amsupq(1)**2
       xmu2=amneut(1)**2/amsupq(1)**2
       xmu3=amw**2/amsupq(1)**2
  
      nx1t  = 128
      ny1t  = 128

       call SD_integ2(SD_3bodbchiW,SD_ax,SD_bx,SD_ay,SD_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum1)

       

      width_bchiw=1d0/32.D0/(2.D0*pi)**3*amsupq(1)*
     .                            sum1




        xmu1=0d0

             call SD_integ2(SD_3bodjchiW,SD_ax,SD_bx,SD_ay,SD_by,xmu1,
     .                     xmu2,xmu3,nx1t,ny1t,sum2)

      width_jchiw=1d0/32.D0/(2.D0*pi)**3*amsupq(1)*
     .                            sum2



      width_tot3bod=width_bchiW+width_jchiw


      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      matrix element squared for stop-->b chi_0 W                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function SD_3bodbchiW(x1, x2)
      implicit none
      integer i, j, k, l
      integer ifavvio
      double precision x1, x2, x3, y1, y2, y3
      double precision sdgf,amz,amw,pi,g2
      double precision amt, amb, amtau
      double precision amneut(4),xmneut(4),amchar(2),xmchar(2)
      double precision amsupq(6), amsdownq(6), amslepton(6),
     . amsneutrino(3)
      double precision sbotexchange, stbchiwtt, stbchiwchichi, 
     . stbchiwchib, stbchiwbt, stbchiwchit
      double precision xmuw, xmut(3), xmub, xmuneut, xmuchar(2),
     . xmusb(6)
      double precision DcharL(3,6,2), EcharR(3,6,2), 
     . cchaneutL(4,2), cchaneutR(4,2), ccharl(3,6,2), fcharr(3,6,2)
      double precision gneutul(3,4,6), 
     . gneutur(3,4,6), gneutdl(3,4,6), gneutdr(3,4,6)
      double precision gsqsqW(6)
      double precision dsb(6), dt(3), dchi(2)
      double precision mfu(3)
      double precision vckm(3,3),msq2(3,3), msu2(3,3), 
     .msd2(3,3), td(3,3), tu(3,3),
     . usqmix(6,6), dsqmix(6,6)
      double precision g2ew, g1ew, alsew
      double precision amwp, amzp
      double precision ratiotopcharg
c----- ramona changed 2/2/15
      double precision widthmfu(3), topwidth
c----- end ramona changed
      
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      COMMON/SD_fermion/amt,amb,amtau
      COMMON/SD_massgino/amneut,xmneut,amchar,xmchar
      COMMON/msfermion/ amsupq, amsdownq, amslepton, amsneutrino
      common/charcoup/DcharL, EcharR, FcharR, CCharL
      common/neutcoup/Gneutul, Gneutur, Gneutdr, Gneutdl
      common/Wcharneutcoup/CchaneutL, CchaneutR
      common/squarkWcoup/gsqsqW
      COMMON/flavviolation/vckm, msq2, msd2, msu2, td, 
     .tu, usqmix, ifavvio, dsqmix
      COMMON/SD_coupewsb/alsew,g2ew,g1ew
      COMMON/SD_mwmzpole/amwp,amzp
      common/ratio3bod/ratiotopcharg
c---- ramona cahnged 2/2/15
      COMMON/widthtop/topwidth
c----- end ramona changed

      mfu(1)=0d0
      mfu(2)=0d0
      mfu(3)=amt

c---- ramona changed 2/2/15
      widthmfu(1)=0d0
      widthmfu(2)=0d0
      widthmfu(3)=topwidth
c---- end ramona changed
      
      xmuw       = amw**2/amsupq(1)**2
      xmub       = amb**2/amsupq(1)**2
      xmuneut    = amneut(1)**2/amsupq(1)**2
      xmuchar(1) = amchar(1)**2/amsupq(1)**2
      xmuchar(2) = amchar(2)**2/amsupq(1)**2

      x3=2.D0-x1-x2
      y1=(1.D0+xmuw-xmuneut-xmub-x3)/2.D0
      y2=(1.D0-xmuw+xmuneut-xmub-x2)/2.D0
      y3=(1.D0-xmuw-xmuneut+xmub-x1)/2.D0

      do i=1,6
      xmusb(i)     = amsdownq(i)**2/amsupq(1)**2
      dsb(i)       = 1.D0-x3+xmuw-xmusb(i)
      enddo

      do i= 1, 3
      dt(i)=1.D0-x2+xmuneut-mfu(i)**2/amsupq(1)**2
      xmut(i)       = mfu(i)**2/amsupq(1)**2
      enddo

      do i=1,2
      dchi(i) = 1.D0-x1+xmub-xmuchar(i)
      enddo

  

c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
	
      sbotexchange=0.D0

c------ note that compared to SD_stbchiw W squark coupling defined differently

      do k=1,6
         do l=1,6
            sbotexchange=sbotexchange+gsqsqW(k)*gsqsqw(l)/dsb(k)/dsb(l)*
     .           8.D0*(
     .           (gneutdl(3,1,k)*gneutdl(3,1,l)
     .           +gneutdr(3,1,k)*gneutdr(3,1,l))*
     .           (-y1*(2.D0*y1+xmuneut+xmub)+y1/xmuw*(y2+y3)**2) +
     .           dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           (gneutdl(3,1,k)*gneutdr(3,1,l)
     .           +gneutdL(3,1,l)*gneutdr(3,1,k))*
     .           (2.D0*y1+xmub+xmuneut-1.D0/xmuw*(y2+y3)**2) )

         enddo
      enddo


c -------------------------------------------------------------------- c
c                             top exchange
c -------------------------------------------------------------------- c

      stbchiwtt=0d0

c---- g2e**2 needs to be added compared to SD_stbchiw for W-boson coupling
      do k=1,3
      do i=1,3
      stbchiwtt=stbchiwtt+g2ew**2*
c---- ramona changed 2/2/15 (topwidth can just be added like this since mfu(1/2)=0)
     .     1d0/2d0/((1.D0-x2+xmuneut-mfu(i)**2/amsupq(1)**2)*
     .     (1.D0-x2+xmuneut-mfu(k)**2/amsupq(1)**2)
     .     +mfu(k)*mfu(i)/amsupq(1)**4*
     .     topwidth**2)*vckm(i,3)*vckm(k,3)
c---- end ramona changed
     .     *((dsqrt(xmut(i))/2d0*gneutur(i,1,1)*gneutul(k,1,1)
     .     +dsqrt(xmut(k))/2d0*gneutur(k,1,1)*gneutul(i,1,1))* 
     .     xmneut(1)/amsupq(1)*
     .     (-4.D0)*(xmub+3.D0*y2+2.D0*y2**2/xmuw) +
     .     2.D0*gneutur(i,1,1)*gneutur(k,1,1)*dsqrt(xmut(i)*xmut(k))
     .     *(y1+2.D0*y2*y3/xmuw) +
     .     2.D0*gneutul(i,1,1)*gneutul(k,1,1)
     .     *(y1*(xmub-xmuw+4.D0*y2)+2.D0*y3*xmub+
     .     4.D0*y2*y3+1.D0/xmuw*(4.D0*y2**2*y1-2.D0*y2*y3*xmub)) )

      enddo
      enddo

c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c

      stbchiwchichi=0.D0

      do k=1,2
         do l=1,2
            stbchiwchichi=stbchiwchichi+2.D0/dchi(k)/dchi(l)*(
     .           (DcharL(3,k, 1)*DcharL(3,l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l)+
     .            EcharR(3,k, 1)*EcharR(3,l, 1)
     .            *cchaneutR(1,k)*cchaneutR(1,l))*
     .           (4.D0*y3*(y1+y2+y1*y3/xmuw)+y1*(xmuneut-xmuw)
     .            +2.D0*xmuneut*y2*(1-y3/xmuw)) +
     .           xmchar(k)*xmchar(l)/amsupq(1)**2*
     .           (DcharL(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,k)*cchaneutR(1,l)+
     .            EcharR(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l))*
     .           (y1+2.D0/xmuw*y2*y3) +
     .           (-3.D0)*(y1+y2)*xmneut(1)/amsupq(1)*(
     .           (DcharL(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,k)*cchaneutL(1,l)+
     .            EcharR(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,k)*cchaneutR(1,l))*
     .           xmchar(k)/amsupq(1) +
     .           (DcharL(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,l)*cchaneutL(1,k)+
     .            EcharR(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,l)*cchaneutR(1,k))*
     .           xmchar(l)/amsupq(1) ) +
     .           (-1.D0)*(xmuneut+2.D0*y3**2/xmuw+3.D0*y3)*
     .           dsqrt(xmub)*(
     .           (DcharL(3, k,1)*EcharR(3, l,1)
     .           *cchaneutR(1,k)*cchaneutR(1,l)+
     .            EcharR(3, k,1)*DcharL(3, l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l))*
     .           xmchar(k)/amsupq(1) + 
     .           (DcharL(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,l)*cchaneutL(1,k)+
     .            EcharR(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,l)*cchaneutR(1,k))*
     .           xmchar(l)/amsupq(1) ) +
     .           3.D0*dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           xmchar(k)/amsupq(1)*xmchar(l)/amsupq(1)*
     .           (DcharL(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,l)*cchaneutR(1,k)+
     .            EcharR(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,l)*cchaneutL(1,k)) +
     .           3.D0*dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           (DcharL(3, k,1)*EcharR(3, l,1)
     .            *cchaneutL(1,k)*cchaneutR(1,l)+
     .            EcharR(3, k,1)*DcharL(3, l,1)
     .            *cchaneutR(1,k)*cchaneutL(1,l))*
     .           (xmuw+xmuneut+2.D0*y3) )


         enddo
      enddo
   
       
c -------------------------------------------------------------------- c
c                    chargino sbottom interference
c -------------------------------------------------------------------- c

      stbchiwchib=0.D0
c---- note that the definitions need to be matched with sdecay---> - sign in
c---  and Sqrt(2) different in squark squark W couplings
c---- a g2ew needs to be added for W boson fermion coupling   
      do k=1,2
         do i=1,6
            stbchiwchib=stbchiwchib-4.D0*gsqsqW(i)/
     .           dsb(i)/dchi(k)*( 
     .           xmchar(k)/amsupq(1)*xmneut(1)/amsupq(1)*
     .           (DcharL(3, k,1)*gneutdl(3,1,1)*cchaneutR(1,k)+
     .            EcharR(3, k,1)*gneutdr(3,1,1)*cchaneutL(1,k))*
     .           (y1-y2/xmuw*(y2+y3)+xmub) +
     .           (DcharL(3, k,1)*gneutdl(3,1,1)*cchaneutL(1,k)+
     .            EcharR(3, k,1)*gneutdr(3,1,1)*cchaneutR(1,k))*
     .           ((y2+y3)*(xmuneut*y2-2.D0*y1*y3)/xmuw+y1*(2.D0*y1+y2
     .            -y3+xmuneut)+xmuneut*y2-xmub*(xmuneut+y3)) +
     .           dsqrt(xmub)*( xmneut(1)/amsupq(1)*
     .           (gneutdl(3,1,1)*EcharR(3, k,1)*cchaneutR(1,k)
     .           +gneutdr(3,1,1)*DcharL(3, k,1)*cchaneutL(1,k)) +
     .           (gneutdl(3,1,1)*EcharR(3, k,1)*cchaneutL(1,k)
     .           +gneutdr(3,1,1)*DcharL(3, k,1)*cchaneutR(1,k))
     .           *xmchar(k)/amsupq(1))*
     .           (1.D0/xmuw*y3*(y2+y3)-xmuneut-y1) )
         end do
      end do

      
c -------------------------------------------------------------------- c
c                       top sbottom interference
c -------------------------------------------------------------------- c

      stbchiwbt=0.d0
c---- note that the definitions need to be matched with sdecay---> - sign in
c---  and Sqrt(2) different in squark squark W couplings
c---- a g2ew needs to be added for W boson fermion coupling   
      do i=1,6
        do j=1,3
         stbchiwbt=stbchiwbt-4.D0*g2ew
     .        *gsqsqw(i)/dsqrt(2.D0)/dt(j)
     .        /dsb(i)*vckm(j,3)*(dsqrt(xmut(j))*xmneut(1)/amsupq(1)
     .        *gneutuR(j,1,1)*gneutdL(3,1,i)*
     .        (y1+xmub-y2/xmuw*(y2+y3)) +
     .        gneutuL(j,1,1)*gneutdL(3,1,i)
     .        *(y1*y2*(1.D0+2.D0*(y2+y3)/xmuw)+
     .        xmuneut*y2-y1*y3-2.D0*y1**2-y1*xmub+xmub*(xmuneut-y3)+
     .        1.D0/xmuw*(-xmub*y2*y3-xmub*y3**2)) +
     .        dsqrt(xmub)*dsqrt(xmut(j))
     .        *gneutuR(j,1,1)*gneutdR(3,1,i)*(-y1
     .        -xmuneut+1.D0/xmuw*y3*(y2+y3)) +
     .        dsqrt(xmub)*xmneut(1)/amsupq(1)
     .        *gneutuL(j,1,1)*gneutdR(3,1,i)*
     .        (y1+xmub-1.D0/xmuw*y2*(y2+y3)) )
      end do
      end do

 

c -------------------------------------------------------------------- c
c                      chargino top interference
c -------------------------------------------------------------------- c

      stbchiwchit=0.D0

      do i=1,2,1
       do j=1,3
         stbchiwchit=stbchiwchit+g2ew/sqrt(2d0)/dt(j)/dchi(i)*vckm(j,3)
     .        *(xmchar(i)/amsupq(1)*xmneut(1)/amsupq(1)*
     .        gneutuL(j,1,1)*DcharL(3, i,1)*cchaneutR(1,i)*
     .        (-6.D0*y2-4.D0*y2**2/xmuw-2.D0*xmub) +
     .        gneutuL(j,1,1)*DcharL(3, i,1)*cchaneutL(1,i)*2.D0*(
     .        y1*(2.D0*y3+2.D0*y2+4.D0*y1-xmuw)+y2*(4.D0*y3+xmuneut)
     .        -2.D0*y2*(2.D0*y1*y3-xmuneut*y2)/xmuw+
     .        2.D0*xmub/xmuw*y3**2-xmub*xmuneut+xmub*y3) +
     .        dsqrt(xmut(j))*xmneut(1)/amsupq(1)*
     .        gneutuR(j,1,1)*DcharL(3, i,1)*cchaneutL(1,i)*
     .        (-6.D0)*(y1+y2) +
     .        dsqrt(xmut(j))*xmchar(i)/amsupq(1)
     .        *DcharL(3, i,1)*gneutuR(j,1,1)*
     .        cchaneutR(1,i)*(2.D0*y1+4.D0/xmuw*y2*y3) +
     .        dsqrt(xmub*xmut(j))*gneutuR(j,1,1)
     .        *EcharR(3, i,1)*cchaneutR(1,i)*
     .        (-6.D0*y3-4.D0/xmuw*y3**2-2.D0*xmuneut) +
     .        6.D0*dsqrt(xmut(j)*xmub)*xmchar(i)/amsupq(1)*xmneut(1)
     .        /amsupq(1)*gneutuR(j,1,1)*EcharR(3, i,1)*cchaneutL(1,i) +
     .        dsqrt(xmub)*xmchar(i)/amsupq(1)*gneutuL(j,1,1)*
     .        EcharR(3, i,1)*cchaneutL(1,i)*(-6.D0)*(y1+y3) +
     .        dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .        gneutuL(j,1,1)*EcharR(3, i,1)*cchaneutR(1,i)
     .        *( 6.D0*xmuw+6.D0*y3+
     .        6.D0*y2+2.D0*y1+4.D0/xmuw*y2*y3) )
       end do
      end do
      

      SD_3bodbchiw=stbchiwtt+sbotexchange+stbchiwchichi
     . +2d0*stbchiwbt+2d0*stbchiwchit+2d0*stbchiwchib



      ratiotopcharg=stbchiwtt/stbchiwchichi

      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      matrix element squared for stop-->j chi_0 W                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function SD_3bodjchiW(x1, x2)
      implicit none
      integer i, j, k, l
      integer ifavvio
      double precision x1, x2, x3, y1, y2, y3
      double precision sdgf,amz,amw,pi,g2
      double precision amt, amb, amtau
      double precision amneut(4),xmneut(4),amchar(2),xmchar(2)
      double precision amsupq(6), amsdownq(6), amslepton(6),
     . amsneutrino(3)
      double precision sbotexchangej, stjchiwtt, stjchiwchichi, 
     . stjchiwchib, stjchiwbt, stjchiwchit
      double precision xmuw, xmut(3), xmub, xmuneut, xmuchar(2),
     . xmusb(6)
      double precision DcharL(3,6,2), EcharR(3,6,2), 
     . cchaneutL(4,2), cchaneutR(4,2), ccharl(3,6,2), fcharr(3,6,2)
      double precision gneutul(3,4,6), 
     . gneutur(3,4,6), gneutdl(3,4,6), gneutdr(3,4,6)
      double precision gsqsqW(6)
      double precision dsb(6), dt(3), dchi(2)
      double precision mfu(3)
      double precision vckm(3,3),msq2(3,3), msu2(3,3), 
     .msd2(3,3), td(3,3), tu(3,3),
     . usqmix(6,6), dsqmix(6,6)
      double precision alsew,g2ew,g1ew
      double precision amwp, amzp
c----- ramona changed 2/2/15
      double precision widthmfu(3), topwidth
c----- end ramona changed

      COMMON/SD_param/sdgf,amz,amw,pi,g2
      COMMON/SD_fermion/amt,amb,amtau
      COMMON/SD_massgino/amneut,xmneut,amchar,xmchar
      COMMON/msfermion/ amsupq, amsdownq, amslepton, amsneutrino
      common/charcoup/DcharL, EcharR, FcharR, CCharL
      common/neutcoup/Gneutul, Gneutur, Gneutdr, Gneutdl
      common/Wcharneutcoup/CchaneutL, CchaneutR
      common/squarkWcoup/gsqsqW
      COMMON/flavviolation/vckm, msq2, msd2, msu2, td, 
     .tu, usqmix, ifavvio, dsqmix
      COMMON/SD_coupewsb/alsew,g2ew,g1ew
      COMMON/SD_mwmzpole/amwp,amzp
c---- ramona cahnged 2/2/15
      COMMON/widthtop/topwidth
c----- end ramona changed

      mfu(1)=0d0
      mfu(2)=0d0
      mfu(3)=amt
      
      xmuw       = amw**2/amsupq(1)**2

      xmub       = 0d0
      xmuneut    = amneut(1)**2/amsupq(1)**2
      xmuchar(1) = amchar(1)**2/amsupq(1)**2
      xmuchar(2) = amchar(2)**2/amsupq(1)**2

      x3=2.D0-x1-x2
      y1=(1.D0+xmuw-xmuneut-xmub-x3)/2.D0
      y2=(1.D0-xmuw+xmuneut-xmub-x2)/2.D0
      y3=(1.D0-xmuw-xmuneut+xmub-x1)/2.D0

      do i=1,6
      xmusb(i)     = amsdownq(i)**2/amsupq(1)**2
      dsb(i)       = 1.D0-x3+xmuw-xmusb(i)
      enddo

      do i= 1, 3
      dt(i)=1.D0-x2+xmuneut-mfu(i)**2/amsupq(1)**2
      enddo

      do i=1,3
      xmut(i)       = mfu(i)**2/amsupq(1)**2
      enddo

      do i=1,2
      dchi(i) = 1.D0-x1+xmub-xmuchar(i)
      enddo

  

c -------------------------------------------------------------------- c
c                           sbottom exchange
c -------------------------------------------------------------------- c
	
      sbotexchangej=0.D0

      do k=1,6
       do l=1,6
        do j=1,2
         
            sbotexchangej=sbotexchangej+gsqsqW(k)*gsqsqw(l)/dsb(k)
     .           /dsb(l)*8.D0*(
     .           (gneutdl(j,1,k)*gneutdl(j,1,l)
     .           +gneutdr(j,1,k)*gneutdr(j,1,l))*
     .           (-y1*(2.D0*y1+xmuneut+xmub)+y1/xmuw*(y2+y3)**2) +
     .           dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           (gneutdl(j,1,k)*gneutdr(j,1,l)
     .           +gneutdL(j,1,l)*gneutdr(j,1,k))*
     .           (2.D0*y1+xmub+xmuneut-1.D0/xmuw*(y2+y3)**2) )
        enddo
       enddo
      enddo



 

c -------------------------------------------------------------------- c
c                             top exchange
c -------------------------------------------------------------------- c

      stjchiwtt=0d0

      
      do i=1,3
      do k=1,3
       do j=1,2
      stjchiwtt=stjchiwtt+g2ew**2
c---- ramona changed 2/2/15 (topwidth can just be added like this since mfu(1/2)=0)
     .     *1d0/2d0/((1.D0-x2+xmuneut-mfu(i)**2/amsupq(1)**2)*
     .     (1.D0-x2+xmuneut-mfu(k)**2/amsupq(1)**2)
     .     +mfu(i)*mfu(k)/amsupq(1)**4*topwidth**2)
c---- end ramona changed
     .     *vckm(i,j)*vckm(k,j)*
     .     ((dsqrt(xmut(i))/2d0*gneutur(i,1,1)*gneutul(k,1,1)
     .     +dsqrt(xmut(k))/2d0*gneutur(k,1,1)*gneutul(i,1,1))* 
     .     xmneut(1)/amsupq(1)*
     .     (-4.D0)*(xmub+3.D0*y2+2.D0*y2**2/xmuw) +
     .     2.D0*gneutur(i,1,1)*gneutur(k,1,1)*dsqrt(xmut(i)*xmut(k))
     .     *(y1+2.D0*y2*y3/xmuw) +
     .     2.D0*gneutul(i,1,1)*gneutul(k,1,1)*
     .     (y1*(xmub-xmuw+4.D0*y2)+2.D0*y3*xmub+
     .     4.D0*y2*y3+1.D0/xmuw*(4.D0*y2**2*y1-2.D0*y2*y3*xmub)) )

	
      
       enddo
      enddo
      enddo


c -------------------------------------------------------------------- c
c                           chargino exchange
c -------------------------------------------------------------------- c

      stjchiwchichi=0.D0

      do k=1,2
       do l=1,2
        do j=1,2
            stjchiwchichi=stjchiwchichi+2.D0/dchi(k)/dchi(l)*(
     .           (DcharL(j,k, 1)*DcharL(j,l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l)+
     .            EcharR(j,k, 1)*EcharR(j,l, 1)
     .            *cchaneutR(1,k)*cchaneutR(1,l))*
     .           (4.D0*y3*(y1+y2+y1*y3/xmuw)+y1*(xmuneut-xmuw)
     .            +2.D0*xmuneut*y2*(1-y3/xmuw)) +
     .           xmchar(k)*xmchar(l)/amsupq(1)**2*
     .           (DcharL(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,k)*cchaneutR(1,l)+
     .            EcharR(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l))*
     .           (y1+2.D0/xmuw*y2*y3) +
     .           (-3.D0)*(y1+y2)*xmneut(1)/amsupq(1)*(
     .           (DcharL(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,k)*cchaneutL(1,l)+
     .            EcharR(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,k)*cchaneutR(1,l))*
     .           xmchar(k)/amsupq(1) +
     .           (DcharL(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,l)*cchaneutL(1,k)+
     .            EcharR(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,l)*cchaneutR(1,k))*
     .           xmchar(l)/amsupq(1) ) +
     .           (-1.D0)*(xmuneut+2.D0*y3**2/xmuw+3.D0*y3)*
     .           dsqrt(xmub)*(
     .           (DcharL(j, k,1)*EcharR(j, l,1)
     .           *cchaneutR(1,k)*cchaneutR(1,l)+
     .            EcharR(j, k,1)*DcharL(j, l,1)
     .            *cchaneutL(1,k)*cchaneutL(1,l))*
     .           xmchar(k)/amsupq(1) + 
     .           (DcharL(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,l)*cchaneutL(1,k)+
     .            EcharR(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,l)*cchaneutR(1,k))*
     .           xmchar(l)/amsupq(1) ) +
     .           3.D0*dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           xmchar(k)/amsupq(1)*xmchar(l)/amsupq(1)*
     .           (DcharL(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,l)*cchaneutR(1,k)+
     .            EcharR(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,l)*cchaneutL(1,k)) +
     .           3.D0*dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .           (DcharL(j, k,1)*EcharR(j, l,1)
     .            *cchaneutL(1,k)*cchaneutR(1,l)+
     .            EcharR(j, k,1)*DcharL(j, l,1)
     .            *cchaneutR(1,k)*cchaneutL(1,l))*
     .           (xmuw+xmuneut+2.D0*y3) )
        enddo
       enddo
      enddo

c -------------------------------------------------------------------- c
c                    chargino sbottom interference
c -------------------------------------------------------------------- c

      stjchiwchib=0.D0

      do k=1,2
       do i=1,6
        do j=1,2
            stjchiwchib=stjchiwchib-4.D0*gsqsqW(i)/
     .           dsb(i)/dchi(k)*( 
     .           xmchar(k)/amsupq(1)*xmneut(1)/amsupq(1)*
     .           (DcharL(j, k,1)*gneutdl(j,1,1)*cchaneutR(1,k)+
     .            EcharR(j, k,1)*gneutdr(j,1,1)*cchaneutL(1,k))*
     .           (y1-y2/xmuw*(y2+y3)+xmub) +
     .           (DcharL(j, k,1)*gneutdl(j,1,1)*cchaneutL(1,k)+
     .            EcharR(j, k,1)*gneutdr(j,1,1)*cchaneutR(1,k))*
     .           ((y2+y3)*(xmuneut*y2-2.D0*y1*y3)/xmuw+y1*(2.D0*y1+y2
     .            -y3+xmuneut)+xmuneut*y2-xmub*(xmuneut+y3)) +
     .           dsqrt(xmub)*( xmneut(1)/amsupq(1)*
     .           (gneutdl(j,1,1)*EcharR(j, k,1)*cchaneutR(1,k)
     .           +gneutdr(j,1,1)*DcharL(j, k,1)*cchaneutL(1,k)) +
     .           (gneutdl(j,1,1)*EcharR(j, k,1)*cchaneutL(1,k)
     .           +gneutdr(j,1,1)*DcharL(j, k,1)*cchaneutR(1,k))
     .           *xmchar(k)/amsupq(1))*
     .           (1.D0/xmuw*y3*(y2+y3)-xmuneut-y1) )
        end do
       end do
      end do


c -------------------------------------------------------------------- c
c                       top sbottom interference
c -------------------------------------------------------------------- c

      stjchiwbt=0.d0
        
      do i=1,6
        do j=1,3
         do k=1,2
         stjchiwbt=stjchiwbt-g2ew*4.D0/dsqrt(2d0)*Vckm(j,k)
     .        *gsqsqw(i)/dt(j)
     .        /dsb(i)*(dsqrt(xmut(j))*xmneut(1)/amsupq(1)
     .        *gneutuR(j,1,1)*gneutdL(k,1,i)*
     .        (y1+xmub-y2/xmuw*(y2+y3)) +
     .        gneutuL(j,1,1)*gneutdL(k,1,i)
     .        *(y1*y2*(1.D0+2.D0*(y2+y3)/xmuw)+
     .        xmuneut*y2-y1*y3-2.D0*y1**2-y1*xmub+xmub*(xmuneut-y3)+
     .        1.D0/xmuw*(-xmub*y2*y3-xmub*y3**2)) +
     .        dsqrt(xmub)*dsqrt(xmut(j))
     .        *gneutuR(j,1,1)*gneutdR(k,1,i)*(-y1
     .        -xmuneut+1.D0/xmuw*y3*(y2+y3)) +
     .        dsqrt(xmub)*xmneut(1)/amsupq(1)
     .        *gneutuL(j,1,1)*gneutdR(k,1,i)*
     .        (y1+xmub-1.D0/xmuw*y2*(y2+y3)) )
        end do
       end do
      end do

c -------------------------------------------------------------------- c
c                      chargino top interference
c -------------------------------------------------------------------- c

      stjchiwchit=0.D0

      do i=1,2,1
       do j=1,3
        do k=1,2
         stjchiwchit=stjchiwchit+g2ew/sqrt(2d0)/dt(j)/dchi(i)*vckm(j,k) 
     .        *(xmchar(i)/amsupq(1)*xmneut(1)/amsupq(1)*
     .        gneutuL(j,1,1)*DcharL(k, i,1)*cchaneutR(1,i)*
     .        (-6.D0*y2-4.D0*y2**2/xmuw-2.D0*xmub) +
     .        gneutuL(j,1,1)*DcharL(k, i,1)*cchaneutL(1,i)*2.D0*(
     .        y1*(2.D0*y3+2.D0*y2+4.D0*y1-xmuw)+y2*(4.D0*y3+xmuneut)
     .        -2.D0*y2*(2.D0*y1*y3-xmuneut*y2)/xmuw+
     .        2.D0*xmub/xmuw*y3**2-xmub*xmuneut+xmub*y3) +
     .        dsqrt(xmut(j))*xmneut(1)/amsupq(1)*
     .        gneutuR(j,1,1)*DcharL(k, i,1)*cchaneutL(1,i)*
     .        (-6.D0)*(y1+y2) +
     .        dsqrt(xmut(j))*xmchar(i)/amsupq(1)
     .        *DcharL(k, i,1)*gneutuR(j,1,1)*
     .        cchaneutR(1,i)*(2.D0*y1+4.D0/xmuw*y2*y3) +
     .        dsqrt(xmub*xmut(j))*gneutuR(j,1,1)
     .        *EcharR(k, i,1)*cchaneutR(1,i)*
     .        (-6.D0*y3-4.D0/xmuw*y3**2-2.D0*xmuneut) +
     .        6.D0*dsqrt(xmut(j)*xmub)*xmchar(i)/amsupq(1)*xmneut(1)
     .        /amsupq(1)*gneutuR(j,1,1)*EcharR(k, i,1)*cchaneutL(1,i) +
     .        dsqrt(xmub)*xmchar(i)/amsupq(1)*gneutuL(j,1,1)*
     .        EcharR(k, i,1)*cchaneutL(1,i)*(-6.D0)*(y1+y3) +
     .        dsqrt(xmub)*xmneut(1)/amsupq(1)*
     .        gneutuL(j,1,1)*EcharR(k, i,1)*cchaneutR(1,i)
     .        *( 6.D0*xmuw+6.D0*y3+
     .        6.D0*y2+2.D0*y1+4.D0/xmuw*y2*y3) )
        end do
       end do
      end do
      

      SD_3bodjchiw=stjchiwtt+stjchiwchichi+sbotexchangej
     .    +2d0*stjchiwchit+2d0*stjchiwbt+2d0*stjchiwchib


      end
