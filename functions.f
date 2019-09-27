*****************************************************************
*****************************************************************
*****************************************************************
      double precision function fun(s)
!     compute number of BH mergers with eccentricity <e that come from binaries
!     produced at  t<td
      USE commonV
      USE modelPar
      USE inputPar
      implicit none            
      real*8 s,aej,agw,agw2,c,eccen,X0,X1,Edot,elcap,epsilon,F,fbh,fL
      real*8 g, Gamma, K1,K2,m1,m2,m3,Mb,Mbh,Mcl,mej,Ncap,Nej,Ngw,Nis,P
      real*8 Pcap,Pej,Pex,Pgw,pi,Pin,psi,q3,rh,Rs,t,ut,vesc,x,trh
      real*8 funm,SevalDouble
      common/ecc/eccen,X0,X1
      pi=2.d0*ASIN(1.d0)      
      ut=58.*24.*60.*60.
      c=1.d4                    !sol
      t=tf-s
      rh=SevalDouble(t/1.d9,tabt,tabr,c7,c8,c9,cont)
      Mcl=SevalDouble(t/1.d9,tabt,tabm,c4,c5,c6,cont)
      Mbh=max(SevalDouble(t/1.d9,tabt,tabmbh,c1,c2,c3,cont),0.)      
      vesc=3.69d-3*sqrt(Mcl/rh) !for fc=1, and in M_sun, 1AU, G=1      

!     BH mass      
      IMFBH: select case(IMF)
      case (1)
         inp: select case(input)
         case (1)
         m1=SevalDouble(t/1.d9,tabt,tabmax,c10,c11,c12,cont)
         case (2)
         m1=SevalDouble(Mbh,tabmbh_p,tabmax,c10,c11,c12,cont2)
         end select inp
      case (0)
         P=max(MBH/SevalDouble(0.d0,tabt,tabmbh,c1,c2,c3,cont),0.)
         m1=funm(P)             !deterime the mass of the BHs
      case DEFAULT
         write(*,*)"value of IMF not valid"
         stop
      end select IMFBH 
      
      
!     compute Edot     
      fbh=Mbh/Mcl
      psi=1.+a1*fbh/1.d-2        !relaxation alla Antonini+Guiles 2019
      trh=2.06d5*sqrt(Mcl)*rh**1.5/psi/m_mean
      Edot=1.53d-7*csi*(Mcl**2/rh)/trh !convert to M_sun, 1AU, G=1
    
!     binary properties: assume all BHs have same mass!
      m2=m1
      m3=m1
      Mb=m1+m2
      q3=m3/Mb 
      aej=dE*m1*m2/(m1+m2+m3)*q3/vesc**2
      Rs=4.*(m1+m2)/c**2         
      
      epsilon=1./(1.+dE)      
      g=1.3*((m1*m2)**2*Mb/c**5/Edot)**(1./7.)    
      agw2=(g**2*7./10./(1.-epsilon))**(7./10.)
      Nis=20.                   !number of intermediate resonant states
      K2=(1.4*Rs**(5./14.))**2*Nis*7./5./(1-epsilon)     
      K1=g**2*7./10./(1-epsilon)          
      agw=((-K2+sqrt(K2**2+4.*K1))/(2.*K1))**(-7./5.)         
      agw=max(agw,aej)
      
      Pgw=7./10./(1.-epsilon)*g**2/agw**(10./7.)
      Pcap=K2/agw**(5./7.)
     
      fL=freq*ut             
      x=eccen
      F=1.*(1./fL)**(2./3.)*(1./x)**(12./19.)
     &     *(1./(1.+121./304.*x**2))**(870./2299.)*(1.+x)**(0.7969)
      Ngw=max((Pgw*(1.-epsilon)-(Mb**(1./3.)/pi**(2./3.)
     &     /agw*F))/(Pgw*(1.-epsilon)),0.)      
      elcap=2.*(Rs/agw)**(5./7.)
      Ncap=max((elcap-(Mb**(1./3.)/pi**(2./3.)
     &     /agw*F)*5./7.)/elcap,0.) 
      
      mej=Mb+m3/(1.-epsilon)*max(log(aej/agw/q3**2),0.)  !mass ejected by binary      
      Gamma=(beta*Mcl/trh)/mej  !formation rate of binaries     
      
      fun=Gamma*(Pgw*Ngw+Pcap*Ncap)      
      if(Mbh.le.10.or.t/1d9.lt.tcc)fun=0.d0
      return
      end   

      double precision function funej(s)
!     compute number of BH mergers with eccentricity <e that come from binaries
!     that are ejected at t<td
      USE commonV
      USE modelPar
      USE inputPar
      implicit none            
      real*8 s,aej,agw,agw2,c,eccen,X0,X1,Edot,elcap,epsilon,F,fbh,fL
      real*8 g, Gamma, K1,K2,m1,m2,m3,Mb,Mbh,Mcl,mej,Ncap,Nej,Ngw,Nis,P
      real*8 Pcap,Pej,Pex,Pgw,pi,Pin,psi,q3,rh,Rs,t,ut,vesc,x,trh
      real*8 funm,SevalDouble
      common/ecc/eccen,X0,X1
      pi=2.d0*ASIN(1.d0)      
      ut=58.*24.*60.*60.
      c=1.d4                    !sol
      t=tf-s
      rh=SevalDouble(t/1.d9,tabt,tabr,c7,c8,c9,cont)
      Mcl=SevalDouble(t/1.d9,tabt,tabm,c4,c5,c6,cont)
      Mbh=max(SevalDouble(t/1.d9,tabt,tabmbh,c1,c2,c3,cont),0.)
      vesc=3.69d-3*sqrt(Mcl/rh) !for fc=0.1, and in M_sun, 1AU, G=1      

!     BH mass
      IMFBH: select case(IMF)
      case (1)
         inp: select case(input)
         case (1)
         m1=SevalDouble(t/1.d9,tabt,tabmax,c10,c11,c12,cont)
         case (2)
         m1=SevalDouble(Mbh,tabmbh_p,tabmax,c10,c11,c12,cont2)
         end select inp
      case (0)
         P=max(MBH/SevalDouble(0.d0,tabt,tabmbh,c1,c2,c3,cont),0.)
         m1=funm(P)             !deterime the mass of the BHs
      case DEFAULT
         write(*,*)"value of IMF not valid"
         stop
      end select IMFBH 

!     compute Edot     
      fbh=Mbh/Mcl
      psi=1.+a1*fbh/1.d-2        !relaxation alla Antonini+Guiles 2019
      trh=2.06d5*sqrt(Mcl)*rh**1.5/psi/m_mean

      Edot=1.53d-7*csi*(Mcl**2/rh)/trh !convert to M_sun, 1AU, G=1
    
!     binary properties
      m2=m1
      m3=m1
      Mb=m1+m2
      q3=m3/Mb 
      aej=dE*m1*m2/(m1+m2+m3)*q3/vesc**2
      Rs=4.*(m1+m2)/c**2         
      
      epsilon=1./(1.+dE)      
      g=1.3*((m1*m2)**2*Mb/c**5/Edot)**(1./7.)    
      agw2=(g**2*7./10./(1.-epsilon))**(7./10.)
      Nis=20.                   !number of intermediate resonant states
      K2=(1.4*Rs**(5./14.))**2*Nis*7./5./(1-epsilon)     
      K1=g**2*7./10./(1-epsilon)          
      agw=((-K2+sqrt(K2**2+4.*K1))/(2.*K1))**(-7./5.)         
      agw=max(agw,aej)
      
      Pgw=7./10./(1.-epsilon)*g**2/agw**(10./7.)
      Pcap=K2/agw**(5./7.)
      Pin=Pgw+Pcap
      Pex=min((0.7d-1/aej)**(8./7.)*(m1*m2*Mb*s/1.e10/1.e3
     &     )**(2./7.),1.)
      Pej=(1.-Pin)*Pex
      if(Pej.lt.0)Pej=0.
    
      fL=freq*ut             
      x=eccen
      F=1.*(1./fL)**(2./3.)*(1./x)**(12./19.)
     &     *(1./(1.+121./304.*x**2))**(870./2299.)*(1.+x)**(0.7969)     
      Nej=max((Pex-(Mb/pi**2)**(1./3.)*F/aej)/Pex,0.)
      
      mej=Mb+m3/(1.-epsilon)*max(log(aej/agw/q3**2),0.) !mass ejected by binary
          
      Gamma=(beta*Mcl/trh)/mej  !formation rate of binaries
      funej=Gamma*(Pej*Nej)
      if(Mbh.le.10.or.t/1d9.lt.tcc)funej=0.d0
      return
      end   
      
      
      double precision function funej2(s)
!     compute number of BH mergers with eccentricity <e that come from binaries
!     that are ejected at t>td but merger at t<td
      USE commonV
      USE modelPar
      USE inputPar
      implicit none            
      real*8 s,aej,agw,agw2,c,eccen,X0,X1,Edot,elcap,epsilon,F,fbh,fL
      real*8 g, Gamma, K1,K2,m1,m2,m3,Mb,Mbh,Mcl,mej,Ncap,Nej,Ngw,Nis,P
      real*8 Pcap,Pej,Pex,Pgw,pi,Pin,psi,q3,rh,Rs,t,ut,vesc,x,trh,F1,F2
      real*8 funm,SevalDouble,Pex_p
      common/ecc/eccen,X0,X1
      pi=2.d0*ASIN(1.d0)      
      ut=58.*24.*60.*60.
      c=1.d4                    !sol
      t=tf-s     
      rh=SevalDouble(t/1.d9,tabt,tabr,c7,c8,c9,cont)
      Mcl=SevalDouble(t/1.d9,tabt,tabm,c4,c5,c6,cont)
      Mbh=max(SevalDouble(t/1.d9,tabt,tabmbh,c1,c2,c3,cont),0.)      
      vesc=3.69d-3*sqrt(Mcl/rh) !for fc=0.1, and in M_sun, 1AU, G=1      
 
!     BH mass
      IMFBH: select case(IMF)
      case (1)
         inp: select case(input)
         case (1)
         m1=SevalDouble(t/1.d9,tabt,tabmax,c10,c11,c12,cont)
         case (2)
         m1=SevalDouble(Mbh,tabmbh_p,tabmax,c10,c11,c12,cont2)
         end select inp
      case (0)
         P=max(MBH/SevalDouble(0.d0,tabt,tabmbh,c1,c2,c3,cont),0.)
         m1=funm(P)             !deterime the mass of the BHs
      case DEFAULT
         write(*,*)"value of IMF not valid"
         stop
      end select IMFBH

      
!     compute Edot
      fbh=Mbh/Mcl
      psi=1.+a1*fbh/1.d-2        !relaxation alla Antonini+Guiles'19
      trh=2.06d5*sqrt(Mcl)*rh**1.5/psi/m_mean
      Edot=1.53d-7*csi*(Mcl**2/rh)/trh !convert to M_sun, 1AU, G=1

!     binary properties
      m2=m1
      m3=m1
      Mb=m1+m2
      q3=m3/Mb 
      aej=dE*m1*m2/(m1+m2+m3)*q3/vesc**2
      Rs=4.*(m1+m2)/c**2         
      
      epsilon=1./(1+dE)      
      g=1.3*((m1*m2)**2*Mb/c**5/Edot)**(1./7.)    
      agw2=(g**2*7./10./(1.-epsilon))**(7./10.)
      Nis=20.                   !number of intermediate resonant states
      K2=(1.4*Rs**(5./14.))**2*Nis*7./5./(1-epsilon)     
      K1=g**2*7./10./(1-epsilon)          
      agw=((-K2+sqrt(K2**2+4.*K1))/(2.*K1))**(-7./5.)         
      agw=max(agw,aej)
      
      Pgw=7./10./(1.-epsilon)*g**2/agw**(10./7.)
      Pcap=K2/agw**(5./7.)
      Pin=Pgw+Pcap
      Pex=min((0.7d-1/aej)**(8./7.)*(m1*m2*Mb*s/1.e10/1.e3
     &     )**(2./7.),1.)
      Pej=(1.-Pin)*Pex
      if(Pej.lt.0)Pej=0.

      fL=freq*ut             
      x=eccen
      F=1.*(1./fL)**(2./3.)*(1./x)**(12./19.)
     &     *(1./(1.+121./304.*x**2))**(870./2299.)*(1.+x)**(0.7969)

      Nej=max((Pex-(Mb/pi**2)**(1./3.)*F/aej)/Pex,0.)        
      F1=Pej*Nej                !all mergers at t'

      Pex=min((0.7d-1/aej)**(8./7.)*(m1*m2*Mb*(s-X0)/1.e10/1.e3
     &     )**(2./7.),1.)      
      Pej=(1.-Pin)*Pex          !all mergers with very short merger time      
      if(Pej.lt.0)Pej=0.
      
      Nej=max((Pex-(Mb/pi**2)**(1./3.)*F/aej)/Pex,0.)  
      F2=Pej*Nej

      mej=Mb+m3/(1.-epsilon)*max(log(aej/agw/q3**2),0.) !mass ejected by binary
      Gamma=(beta*Mcl/trh)/mej  !formation rate of binaries      
      funej2=(F1-F2)*Gamma

      if(Mbh.le.10.or.t/1.d9.lt.tcc.or.funej2.lt.0)funej2=0.d0
      return
      end   

      double precision function funm(P)
!     mass of the BHs obtained assuming the mass function is progressively depleated
!     from top to bottom and power law initial mass function
      USE inputPar
      implicit none
      real*8 mu,ml,P
      parameter ( mu = 33.d0 )	
      parameter ( ml = 3.d0 )
      funm=(P*(mu**(2.+a)-ml**(2.+a))+ml**(2.+a))**(1./(2.+a))
      return
      end

      subroutine Tbalanced(Mbh,Mcl,rh,tb)
      USE modelPar
      USE inputPar
      implicit none
      real*8 Mbh,Mcl,rh,tb,fbh,psi
      fbh=Mbh/Mcl
      psi=1.+a1*fbh/0.01
      tb=6.22d0*2.06d5*sqrt(Mcl)*(rh)**1.5
     &     /psi/m_mean/1.e9     !start of balanced evolution (alla Antonini+Gieles'19)

      if(tb.gt.tf/1.e9)then
         print*,
     &  'no mergers are produced because the core-collapse time is >t_f'
         stop
      end if

      
      return
      end subroutine
      
C---------------------------------------------------------------------
      SUBROUTINE RK78(IR,T,DT,X,N,TOL,DER)
C---------------------------------------------------------------------
C     Variable step-size automatic one-step integrator for a system of
C     N firts order ordinary differential equations with initial values
C     The Runge-Kutta-Fehlberg formula 7(8) is used
C     REF   E. FEHLBERG, NASA TECHNICAL REPORT TR R-287, 1968
C     Description of parameters list
C     (All floating variables in DOUBLE PRECISION)
C     IR      O    NUMBER OF REJECTIONS OF THE LAST STEP
C     (IN CASE  DT WAS TOO LARGE)
C     T       I-O  INDEPENDENT VARIABLE
C     DT      I-O  STEP SIZE
C     A RECOMMENDED VALUE FOR THE NEXT STEP IS OUTPUT
C     X(N)    I-O  DEPENDENT VARIABLES
C     Fm(N)        AUXILIARY ARRAYS   WITH m = 0 TO 6
C     F7(N)   O    ABSOLUTE ESTIMATED TRUNCATION ERROR ON EACH COMPONENT
C     N       I    ORDER OF THE DIFFERENTIAL EQUATIONS SYSTEM
C     TOL(N)  I    RELATIVE TOLERATED ERROR ON EACH COMPONENT
C     DER     I    NAME OF THE SUBROUTINE COMPUTING THE DERIVATIVES. THIS
C     SUBROUTINE HAS TO HAVE THE STANDARD CALLING SEQUENCE
C     CALL DER(T,X,F0)
C---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      PARAMETER ( CH1 = 34D0/105D0, CH2 = 9D0/35D0, CH3 = 9D0/280D0,
     &     CH4 = 41D0/840D0, AL2 = 2D0/27D0, AL3 = 1D0/9D0,
     &     AL4 = 1D0/6D0,    AL5 = 5D0/12D0, AL6 = 5D-1,
     &     AL7 = 5D0/6D0,    AL9 = 2D0/3D0,  ALA = 1D0/3D0,
     &     B21 = 2D0/27D0,     B31 = 1D0/36D0,    B41 = 1D0/24D0,
     &     B51 = 5D0/12D0,     B61 = 5D-2,        B71 = -25D0/108D0,
     &     B81 = 31D0/3D2,     B101= -91D0/108D0, B111= 2383D0/41D2,
     &     B121= 3D0/205D0,    B131= -1777D0/41D2,B32 = 1D0/12D0,
     &     B43 = .125D0,       B53 = -25D0/16D0,  B64 = 25D-2,
     &     B74 = 125D0/108D0,  B94 = -53D0/6D0,   B104= 23D0/108D0,
     &     B114= -341D0/164D0, B65 = 2D-1,        B75 = -65D0/27D0,
     &     B85 = 61D0/225D0,   B95 = 704D0/45D0,  B105= -976D0/135D0,
     &     B115= 4496D0/1025D0,B76 = 125D0/54D0,  B86 = -2D0/9D0,
     &     B96 = -107D0/9D0,   B106= 311D0/54D0,  B116= -301D0/82D0,
     &     B126= -6D0/41D0,    B136= -289D0/82D0, B87 = 13D0/9D2,
     &     B97 = 67D0/9D1,     B107= -19D0/6D1,   B117= 2133D0/41D2,
     &     B127= -3D0/205D0,   B137= 2193D0/41D2, B108= 17D0/6D0,
     &     B118= 45D0/82D0,    B128= -3D0/41D0,   B138= 51D0/82D0,
     &     B119= 45D0/164D0,   B139= 33D0/164D0,  B1110= 18D0/41D0,
     &     B1310= 12D0/41D0)
C     
      integer IR,n,i
      DIMENSION X(N), TOL(N)
      DIMENSION F0(42),F1(42),F2(42),F3(42),F4(42),F5(42),F6(42),F7(42)
C     

      IF (N .GT. 42) STOP 'N > 42'
C     
      IR = 0
      CALL DER(T, X, F1)
C     
 104  DO I = 1, N
         F0(I) = X(I) + DT*B21*F1(I)
      ENDDO
      CALL DER(T + AL2*DT, F0, F2)
      DO I = 1, N
         F0(I) = X(I) + DT*(B31*F1(I) + B32*F2(I))
      ENDDO
      CALL DER(T + AL3*DT, F0, F3)
      DO I = 1, N
         F0(I) = X(I) + DT*(B41*F1(I) + B43*F3(I))
      ENDDO
      CALL DER(T + AL4*DT, F0, F4)
      DO I = 1, N
         F0(I) = X(I) + DT*(B51*F1(I) + B53*(F3(I) - F4(I)))
      ENDDO
      CALL DER(T + AL5*DT, F0, F5)
      DO I = 1, N
         F0(I) = X(I) + DT*(B61*F1(I) + B64*F4(I) + B65*F5(I))
      ENDDO
      CALL DER(T + AL6*DT, F0, F6)
      DO I = 1, N
         F0(I) = X(I) + DT*(B71*F1(I) + B74*F4(I) + B75*F5(I) +
     &        B76*F6(I))
      ENDDO
      CALL DER(T + AL7*DT, F0, F7)
      DO I = 1, N
         F0(I) = X(I) + DT*(B81*F1(I) + B85*F5(I) + B86*F6(I) +
     &        B87*F7(I))
      ENDDO
      CALL DER(T + AL4*DT, F0, F2)
      DO I = 1, N
         F0(I) = X(I) + DT*(2D0*F1(I) + B94*F4(I) + B95*F5(I) +
     &        B96*F6(I) + B97*F7(I) + 3D0*F2(I))
      ENDDO
      CALL DER(T + AL9*DT, F0, F3)
      DO I = 1, N
         X1 = F1(I)
         X4 = F4(I)
         X5 = F5(I)
         X6 = F6(I)
         X7 = F7(I)
         X8 = F2(I)
         X9 = F3(I)
         F2(I) = CH1*X6 + CH2*(X7 + X8) + CH3*X9
         F0(I) = X(I) + DT*(B101*X1 + B104*X4 + B105*X5 + B106*X6 +
     &        B107*X7 + B108*X8 - B32*X9)
         F4(I) = B111*X1 + B114*X4 + B115*X5 + B116*X6 + B117*X7 +
     &        B118*X8 + B119*X9
         F5(I) = B121*X1 + B126*X6 + B127*X7 + B128*(X8 - X9)
         F6(I) = B131*X1 + B114*X4 + B115*X5 + B136*X6 + B137*X7 +
     &        B138*X8 + B139*X9
      ENDDO
      CALL DER(T + ALA*DT, F0, F3)
      DO I = 1, N
         F7(I) = X(I) + DT*(F4(I) + B1110*F3(I))
         F0(I) = X(I) + DT*(F5(I) - B126*F3(I))
      ENDDO
      CALL DER(T + DT, F7, F4)
      CALL DER(T,      F0, F5)
      DO I = 1, N
         F0(I) = X(I) + DT*(F6(I) + B1310*F3(I) + F5(I))
      ENDDO
      CALL DER(T + DT, F0, F6)
      X7 = 1D-30
      DO I = 1, N
         F0(I) = X(I)
         X(I) = X(I) + DT*(CH3*F3(I) + CH4*(F5(I) + F6(I)) + F2(I))
         F7(I) = DT*(F1(I) + F4(I) - F5(I) - F6(I))*CH4
         X7 = X7 + (F7(I)/TOL(I))**2
      ENDDO
      X9 = DT
      DT = DT*(25D-4/X7)**625D-4
      IF (X7 .GT. 1D0) THEN
         DO I = 1, N
            X(I) = F0(I)
         ENDDO
         IR = IR + 1
         GOTO 104
      ENDIF
      T = T + X9
      RETURN
      END



      SUBROUTINE FMMsplineDouble(x, y, b, c, d,num)
!     ---------------------------------------------------------------------------
!     PURPOSE - Compute the coefficients b,c,d for a cubic interpolating spline
!     so that the interpolated value is given by
!     s(x) = y(k) + b(k)*(x-x(k)) + c(k)*(x-x(k))**2 + d(k)*(x-x(k))**3
!     when x(k) <= x <= x(k+1)
!     The end conditions match the third derivatives of the interpolated curve to
!     the third derivatives of the unique polynomials thru the first four and
!     last four points.
!     Use Seval or Seval3 to evaluate the spline.

      
      REAL(8), INTENT(IN)  :: x(num) ! abscissas of knots
      REAL(8), INTENT(IN)  :: y(num) ! ordinates of knots
      REAL(8), INTENT(OUT) :: b(num) ! linear coeff
      REAL(8), INTENT(OUT) :: c(num) ! quadratic coeff.
      REAL(8), INTENT(OUT) :: d(num) ! cubic coeff.
      
      INTEGER:: k,n
      REAL(8):: t
      REAL(8),PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
!----------------------------------------------------------------------------
      n=SIZE(x)
c     print*,n
      IF (n < 3) THEN		! Straight line - special case for n < 3
         b(1)=ZERO
         IF (n == 2) b(1)=(y(2)-y(1))/(x(2)-x(1))
         c(1)=ZERO
         d(1)=ZERO
         IF (n < 2) RETURN
         b(2)=b(1)
         c(2)=ZERO
         d(2)=ZERO
         RETURN
      END IF

!.....Set up tridiagonal system.........................................
!     .    b=diagonal, d=offdiagonal, c=right-hand side
      d(1)=x(2)-x(1)
      c(2)=(y(2)-y(1))/d(1)
      DO k=2,n-1
         d(k)=x(k+1)-x(k)
         b(k)=TWO*(d(k-1)+d(k))
         c(k+1)=(y(k+1)-y(k))/d(k)
         c(k)=c(k+1)-c(k)
      END DO

!.....End conditions.  third derivatives at x(1) and x(n) obtained
!     .       from divided
!     differences.......................................
      b(1)=-d(1)
      b(n)=-d(n-1)
      c(1)=ZERO
      c(n)=ZERO
      IF (n > 3) THEN
         c(1)=c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
         c(n)=c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
         c(1)=c(1)*d(1)*d(1)/(x(4)-x(1))
         c(n)=-c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF

      DO k=2,n                  ! forward elimination
         t=d(k-1)/b(k-1)
         b(k)=b(k)-t*d(k-1)
         c(k)=c(k)-t*c(k-1)
      END DO

      c(n)=c(n)/b(n)		! back substitution ( makes c the sigma of text)
      DO k=n-1,1,-1
         c(k)=(c(k)-d(k)*c(k+1))/b(k)
      END DO

!.....Compute polynomial coefficients...................................
      b(n)=(y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
         b(k)=(y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
         d(k)=(c(k+1)-c(k))/d(k)
         c(k)=THREE*c(k)
      END DO
      c(n)=THREE*c(n)
      d(n)=d(n-1)

      RETURN
      END Subroutine FMMsplineDouble   

      
      FUNCTION SevalDouble(u, x,y, b,c,d,num) RESULT(SevalResult)
!     ---------------------------------------------------------------------------
!     PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!     where  x(i) <= u < x(i+1)

!     NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

      REAL(8),INTENT(IN) :: u	! abscissa at which the spline is to be evaluated
      REAL(8),INTENT(IN):: x(num) ! abscissas of knots
      REAL(8),INTENT(IN):: y(num) ! ordinates of knots
      REAL(8),INTENT(IN):: b(num),c(num),d(num) ! linear,quadratic,cubic coeff
      REAL(8):: SevalResult
      integer, intent(in)::num
      INTEGER, SAVE :: i=1
      INTEGER :: j, k, n
      REAL(8):: dx
!----------------------------------------------------------------------------
      n=SIZE(x)
c     print*,num,n

!.....First check if u is in the same interval found on the
!     last call to Seval.............................................
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
         i=1			! binary search
         j=n+1

         DO
            k=(i+j)/2
            IF (u < x(k)) THEN
               j=k
            ELSE
               i=k
            END IF
            IF (j <= i+1) EXIT
         END DO
      END IF

      dx=u-x(i)                 ! evaluate the spline
      SevalResult=(y(i)+dx*(b(i)+dx*(c(i)+dx*d(i))))
      RETURN
      END     

      
      SUBROUTINE Seval3Double(u, x,y, b,c,d, f,fp,fpp,fppp,num)
!     ---------------------------------------------------------------------------
!     PURPOSE - Evaluate the cubic spline function
!     Seval=y(i)+b(i)!(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
!     where  x(i) <= u < x(i+1)

!     NOTES- if u<x(1), i=1 is used;if u>x(n), i=n is used

      
      REAL(8),INTENT(IN) :: u	! abscissa at which the spline is to be evaluated
      REAL(8),INTENT(IN):: x(num) ! abscissas of knots
      REAL(8),INTENT(IN):: y(num) ! ordinates of knots
      REAL(8),INTENT(IN):: b(num),c(num),d(num) ! linear,quadratic,cubic coeff
      REAL(8),INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp ! function, 1st,2nd,3rd deriv
      integer, intent(in)::num
      INTEGER, SAVE :: i=1
      INTEGER :: j, k, n
      REAL(8)    :: dx
      REAL(8),PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
!----------------------------------------------------------------------------
      n=SIZE(x)

!.....First check if u is in the same interval found on the
!     last call to Seval.............................................
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
         i=1			! binary search
         j=n+1

         DO
            k=(i+j)/2
            IF (u < x(k)) THEN
               j=k
            ELSE
               i=k
            END IF
            IF (j <= i+1) EXIT
         END DO
      END IF

      dx=u-x(i)                 ! evaluate the spline
      IF (Present(f))f=y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))fp=b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))fpp=TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp=SIX*d(i)

      RETURN
      END Subroutine Seval3Double
      

      

      
      SUBROUTINE SORTRX(N,DATA,INDEX)
c     ordinamento veloce con quicksort  
c     piu' piccolo DATA(index(1))
c     piu' grande  DATA(index(N))
      INTEGER   N,INDEX(N)
      REAL      DATA(N)
      
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      REAL      DATAP
      
C     QuickSort Cutoff
C     

      INTEGER   M
      PARAMETER (M=9)
      
C===================================================================
C     
C     fai una stima iniziale dell'indice
      
      DO 50 I=1,N
         INDEX(I)=I
 50   CONTINUE
      
C     If array is short, skip QuickSort and go directly to
C     the straight insertion sort.
      
      IF (N.LE.M) GOTO 900
      
      ISTK=0
      L=1
      R=N
      
 200  CONTINUE
      
      
      I=L
      J=R
      
      
      P=(L+R)/2
      INDEXP=INDEX(P)
      DATAP=DATA(INDEXP)
      
      IF (DATA(INDEX(L)) .GT. DATAP) THEN
         INDEX(P)=INDEX(L)
         INDEX(L)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
      
      IF (DATAP .GT. DATA(INDEX(R))) THEN
         IF (DATA(INDEX(L)) .GT. DATA(INDEX(R))) THEN
            INDEX(P)=INDEX(L)
            INDEX(L)=INDEX(R)
         ELSE
            INDEX(P)=INDEX(R)
         ENDIF
         INDEX(R)=INDEXP
         INDEXP=INDEX(P)
         DATAP=DATA(INDEXP)
      ENDIF
      
      
 300  CONTINUE
      
      
      I=I+1
      IF (DATA(INDEX(I)).LT.DATAP) GOTO 300
      
 400  CONTINUE
      
      J=J-1
      IF (DATA(INDEX(J)).GT.DATAP) GOTO 400
      
C     Q5: collisione?
      
      IF (I.LT.J) THEN
         
C     Q6: interscambio DATA[I] <--> DATA[J] 
         
         INDEXT=INDEX(I)
         INDEX(I)=INDEX(J)
         INDEX(J)=INDEXT
         GOTO 300
      ELSE
         
C     Q7: Yes, select next subsequence to sort
C     
C     A questo punto, I >= J ae DATA[l] <= DATA[I] == DATAP <= DATA[r],
C     per tutti L <= l < I e J < r <= R.  
         
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
C     Q8: Pop the stack, or terminate QuickSort if empty
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
      
 900  CONTINUE
      
C===================================================================
C     
C     Q9: Straight Insertion sort
      
      DO 950 I=2,N
         IF (DATA(INDEX(I-1)) .GT. DATA(INDEX(I))) THEN
            INDEXP=INDEX(I)
            DATAP=DATA(INDEXP)
            P=I-1
 920        CONTINUE
            INDEX(P+1) = INDEX(P)
            P=P-1
            IF (P.GT.0) THEN
               IF (DATA(INDEX(P)).GT.DATAP) GOTO 920
            ENDIF
            INDEX(P+1) = INDEXP
         ENDIF
 950  CONTINUE
      
C===================================================================
C     
C     All done
      
      END
