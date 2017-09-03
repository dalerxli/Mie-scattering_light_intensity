! *********************************************************************************************************
! * Calculating the scattered light by a spherical particle using Bohren-Huffman Mie scattering algorithm *
! * Measuring the intensity of collected light by a mirror/lens/sensor with an arbitrary size             *
! * Developed by A. F. Forughi (Aug. 2012), Mech. Eng. Dept., Sharif University of Technology             *
! *********************************************************************************************************

implicit none
INTEGER MXNANG; PARAMETER(MXNANG=1000)
INTEGER IREADEP,J,NAN,NANG,NANG0,k
REAL AJ,ANG,DANG,GSCA,PI,POL,QABS,QBACK,QEXT,QSCA,DIA,REFMED,S11,S12,S33,S34,WAVEL,X,ih(2*MXNANG-1),iv(2*MXNANG-1),coeffd
COMPLEX REFREL,CXEPS,S1(2*MXNANG-1),S2(2*MXNANG-1)
REAL tti,tmin,tmax,Sumih,Sumiv,r,dmin,dmax,coeff
INTEGER jmin,jmax,nmirr,nd
PI=4.E0*atan(1.E0)
open(UNIT=7,FILE='Result.plt')

6012 format(' Complex refractive index=',1P2E10.3)
6013 format(' diameter=',1PE11.4,' lambda=',E11.4,' x=',E11.4)
6017 format(2X,'theta',7X,'S11',11X,'POL',11X,'S33',11X,'S34')
6065 format(/,'Qext=',1PE11.4,' Qsca=',E11.4,' Qabs=',E11.4,' <cos>=',E11.4,/,17X,'Qbk =',E11.4)
6075 format(1X,F6.2,2X,1PE12.5,2X,E12.5,2X,E12.5,2X,E12.5)

REFMED=1.0  !Real refractive index of the surrounding medium'
IREADEP=0   !Wish to enter refr.index or epsilon? (0 or 1)'
if(IREADEP.LE.0)then !complex refractive index of sphere, Imaganary part must be Positive:
	!REFREL=(1.75,0.43)       !SOOT
	!REFREL=(1.53,0.008)      !DUST
	!REFREL=(1.332,1.46e-8)   !WATER
	!REFREL=(1.49,2.0e-8)     !Sea Salt
	!REFREL=(1.59,0.0)        !PSL (Duck Sci.)
	REFREL=(1.2,0.0)
else
    read(*,*)CXEPS !complex epsilon of sphere'
    REFREL=sqrt(CXEPS)
endif

DIA=1.0e-6     !particle diameter
WAVEL=650.0e-9 !wavelength'
NANG0=500      !number of angles between 0 and 90(incl. 0 and 90)

REFREL=REFREL/REFMED !Relative n
write(*,6012)REFREL


if(NANG0.GT.MXNANG)STOP'*Error: NANG > MXNANG'
NANG=NANG0
if(NANG0.LT.2)NANG=2
X=PI*DIA*REFMED/WAVEL !Shape Factor!
write(*,6013)DIA,WAVEL,X

if(NANG.GT.1)DANG=0.5E0*PI/float(NANG-1)
r=1.0d0 !Mirror distance

tmin=5.*PI/180.   !Teta Minimum
tmax=38.*PI/180.   !Teta Maximum

NAN=2*NANG-1
jmin=aint(real(NAN)*tmin/PI+1.0)+1
jmax=aint(real(NAN)*tmax/PI+1.0)
nmirr=jmax-jmin
print*,"No. of elements on the mirror= ",nmirr

dmin=0.1e-6
dmax=10.0e-6
nd=3000      !No. of ds between dmin and dmax

Do k=0,nd !Diameters
	X=PI*(dmin+real(k)*(dmax-dmin)/real(nd))*REFMED/WAVEL !Shape Factor!
	call BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
    !QABS=QEXT-QSCA
    !WRITE(*,6065)QEXT,QSCA,QABS,GSCA,QBACK
    if(NANG0.GT.1)then
        NAN=2*NANG-1
		do J=1,NAN
			iv(j)=CABS(S1(J))**2/((2*PI/WAVEL)*(2*PI/WAVEL))
			ih(j)=CABS(S2(J))**2/((2*PI/WAVEL)*(2*PI/WAVEL))
			!ANG=DANG*(J-1.E0)*180.E0/PI
			!WRITE(7,6075)ANG,iv,ih,(iv+ih)/2.0
        enddo
    endif

	Sumiv=0.0d0;Sumih=0.0d0;coeffd=0.0d0
	do j=jmin,jmax;
		tti=real(j-1)*PI/real(NAN)
		coeff=2.0*r*sqrt((sin((tmax-tmin)/2.0))**2.0-(cos((tmax-tmin)/2.0)*tan((tmax+tmin)/2.0-tti))**2.0)
		coeff=coeff*r*(1.0+(tan((tmax+tmin)/2.0-tti))**2.0)*((cos((tmax+tmin)/2.0-tti))**2.0)/cos((tmax-tmin)/2.0)/(r*r)
		coeff=coeff*(tmax-tmin)/real(nmirr)
		Sumiv=Sumiv+iv(j)*coeff
		Sumih=Sumih+ih(j)*coeff
		coeffd=coeffd+coeff
	enddo
	write(7,*) (dmin+real(k)*(dmax-dmin)/real(nd)),(Sumiv+Sumih)/2.0
enddo

close(7)

STOP
END



SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
! This Subroutine Adapted by B.T.Draine, Princeton Univ. Obs.
! Declare parameters:
! Note: important that MXNANG be consistent with dimension of S1 and S2 in calling routine!
      INTEGER MXNANG,NMXX
! PARAMETER(MXNANG=1000,NMXX=15000)
      PARAMETER(MXNANG=1000,NMXX=150000)
! Arguments:
      INTEGER NANG
      REAL GSCA,QBACK,QEXT,QSCA,X
      COMPLEX REFREL
      COMPLEX S1(2*MXNANG-1),S2(2*MXNANG-1)
! Local variables:
      INTEGER J,JJ,N,NSTOP,NMX,NN
      DOUBLE PRECISION CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,THETA,XSTOP,YMOD
      DOUBLE PRECISION AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG),TAU(MXNANG)
      DOUBLE COMPLEX AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
      DOUBLE COMPLEX D(NMXX)
!***********************************************************************
! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
!    to calculate scattering and absorption by a homogenous isotropic
!    sphere.
! Given:
!    X = 2*pi*a/lambda
!    REFREL = (complex refr. index of sphere)/(real index of medium)
!    NANG = number of angles between 0 and 90 degrees
!           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
!           if called with NANG<2, will set NANG=2 and will compute
!           scattering for theta=0,90,180.
! Returns:
!    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
!                                scatt. E perp. to scatt. plane)
!    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
!                                scatt. E parr. to scatt. plane)
!    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
!    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
!    QBACK = (dC_sca/domega)/pi*a**2
!          = backscattering efficiency [NB: this is (1/4*pi) smaller
!            than the "radar backscattering efficiency"; see Bohren &
!            Huffman 1983 pp. 120-123]
!    GSCA = <cos(theta)> for scattering
!
! Original program taken from Bohren and Huffman (1983), Appendix A
! Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
! in order to compute <cos(theta)>
! 91/05/07 (BTD): Modified to allow NANG=1
! 91/08/15 (BTD): Corrected error (failure to initialize P)
! 91/08/15 (BTD): Modified to enhance vectorizability.
! 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
! 91/08/15 (BTD): Changed definition of QBACK.
! 92/01/08 (BTD): Converted to full double precision and double complex
!                 eliminated 2 unneed lines of code
!                 eliminated redundant variables (e.g. APSI,APSI0)
!                 renamed RN -> EN = double precision N
!                 Note that DOUBLE COMPLEX and DCMPLX are not part
!                 of f77 standard, so this version may not be fully
!                 portable.  In event that portable version is
!                 needed, use src/bhmie_f77.f
! 93/06/01 (BTD): Changed AMAX1 to generic function MAX
!***********************************************************************
!*** Safety checks
      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2
!*** Obtain pi:
      PII=4.*ATAN(1.D0)
      DX=X
      DREFRL=REFREL
      Y=X*DREFRL
      YMOD=ABS(Y)
!
!*** Series expansion terminated after NSTOP terms Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.*X**0.3333+2.
      NMX=MAX(XSTOP,YMOD)+15
! BTD experiment 91/1/15: add one more term to series and compare results
! NMX=AMAX1(XSTOP,YMOD)+16
! test: compute 7001 wavelengths between .0001 and 1000 micron
! for a=1.0micron SiC grain.  When NMX increased by 1, only a single
! computed number changed (out of 4*7001) and it only changed by 1/8387
! conclusion: we are indeed retaining enough terms in series!
      NSTOP=XSTOP

      IF(NMX.GT.NMXX)THEN
          WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
          STOP
      ENDIF
!*** Require NANG.GE.1 in order to calculate scattering intensities
      DANG=0.
      IF(NANG.GT.1)DANG=.5*PII/DBLE(NANG-1)
      DO J=1,NANG
          THETA=DBLE(J-1)*DANG
          AMU(J)=COS(THETA)
      ENDDO
      DO J=1,NANG
          PI0(J)=0.
          PI1(J)=1.
      ENDDO
      NN=2*NANG-1
      DO J=1,NN
          S1(J)=(0.,0.)
          S2(J)=(0.,0.)
      ENDDO
!
!*** Logarithmic derivative D(J) calculated by downward recurrence
!    beginning with initial value (0.,0.) at J=NMX
!
      D(NMX)=(0.,0.)
      NN=NMX-1
      DO  N=1,NN
          EN=NMX-N+1
          D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
      ENDDO
!
!*** Riccati-Bessel functions with real argument X
!    calculated by upward recurrence
!
      PSI0=COS(DX)
      PSI1=SIN(DX)
      CHI0=-SIN(DX)
      CHI1=COS(DX)
      XI1=DCMPLX(PSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO N=1,NSTOP
          EN=N
          FN=(2.E0*EN+1.)/(EN*(EN+1.))
! for given N, PSI  = psi_n        CHI  = chi_n
!              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
!              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
! Calculate psi_n and chi_n
          PSI=(2.E0*EN-1.)*PSI1/DX-PSI0
          CHI=(2.E0*EN-1.)*CHI1/DX-CHI0
          XI=DCMPLX(PSI,-CHI)
!
!*** Store previous values of AN and BN for use
!    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
              AN1=AN
              BN1=BN
          ENDIF
!
!*** Compute AN and BN:
          AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
          AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
          BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
          BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)
!
!*** Augment sums for Qsca and g=<cos(theta)>
          QSCA=QSCA+(2.*EN+1.)*(ABS(AN)**2+ABS(BN)**2)
          GSCA=GSCA+((2.*EN+1.)/(EN*(EN+1.)))*(REAL(AN)*REAL(BN)+IMAG(AN)*IMAG(BN))
          IF(N.GT.1)THEN
              GSCA=GSCA+((EN-1.)*(EN+1.)/EN)*(REAL(AN1)*REAL(AN)+IMAG(AN1)*IMAG(AN)+REAL(BN1)*REAL(BN)+IMAG(BN1)*IMAG(BN))
          ENDIF
!
!*** Now calculate scattering intensity pattern
!    First do angles from 0 to 90
          DO J=1,NANG
              JJ=2*NANG-J
              PI(J)=PI1(J)
              TAU(J)=EN*AMU(J)*PI(J)-(EN+1.)*PI0(J)
              S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
              S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
          ENDDO
!
!*** Now do angles greater than 90 using PI and TAU from
!    angles less than 90.
!    P=1 for N=1,3,...; P=-1 for N=2,4,...
          P=-P
          DO J=1,NANG-1
              JJ=2*NANG-J
              S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
              S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
          ENDDO
          PSI0=PSI1
          PSI1=PSI
          CHI0=CHI1
          CHI1=CHI
          XI1=DCMPLX(PSI1,-CHI1)
!
!*** Compute pi_n for next value of n
!    For each angle J, compute pi_n+1
!    from PI = pi_n , PI0 = pi_n-1
          DO J=1,NANG
              PI1(J)=((2.*EN+1.)*AMU(J)*PI(J)-(EN+1.)*PI0(J))/EN
              PI0(J)=PI(J)
         ENDDO
      ENDDO
!
!*** Have summed sufficient terms.
!    Now compute QSCA,QEXT,QBACK,and GSCA
      GSCA=2.*GSCA/QSCA
      QSCA=(2./(DX*DX))*QSCA
      QEXT=(4./(DX*DX))*REAL(S1(1))
      QBACK=(ABS(S1(2*NANG-1))/DX)**2/PII
      RETURN
END
