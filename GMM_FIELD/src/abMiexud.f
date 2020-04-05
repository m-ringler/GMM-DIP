c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      SUBROUTINE abMiexud(X,REFREL,NP,NMAX,NM,AN,BN,NADD,
     +                    RSR,RSI,RSX,PX,AR,AI,BR,BI,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NMAX,NM,NADD,NSTOP,NX,K,N
      COMPLEX*16 REFREL,AN(NP),BN(NP)
      DOUBLE PRECISION AR(NP),AI(NP),BR(NP),BI(NP),
     +   RSR(NMAX),RSI(NMAX),RSX(NMAX),PX(NMAX)
      IF(EPS.GT.1.D0.OR.EPS.LT.0.D0) EPS=1.0D-20
      IF(NADD.NE.0) EPS=0.D0
      CTC=EPS
      XM=DBLE(REFREL)
      YM=DIMAG(REFREL)
      XMX=X*XM
      YMX=X*YM
      RP2=XMX*XMX+YMX*YMX
      NSTOP=X+4.D0*X**.3333D0
      NSTOP=NSTOP+2
      NM=NSTOP+NADD
      XN=DSQRT(XM**2+YM**2)*X
      NX=1.1D0*XN+10.D0
      if(NX-NM.lt.10) NX=NM+10
      write(6,*) 'Wiscombe criterion:',NSTOP
      write(6,*) 'NADD:',NADD
      write(6,*) 'NX:',NX
      IF(NX.GT.NMAX) THEN
         WRITE(6,*) 'Parameter NXMAX too small'
         WRITE(6,*) '  NXMAX must be greater than', NX
         WRITE(6,*) 'Please correct NXMAX in main code,'
         WRITE(6,*) '  recompile, then try again'
         STOP
      ENDIF
      IF(NM.GT.NP) THEN
         WRITE(6,*) 'Parameter np too small'
         WRITE(6,*) '  np must be greater than', NM
         WRITE(6,*) 'Please correct np in gmm01f.par,'
         WRITE(6,*) '  recompile the code, then try again'
         STOP
      ENDIF
C  DOWN RECURSION FOR RATIOS RSR,RSI,RSX,PNR,PNI,PX
      PNX=X/DBLE(2*NX+3)
      PNR=XMX/DBLE(2*NX+3)
      PNI=YMX/DBLE(2*NX+3)
      DO 5 K=1,NX
         N=NX-K+1
         CN=DBLE(N)
         ALN=(2.D0*CN+1.D0)*XMX/RP2-PNR
         BEN=(2.D0*CN+1.D0)*YMX/RP2+PNI
         RSR(N)=-CN*XMX/RP2+ALN
         RSI(N)=CN*YMX/RP2-BEN
         PZD=ALN*ALN+BEN*BEN
         PNR=ALN/PZD
         PNI=BEN/PZD
         RSX(N)=(CN+1.D0)/X-PNX
         IF(N.EQ.1) GO TO 20
         PNX=X/(2.D0*CN+1.D0-PNX*X)
         PX(N)=PNX
    5    CONTINUE
   20 SNM1X=DSIN(X)
      CNM1X=DCOS(X)
      IF(X-0.1D0)21,22,22
   21 SNX=X**2./3.D0-X**4./30.D0+X**6./840.D0-X**8./45360.D0
      GO TO 23
   22 SNX=SNM1X/X-CNM1X
   23 CNX=CNM1X/X+SNM1X
      DO 10 N=1,NX
         PX(N)=SNX      !preparing for the calculation of Cabs
         C=DBLE(N)
         DCNX=CNM1X-C*CNX/X
         DSNX=RSX(N)*SNX
C  CALCULATION OF EXTERIOR COEFFICIENTS AN AND BN
         ANNR=RSR(N)*SNX-XM*DSNX
         ANNI=RSI(N)*SNX-YM*DSNX
         TA1=RSR(N)*SNX-RSI(N)*CNX
         TA2=RSI(N)*SNX+RSR(N)*CNX
         ANDR=TA1-XM*DSNX+YM*DCNX
         ANDI=TA2-XM*DCNX-YM*DSNX
         AND=ANDR*ANDR+ANDI*ANDI
         BNNR=(XM*RSR(N)-YM*RSI(N))*SNX-DSNX
         BNNI=(XM*RSI(N)+YM*RSR(N))*SNX
         TB1=RSR(N)*SNX-RSI(N)*CNX
         TB2=RSR(N)*CNX+RSI(N)*SNX
         BNDR=XM*TB1-YM*TB2-DSNX
         BNDI=XM*TB2+YM*TB1-DCNX
         BND=BNDR*BNDR+BNDI*BNDI
         AR(N)=(ANNR*ANDR+ANNI*ANDI)/AND
         AI(N)=(ANNI*ANDR-ANNR*ANDI)/AND
         BR(N)=(BNNR*BNDR+BNNI*BNDI)/BND
         BI(N)=(BNNI*BNDR-BNNR*BNDI)/BND
C  MIE SERIES CONVERGENCE TEST IS MADE BY TESTING AN'S AND BN'S
         TI=AR(N)*AR(N)+AI(N)*AI(N)+BR(N)*BR(N)+BI(N)*BI(N)
         TI=TI/(AR(1)*AR(1)+AI(1)*AI(1)+BR(1)*BR(1)+BI(1)*BI(1))
         IF(TI-CTC) 16,18,18
   18    IF(NM-N) 15,15,6
    6    IF(N-NX)7,15,15
    7    M=N+1
         SNX=PX(M)*SNX
         CNM2X=CNM1X
         CNM1X=CNX
         CNX=(2.D0*C+1.D0)*CNM1X/X-CNM2X
   10 CONTINUE
      GO TO 15
   16 WRITE(6,*) '*** NOTE THAT THE FIELD-EXPANSION TRANCATION'
      WRITE(6,*) '*** IS DETERMINED BY eps GIVEN IN THE INPUT'
      WRITE(6,*) '*** FILE gmm01f.in'
      WRITE(6,*) '*** IN CASE YOU NEED A HIGHER ORDER, eps MUST'
      WRITE(6,'(a,e9.1)')
     +          ' *** BE SMALLER THAN THE CURRENT VALUE',EPS
   15 NM=N
      DO I=1,NM
         AN(I)=DCMPLX(AR(I),-AI(I))
         BN(I)=DCMPLX(BR(I),-BI(I))
      ENDDO
      RETURN
      END


