c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      SUBROUTINE orientcd(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
     &                  MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI,
     &                  BETA,THETA,PHI)
C Arguments:
      INTEGER MXBETA,MXTHET,MXPHI,NBETA,NTHETA,NPHI
      double precision BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX
      double precision BETA(MXBETA),THETA(MXTHET),PHI(MXPHI)
C Local variables:
      INTEGER J
      double precision DELTA
C***********************************************************************
C Given: BETAMI=minimum value of beta (radians)
C        BETAMX=maximum value of beta (radians)
C        THETMI=minimum value of theta (radians)
C        THETMX=maximum value of theta (radians)
C        PHIMIN=minimum value of phi (radians)
C        PHIMAX=maximum value of phi (radians)
C        MXBETA,MXTHET,MXPHI=dimensions of the arrays BETA,THETA,PHI
C        NBETA=number of values of beta
C        NTHETA=number of values of theta
C        NPHI=number of values of PHI
C Returns: BETA(1-NBETA)=beta values (radians)
C          THETA(1-NTHETA)=theta values (radians)
C          PHI(1-NPHI)=phi values (radians)
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        cos(theta) to be uniformly distributed between cos(THETMI) and
C                   cos(THETMX)
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C        If NPHI=1, first angle is THETMI, last angle is THETMX
C        If NPHI>1, then angles are midpoints of intervals of equal
C            range in theta subdividing range from THETMI to THETMX
C***********************************************************************
      BETA(1)=BETAMI
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      IF(NPHI.EQ.1.AND.NTHETA.GT.1)THEN
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA-1)
         THETA(1)=THETMI
      ELSE
         DELTA=(DCOS(THETMX)-DCOS(THETMI))/DBLE(NTHETA)
         THETA(1)=DACOS(DCOS(THETMI)+0.5d0*DELTA)
      ENDIF
      IF(NTHETA.GT.1)THEN
         DO 2000 J=2,NTHETA
            THETA(J)=DACOS(DCOS(THETA(1))+DELTA*DBLE(J-1))
 2000    CONTINUE
      ENDIF
c      DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI)
c      PHI(1)=PHIMIN+0.5D0*DELTA
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END


