c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      SUBROUTINE orientud(BETAMI,BETAMX,THETMI,THETMX,PHIMIN,PHIMAX,
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
C Returns:  BETA(1-NBETA)=beta values (radians)
C           THETA(1-NTHETA)=theta values (radians)
C           PHI(1-NPHI)=phi values (radians)
C Note: it is assumed that target orientation weight function
C       can be factored into WGTA*WGTB -- i.e., that rotations
C       around a1 are decoupled from orientation of a1.
C Purpose: to generate a sequence of desired target orientations
C Present version assumes:
C        beta to be uniformly distributed between BETAMI and BETAMX
C        theta to be uniformly distributed between THETMI and THETMX
C        phi to be uniformly distributed between PHIMIN and PHIMAX
C            first angle is THETMI, last angle is THETMX
C***********************************************************************
      BETA(1)=BETAMI
      IF(NBETA.GT.1)THEN
         DELTA=(BETAMX-BETAMI)/DBLE(NBETA-1)
         DO 1000 J=2,NBETA
           BETA(J)=BETA(1)+DELTA*DBLE(J-1)
 1000    CONTINUE
      ENDIF
      THETA(1)=THETMI
      IF(NTHETA.GT.1)THEN
         DELTA=(THETMX-THETMI)/DBLE(NTHETA-1)
         DO 2000 J=2,NTHETA
            THETA(J)=THETA(1)+DELTA*DBLE(J-1)
 2000    CONTINUE
      ENDIF
      PHI(1)=PHIMIN
      IF(NPHI.GT.1)THEN
         DELTA=(PHIMAX-PHIMIN)/DBLE(NPHI-1)
         DO 3000 J=2,NPHI
            PHI(J)=PHI(1)+DELTA*DBLE(J-1)
 3000    CONTINUE
      ENDIF
      RETURN
      END



