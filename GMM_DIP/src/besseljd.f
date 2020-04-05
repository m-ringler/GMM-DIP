c gmm01f
c  subroutine besseljd.f  (in double precision arithmetic)
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  returns an array of the spherical Bessel function of the
c  first kind with a real argument: j_0,j_1,j_2,...,j_{NC}
c  uses Ru Wang's ratio method for the downward recursive
c  calculation of the Riccati-Bessel function Psi_n(z)=z j_n(z)
c  [see Xu et al., Physical Review E, v.60, 2347-2365 (1999)]
      SUBROUTINE besseljd(NC,X,BESJ)
      INTEGER NC,NX,K,N
      DOUBLE PRECISION X,BESJ(0:NC),PN,CN,X2
      DO K=1,NC
         BESJ(K)=0.D0
      ENDDO
      IF(DABS(X).LT.1.D-12) THEN
         BESJ(0)=1.D0
         BESJ(1)=1.D0/3.D0*X
         RETURN
      ENDIF
c  down-recursively calculating an array of the ratio functions
c  P_{NC},P_{NC-1},...,P(2) stored in the same array for j_n,
c  starting with an asymptotic value P_{NX+1}=X/(2 NX+3) at the
c  highest order NX+1, where NX=NC+1.1X+10
      NX=1.1D0*X+10.D0
      NX=NC+NX
      PN=X/DBLE(2*NX+3)
      DO 5 K=1,NX-1
         N=NX-K+1
         CN=DBLE(N)
         PN=X/(DBLE(2*N+1)-PN*X)
         IF(N.GT.NC) GOTO 5
         BESJ(N)=PN
    5 CONTINUE
C  calculating j_0(x) and j_1(x)
      IF(DABS(X)-0.1D0) 10,11,11
   10 X2=X*X
      BESJ(0)=1.D0-X2/72.D0
      BESJ(0)=1.D0-X2/42.D0*BESJ(0)
      BESJ(0)=1.D0-X2/20.D0*BESJ(0)
      BESJ(0)=1.D0-X2/6.D0*BESJ(0)
      BESJ(1)=1.D0/45360.D0-X2/3991680.D0
      BESJ(1)=1.D0/840.D0-X2*BESJ(1)
      BESJ(1)=1.D0/30.0d0-X2*BESJ(1)
      BESJ(1)=X*(1.D0/3.0d0-X2*BESJ(1))
      GOTO 12
   11 BESJ(0)=DSIN(X)/X
      BESJ(1)=(DSIN(X)/X-DCOS(X))/X
c  calculating j_2,...,j_{NC}
   12 DO 20 N=2,NC
         BESJ(N)=BESJ(N)*BESJ(N-1)
   20 CONTINUE
      RETURN
      END
      
