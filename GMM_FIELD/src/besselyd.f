c gmm01f
c  sub. besselyd.f  (in double precision arithmetic)
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  returns an array of the spherical Bessel function of
c  the second kind with a real argument: y_0,y_1,...,y_n
      subroutine besselyd(n,x,besy)
      implicit logical (a-z)
      
C SUBROUTINE PARAMETERS
      integer n
      double precision x,besy(0:n)
      
C LOCAL VARIABLES
      integer i
      double precision besyn,x2
      
      if(x .eq. 0.D0) then
        write(6,*) 'bad argument in sub. besselyd'
        stop
      endif
      
      if(dabs(x)-0.1d0)10,11,11
  10  x2=x*x
      besyn=1.d0-x2/72.d0
      besyn=1.d0-x2/42.d0*besyn
      besyn=1.d0-x2/20.d0*besyn
      besyn=1.d0-x2/6.d0*besyn
      besy(0)=1.d0-x2/56.d0
      besy(0)=1.d0-x2/30.d0*besy(0)
      besy(0)=1.d0-x2/12.d0*besy(0)
      besy(0)=1.d0-0.5d0*x2*besy(0)
      besy(0)=-besy(0)/x
      goto 12
  11  besyn=dsin(x)/x
      besy(0)=-dcos(x)/x
  12  besy(1)=besy(0)/x-besyn
      do i=2,n
        besy(i)=dble(2*i-1)/x*besy(i-1)-besy(i-2)
      enddo
      return
      end
