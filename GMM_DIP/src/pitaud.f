c    subroutine pitaud.f
c    modified tipitaud that does not use COMMON for pi & tau
c    fnr is initialized in gmm01f.f
c    calculates pi(m,n) & tau(m,n) up to a specified nmax for all
c    m=0,1,...n at a given x
c    pi(m,n) & tau(m,n) calculated are normalized by
c           C_mn=[(2n+1)(n-m)!/n/(n+1)/(n+m)!]^(1/2)
c    Yu-lin Xu    12/2000

      subroutine pitaud(nmax, x, pi, tau)

      implicit logical (a-z)

      integer np, nLp
      include 'gmm01f.par'
      integer nmp0
      parameter (nmp0=(np+1)*(np+4)/2)

      integer nmax
      double precision x
      double precision pi(nmp0),tau(nmp0)

      double precision fnr(0:2*(np+2))
      common/fnr/fnr

      integer nt
      double precision sx, t, t125, tx
      integer imn, m, n, i, j, i1, i2
      
      nt=(nmax+1)*(nmax+4)/2          ! calculates pi up to nmax+1
      if(nt.gt.nmp0.or.dabs(x).gt.1.d0) then
         write(6,*) 'dimension or argument wrong in sub. tipitaud'
         write(6,*) 'argument: ',x
         stop
      endif
      sx=dsqrt(1.d0-x*x)
      pi(1)=0.d0                     ! pi(0,1)  pi(0,n)=0 when m=0
      pi(2)=dsqrt(.75d0)             ! pi(1,1)
      pi(3)=0.d0                     ! pi(0,2)
      t125=dsqrt(1.25d0)
      pi(4)=t125*x                      ! pi(1,2)
      pi(5)=t125*sx                     ! pi(2,2)
      imn=5
      do i=3,nmax+1
         n=i
         imn=imn+1
         pi(imn)=0.d0                ! pi(0,n)=0
         do 11 j=2,n
            m=j-1
            imn=imn+1
            i1=(n-2)*(n+1)/2+m+1
            if(m.eq.n-1) then
               pi(imn)=fnr(n-1)*fnr(2*n+1)/fnr(n+1)*x*pi(i1)
               goto 11
            endif
            i2=(n-3)*n/2+m+1
            t=fnr(n)*fnr(2*n-3)
            t=fnr(n-2)*fnr(n-m-1)*fnr(n+m-1)/t
            pi(imn)=fnr(2*n-1)*x*pi(i1)-t*pi(i2)
            t=fnr(n+1)*fnr(n-m)*fnr(n+m)
            t=fnr(n-1)*fnr(2*n+1)/t
            pi(imn)=t*pi(imn)
 11	     continue
         imn=imn+1
         i1=(n-2)*(n+1)/2+n
         t=fnr(n-1)*fnr(n+1)
         t=dsqrt(.5d0)*fnr(n)*fnr(2*n+1)/t
         pi(imn)=t*sx*pi(i1)
      enddo
      tx=x*sx
      tau(1)=-dsqrt(1.5d0)*sx          ! tau(0,1)
      tau(2)=pi(2)*x                   ! tau(1,1)
      tau(3)=-dsqrt(7.5d0)*tx          ! tau(0,2)
      tau(4)=t125*(2.d0*x*x-1.d0)      ! tau(1,2)
      tau(5)=t125*tx                   ! tau(2,2)
      imn=5
      do i=3,nmax
         n=i
         do 21 j=1,n+1
            m=j-1
            imn=imn+1
            if(m.eq.0) then
               i1=(n-2)*(n+1)/2+1
               i2=(n-3)*n/2+1
               t=fnr(2*n-3)
               t=fnr(n-2)*fnr(n)/t
               tau(imn)=fnr(2*n-1)*x*tau(i1)-t*tau(i2)
               t=fnr(n-1)*fnr(n+1)
               t=fnr(2*n+1)/t
               tau(imn)=t*tau(imn)
               goto 21
            endif
            i1=n*(n+3)/2+m+1
            t=fnr(n)*fnr(2*n+3)
            t=fnr(n+2)*fnr(2*n+1)*fnr(n-m+1)*fnr(n+m+1)/t
            tau(imn)=t*pi(i1)-dble(n+1)*x*pi(imn)
            tau(imn)=tau(imn)/dble(m)
 21	     continue
      enddo
      return
      end

