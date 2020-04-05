c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
C  subroutine gau0.f generates tabulated values for
C  Gaunt coefficients up to n=v=n_max
      subroutine gau0(nmax)
      include 'gmm01f.par'
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
      parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
      integer v,qmax,uvmax,iga0(ni0)
      double precision ga0(ng0)
      common/g0/ga0
      common/ig0/iga0
          na=0
          uvmax=nmax*(nmax+2)
          i=0
      do m=-nmax,nmax
         ns=max(1,iabs(m))
         do n=ns,nmax
            do v=ns,nmax
               call gxurcd0(-m,n,v,qmax,na)
               i=i+1
                   iga0(i)=na
               na=na+qmax+1
                enddo
         enddo
          enddo
      return
      end

