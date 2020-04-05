c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      subroutine cofsrd(nmax)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      double precision cofsr(nmp),lnfacd,c
      common/crot/cofsr
          i=0
      do n=1,nmax
         do m=-n,n
            i=i+1
            c=lnfacd(dble(n-m))-lnfacd(dble(n+m))
            cofsr(i)=0.5d0*c
c	        c=0.5d0*c
c	        cofsr(i)=dexp(c)
             enddo
          enddo
      return
      end  

