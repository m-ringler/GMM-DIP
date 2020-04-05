c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      subroutine cofnv0(nmax)
      include 'gmm01f.par'
      integer n,v
      double precision c1,lnfacd,cnv(np,np)
      common/cfnv/cnv
      do n=1,nmax
         do v=n,nmax
            c1=lnfacd(dble(2*n))+lnfacd(dble(2*v))
            c1=c1-lnfacd(dble(2*n+2*v))
            c1=c1+2.d0*lnfacd(dble(n+v))
            c1=c1-lnfacd(dble(n))-lnfacd(dble(v))
            cnv(n,v)=c1
         enddo
      enddo
      return
      end
