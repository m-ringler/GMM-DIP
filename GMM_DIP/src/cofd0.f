c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      subroutine cofd0(nmax)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
      integer v
      double precision lnfacd
      common/cofmnv0/cof0(ni0)
      common/crot/cofsr(nmp)
          common/fnr/fnr(0:2*(np+2))
      i=0
      sm=-0.5d0*dble((-1)**nmax)
      do m=-nmax,nmax
         ns=max(1,iabs(m))
         sm=-sm
         do n=ns,nmax
            inm=n*(n+1)-m
            do v=ns,nmax
               i=i+1
               ivm=v*(v+1)+m
               c=cofsr(inm)+cofsr(ivm)
               c=sm*dexp(c)
               c0=fnr(2*n+1)*fnr(2*v+1)
               c1=fnr(n)*fnr(v)*fnr(n+1)*fnr(v+1)
               c0=c0/c1
                   cof0(i)=c*c0
                enddo
             enddo
          enddo
      return
      end
 
