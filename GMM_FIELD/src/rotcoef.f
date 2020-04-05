c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c
c  This subroutine is originally written by D.W. Mackowski
c  (taken from Mackowski's multisphere-scattering code scsmfo1b.for
c  released to public by the author at 8/1999).
c  Very slightly modified for fitting into this code.
c  Yu-lin Xu   12/2000
c
      subroutine rotcoef(cbe,nmax)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      implicit double precision (a-h,o-z)
      double precision dk0(-2*np:2*np),dk01(-2*np:2*np)
      common/rot/bcof(0:np+2),dc(-np:np,0:nmp)
      common/fnr/fnr(0:2*(np+2))
      sbe=dsqrt((1.d0+cbe)*(1.d0-cbe))
      cbe2=.5d0*(1.d0+cbe)
      sbe2=.5d0*(1.d0-cbe)
      in=1
      dk0(0)=1.d0
      sben=1.d0
      dc(0,0)=1.d0
      dk01(0)=0.d0
      do n=1,nmax
         nn1=n*(n+1)
         in=-in
         sben=sben*sbe/2.d0
         dk0(n)=dble(in)*sben*bcof(n)
         dk0(-n)=dble(in)*dk0(n)
         dk01(n)=0.d0
         dk01(-n)=0.d0
         dc(0,nn1+n)=dk0(n)
         dc(0,nn1-n)=dk0(-n)
         do k=-n+1,n-1
            kn=nn1+k
            dkt=dk01(k)
            dk01(k)=dk0(k)
            dk0(k)=(cbe*dble(n+n-1)*dk01(k)-fnr(n-k-1)*fnr(n+k-1)*dkt)
     1             /(fnr(n+k)*fnr(n-k))
            dc(0,kn)=dk0(k)
         enddo
         im=1
         do m=1,n
            im=-im
            fmn=1.d0/fnr(n-m+1)/fnr(n+m)
            m1=m-1
            dkm0=0.d0
            do k=-n,n
               kn=nn1+k
               dkm1=dkm0
               dkm0=dc(m1,kn)
               if(k.eq.n) then
                  dkn1=0.d0
               else
                  dkn1=dc(m1,kn+1)
               endif
               dc(m,kn)=(fnr(n+k)*fnr(n-k+1)*cbe2*dkm1
     1           -fnr(n-k)*fnr(n+k+1)*sbe2*dkn1
     1              -dble(k)*sbe*dc(m1,kn))*fmn
               dc(-m,nn1-k)=dble((-1)**(k)*im)*dc(m,kn)
            enddo
         enddo
      enddo
      return
      end
      
