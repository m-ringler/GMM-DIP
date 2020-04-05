c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      subroutine trv(anpt,nodrj,nodri)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      complex*16 anpt(2,*),ant(2,2*np),a,b,atr(2,np,nmp)
      common/tran/atr
      mmax=min(nodrj,nodri)
      do m=-mmax,mmax
         n1=max(1,iabs(m))
         do n=n1,nodrj
            imn=n*(n+1)+m
            do ip=1,2
               ant(ip,n)=anpt(ip,imn)
            enddo
         enddo
         do n=n1,nodri
            imn=n*(n+1)+m
            a=0.d0
            b=0.d0
            do l=n1,nodrj
               ml=l*(l+1)+m
               a=a+atr(1,n,ml)*ant(1,l)
     1            +atr(2,n,ml)*ant(2,l)
               b=b+atr(1,n,ml)*ant(2,l)
     1            +atr(2,n,ml)*ant(1,l)
            enddo
            anpt(1,imn) = a
            anpt(2,imn) = b
         enddo
      enddo
      return
      end
