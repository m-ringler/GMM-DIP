c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c
c  This subroutine is the implementation of Mackowski's three-step
c  method for decomposition of translation matrix into rotational
c  and axial translational parts.
c  The rotation-translation-rotation formulation is described
c  in Mackowski, Proc. R. Soc. Lond. A 433, pp.611-612 (1991).
c  Yu-lin Xu   12/2000
c
      subroutine rtr(anpt,nodrj,nodri,ekt,drot)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter(nmp=np*(np+2))
      double precision drot(*)
      complex*16 anpt(2,*),ant(2,2*np),amt(2,-np:np),a,b,
     +           ekt(*),ek(-np:np),atr(2,np,nmp)
      common/tran/atr
c
      ek(0)=1.d0
      nmax=max(nodrj,nodri)
      do m=1,nmax
         ek(m)=ekt(m)
         ek(-m)=dconjg(ek(m))
      enddo
      irc=0
      do n=1,nodrj
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         do k=-n,n
            kn=n1+k
            a=ek(k)*anpt(1,kn)
            b=ek(k)*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         do m=-n,n
            imn=n1+m
            anpt(1,imn)=amt(1,m)
            anpt(2,imn)=amt(2,m)
         enddo
      enddo
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
      in=1
      irc=0
      do n=1,nodri
         in=-in
         n1=n*(n+1)
         do m=-n,n
            amt(1,m)=0.d0
            amt(2,m)=0.d0
         enddo
         sik=-in
         do k=-n,n
        sik=-sik
        kn=n1+k
        a=sik*anpt(1,kn)
            b=sik*anpt(2,kn)
            do m=-n,n
               irc=irc+1
               amt(1,m)=amt(1,m)+a*drot(irc)
               amt(2,m)=amt(2,m)+b*drot(irc)
            enddo
         enddo
         sik=-in
         do m=-n,n
            sik=-sik
            imn=n1+m
            anpt(1,imn)=amt(1,m)*ek(-m)*sik
            anpt(2,imn)=amt(2,m)*ek(-m)*sik
         enddo
      enddo
      return
      end

