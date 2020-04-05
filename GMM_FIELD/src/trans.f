c gmm01f
c For each sphere j, calculates the expansion
c a1(m n j) and b1(m n j)
c of the field originating from the other spheres l != j
c in the basis of VSWF M1(m n j) and N1(m n j)
c in the j-centered coordinate system.
c The input fields are given as expansions
c a(mu nu l) and b(mu nu l)
c in the basis of VSWF M3(mu nu l) and N3(mu nu l)
c in the l-centered coordinate systems.
c The formulas (xuy95, part of eq. 30) are:
c a1(m n j)) = Sum_{l != j, mu; nu} [
c       a(mu nu l) A(mu nu m n l j) +
c       b(mu nu l) B(mu nu m n l j)]
c and
c b1(m n j) = Sum_{l:l!=j, mu, nu} [
c       a(mu nu l) B(mu nu m n l j) +
c       b(mu nu l) A(mu nu m n l j)]
c @param nL number of particles
c @param r0
c @param nmax maximum order
c @param uvmax maximum value for uv index
c @param fint interaction index; from input file; see comment in gmm01f
c @param atr0 the translation coefficients A
c @param btr0 the translation coefficients B
c @param ek exp(i m phi); which phi?
c @param drot rotation coefficients (?)
c @param as the input scattering coefficients a(mu nu l)
c @param bs the scattering coefficients a(mu nu l)
c @param as1 the result coefficients a1(m n l)
c @param bs1 the result coefficients b1(m n l)
c @param ind polarization
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
      subroutine trans(nL,r0,nmax,uvmax,fint,atr0,btr0,ek,
     +                   drot,as,bs,as1,bs1,ind)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter (nmp=np*(np+2))
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
      parameter (nrc=4*np*(np+1)*(np+2)/3+np)
      parameter (nij=nLp*(nLp-1)/2)
      integer v,nmax(nLp),uvmax(nLp),ind(nL)
      double precision r0(6,nLp),drot(nrc,nij)
      complex*16 atr(2,np,nmp),atr0(ni0,nij),btr0(ni0,nij),
     +     ek(np,nij),as(nLp,nmp),bs(nLp,nmp),as1(nLp,nmp),
     +     bs1(nLp,nmp),at1(2,nmp),bt1(2,nmp)
          common/tran/atr
      do i=1,nL
         if(ind(i).gt.0) goto 11
             do imn=1,uvmax(i)
                as1(i,imn)=dcmplx(0.d0,0.d0)
            bs1(i,imn)=dcmplx(0.d0,0.d0)
         enddo
         do 10 j=1,nL
            if(j.eq.i) goto 10
            x0=r0(1,i)-r0(1,j)
            y0=r0(2,i)-r0(2,j)
            z0=r0(3,i)-r0(3,j)
            d=dsqrt(x0*x0+y0*y0+z0*z0)
            temp=(r0(4,i)+r0(4,j))/d
            if(temp.le.fint) goto 10
            if(i.lt.j) then
                   ij=(j-1)*(j-2)/2+j-i
                else
                   ij=(i-1)*(i-2)/2+i-j
                endif
                nlarge=max(nmax(i),nmax(j))
                itrc=0
                nsmall=min(nmax(i),nmax(j))
                do m=-nsmall,nsmall
                   n1=max(1,iabs(m))
                   do n=n1,nlarge
                      do v=n1,nlarge
                         itrc=itrc+1
                         iuv=v*(v+1)+m
                         atr(1,n,iuv)=atr0(itrc,ij)
                         atr(2,n,iuv)=btr0(itrc,ij)
                         if(x0.eq.0.d0.and.y0.eq.0.d0) then
                            if(z0.lt.0.d0) goto 20
                         else
                            if(j.lt.i) goto 20
                         endif
                         goto 21
 20                      sic=dble((-1)**(n+v))
                         atr(1,n,iuv)=sic*atr(1,n,iuv)
                         atr(2,n,iuv)=-sic*atr(2,n,iuv)
 21                      continue
                      enddo
                   enddo
                enddo
                do iuv=1,uvmax(j)
                   at1(1,iuv)=as(j,iuv)
                   at1(2,iuv)=bs(j,iuv)
                enddo
                if(x0.eq.0.d0.and.y0.eq.0.d0) then
                   call trv(at1,nmax(j),nmax(i))
                else
                   call rtr(at1,nmax(j),nmax(i),ek(1,ij),drot(1,ij))
                endif
                do imn=1,uvmax(i)
                   as1(i,imn)=as1(i,imn)+at1(1,imn)
                   bs1(i,imn)=bs1(i,imn)+at1(2,imn)
                enddo
 10	     continue
 11	     continue
      enddo
      return
      end

