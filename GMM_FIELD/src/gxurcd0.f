c gmm01f
c  gxurcd0.f to compute Gaunt coefficients a(-m,n,m,v,p)
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
cu uses lnfacd.f to compute ln(z!)
      subroutine gxurcd0(m,n,v,qmax,na)
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
      integer v,qmax,p
      double precision lnfacd,cnv(np,np),ga0(ng0)
      common/cfnv/cnv
      common/g0/ga0
      fb(n,v,p)=dble(p-(n+v+1))*dble(p+(n+v+1))
     +           *dble(p-(n-v))*dble(p+(n-v))
     +           /(dble(2*p+1)*dble(2*p-1))
      if(iabs(m).gt.n.or.iabs(m).gt.v) then
         write(6,*) 'warning: |m|>n or v in gxurcd0'
         qmax=-1
         return
      endif
      qmax=min(n,v)
      nq=qmax+1
      if(n.le.v) then
         c1=cnv(n,v)
      else
         c1=cnv(v,n)
      endif
      c1=c1-lnfacd(dble(n-m))-lnfacd(dble(v+m))
      ga0(na+1)=dexp(c1)
      if(qmax.lt.1) return
      p=n+v
      do 8 i=2,nq
         p=p-2
         if(m.eq.0) then
            c1=fb(n,v,p+1)
            c2=fb(n,v,p+2)
            goto 2
         endif
         c1=fb(n,v,p+1)
         c2=dble(4*m*m)+fb(n,v,p+2)+fb(n,v,p+3)
         if(i.eq.2) goto 2
         c3=-fb(n,v,p+4)
         goto 4
  2	     ga0(na+i)=c2*ga0(na+i-1)/c1
         goto 8
  4	     ga0(na+i)=(c2*ga0(na+i-1)+c3*ga0(na+i-2))/c1
  8       continue
      return
      end
      
      subroutine gid0(nmax,m,n,iv,id)
      nt=nmax*(nmax+1)*(2*nmax+1)/3+nmax*nmax
      ns=max(1,iabs(m))
      nc0=nmax-iabs(m)
      id=nc0*(nc0+1)*(2*nc0+1)/6
      if(m) 10,11,12
 10	  id=id+(n-ns)*(nc0+1)+iv-ns+1
      return
 11	  id=id+(n-ns)*nmax+iv
      return
 12	  id=id+(nmax-n)*(nc0+1)+nmax-iv
      id=nt-id
      return
      end

