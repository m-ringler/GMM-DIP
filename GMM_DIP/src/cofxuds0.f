c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  subroutine cofxuds0.f returns the two classes of vector
c  (axial) translation coefficients for a given combination of
c  (m,n,m,v) and a given dimensionless translation distance kd
cu uses subroutine gid0.f
      subroutine cofxuds0(nmax,m,n,v,sja,sya,A,B,Aj,Bj)
      
      implicit double precision (a-h,o-z)
      include 'gmm01f.par'
      parameter (ni0=np*(np+1)*(2*np+1)/3+np*np)
      parameter (ng0=np*(2*np**3+10*np**2+19*np+5)/6)
      
      integer v,p,qmax
      double precision sja(0:n+v+1),sya(0:n+v+1)
      complex*16 A,B,Aj,Bj,signz
      common/ig0/iga0(ni0)
      common/g0/ga0(ng0)
      common/cofmnv0/cof0(ni0)
      fa(m,p)=dble(-2*m*p*(p-1))
      fb(n,v,p)=dble(p*p-(n+v+1)*(n+v+1))
     +            *dble(p*p-(n-v)*(n-v))/dble(4*p*p-1)
      if(iabs(m).gt.n.or.iabs(m).gt.v) then
         write(6,*) '|m|>n or v in subroutine cofxuds0.f'
         stop
      endif
      A=0.d0
      B=0.d0
      Aj=0.d0
      Bj=0.d0
      call gid0(nmax,m,n,v,id)
      c=cof0(id)
      ig=iga0(id)
      nv2=v*(v+1)+n*(n+1)
      signz=dcmplx(0.d0,1.d0)**(n+v)
      p=n+v+2
      qmax=min(n,v)
      do i=1,qmax+1
         p=p-2
         cp=dble(nv2-p*(p+1))*ga0(ig+i)
         sj=sja(p)
         sy=sya(p)
         A=A+dcmplx(sj,sy)*signz*cp
         Aj=Aj+sj*signz*cp
         signz=-signz
      enddo
      A=A*c
      Aj=Aj*c
      if(m.eq.0) return
      signz=dcmplx(0.d0,1.d0)**(n+v+1)
      p=n+v
      do 20 i=1,qmax
         p=p-2
         signz=-signz
         if(i.eq.1) then
                cp=dble(2*p+3)*fa(m,p+3)
                cp=cp*ga0(ig+1)/dble((p+3)*(p+2))
                goto 21
             endif
             if(i.eq.qmax) then
                if(p.eq.0) goto 22
                nv2=p*(p+1)
                cp=dble(2*p+3)*fa(m,p+1)
                cp=-cp*ga0(ig+i+1)/dble(nv2)
                goto 21
             endif
 22      c4=fa(m,p+2)
             cp=-dble((p+1)*(p+2))*fb(n,v,p+2)*ga0(ig+i)
         cp=cp+dble((p+2)*(p+1))*fb(n,v,p+1)*ga0(ig+i+1)
         cp=cp*dble(2*p+3)/c4
 21      sj=sja(p+1)
         sy=sya(p+1)
         B=B+dcmplx(sj,sy)*signz*cp
         Bj=Bj+sj*signz*cp
 20   continue
      B=B*c
      Bj=Bj*c
      return
      end

