c gmm01f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  transforms the rectangular coordinates (x,y,z)
c  to spherical coordinates (r,theta,phi)
      subroutine carsphd(x,y,z,r,xt,sphi,cphi)
      double precision x,y,z,r,xt,sphi,cphi
      r=dsqrt(x*x+y*y+z*z)
      if(r.eq.0.d0) then
         xt=1.d0
         sphi=0.d0
         cphi=1.d0
         return
      endif
      xt=z/r
      if(y.eq.0.d0.and.x.eq.0.d0) then
         sphi=0.d0
         cphi=1.d0
         return
      endif
      sphi=dsqrt(x*x+y*y)
      cphi=x/sphi
      sphi=y/sphi
      return
      end   

