c gmm01f
c  lnfacd.f  (double precision function)
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  returns ln(z!)  z>-1.0
c  based on Lanczos' method [see Xu, Journal of Computational
c  Physics, v.139, 137-165 (1998)]
      double precision function lnfacd(z)
      integer i
      double precision z,a,b,cp,c0(11)
      data c0/0.16427423239836267d5, -0.48589401600331902d5,
     +        0.55557391003815523d5, -0.30964901015912058d5,
     +        0.87287202992571788d4, -0.11714474574532352d4,
     +        0.63103078123601037d2, -0.93060589791758878d0,
     +        0.13919002438227877d-2,-0.45006835613027859d-8,
     +        0.13069587914063262d-9/
      a=1.d0
      cp=2.5066282746310005d0
      b=z+10.5d0
      b=(z+0.5d0)*dlog(b)-b
      do i=1,11
        z=z+1.d0
        a=a+c0(i)/z
      enddo
      lnfacd=b+dlog(cp*a)
      return
      end
