c gmm01f
c  subroutine mueller.f
c  THIS FILE IS PART OF THE ORIGINAL GMM01F DISTRIBUTION
c  returns the values of 16 elements of the 4x4
c  Mueller matrix calculated from the known 2x2
c  amplitude scattering matrix
c  using Bohren & Huffman's formulas (p.65)
      subroutine mueller(s1,s2,s3,s4,s)
      double precision s(4,4),s1s,s2s,s3s,s4s
      complex*16 s1,s2,s3,s4
      complex*16 s2s3c,s1s4c,s2s4c,s1s3c,s1s2c
      complex*16 s3s4c,s2s1c,s4s3c,s2cs4,s3cs1
      s1s=cdabs(s1)**2
      s2s=cdabs(s2)**2
      s3s=cdabs(s3)**2
      s4s=cdabs(s4)**2
      s2s3c=s2*dconjg(s3)
      s1s4c=s1*dconjg(s4)
      s2s4c=s2*dconjg(s4)
      s1s3c=s1*dconjg(s3)
      s1s2c=s1*dconjg(s2)
      s3s4c=s3*dconjg(s4)
      s2s1c=s2*dconjg(s1)
      s4s3c=s4*dconjg(s3)
      s2cs4=dconjg(s2)*s4
      s3cs1=dconjg(s3)*s1
      s(1,1)=0.5d0*(s1s+s2s+s3s+s4s)
      s(1,2)=0.5d0*(s2s-s1s+s4s-s3s)
      s(1,3)=s2s3c+s1s4c
      s(1,4)=dimag(s2s3c-s1s4c)
      s(2,1)=0.5d0*(s2s-s1s-s4s+s3s)
      s(2,2)=0.5d0*(s2s+s1s-s4s-s3s)
      s(2,3)=s2s3c-s1s4c
      s(2,4)=dimag(s2s3c+s1s4c)
      s(3,1)=s2s4c+s1s3c
      s(3,2)=s2s4c-s1s3c
      s(3,3)=s1s2c+s3s4c
      s(3,4)=dimag(s2s1c+s4s3c)
      s(4,1)=dimag(s2cs4+s3cs1)
      s(4,2)=dimag(s2cs4-s3cs1)
      s(4,3)=dimag(s1s2c-s3s4c)
      s(4,4)=s1s2c-s3s4c
      return
      end
