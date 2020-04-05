c field.f
c subroutine field.f
c
c     (c) Photonics & Optoelectronics Group, LMU München, 2007-2008
c     @author Moritz Ringler
c
c     This file is part of the field extension to Yu-lin Xu's GMM
c     by Moritz Ringler.
c     This program is free software: you can redistribute it and/or
c     modify it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3 of
c     the License, or (at your option) any later version.
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c     You should have received a copy of the GNU General Public License
c     along with this program.  If not, see
c     <http://www.gnu.org/licenses/>.
c
c
c die Idee: benutze die partiellen Streukoeffizienten, die von
c gmm01f ausgespuckt werden, um die gestreuten Felder der einzelnen
c Kugeln an jedem Punkt zu berechnen
c Diese Felder sind zunaechst als Komponenten im sphaerischen
c Koordinatensystem mit Mittelpunkt in der betreffenden Kugel gegeben.
c Sie werden dann ins kartesische Laborsystem umgerechnet und zum
c einfallenden Feld dazuaddiert. Fertig ist das ganze Feld.
c
c  #L index of sphere
c  a[m,n,L], b[m,n,L] = calculateScatteringCoefficients(N_max)
c  e[m, n] = calculateNormalizationFactor(N_max)
c
c  for x, y, z
c      ETotal[] = incidentField(x,y,z)
c      for L < NumSpheres
c          r,phi,theta = toSphericalCoordinates(x,y,z, origin[L])
c          N3[m,n],M3[m,n] = calculateVectorSphericalHarmonics(N_max,
c                           r, phi, theta)
c          for n, m, i
c              EScaL[i] += I  e[m,n]
c                    ( a[m,n,L] N3[m,n,i] + b[m,n,L] M3[m,n,i] )
c          ETotal += toCartesianVector(EScaL[], origin[L])
c      write x, y, z, ETotal[], |ETotal[]|^2
c
c START OF PROGRAM
c This program assumes that the incident field vector points
c in the x direction (beta == theta == 0).
c Output is in units of the incident field amplitude.
c @param nL
c    an integer holding the number of spheres
c @param r0
c    a double array holding six coordinates for each sphere
c    r0(1-3, i) hold the origin of the ith sphere (x, y, z)
c    r0(4, i) is the radius of the ith sphere
c    r0(5, i) is the real part of the refr. index of the ith sphere
c    r0(6, i) is the imaginary part
c @param k
c    the wavenumber of the incident wave = 2 Pi/ lambda
c @param nmax
c    an integer array. nmax(i) holds the highest scattering
c    order for the ith sphere (cannot be larger than np)
c @param as
c    a complex array holding the partial scattering coefficients:
c    as(i, j) is the scattering coefficient A_m,n,i of the ith sphere
c    where j is calculated from m, n as j = (m + n) + n^2; here
c    n runs from 1 to nmax(i), and m from -n to +n for each n
c @param bs
c    the partial scattering coefficients B_m,n,i
      subroutine field(nL,r0,k,nmax,as,bs)
      implicit logical (a-z)

c COMPILE TIME PARAMETERS
      integer np, nLp
      include 'gmm01f.par'

c     the # of points along each grid dimension;
c     must be odd and gt 2
      integer ngrd(3)
c     the minimum value of x, y, and z;
c     specified in the same unit as wavelength (micron)
      double precision grdmin(3)
c     the spacing between two gridpoints in one dimension;
c     same unit as grdmin
      double precision grdstp(3)

      integer nmp, nmp0
      parameter (nmp=np*(np+2),nmp0=(np+1)*(np+4)/2)

c     the complex numbers i and 0
      complex*16 cplxi, cplx0
      parameter (cplxi = (0.0D0, 1.0D0))
      parameter (cplx0 = (0.0D0, 0.0D0))

c SUBROUTINE ARGUMENTS
c see @param comments above
      integer nL
      integer nmax(nLp)
      double precision r0(6,nLp)
      double precision k
      complex*16 as(nLp,nmp),bs(nLp,nmp)

c FUNCTIONS
      double precision vnorm

c LOCAL VARIABLES
c     the index of the current sphere
      integer i

c     indexes for pairs (m,n)
      integer m, n, imn, j

c     the current integer grid position
      integer ix, iy, iz

c     the current cartesian grid position
      double precision x, y, z

c     the current grid position in spherical coordinates
c     r is the radius sqrt(x^2 + y^2 + z^2)
c     xt = cos(theta) = z/r
c     sphi = sin(phi) = y/rho; cphi = cos(phi) = z/rho
c     with rho = sqrt(x^2 + y^2)
      double precision r, xt, sphi, cphi
      complex*16 eiphi


c     the total field at x,y,z (cartesian components)
c     this is what we want to calculate
      complex*16 etot(3)

c     the scattered field of a single sphere (spherical components)
      complex*16 escati(3)
c     the scattered field of a single sphere (cartesian components)
      complex*16 escatc(3)

c     the normalized spherical harmonic wave functions
c     as vectors in spherical coordinates
c     for m >= 0 at the current grid point
c     Nmn3(1, imn) = <Nmn3, er>
c     Nmn3(2, imn) = <Nmn3, ephi>
c     Nmn3(3, imn) = <Nmn3, etheta>
c     imn is calculated from m, n as
c     imn = (n-1)(n+2)/2 + m + 1; here
c     n runs from 1 to nmax, and m from 0 to +n for each n,
c     the first (radial) component of Mmn3 is always zero
      complex*16 Mmn3(2:3,nmp0), Nmn3(3,nmp0)

c    the normalization factors Emn (eq. 3 from the gmm01f.f manual)
      complex*16 Emn(nmp0)

c    the factor i Emn a^i_mn or
c               i Emn b^i_mn
      complex*16 cimn

c    calculate normalization factors
      call normlz(np, Emn)

c    read grid
      open(UNIT=14,FILE='grid.in',STATUS='OLD')
      read(14,*) ngrd(1), ngrd(2), ngrd(3)
      read(14,*) grdmin(1), grdmin(2), grdmin(3)
      do i = 1,3
        grdstp(i) = -grdmin(i) * 2.0D0 / (ngrd(i) - 1.0D0)
      enddo
      close(14)

      open(UNIT=14,FILE='grid.out',STATUS='unknown')
      write (14,*) "number of grid points (nx, ny, nz): ",
     +             ngrd(1), ngrd(2), ngrd(3)
      write (14,*) "grid corner (x0, y0, z0): ",
     +             grdmin(1), grdmin(2), grdmin(3)
      write (14,*) "grid step (dx, dy, dz): ",
     +               grdstp(1), grdstp(2), grdstp(3)
      close(14)

c    open output file
      open(67,file='field.dat',status='unknown')

c     loop over grid positions
      z = grdmin(3)
      do iz=1,ngrd(3)
         write (*,*) iz, "/", ngrd(3)
         y = grdmin(2)
         do iy=1,ngrd(2)
             x = grdmin(1)
             do ix=1,ngrd(1)




c                 calculate incident field as
c                 Einc = \hat(ex) exp(i k z) =
c                 { cos(k z) + i sin(k, z), 0, 0 }
c                  etot(1) = dcmplx( dcos(k * z), dsin(k * z) )
                  etot(1) = cplx0
                  etot(2) = cplx0
                  etot(3) = cplx0

c                 loop over spheres
                  do i = 1, nL

c                     calculate the current grid position
c                     in a spherical coordinate system centered at the
c                     ith sphere using subroutine carsphd from gmm01f
                      call carsphd(x - r0(1, i),
     +                             y - r0(2, i),
     +                             z - r0(3, i), r, xt, sphi, cphi)

c                     calculate exp(i phi)
                      eiphi = dcmplx(cphi, sphi)

c                     if we are in the interior of one of our spheres,
c                     set the field to zero
c                     break the innermost loop (over spheres),
c                     and proceed to the next point
                      if(r < r0(4, i)) then
                          etot(1) = cplx0
                          etot(2) = cplx0
                          etot(3) = cplx0
                          goto 11
                      endif

c                     calculate the vector spherical harmonics
c                     Nmn3 and Mmn3 in spherical coordinates
c                     centered at the ith sphere
                      call vswf(nmax(i), k*r, xt, sphi, cphi,
     +                          Nmn3, Mmn3)

c                     perform the sum in eq. 11a of the manual
                      escati(1) = cplx0
                      escati(2) = cplx0
                      escati(3) = cplx0
                      imn = 1
                      j = 0
c                     loop over degree of VSWF
                      do n = 1,nmax(i)
c                         loop over negative order m
                          imn = imn + n
                          do m=-n,-1
                          j = j + 1
c                         from identities A3 (2) and A4 (2)
c                         of Xu Appl Opt 1995
c                         and eq 6a of the manual we have
c                         pi[-m,n] = (-1)^m+1 (n - m)!/(n + m)! pi[m,n]
c                         P[-m,n] = (-1)^m (n - m)!/(n + m)! P[m,n]
c                         tau[-m,n] = (-1)^m (n - m)!/(n + m)! tau[m,n]
c
c                         and from eq 3 of the manual
c                         E[-m,n] = E[m, n] (n + m)!/(n - m)!
c
c                        and so we get
c                        pi[-m,n] E[-m,n] = (-1)^(m+1) E[m,n] pi[m,n]
c                        tau[-m,n] E[-m,n] = (-1)^m E[m,n] tau[m,n]
c                        P[-m,n] E[-m,n] = (-1)^m E[m,n] P[m,n]
c
c                        plugging this into eq 5a of the manual
c                        we have
c                        E[m,n] M3[n,m] =
c                            E[-m, n] (-1)^m exp(i 2 m phi)
c                            {0, M3[-m,n]_phi, -M3[-m,n]_theta}
c                        and hence
c                        nM3[n, m] = (-1)^m exp(i 2 m phi)
c                            {0, nM3[-m,n]_phi, -nM3[-m,n]_theta}
c
c                        and with eq. 5b of the manual we have
c                        E[m,n] N3[m,n] =
c                           E[-m,n] (-1)^m exp(i 2 m phi)
c                           {N3[-m,n]_r, - N3[-m,n]_phi, N3[m,n]_theta}
c
                          cimn = cplxi * as(i, j) * Emn(imn) * (-1)**m *
     +                               eiphi**(2 * m)
                          escati(1) = escati(1) + cimn * Nmn3(1,imn)
                          escati(2) = escati(2) - cimn * Nmn3(2,imn)
                          escati(3) = escati(3) + cimn * Nmn3(3,imn)

                          cimn = cplxi * bs(i, j) * Emn(imn) * (-1)**m *
     +                               eiphi**(2 * m)
c                         radial component of Mmn3 is zero
c                         escati(2) = escati(1) + cimn * Mmn3(1,imn)
                          escati(2) = escati(2) + cimn * Mmn3(2,imn)
                          escati(3) = escati(3) - cimn * Mmn3(3,imn)

                          imn = imn - 1
                          enddo

c                         loop over positive order m
                          do m=0,n
                          j = j + 1

                          cimn = cplxi * Emn(imn) * as(i, j)
                          escati(1) = escati(1) + cimn * Nmn3(1,imn)
                          escati(2) = escati(2) + cimn * Nmn3(2,imn)
                          escati(3) = escati(3) + cimn * Nmn3(3,imn)

                          cimn = cplxi * Emn(imn) * bs(i, j)
c                         radial component of Mmn3 is zero
c                         escati(1) = escati(1) + cimn * Mmn3(1,imn)
                          escati(2) = escati(2) + cimn * Mmn3(2,imn)
                          escati(3) = escati(3) + cimn * Mmn3(3,imn)

                          imn = imn + 1
                          enddo
                      enddo

c                     at this point we have the scattered field
c                     from the ith sphere in escati
c                     we need to build its cartesian representation
                      call sphcrtv(r, xt, sphi, cphi, escati, escatc)

c                     and add it to etot
                      etot(1) = etot(1) + escatc(1)
                      etot(2) = etot(2) + escatc(2)
                      etot(3) = etot(3) + escatc(3)

                  enddo
   11             continue
c                 OUTPUT
                  write(67,'(3e13.5, 7e13.5)') x, y, z, dble(etot(1)),
     +                                 dimag(etot(1)),
     +                                 dble(etot(2)),dimag(etot(2)),
     +                                 dble(etot(3)),dimag(etot(3)),
     +                                 vnorm(etot)

                  x = x + grdstp(1)
              enddo
              y = y + grdstp(2)
          enddo
          z = z + grdstp(3)
      enddo
      close(67)
      return
      end

      double precision function vnorm(vec)
      complex*16 vec(3)
      double precision tmp
      integer i

      vnorm = 0
      do i=1,3
        tmp = abs(vec(i))
        vnorm = vnorm + tmp * tmp
      enddo
      vnorm = dsqrt(vnorm)
      return
      end
