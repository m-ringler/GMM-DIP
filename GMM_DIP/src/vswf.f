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
c     Calculates the vector spherical harmonics Nmn3 and Mmn3
c     in spherical coordinates at a given point r, phi, theta
c     for all non-negative m up to degree nmax
c     @param nmax
c        the maximum degree of the VSWF to calclate
c     @param kr
c        the dimensionless distance from the origin calculated
c        as k * r = 2 Pi r / lambda
c     @param xt
c        cos(theta)
c     @param sphi
c        sin(phi)
c     @param cphi
c        cos(phi)
c     @param Nmn3, Mmn3
c        output arrays for
c        the spherical harmonic wave functions in spherical coordinates
c        for m >= 0 at the specified point
c        Mmn3(1, imn) = <Mmn3, er>
c        Mmn3(2, imn) = <Mmn3, ephi>
c        Mmn3(3, imn) = <Mmn3, etheta>
c        imn is calculated from m, n as
c        imn = (n-1)(n+2)/2 + m + 1; here
c        n runs from 1 to nmax, and m from 0 to +n for each n,
      subroutine vswf(nmax, kr, xt, sphi, cphi, Nmn3, Mmn3)
      implicit logical (a-z)

      integer np, nLp
      include 'gmm01f.par'
      integer nmp, nmp0
      parameter (nmp=np*(np+2),nmp0=(np+1)*(np+4)/2)
      complex*16 cplxi
      parameter (cplxi = (0.0D0, 1.0D0))

c SUBROUTINE ARGUMENTS
c INPUTS
      double precision kr, xt, sphi, cphi
      integer nmax
c OUTPUTS
      complex*16 Mmn3(2:3,nmp0), Nmn3(3,nmp0)

c LOCAL VARIABLES
c     arrays holding pi_n,m(cos theta) and tau_n,m(cos theta)
c     here, n = 1, ..., nmax + 1
c     and   m = 0, ..., n for each n
c     for a given n,m:
c       pi(n,m) = pi(imn)
c     with
c       imn = (n-1)(n+2)/2 + m + 1
      double precision pi(nmp0),tau(nmp0)

c     single values of pi and tau for a given n, m
      double precision pimn, taumn

c     values of exp(i phi) and exp(i m phi)
      complex*16 eiphi
      complex*16 eimphi

c     an array holding the spherical bessel function of the first kind j_n(x)
c     at x = kr, for n between 0 and nmax
      double precision besj(0:np)

c     an array holding the spherical bessel function of the second kind y_n(x)
c     at x = kr, for n between 0 and nmax
      double precision besy(0:np)

c     the spherical bessel function of the 3rd kind or
c     sph. Hankel function of  the 1st kind
c     h1[n] = jn(x) + yn(x)
      complex*16 hankln

c     the derivative psi_n' = d/dz( z h1[n](z) ) at z = kr
      complex*16 psinpr

c     the value of the associated legendre function P^m_n(xt)
      double precision pmn

c     the legendre function P_n(z) at given cos(theta) = xt
c     for n = 0 ... nmax
      double precision p(0:np)

c     loop variables
      integer n, m, imn

c END OF DECLARATIONS

c     calculate pi and tau
      call pitaud(nmax, xt, pi, tau)

c     calculate legendre polynomials P_n(xt)
      call legdre(nmax, xt, p)

c     calculate bessel functions
      call besseljd(nmax, kr, besj)
      call besselyd(nmax, kr, besy)

      eiphi = dcmplx(cphi, sphi)

c     calculate vector spherical harmonic Mmn3, Nmn3
      imn = 0

c     loop over degree
      do n=1,nmax

c        calculate the hankel function of the first kind
         hankln = dcmplx(besj(n), besy(n))
c         write (*,*) "hankl", n, hankln

c        calculate the derivative psi_n'(z) = d/dz psi_n(z)
c        of psi_n(z) = z h1[n](z)
c        using Abramowitz 10.1.21
c       (z h1[n](z))' = z h1[n-1](z) - n h1[n](z)
        psinpr = dcmplx(kr * besj(n - 1) - n * besj(n),
     +                   kr * besy(n - 1) - n * besy(n))

c        write (*,*) "psi'", n, psinpr

c        loop over order
         do m=0,n
c           calculate exp(i m phi)
            eimphi = eiphi**m
c           increment n,m array index
            imn = imn + 1

c           assign values of pi_nm(xt) and tau_nm(xt) for this m and n
            pimn = pi(imn)
            taumn = tau(imn)

c           calculate the associated legendre polynomial P^m_n(xt)
c           needed for Nmn3
c           using equation (4) from the gmm01f.f manual if m != 0
c           sin(theta) = sqrt(1 - cos^2(theta))
c           or the previously calculated value of P_n(xt) if m == 0
            if(m.eq.0) then
              pmn = p(n)
            else
              pmn = pimn * dsqrt(1 - xt*xt)/m
            endif

c           calculate Mnm3(imn) using eq. 5a from the manual
c           with h1[n] instead of j[n]
c           Mmn3(1, imn) = 0
            Mmn3(2, imn) = - taumn * hankln * eimphi
            Mmn3(3, imn) = cplxi * pimn * hankln * eimphi

c           calculate Nnm3(imn) using eq. 5b from the manual
c           with h1[n] instead of j[n]
            Nmn3(1, imn) = n * (n + 1) * pmn * hankln * eimphi/kr
            Nmn3(2, imn) = cplxi * pimn * psinpr * eimphi/kr
            Nmn3(3, imn) = taumn * psinpr * eimphi/kr

         enddo
      enddo

      return
      end


