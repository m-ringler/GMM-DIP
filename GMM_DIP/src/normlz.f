c
c     (c) Photonics & Optoelectronics Group, LMU M�nchen, 2007-2008
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
c  calculate Emn normalization factors (eq 3 from the manual
c  divided by E0)
c  for n = 1 ... nmax, m = 0...n
c  tested OK with tstnml.f against testnormlz.nb
      subroutine normlz(nmax, Emn)
      implicit logical (a-z)

      complex*16 cplxi
      parameter (cplxi = (0.0D+0, 1.0D+0))

      integer nmax
      integer n, m, imn
      complex*16 Emn(*)
      complex*16 cin
c      double precision en, emnd, lnfacd

      imn = 0
      cin = dcmplx(1.0D+0, 0.0D+0)
      do n=1,nmax
          cin = cin * cplxi
c          en = dsqrt(dble(2 * n + 1)/dble(n * (n + 1)))
          do m=0,n
            imn = imn + 1
c            emnd = en * dexp(0.5D+0 *
c     +                (lnfacd(dble(n - m)) - lnfacd(dble(n + m))))
c Matthew Arnold - corrected Emn based on normalization of tau & pi
c            emnd = 1d0
c            Emn(imn) = dcmplx(emnd, 0) * cin
            Emn(imn) = cin
          enddo
      enddo

      return
      end
