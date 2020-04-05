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
c   Calculates the legendre polynomials P_n(x) = P^0_n(x) for n <= nmax
c   using Abramowitz 8.5.3 and the normalization
c   sqrt((2n+1)/(n(n+1)))
c   This normalization is consistent with pi and tau as used in gmm01f.
c   See legdre.f.nb and ../pitaud/pitaud.nb.
c
c   TESTED OK WITH tstlgd.f AGAINST legdre.f.nb
c
      subroutine legdre(nmax, xt, p)

c     defines np
      include 'gmm01f.par'

c SUBROUTINE PARAMETERS
      integer nmax
      double precision xt
      double precision p(0:np)

c LOCAL VARIABLES
      integer n, npp, nm, n2p
      double precision en

c COMPUTATION
      p(0) = 1
      p(1) = xt
c     npp = n + 1
      npp = 2
c     nm = n - 1
      nm = 0
c     n2p = 2n + 1
      n2p = 3
      do n = 1, nmax-1
          p(npp) = dble(n2p * xt * p(n) - n * p(nm)) / dble(npp)
          npp = npp + 1
          nm = nm + 1
          n2p = n2p + 2
      enddo
      do n = 1, nmax-1
          en = dsqrt(dble(2 * n + 1)/dble(n * (n + 1)))
          p(n) = p(n) * en
      enddo
      return
      end
