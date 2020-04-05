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
c     transform spherical to cartesian vector components
      subroutine sphcrtv(r, xt, sphi, cphi, vsph, vcrt)
        double precision r, xt, sphi, cphi
        complex*16 vsph(3), vcrt(3)

        double precision st

        st = dsqrt(1-xt*xt)

        vcrt(1) =   vsph(1) * cphi * st
     +            - vsph(2) * sphi
     +            + vsph(3) * cphi * xt

        vcrt(2) =   vsph(1) * sphi * st
     +            + vsph(2) * cphi
     +            + vsph(3) * sphi * xt

        vcrt(3) =   vsph(1) * xt
     +            - vsph(3) * st

      return
      end
