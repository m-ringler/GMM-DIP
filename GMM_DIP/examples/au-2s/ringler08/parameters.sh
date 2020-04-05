#!/bin/bash
# $Id: parameters.sh,v 1.10 2008/02/13 20:45:05 Moritz.Ringler Exp $
# MQMie NK file for the particle material
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
PARTICLENK=au_jc_smooth.nk
#PARTICLENK=H2O_b.nk
# MQMie NK file for the medium (empty for none)
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
MEDIUMNK=h2o_b.nk
# start wavelength in nm (int)
STARTWL=540
# stop wavelength in nm (int)
STOPWL=721
# wavelength step in nm (int)
WLSTEP=1

#distances (micron)
DISTANCES="0.006422432 0.003620182 0.00352208 0.002329881 0.002122791 0.002030293 0.001825797 0.001649431 0.001617513 0.001448625 0.0013035 0.001260076 0.001129548 0.000987448 0.000920278 0.00079448"

# radius of particles (micron)
R1=0.02
R2=0.02

# coordinates of particles (micron)
#XYZ1="-0.025 0. 0."
#XYZ2="0.025 0. 0."

# coordinate of dipole
XYZD="0 0 0"

# log file
#LOG=debug.log
LOG=/dev/null
#GMMDIPLOG=gmmdip.log
GMMDIPLOG=/dev/null
