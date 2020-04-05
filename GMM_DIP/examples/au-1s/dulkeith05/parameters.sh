#!/bin/bash
# $Id: parameters.sh,v 1.3 2008/02/13 20:44:48 Moritz.Ringler Exp $
# MQMie NK file for the particle material
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
PARTICLENK=AU_JC_SMOOTH.nk
#PARTICLENK=H2O_b.nk
# MQMie NK file for the medium (empty for none)
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
MEDIUMNK=H2O_b.NK
# start wavelength in nm (int)
STARTWL=668
# stop wavelength in nm (int)
STOPWL=669
#STOPWL=441
# wavelength step in nm (int)
WLSTEP=1

#Angstroem
STARTDIST=1
STOPDIST=171
DISTSTEP=1

# radius of particles (micron)
R1=0.008

# coordinate of dipole
XYZD="0 0 0"

# log file
LOG=debug.log
LOG=/dev/null
GMMDIPLOG=gmmdip.log
#GMMDIPLOG=/dev/null
