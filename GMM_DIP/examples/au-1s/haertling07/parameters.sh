#!/bin/bash
# $Id: parameters.sh,v 1.1 2008/02/13 20:44:52 Moritz.Ringler Exp $
# MQMie NK file for the particle material
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
PARTICLENK=AU_JC_SMOOTH.nk
#PARTICLENK=H2O_b.nk
# MQMie NK file for the medium (empty for none)
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
MEDIUMNK=H2O_b.NK
# start wavelength in nm (int)
STARTWL=560
# stop wavelength in nm (int)
STOPWL=561
#STOPWL=441
# wavelength step in nm (int)
WLSTEP=1

#distances (micron)
#DISTANCES="0.0001 0.0004 0.0008 0.0010 0.0016 0.0020 0.0040 0.0080 0.012 0.016 0.017"
STARTDIST=1
STOPDIST=201
DISTSTEP=1

# radius of particles (micron)
R1=0.04

# coordinate of dipole
XYZD="0 0 0"

# log file
LOG=debug.log
LOG=/dev/null
GMMDIPLOG=gmmdip.log
#GMMDIPLOG=/dev/null
