#!/bin/bash
# $Id: parameters.sh,v 1.1 2007/12/12 17:43:30 Moritz.Ringler Exp $
# MQMie NK file for the particle material
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
PARTICLENK=au_jc_smooth.nk
# MQMie NK file for the medium (empty for none)
# ALL wavelengths in this file MUST be specified to 3-digit precision (0.xxx)
#MEDIUMNK=h2o_b.nk
# start wavelength in nm (int)
STARTWL=532
# maximum distance of each particle from origin (nm)
DMAX=200
# minimum distance of each particle from origin (nm)
# must be larger than radius
DMIN=1
# distance step (nm)
DSTEP=5

# radius of particles (micron)
R1=0.04

# log file
LOG=debug.log

