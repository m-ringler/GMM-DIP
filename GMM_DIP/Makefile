
#******************************************************************************
#
# 1. Compiler Definition
#    ===================
#    FC specifies the fortran compiler
#    FFLAGS is used to compile the modules in "OBJS"
#    LDFLAGS are flags for linking
#
#	timeit_titan
#
#----------------------Linux with g77 -------- --------------------------------

FC		= gfortran
FFLAGS		= -c -u -O
LDFLAGS 	=

#******************************************************************************
#
#                  End of option specifications.
#      User should not need to edit the Makefile below this point
#
#******************************************************************************

SRCS	      = src/gmmdip.f\
		src/field.f\
		src/abMiexud.f\
		src/besseljd.f\
		src/besselyd.f\
		src/carsphd.f\
		src/cofd0.f\
		src/cofnv0.f\
		src/cofsrd.f\
		src/cofxuds0.f\
		src/gau0.f\
		src/gxurcd0.f\
		src/legdre.f\
		src/lnfacd.f\
		src/mueller.f\
		src/normlz.f\
		src/orientcd.f\
		src/orientud.f\
		src/pitaud.f\
		src/ran1d.f\
		src/rotcoef.f\
		src/rtr.f\
		src/sphcrtv.f\
		src/tipitaud.f\
		src/trans.f\
		src/trv.f\
		src/vswf.f

OBJS = $(SRCS:src/%.f=%.o)

%.o: src/%.f
		$(FC) $(FFLAGS) $<

all:		gmm01f

clean:
		rm -f *.o

gmm01f:		$(OBJS) Makefile
		$(FC) $(OBJS) $(LDFLAGS)  -o gmmdip

