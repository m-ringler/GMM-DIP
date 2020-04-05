Calculations of field enhancement for

Bek, Alpan; Jansen, Reiner; Ringler, Moritz; Mayilo, Sergiy; Klar, Thomas A.; Feldmann, Jochen.
Fluorescence enhancement in hot spots of AFM-designed gold nanoparticle sandwiches.
Nano Lett. 8 (2008).
http://dx.doi.org/10.1021/nl072602n

Just run gmmfield in the longitudinal/ or transverse/ subdirectory.

AVERAGE ENHANCEMENTS
can be obtained from field.dat by running

awk 'BEGIN {n=0; inte=0 ; inte2=0};
     $1*$1+$2*$2+$3*$3<=0.0004 {n++; inte+=$10; inte2+=$10*$10};
     END {print n, inte/n, inte2/n}' \
field.dat

from a (bash) shell.