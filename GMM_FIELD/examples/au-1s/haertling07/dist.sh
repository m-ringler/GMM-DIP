#!/bin/bash
# Calculates field and cross section of particles in a dielectric medium
# with a real refractive index
export LC_ALL=C
. parameters.sh
mv dist-field.dat dist-field.dat~
echo "d/nm |E_x/E_in|^2" > dist-field.dat
mv dist-cro.dat dist-cro.dat~
echo "d/nm Cext Cabs Csca" > dist-cro.dat

#if test -z "$XYZ1"
#then
#    echo particle 1 coordinates missing
#    exit 1
#elif test -z "$XYZ2"
#then
#    echo particle 2 coordinates missing
#    exit 1
if test -z "$R1"
then
    echo particle 1 radius missing
    exit 1
fi

# convert wavelength to micron
wavelength="0.$STARTWL"

# print progress
echo -ne "calculating field and cross sections at $wavelength micron"
echo wavelength $wavelength >> $LOG

# get the particle refractive index (real and imaginary part)
nk=(`grep -e '^ *'$wavelength $PARTICLENK | awk '{print $2, $3}'`)
particleN=${nk[0]}
particleK=${nk[1]}
echo n k $particleN $particleK >> $LOG

# initialize effective wavelength to vacuum wavelength
wlEff=$wavelength

if test -f "$MEDIUMNK"
    then
    # get the refractive index of the medium (real part)
    mediumN=`grep -e '^ *'$wavelength $MEDIUMNK | awk '{print $2}'`
    echo n_medium $mediumN >> $LOG

    # calculate the wavelength in the dielectric
    # wl_medium = wl_vac / n_medium
    wlEff=$(echo $wavelength $mediumN | awk '{print $1/$2}')
    echo wl_medium $wlEff >> $LOG

    # calculate the effective complex refractive index of the particle
    # n = n/n_medium
    particleN=$(echo $particleN $mediumN | awk '{print $1/$2}')
    particleK=$(echo $particleK $mediumN | awk '{print $1/$2}')
elif test -n "$MEDIUMNK"
then
    echo File not found $MEDIUMNK
    exit 1
fi


echo n_eff $particleN >> $LOG
echo k_eff $particleK >> $LOG

if test -z "$particleN"
then
    echo n missing
    exit 1
elif test -z "$particleK"
then
    echo k missing
    exit 1
elif test -z "$wlEff"
then
    echo wavelength missing
    exit 1
fi


for ((dnm=DMIN; dnm<=DMAX; dnm+=DSTEP))
do

    # print progress
    distance=$(echo "$dnm / 1000." | bc -l)
    echo -ne "\rcalculating field and cross sections at distance $distance micron"
    #dnm=$(echo "$distance * 1000" | bc -l)
    echo distance $distance >> $LOG
    #calculate y
    D1=$(echo "$R1 $distance" | awk '{print ($1 + $2)}')

    #create the .k file for gmmf01
    echo "$wlEff = effective wavelength (vacuum wavelength: $wavelength)" > tmp.k
    echo "1 = number of spheres" >> tmp.k
    echo "$D1 0 0 $R1 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
    cat tmp.k >> $LOG

    ./gmmfield > gmm01f.log || exit $?
    echo -n "$dnm " >> dist-field.dat
    awk '{print $4 * $4 + $5 * $5}' field.dat  >>  dist-field.dat
    echo -n "$dnm " >> dist-cro.dat
    #grep -m 1 -e '^ *x' gmm01f.out  | awk '{print $2, $3, $4}' >> dist-cro.dat
done
echo -e "\nCleaning up..."
rm field.dat
rm gmm01f.log
rm *out
rm tmp.k
echo Done

