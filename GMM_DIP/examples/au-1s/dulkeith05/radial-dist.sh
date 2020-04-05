#!/bin/bash
# Calculates emitted power of particles in a dielectric medium
# with a real refractive index illuminated by a dipole
export LC_ALL=C
. parameters.sh

if test -z "$XYZD"
then
    echo dipole coordinates missing
    exit 1
elif test -z "$R1"
then
    echo particle 1 radius missing
    exit 1
fi

# convert wavelength to micron
wavelength="0.$STARTWL"

# print progress
echo -ne "\rlambda = $wavelength micron"

# get the particle refractive index (real and imaginary part)
#nk=(`grep -e '^ *'$wavelength $PARTICLENK | awk '{print $2, $3}'`)
#particleN=${nk[0]}
#particleK=${nk[1]}
particleN=0.146405
particleK=3.77378
echo n k $particleN $particleK >> $LOG

# initialize effective wavelength to vacuum wavelength
wlEff=$wavelength

# get the refractive index of the medium (real part)
mediumN=1.33475
echo n_medium $mediumN >> $LOG

# calculate the wavelength in the dielectric
# wl_medium = wl_vac / n_medium
wlEff=$(echo $wavelength $mediumN | awk '{print $1/$2}')
echo wl_medium $wlEff >> $LOG

# calculate the effective complex refractive index of the particle
# n = n/n_medium
particleN=$(echo $particleN $mediumN | awk '{print $1/$2}')
particleK=$(echo $particleK $mediumN | awk '{print $1/$2}')


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


output=radial.dat
echo "d/nm grad get" > $output
for ((da=$STARTDIST;da<$STOPDIST;da+=$DISTSTEP))
do
    dnm=$(echo "$da / 10." | bc -l)
    distance=$(echo "$dnm / 1000." | bc -l)
    echo -ne "\rcalculating rates at d = $distance micron"
    echo -n "${dnm} " >> $output
    #calculate x
    x1=$(echo "$R1 $distance" | awk '{print ($1 + $2)}')

    #create the .k file for gmmf01
    # with nanoparticle dimer
    echo "$wlEff = effective wavelength (vacuum wavelength: $wavelength)" > tmp.k
    echo "2 = #spheres + 1 (dipole)" >> tmp.k
    echo "$x1 0 0 $R1 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
    echo "$XYZD 0.00001 1 0 = x, y, z, r, Re(m), Im(m)" >> tmp.k
    cat tmp.k >> $LOG
    ./gmmdip > $GMMDIPLOG || exit $?
    grad=$(tail -n 1 gmmdip.out | awk '{print $1}')
    get=$(tail -n 1 gmmdip.out | awk '{print $2}')

    #write rates to rates.dat
    echo "$grad $get" >> $output
done

echo -e "\nCleaning up..."
rm field.dat
rm *out
rm tmp.k
echo Done

