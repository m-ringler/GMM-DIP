#!/bin/bash
# Calculates purcell factor and energy transfer factor for a dipole
# in a nanoparticle dimer embedded in a medium
# with a purely real refractive index
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
elif test -z "$R2"
then
    echo particle 1 radius missing
    exit 1
fi

for distance in $DISTANCES
do
    outputrad=wl-gr-${distance}mu.dat
    outputet=wl-get-${distance}mu.dat
    echo -e "\rcalculating rates at d = $distance micron"
    dnm=$(echo "$distance * 1000" | bc -l)
    echo "wl/nm d${dnm}" > $outputet
    echo "wl/nm d${dnm}" > $outputrad
    #calculate x
    x1=$(echo "$R1 $distance" | awk '{print ($1 + $2/2.)}')
    x2=$(echo "$R2 $distance" | awk '{print -($1 + $2/2.)}')
    for((i=$STARTWL;i<$STOPWL;i+=$WLSTEP))
        do
        # convert wavelength to micron
        wavelength="0.$i"

        # print progress
        echo -ne "\rlambda = $wavelength micron"
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

        #create the .k file for gmmf01
        # with nanoparticle dimer
        echo "$wlEff = effective wavelength (vacuum wavelength: $wavelength)" > tmp.k
        echo "3 = #spheres + 1 (dipole)" >> tmp.k
        echo "$x1 0 0 $R1 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
        echo "$x2 0 0 $R2 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
        echo "$XYZD 0.00001 1 0 = x, y, z, r, Re(m), Im(m)" >> tmp.k
        cat tmp.k >> $LOG
        ./gmmdip > $GMMDIPLOG || exit $?
        grad=$(tail -n 1 gmmdip.out | awk '{print $1}')
        get=$(tail -n 1 gmmdip.out | awk '{print $2}')

        echo "$wavelength $grad" >> $outputrad
        echo "$wavelength $get" >> $outputet
    done
done

echo -e "\njoining individual files to new file wl-gr.dat wl-get.dat"
. join.sh wl-gr
. join.sh wl-get
echo -e "\nCleaning up..."
rm field.dat
rm gmmdip.log
rm *out
rm tmp.k
echo Done

