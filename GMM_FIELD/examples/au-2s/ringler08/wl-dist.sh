#!/bin/bash
# Calculates cross sections, field enhancement, phase of total field
# for a dimer.
export LC_ALL=C
. parameters.sh

if test -z "$R1"
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
    outputcsca=wl-csca-${distance}mu.dat
    outputcabs=wl-cabs-${distance}mu.dat
    outputcext=wl-cext-${distance}mu.dat
    outputfield=wl-absfield-${distance}mu.dat
    outputphase=wl-phase-${distance}mu.dat
    echo -e "\rcalculating cross sections and field at d = $distance micron"
    dnm=$(echo "$distance * 1000" | bc -l)
    echo "wl/nm d${dnm}" > $outputcsca
    echo "wl/nm d${dnm}" > $outputcabs
    echo "wl/nm d${dnm}" > $outputcext
    echo "wl/nm d${dnm}" > $outputfield
    echo "wl/nm d${dnm}" > $outputphase
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
        fi

        #create the .k file for gmm01f
        # with nanoparticle dimer
        echo "$wlEff = effective wavelength (vacuum wavelength: $wavelength)" > tmp.k
        echo "2 = #spheres" >> tmp.k
        echo "$x1 0 0 $R1 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
        echo "$x2 0 0 $R2 $particleN $particleK = x, y, z, r, Re(m), Im(m)" >> tmp.k
        cat tmp.k >> $LOG
        ./gmmfield > $GMMLOG || exit $?

        #extract rates
        cro=$(grep -m 1 -e '^ *x' gmm01f.out)
        cext=$(echo $cro  | awk '{print $2}')
        cabs=$(echo $cro  | awk '{print $3}')
        csca=$(echo $cro  | awk '{print $4}')

        #write rates to wl-rates.dat
        echo "$wavelength $csca" >> $outputcsca
        echo "$wavelength $cext" >> $outputcext
        echo "$wavelength $cabs" >> $outputcabs

        #write fields to wl-field.dat
        norme=$(awk '{print $10}' field.dat)
        phase=$(awk '{print atan2($5, $4)}' field.dat)
        echo "$wavelength $norme" >> $outputfield
        echo "$wavelength $phase" >> $outputphase
    done
done

echo -e "\njoining cro files to new file wl-csca.dat wl-cext.dat wl-cabs.dat"
. join.sh wl-csca
. join.sh wl-cext
. join.sh wl-cabs
echo -e "\njoining field files to new file wl-absfield.dat"
. join.sh wl-absfield
. join.sh wl-phase
echo -e "\nCleaning up..."
rm field.dat
#rm gmmdip.log
rm *mu.dat
rm *out
rm tmp.k
echo Done

