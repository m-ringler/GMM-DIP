. ./parameters.sh
d0=$(echo $DISTANCES | awk '{print $1}')
prefix=$1
joint=$prefix.dat
if [[ -e "$joint" ]]
then
    mv $joint ${joint}~
fi
cp $prefix-${d0}mu.dat $joint
for distance in $DISTANCES
do
    if [[ "$distance" != "$d0" ]]
    then
        f=$prefix-${distance}mu.dat
        if [[ -e "$f" ]]
        then
            join $joint $f > tmp.dat
            mv tmp.dat $joint
        else
            echo Missing file: $f >&2
            exit 1
        fi
    fi
done