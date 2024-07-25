for db in $(seq 23); do
    grep real time_*_${db}_* | sed "s/^.*real\s*//" | awk -Fm '{ print ($1 * 60) + $2 }' | sed "s/.*/${db},&/"
done
