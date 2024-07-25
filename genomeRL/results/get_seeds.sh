for db in $(seq 23); do
     cat seeds_*_${db}_*.txt | sed "s/.*/${db},&/"
done
