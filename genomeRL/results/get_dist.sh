for db in $(seq 23); do
    tail results_*_${db}_*.txt |  grep -o dist..[0-9][0-9]* | sed "s/^.*dist../${db},/"
done
