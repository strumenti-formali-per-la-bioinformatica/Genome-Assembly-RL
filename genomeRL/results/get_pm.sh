for db in $(seq 23); do
    tail results_*_${db}_*.txt |  grep dist..[0-9][0-9]* | grep -o test.rw..[0-9\.]* | sed "s/^.*test.rw../${db},/"
done
