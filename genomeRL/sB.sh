for iter in $(seq 20); do
     { time ./sA.sh $1 $2 $iter >/dev/null; } &> time_$1_$2_$iter
done;
