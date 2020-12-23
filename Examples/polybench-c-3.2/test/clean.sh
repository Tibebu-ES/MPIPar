#find outer loops detected in each bench output

for i in $(cat benchmark_list); do
	input=$i;
	inputDir=$(dirname -- "$i");
	benchname="$(basename "$input" .c)";
	rm -r $inputDir"/mpi"
			
done
