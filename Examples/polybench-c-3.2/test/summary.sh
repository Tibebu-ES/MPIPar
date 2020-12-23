#find outer loops detected in each bench output

for i in $(cat benchmark_list); do
	input=$i;
	inputDir=$(dirname -- "$i");
	benchname="$(basename "$input" .c)";
	echo $benchname" ==============================================="
	head $inputDir"/mpi/"$benchname".tm"
			
done
