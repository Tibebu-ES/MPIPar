#for each polybench - do the follwing
#	--1 parallelize and put the parallelized code in 'inputdir/mpi/benchname.c
#	--2 compile and run both serial and parallelized code ... put the bin files in - inputdir/ and /inputdir/mpi/  
#	-- 3 check if the ouputs of the serial and parallel codes are identical


for i in $(cat benchmark_list); do
	input=$i;
	inputDir=$(dirname -- "$i");
	outputDir=$inputDir/mpi
	benchname="$(basename "$input" .c)";

#compiler options 
	#dataset="-DMINI_DATASET"
	dataset="-DSMALL_DATASET"
	#related to parallelization
	# -for inlining kernel_input function in the main fuction 
	inline="-tinline=mode=2:functions=kernel_$benchname"
	alias="-alias=2" #assuming no alias
	range="-range=2"
	#include '-TIME' - in the mpipar_options to insert time instrumenting codes
	#MPI flag is included for making instrumentation by master node -DPOLYBENCH_TIME -DMPI  $dataset
	#-MPIGen
	mpipar_options="-MPIGen -TIME -outdir=$outputDir -preprocessor=cpp -C -I../utilities/ -I$inputDir -DPOLYBENCH_TIME $dataset" 
#

#Step 1 - parallelize 
#1>&- 2>&-

#java -jar MPIPar.jar $mpipar_options $input 1>&- 2>&-	

#Step 2 compile    To generate the reference output of a benchmark:

#compile serial version
#gcc -O0 -I ../utilities -I $inputDir ../utilities/polybench.c $input '-lm'  $dataset -o $inputDir/$benchname

#compile parallel version
mpicc -I ../utilities -I $inputDir ../utilities/polybench.c $outputDir/$benchname'.c' '-lm'   -o $outputDir/$benchname

#Step3  - run 
#run serial version  - dump the output in $inputDir/$benchname.out file
#$inputDir/$benchname > $inputDir/$benchname'.tm'

#run parallel version  - dump the output in $outputDir/$benchname.out file
#mpiexec -n 4 $outputDir/$benchname > $outputDir/$benchname'.tm'

	
	
			
done
