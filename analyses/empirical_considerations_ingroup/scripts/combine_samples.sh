#!/usr/bin/bash

# specify default parameters
output=${output:-./output}
ncores=${ncores:-1}
overwrite=${overwrite:-0}

# find named parameters
while getopts ":o:n:f:" opt; do
	case $opt in
		o ) # output directory
			output=$OPTARG;;
		n ) # number of cores
			ncores=$OPTARG;;
		f ) # force recomputation
			overwrite=$OPTARG;;
		\? ) echo "Usage: cmd [-o] [-n] [-f]";;
	esac
done

# combine post-burnin samples
if [ "$overwrite" = "1" ]; then
	Rscript src/mcmc_diagnosis_combine_samples.R --output $output --ncores $ncores --reset;
else
	Rscript src/mcmc_diagnosis_combine_samples.R --output $output --ncores $ncores;
fi
