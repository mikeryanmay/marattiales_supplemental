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

# MCMC diagnosis for parameters
if [ "$overwrite" = "1" ]; then
	Rscript src/mcmc_diagnosis_params.R --output $output --ncores $ncores --reset;
else
	Rscript src/mcmc_diagnosis_params.R --output $output --ncores $ncores;	
fi

# MCMC diagnosis for trees
if [ "$overwrite" = "1" ]; then
	Rscript src/mcmc_diagnosis_trees.R --output $output --ncores $ncores --reset;
else
	Rscript src/mcmc_diagnosis_trees.R --output $output --ncores $ncores;
fi

# compile a table of ess for each run
Rscript src/mcmc_diagnosis_compile_diagnostics.R --output $output;

# create LTT plots
if [ "$overwrite" = "1" ]; then
	Rscript src/mcmc_diagnosis_ltt.R --output $output --ncores $ncores --reset;
else
	Rscript src/mcmc_diagnosis_ltt.R --output $output --ncores $ncores;
fi

# create MDS plots
if [ "$overwrite" = "1" ]; then
	Rscript src/mcmc_diagnosis_mds.R --output $output --ncores $ncores --reset;
else
	Rscript src/mcmc_diagnosis_mds.R --output $output --ncores $ncores;
fi
