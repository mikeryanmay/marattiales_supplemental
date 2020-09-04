#!/bin/bash

# specify default parameters
data=${data:-./data}
samples=${samples:-./output}
output=${output:-./output}
ncores=${ncores:-1}
nsim=${nsim:-1000}
overwrite=${overwrite:-false}

# read arguments
while getopts ":d:s:o:n:i:r:" opt; do
	case $opt in
		d ) # data directory
			data=$OPTARG;;
		s ) # sample directory
			samples=$OPTARG;;
		o ) # output directory
			output=$OPTARG;;
		n ) # the number of cores
			ncores=$OPTARG;;
		i ) # the number of simulation
			nsim=$OPTARG;;
		r ) # whether to overwrite
			overwrite=$OPTARG;;
		\? ) echo "Usage: cmd [-d] [-n] [-s] [-o] [-r]";;
	esac
done

##################
# check the data #
##################

echo -e "\nLocating data files."
if [ ! -d $data ]; then
	echo -e "\tProvided data directory does not exist."
	exit 1
fi

if [ ! -f "$data/morpho.nex" ]; then
	echo -e "\tNo morphological data found. Please provide a $data/morpho.nex file."
	exit 1
else
	echo -e "\tMorphological data located at $data/morpho.nex"
fi

#########################
# check the input files #
#########################

echo -e "\nLocating posterior samples."
if [ ! -d $samples ]; then
	echo -e "\tProvided sample directory does not exist."
	exit 1
fi

# get all the directories
sample_dirs=(${samples}/*/)

# check for combined samples in each directory
# only include outputs that are combined
sample_dirs_combined=()
output_dirs_combined=()
for sample_dir in ${sample_dirs[@]}; do
	if [ -f "${sample_dir}/tree_combined.nex" ] && [ -f "${sample_dir}/params_combined.log" ]; then
		sample_dirs_combined+=($sample_dir)
		output_dirs_combined+=("./${output}/$(basename $sample_dir)")
	else
		echo -e "\t Warning: no combined samples found in ${sample_dir}."
	fi
done

##########################
# perform the simulation #
##########################

parallel --eta --jobs $ncores "Rscript src/pps_sim.R --data ${data} --samples {1} --output {2} --nsim ${nsim} --overwrite ${overwrite}" ::: ${sample_dirs_combined[@]} :::+ ${output_dirs_combined[@]}
