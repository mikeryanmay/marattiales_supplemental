#!/bin/bash

# specify default parameters
data=${data:-./data}
modules=${modules:-./modules}
analysis=${analysis:-/modules/analysis/MCMC.Rev}
template=${template:-./modules/templates/template.Rev}
output=${output:-./output}
nruns=${nruns:-4}
jobdir=${jobdir:-./jobs}

# read arguments
while getopts ":d:m:a:t:o:n:j:" opt; do
	case $opt in
		d ) # data directory
			data=$OPTARG;;
		m ) # modules directory
			modules=$OPTARG;;
		a ) # analysis file
			analysis=$OPTARG;;
		t ) # template file
			template=$OPTARG;;
		o ) # output directory
			output=$OPTARG;;
		n ) # number of runs per analysis
			nruns=$OPTARG;;
		j ) # directory for the jobs
			jobdir=$OPTARG;;
		\? ) echo "Usage: cmd [-d] [-m] [-a] [-t] [-o] [-n] [-j]";;
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

if [ ! -f "$data/molecular.nex" ]; then
	echo -e "\tNo molecular data found. Please provide a $data/molecular.nex file."
	exit 1
else
	echo -e "\tMolecular data located at $data/molecular.nex"
fi

if [ ! -f "$data/morpho.nex" ]; then
	echo -e "\tNo morphological data found. Please provide a $data/morpho.nex file."
	exit 1
else
	echo -e "\tMorphological data located at $data/morpho.nex"
fi

if [ ! -f "$data/taxa.tsv" ]; then
	echo -e "\tNo taxon data found. Please provide a $data/taxa.tsv file."
	exit 1
else
	echo -e "\tTaxon data located at $data/taxa.tsv"
fi

if [ ! -f "$data/epochs.csv" ]; then
	echo -e "\tWarning: No epochs provided. Some tree models need epochs. Consider providing a $data/epochs.csv file."else
else
	echo -e "\tEpochs located at $data/epochs.csv"
fi

starting_trees=true
starting_tree_file="$data/starting_trees.nex"
if [ ! -f $starting_tree_file ]; then
	starting_trees=false
	echo -e "\tWarning: No starting trees provided. Some tree models need starting trees. Consider providing a $data/starting_trees.nex file."
else
	echo -e "\tStarting trees located at $data/starting_trees.nex"
fi

###############
# get modules #
###############

echo -e "\nLocating module files."
if [ ! -d $modules ]; then
	echo -e "\tProvided module directory does not exist."
	exit 1
fi

# find the tree models
if [ ! -d "$modules/tree_models" ]; then
	echo -e "\tno tree_model directory in $modules."
	exit 1
fi
tree_model_count=`ls -1 $modules/tree_models/*.Rev 2>/dev/null | wc -l`
if [ $tree_model_count = 0 ]; then
	echo -e "\tno tree models found in $modules/tree_models"
	exit 1
fi
tree_models=(`ls $modules/tree_models/*.Rev`)
echo -e "\tFound $tree_model_count tree model(s):"
for tree_model in ${tree_models[@]}; do
	echo -e "\t\t$tree_model"
done

# find the substitution models
if [ ! -d "$modules/substitution_models" ]; then
	echo -e "\tno substitution_models directory in $modules."
	exit 1
fi
substitution_models_count=`ls -1 $modules/substitution_models/*.Rev 2>/dev/null | wc -l`
if [ $substitution_models_count = 0 ]; then
	echo -e "\tno substitution models found in $modules/substitution_models"
	exit 1
fi
substitution_models=(`ls $modules/substitution_models/*.Rev`)
echo -e "\tFound $substitution_models_count substitution model(s):"
for substitution_model in ${substitution_models[@]}; do
	echo -e "\t\t$substitution_model"
done

# find the molecular-clock models
if [ ! -d "$modules/molecular_clock_models" ]; then
	echo -e "\tno molecular_clock_models directory in $modules."
	exit 1
fi
molecular_clock_models_count=`ls -1 $modules/molecular_clock_models/*.Rev 2>/dev/null | wc -l`
if [ $molecular_clock_models_count = 0 ]; then
	echo -e "\tno molecular clock models found in $modules/molecular_clock_models"
	exit 1
fi
molecular_clock_models=(`ls $modules/molecular_clock_models/*.Rev`)
echo -e "\tFound $molecular_clock_models_count molecular clock model(s):"
for molecular_clock_model in ${molecular_clock_models[@]}; do
	echo -e "\t\t$molecular_clock_model"
done

# find the morph matrix models
if [ ! -d "$modules/morph_matrix_models" ]; then
	echo -e "\tno morph_matrix_models directory in $modules."
	exit 1
fi
morph_matrix_models_count=`ls -1 $modules/morph_matrix_models/*.Rev 2>/dev/null | wc -l`
if [ $morph_matrix_models_count = 0 ]; then
	echo -e "\tno morphological matrix models found in $modules/morph_matrix_models"
	exit 1
fi
morph_matrix_models=(`ls $modules/morph_matrix_models/*.Rev`)
echo -e "\tFound $morph_matrix_models_count morphological matrix model(s):"
for morph_matrix_model in ${morph_matrix_models[@]}; do
	echo -e "\t\t$morph_matrix_model"
done

# find the morph clock models
if [ ! -d "$modules/morph_clock_models" ]; then
	echo -e "\tno morph_clock_models directory in $modules."
	exit 1
fi
morph_clock_models_count=`ls -1 $modules/morph_clock_models/*.Rev 2>/dev/null | wc -l`
if [ $morph_clock_models_count = 0 ]; then
	echo -e "\tno morphological clock models found in $modules/morph_clock_models"
	exit 1
fi
morph_clock_models=(`ls $modules/morph_clock_models/*.Rev`)
echo -e "\tFound $morph_clock_models_count morphological clock model(s):"
for morph_clock_model in ${morph_clock_models[@]}; do
	echo -e "\t\t$morph_clock_model"
done

# find the morph relative-rates models
if [ ! -d "$modules/morph_rel_rate_models" ]; then
	echo -e "\tno morph_rel_rate_models directory in $modules."
	exit 1
fi
morph_rel_rate_models_count=`ls -1 $modules/morph_rel_rate_models/*.Rev 2>/dev/null | wc -l`
if [ $morph_rel_rate_models_count = 0 ]; then
	echo -e "\tno morphological relative-rate models found in $modules/morph_rel_rate_models"
	exit 1
fi
morph_rel_rate_models=(`ls $modules/morph_rel_rate_models/*.Rev`)
echo -e "\tFound $morph_rel_rate_models_count morphological rel_rate model(s):"
for morph_rel_rate_model in ${morph_rel_rate_models[@]}; do
	echo -e "\t\t$morph_rel_rate_model"
done

# compute the total number of jobs
(( num_jobs = $tree_model_count * $substitution_models_count * $molecular_clock_models_count * $morph_matrix_models_count * $morph_rel_rate_models_count * $morph_clock_models_count * $nruns  ))

######################
# check reader files #
######################

echo -e "\nLocating data readers."
if [ ! -f "$modules/data_readers/read_molecular_data.Rev" ]; then
	echo -e "\tNo molecular data reader found."
	exit 1
fi

if [ ! -f "$modules/data_readers/read_morphological_data.Rev" ]; then
	echo -e "\tNo morphological data reader found."
	exit 1
fi

if [ ! -f "$modules/data_readers/read_taxon_data.Rev" ]; then
	echo -e "\tNo taxon data reader found."
	exit 1
fi

if [ $starting_trees ]; then
	if [ ! -f "$modules/data_readers/read_starting_trees.Rev" ]; then
		echo -e "\tNo starting-tree reader found."
		exit 1
	fi	
fi

#######################
# locate the template #
#######################

echo -e "\nLocating template file."
if [ ! -f $template ]; then
	echo -e "\tNo template found. Please provide a template file."
	exit 1
else
	echo -e "\tTemplate located at $template."
fi

###################
# locate analysis #
###################

echo -e "\nLocating analysis file."
if [ ! -f "$analysis" ]; then
	echo -e "\tNo analysis found. Please provide an analysis file."
	exit 1
else
	echo -e "\tAnalysis located at $analysis."
fi

#########################
# locate the tree moves #
#########################

echo -e "\nLocating tree moves file."
if [ ! -f "$modules/tree_moves/tree_moves.Rev" ]; then
	echo -e "\tNo tree moves found. Please provide a $modules/tree_moves/tree_moves.Rev file."
	exit 1
else
	echo -e "\tTree moves located at $modules/tree_moves/tree_moves.Rev"
fi

###################
# create the jobs #
###################

echo -e "\nCreating scripts."
echo -e "\tThere are $tree_model_count x $substitution_models_count x $molecular_clock_models_count x $morph_matrix_models_count x $morph_rel_rate_models_count x $morph_clock_models_count x $nruns = $num_jobs scripts."

if [ -d $jobdir ]; then
	rm -r $jobdir
fi
mkdir $jobdir

# create the output directory
mkdir -p $output
if [ -f "${output}/model_table.csv" ]; then
	rm ${output}/model_table.csv
fi
echo -e "number,tree_model,sub_model,mol_clock_model,morph_matrix_model,morph_rel_rate_model,morph_clock_model,run,status" >> ${output}/model_table.csv

# loop over tree models
job_id=0;
for tree_model in ${tree_models[@]}; do
	
	# loop over substitution models
	for substitution_model in ${substitution_models[@]}; do
		
		# loop over molecular clock models
		for molecular_clock_model in ${molecular_clock_models[@]}; do
			
			# loop over morphological matrix models
			for morph_matrix_model in ${morph_matrix_models[@]}; do
				
				# loop over morphological relative-rate models
				for morph_rel_rate_model in ${morph_rel_rate_models[@]}; do
					
					# loop over morphological clock models
					for morph_clock_model in ${morph_clock_models[@]}; do
						
						# loop over runs
						for run in $(seq 1 $nruns); do

							# iterate the job id
							((job_id++))
							
							# create the file
							jobfile="$jobdir/job_${job_id}.Rev"
							touch $jobfile

							# specify the data variables
							echo -e "# data files" >> $jobfile
							echo -e "mol_data_file      = \"${data}/molecular.nex\"" >> $jobfile
							echo -e "morph_data_file    = \"${data}/morpho.nex\"" >> $jobfile
							echo -e "taxa_file          = \"${data}/taxa.tsv\"" >> $jobfile
							echo -e "epoch_file         = \"${data}/epochs.csv\"" >> $jobfile
							echo -e "use_starting_tree  = ${starting_trees}" >> $jobfile
							if [ ${starting_trees} ]; then
								echo -e "starting_tree_file = \"${starting_tree_file}\"" >> $jobfile
							fi
							
							# specify data readers
							echo -e "\n# data readers" >> $jobfile
							echo -e "mol_data_reader      = \"${modules}/data_readers/read_molecular_data.Rev\"" >> $jobfile
							echo -e "morph_data_reader    = \"${modules}/data_readers/read_morphological_data.Rev\"" >> $jobfile
							echo -e "taxon_data_reader    = \"${modules}/data_readers/read_taxon_data.Rev\"" >> $jobfile
							if [ ${starting_trees} ]; then
								echo -e "starting_tree_reader = \"${modules}/data_readers/read_starting_trees.Rev\"" >> $jobfile
							fi
							
							# specify the variables
							echo -e "\n# model files" >> $jobfile
							echo -e "tree_model_file     = \"${tree_model}\"" >> $jobfile
							echo -e "sub_mode_file       = \"${substitution_model}\"" >> $jobfile
							echo -e "mol_clock_file      = \"${molecular_clock_model}\"" >> $jobfile
							echo -e "morph_matrix_file   = \"${morph_matrix_model}\"" >> $jobfile
							echo -e "morph_rel_rate_file = \"${morph_rel_rate_model}\"" >> $jobfile
							echo -e "morph_clock_file    = \"${morph_clock_model}\"" >> $jobfile
							echo -e "analysis_file       = \"${analysis}\"" >> $jobfile
							echo -e "run_ID              = ${run}" >> $jobfile
								
							# specify the tree-move file
							echo -e "\n# tree moves" >> $jobfile
							echo -e "tree_moves_file = \"${modules}/tree_moves/tree_moves.Rev\"" >> $jobfile
																		
							# name the output files		
							tree_name=${tree_model##*/}
							sub_model=${substitution_model##*/}
							mol_clock=${molecular_clock_model##*/}
							morph_mat=${morph_matrix_model##*/}
							morph_rel=${morph_rel_rate_model##*/}
							morph_clock=${morph_clock_model##*/}
							analysis_type=${analysis##*/}
							
							echo -e "\n# output path" >> $jobfile
							echo -e "output_dir = \"${output}/tree_model_${tree_name%.Rev}_sub_model_${sub_model%.Rev}_mol_clock_${mol_clock%.Rev}_morph_mat_${morph_mat%.Rev}_morph_rel_${morph_rel%.Rev}_morph_clock_${morph_clock%.Rev}/\"" >> $jobfile

							# source the template
							echo -e "\n# source the template" >> $jobfile
							echo -e "source(\"${template}\")" >> $jobfile
					
							# add the model to the model table
							echo -e "${job_id},${tree_name%.Rev},${sub_model%.Rev},${mol_clock%.Rev},${morph_mat%.Rev},${morph_rel%.Rev},${morph_clock%.Rev},${run},0" >> ${output}/model_table.csv		
							
						done
					done
				done
			done
		done
	done
done

echo -e "\nDone creating scripts."
