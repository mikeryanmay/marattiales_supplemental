# Supplementary Archive

Welcome to the supplementary archive for "Inferring the total-evidence timescale of marattialean fern evolution in the face of extreme model sensitivity and the absence of close extant relatives".
This archive is divided into three subdirectories: data, analyses, and scripts, each described below.

The code presented here is suitable for use with `RevBayes` compiled from the `ssbdp_fixes` branch, commit `ff39b2f1cf5604870b185ed0b85fd074cea517da`.

## Data

We provide our raw morphological and molecular datasets in `data/raw`, and corresponding curated datasets in `data/curated`.
For morphological data, we provide character data for all (including the ancient plant taxa we used in one set of analyses).
We used subsets of these morphological and molecular datasets for most of our analyses; we also provide these subsets in the analysis-specific directories (described below).

## Analyses

This directory contains five subdirectories, corresponding to the primary analyses performed in our study.
Each subdirectory contains four folders: `/data`, `/modules`, `/scripts`, and `/src`, described next.
**The scripts folder is the most important**: the scripts in this folder can be used to generate `RevBayes` analyses using the files in the `data` and the code in `modules`.

### `data`

The directory `/analyses/\*/data` contains data subsets specific to the analysis:
- `molecular.nex`: molecular alignments.
- `morpho.nex`: morphological alignments.
- `taxa.tsv`: age data (min and max age) for each taxon.
- `epochs.csv`: age information for epochs used by the episodic fossilized birth-death models.
- `char_table.tsv`: a table of morphological metadata (number of states, state labels, etc) used mainly for plotting.
- `starting_trees.nex`: a set of starting trees needed to initialize the epochal fossilized birth-death analyses. We generated these trees under a uniform time-tree model.

### `modules`

The `/modules` directory contains code chunks that correspond to the different components of the total-evidence model, organized by component, as well as code chunks for performing the MCMC analyses.
These are mostly used by the script `make_analyses.sh` (described in the next section) to assemble functional `RevBayes` scripts.

### `scripts`

This directory contains bash scripts that can be used to create the `RevBayes` code we used in our study.
In particular, `make_analyses.sh` generates a directory containing RevBayes scripts that are built from the code in modules.
The remaining scripts perform MCMC diagnosis, post-processing, posterior-predictive simulation, and stochastic mapping.
We describe each script in the order they should be used.
All of these scripts are intended to be run from the corresponding /analyses directory

##### Using `make_analyses.sh`

We assume the TED model has five components: the substitution model, the molecular clock model, the morphological transition model, the morphological clock model, and the tree model.
There is one directory in modules for each of these components.
Running `make_analyses.sh` from the /analyses/* directory will assemble one script per combination of model components.
This script has arguments (and defaults when applicable):
- `-d`: the directory to locate data (default /data)
- `-m`: the directory to locate module files (default /modules)
- `-a`: the type of analysis to perform (e.g., MCMC or MCMCMC), (default /modules/analysis/MCMC.Rev)
- `-t`: the RevBayes template used to stitch modules together (default /modules/templates/template.Rev)
- `-o`: the directory in which to place output (default /output)
- `-n`: the number of replicate analyses to perform (default 4)
- `-j`: the directory to place the RevBayes scripts (default /jobs)

The default values largely correspond to those we used in our analyses, so that executing the following code with populate the jobs folder with scripts necessary to recreate our analyses:

	bash scripts/make_analyses.sh -a modules/analysis/MCMCMC.Rev -t modules/templates/template.Rev

The resulting RevBayes scripts are intended to be run from the /analyses/* directory.

`make_analyses.sh` is agnostic to the contents of the model component directories in /modules, as long as there is at least one model per module.
Additional models can be added to this workflow by placing them in the correct model component within /modules.

##### Using `mcmc_diagnosis.sh`

This script reads output files and checks MCMC behavior for continuous and tree parameters as described in our supplemental material. It is invoked with:

	bash scripts/mcmc_diagnosis.sh -o output -n 4 -f 0

The argument ncores (`-n`) can be any suitable value, and the overwrite argument (`-f`) determines whether output that has already been diagnosed is re-diagnosed (default false).

This script creates summaries of average burnin over time for numerical parameters and tree distance metrics, as well as MDS plots among runs in the output directory for the set of analyses.
Users should inspect these summaries to determine whether each run behaved adequately, and that the runs converged to the same joint posterior distribution.

##### Using `combine_samples.sh`

After performing analyses and mcmc diagnosis, this script can be used to combine (post-burnin) samples.
The `make_scripts.sh` script creates a table, /output/model_table.csv.
The last column in this table, "status" (default 0), indicates whether the corresponding analysis should be included in the combined sample.
After performing MCMC diagnosis (described above), users can set the "status" of successful runs to 1.

Invoking the command:

	bash scripts/combine_samples.sh -o output -n 4

will concatenate post-burnin samples from successful runs (determined by the status value in the model table) into combined log files for numerical parameters and tree samples.
The argument -n determines the number of cores.

**Important**: After combining samples, the user will need to perform an additional analysis to create summary trees.
If desired, the user should use:

	bash scripts/make_analyses.sh -a modules/analysis/summarize_combined.Rev -j jobs_summarize

and run the resulting scripts to create summary trees (which can be time consuming).
The resulting MCC trees are also needed for stochastic mapping, described below.

##### Using `pps_simulations.sh`

This script uses combined samples (create in the previous step) to simulate new morphological datasets.
We then compute the summary statistics described in the main text for each simulated dataset, and write them into a .csv for post-processing.

This script has arguments (and defaults when applicable):
- `-d`: the directory to locate data (default /data)
- `-s`: the directory with the MCMC output (default /output)
- `-o`: the directory in which to place the PPS statistics (default /output)
- `-n`: the number of cores to use (default 1)
- `-i`: the number of posterior-predictive simulation replicates (default 1000)
- `-r`: whether to overwrite simulations (0) or to skip analyses that have already been simulated (1) (default 0)

Invoking the command:
	
	bash scripts/pps_simulations.sh

will create a pps.csv in each output file that contains one row per simulated dataset, with the first two columns corresponding to the S and V statistics, respectively.
(The remaining columns correspond to the statistic calculated for each morphological partition and are reported for completeness.)

##### Using `stochastic_mapping.sh`

This script performs stochastic mapping with the same arguments as `pps_simulations.sh`:

	bash scripts/stochastic_mapping.sh -n 4

### `src`

This directory contains code used by the bash scripts in the scripts directory.

## Scripts

This directory contains miscellaneous scripts that either do not belong to individual analyses, or that may be useful outside of the individual analyses.
- `geoscale_axis.R`: for plotting geological timescales
- `mds.R`: for performing MDS analysis with posterior samples of trees. Examples of use can be found in `analyses/src/mcmc_diagnosis_mds.R`
- `plot_tree.R`: for plotting nice serial-sampled trees with sampled ancestors.
- `simulate_epochs.R` and `simulate_epochs.cpp`: for simulating implied diversity over time.
