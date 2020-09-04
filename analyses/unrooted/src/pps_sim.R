# this script will simulate one dataset per posterior sample,
# compute summary statistics based on the simulation and the observed data,
# and create a csv file that includes one column per summary statistic.
# each row corresponds to a simulated dataset

# get the arguments
args = commandArgs(trailingOnly = TRUE)

# args = c("--data","data",
#          "--samples","output/tree_model_variable_rate_fbd_mixture_all_outgroup_sub_model_GTR_I_G_mol_clock_UCLN_morph_mat_JC_G_morph_rel_unlinked_morph_clock_UCLN",
#          "--output","output_pps/tree_model_variable_rate_fbd_mixture_all_outgroup_sub_model_GTR_I_G_mol_clock_UCLN_morph_mat_JC_G_morph_rel_unlinked_morph_clock_UCLN",
#          "--overwrite","true",
#          "--ncores","1")

# get the data directory
if ( "--data" %in% args ) {
    data_dir = args[which(args == "--data") + 1]
} else {
    stop("Must provide an --data argument!")
}

# get the sample directory
if ( "--samples" %in% args ) {
    sample_dir = args[which(args == "--samples") + 1]
} else {
    stop("Must provide an --sample argument!")
}

# get the output directory
if ( "--output" %in% args ) {
    output_dir = args[which(args == "--output") + 1]
} else {
    stop("Must provide an --output argument!")
}

# the number of simulations
nsim = 1000
if ( "--nsim" %in% args ) {
    nsim = as.numeric(args[which(args == "--nsim") + 1])
}

# get the overwrite
overwrite = FALSE
if ( "--overwrite" %in% args ) {
    if ( tolower(args[which(args == "--overwrite") + 1]) == "true" ) {
        overwrite = TRUE
    } else if ( tolower(args[which(args == "--overwrite") + 1]) == "false" ) {
        overwrite = FALSE
    } else {
        stop("Invalided --overwrite value!")
    }
}

cat("Posterior predictive simulation for:\n  ", sample_dir,"\n", sep="")

# cat(overwrite, "\n")
# q()

# create filenames
data_fn   = "data/morpho.nex"
sample_fn = paste0(sample_dir,"/params_combined.log")
tree_fn   = paste0(sample_dir,"/morph_phylogram_combined.nex")
output_fn = paste0(output_dir,"/pps.csv")

# check that files exist
if ( file.exists(data_fn) == FALSE ) {
    cat("Data file could not be located.\n")
    q()
}

if ( file.exists(sample_fn) == FALSE ) {
    cat("Parameter file could not be located.\n")
    q()
}

if ( file.exists(tree_fn) == FALSE ) {
    cat("Tree file could not be located.\n")
    q()
}

if ( file.exists(output_fn) == TRUE & overwrite == FALSE ) {
    cat("Already complete.\n")
    q()
}

# install packages
cat("Checking for packages.\n")
if ( "ape" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("ape")
}
library(ape)

if ( "phytools" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("phytools")
}
library(phytools)

if ( "phangorn" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("phangorn")
}
library(phangorn)

if ( "rncl" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("rncl")
}
library(rncl)

if ( "stringr" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("stringr")
}

if ( "parallel" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("parallel")
}
library(parallel)

# source the scripts
source("src/pps_functions.R")

# read the observed data
obs_data  = readParsedData(data_fn)

# compute the number of partitions
num_partitions = length(obs_data)
num_states     = as.numeric(names(obs_data))

# make sure data is properly formatted
for(i in 1:num_partitions) {

    # get this partition
    this_part = obs_data[[i]]

    # replace - with ?
    this_part[this_part == "-"] = "?"

    # remove invariant characters
    is_invariant = apply(this_part, 2, function(x) {
        obs_states = unique(x)
        sum(obs_states != "?") == 1
    })
    this_part = this_part[,is_invariant == FALSE]

    # make sure it's a matrix
    if ( class(this_part) == "character" ) {
        this_part = as.matrix(this_part)
    }

    # replace the data
    obs_data[[i]] = this_part

}

# read the samples
cat("Reading parameters.\n")
param_samples = read.table(sample_fn, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")

cat("Reading trees.\n")
tree_samples  = read.nexus(tree_fn)
tree_samples  = lapply(tree_samples, make_tree_bifurcating)

# do the simulation
num_samples = nsim

# can't parallelize this because it takes too much memory
cat("Simulating.\n")
pps_stats = vector("list", num_samples)
# bar = txtProgressBar(style=1, width=40)
for(i in 1:num_samples) {

    # get this sample
    this_sample = sample.int(nrow(param_samples), size=1)

    # get this param
    these_params = param_samples[this_sample,]

    # get this tree
    this_tree = tree_samples[[this_sample]]

    cat(i, sep="")

    # simulate a dataset
    these_sims = vector("list", num_partitions)
    for(j in 1:num_partitions) {
        cat("\t", j, sep="")
        these_sims[[j]] = simulate_pps_dataset(obs_data[[j]], this_tree, these_params, num_states[j], j)
    }
    cat("\n")

    # compute statistics
    pps_stats[[i]] = pps_statistics(obs_data, these_sims, this_tree, num_states)

    # setTxtProgressBar(bar, i / num_samples)

}
cat("\n")
pps_stats = do.call(rbind, pps_stats)

# write the statistics to file
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

pdf(paste0(output_dir,"/pps.pdf"), height=4)
par(mar=c(6,3,0,0)+0.1, las=2)
boxplot(pps_stats, outline=FALSE)
abline(h=0, lty=2)
dev.off()

write.table(pps_stats, file=output_fn, quote=FALSE, row.names=FALSE, sep=",")

# exit
q()










