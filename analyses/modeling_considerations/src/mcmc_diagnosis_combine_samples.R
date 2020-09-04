# get the arguments
args = commandArgs(trailingOnly = TRUE)

# args = c("--output", "output")

# get the output directory
if ( "--output" %in% args ) {
    output_dir = args[which(args == "--output") + 1]
} else {
    stop("Must provide an --output argument!")
}

# get the number of cores
ncores = 1
if ( "--ncores" %in% args ) {
    ncores = as.numeric(args[which(args == "--ncores") + 1])
}

# get whether to reset
redo = "--reset" %in% args

# found the arguments
cat("Combining MCMC samples from ", output_dir, " with ", ncores, " cores.\n", sep="")

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

if ( "coda" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("coda")
}
library(coda)

if ( "parallel" %in% rownames(installed.packages()) == FALSE ) {
    install.packages("parallel")
}
library(parallel)

# list output files
dirs = list.dirs(output_dir, recursive=FALSE)

# read the model table
model_table = read.csv("output/model_table.csv", header=TRUE, stringsAsFactors=FALSE)

# get models w/o respect to run
sub_model_table = model_table[model_table$run == 1,]

# loop over models
cat("Combining samples.\n")
x = mclapply(1:nrow(sub_model_table), function(i) {

    # get this model
    this_model = sub_model_table[i,]

    # get the model settings
    this_tree_model         = this_model$tree_model
    this_sub_model          = this_model$sub_model
    this_clock_model        = this_model$mol_clock_model
    this_morph_matrix_model = this_model$morph_matrix_model
    this_morph_rel_model    = this_model$morph_rel_rate_model
    this_morph_clock_model  = this_model$morph_clock_model

    # get the runs to combine
    which_runs = model_table$tree_model == this_tree_model &
                 model_table$sub_model == this_sub_model &
                 model_table$mol_clock_model == this_clock_model &
                 model_table$morph_matrix_model == this_morph_matrix_model &
                 model_table$morph_rel_rate_model == this_morph_rel_model &
                 model_table$morph_clock_model == this_morph_clock_model

    these_runs = model_table[which_runs,]

    # remove analysis for which status != 1
    these_runs = these_runs[these_runs$status == 1,]

    # construct the directory name
    this_dir = paste0(output_dir, "/",
                      "tree_model_",   this_tree_model,
                      "_sub_model_",   this_sub_model,
                      "_mol_clock_",   this_clock_model,
                      "_morph_mat_",   this_morph_matrix_model,
                      "_morph_rel_",   this_morph_rel_model,
                      "_morph_clock_", this_morph_clock_model)

    # check that the directory exists
    if ( dir.exists(this_dir) == FALSE ) {
        warning("Could not find directory: ", this_dir)
        return()
    }

    # check if the combined samples already exist

    # check for combined trees
    if ( "tree_combined.nex" %in% list.files(this_dir) & redo == FALSE ) {
        return()
    }

    # screen log
    cat(this_dir, "\n", sep="")

    ##############
    # parameters #
    ##############

    files = paste0(this_dir, "/params_", these_runs$run, ".log")
    param_samples = vector("list", length(files))
    rows_to_drop  = vector("list", length(files))
    for(i in 1:length(param_samples)) {

        # get the file
        this_file = files[i]

        # name the tree burnin file
        tree_ess_file    = gsub(".log", ".txt", paste0(this_dir, gsub("params", "tree_ess" , gsub(this_dir, "", this_file))))
        tree_burnin_file = gsub("ess", "burnin", tree_ess_file)
        if ( file.exists(tree_ess_file) == FALSE ) {
            next
        }
        tree_ess = as.numeric(readLines(tree_ess_file))

        # name the param ess file
        param_ess_file    = gsub(".log", ".txt", paste0(this_dir, gsub("params", "ess" , gsub(this_dir, "", this_file))))
        param_burnin_file = gsub("ess", "burnin", param_ess_file)
        if ( file.exists(param_ess_file) == FALSE ) {
            next
        }
        param_ess = as.numeric(readLines(param_ess_file))

        # compute the burnin
        tree_burnin  = as.numeric(readLines(tree_burnin_file))
        param_burnin = as.numeric(readLines(param_burnin_file))
        burnin = pmax(tree_burnin, param_burnin)

        # read the parameters
        these_samples = read.table(this_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

        # decide which samples to drop
        rows_to_drop[[i]] = which(apply(these_samples, 1, function(x) any(is.nan(x)) ))
        rows_to_drop[[i]] = sort(c(1:burnin, rows_to_drop[[i]]))

        # drop NA and burnin
        these_samples = these_samples[-rows_to_drop[[i]],]
        # these_samples = these_samples[-1 * 1:burnin,]
        param_samples[[i]] = these_samples

    }
    param_samples = do.call(rbind, param_samples)
    param_samples$Iteration = 1:nrow(param_samples)
    write.table(param_samples, file=paste0(this_dir,"/params_combined.log"), row.names=FALSE, sep="\t", quote=FALSE)

    ##############
    # full trees #
    ##############

    files = paste0(this_dir, "/tree_", these_runs$run, ".nex")
    nexus_samples = vector("list", length(files))
    tree_samples  = vector("list", length(files))
    for(i in 1:length(tree_samples)) {

        # get the file
        this_file = files[i]

        # name the tree burnin file
        tree_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("tree", "tree_ess" , gsub(this_dir, "", this_file))))
        tree_burnin_file = gsub("ess", "burnin", tree_ess_file)
        if ( file.exists(tree_ess_file) == FALSE ) {
            next
        }
        tree_ess = as.numeric(readLines(tree_ess_file))

        # name the param ess file
        param_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("tree", "ess" , gsub(this_dir, "", this_file))))
        param_burnin_file = gsub("ess", "burnin", param_ess_file)
        if ( file.exists(param_ess_file) == FALSE ) {
            next
        }
        param_ess = as.numeric(readLines(param_ess_file))

        # compute the burnin
        tree_burnin  = as.numeric(readLines(tree_burnin_file))
        param_burnin = as.numeric(readLines(param_burnin_file))
        burnin = pmax(tree_burnin, param_burnin)

        # read the nexus trees
        these_nexus_samples = read.nexus(this_file)
        these_nexus_samples = these_nexus_samples[-rows_to_drop[[i]]]
        nexus_samples[[i]] = these_nexus_samples

        # read the trees
        tree_file = gsub(".nex", ".trees", this_file)
        these_tree_samples = read.table(tree_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
        these_tree_samples = these_tree_samples[-rows_to_drop[[i]],]
        tree_samples[[i]]  = these_tree_samples

    }
    nexus_samples = nexus_samples[lengths(nexus_samples) > 0]
    nexus_samples = do.call("c", nexus_samples)
    tree_samples  = do.call("rbind", tree_samples)
    tree_samples$Iteration  = 1:nrow(tree_samples)
    write.nexus(nexus_samples, file=paste0(this_dir, "/tree_combined.nex"), translate=FALSE)
    write.table(tree_samples,  file=paste0(this_dir, "/tree_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)

    ################
    # extant trees #
    ################

    files = paste0(this_dir, "/extant_tree_", these_runs$run, ".nex")
    nexus_samples = vector("list", length(files))
    tree_samples  = vector("list", length(files))
    for(i in 1:length(tree_samples)) {

        # get the file
        this_file = files[i]

        # name the tree burnin file
        tree_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("extant_tree", "tree_ess" , gsub(this_dir, "", this_file))))
        tree_burnin_file = gsub("ess", "burnin", tree_ess_file)
        if ( file.exists(tree_ess_file) == FALSE ) {
            next
        }
        tree_ess = as.numeric(readLines(tree_ess_file))

        # name the param ess file
        param_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("extant_tree", "ess" , gsub(this_dir, "", this_file))))
        param_burnin_file = gsub("ess", "burnin", param_ess_file)
        if ( file.exists(param_ess_file) == FALSE ) {
            next
        }
        param_ess = as.numeric(readLines(param_ess_file))

        # compute the burnin
        tree_burnin  = as.numeric(readLines(tree_burnin_file))
        param_burnin = as.numeric(readLines(param_burnin_file))
        burnin = pmax(tree_burnin, param_burnin)

        # read the nexus trees
        these_nexus_samples = read.nexus(this_file)
        these_nexus_samples = these_nexus_samples[-rows_to_drop[[i]]]
        nexus_samples[[i]] = these_nexus_samples

        # read the trees
        tree_file = gsub(".nex", ".trees", this_file)
        these_tree_samples = read.table(tree_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
        these_tree_samples = these_tree_samples[-rows_to_drop[[i]],]
        tree_samples[[i]]  = these_tree_samples

    }

    nexus_samples = do.call("c", nexus_samples)
    tree_samples  = do.call("rbind", tree_samples)
    tree_samples$Iteration  = 1:nrow(tree_samples)
    write.nexus(nexus_samples, file=paste0(this_dir, "/extant_tree_combined.nex"), translate=FALSE)
    write.table(tree_samples,  file=paste0(this_dir, "/extant_tree_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)

    ##############
    # phylograms #
    ##############

    files = paste0(this_dir, "/mole_phylogram_", these_runs$run, ".nex")
    nexus_samples = vector("list", length(files))
    tree_samples  = vector("list", length(files))
    for(i in 1:length(tree_samples)) {

        # get the file
        this_file = files[i]

        # name the tree burnin file
        tree_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("mole_phylogram", "tree_ess" , gsub(this_dir, "", this_file))))
        tree_burnin_file = gsub("ess", "burnin", tree_ess_file)
        if ( file.exists(tree_ess_file) == FALSE ) {
            next
        }
        tree_ess = as.numeric(readLines(tree_ess_file))

        # name the param ess file
        param_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("mole_phylogram", "ess" , gsub(this_dir, "", this_file))))
        param_burnin_file = gsub("ess", "burnin", param_ess_file)
        if ( file.exists(param_ess_file) == FALSE ) {
            next
        }
        param_ess = as.numeric(readLines(param_ess_file))

        # compute the burnin
        tree_burnin  = as.numeric(readLines(tree_burnin_file))
        param_burnin = as.numeric(readLines(param_burnin_file))
        burnin = pmax(tree_burnin, param_burnin)

        # read the nexus trees
        these_nexus_samples = read.nexus(this_file)
        these_nexus_samples = these_nexus_samples[-1 * 1:burnin]
        nexus_samples[[i]] = these_nexus_samples

        # read the nexus trees
        these_nexus_samples = read.nexus(this_file)
        these_nexus_samples = these_nexus_samples[-rows_to_drop[[i]]]
        nexus_samples[[i]] = these_nexus_samples

        # read the trees
        tree_file = gsub(".nex", ".trees", this_file)
        these_tree_samples = read.table(tree_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
        these_tree_samples = these_tree_samples[-rows_to_drop[[i]],]
        tree_samples[[i]]  = these_tree_samples

    }
    nexus_samples = do.call("c", nexus_samples)
    tree_samples  = do.call("rbind", tree_samples)
    tree_samples$Iteration  = 1:nrow(tree_samples)
    write.nexus(nexus_samples, file=paste0(this_dir, "/mole_phylogram_combined.nex"),   translate=FALSE)
    write.table(tree_samples,  file=paste0(this_dir, "/mole_phylogram_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)

    ###############
    # morphograms #
    ###############

    files = paste0(this_dir, "/morph_phylogram_", these_runs$run, ".nex")
    nexus_samples = vector("list", length(files))
    tree_samples  = vector("list", length(files))
    for(i in 1:length(tree_samples)) {

        # get the file
        this_file = files[i]

        # name the tree burnin file
        tree_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("morph_phylogram", "tree_ess" , gsub(this_dir, "", this_file))))
        tree_burnin_file = gsub("ess", "burnin", tree_ess_file)
        if ( file.exists(tree_ess_file) == FALSE ) {
            next
        }
        tree_ess = as.numeric(readLines(tree_ess_file))

        # name the param ess file
        param_ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("morph_phylogram", "ess" , gsub(this_dir, "", this_file))))
        param_burnin_file = gsub("ess", "burnin", param_ess_file)
        if ( file.exists(param_ess_file) == FALSE ) {
            next
        }
        param_ess = as.numeric(readLines(param_ess_file))

        # compute the burnin
        tree_burnin  = as.numeric(readLines(tree_burnin_file))
        param_burnin = as.numeric(readLines(param_burnin_file))
        burnin = pmax(tree_burnin, param_burnin)

        # read the nexus trees
        these_nexus_samples = read.nexus(this_file)
        these_nexus_samples = these_nexus_samples[-rows_to_drop[[i]]]
        nexus_samples[[i]] = these_nexus_samples

        # read the trees
        tree_file = gsub(".nex", ".trees", this_file)
        these_tree_samples = read.table(tree_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
        these_tree_samples = these_tree_samples[-rows_to_drop[[i]],]
        tree_samples[[i]]  = these_tree_samples

    }
    nexus_samples = do.call("c", nexus_samples)
    tree_samples  = do.call("rbind", tree_samples)
    tree_samples$Iteration  = 1:nrow(tree_samples)
    write.nexus(nexus_samples, file=paste0(this_dir, "/morph_phylogram_combined.nex"), translate=FALSE)
    write.table(tree_samples,  file=paste0(this_dir, "/morph_phylogram_combined.trees"), row.names=FALSE, sep="\t", quote=FALSE)

}, mc.cores=ncores, mc.preschedule=FALSE)


# quit
q()






















