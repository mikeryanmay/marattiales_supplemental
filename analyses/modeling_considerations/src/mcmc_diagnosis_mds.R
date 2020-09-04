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

library(RColorBrewer)
source("src/tsv.R")
source("src/pps_scripts.R")
source("src/lttPlot.R")

cols = brewer.pal(4, "Set1")

# list output files
dirs = list.dirs(output_dir, recursive=FALSE)

# loop over models
cat("Computing MDS plots.\n")
x = mclapply(dirs, function(this_dir) {

    # screen log
    cat(this_dir, "\n", sep="")

    ##############
    # full trees #
    ##############

    files     = list.files(this_dir, full.names=TRUE, pattern="^tree_\\d+.nex")
    names     = numeric(length(files))
    samples   = vector("list", length(files))
    map_trees = vector("list", length(files))
    mcc_trees = vector("list", length(files))
    mrc_trees = vector("list", length(files))
    for(i in 1:length(samples)) {

        # get the file
        this_file = files[i]

        name = gsub(".nex", "", basename(this_file))
        name = gsub("tree_", "run ", name)
        names[i] = name

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
        these_nexus_samples = these_nexus_samples[-1 * 1:burnin]

        # make the trees bifurcating
        these_nexus_samples = lapply(these_nexus_samples, make_tree_bifurcating)
        class(these_nexus_samples) = "multiPhylo"

        # store
        samples[[i]] = these_nexus_samples

        # read the consensus trees
        this_consensus_tree = list(make_tree_bifurcating(multi2di(read.nexus(gsub(".nex", "_MAP.tre", this_file)))))
        class(this_consensus_tree) = "multiPhylo"
        map_trees[[i]] = this_consensus_tree

        this_consensus_tree = list(make_tree_bifurcating(multi2di(read.nexus(gsub(".nex", "_MCC.tre", this_file)))))
        class(this_consensus_tree) = "multiPhylo"
        mcc_trees[[i]] = this_consensus_tree

        this_consensus_tree = list(make_tree_bifurcating(multi2di(read.nexus(gsub(".nex", "_MRC.tre", this_file)))))
        class(this_consensus_tree) = "multiPhylo"
        mrc_trees[[i]] = this_consensus_tree

    }

    # drop samples without burnin
    names = names[sapply(samples, is.null) == FALSE]
    samples = samples[sapply(samples, is.null) == FALSE]
    map_trees = map_trees[sapply(map_trees, is.null) == FALSE]
    mcc_trees = mcc_trees[sapply(mcc_trees, is.null) == FALSE]
    mrc_trees = mrc_trees[sapply(mrc_trees, is.null) == FALSE]

    # combine the lists
    combined_samples = c(samples, map_trees, mcc_trees, mrc_trees)

    # define the colors
    colors = rep(brewer.pal(length(samples), "Set1") , times=4)
    pch    = rep(c(19, 3, 4, 8), each=length(samples))
    cex    = rep(1.0, length(pch))
    cex[1:length(samples)] = 0.2

    # do MDS
    this_mds = treeSetVis(combined_samples, num_samples=400)

    mds_dir = paste0(this_dir, "/mds.pdf")
    pdf(mds_dir, height=3)
    par(mar=c(0,0,0,0)+0.1, mfrow=c(1,2))
    this_mds$plot(type="KF", bty="n", xaxt="n", yaxt="n", plot_lines=FALSE, cex=cex, pch=pch, colors=colors)
    legend("bottomright", legend=names, pch=19, col=colors[1:4], bg="white", ncol=2)
    this_mds$plot(type="RF", bty="n", xaxt="n", yaxt="n", plot_lines=FALSE, cex=cex, pch=pch, colors=colors)
    dev.off()

}, mc.cores=ncores, mc.preschedule=FALSE)


# quit
q()

