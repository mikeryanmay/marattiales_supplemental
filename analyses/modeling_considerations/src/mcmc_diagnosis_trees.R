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
cat("MCMC tree diagnosis on ", output_dir, " with ", ncores, " cores.\n", sep="")

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

make_tree_bifurcating = function(tree) {

    # check if there are sampled ancestors
    if ( is.binary(tree) == TRUE ) {
        return(tree)
    }

    # find the sampled ancestors
    sampled_ancestors = tree$node.label[tree$node.label != ""]
    num_sampled_ancestors = length(sampled_ancestors)

    # add any sampled ancestors as zero-length branches
    tmp_tree = tree
    for(i in 1:num_sampled_ancestors) {

        # get this sampled ancestor
        this_sampled_ancestor = sampled_ancestors[i]

        # find the node to attach it to
        attachment_point = which(tmp_tree$node.label == this_sampled_ancestor) + length(tmp_tree$tip.label)

        # attach the tip
        tmp_tree = bind.tip(tmp_tree, this_sampled_ancestor, edge.length=0.0, where=attachment_point)

    }

    return(collapse.singles(tmp_tree))

}

# list output files
dirs = list.dirs(output_dir, recursive=FALSE)

# loop over models
cat("Diagnosing trees.\n")
x = mclapply(dirs, function(this_dir) {

    # locate parameter output files
    param_files = list.files(this_dir, full.names=TRUE, pattern="^tree_[0-9].nex")

    # read each parameter file and do
    # MCMC diagnosis
    for(this_file in param_files) {

        ess_file    = gsub(".nex", ".txt", paste0(this_dir, gsub("tree", "tree_ess" , gsub(this_dir, "", this_file))))
        burnin_file = gsub("ess", "burnin", ess_file)

        if ( file.exists(ess_file) & redo == FALSE ) {
            next
        }

        cat(this_file, "\n", sep="")

        # read the samples
        samples = read.nexus(this_file)

        # make the trees bifurcating
        samples = lapply(samples, make_tree_bifurcating)
        class(samples) = "multiPhylo"

        # read the summary trees
        summary_tree_files = c(gsub(".nex", "_MCC.tre", this_file), gsub(".nex", "_MAP.tre", this_file))
        # summary_tree_files = c(gsub(".nex", "_MCC.tre", this_file), gsub(".nex", "_MAP.tre", this_file))
        focal_trees = lapply(summary_tree_files, function(x) {
            make_tree_bifurcating( read.nexus(x) )
        })

        # compute the distance to the summary trees
        distance_samples = do.call(cbind, lapply(focal_trees, function(x) KF.dist(samples, x, rooted=TRUE) ))

        # compute the burnin
        num_samples = length(samples)
        candidate_burnins = floor(seq(0.01, 0.90, 0.01) * num_samples)
        ess_per_burnin = numeric(length(candidate_burnins))
        # cat("Computing ESS for various burnin values.\n")
        # bar = txtProgressBar(style=3, width=40)
        for(i in 1:length(candidate_burnins)) {

            # get the burnin fraction
            this_burnin = candidate_burnins[i]

            # compute the ESS for each parameter
            these_ess   = as.numeric(effectiveSize(distance_samples[-1 * 1:this_burnin,]))

            # compute the harmonic mean ESS
            ess_per_burnin[i] = 1 / mean(1 / these_ess)

            # setTxtProgressBar(bar, i / length(candidate_burnins))

        }
        # cat("\n")

        pdf_file = gsub(".txt", ".pdf", burnin_file)
        pdf(pdf_file, height=3)
        par(mar=c(2,2,0,0)+0.1, mfrow=c(1,3))
        matplot(distance_samples, type="l", lty=1, xlim=c(0,length(samples)))
        plot(candidate_burnins, ess_per_burnin, type="l")

        max_candidate = max(which(candidate_burnins < length(samples) * 0.5))

        # max_ess = max(ess_per_burnin[1:max_candidate])
        max_ess = max(ess_per_burnin)
        which_max_ess = which(ess_per_burnin == max_ess)
        if (which_max_ess == 1) {
            min_burnin = candidate_burnins[which_max_ess]
            max_burnin = candidate_burnins[which_max_ess + 1]
        } else if (which_max_ess == length(candidate_burnins) ) {
            min_burnin = candidate_burnins[which_max_ess]
            max_burnin = nrow(samples) - 1
        } else {
            min_burnin = candidate_burnins[which_max_ess - 1]
            max_burnin = candidate_burnins[which_max_ess + 1]
        }

        candidate_burnins = min_burnin : max_burnin
        ess_per_burnin = numeric(length(candidate_burnins))
        # cat("Refining burnin.\n")
        # bar = txtProgressBar(style=3, width=40)
        for(i in 1:length(candidate_burnins)) {

            # get the burnin fraction
            this_burnin = candidate_burnins[i]

            # compute the ESS for each parameter
            these_ess   = as.numeric(effectiveSize(distance_samples[-1 * 1:this_burnin,]))

            # compute the harmonic mean ESS
            ess_per_burnin[i] = 1 / mean(1 / these_ess)

            # setTxtProgressBar(bar, i / length(candidate_burnins))

        }
        # cat("\n")

        # get the final burnin
        sample_index = which.max(ess_per_burnin)
        burnin = candidate_burnins[sample_index]
        ess    = ess_per_burnin[sample_index]

        plot(candidate_burnins, ess_per_burnin, type="l")
        abline(v=burnin, lty=2)
        dev.off()

        # record the ESS
        cat(ess, "\n", sep="", file=ess_file)
        cat(burnin, "\n", sep="", file=burnin_file)
        # cat(burnin, "\t", ess, "\n", sep="")

    }

}, mc.cores=ncores, mc.preschedule=FALSE)

# quit
q()
