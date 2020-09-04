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
cat("Diagnosing LTT for MCMC from ", output_dir, " with ", ncores, " cores.\n", sep="")

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
source("src/lttPlot.R")

cols = brewer.pal(4, "Set1")

# list output files
dirs = list.dirs(output_dir, recursive=FALSE)

# loop over models
cat("Computing LTTs.\n")
x = mclapply(dirs, function(this_dir) {

    # screen log
    cat(this_dir, "\n", sep="")

    ##############
    # full trees #
    ##############

    files = list.files(this_dir, full.names=TRUE, pattern="^tree_\\d+.nex")
    samples = vector("list", length(files))
    for(i in 1:length(samples)) {

        # get the file
        this_file = files[i]

        name = gsub(".nex", "", basename(this_file))
        name = gsub("tree_", "run ", name)
        names(samples)[i] = name

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
        these_nexus_samples = lapply(these_nexus_samples, collapse.singles)
        class(these_nexus_samples) = "multiPhylo"

        # store
        samples[[i]] = these_nexus_samples

    }

    samples = samples[sapply(samples, is.null) == FALSE]

    # compute LTTs
    ltts         = lapply(samples, lttCurve, verbose=FALSE)
    combined_ltt = lttCurve(do.call(c, samples), verbose=FALSE)

    # get the name for each LTT
    names = names(samples)
    names = paste0(names, " (ESS = ", sprintf("%.3f",sapply(ltts, function(x) x$ess)), ")")

    # compute the range
    ylim = range(combined_ltt$quantiles)

    # make the plot
    ltt_dir = paste0(this_dir, "/ltt.pdf")
    pdf(ltt_dir, height=3)
    par(mar=c(2, 2, 0, 0) + 0.1, mfrow=c(1,1))
    plot(1, type="n", xlab="time", ylab="number of lineages", ylim=ylim, xlim=c(max(combined_ltt$times),0), xaxt="n", yaxt="n", log="y")
    for(i in 1:length(ltts)) {
        this_ltt = ltts[[i]]
        polygon(x = c(this_ltt$times, rev(this_ltt$times)),
                y = c(this_ltt$quantiles[1,], rev(this_ltt$quantiles[2,])),
                col = paste0(cols[i],"70"), border=NA)
    }
    for(i in 1:length(ltts)) {
        this_ltt = ltts[[i]]
        lines(this_ltt$times, this_ltt$mean, col=cols[i])
    }
    polygon(x = c(combined_ltt$times, rev(combined_ltt$times)),
            y = c(combined_ltt$quantiles[1,], rev(combined_ltt$quantiles[2,])),
            lty=2)
    lines(combined_ltt$times, combined_ltt$mean, lty=2)
    axis(1, lwd=0, lwd.tick=1)
    axis(2, lwd=0, lwd.tick=1)
    legend("bottomright", legend=names, fill=cols, bty="n", ncol=2)
    dev.off()

}, mc.cores=ncores, mc.preschedule=FALSE)


# quit
q()

