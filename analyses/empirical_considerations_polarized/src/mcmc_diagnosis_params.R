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
cat("MCMC diagnosis on ", output_dir, " with ", ncores, " cores.\n", sep="")

# install packages
cat("Checking for packages.\n")
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

# loop over models
cat("Diagnosing parameters.\n")
x = mclapply(dirs, function(this_dir) {

    # locate parameter output files
    param_files = list.files(this_dir, full.names=TRUE, pattern="params_[0-9].log")
    param_files = param_files[grepl("filtered", param_files) == FALSE]

    # read each parameter file and do
    # MCMC diagnosis
    for(this_file in param_files) {

        ess_file    = gsub(".log", ".txt", gsub("params", "ess", this_file))
        burnin_file = gsub("ess", "burnin", ess_file)
        remove_file = gsub("ess", "remove", ess_file)

        if ( file.exists(ess_file) & redo == FALSE ) {
            next
        }

        # screen log
        cat(gsub(paste0(output_dir,"/"), "", this_file), "\n", sep="")

        # read the samples
        samples = read.table(this_file, header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)

        # discard the iteration
        iterations = samples[,1]
        samples = samples[,-1]

        # discard branch rates
        samples = samples[,grepl("^branch_rate_mole", colnames(samples)) == FALSE & grepl("^branch_rate_morph", colnames(samples)) == FALSE]

        # discard mixture model parameters (non-identifiable)
        samples = samples[,grepl("^morph_matrix_weights", colnames(samples)) == FALSE]
        samples = samples[,grepl("^pi_morpho", colnames(samples)) == FALSE]

        samples = samples[,grepl("^speciation_mixture_rates", colnames(samples)) == FALSE]
        samples = samples[,grepl("^speciation_mixture_weights", colnames(samples)) == FALSE]

        samples = samples[,grepl("^extinction_mixture_rates", colnames(samples)) == FALSE]
        samples = samples[,grepl("^extinction_mixture_weights", colnames(samples)) == FALSE]

        samples = samples[,grepl("^fossilization_mixture_rates", colnames(samples)) == FALSE]
        samples = samples[,grepl("^fossilization_mixture_weights", colnames(samples)) == FALSE]

        # discard tree parameters
        samples = samples[,grepl("^t\\[", colnames(samples)) == FALSE]
        samples = samples[,grepl("stem_age", colnames(samples)) == FALSE]
        samples = samples[,grepl("root_age", colnames(samples)) == FALSE]
        samples = samples[,grepl("mole_tree_length", colnames(samples)) == FALSE]
        samples = samples[,grepl("morph_tree_length", colnames(samples)) == FALSE]
        samples = samples[,grepl("mean_branch_rate_mole", colnames(samples)) == FALSE]
        samples = samples[,grepl("mean_branch_rate_morph", colnames(samples)) == FALSE]

        # discard rows with NaN
        rows_to_remove = which(colSums(apply(samples, 1, is.nan)) > 0)
        if ( length(rows_to_remove) > 0 ) {
            iterations = iterations[-rows_to_remove]
            samples    = samples[-rows_to_remove,]
        }

        # compute the burnin
        num_samples = nrow(samples)
        candidate_burnins = floor(seq(0.1, 0.95, 0.01) * num_samples)
        ess_per_burnin = numeric(length(candidate_burnins))
        # cat("Computing ESS for various burnin values.\n")
        # bar = txtProgressBar(style=3, width=40)
        for(i in 1:length(candidate_burnins)) {

            # get the burnin fraction
            this_burnin = candidate_burnins[i]

            # compute the ESS for each parameter
            these_ess   = as.numeric(effectiveSize(samples[-1 * 1:this_burnin,]))

            # compute the harmonic mean ESS
            ess_per_burnin[i] = 1 / mean(1 / these_ess)

            # setTxtProgressBar(bar, i / length(candidate_burnins))

        }
        # cat("\n")

        pdf_file = gsub(".txt", ".pdf",  gsub("^tree", "^burnin", burnin_file))
        pdf(pdf_file, height=3)
        par(mar=c(2,2,0,0)+0.1, mfrow=c(1,2))
        plot(candidate_burnins, ess_per_burnin, type="l")

        # refine the burnin
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
            these_ess   = as.numeric(effectiveSize(samples[-1 * 1:this_burnin,]))

            # compute the harmonic mean ESS
            ess_per_burnin[i] = 1 / mean(1 / these_ess)

            # setTxtProgressBar(bar, i / length(candidate_burnins))

        }
        # cat("\n")

        # get the final burnin
        burnin = candidate_burnins[which.max(ess_per_burnin)]
        ess    = ess_per_burnin[which.max(ess_per_burnin)]
        iteration_burnin = iterations[burnin]

        plot(candidate_burnins, ess_per_burnin, type="l")
        abline(v=burnin, lty=2)
        dev.off()

        # write new samples
        new_file = gsub("params_", "filtered_params_", this_file)
        filtered_samples = read.table(this_file, header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
        if ( length(rows_to_remove) > 0 ) {
            filtered_samples = filtered_samples[-rows_to_remove,]
        }
        filtered_samples$Iteration = 1:nrow(filtered_samples)
        write.table(filtered_samples, file=new_file, quote=FALSE, sep="\t", row.names=FALSE)

        # record the ESS
        cat(ess, "\n", sep="", file=ess_file)
        cat(burnin, "\n", sep="", file=burnin_file)
        cat(iterations[rows_to_remove], "\n", sep="\t", file=remove_file)
        # cat(this_file, "\n", burnin, "\t", ess, "\n", sep="")

    }

}, mc.cores=ncores, mc.preschedule=FALSE)

# quit
q()
