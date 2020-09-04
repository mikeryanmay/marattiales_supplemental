# get the arguments
args = commandArgs(trailingOnly = TRUE)

# args = c("--output", "output/full_test")

# get the output directory
if ( "--output" %in% args ) {
    output_dir = args[which(args == "--output") + 1]
} else {
    stop("Must provide an --output argument!")
}

# found the arguments
cat("Compiling MCMC diagnostic statistics for ", output_dir, ".\n", sep="")

# list output files
dirs = list.dirs(output_dir, recursive=FALSE)

# loop over models
ess_df = vector("list", length(dirs))
i = 0
for(this_dir in dirs) {

    # check for files
    tree_ess_files = list.files(this_dir, full.names=TRUE, pattern="^tree_ess_[0-9].txt")

    # create containers
    models        = c()
    runs          = c()
    tree_ess_vals = c()
    ess_vals      = c()

    # first check if tree_ess is done
    for(tree_ess_file in tree_ess_files) {

        # check for the ess file
        ess_file = gsub("tree_ess", "ess", tree_ess_file)
        if ( file.exists(ess_file) == TRUE ) {

            # read the files
            tree_ess = as.numeric(readLines(tree_ess_file))
            ess      = as.numeric(readLines(ess_file))

            # get the model name
            model_name = tail(strsplit(this_dir, "/")[[1]], 1)

            # get the run number
            run_ID = as.numeric(gsub(".txt", "", strsplit(tree_ess_file, "tree_ess_")[[1]][2]))

            # create the identifier
            this_fn = paste0(model_name, "/run_", run_ID)

            # store
            models        = c(models, model_name)
            runs          = c(runs, run_ID)
            tree_ess_vals = c(tree_ess_vals, tree_ess)
            ess_vals      = c(ess_vals, ess)

        }

    } # end loop over runs

    # create a data frame
    i = i + 1
    ess_df[[i]] = data.frame(model=models, run=runs, ess=ess_vals, tree_ess=tree_ess_vals)


} # end loop over models
ess_df = do.call(rbind, ess_df)

# write the table
write.table(ess_df, file=paste0(output_dir,"/ess_summary.csv"), quote=FALSE, row.names=FALSE, sep=",")

# quit
q()