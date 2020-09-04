library(ape)
library(phangorn)
library(phytools)

read_pps_summary = function(sim_dir) {

    # read the file
    stats = read.csv(sim_dir, header=FALSE)

    # get the statistics
    scores = stats[,2]
    names(scores) = stats[,1]

    return(scores)

}

read_sim_data = function(filename) {

    # read the lines
    lines = readLines(filename)

    # drop any comments
    lines = lines[grepl("\\[", lines) == FALSE]

    # find all the blocks
    start_lines = which(toupper(lines) == "MATRIX") + 1
    end_lines   = which(toupper(lines) == "END;") - 2

    # for each block, extract the data
    alns = vector("list", length(start_lines))
    for(i in 1:length(start_lines)) {

        # get the lines for this block
        these_lines = lines[start_lines[i]:end_lines[i]]

        # split the taxa and the data
        these_lines = strsplit(these_lines, "\t")

        # get the taxa
        taxa = gsub(" ","", sapply(these_lines, function(x) x[[1]]))

        # get the data
        data = do.call(rbind, lapply(these_lines, function(x) {
            tmp = gsub(" ","",x[[2]])
            tmp = gsub("\\(01\\)","?",tmp)
            return(strsplit(tmp,"")[[1]])
        }))

        rownames(data) = taxa

        alns[[i]] = data

    }

    return(alns)

}

pps_summary = function(obs_data, sim_dir) {

    missing = colMeans(obs_data == "?")

    # read the tree
    tree = read.tree(paste0(sim_dir,"/tree.tre"))
    taxa = tree$tip.label

    # read the simulated data
    sim_data = read_sim_data(paste0(sim_dir,"/sim_char.nex"))[[1]]

    # drop missing taxa from both alignments
    obs_data = obs_data[taxa,]
    sim_data = sim_data[taxa,]

    # convert the data to phangorn format
    obs_data = phyDat(obs_data, type="USER", levels=c("0","1"))
    sim_data = phyDat(sim_data, type="USER", levels=c("0","1"))

    # calculate the statistics
    obs_pars = parsimony(tree, obs_data, site="site")
    sim_pars = parsimony(tree, sim_data, site="site")

    # drop characters with score 0
    sim_pars = sim_pars[obs_pars != 0]
    obs_pars = obs_pars[obs_pars != 0]

    # the mean parsimony scores
    sum_obs_pars = sum(obs_pars)
    sum_sim_pars = sum(sim_pars)

    # the variance in parsimony scores
    var_obs_pars = var(obs_pars)
    var_sim_pars = var(sim_pars)

    # compute the statistics
    sum_stat = sum_sim_pars - sum_obs_pars
    var_stat = var_sim_pars - var_obs_pars

    # TODO: more stats here

    # record the statistics
    stat_fn = paste0(sim_dir, "/stats.csv")
    cat("sum_stat,", sum_stat, "\n", sep="", file=stat_fn)
    cat("var_stat,", var_stat, "\n", sep="", file=stat_fn, append=TRUE)

    return(invisible())

    # plot(x=missing, y=sim_pars - obs_pars)
    # plot(x=missing, y=sim_pars)

    # plot per character
    dummy_tree = tree
    # dummy_tree$edge.length[] = 1
    num_chars = attr(obs_data, "nr")
    obs_anc_states = ancestral.pars(tree, obs_data, type="ACCTRAN")
    sim_anc_states = ancestral.pars(tree, sim_data, type="ACCTRAN")
    pdf("test/dump.pdf", height=10, width=10)
    for(i in 1:num_chars) {
        par(mfrow=c(1,2), oma=c(0,0,2,0))
        plotAnc(dummy_tree, obs_anc_states, i=i, show.tip.label=TRUE, cex=0.5, cex.pie=0.4, pos=NULL, no.margin=TRUE)
        mtext(paste0("N = ", obs_pars[i]), outer=FALSE)
        plotAnc(dummy_tree, sim_anc_states, i=i, show.tip.label=TRUE, cex=0.5, cex.pie=0.4, pos=NULL, no.margin=TRUE)
        mtext(paste0("N = ", sim_pars[i]), outer=FALSE)
        mtext(paste0("character ", i), outer=TRUE)
    }
    dev.off()

}

simulate_pps_model = function(
    obs_data,
    taxon_data,
    output_dir,
    sim_dir,
    run_ID,
    nsim           = 1000,
    num_gamma_cats = 4,
    num_mix_cats   = 5,
    num_states     = 2,
    part           = 1,
    verbose        = TRUE
) {

    # read the observed data
    if (verbose) cat("Reading observed data.\n")
    obs_data = readParsedData(obs_data)[[part]]
    taxa = read.table(taxon_data, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t")
    missing_taxa = taxa$taxon[taxa$taxon %in% rownames(obs_data) == FALSE]

    # read the trees
    if (verbose) cat("Reading trees.\n")

    # prefer nexus files, if they exist
    tree_fn = paste0(output_dir,"/morph_phylogram_", run_ID, ".nex")
    if (file.exists(tree_fn)) {
        trees = read.nexus(tree_fn)
    } else {
        tree_fn = paste0(output_dir,"/morph_phylogram_", run_ID, ".trees")
        trees   = read.tree(tree_fn)
    }
    num_trees = length(trees)

    # read the samples
    if (verbose) cat("Reading samples.\n")
    param_fn = paste0(output_dir,"/params_",run_ID,".log")
    samples  = read.table(param_fn, header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE)
    num_samples = nrow(samples)

    # determine the samples to use
    max_sample    = pmin(num_trees, num_samples, na.rm=TRUE)
    these_samples = sample.int(max_sample, size=nsim, replace=TRUE)

    # create the directory
    pps_dir = paste0(sim_dir,"/pps_sim_",run_ID)
    dir.create(pps_dir, recursive=TRUE, showWarnings=FALSE)

    # do the simulation
    if (verbose) cat("Simulating data.\n")
    if (verbose) bar = txtProgressBar(style=3, width=40)
    for(i in 1:nsim) {

        # get this sample
        this_sample = these_samples[i]

        # get the relevant things
        this_tree   = trees[[this_sample]]
        this_params = samples[this_sample,]

        # make the tree bifurcating, drop missing taxa, etc.
        this_tree = make_tree_bifurcating(this_tree)
        this_tree = drop.tip(this_tree, missing_taxa)

        # simulate the dataset
        this_sim = list(simulate_pps_dataset(obs_data, this_tree, this_params, num_gamma_cats=num_gamma_cats))
        names(this_sim) = "1"

        # write the simulation (including the tree!)
        this_sim_dir = paste0(pps_dir,"/sim_",i)
        dir.create(this_sim_dir, showWarnings=FALSE)

        # write the tree
        write.tree(this_tree, file=paste0(this_sim_dir,"/tree.tre"))

        # write the alignment
        writeData(this_sim, paste0(this_sim_dir,"/sim_char.nex"))

        # write metadata
        cat(this_sample, sep="", file=paste0(this_sim_dir,"/metadata.txt"))
        cat("\n", sep="", file=paste0(this_sim_dir,"/metadata.txt"), append=TRUE)

        if (verbose) setTxtProgressBar(bar, i / nsim)

    } # end loop over simulations

}

simulate_pps_dataset = function(obs, tree, params, part=1, num_gamma_cats=4, num_mix_cats=5, num_states = 2) {

    recover()

    # get some statistics
    taxa  = tree$tip.label
    obs   = obs[taxa,]
    nchar = ncol(obs)
    ntax  = nrow(obs)
    obs[obs == "-"] = "?"

    # determine if this model has among-character rate variation
    has_morph_gamma = "gamma_alpha_morph[1]" %in% colnames(params)
    if ( has_morph_gamma == TRUE ) {

        # get the parameter
        alpha = params[[paste0("gamma_alpha_morph[",part,"]")]]

        # get the rate categories
        discrete_rates = discrete.gamma(alpha, num_gamma_cats)

    } else {

        discrete_rates = 1.0

    }

    # determine the morphology partition rate-multiplier (if present)
    has_prop_rates = "morph_prop_rates[1]" %in% colnames(params)
    if ( has_prop_rates == TRUE ) {
        discrete_rates = discrete_rates * params[[paste0("morph_prop_rates[",part,"]")]]
    }

    # determine the Q matrix
    has_mix_cats = "var_mix" %in% colnames(params)
    # mat = make_Q_matrix(rep(1, num_states))
    mat = matrix(1, 2, 2)
    diag(mat) = -1
    if ( has_mix_cats == TRUE ) {

        # get the parameter
        alpha = params[["alpha_mix"]]

        # discretize the beta distribution
        disc_rates     = discretize_beta(alpha, alpha, num_mix_cats)
        discrete_pis   = lapply(disc_rates, function(x) c(x, 1-x) )

    } else {

        discrete_pis = list(rep(1 / num_states, num_states))

    }

    # simulate one realization for each character
    sim_mat = matrix(NA, nrow=ntax, ncol=nchar)
    for(i in 1:nchar) {

        # create the data phyDat object
        dat = phyDat(obs[,i], type="USER", levels=c("0","1"))

        # get the gap pattern for this character
        missing = obs[,i] %in% c("-","?")

        # compute the log-likelihood in each mixture category
        lik_mat = matrix(NA, nrow=length(discrete_rates), ncol=length(discrete_pis))
        for(j in 1:length(discrete_rates)) {
            local_rate = discrete_rates[j]
            for(k in 1:length(discrete_pis)) {
                local_pi = discrete_pis[[k]]
                fit = pml(tree, dat, bf=local_pi, rate=local_rate)
                lik_mat[j,k] = fit$logLik
            }
        }

        # relativize and exponentiate
        lik_mat = exp(lik_mat - max(lik_mat))

        # choose the mixture category
        this_rate = sample.int(length(discrete_rates), size=1, prob=rowSums(lik_mat))    # marginalize over pi to choose the rate
        this_pi   = sample.int(length(discrete_pis),   size=1, prob=lik_mat[this_rate,]) # choose pi conditional on rate

        # simulate this character
        sim_mat[,i] = simulate_character(tree, discrete_rates[this_rate], mat, discrete_pis[[this_pi]], missing, taxa)

    }

    # name the rows by taxon
    rownames(sim_mat) = taxa

    # return the simulated matrix
    return(sim_mat)

}

simulate_character = function(tree, rate, matrix, pi, gap_pattern, taxa) {

    # repeat until variable
    repeat {

        # simulate the sequence
        sim = unlist(simSeq(tree, l=1, bf=pi, type="USER", levels=as.character(1:length(pi)-1), rate=rate))

        # re-order
        sim = sim[taxa]

        # check that it is variable
        if (length(unique(sim[gap_pattern == FALSE])) > 1) {
            break
        }

    }

    # drop missing data
    sim = sim - 1
    sim[gap_pattern] = "?"

    # return simulation
    return(sim)

}




simulate_stochastic_map = function(obs, tree, params, part=1, num_gamma_cats=4, num_mix_cats=5, num_states=2, nsims=1) {

    # get some statistics
    taxa  = tree$tip.label
    obs   = obs[taxa,]
    nchar = ncol(obs)
    ntax  = nrow(obs)
    obs[obs == "-"] = "?"

    # determine if this model has among-character rate variation
    has_morph_gamma = "gamma_alpha_morph[1]" %in% colnames(params)
    if ( has_morph_gamma == TRUE ) {

        # get the parameter
        alpha = params[[paste0("gamma_alpha_morph[",part,"]")]]

        # get the rate categories
        discrete_rates = discrete.gamma(alpha, num_gamma_cats)

        stop("implement gamma rates")

    } else {

        discrete_rates = 1.0

    }

    # determine the morphology partition rate-multiplier (if present)
    has_prop_rates = "morph_prop_rates[1]" %in% colnames(params)
    if ( has_prop_rates == TRUE ) {
        discrete_rates = discrete_rates * mean(params[[paste0("morph_prop_rates[",part,"]")]])
    }

    # determine the Q matrix
    has_mix_cats = "var_mix" %in% colnames(params)
    # mat = make_Q_matrix(rep(1, num_states))
    mat = matrix(1, 2, 2)
    diag(mat) = -1
    if ( has_mix_cats == TRUE ) {

        # get the parameter
        alpha = mean(params[["alpha_mix"]])

        # discretize the beta distribution
        disc_rates     = discretize_beta(alpha, alpha, num_mix_cats)
        discrete_pis   = lapply(disc_rates, function(x) c(x, 1-x) )

    } else {

        discrete_pis = list(rep(1 / num_states, num_states))

    }

    # simulate one map for each character
    sims = vector("list", nchar)
    for(i in 1:nchar) {

        cat("Simulating character ", i, "\n", sep="")

        # create the data phyDat object
        dat = phyDat(obs[,i], type="USER", levels=c("0","1"))

        # get the gap pattern for this character
        missing = obs[,i] %in% c("-","?")

        # get the data
        this_dat = obs[,i]

        # compute probabilities for each state
        this_dat_mat = matrix(0, nrow=nrow(obs), ncol=part + 1)
        rownames(this_dat_mat) = rownames(obs)
        colnames(this_dat_mat) = 0:part
        for(j in 1:length(this_dat)) {
            this_datum = this_dat[j]
            this_taxon = names(this_datum)
            if ( this_datum == "?" ) {
                this_dat_mat[j,] = 1 / (part + 1)
            } else {
                this_dat_mat[this_taxon,this_datum] = 1
            }
        }

        # if the data are invariant, skip
        # if (length(table(this_dat)) < 2) {
        #     next
        # }

        # compute the log-likelihood in each mixture category
        lik_mat = matrix(NA, nrow=length(discrete_rates), ncol=length(discrete_pis))
        for(j in 1:length(discrete_rates)) {
            local_rate = discrete_rates[j]
            for(k in 1:length(discrete_pis)) {
                local_pi = discrete_pis[[k]]
                fit = pml(tree, dat, bf=local_pi, rate=local_rate)
                lik_mat[j,k] = fit$logLik
            }
        }

        # relativize and exponentiate
        lik_mat = exp(lik_mat - max(lik_mat))

        these_sims = vector("list", nsims)
        bar = txtProgressBar(style=3, width=40)
        for(j in 1:nsims) {

            # choose the mixture category
            this_rate_cat = sample.int(length(discrete_rates), size=1, prob=rowSums(lik_mat))        # marginalize over pi to choose the rate
            this_pi_cat   = sample.int(length(discrete_pis),   size=1, prob=lik_mat[this_rate_cat,]) # choose pi conditional on rate

            # get the parameters
            this_rate = discrete_rates[this_rate_cat]
            this_pi   = discrete_pis[[this_pi_cat]]
            this_Q    = make_Q_matrix(this_pi)

            rownames(this_Q) = colnames(this_Q) = as.character(0:part)
            names(this_pi)   = as.character(0:part)

            # simulate this character
            these_sims[[j]] = make.simmap(tree, this_dat_mat, model=this_rate * this_Q, pi=this_pi, Q=this_rate * this_Q, message=FALSE)

            setTxtProgressBar(bar, j / nsims)

        }
        sims[[i]] = these_sims

        cat("\n")

    }

    # return the simulated matrix
    return(sims)

}

discretize_beta = function(alpha, beta, num_cats) {

    factor = (alpha / (alpha + beta)) * num_cats
    q = qbeta(1:(num_cats-1) / num_cats, alpha, beta)
    p = pbeta( q, alpha + 1, beta )
    p[num_cats] = 1.0
    p[-1] = p[-1] - p[-num_cats]
    p = p * factor
    return(p)

}

make_Q_matrix = function(pi) {

    # make the matrix
    matrix = matrix(pi, length(pi), length(pi), byrow=TRUE)
    diag(matrix) = 0
    diag(matrix) = -rowSums(matrix)

    # normalize the matrix
    matrix = matrix / -sum(diag(matrix) * pi)

    return(matrix)

}







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
        # attachment_point = tmp_tree$Nnode + which(tmp_tree$node.label == this_sampled_ancestor) - 1

        # attach the tip
        tmp_tree = bind.tip(tmp_tree, this_sampled_ancestor, edge.length=0.0, where=attachment_point)

    }

    # recover()
    #
    # par(mfrow=1:2)
    # plot(tree, no.margin=TRUE, cex=0.5)
    # nodelabels(pch=19, cex=0.5)
    # # plot(multi2di(tmp_tree), no.margin=TRUE, cex=0.5)
    # plot(collapse.singles(tmp_tree), no.margin=TRUE, cex=0.5)
    # nodelabels(pch=19, cex=0.5)
    #
    # drop.tip(tmp_tree, tip="", collapse.singles=TRUE)
    #
    # has.singles(tmp_tree)
    # collapse.singles(tmp_tree)

    #  return the tree
    return(collapse.singles(tmp_tree))

}



