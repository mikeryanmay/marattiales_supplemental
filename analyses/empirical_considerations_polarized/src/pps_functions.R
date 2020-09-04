library(rncl)
library(stringr)

pps_statistics = function(obs, sim, tree, nstates) {

    # compute the total number of characters
    num_chars = sum(sapply(obs, ncol))

    # container for scores
    obs_scores = numeric(num_chars)
    sim_scores = numeric(num_chars)

    # iterate over partitions
    num_partitions = length(obs)
    part_obs_scores = vector("list", num_partitions)
    part_sim_scores = vector("list", num_partitions)
    for(i in 1:num_partitions) {

        # get the data
        this_obs_data = phyDat(obs[[i]], type="USER", levels=paste0(1:nstates[i]-1))
        this_sim_data = phyDat(sim[[i]], type="USER", levels=paste0(1:nstates[i]-1))

        # calculate parsimony scores
        part_obs_scores[[i]] = parsimony(tree, this_obs_data, site="site")
        part_sim_scores[[i]] = parsimony(tree, this_sim_data, site="site")

    }

    # for each partititon
    sum_stats = numeric(num_partitions)
    var_stats = numeric(num_partitions)
    for(i in 1:num_partitions) {
        sum_stats[i] = sum(part_sim_scores[[i]]) - sum(part_obs_scores[[i]])
        var_stats[i] = var(part_sim_scores[[i]]) - var(part_obs_scores[[i]])
    }

    # concatenate among partitions
    part_obs_scores = unlist(part_obs_scores)
    part_sim_scores = unlist(part_sim_scores)

    # compute statistics
    sum_stat_cat = sum(part_sim_scores) - sum(part_obs_scores)
    var_stat_cat = var(part_sim_scores) - var(part_obs_scores)

    # create the result
    res = data.frame(sum_stat_cat = sum_stat_cat, var_stat_cat = var_stat_cat)

    # append the stat for each partition
    for(i in 1:num_partitions) {
        res[paste0("sum_stat_",i)] = sum_stats[i]
        res[paste0("var_stat_",i)] = var_stats[i]
    }

    return(res)

}


readParsedData = function(filename) {

    # read the lines
    lines = readLines(filename)

    # drop any comments
    lines = lines[grepl("\\[", lines) == FALSE]

    # find all the blocks
    start_lines  = which(toupper(lines) == "MATRIX") + 1
    end_lines    = which(toupper(lines) == "END;") - 2
    symbol_lines = which(grepl("SYMBOLS=", toupper(lines)))
    nchars_lines = which(grepl("NCHAR=", toupper(lines)))

    # for each block, extract the data
    alns = vector("list", length(start_lines))
    nstates = numeric(length(start_lines))
    for(i in 1:length(start_lines)) {

        # get the number of states
        these_symbol_lines = lines[symbol_lines[i]]
        nstates[i] = nchar(strsplit(these_symbol_lines, "\\\"")[[1]][2])

        # get the lines for this block
        these_lines = lines[start_lines[i]:end_lines[i]]

        # split the taxa and the data
        these_lines = strsplit(these_lines, "   ")

        # get the taxa
        taxa = gsub(" ","", sapply(these_lines, function(x) x[[1]]))

        # get the data
        data = do.call(rbind, lapply(these_lines, function(x) {
            tmp = gsub(" ","",x[[2]])
            # tmp = gsub("\\(01\\)","?",tmp)
            tmp = gsub("\\(.[0-9]*\\)","?",tmp)
            return(strsplit(tmp,"")[[1]])
        }))

        rownames(data) = taxa

        alns[[i]] = data

    }

    names(alns) = nstates

    return(alns)

}

readDataAsProbs = function(filename) {

    # read the lines
    lines = readLines(filename)

    # drop any comments
    lines = lines[grepl("\\[", lines) == FALSE]

    # find all the blocks
    start_lines  = which(toupper(lines) == "MATRIX") + 1
    end_lines    = which(toupper(lines) == "END;") - 2
    symbol_lines = which(grepl("SYMBOLS=", toupper(lines)))
    nchars_lines = which(grepl("NCHAR=", toupper(lines)))

    # for each block, extract the data
    alns = vector("list", length(start_lines))
    nstates = numeric(length(start_lines))
    for(i in 1:length(start_lines)) {

        # get the number of states
        these_symbol_lines = lines[symbol_lines[i]]
        nstates[i] = nchar(str_match(toupper(these_symbol_lines), "SYMBOLS=\\s*\\\"(\\d+)\\\"")[1,2])
        symbols = strsplit(str_match(toupper(these_symbol_lines), "SYMBOLS=\\s*\\\"(\\d+)\\\"")[1,2],"")[[1]]

        # get the number of characters
        these_nchar_lines = lines[nchars_lines[i]]
        nchars = as.numeric(str_match(toupper(these_nchar_lines), "NCHAR=\\s*(\\d+)")[,2])
        ntaxa  = as.numeric(str_match(toupper(these_nchar_lines), "NTAX=\\s*(\\d+)")[,2])

        # get the lines for this block
        these_lines = lines[start_lines[i]:end_lines[i]]

        # split the taxa and the data
        split_data = str_split_fixed(these_lines, " ", 2)
        taxa  = split_data[,1]
        chars = split_data[,2]

        # strip white space from characters
        chars = gsub(" ", "", chars)
        chars = gsub("-","?", chars)

        # create the data array
        data = array(0, dim=c(ntaxa, nchars, nstates[i]), dimnames=list(taxa, 1:nchars, symbols) )

        # extract the data
        for(t in 1:ntaxa) {

            # get the taxon data
            this_data = chars[t]
            for(c in 1:nchars) {

                # if the next character is a parens, ambiguous
                nextchar = substr(this_data, 1, 1)
                if ( nextchar == "(" ) {

                    # take off the first character
                    this_data = sub(".", "", this_data)

                    # extract up to the next close parens
                    ambig_string = str_match(this_data, "(\\d+)\\)")[1,2]
                    ambig_string = strsplit(ambig_string,"")[[1]]
                    num_ambig_values = length(ambig_string)
                    data[t,c,ambig_string] = 1 / num_ambig_values

                    # chop off everything up to the next close parens
                    this_data = sub(paste0("[A-Za-z0-9_]{",num_ambig_values,"}\\)"), "", this_data)

                } else {

                    # assign the probability
                    if ( nextchar %in% c("?","-") ) {
                        data[t,c,] = 1 / nstates[i]
                    } else {
                        data[t,c,nextchar] = 1
                    }

                    # take off the first character
                    this_data = sub(".", "", this_data)

                }


            }

        }

        # store the data array
        alns[[i]] = data

    }

    names(alns) = nstates

    return(alns)

}

readDataAsPhyDat = function(filename) {

    # read the lines
    lines = readLines(filename)

    # drop any comments
    lines = lines[grepl("\\[", lines) == FALSE]

    # find all the blocks
    start_lines  = which(toupper(lines) == "MATRIX") + 1
    end_lines    = which(toupper(lines) == "END;") - 2
    symbol_lines = which(grepl("SYMBOLS=", toupper(lines)))
    nchars_lines = which(grepl("NCHAR=", toupper(lines)))

    # for each block, extract the data
    alns = vector("list", length(start_lines))
    nstates = numeric(length(start_lines))
    for(i in 1:length(start_lines)) {

        # get the number of states
        these_symbol_lines = lines[symbol_lines[i]]
        nstates[i] = nchar(str_match(toupper(these_symbol_lines), "SYMBOLS=\\s*\\\"(\\d+)\\\"")[1,2])
        symbols = strsplit(str_match(toupper(these_symbol_lines), "SYMBOLS=\\s*\\\"(\\d+)\\\"")[1,2],"")[[1]]

        # get the number of characters
        these_nchar_lines = lines[nchars_lines[i]]
        nchars = as.numeric(str_match(toupper(these_nchar_lines), "NCHAR=\\s*(\\d+)")[,2])
        ntaxa  = as.numeric(str_match(toupper(these_nchar_lines), "NTAX=\\s*(\\d+)")[,2])

        # get the lines for this block
        these_lines = lines[start_lines[i]:end_lines[i]]

        # split the taxa and the data
        split_data = str_split_fixed(these_lines, " ", 2)
        taxa  = split_data[,1]
        chars = split_data[,2]

        # strip white space from characters
        chars = gsub(" ", "", chars)
        chars = gsub("-","?", chars)

        # create the contrast matrix
        total_num_states = 2^nstates[i] - 1
        contrast = do.call(rbind, lapply(binSeq(1:total_num_states, n=nstates[i]), function(x) as.integer(as.logical(x)) ))
        contrast = contrast[order(rowSums(contrast)),nstates[i]:1]

        if ( nstates[i] == 2 ) {
            total_symbols = c(symbols,"?")
        } else {
            total_symbols = c(symbols, (nstates[i] + 1):(nrow(contrast)-1),"?")
        }

        rownames(contrast) = total_symbols
        colnames(contrast) = symbols

        # create the data array
        data = matrix("", nrow=ntaxa, ncol=nchars, dimnames=list(taxa, 1:nchars) )

        # extract the data
        for(t in 1:ntaxa) {

            # get the taxon data
            this_data = chars[t]
            for(c in 1:nchars) {

                # if the next character is a parens, ambiguous
                nextchar = substr(this_data, 1, 1)
                if ( nextchar == "(" ) {

                    # take off the first character
                    this_data = sub(".", "", this_data)

                    # extract up to the next close parens
                    ambig_string = str_match(this_data, "(\\d+)\\)")[1,2]
                    ambig_string = strsplit(ambig_string,"")[[1]]
                    num_ambig_values = length(ambig_string)

                    # figure out the symbolic state
                    this_symbol = rownames(contrast)[min(which(rowSums(contrast[,ambig_string]) == num_ambig_values))]

                    # record this symbol
                    data[t,c] = this_symbol

                    # chop off everything up to the next close parens
                    this_data = sub(paste0("[A-Za-z0-9_]{",num_ambig_values,"}\\)"), "", this_data)

                } else {

                    # assign the character
                    data[t,c] = nextchar

                    # take off the first character
                    this_data = sub(".", "", this_data)

                }

            }

        }

        # convert to phyDat
        phy_data = phyDat(data, type="USER", contrast=contrast)

        # store the data array
        alns[[i]] = phy_data

    }

    names(alns) = nstates

    return(alns)

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

    #  return the tree
    return(collapse.singles(tmp_tree))

}





simulate_pps_dataset = function(obs, tree, params, num_states, part, num_gamma_cats=4, num_mix_cats=5) {

    # make sure taxa are consistent
    tree_taxa = tree$tip.label
    data_taxa = rownames(obs)

    # drop missing taxa from the tree
    missing_taxa = tree_taxa[tree_taxa %in% data_taxa == FALSE]
    tree = drop.tip(tree, missing_taxa)

    # get some statistics
    taxa  = tree$tip.label
    obs   = obs[taxa,]
    if ( class(obs) == "character" ) {
        obs = as.matrix(obs)
    }
    nchar = ncol(obs)
    ntax  = nrow(obs)

    # determine if this model has among-character rate variation
    has_morph_gamma = paste0("gamma_alpha_morph[",part,"]") %in% colnames(params)
    if ( has_morph_gamma == TRUE ) {

        # get the parameter
        alpha = params[[paste0("gamma_alpha_morph[",part,"]")]]

        # get the rate categories
        discrete_rates = discrete.gamma(alpha, num_gamma_cats)

    } else {

        discrete_rates = 1.0

    }

    # determine the morphology partition rate-multiplier (if present)
    has_prop_rates = paste0("morph_prop_rates[",part,"]") %in% colnames(params)
    if ( has_prop_rates == TRUE ) {
        discrete_rates = discrete_rates * params[[paste0("morph_prop_rates[",part,"]")]]
    }

    # determine the Q matrix
    # case 1: no mixture model
    # case 2: mixture model with pis and mixture weights
    # case 3: mixture model with discretized beta
    has_mix_cats = "var_mix" %in% colnames(params)
    first_mix_cat_index = num_mix_cats * (part - 1) - num_mix_cats + 1
    if ( has_mix_cats == FALSE ) {
        discrete_pis = list(rep(1 / num_states, num_states))
        discrete_pi_weights = 1
    } else if ( paste0("pi_morpho[",first_mix_cat_index,"][1]") %in% colnames(params) ) {

        # get the stationary frequencies
        col_prefix = paste0("pi_morpho[",first_mix_cat_index+1:num_mix_cats-1,"]")
        discrete_pis = vector("list", num_mix_cats)
        for(i in 1:num_mix_cats) {
            col_names = paste0(col_prefix[i],"[",1:num_states,"]")
            discrete_pis[[i]] = as.numeric(params[col_names])
        }

        # get the weights
        discrete_pi_weights = as.numeric(params[paste0("morph_matrix_weights[",part,"][",1:num_mix_cats,"]")])
        discrete_pi_weights = discrete_pi_weights / sum(discrete_pi_weights)

    } else {

        # get the parameter
        alpha = params[["alpha_mix"]]

        # discretize the beta distribution
        disc_pi_0    = discretize_beta(alpha, alpha, num_mix_cats)
        discrete_pis = lapply(disc_pi_0, function(x) c(x, 1-x) )
        discrete_pi_weights = rep(1, num_mix_cats) / num_mix_cats

    }

    # simulate one realization for each character
    sim_mat = matrix(NA, nrow=ntax, ncol=nchar)
    for(i in 1:nchar) {

        # create the data phyDat object
        dat = phyDat(obs[,i], type="USER", levels=as.character(1:num_states-1))

        # get the gap pattern for this character
        missing = obs[,i] %in% c("-","?")

        # compute the posterior probability of each mixture category
        pp_mat = matrix(NA, nrow=length(discrete_rates), ncol=length(discrete_pis))
        for(j in 1:length(discrete_rates)) {
            local_rate = discrete_rates[j]
            for(k in 1:length(discrete_pis)) {
                local_pi = discrete_pis[[k]]
                fit = pml(tree, dat, bf=local_pi, rate=local_rate)
                pp_mat[j,k] = fit$logLik + log(discrete_pi_weights[k])
            }
        }

        # relativize and exponentiate
        pp_mat = exp(pp_mat - max(pp_mat))

        # choose the mixture category
        this_rate = sample.int(length(discrete_rates), size=1, prob=rowSums(pp_mat))    # marginalize over pi to choose the rate
        this_pi   = sample.int(length(discrete_pis),   size=1, prob=pp_mat[this_rate,]) # choose pi conditional on rate

        # simulate this character
        sim_mat[,i] = simulate_character(tree, discrete_rates[this_rate], discrete_pis[[this_pi]], missing, taxa)

    }

    # name the rows by taxon
    rownames(sim_mat) = taxa

    # return the simulated matrix
    return(sim_mat)

}

simulate_character = function(tree, rate, pi, gap_pattern, taxa) {

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

simulate_stochastic_map = function(probdat, tree, params, part=1, num_gamma_cats=4, num_mix_cats=5, num_states=2, nsims=1) {

    # make sure taxa are consistent
    tree_taxa = tree$tip.label
    data_taxa = rownames(probdat)

    # drop missing taxa from the tree
    missing_taxa = tree_taxa[tree_taxa %in% data_taxa == FALSE]
    tree = drop.tip(tree, missing_taxa)

    # get some statistics
    taxa  = tree$tip.label
    nchar = ncol(probdat)
    ntax  = nrow(probdat)

    # determine if this model has among-character rate variation
    has_morph_gamma = paste0("gamma_alpha_morph[",part,"]") %in% colnames(params)
    if ( has_morph_gamma == TRUE ) {

        # get the parameter
        alpha = params[[paste0("gamma_alpha_morph[",part,"]")]]

        # get the rate categories
        discrete_rates = discrete.gamma(alpha, num_gamma_cats)

    } else {

        discrete_rates = 1.0

    }

    # determine the morphology partition rate-multiplier (if present)
    has_prop_rates = paste0("morph_prop_rates[",part,"]") %in% colnames(params)
    if ( has_prop_rates == TRUE ) {
        discrete_rates = discrete_rates * params[[paste0("morph_prop_rates[",part,"]")]]
    }

    # determine the Q matrix
    # case 1: no mixture model
    # case 2: mixture model with pis and mixture weights
    # case 3: mixture model with discretized beta
    has_mix_cats = "var_mix" %in% colnames(params)
    first_mix_cat_index = num_mix_cats * (part - 1) - num_mix_cats + 1
    if ( has_mix_cats == FALSE ) {
        discrete_pis = list(rep(1 / num_states, num_states))
        discrete_pi_weights = 1
    } else if ( paste0("pi_morpho[",first_mix_cat_index,"][1]") %in% colnames(params) ) {

        # get the stationary frequencies
        col_prefix = paste0("pi_morpho[",first_mix_cat_index+1:num_mix_cats-1,"]")
        discrete_pis = vector("list", num_mix_cats)
        for(i in 1:num_mix_cats) {
            col_names = paste0(col_prefix[i],"[",1:num_states,"]")
            discrete_pis[[i]] = as.numeric(params[col_names])
        }

        # get the weights
        discrete_pi_weights = as.numeric(params[paste0("morph_matrix_weights[",part,"][",1:num_mix_cats,"]")])
        discrete_pi_weights = discrete_pi_weights / sum(discrete_pi_weights)

    } else {

        # get the parameter
        alpha = params[["alpha_mix"]]

        # discretize the beta distribution
        disc_pi_0    = discretize_beta(alpha, alpha, num_mix_cats)
        discrete_pis = lapply(disc_pi_0, function(x) c(x, 1-x) )
        discrete_pi_weights = rep(1, num_mix_cats) / num_mix_cats

    }

    # simulate one map for each character
    sims = vector("list", nchar)
    for(i in 1:nchar) {

        cat("Simulating character ", i, "\n", sep="")

        # get the data as probabilities
        this_dat = probdat[,i,]

        # compute the log-likelihood in each mixture category
        lik_mat = matrix(NA, nrow=length(discrete_rates), ncol=length(discrete_pis))
        for(j in 1:length(discrete_rates)) {
            local_rate = discrete_rates[j]
            for(k in 1:length(discrete_pis)) {
                local_pi = discrete_pis[[k]]
                local_Q  = make_Q_matrix(local_pi)
                fit = fitMk(tree, this_dat, model=local_rate * local_Q, pi=local_pi, fixedQ=local_rate * local_Q)
                lik_mat[j,k] = fit$logLik * discrete_pi_weights[k]
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

            rownames(this_Q) = colnames(this_Q) = as.character(1:num_states-1)
            names(this_pi)   = as.character(1:num_states-1)

            # simulate this character
            these_sims[[j]] = make.simmap(tree, this_dat, model=this_rate * this_Q, pi=this_pi, Q=this_rate * this_Q, message=FALSE)

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
