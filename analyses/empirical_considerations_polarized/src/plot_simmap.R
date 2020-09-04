plot_simmap = function(time_tree, tree, data, simmaps, states, colors, nt=1001, show.tip.label=FALSE, edge.width=0, plot_pie=TRUE, pie_size=1.5, label.offset=0, label.cex=1, lwd=1, ...) {

    # compute dt
    dt = sum(time_tree$edge.length) / nt

    # for each branch, compute the probability at each time slice
    num_branches = nrow(tree$edge)
    num_sims     = length(simmaps)
    num_states   = length(states)

    if ( missing(colors) ) {
        colors = 1:num_states
    }
    col_rgb = col2rgb(colors)

    branch_x0 = vector("list", num_branches)
    branch_x1 = vector("list", num_branches)
    branch_sample_colors = vector("list", num_branches)
    for(i in 1:num_branches) {

        # compute the scale factor
        factor = time_tree$edge.length[i] / tree$edge.length[i]

        # compute the time points
        sampled_times = seq(0, time_tree$edge.length[i], by=dt)

        # compute the states per time slice
        sampled_states = vector("list", num_sims)
        for(j in 1:num_sims) {
            this_map = simmaps[[j]]$maps[[i]] * factor
            this_map_cum = cumsum(this_map)
            this_map_states = names(this_map_cum)
            sampled_states[[j]] = this_map_states[findInterval(sampled_times, this_map_cum, left.open=TRUE)+1]
        }
        sampled_states = do.call(rbind, sampled_states)

        # compute the probabilities per time slice
        state_probs = matrix(NA, nrow=num_states, ncol=length(sampled_times))
        for(j in 1:num_states) {
            state_probs[j,] = colMeans(sampled_states == states[j])
        }

        # compute the colors per time slice
        segment_colors = numeric(length(sampled_times))
        for(j in 1:length(sampled_times)) {
            these_cols = col_rgb %*% state_probs[,j]
            segment_colors[j] = rgb(these_cols[1,1], these_cols[2,1], these_cols[3,1], maxColorValue = 255)
        }

        # store computed values
        branch_x0[[i]] = sampled_times
        branch_x1[[i]] = c(sampled_times[-1], time_tree$edge.length[i])
        branch_sample_colors[[i]] = segment_colors

    }

    # get the coordinates for the tree
    yy = node.height(time_tree)
    xx = node.depth.edgelength(time_tree)

    # plot the tree
    pp = plot(time_tree, type = "phylogram",
              use.edge.length = TRUE,
              show.tip.label  = show.tip.label,
              edge.width      = edge.width,
              direction       = "rightwards",
              cex             = label.cex,
              label.offset    = label.offset)

    # plot the segments
    for(i in 1:num_branches) {

        # get the horiztonal segments
        this_x0  = branch_x0[[i]] + xx[ time_tree$edge[i,1] ]
        this_x1  = branch_x1[[i]] + xx[ time_tree$edge[i,1] ]
        this_y0  = yy[ time_tree$edge[i,2] ]
        this_col = branch_sample_colors[[i]]
        segments(x0=this_x0, x1=this_x1, y0=this_y0, col=this_col, ...)

        # get the vertical segment
        this_x0  = xx[ time_tree$edge[i,1] ]
        this_y0  = yy[ time_tree$edge[i,1] ]
        this_y1  = yy[ time_tree$edge[i,2] ]

        this_col = branch_sample_colors[[i]][1]
        segments(x0=this_x0, y0=this_y0, y1=this_y1, col=this_col, ...)

    }

    # plot data pie charts
    if ( plot_pie == FALSE ) {
        return()
    }

    for(i in 1:nrow(data)) {

        this_taxon = rownames(data)[i]
        this_data  = data[i,]
        this_index = which(tree$tip.label == this_taxon)

        this_x = xx[ this_index ]
        this_y = yy[ this_index ]

        floating.pie(xpos = this_x, ypos = this_y, x=this_data+1e-10, col=colors, radius=pie_size, lwd=ifelse(edge.width > lwd, edge.width - lwd, lwd)  )

    }


}













