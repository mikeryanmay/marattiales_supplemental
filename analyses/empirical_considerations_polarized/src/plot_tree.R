library(RevGadgets)
library(scales)

matchNodes = function(treedata, phy) {

  # get some useful info
  num_sampled_anc = sum(phy$node.label != "")
  num_tips        = length(phy$tip.label)
  num_nodes       = phy$Nnode
  sampled_ancs    = which(tabulate(phy$edge[,1]) == 1)
  tip_indexes     = 1:(num_tips + num_sampled_anc)
  node_indexes    = (num_tips + num_sampled_anc) + num_nodes:1

  node_map     = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  # current_node = phy$Nnode + 2 - num_sampled_anc
  current_node = num_tips + 1
  k = 1
  t = 1

  while(TRUE) {

    # compute the number of descendants of this tip
    current_num_descendants = sum(phy$edge[,1] == current_node)

    if ( current_node <= num_tips ) {

      treedata_node = which(as.character(treedata@data$node) == current_node)
      node_map$Rev[node_map$R == current_node] = as.numeric(treedata@data[treedata_node,]$index)
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1

    } else if ( current_node %in% sampled_ancs ) {

      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1

      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }

    } else {

      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1

      # if (current_num_descendants > 2) {
      #     cat("POLYTOMY")
      #     break
      # }

      num_visits = node_map$visits[node_map$R == current_node]

      if ( num_visits <= current_num_descendants ) {
          # go to next descendant
          current_node = phy$edge[phy$edge[,1] == current_node,2][current_num_descendants - num_visits + 1]
      } else if ( num_visits > current_num_descendants ) {
          # go down
          if (current_node == num_tips + 1) {
              break
          } else {
              current_node = phy$edge[phy$edge[,2] == current_node,1]
          }
      }

      # if ( node_map$visits[node_map$R == current_node] == 1 ) {
      #   # go right
      #   current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      # } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
      #   # go left
      #   current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      # } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
      #   # go down
      #   if (current_node == num_tips + 1) {
      #     break
      #   } else {
      #     current_node = phy$edge[phy$edge[,2] == current_node,1]
      #   }
      # }

    }

  }

  return(node_map[,1:2])

}

plot_tree = function(tree, bar_colors, bar_cuts, bar_alpha = 0.9, timescale = NULL, label_offset=1, label_cex=0.5, bar_lwd=2, anc_pch=20, anc_cex=0.5, anc_lab_adj=0.75, anc_lab_cex=0.5, anc_leg_pos = "bottomleft", legend_inset=0, ...) {

  # get the phylo
  phylo    = tree@phylo
  node_map = matchNodes(tree, phylo)
  tree@data$node = as.character(node_map[match(as.numeric(tree@data$index), node_map$Rev),]$R)

  # compute the bars for each node
  nodes = as.numeric(tree@data$node)
  ages  = tree@data$age_0.95_HPD

  # get the types of nodes
  sampled_ancestor_nodes = which(tabulate(phylo$edge[,1]) == 1)
  internal_nodes         = sort(unique(phylo$edge[phylo$edge[,1] %in% sampled_ancestor_nodes == FALSE,1]))
  tips                   = phylo$edge[phylo$edge[,2] %in% phylo$edge[,1] == FALSE,2]
  extant_tips            = tips[sapply(tree@data[nodes %in% tips,]$age_0.95_HPD, function(x) any(is.na(x)))]
  extinct_tips           = tips[sapply(tree@data[nodes %in% tips,]$age_0.95_HPD, function(x) any(is.na(x)) == FALSE)]

  # internal node bars
  internal_node_bars = vector("list", length(internal_nodes))
  internal_node_prob = numeric(length(internal_nodes))
  for(i in 1:length(internal_nodes)) {

    # get the node
    this_node = internal_nodes[i]
    this_idx  = which(nodes == this_node)

    # get the values
    this_bar = as.numeric(ages[this_idx][[1]])
    this_pp  = as.numeric(tree@data$posterior[this_idx])

    # store
    internal_node_bars[[i]] = this_bar
    internal_node_prob[i]   = this_pp

  }

  # tip bars
  if ( length(extinct_tips) > 0 ) {

      extinct_tip_bars = vector("list", length(extinct_tips))
      extinct_tip_prob = numeric(length(extinct_tips))
      for(i in 1:length(extinct_tips)) {

          # get the node
          this_node = extinct_tips[i]
          this_idx  = which(nodes == this_node)

          # get the values
          this_bar = as.numeric(ages[this_idx][[1]])
          if ( "sampled_ancestor" %in% colnames(tree@data) == FALSE ) {
              this_pp = 0.0
          } else {
              this_pp  = as.numeric(tree@data$sampled_ancestor[this_idx])
              if ( is.na(this_pp) ) {
                  this_pp = 0.0
              }
          }

          # store
          extinct_tip_bars[[i]] = this_bar
          extinct_tip_prob[i]   = this_pp

      }

  }

  # sampled ancestors
  if ( length(sampled_ancestor_nodes) > 0 ) {

      sampled_ancestor_bars = vector("list", length(sampled_ancestor_nodes))
      sampled_ancestor_prob = numeric(length(sampled_ancestor_nodes))
      for(i in 1:length(sampled_ancestor_nodes)) {

          # get the node
          this_node = sampled_ancestor_nodes[i]
          this_idx  = which(nodes == this_node)

          # get the values
          this_bar = as.numeric(ages[this_idx][[1]])
          this_pp  = as.numeric(tree@data$sampled_ancestor[this_idx])
          if ( is.na(this_pp) ) {
              this_pp = 0.0
          }

          # store
          sampled_ancestor_bars[[i]] = this_bar
          sampled_ancestor_prob[i]   = this_pp

      }

  }


  # get the coordinates for the tree
  phylo = ladderize(phylo)
  p = plot.phylo(phylo, plot=FALSE, label.offset=label_offset, show.tip.label=TRUE, ...)
  yy = node.height(phylo)
  xx = node.depth.edgelength(phylo)
  x_max = internal_node_bars[[1]][2] - max(xx)
  x_lim = p$x.lim
  x_lim[1] = -x_max

  # plot bars for timescale
  if ( is.null(timescale) == FALSE ) {

    # reverse the axis
    new_lim = c(internal_node_bars[[1]][2], -(diff(x_lim) - internal_node_bars[[1]][2]))

    par(new=TRUE)
    plot(rev(x_lim), type="n", bty="n", xlab=NA, ylab=NA, xaxt="n", yaxt="n", xlim=new_lim)

    # plot the epoch
    for(i in seq(2, nrow(timescale), by=2)) {

      this_epoch = timescale[i,]
      polygon( x = c(this_epoch$Start, this_epoch$End, this_epoch$End, this_epoch$Start),
               y = c(-1000, -1000, 1000, 1000),
               border=NA, col="grey90")

    }

  }

  # plot the tree
  par(new=TRUE)
  p = plot(phylo, plot=TRUE, label.offset=label_offset, show.tip.label=FALSE, x.lim=x_lim, ...)

  # plot internal node bars
  node_y0      = numeric(length(internal_nodes))
  node_x0      = numeric(length(internal_nodes))
  node_x1      = numeric(length(internal_nodes))
  node_bar_col = character(length(internal_nodes))
  for(i in 1:length(internal_nodes)) {

    # get the node
    this_node = internal_nodes[i]

    # get the bar
    this_bar = max(xx) - internal_node_bars[[i]]
    this_pp  = internal_node_prob[i]

    # get the coordinates
    node_y0[i] = yy[this_node]
    node_x0[i] = this_bar[1]
    node_x1[i] = this_bar[2]

    # get the color
    node_bar_col[i] = bar_colors[findInterval(this_pp, bar_cuts, left.open = TRUE) + 1]

  }

  segments(x0 = node_x0 + 1 - bar_lwd, x1 = node_x1 - 1 + bar_lwd, y0 = node_y0, lwd=bar_lwd, col=alpha(node_bar_col, bar_alpha))
  segments(x0 = node_x0, y0 = node_y0 - 0.25, y1 = node_y0 + 0.25, col=node_bar_col, lwd=1)
  segments(x0 = node_x1, y0 = node_y0 - 0.25, y1 = node_y0 + 0.25, col=node_bar_col, lwd=1)
  # arrow_width = 0.25 * diff(grconvertY(c(1,2), "user", "inches"))
  # arrows(x0 = node_x0, x1 = node_x1, y0 = node_y0, lwd=bar_lwd, col=alpha(node_bar_col, bar_alpha), length=arrow_width, code=3, angle=90)

  # plot tip bars
  if ( length(extinct_tips) > 0 ) {

      tip_y0      = numeric(length(extinct_tips))
      tip_x0      = numeric(length(extinct_tips))
      tip_x1      = numeric(length(extinct_tips))
      tip_mean    = numeric(length(extinct_tips))
      tip_name    = character(length(extinct_tips))
      tip_bar_col = character(length(extinct_tips))
      for(i in 1:length(extinct_tips)) {

          # get the tip
          this_tip = extinct_tips[i]

          # get the bar
          this_bar = max(xx) - extinct_tip_bars[[i]]
          this_pp  = 1 - extinct_tip_prob[i]

          # get the coordinates
          tip_y0[i]   = yy[this_tip]
          tip_x0[i]   = this_bar[1]
          tip_x1[i]   = this_bar[2]
          tip_mean[i] = xx[this_tip]
          tip_name[i] = gsub("_", " ",phylo$tip.label[this_tip])

          # get the color
          tip_bar_col[i] = bar_colors[findInterval(this_pp, bar_cuts, left.open = TRUE) + 1]

      }

      points(x = tip_mean, y=tip_y0, pch=anc_pch, cex=anc_cex, col=tip_bar_col)
      segments(x0 = tip_x0 + 1 - bar_lwd, x1 = tip_x1 - 1 + bar_lwd, y0 = tip_y0, lwd=bar_lwd, col=alpha(tip_bar_col, bar_alpha))
      segments(x0 = tip_x0, y0 = tip_y0 - 0.25, y1 = tip_y0 + 0.25, col=tip_bar_col, lwd=1)
      segments(x0 = tip_x1, y0 = tip_y0 - 0.25, y1 = tip_y0 + 0.25, col=tip_bar_col, lwd=1)

      # arrows(x0 = tip_x0, x1 = tip_x1, y0 = tip_y0, lwd=bar_lwd, col=alpha(tip_bar_col, bar_alpha), length=arrow_width, code=3, angle=90)
      text(x = tip_x0 + label_offset, y = tip_y0, label=tip_name, adj=0, cex=label_cex, xpd=NA, font=3)

  }

  # plot the tip labels
  tip_y0   = numeric(length(extant_tips))
  tip_x0   = numeric(length(extant_tips))
  tip_name = character(length(extant_tips))
  for(i in 1:length(extant_tips)) {

    # get the tip
    this_tip = extant_tips[i]

    # get the coordinates
    tip_y0[i]   = yy[this_tip]
    tip_x0[i]   = xx[this_tip]
    tip_name[i] = gsub("_", " ",phylo$tip.label[this_tip])

  }

  text(x = tip_x0 + label_offset, y = tip_y0, label=tip_name, adj=0, cex=label_cex, xpd=NA, font=3)

  # plot ancestor bars
  if ( length(sampled_ancestor_nodes) > 0 ) {

      anc_y0      = numeric(length(sampled_ancestor_nodes))
      anc_x0      = numeric(length(sampled_ancestor_nodes))
      anc_x1      = numeric(length(sampled_ancestor_nodes))
      anc_mean    = numeric(length(sampled_ancestor_nodes))
      anc_names   = character(length(sampled_ancestor_nodes))
      anc_bar_col = character(length(sampled_ancestor_nodes))
      for(i in 1:length(sampled_ancestor_nodes)) {

        # get the anc
        this_anc = sampled_ancestor_nodes[i]

        # get the bar
        this_bar = max(xx) - sampled_ancestor_bars[[i]]
        this_pp  = sampled_ancestor_prob[i]

        # get the coordinates
        anc_y0[i]    = yy[this_anc]
        anc_x0[i]    = this_bar[1]
        anc_x1[i]    = this_bar[2]
        anc_mean[i]  = xx[this_anc]
        anc_names[i] = gsub("_", " ",phylo$node.label[this_anc - length(phylo$tip.label)])

        # get the color
        anc_bar_col[i] = bar_colors[findInterval(this_pp, bar_cuts, left.open = TRUE) + 1]

      }

      points(x = anc_mean, y = anc_y0, pch=anc_pch, cex=anc_cex, col=anc_bar_col)
      segments(x0 = anc_x0 + 1 - bar_lwd, x1 = anc_x1 - 1 + bar_lwd, y0 = anc_y0, lwd=bar_lwd, col=alpha(anc_bar_col, bar_alpha))
      segments(x0 = anc_x0, y0 = anc_y0 - 0.25, y1 = anc_y0 + 0.25, col=anc_bar_col, lwd=1)
      segments(x0 = anc_x1, y0 = anc_y0 - 0.25, y1 = anc_y0 + 0.25, col=anc_bar_col, lwd=1)

      # segments(x0 = tip_x0 + 1 - bar_lwd, x1 = tip_x1 - 1 + bar_lwd, y0 = tip_y0, lwd=bar_lwd, col=alpha(tip_bar_col, bar_alpha))
      # segments(x0 = tip_x0, y0 = tip_y0 - 0.25, y1 = tip_y0 + 0.25, col=tip_bar_col, lwd=1)
      # segments(x0 = tip_x1, y0 = tip_y0 - 0.25, y1 = tip_y0 + 0.25, col=tip_bar_col, lwd=1)

      text(x = anc_mean, y = anc_y0 + anc_lab_adj, label=1:length(sampled_ancestor_nodes), cex=anc_lab_cex)
      if ( is.null(timescale)  ) {
        legend(anc_leg_pos, legend=paste0(1:length(anc_names), ": ", c(anc_names)), bty="n", text.font=3, cex=label_cex, inset=legend_inset)
      } else {
        legend(anc_leg_pos, legend=paste0(1:length(anc_names), ": ", c(anc_names)), bg="white", text.font=3, cex=label_cex, inset=legend_inset)
      }

  }

  # re-orient the axes
  if ( is.null(timescale) == FALSE ) {

    # reverse the axis
    new_lim = c(internal_node_bars[[1]][2], -(diff(x_lim) - internal_node_bars[[1]][2]))

    par(new=TRUE)
    plot(rev(x_lim), type="n", bty="n", xlab=NA, ylab=NA, xaxt="n", yaxt="n", xlim=new_lim)

  }

}






