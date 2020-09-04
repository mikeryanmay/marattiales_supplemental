library(gplots)
library(coda)

lttPlot = function(trees,
                   CI=0.95, CI_lty=2, CI_type="polygon",
                   num_bins=10001,
                   log=TRUE, add=FALSE, verbose=TRUE, ...) {

  num_trees = length(trees)

  # compute the node heights
  node_heights = do.call(rbind, lapply(trees, function(x) {

    # compute the heights
    heights = nodeHeights(x)
    heights = abs(heights - max(heights))
    heights = round(heights, digits=10)

    # node and tip indexes
    tip_idx  = 1:length(x$tip.label)
    node_idx = 1:x$Nnode + length(x$tip.label)

    # compute the branching times
    branching_times = sort(unique(heights[,1]), decreasing=FALSE)

    # compute the tip ages
    tip_ages = sort(unique(heights[x$edge[,2] %in% tip_idx,2]), decreasing=FALSE)

    return(data.frame(na = I(list(branching_times)), ta = I(list(tip_ages))))

  }))

  # get the points to compute times at
  max_time = max(unlist(node_heights))
  times = seq(0, max_time, length.out=num_bins)

  # # re-orient the branching times
  # branching_times = lapply(branching_times, function(x) round(sort(max_time - x),15))

  # for each tree, compute the number of lineages at each time
  num_lineages = vector("list", length(trees))
  if(verbose) bar = txtProgressBar(style=3, width=40)
  for(i in 1:num_trees) {

    num_extinct = length(node_heights$ta[[i]]) - 1
    num_extant  = length(trees[[i]]$tip.label) - num_extinct

    total_intervals = findInterval(times, node_heights$na[[i]]) + 1
    total_num_per_intervals = ((num_extinct + num_extant + 1):1)[total_intervals]

    extinct_intervals = findInterval(times, node_heights$ta[[i]])
    extinct_num_per_intervals = (num_extinct:0)[extinct_intervals]

    # findInterval(times[9774], node_heights$na[[i]])
    # findInterval(times[9774], node_heights$ta[[i]])
    # times[8779]

    num_lineages[[i]] = total_num_per_intervals - extinct_num_per_intervals - 1

    if(verbose) setTxtProgressBar(bar, i / num_trees)

  }
  num_lineages = do.call(rbind, num_lineages)

  # compute the mean and confidence interval
  mean      = colMeans(num_lineages)
  alpha     = (1 - CI) / 2
  probs     = c(alpha, 1 - alpha)
  quantiles = apply(num_lineages, 2, quantile, prob=probs)

  # plot the curve
  if(add == FALSE) {
    plot(x=times, y=mean, xlim=rev(range(times)), log=ifelse(log,"y",""), ...)
  } else {
    lines(x=times, y=mean, ...)
  }

  if(CI_type == "lines") {
    lines(x=times, y=quantiles[1,], lty=CI_lty, ...)
    lines(x=times, y=quantiles[2,], lty=CI_lty, ...)
  } else if (CI_type == "polygon") {
    color = paste0(col2hex(eval(match.call()$col)),"50")
    polygon(x=c(times,rev(times)), y=c((quantiles[1,]),rev(quantiles[2,])), border=NA, col=color)
  }

  # matplot(t(num_lineages), type="l", lty=1, col="red")

}

lttCurve = function(trees,
                    CI=0.95, CI_lty=2, CI_type="polygon",
                    num_bins=10001,
                    log=TRUE, add=FALSE, verbose=TRUE, ...) {

  num_trees = length(trees)

  # compute the node heights
  node_heights = do.call(rbind, lapply(trees, function(x) {

    # compute the heights
    heights = nodeHeights(x)
    heights = abs(heights - max(heights))
    heights = round(heights, digits=10)

    # node and tip indexes
    tip_idx  = 1:length(x$tip.label)
    node_idx = 1:x$Nnode + length(x$tip.label)

    # compute the branching times
    branching_times = sort(unique(heights[,1]), decreasing=FALSE)

    # compute the tip ages
    tip_ages = sort(unique(heights[x$edge[,2] %in% tip_idx,2]), decreasing=FALSE)

    return(data.frame(na = I(list(branching_times)), ta = I(list(tip_ages))))

  }))

  # get the points to compute times at
  max_time = max(unlist(node_heights))
  times = seq(0, max_time, length.out=num_bins)

  # # re-orient the branching times
  # branching_times = lapply(branching_times, function(x) round(sort(max_time - x),15))

  # for each tree, compute the number of lineages at each time
  num_lineages = vector("list", length(trees))
  if(verbose) bar = txtProgressBar(style=3, width=40)
  for(i in 1:num_trees) {

    num_extinct = length(node_heights$ta[[i]]) - 1
    num_extant  = length(trees[[i]]$tip.label) - num_extinct

    total_intervals = findInterval(times, node_heights$na[[i]]) + 1
    total_num_per_intervals = ((num_extinct + num_extant + 1):1)[total_intervals]

    extinct_intervals = findInterval(times, node_heights$ta[[i]])
    extinct_num_per_intervals = (num_extinct:0)[extinct_intervals]

    # findInterval(times[9774], node_heights$na[[i]])
    # findInterval(times[9774], node_heights$ta[[i]])
    # times[8779]

    num_lineages[[i]] = total_num_per_intervals - extinct_num_per_intervals - 1

    if(verbose) setTxtProgressBar(bar, i / num_trees)

  }
  num_lineages = do.call(rbind, num_lineages)

  # compute the mean and confidence interval
  mean      = colMeans(num_lineages)
  median    = apply(num_lineages, 2, median)
  alpha     = (1 - CI) / 2
  probs     = c(alpha, 1 - alpha)
  quantiles = apply(num_lineages, 2, quantile, prob=probs)

  # compute the ESS!
  diffs = colSums(abs((t(num_lineages) - mean))) / num_bins
  ess   = effectiveSize(diffs)

  res = list(mean=mean, median=median, quantiles=quantiles, times=times, ess=ess)
  return(res)

}