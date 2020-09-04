library(smacof)
library(vegan)

treeMDS <- setRefClass(

  Class = "treeMDS",

  fields = c(

    ####################
    # basic properties #
    ####################

    "samples",
    "samples_per_chain",
    "max_num_samples",
    "samples_pre_burnin",
    "total_num_samples",
    "names",
    "burnin",
    "mds_method",

    ##################
    # plotting stuff #
    ##################

    "plot_colors",
    "plot_pch",
    "plot_coords",

    ###################
    # the RF distance #
    ###################

    "RF_dist_matrix_clean",
    "RF_dist_matrix",
    "RF_mds",

    ############################
    # the weighted RF distance #
    ############################

    "wRF_dist_matrix_clean",
    "wRF_dist_matrix",
    "wRF_mds",

    ###################
    # the KF distance #
    ###################

    "KF_dist_matrix_clean",
    "KF_dist_matrix",
    "KF_mds"


  ),

  methods = list(

    ###################
    ### Initializer ###
    ###################

    initialize = function(samples,
                          num_samples = 1000,
                          burn_in     = 0,
                          mds_method  = "smacof")
    {

      if ( is.null(names(samples)) ) {
        names <<- paste0("run_", 1:length(samples))
      } else {
        names <<- names(samples)
      }

      mds_method <<- mds_method

      #######################
      # discard some burnin #
      #######################

      max_num_samples <<- num_samples
      burnin <<- burn_in
      samples_pre_burnin <<- samples
      # .self$setBurnin(burn_in)
      .self$processSamples()

    },

    setBurnin = function(b) {

      if ( b != burnin ) {
        burnin <<- b
        .self$processSamples()
      }

    },

    setNumSamples = function(n) {

      if ( n != max_num_samples ) {
        max_num_samples <<- n
        .self$processSamples()
      }

    },

    processSamples = function() {

      ###################
      # toss the burnin #
      ###################

      post_burnin <- vector("list", length(samples_pre_burnin))
      for(i in 1:length(post_burnin)) {
        these_samples <- samples_pre_burnin[[i]]
        if(burnin > 0) {
          these_samples <- these_samples[-c(1:floor(length(these_samples) * burnin))]
        }
        post_burnin[[i]] <- these_samples
      }

      ##################
      # get subsamples #
      ##################

      sub_samples <- vector("list", length(post_burnin))
      for(i in 1:length(post_burnin)) {
        these_samples <- post_burnin[[i]]
        if (length(these_samples) > max_num_samples) {
          these_samples <- these_samples[round(seq(1,length(these_samples),length.out=max_num_samples))]
        }
        sub_samples[[i]] <- these_samples
      }

      ####################
      # save the samples #
      ####################

      samples_per_chain <<- lengths(sub_samples)
      total_num_samples <<- sum(samples_per_chain)
      samples           <<- do.call(c, sub_samples)
      class(samples)    <<- "multiPhylo"

      ######################
      # set flags to false #
      ######################

      RF_dist_matrix_clean  <<- FALSE
      wRF_dist_matrix_clean <<- FALSE
      KF_dist_matrix_clean  <<- FALSE

    },

    computeRFDist = function(dim=2, log=FALSE, method="ratio") {

      if(RF_dist_matrix_clean == FALSE) {

        # compute the distance matrix
        RF_dist_matrix <<- RF.dist(samples, rooted=is.rooted(samples[[1]]))

        # if ( method == "ordinal" ) {
        #   # get the upper diagonal
        #   k <- as.matrix(RF_dist_matrix)
        #   z <- k[upper.tri(k)]
        #   o <- order(z)
        #   a <- numeric(length(o))
        #   a[o] = 1:length(z)
        #   tmp_matrix <- matrix(0, nrow(k), ncol(k))
        #   diag(tmp_matrix) <- 0
        #   tmp_matrix[upper.tri(tmp_matrix)] <- a + 0
        #   tmp_matrix[lower.tri(tmp_matrix)] <- a + 0
        #   RF_dist_matrix <<- as.dist(tmp_matrix)
        # }
        
        if ( log == TRUE ) {
          # RF_dist_matrix <<- log(RF_dist_matrix)
          RF_dist_matrix <<- sqrt(RF_dist_matrix)
        }
        
        # compute the MDS
        if( mds_method == "smacof" & method != "ordinal" ) {
          init_matrix <- cmdscale(RF_dist_matrix, dim)
          RF_mds <<- smacofSym(RF_dist_matrix, dim, init=init_matrix, verbose=TRUE, type=method)$conf
          # RF_mds <<- smacofSphere(RF_dist_matrix, dim, init=init_matrix, itmax=2000, verbose=TRUE)$conf
        } else if ( method == "ordinal" ) {
          init_matrix <- cmdscale(RF_dist_matrix, dim)
          RF_mds <<- metaMDS(RF_dist_matrix, model="hybrid", zerodist="add", autotransform=TRUE, maxit=1000)$points
        } else if ( mds_method == "cmdscale" ) {
          RF_mds <<- cmdscale(RF_dist_matrix, dim)
        } else if ( mds_method == "isoMDS" ) {
          RF_mds <<- isoMDS(as.matrix(RF_dist_matrix) + 1e-16, k=dim, maxit=200)$points
        }

        # set the metric as clean
        RF_dist_matrix_clean <<- TRUE

      }

    },

    computeWeightedRFDist = function(dim=2, log=FALSE, method="ratio") {

      if(wRF_dist_matrix_clean == FALSE) {

        # compute the distance matrix
        wRF_dist_matrix <<- wRF.dist(samples, rooted=is.rooted(samples[[1]]))

        if ( log == TRUE ) {
          # wRF_dist_matrix <<- log(wRF_dist_matrix)
          wRF_dist_matrix <<- sqrt(wRF_dist_matrix)
        }
        
        # compute the MDS
        if( mds_method == "smacof" ) {
          init_matrix <- cmdscale(wRF_dist_matrix, dim)
          wRF_mds <<- smacofSym(wRF_dist_matrix, dim, init=init_matrix, itmax=2000, verbose=TRUE, type=method)$conf
        } else if ( mds_method == "cmdscale" ) {
          wRF_mds <<- cmdscale(wRF_dist_matrix, dim)
        } else if ( mds_method == "isoMDS" ) {
          wRF_mds <<- isoMDS(wRF_dist_matrix, k=dim, trace=FALSE, maxit=200)$points
        }

        # set the metric as clean
        wRF_dist_matrix_clean <<- TRUE

      }

    },

    computeKFDist = function(dim=2, log=FALSE, method="ratio") {

      if(KF_dist_matrix_clean == FALSE) {

        # compute the distance matrix
        KF_dist_matrix <<- KF.dist(samples, rooted=is.rooted(samples[[1]]))

        # if ( method == "ordinal" ) {
        #   # get the upper diagonal
        #   k <- as.matrix(KF_dist_matrix)
        #   z <- k[upper.tri(k)]
        #   o <- order(z)
        #   a <- numeric(length(o))
        #   a[o] = 1:length(z)
        #   tmp_matrix <- matrix(0, nrow(k), ncol(k))
        #   diag(tmp_matrix) <- 0
        #   tmp_matrix[upper.tri(tmp_matrix)] <- a + 0
        #   tmp_matrix[lower.tri(tmp_matrix)] <- a + 0
        #   KF_dist_matrix <<- as.dist(tmp_matrix)
        # }

        if ( log == TRUE ) {
          # KF_dist_matrix <<- log(KF_dist_matrix)
          KF_dist_matrix <<- sqrt(KF_dist_matrix)
        }
        
        # plot(monoMDS(KF_dist_matrix, model="global"))
        
        # compute the MDS
        if( mds_method == "smacof" & method != "ordinal" ) {
          init_matrix <- cmdscale(KF_dist_matrix, dim)
          KF_mds <<- smacofSym(KF_dist_matrix, dim, init=init_matrix, verbose=TRUE, type=method)$conf
        } else if ( method == "ordinal" ) {
          init_matrix <- cmdscale(KF_dist_matrix, dim)
          KF_mds <<- metaMDS(KF_dist_matrix, model="hybrid", zerodist="add", autotransform=TRUE, maxit=1000)$points
        } else if ( mds_method == "cmdscale" ) {
          KF_mds <<- cmdscale(KF_dist_matrix, dim)
        } else if ( mds_method == "isoMDS" ) {
          KF_mds <<- isoMDS(KF_dist_matrix, k=dim, maxit=200)$points
        }

        # set the metric as clean
        KF_dist_matrix_clean <<- TRUE

      }

    },

    getCoordinates = function(type="RF", log=FALSE, method="ratio") {

      if ( type == "RF" ) {
        computeRFDist(log=log, method=method)
        xy = RF_mds
      } else if ( type == "wRF" ) {
        computeWeightedRFDist(log=log, method=method)
        xy = wRF_mds
      } else if ( type == "KF" ) {
        computeKFDist(log=log, method=method)
        xy = KF_mds
      } else {
        stop("Invalid type option. Choose RF, wRF, or KF.")
      }

      return(xy)

    },

    plot = function(type="RF", method="ratio", colors, cex, pch, main=NA, plot_lines=TRUE, xlab=NA, ylab=NA, log=FALSE, ...) {

      # get color palette
      if (missing(colors)) {
        if(length(samples_per_chain) <= 9) {
          colors <- brewer.pal(pmax(3,length(samples_per_chain)), "Set1")[1:length(samples_per_chain)]
        } else {
          colors <- rev(rainbow(length(samples_per_chain), end=0.7))
        }
      }
      plot_colors <<- colors

      # get point types
      if (missing(pch)) {
        pch <- 1:length(samples_per_chain)
      }
      plot_pch <<- pch

      # get the cex
      if ( missing(cex) ) {
        cex <- rep(1, length(samples_per_chain))
      }
      plot_cex <- cex
      if ( length(plot_cex) == 1 ) {
        plot_cex <- rep(plot_cex, length(samples_per_chain))
      }

      # get the coordinates
      xy = getCoordinates(type, log, method)
      plot_coords <<- xy

      # plot the coordinates
      plot.default(xy, type="n", main=main, xlab=xlab, ylab=ylab, ...)

      if (plot_lines) {
        n = 0
        for(i in 1:length(samples_per_chain)) {
          ns = samples_per_chain[i]
          lines(plot_coords[1:ns + n,], col=colors[i], ...)
          n = n + ns
        }
      }

      # randomize the points
      order = sample(nrow(plot_coords))
      sample_colors = rep(plot_colors, times=samples_per_chain)
      sample_pch    = rep(plot_pch,    times=samples_per_chain)
      sample_cex    = rep(plot_cex,    times=samples_per_chain)
      # points(plot_coords[order,], cex=sample_cex[order] + 0.1, col="black", pch=sample_pch[order], ...)
      points(plot_coords[order,], cex=sample_cex[order], col=sample_colors[order], pch=sample_pch[order], ...)

    },

    addLegend = function(x, bty="n", ...) {
      legend(x, legend=names, bty=bty, col=plot_colors, pch=plot_pch, ...)
    },

    getClosestTree = function(coords) {

      closest_tree_index = which.min(sqrt((plot_coords[,1] - coords[1])^2 + (plot_coords[,2] - coords[2])^2))
      tree = samples[[closest_tree_index]]

      return(tree)

    },

    getTree = function(n=1, index=FALSE, ...) {

      # identify the points
      these_points = identify(plot_coords, n=n, plot=FALSE)

      # get the samples for these points
      if(index == TRUE) {
        these_samples = these_points
      } else {
        these_samples = samples[these_points]
      }

      return(these_samples)

    },

    plotAndIdentify = function(type="RF", colors, pch, plot_lines=TRUE, xlab=NA, ylab=NA,
                               show.tip.label=TRUE, tree_type="unrooted",
                               select_color="gold", select_lwd=2, ...) {

      layout_mat = matrix(2:1, nrow=1)
      layout(layout_mat)
      plot.default(1:10, type="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA, bty="n")
      plot(type=type, colors=colors, pch=pch, plot_lines=plot_lines, xlab=xlab, ylab=ylab, ...)
      repeat {
        cat("\nChoose a point.")
        this_point   = getTree(index=TRUE)
        this_tree    = samples[[this_point]]
        point_coords = plot_coords[this_point,]
        from_chain   = findInterval(this_point, cumsum(samples_per_chain)) + 1
        plot.phylo(this_tree, type=tree_type, show.tip.label=show.tip.label)
        plot(type=type, colors=colors, pch=pch, plot_lines=plot_lines, xlab=xlab, ylab=ylab, ...)
        points(x=point_coords[1], y=point_coords[2], pch=plot_pch[from_chain], col=select_color, lwd=select_lwd)
      }

    }

  )

)
