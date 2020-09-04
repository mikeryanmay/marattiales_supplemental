epochLTT = function(sim, age, num_bins=10001, ...) {
  
  # get the points to compute times at
  times       = seq(0, age, length.out=num_bins)
  event_times = rev(sim$event_times)
  num_taxa    = rev(sim$num_taxa_per_time)
  
  # compute the number of taxa per time
  num_taxa_per_time = num_taxa[findInterval(times, event_times) + 1]
  num_taxa_per_time[is.na(num_taxa_per_time)] = 0
    
  return(num_taxa_per_time)
  
}

harmonicMeanPValue = function(obs, sims) {
  
  # compute the number of epochs
  num_epochs  = length(obs)
  num_samples = nrow(sims)
  
  # compute the two-sided p-value for each interval
  p_values = numeric(num_epochs)
  for(i in 1:num_epochs) {
    this_obs     = round(obs[i])
    this_table   = table(sims[,i]) / num_samples
    these_values = as.numeric(names(this_table))
    if ( this_obs %in% these_values == FALSE ) {
      p_values[i] = 0.0
    } else {
      this_prob    = this_table[these_values == this_obs]
      p_values[i]  = sum(this_table[this_table <= this_prob])
    }
  }
  
  # compute the one-sided p-value
  # p_values = rowMeans(t(sims) > obs)
  
  # calculate the harmonic mean
  w = 1 / num_epochs
  harmonic_p = 1 / sum(w / p_values)
  
  # if (harmonic_p == 0) {
  #     recover()
  # }
  
  return(harmonic_p)
  
}

prodPValue = function(obs, sims) {
  
  # compute the number of epochs
  num_epochs = length(obs)
  
  # compute the p-value for each interval
  p_values = rowMeans(t(sims) > obs)
  
  # calculate the product
  prod_p = prod(p_values)
  
  return(prod_p)
  
}

calculateObservedEpochNums = function(data, epoch_times) {
  
  # make the container
  num_taxa_per_epoch = numeric(length(epoch_times))
  
  # split extant and extinct
  extant_data  = data[data$min == 0,]
  extinct_data = data[data$min > 0,]
  
  # loop over extinct taxa
  for(i in 1:nrow(extinct_data)) {
    
    # get this taxon
    this_taxon = extinct_data[i,]
    
    # compute the total duration
    duration = this_taxon$max - this_taxon$min
    
    # compute the start and end intervals
    start_interval = findInterval(this_taxon$min, epoch_times)
    end_interval   = findInterval(this_taxon$max, epoch_times)
    
    if ( start_interval != end_interval ) {
      
      for(j in start_interval:end_interval) {
        
        if (j == 0) {
          x0 = 0
        } else {
          x0 = epoch_times[j]
        }
        x0 = pmax(x0, this_taxon$min)
        x1 = epoch_times[j + 1]
        x1 = pmin(x1, this_taxon$max)
        x  = (x1 - x0) / duration
        num_taxa_per_epoch[j + 1] = num_taxa_per_epoch[j + 1] + x
        
      }
      
    } else {
      num_taxa_per_epoch[start_interval + 1] = num_taxa_per_epoch[start_interval + 1] + 1
    }
    
  }
  
  # add extant data
  num_taxa_per_epoch = c(nrow(extant_data), num_taxa_per_epoch)
  
  return(num_taxa_per_epoch)
  
}

makeEpochRates = function(samples, epoch_times) {
  
  num_epochs = length(epoch_times)
  
  # speciation rates
  if ( "speciation_rate[1]" %in% colnames(samples) == TRUE ) {
    # variable rate
    speciation_rate = samples[paste0("speciation_rate[",(num_epochs+1):1,"]")]
    speciation_rate = lapply(1:nrow(speciation_rate), function(x) return(as.numeric(speciation_rate[x,])))
  } else {
    # constant rate
    speciation_rate = lapply(samples$speciation_rate, rep, times=nrow(epochs) + 1)
  }
  
  # extinction rates
  if ( "extinction_rate[1]" %in% colnames(samples) == TRUE ) {
    # variable rate
    extinction_rate = samples[paste0("extinction_rate[",(num_epochs+1):1,"]")]
    extinction_rate = lapply(1:nrow(extinction_rate), function(x) return(as.numeric(extinction_rate[x,])))
  } else {
    # constant rate
    extinction_rate = lapply(samples$extinction_rate, rep, times=nrow(epochs) + 1)
  }
  
  # fossilization rates
  if ( "fossilization_rate[1]" %in% colnames(samples) == TRUE ) {
    # variable rate
    fossilization_rate = samples[paste0("fossilization_rate[",(num_epochs+1):1,"]")]
    fossilization_rate = lapply(1:nrow(fossilization_rate), function(x) return(as.numeric(fossilization_rate[x,])))
  } else {
    # constant rate
    fossilization_rate = lapply(samples$fossilization_rate, rep, times=nrow(epochs) + 1)
  }
  
  # bundle the rates
  res = list(speciation_rate    = speciation_rate,
             extinction_rate    = extinction_rate,
             fossilization_rate = fossilization_rate)
  
  return(res)
  
}

simulateEpochsConditional = function(origin, epoch_times, lambda, mu, phi, rho,
                                     target = NULL) {
  
  i = 1
  repeat {
    # sim = simulateEpochs(origin, epoch_times, lambda, mu, phi, rho)
    sim = simulateEpochsCPP(origin, epoch_times, lambda, mu, phi, rho)
    if ( is.null(target) ) {
      if ( sim$num_taxa_per_epoch[1] > 0 ) {
        break
      }
    } else {
      if ( all(sim$num_taxa_per_epoch == target) ) {
        break
      }
    }
    i = i + 1
    # cat(i, ":\t", sim, "\n", sep="")
    if ( i > 100000 ) {
      return(NULL)
    }
  }
  
  return(sim)
  
  
}

simulateEpochs = function(origin, epoch_times, lambda, mu, phi, rho) {
  
  # initialize the containers
  current_num_taxa   = 1
  num_taxa_per_epoch = numeric(length(epoch_times))
  
  # get the current time
  current_time = pmin(max(epoch_times), origin)
  # current_time = origin
  
  # get the current epoch index
  current_epoch   = findInterval(current_time, epoch_times) + 1
  next_epoch      = current_epoch - 1
  next_epoch_time = epoch_times[next_epoch]
  
  # get the rates
  current_lambda = lambda[current_epoch]
  current_mu     = mu[current_epoch]
  current_phi    = phi[current_epoch]
  
  # forward simulate
  repeat {
    
    # compute the rate
    total_rate = current_lambda + current_mu + current_phi
    
    # draw the next time
    next_time = current_time - rexp(1, current_num_taxa * total_rate)
    
    if ( next_time < next_epoch_time ) {
      
      # terminate if the next epoch is the present
      if ( next_epoch_time == 0 ) {
        break
      }
      
      # we moved into the next epoch
      current_time    = next_epoch_time
      current_epoch   = next_epoch
      next_epoch      = current_epoch - 1
      next_epoch_time = epoch_times[next_epoch]
      
      # compute the rates
      current_lambda = lambda[current_epoch]
      current_mu     = mu[current_epoch]
      current_phi    = phi[current_epoch]
      
      if ( next_epoch == 0 ) {
        next_epoch_time = 0
      }
      
    } else {
      
      # increment time
      current_time = next_time
      
      # simulate an event
      next_event_type = sample.int(3, 1, prob = c(current_lambda, current_mu, current_phi))
      
      if ( next_event_type == 1 ) {
        # speciation
        current_num_taxa = current_num_taxa + 1
      } else if ( next_event_type == 2 ) {
        # extinction
        current_num_taxa = current_num_taxa - 1
      } else {
        # fossilization
        num_taxa_per_epoch[current_epoch] = num_taxa_per_epoch[current_epoch] + 1
      }
      
    }
    
    # cat(current_time, current_num_taxa, "\n")
    
    # check for extinction
    if (current_num_taxa == 0) {
      break
    }
    
    # if ( current_num_taxa * rho > 1000 ) {
    #     current_num_taxa = 0
    #     break
    # }
    
  }
  
  # add the extant taxa
  current_num_taxa = rbinom(1, current_num_taxa, rho)
  # cat(mean(rbinom(1000, 1000, rho)))
  num_taxa_per_epoch = c(current_num_taxa, num_taxa_per_epoch)
  
  return(num_taxa_per_epoch)
  
}