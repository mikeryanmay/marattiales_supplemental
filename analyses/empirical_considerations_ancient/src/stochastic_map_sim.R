# get the arguments
args = commandArgs(trailingOnly = TRUE)

# args = c("--data", "data",
#          "--samples", "output/tree_model_variable_rate_fbd_mixture_diversification_ingroup_sub_model_GTR_I_G_mol_clock_UCLN_morph_mat_F81_MIX_G_morph_rel_unlinked_morph_clock_linked/",
#          "--output", "output_stoch_map/tree_model_variable_rate_fbd_mixture_diversification_ingroup_sub_model_GTR_I_G_mol_clock_UCLN_morph_mat_F81_MIX_G_morph_rel_unlinked_morph_clock_linked/",
#          "--nsim", 1000)

# get the data directory
if ( "--data" %in% args ) {
  data_dir = args[which(args == "--data") + 1]
} else {
  stop("Must provide an --data argument!")
}

# get the sample directory
if ( "--samples" %in% args ) {
  sample_dir = args[which(args == "--samples") + 1]
} else {
  stop("Must provide an --sample argument!")
}

# get the output directory
if ( "--output" %in% args ) {
  output_dir = args[which(args == "--output") + 1]
} else {
  stop("Must provide an --output argument!")
}

# the number of simulations
nsim = 1000
if ( "--nsim" %in% args ) {
  nsim = as.numeric(args[which(args == "--nsim") + 1])
}

ncores = 1
if ( "--ncores" %in% args ) {
  ncores = as.numeric(args[which(args == "--ncores") + 1])
}

# get the overwrite
overwrite = FALSE
if ( "--overwrite" %in% args ) {
  if ( tolower(args[which(args == "--overwrite") + 1]) == "true" ) {
    overwrite = TRUE
  } else if ( tolower(args[which(args == "--overwrite") + 1]) == "false" ) {
    overwrite = FALSE
  } else {
    stop("Invalided --overwrite value!")
  }
}

cat("Stochastic mapping for:\n  ", sample_dir,"\n", sep="")

# check for tree files
tree_fn = paste0(sample_dir, "/morph_phylogram_combined_MCC.tre")
tt_fn   = paste0(sample_dir, "/tree_combined_MCC.tre")

if ( file.exists(tree_fn) == FALSE | file.exists(tt_fn) == FALSE ) {
  cat("Tree files do not exist.")
  q()
}

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

if ( "rncl" %in% rownames(installed.packages()) == FALSE ) {
  install.packages("rncl")
}
library(rncl)

if ( "stringr" %in% rownames(installed.packages()) == FALSE ) {
  install.packages("stringr")
}

if ( "parallel" %in% rownames(installed.packages()) == FALSE ) {
  install.packages("parallel")
}
library(parallel)

if ( "plotrix" %in% rownames(installed.packages()) == FALSE ) {
  install.packages("plotrix")
}
library(plotrix)

if ( "RColorBrewer" %in% rownames(installed.packages()) == FALSE ) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999", "#D55E00", "#0072B2", "#CC79A7")

# source the appropriate scripts
source("src/pps_functions.R")
source("src/plot_simmap.R")

rename = c(
  "Marattia_asiatica"           = "Marattiopsis_asiatica",
  "Marattia_anglica"            = "Marattiopsis_anglica",
  "Marattia_aganzhenensis"      = "Marattiopsis_aganzhenensis",
  "Scolecopteris_alta_A"        = "Scolecopteris_alta",
  "Scolecopteris_antarctica_L"  = "Scolecopteris_antarctica",
  "Scolecopteris_calicifolia_L" = "Scolecopteris_calicifolia",
  "Scolecopteris_charma_O"      = "Scolecopteris_charma",
  "Scolecopteris_fragilis_L"    = "Scolecopteris_fragilis",
  "Scolecopteris_incisifolia_L" = "Scolecopteris_incisifolia",
  "Scolecopteris_iowensis_O"    = "Scolecopteris_iowensis",
  "Scolecopteris_latifolia_L"   = "Scolecopteris_latifolia",
  "Scolecopteris_majopsis_O"    = "Scolecopteris_majopsis",
  "Scolecopteris_mamayi_L"      = "Scolecopteris_mamayi",
  "Scolecopteris_minor_M"       = "Scolecopteris_minor",
  "Scolecopteris_monothrix_L"   = "Scolecopteris_monothrix",
  "Scolecopteris_nigra_A"       = "Scolecopteris_nigra",
  "Scolecopteris_oliveri_O"     = "Scolecopteris_oliveri",
  "Scolecopteris_parkerensis_L" = "Scolecopteris_parkerensis",
  "Scolecopteris_saharaensis_M" = "Scolecopteris_saharaensis",
  "Scolecopteris_vallumii_L"    = "Scolecopteris_vallumii"
)

# read tree and samples
tree      = make_tree_bifurcating(read.nexus(tree_fn))
time_tree = make_tree_bifurcating(read.nexus(tt_fn))

# drop sampled ancestors
time_tree = drop.tip(time_tree, tree$tip.label[tree$edge[tree$edge.length == 0,2]])
tree      = drop.tip(tree, tree$tip.label[tree$edge[tree$edge.length == 0,2]])

# rename the time tree
for(i in 1:length(rename)) {
  this_rename = rename[i]
  time_tree$tip.label[time_tree$tip.label == names(this_rename)] = this_rename
}

# read observed data
obs_data = readDataAsProbs(paste0(data_dir, "/morpho.nex"))

# read the character information
char_data = read.csv(paste0(data_dir, "/char_table.tsv"), check.names=FALSE, sep="\t", stringsAsFactors=FALSE)

# read the parameters
param_files  = paste0(paste0(sample_dir, "/params_combined.log"))
samples      = do.call(rbind, lapply(param_files, read.table, header=TRUE, sep="\t", stringsAsFactors=TRUE, check.names=FALSE))
samples_mean = colMeans(samples)

# compute the number of partitions
num_partitions = length(obs_data)
num_states     = as.numeric(names(obs_data))

# create the output directory
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

mclapply(1:num_partitions, function(i){
  
  # simulate maps
  sims = simulate_stochastic_map(obs_data[[i]], tree, samples_mean, part=i, num_states=num_states[i], nsims=nsim)
  
  # get char info
  this_char_data = char_data[char_data$num_states == num_states[i],]
  
  # make colors
  labels = as.character(1:num_states[i]-1)
  # colors = brewer.pal(9,"Set1")[1:num_states[i]]
  colors = cbPalette[1:num_states[i]]
  
  for(j in 1:nrow(this_char_data)) {
    
    # get this character
    this_char_name  = this_char_data$names[j]
    this_char_title = gsub("_"," ", this_char_name)
    this_char_id    = this_char_data$index[j]
    these_sims      = sims[[j]]
    state_labels    = strsplit(this_char_data$states[j],",")[[1]]
    
    # plot this character
    this_fig_name = paste0(output_dir,"/stoch_map_", this_char_id,"_",this_char_name,".pdf")
    
    pdf(this_fig_name)
    par(mar=c(2,0,1,0), lend=2)
    plot_simmap(time_tree, tree, obs_data[[i]][,j,], these_sims, labels, colors=colors, nt=2001, show.tip.label=TRUE, lwd=2, edge.width=3, lend=2, pie_size=3.5, label.offset=10, label.cex=0.5)
    legend("bottomleft", legend=gsub("_"," ",state_labels), fill=colors, bty="n")
    mtext(this_char_title)
    axisPhylo()
    dev.off()
    
  }
  
}, mc.cores=ncores, mc.preschedule=FALSE)

# for(i in 1:num_partitions) {
#   
#   # simulate maps
#   sims = simulate_stochastic_map(obs_data[[i]], tree, samples_mean, part=i, num_states=num_states[i], nsims=nsim)
#   
#   # get char info
#   this_char_data = char_data[char_data$num_states == num_states[i],]
#   
#   # make colors
#   labels = as.character(1:num_states[i]-1)
#   colors = brewer.pal(9,"Set1")[1:num_states[i]]
#   
#   for(j in 1:nrow(this_char_data)) {
#     
#     # get this character
#     this_char_name  = this_char_data$names[j]
#     this_char_title = gsub("_"," ", this_char_name)
#     this_char_id    = this_char_data$index[j]
#     these_sims      = sims[[j]]
#     state_labels    = strsplit(this_char_data$states[j],",")[[1]]
#     
#     # plot this character
#     this_fig_name = paste0(output_dir,"/stoch_map_", this_char_id,"_",this_char_name,".pdf")
#     
#     pdf(this_fig_name)
#     par(mar=c(2,0,1,0), lend=2)
#     plot_simmap(time_tree, tree, obs_data[[i]][,j,], these_sims, labels, colors=colors, nt=2001, show.tip.label=TRUE, lwd=2, edge.width=3, lend=2, pie_size=3.5, label.offset=10, label.cex=0.5)
#     legend("bottomleft", legend=gsub("_"," ",state_labels), fill=colors, bty="n")
#     mtext(this_char_title)
#     axisPhylo()
#     dev.off()
#     
#   }
#   
# }
