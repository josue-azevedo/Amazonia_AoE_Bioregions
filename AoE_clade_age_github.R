library(ape)
library(phytools)
library(ggplot2)


#fully_sampled_rep = read.tree("phylogenies/squam_shl_new_Posterior_9755.1000.trees")

fully_sampled_rep[[1]]

#species_by_aoe = read.csv("amazniarecorte/species_tonini_AoE.csv")
species_by_aoe = read.csv("species_tonini.csv")
head(species_by_aoe)
str(species_by_aoe)
summary(species_by_aoe)


library(ape)
library(phytools)
library(picante) # For phylogenetic structure metrics

# Step 1: Select a tree to work with
tree <- fully_sampled_rep[[1]]

# Step 2: Get list of species by AoE
aoe_species_list <- split(species_by_aoe$speciesnames, species_by_aoe$AoE)

# Step 3: Create community matrix for all AoEs (for phylogenetic structure calculation)
all_species <- unique(unlist(aoe_species_list))
all_species <- all_species[all_species %in% tree$tip.label]

comm_matrix <- matrix(0, ncol=length(aoe_species_list), nrow=length(all_species))
rownames(comm_matrix) <- all_species
colnames(comm_matrix) <- names(aoe_species_list)

# Fill community matrix
for(aoe in names(aoe_species_list)) {
  species <- aoe_species_list[[aoe]][aoe_species_list[[aoe]] %in% tree$tip.label]
  if(length(species) > 0) {
    comm_matrix[species, aoe] <- 1
  }
}

# Verify community matrix structure
cat("Community matrix dimensions:", dim(comm_matrix), "\n")
cat("First few row sums (species):", head(rowSums(comm_matrix)), "\n")
cat("First few column sums (AoEs):", head(colSums(comm_matrix)), "\n")

# Step 4: Calculate cophenetic distance matrix
phy_dist <- cophenetic.phylo(tree)
cat("Phylogenetic distance matrix dimensions:", dim(phy_dist), "\n")

# Step 5: Create results dataframe
results <- data.frame(
  AoE = character(),
  n_species = integer(),
  is_monophyletic = logical(),
  mrca_age = numeric(),
  mean_terminal_length = numeric(),
  mpd_obs = numeric(),
  mpd_rand_mean = numeric(),
  mpd_rand_sd = numeric(),
  mpd_z_score = numeric(),
  mntd_obs = numeric(),
  mntd_rand_mean = numeric(),
  mntd_rand_sd = numeric(),
  mntd_z_score = numeric(),
  phylo_structure = character(),
  diversification_metric = numeric(),
  stringsAsFactors = FALSE
)

# Step 6: Process each AoE individually
# Let's start with just a few AoEs for debugging
test_aoes <- head(names(aoe_species_list), 3)
cat("Testing with the following AoEs:", paste(test_aoes, collapse=", "), "\n\n")

for(aoe in test_aoes) {
  cat("=== Processing AoE:", aoe, "===\n")
  species <- aoe_species_list[[aoe]][aoe_species_list[[aoe]] %in% tree$tip.label]
  n_species <- length(species)
  
  cat("Number of species in tree:", n_species, "\n")
  if(n_species > 0) {
    # Basic metrics
    is_mono <- NA
    if(n_species > 1) {
      is_mono <- is.monophyletic(tree, species)
      cat("Is monophyletic:", is_mono, "\n")
    }
    
    # Calculate MRCA age
    if(n_species > 1) {
      mrca_node <- findMRCA(tree, species)
      cat("MRCA node:", mrca_node, "\n")
      mrca_age <- nodeheight(tree, mrca_node)
      cat("MRCA age:", mrca_age, "\n")
    } else if(n_species == 1) {
      # For singletons, use the tip's parent node age
      tip_idx <- which(tree$tip.label == species[1])
      parent_node <- tree$edge[which(tree$edge[,2] == tip_idx), 1]
      mrca_age <- nodeheight(tree, parent_node)
      cat("Single species with parent node age:", mrca_age, "\n")
    } else {
      mrca_age <- NA
    }
    
    # Calculate terminal branch lengths (true tip ages)
    if(n_species > 0) {
      terminal_branches <- numeric(n_species)
      for(i in 1:n_species) {
        tip_idx <- which(tree$tip.label == species[i])
        edge_idx <- which(tree$edge[,2] == tip_idx)
        terminal_branches[i] <- tree$edge.length[edge_idx]
      }
      cat("Terminal branch lengths:", terminal_branches, "\n")
      mean_terminal_length <- mean(terminal_branches)
      cat("Mean terminal branch length:", mean_terminal_length, "\n")
    } else {
      mean_terminal_length <- NA
    }
    
    # Calculate phylogenetic metrics for AoEs with >1 species
    mpd_obs <- NA
    mpd_rand_mean <- NA
    mpd_rand_sd <- NA
    mpd_z <- NA
    mntd_obs <- NA
    mntd_rand_mean <- NA
    mntd_rand_sd <- NA
    mntd_z <- NA
    phylo_structure <- "NA"
    
    if(n_species > 1) {
      # Test if we have enough species for meaningful analysis
      if(n_species >= 3) {
        cat("Calculating phylogenetic structure metrics...\n")
        
        # Get only the distance matrix entries for our species
        species_dist <- phy_dist[species, species]
        cat("Pairwise distance matrix for this AoE size:", dim(species_dist), "\n")
        
        # Calculate MPD (Mean Pairwise Distance) manually
        species_pairs <- combn(n_species, 2)
        pair_distances <- numeric(ncol(species_pairs))
        for(i in 1:ncol(species_pairs)) {
          sp1 <- species[species_pairs[1,i]]
          sp2 <- species[species_pairs[2,i]]
          pair_distances[i] <- phy_dist[sp1, sp2]
        }
        mpd_obs <- mean(pair_distances)
        cat("Observed MPD:", mpd_obs, "\n")
        
        # Calculate MNTD (Mean Nearest Taxon Distance) manually
        nearest_distances <- numeric(n_species)
        for(i in 1:n_species) {
          sp <- species[i]
          other_spp <- species[-i]
          nearest_distances[i] <- min(phy_dist[sp, other_spp])
        }
        mntd_obs <- mean(nearest_distances)
        cat("Observed MNTD:", mntd_obs, "\n")
        
        # Run null model to get z-scores (do fewer runs for debugging)
        n_rand <- 99 # Fewer runs for debugging
        
        # For MPD null distribution
        rand_mpds <- numeric(n_rand)
        for(i in 1:n_rand) {
          # Randomly sample n_species from all species
          rand_species <- sample(all_species, n_species)
          # Calculate MPD for random sample
          rand_dist <- phy_dist[rand_species, rand_species]
          rand_pairs <- combn(n_species, 2)
          rand_pair_dists <- numeric(ncol(rand_pairs))
          for(j in 1:ncol(rand_pairs)) {
            rs1 <- rand_species[rand_pairs[1,j]]
            rs2 <- rand_species[rand_pairs[2,j]]
            rand_pair_dists[j] <- phy_dist[rs1, rs2]
          }
          rand_mpds[i] <- mean(rand_pair_dists)
        }
        
        # Calculate z-score for MPD
        mpd_rand_mean <- mean(rand_mpds)
        mpd_rand_sd <- sd(rand_mpds)
        mpd_z <- (mpd_obs - mpd_rand_mean) / mpd_rand_sd
        cat("MPD random mean:", mpd_rand_mean, "\n")
        cat("MPD random SD:", mpd_rand_sd, "\n")
        cat("MPD Z-score:", mpd_z, "\n")
        
        # For MNTD null distribution (similar process)
        rand_mntds <- numeric(n_rand)
        for(i in 1:n_rand) {
          rand_species <- sample(all_species, n_species)
          rand_nearest_dists <- numeric(n_species)
          for(j in 1:n_species) {
            rs <- rand_species[j]
            other_rs <- rand_species[-j]
            rand_nearest_dists[j] <- min(phy_dist[rs, other_rs])
          }
          rand_mntds[i] <- mean(rand_nearest_dists)
        }
        
        # Calculate z-score for MNTD
        mntd_rand_mean <- mean(rand_mntds)
        mntd_rand_sd <- sd(rand_mntds)
        mntd_z <- (mntd_obs - mntd_rand_mean) / mntd_rand_sd
        cat("MNTD random mean:", mntd_rand_mean, "\n")
        cat("MNTD random SD:", mntd_rand_sd, "\n")
        cat("MNTD Z-score:", mntd_z, "\n")
        
        # Interpret phylogenetic structure (using just a single z-score to simplify)
        if(mpd_z < -1.96) {
          phylo_structure <- "clustered_deep"
          cat("Phylogenetic structure: CLUSTERED at deep nodes\n")
        } else if(mpd_z > 1.96) {
          phylo_structure <- "overdispersed_deep"
          cat("Phylogenetic structure: OVERDISPERSED at deep nodes\n")
        } else {
          phylo_structure <- "random"
          cat("Phylogenetic structure: RANDOM\n")
        }
      } else {
        cat("Too few species for phylogenetic structure analysis\n")
        phylo_structure <- "too_few_species"
      }
    }
    
    # Calculate diversification metric for all multi-species AoEs
    div_metric <- NA
    if(n_species > 1) {
      # Simple metric: ln(species) / MRCA age
      div_metric <- log(n_species) / mrca_age
      cat("Diversification metric (ln(species)/MRCA age):", div_metric, "\n")
    }
    
    # Create result row
    aoe_row <- data.frame(
      AoE = aoe,
      n_species = n_species,
      is_monophyletic = is_mono,
      mrca_age = mrca_age,
      mean_terminal_length = mean_terminal_length,
      mpd_obs = mpd_obs,
      mpd_rand_mean = mpd_rand_mean,
      mpd_rand_sd = mpd_rand_sd,
      mpd_z_score = mpd_z,
      mntd_obs = mntd_obs,
      mntd_rand_mean = mntd_rand_mean,
      mntd_rand_sd = mntd_rand_sd,
      mntd_z_score = mntd_z,
      phylo_structure = phylo_structure,
      diversification_metric = div_metric,
      stringsAsFactors = FALSE
    )
    
    # Add to results
    results <- rbind(results, aoe_row)
    cat("\n")
  }
}

# Step 7: Display results for the test AoEs
cat("\nResults for test AoEs:\n")
print(results)

# Step 8: Run for all AoEs
cat("\nNow processing all AoEs...\n")
all_results <- data.frame(
  AoE = character(),
  n_species = integer(),
  is_monophyletic = logical(),
  mrca_age = numeric(),
  mean_terminal_length = numeric(),
  mpd_obs = numeric(),
  mpd_rand_mean = numeric(),
  mpd_rand_sd = numeric(),
  mpd_z_score = numeric(),
  mntd_obs = numeric(),
  mntd_rand_mean = numeric(),
  mntd_rand_sd = numeric(),
  mntd_z_score = numeric(),
  phylo_structure = character(),
  diversification_metric = numeric(),
  stringsAsFactors = FALSE
)

# This is just a demonstration - you can comment this out if you want to run only the test AoEs
for(aoe in setdiff(names(aoe_species_list), test_aoes)) {
  cat("Processing AoE:", aoe, "\n")
  species <- aoe_species_list[[aoe]][aoe_species_list[[aoe]] %in% tree$tip.label]
  n_species <- length(species)
  
  if(n_species > 0) {
    # Basic metrics
    is_mono <- NA
    if(n_species > 1) {
      is_mono <- is.monophyletic(tree, species)
    }
    
    # Calculate MRCA age
    if(n_species > 1) {
      mrca_node <- findMRCA(tree, species)
      mrca_age <- nodeheight(tree, mrca_node)
    } else if(n_species == 1) {
      # For singletons, use the tip's parent node age
      tip_idx <- which(tree$tip.label == species[1])
      parent_node <- tree$edge[which(tree$edge[,2] == tip_idx), 1]
      mrca_age <- nodeheight(tree, parent_node)
    } else {
      mrca_age <- NA
    }
    
    # Calculate terminal branch lengths
    if(n_species > 0) {
      terminal_branches <- numeric(n_species)
      for(i in 1:n_species) {
        tip_idx <- which(tree$tip.label == species[i])
        edge_idx <- which(tree$edge[,2] == tip_idx)
        terminal_branches[i] <- tree$edge.length[edge_idx]
      }
      mean_terminal_length <- mean(terminal_branches)
    } else {
      mean_terminal_length <- NA
    }
    
    # Calculate phylogenetic metrics for AoEs with >1 species
    mpd_obs <- NA
    mpd_rand_mean <- NA
    mpd_rand_sd <- NA
    mpd_z <- NA
    mntd_obs <- NA
    mntd_rand_mean <- NA
    mntd_rand_sd <- NA
    mntd_z <- NA
    phylo_structure <- "NA"
    
    # Only calculate for AoEs with at least 3 species
    if(n_species >= 3) {
      # Calculate observed MPD
      species_pairs <- combn(n_species, 2)
      pair_distances <- numeric(ncol(species_pairs))
      for(i in 1:ncol(species_pairs)) {
        sp1 <- species[species_pairs[1,i]]
        sp2 <- species[species_pairs[2,i]]
        pair_distances[i] <- phy_dist[sp1, sp2]
      }
      mpd_obs <- mean(pair_distances)
      
      # Calculate observed MNTD
      nearest_distances <- numeric(n_species)
      for(i in 1:n_species) {
        sp <- species[i]
        other_spp <- species[-i]
        nearest_distances[i] <- min(phy_dist[sp, other_spp])
      }
      mntd_obs <- mean(nearest_distances)
      
      # Run null model (fewer runs for all AoEs)
      n_rand <- 99
      
      # For MPD null distribution
      rand_mpds <- numeric(n_rand)
      for(i in 1:n_rand) {
        rand_species <- sample(all_species, n_species)
        rand_dist <- phy_dist[rand_species, rand_species]
        rand_pairs <- combn(n_species, 2)
        rand_pair_dists <- numeric(ncol(rand_pairs))
        for(j in 1:ncol(rand_pairs)) {
          rs1 <- rand_species[rand_pairs[1,j]]
          rs2 <- rand_species[rand_pairs[2,j]]
          rand_pair_dists[j] <- phy_dist[rs1, rs2]
        }
        rand_mpds[i] <- mean(rand_pair_dists)
      }
      
      # Calculate z-score for MPD
      mpd_rand_mean <- mean(rand_mpds)
      mpd_rand_sd <- sd(rand_mpds)
      mpd_z <- (mpd_obs - mpd_rand_mean) / mpd_rand_sd
      
      # For MNTD null distribution
      rand_mntds <- numeric(n_rand)
      for(i in 1:n_rand) {
        rand_species <- sample(all_species, n_species)
        rand_nearest_dists <- numeric(n_species)
        for(j in 1:n_species) {
          rs <- rand_species[j]
          other_rs <- rand_species[-j]
          rand_nearest_dists[j] <- min(phy_dist[rs, other_rs])
        }
        rand_mntds[i] <- mean(rand_nearest_dists)
      }
      
      # Calculate z-score for MNTD
      mntd_rand_mean <- mean(rand_mntds)
      mntd_rand_sd <- sd(rand_mntds)
      mntd_z <- (mntd_obs - mntd_rand_mean) / mntd_rand_sd
      
      # Interpret phylogenetic structure
      if(mpd_z < -1.96) {
        phylo_structure <- "clustered_deep"
      } else if(mpd_z > 1.96) {
        phylo_structure <- "overdispersed_deep"
      } else {
        phylo_structure <- "random"
      }
    } else if(n_species == 2) {
      phylo_structure <- "too_few_species"
    } else {
      phylo_structure <- "singleton"
    }
    
    # Calculate diversification metric
    div_metric <- NA
    if(n_species > 1) {
      div_metric <- log(n_species) / mrca_age
    }
    
    # Create result row
    aoe_row <- data.frame(
      AoE = aoe,
      n_species = n_species,
      is_monophyletic = is_mono,
      mrca_age = mrca_age,
      mean_terminal_length = mean_terminal_length,
      mpd_obs = mpd_obs,
      mpd_rand_mean = mpd_rand_mean,
      mpd_rand_sd = mpd_rand_sd,
      mpd_z_score = mpd_z,
      mntd_obs = mntd_obs,
      mntd_rand_mean = mntd_rand_mean,
      mntd_rand_sd = mntd_rand_sd,
      mntd_z_score = mntd_z,
      phylo_structure = phylo_structure,
      diversification_metric = div_metric,
      stringsAsFactors = FALSE
    )
    
    # Add to results
    all_results <- rbind(all_results, aoe_row)
  }
}

# Combine test and all results
final_results <- rbind(results, all_results)

# Step 9: Visualize results
cat("\nCreating visualizations...\n")
library(ggplot2)

# Plot 1: Diversification metric by AoE
p1 <- ggplot(subset(final_results, !is.na(diversification_metric)), 
             aes(x=reorder(AoE, diversification_metric), y=diversification_metric)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Species Accumulation Rate by Area of Endemism", 
       x="Area of Endemism", 
       y="ln(species)/MRCA Age")

# Plot 2: Mean terminal branch length by AoE
p2 <- ggplot(final_results, 
             aes(x=reorder(AoE, mean_terminal_length), y=mean_terminal_length)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Mean Terminal Branch Length by Area of Endemism", 
       x="Area of Endemism", 
       y="Mean Terminal Branch Length")

# Plot 3: Phylogenetic structure (MPD Z-score)
p3 <- ggplot(subset(final_results, !is.na(mpd_z_score)), 
             aes(x=reorder(AoE, mpd_z_score), y=mpd_z_score, fill=mpd_z_score<(-1.96) | mpd_z_score>1.96)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=c(-1.96, 1.96), linetype="dashed", color="red") +
  scale_fill_manual(values=c("grey", "darkred"), name="Significant") +
  coord_flip() +
  labs(title="Phylogenetic Structure by Area of Endemism", 
       x="Area of Endemism", 
       y="MPD Z-Score (Negative = Clustered, Positive = Overdispersed)")

# Display the plots
print(p1)
print(p2)
print(p3)

# Step 10: Summarize findings for manuscript
cat("\nSummary of phylogenetic structure patterns:\n")
structure_table <- table(final_results$phylo_structure)
print(structure_table)

cat("\nAoEs with evidence of in-situ diversification (clustering):\n")
clustered_aoes <- final_results$AoE[final_results$phylo_structure == "clustered_deep"]
print(clustered_aoes)

cat("\nAoEs with evidence of multiple colonization events (overdispersion):\n")
overdispersed_aoes <- final_results$AoE[final_results$phylo_structure == "overdispersed_deep"]
print(overdispersed_aoes)

# Create a summary table for manuscript
summary_table <- final_results[, c("AoE", "n_species", "is_monophyletic", 
                                   "mrca_age", "mean_terminal_length", 
                                   "mpd_z_score", "phylo_structure", 
                                   "diversification_metric")]

# Round numeric columns
summary_table[, c("mrca_age", "mean_terminal_length", "mpd_z_score", "diversification_metric")] <- 
  round(summary_table[, c("mrca_age", "mean_terminal_length", "mpd_z_score", "diversification_metric")], 3)

# Sort by phylogenetic structure and number of species
summary_table <- summary_table[order(summary_table$phylo_structure, -summary_table$n_species), ]

cat("\nSummary table for manuscript:\n")
print(summary_table)

write.csv(summary_table,"amazniarecorte/AoE_clade_age.csv")


########## Now for multiple phylogenies











library(ape)
library(phytools)
library(ggplot2)

# Function to analyze a single tree - based on our validated code
analyze_single_tree <- function(tree, aoe_species_list, tree_id) {
  # Get all species that exist in the tree
  all_species <- unique(unlist(aoe_species_list))
  all_species <- all_species[all_species %in% tree$tip.label]
  
  # Calculate cophenetic distance matrix
  phy_dist <- cophenetic.phylo(tree)
  
  # Set up results dataframe
  results <- data.frame(
    tree_id = integer(),
    AoE = character(),
    n_species = integer(),
    is_monophyletic = logical(),
    mrca_age = numeric(),
    mean_terminal_length = numeric(),
    mpd_obs = numeric(),
    mpd_rand_mean = numeric(),
    mpd_rand_sd = numeric(),
    mpd_z_score = numeric(),
    phylo_structure = character(),
    diversification_metric = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each AoE
  for(aoe in names(aoe_species_list)) {
    species <- aoe_species_list[[aoe]][aoe_species_list[[aoe]] %in% tree$tip.label]
    n_species <- length(species)
    
    if(n_species > 0) {
      # Basic metrics
      is_mono <- NA
      if(n_species > 1) {
        is_mono <- is.monophyletic(tree, species)
      }
      
      # Calculate MRCA age
      if(n_species > 1) {
        mrca_node <- findMRCA(tree, species)
        mrca_age <- nodeheight(tree, mrca_node)
      } else if(n_species == 1) {
        # For singletons, use the tip's parent node age
        tip_idx <- which(tree$tip.label == species[1])
        parent_node <- tree$edge[which(tree$edge[,2] == tip_idx), 1]
        mrca_age <- nodeheight(tree, parent_node)
      } else {
        mrca_age <- NA
      }
      
      # Calculate terminal branch lengths
      if(n_species > 0) {
        terminal_branches <- numeric(n_species)
        for(i in 1:n_species) {
          tip_idx <- which(tree$tip.label == species[i])
          edge_idx <- which(tree$edge[,2] == tip_idx)
          terminal_branches[i] <- tree$edge.length[edge_idx]
        }
        mean_terminal_length <- mean(terminal_branches)
      } else {
        mean_terminal_length <- NA
      }
      
      # Calculate phylogenetic metrics for AoEs with >1 species
      mpd_obs <- NA
      mpd_rand_mean <- NA
      mpd_rand_sd <- NA
      mpd_z <- NA
      phylo_structure <- "NA"
      
      # Only calculate for AoEs with at least 3 species
      if(n_species >= 3) {
        # Calculate observed MPD
        species_pairs <- combn(n_species, 2)
        pair_distances <- numeric(ncol(species_pairs))
        for(i in 1:ncol(species_pairs)) {
          sp1 <- species[species_pairs[1,i]]
          sp2 <- species[species_pairs[2,i]]
          pair_distances[i] <- phy_dist[sp1, sp2]
        }
        mpd_obs <- mean(pair_distances)
        
        # Run null model (fewer runs for multiple trees)
        n_rand <- 99
        
        # For MPD null distribution
        rand_mpds <- numeric(n_rand)
        for(i in 1:n_rand) {
          rand_species <- sample(all_species, n_species)
          rand_pairs <- combn(n_species, 2)
          rand_pair_dists <- numeric(ncol(rand_pairs))
          for(j in 1:ncol(rand_pairs)) {
            rs1 <- rand_species[rand_pairs[1,j]]
            rs2 <- rand_species[rand_pairs[2,j]]
            rand_pair_dists[j] <- phy_dist[rs1, rs2]
          }
          rand_mpds[i] <- mean(rand_pair_dists)
        }
        
        # Calculate z-score for MPD
        mpd_rand_mean <- mean(rand_mpds)
        mpd_rand_sd <- sd(rand_mpds)
        mpd_z <- (mpd_obs - mpd_rand_mean) / mpd_rand_sd
        
        # Interpret phylogenetic structure
        if(mpd_z < -1.96) {
          phylo_structure <- "clustered_deep"
        } else if(mpd_z > 1.96) {
          phylo_structure <- "overdispersed_deep"
        } else {
          phylo_structure <- "random"
        }
      } else if(n_species == 2) {
        phylo_structure <- "too_few_species"
      } else {
        phylo_structure <- "singleton"
      }
      
      # Calculate diversification metric
      div_metric <- NA
      if(n_species > 1) {
        div_metric <- log(n_species) / mrca_age
      }
      
      # Create result row
      aoe_row <- data.frame(
        tree_id = tree_id,
        AoE = aoe,
        n_species = n_species,
        is_monophyletic = is_mono,
        mrca_age = mrca_age,
        mean_terminal_length = mean_terminal_length,
        mpd_obs = mpd_obs,
        mpd_rand_mean = mpd_rand_mean,
        mpd_rand_sd = mpd_rand_sd,
        mpd_z_score = mpd_z,
        phylo_structure = phylo_structure,
        diversification_metric = div_metric,
        stringsAsFactors = FALSE
      )
      
      # Add to results
      results <- rbind(results, aoe_row)
    }
  }
  
  return(results)
}

# Function to analyze multiple trees and aggregate results
analyze_multiple_trees <- function(trees, aoe_species_list, n_trees = 100) {
  # Results from all trees
  all_tree_results <- list()
  
  # Process each tree
  for(i in 1:min(n_trees, length(trees))) {
    cat("Processing tree", i, "of", min(n_trees, length(trees)), "\n")
    tree_results <- analyze_single_tree(trees[[i]], aoe_species_list, i)
    all_tree_results[[i]] <- tree_results
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_tree_results)
  
  # Aggregate numeric metrics across trees
  numeric_cols <- c("mrca_age", "mean_terminal_length", "mpd_z_score", "diversification_metric")
  aggregated_numeric <- aggregate(
    combined_results[, numeric_cols], 
    by = list(AoE = combined_results$AoE), 
    FUN = function(x) {
      c(mean = mean(x, na.rm = TRUE), 
        sd = sd(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE))
    }
  )
  
  # Reshape the aggregated data for easier use
  reshaped_numeric <- data.frame(AoE = aggregated_numeric$AoE)
  for(col in numeric_cols) {
    reshaped_numeric[[paste0(col, "_mean")]] <- aggregated_numeric[[col]][,"mean"]
    reshaped_numeric[[paste0(col, "_sd")]] <- aggregated_numeric[[col]][,"sd"]
    reshaped_numeric[[paste0(col, "_min")]] <- aggregated_numeric[[col]][,"min"]
    reshaped_numeric[[paste0(col, "_max")]] <- aggregated_numeric[[col]][,"max"]
  }
  
  # Count frequency of each phylogenetic structure type
  structure_counts <- table(combined_results$AoE, combined_results$phylo_structure)
  structure_proportions <- prop.table(structure_counts, margin = 1)
  
  # Determine the dominant structure type for each AoE
  dominant_structure <- apply(structure_counts, 1, function(x) {
    names(which.max(x))
  })
  
  # Calculate dominant structure proportion
  dominant_proportion <- numeric(length(dominant_structure))
  for(i in 1:length(dominant_structure)) {
    aoe_name <- names(dominant_structure)[i]
    struct_type <- dominant_structure[i]
    if(!is.na(struct_type)) {
      dominant_proportion[i] <- structure_proportions[aoe_name, struct_type]
    } else {
      dominant_proportion[i] <- NA
    }
  }
  
  # Create structure summary
  structure_summary <- data.frame(
    AoE = names(dominant_structure),
    dominant_structure = unname(dominant_structure),
    structure_proportion = dominant_proportion,
    n_trees = apply(structure_counts, 1, sum)
  )
  
  # Get species counts (should be constant across trees)
  species_counts <- aggregate(n_species ~ AoE, data = combined_results, FUN = function(x) x[1])
  
  # Get monophyly frequencies
  monophyly_counts <- aggregate(
    is_monophyletic ~ AoE, 
    data = subset(combined_results, !is.na(is_monophyletic)), 
    FUN = function(x) mean(as.numeric(x))
  )
  names(monophyly_counts)[2] <- "monophyly_freq"
  
  # Combine all summaries
  final_summary <- merge(species_counts, structure_summary, by = "AoE")
  if(nrow(monophyly_counts) > 0) {
    final_summary <- merge(final_summary, monophyly_counts, by = "AoE", all.x = TRUE)
  } else {
    final_summary$monophyly_freq <- NA
  }
  final_summary <- merge(final_summary, reshaped_numeric, by = "AoE")
  
  # Order by number of species
  final_summary <- final_summary[order(-final_summary$n_species), ]
  
  return(list(
    all_results = combined_results,
    summary = final_summary
  ))
}

# Run the analysis on 10 trees
set.seed(123) # For reproducibility
multi_tree_results <- analyze_multiple_trees(fully_sampled_rep, aoe_species_list, n_trees = 100)

# Display summary
cat("\nSummary across 10 trees:\n")
summary_table <- multi_tree_results$summary

# Format the summary table for manuscript
manuscript_table <- data.frame(
  AoE = summary_table$AoE,
  n_species = summary_table$n_species,
  monophyly_freq = round(summary_table$monophyly_freq, 2),
  mrca_age = round(summary_table$mrca_age_mean, 2),
  mrca_age_sd = round(summary_table$mrca_age_sd, 2),
  terminal_branch = round(summary_table$mean_terminal_length_mean, 2),
  terminal_branch_sd = round(summary_table$mean_terminal_length_sd, 2),
  mpd_z = round(summary_table$mpd_z_score_mean, 2),
  mpd_z_sd = round(summary_table$mpd_z_score_sd, 2),
  diversification = round(summary_table$diversification_metric_mean, 3),
  diversification_sd = round(summary_table$diversification_metric_sd, 3),
  dominant_structure = summary_table$dominant_structure,
  structure_support = paste0(round(summary_table$structure_proportion * 100), "%")
)

# Display the formatted table
print(manuscript_table)

# Create visualizations 
# 1. Mean MRCA age with error bars
p1 <- ggplot(summary_table, aes(x=reorder(AoE, mrca_age_mean), y=mrca_age_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mrca_age_mean-mrca_age_sd, ymax=mrca_age_mean+mrca_age_sd), width=0.2) +
  coord_flip() +
  labs(title="Mean MRCA Age by Area of Endemism (across 10 trees)", 
       x="Area of Endemism", 
       y="MRCA Age (Million Years)")

# 2. Mean MPD Z-score with error bars
p2 <- ggplot(summary_table, aes(x=reorder(AoE, mpd_z_score_mean), y=mpd_z_score_mean)) +
  geom_bar(stat="identity", aes(fill=mpd_z_score_mean < -1.96 | mpd_z_score_mean > 1.96)) +
  geom_errorbar(aes(ymin=mpd_z_score_mean-mpd_z_score_sd, ymax=mpd_z_score_mean+mpd_z_score_sd), width=0.2) +
  geom_hline(yintercept=c(-1.96, 1.96), linetype="dashed", color="red") +
  scale_fill_manual(values=c("grey", "darkred"), name="Significant") +
  coord_flip() +
  labs(title="Phylogenetic Structure by Area of Endemism (across 10 trees)", 
       x="Area of Endemism", 
       y="MPD Z-Score (Negative = Clustered, Positive = Overdispersed)")

# 3. Diversification metric with error bars
p3 <- ggplot(subset(summary_table, !is.na(diversification_metric_mean)), 
             aes(x=reorder(AoE, diversification_metric_mean), y=diversification_metric_mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=diversification_metric_mean-diversification_metric_sd, 
                    ymax=diversification_metric_mean+diversification_metric_sd), width=0.2) +
  coord_flip() +
  labs(title="Species Accumulation Rate by Area of Endemism (across 100 trees)", 
       x="Area of Endemism", 
       y="ln(species)/MRCA Age")

# Display plots
print(p1)
print(p2)
print(p3)

library(dplyr)
# Calculate proportion of trees showing each structure type
structure_summary <- multi_tree_results$all_results %>% 
  group_by(AoE, phylo_structure) %>%
  summarize(count = n()) %>%
  group_by(AoE) %>%
  mutate(proportion = count / sum(count))

# Create stacked bar chart showing structure type proportions
p4 <- ggplot(structure_summary, aes(x=AoE, y=proportion, fill=phylo_structure)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette="Set2") +
  coord_flip() +
  labs(title="Proportion of Trees Supporting Each Phylogenetic Structure Type", 
       x="Area of Endemism", 
       y="Proportion of Trees",
       fill="Phylogenetic Structure")

# Display the chart
print(p4)

# Create a map of diversification processes for manuscript
# This is a suggested visualization structure that would need to be implemented
# with your actual spatial data
cat("\nSuggested visualization for manuscript:\n")
cat("Create a map of Amazonia showing each AoE colored by:\n")
cat("1. Dominant phylogenetic structure (clustered vs. random vs. overdispersed)\n")
cat("2. Mean age of endemic lineages\n")
cat("3. Diversification rate\n")

# Summary of findings to address reviewer comments
cat("\nSummary for addressing reviewer comments:\n")
cat("This analysis reveals that:\n")

clustered_aoes <- manuscript_table$AoE[manuscript_table$dominant_structure == "clustered_deep"]
if(length(clustered_aoes) > 0) {
  cat("- The following AoEs show evidence of in-situ diversification (phylogenetic clustering):\n  ", 
      paste(clustered_aoes, collapse=", "), "\n")
}

overdispersed_aoes <- manuscript_table$AoE[manuscript_table$dominant_structure == "overdispersed_deep"]
if(length(overdispersed_aoes) > 0) {
  cat("- The following AoEs show evidence of multiple colonization events (phylogenetic overdispersion):\n  ", 
      paste(overdispersed_aoes, collapse=", "), "\n")
}

# Oldest AoEs
oldest_aoes <- manuscript_table$AoE[order(-manuscript_table$mrca_age)][1:3]
cat("- The oldest endemic assemblages are found in:\n  ", paste(oldest_aoes, collapse=", "), "\n")

# Highest diversification rates
highest_div_aoes <- manuscript_table$AoE[order(-manuscript_table$diversification)][1:3]
cat("- The highest diversification rates are found in:\n  ", paste(highest_div_aoes, collapse=", "), "\n")

