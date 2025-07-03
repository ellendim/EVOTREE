
# SCRIPT FOR IDENTIFYING CO-EXPRESSOLOGS USING COMPLEX

library(dplyr)
library(tidyr)
library(gdata)


run_version <- "5"
ortholog_group_file <- "DATA/OrthoGroups/Orthogroups_20240823_clean_N1.tsv.gz"
arabidopsis_annot <- "DATA/Annotations/gene_aliases_20140331.txt"

species_list <- c("Lodge", "Asp", "Nor","Scots","Birch","Cher")
pairCombinations <- data.frame(t(combn(species_list, 2)))


# Parameters
cor_method <- "pearson" # pearson spearman
cor_sign <- "" # abs
norm_method <- "MR" # CLR MR
density_thr <- 0.03
randomize <- "no" # yes no

# For test-runs
test_run <- "no" # "yes" "no"  
numb_of_rows <- 1000

# -------------Run ComPlEx --------------
for(row in 1:nrow(pairCombinations)){
  
  species1_keyword <- pairCombinations[row,1]
  species2_keyword <- pairCombinations[row,2] 
  
  outfile <- paste0("DATA/comparisonFiles/orthologs/orthologs-",species1_keyword,"-",species2_keyword,"-table-V",run_version,".RData")
  
  if (!file.exists(outfile)){ 
    
    ortho_original <- read.delim2(ortholog_group_file, header = TRUE, sep = "\t")
    
    ortho_general_filtering <- ortho_original %>%
      mutate(Pinus_sylvestris_cp = Pinus_sylvestris) %>% 
      rename(
        Asp = Populus_tremula,
        Birch = Betula_pendula,
        Nor = Picea_abies,
        Scots = Pinus_sylvestris,
        Cher = Prunus_avium,
        Lodge = Pinus_sylvestris_cp)%>%
      select(OrthoGroup, Asp, Nor, Scots, Birch, Cher, Lodge)
    
    ortho <- ortho_general_filtering %>%
      rename(Species1 = species1_keyword, Species2 = species2_keyword) %>%
      select(OrthoGroup, Species1, Species2) %>%
      filter(Species1 != "", Species2 != "") %>%
      separate_rows(Species1, sep = ", ", convert = FALSE) %>%
      separate_rows(Species2, sep = ", ", convert = FALSE)
    
    # Additional species-specific alterations 
    if(species1_keyword == "Asp"){
      
      ortho <- ortho %>% 
        mutate(Species1 = gsub("\\.\\d\\.p\\d$", "", Species1))
    }
    
    if(species2_keyword == "Asp"){
      ortho <- ortho %>% 
        mutate(Species2 = gsub("\\.\\d\\.p\\d$", "", Species2))
      
    }
    
    if(species1_keyword == "Nor"){
      
      ortho <- ortho %>% 
        mutate(Species1 = gsub("\\.p\\d$", "", Species1))
    }
    
    if(species2_keyword == "Nor"){
      ortho <- ortho %>% 
        mutate(Species2 = gsub("\\.p\\d$", "", Species2))
      
    }
    
    if(species1_keyword == "Cher"){
      
      ortho <- ortho %>% 
        mutate(Species1 = gsub("\\.p\\d$", "", Species1))
    }
    
    if(species2_keyword == "Cher"){
      ortho <- ortho %>% 
        mutate(Species2 = gsub("\\.p\\d$", "", Species2))
      
    }
    ortho <- ortho %>% 
      group_by(Species1, Species2, OrthoGroup) %>% 
      slice(1)
    
    # Add annotations from Arabidopsis
    symbols <- read.delim(arabidopsis_annot, sep = "\t") %>%
      rename(Arab = locus_name,
             Symbol = symbol,
             Name = full_name)
    
    annot <-
      read.delim2(ortholog_group_file, header = TRUE, sep = "\t") %>%
      rename(Arab = Arabidopsis_thaliana) %>%
      select(Arab, OrthoGroup) %>%
      filter(Arab != "") %>%
      separate_rows(Arab, sep = ", ", convert = FALSE) %>%
      mutate(Arab = gsub("\\.\\d.p\\d$", "", Arab)) %>%
      mutate(Arab = gsub("Aratha_", "", Arab)) %>%
      left_join(symbols, by = "Arab", relationship = "many-to-many") %>%
      group_by(OrthoGroup) %>%
      summarise(
        Arab = paste0(unique(Arab), collapse = "; "),
        Symbol = paste0(unique(Symbol), collapse = "; "),
        Name = paste0(unique(Name), collapse = "; ")
      ) %>%
      mutate(Symbol = gsub("NA", "", Symbol),
             Name = gsub("NA", "", Name))
    
    save(ortho, annot, file = outfile)
    
    # Load transcriptomic files for both species
    species1_transcription_txt <- paste0(species1_keyword,"Wood/", species1_keyword,"Wood_transcriptomics.txt.gz")
    species2_transcription_txt <- paste0(species2_keyword,"Wood/", species2_keyword,"Wood_transcriptomics.txt.gz")
    
    if(test_run == "yes"){
      # Only a subset of genes and samples is run
      species1_expr <- read.delim(species1_transcription_txt, sep = "\t", header = TRUE)[1:numb_of_rows, ]
      species2_expr <- read.delim(species2_transcription_txt, sep = "\t", header = TRUE)[1:numb_of_rows, ]
      
    } else{
      # Full expression set
      species1_expr <- read.delim(species1_transcription_txt, sep = "\t", header = TRUE)
      species2_expr <- read.delim(species2_transcription_txt, sep = "\t", header = TRUE)
    }
    
    ortho <- ortho %>%
      filter(Species1 %in% species1_expr$Genes & Species2 %in% species2_expr$Genes)
    
    # Setting up for species comparison
    comparison_RData <- paste0("DATA/comparisonFiles/comparison-", species1_keyword, "_", species2_keyword, "-", 
                               cor_sign, cor_method, norm_method, density_thr, randomize, "-table-V",run_version,".RData")
    
    if (!file.exists(comparison_RData)) {
      if (randomize == "yes") {
        species1_expr$Genes <-
          sample(species1_expr$Genes, nrow(species1_expr), FALSE)
        species2_expr$Genes <-
          sample(species2_expr$Genes, nrow(species2_expr), FALSE)
      }
      
      species1_net <- cor(t(species1_expr[, -1]), method = cor_method) 
      dimnames(species1_net) <- list(species1_expr$Genes, species1_expr$Genes)               
      
      species2_net <- cor(t(species2_expr[, -1]), method = cor_method)  
      dimnames(species2_net) <- list(species2_expr$Genes, species2_expr$Genes)                 
      
      if (cor_sign == "abs") {
        species1_net <- abs(species1_net)
        species2_net <- abs(species2_net)
      }
      
      # Apply normalisation method
      if (norm_method == "CLR") {
        z <- scale(species1_net)
        z[z < 0] <- 0
        species1_net <- sqrt(t(z) ** 2 + z ** 2)
        
        z <- scale(species2_net)
        z[z < 0] <- 0
        species2_net <- sqrt(t(z) ** 2 + z ** 2)
        
      } else if (norm_method == "MR") {
        R <- t(apply(species1_net, 1, rank)) 
        species1_net <- sqrt(R * t(R)) 
        
        R <- t(apply(species2_net, 1, rank))
        species2_net <- sqrt(R * t(R))
      }
      
      diag(species1_net) <- 0
      diag(species2_net) <- 0
      
      R <- sort(species1_net[upper.tri(species1_net, diag = FALSE)], decreasing = TRUE) 
      species1_thr <- R[round(density_thr * length(R))] 
      
      R <-sort(species2_net[upper.tri(species2_net, diag = FALSE)], decreasing = TRUE)
      species2_thr <- R[round(density_thr * length(R))]
      
      comparison <- ortho
      
      comparison$Species1.neigh <- c(NA)
      comparison$Species1.ortho.neigh <- c(NA)
      comparison$Species1.neigh.overlap <- c(NA)
      comparison$Species1.p.val <- c(NA)
      
      comparison$Species2.neigh <- c(NA)
      comparison$Species2.ortho.neigh <- c(NA)
      comparison$Species2.neigh.overlap <- c(NA)
      comparison$Species2.p.val <- c(NA)
      
      # Perform two-way co-expression network comparison
      for (i in 1:nrow(ortho)) {
        
        if (i %% 100 == 0) {
          cat(i, "\n")
        }
        
        # Species 1 -> Species 2
        neigh <- species1_net[ortho$Species1[i],]  
        neigh <- names(neigh[neigh >= species1_thr]) 
        
        ortho_neigh <- species2_net[ortho$Species2[i],] 
        ortho_neigh <- names(ortho_neigh[ortho_neigh >= species2_thr])
        ortho_neigh <- ortho$Species1[ortho$Species2 %in% ortho_neigh]
        
        # Hypergeometric test 1
        N <- nrow(species1_expr)
        m <- length(neigh)
        n <- N-m
        k <- length(unique(ortho_neigh))
        x <- length(unique(intersect(neigh, ortho_neigh)))
        p_val <- 1
        if (x > 1) {
          p_val <- phyper(x-1, m, n, k, lower.tail = FALSE)
        }
        
        comparison$Species1.neigh[i] <- m
        comparison$Species1.ortho.neigh[i] <- k
        comparison$Species1.neigh.overlap[i] <- x
        comparison$Species1.p.val[i] <- p_val
        
        # Species 2 -> Species 1
        neigh <- species2_net[ortho$Species2[i],]
        neigh <- names(neigh[neigh >= species2_thr])
        
        ortho_neigh <- species1_net[ortho$Species1[i],]
        ortho_neigh <- names(ortho_neigh[ortho_neigh >= species1_thr])
        ortho_neigh <- ortho$Species2[ortho$Species1 %in% ortho_neigh]
        
        # Hypergeometric test 2
        N <- nrow(species2_expr)
        m <- length(neigh)
        n <- N-m
        k <- length(unique(ortho_neigh))
        x <- length(unique(intersect(neigh, ortho_neigh)))
        p_val <- 1
        
        if (x > 1) {
          p_val <- phyper(x-1, m, n, k, lower.tail = FALSE)
        }
        
        comparison$Species2.neigh[i] <- m
        comparison$Species2.ortho.neigh[i] <- k
        comparison$Species2.neigh.overlap[i] <- x
        comparison$Species2.p.val[i] <- p_val
      }
      
      
      # Remove orthologs not present in either species network
      comparison <- comparison %>%
        filter(Species1.neigh.overlap > 0 & Species2.neigh.overlap > 0)
      
      # Correct for multiple testing (FDR) and compile comparison table
      comparison$Species1.p.val <- p.adjust(comparison$Species1.p.val, method = "fdr")
      comparison$Species2.p.val <- p.adjust(comparison$Species2.p.val, method = "fdr")
      
      
      comparison_table <- comparison %>%
        rowwise() %>%
        mutate(Max.p.val = max(Species1.p.val, Species2.p.val)) %>%
        left_join(annot, by = "OrthoGroup") %>%
        select(-c("Species1.neigh", "Species1.ortho.neigh", "Species2.neigh", "Species2.ortho.neigh")) %>%
        arrange(Max.p.val)
      
      comparison_table$Species1.p.val <- format(comparison_table$Species1.p.val, digits = 3, scientific = TRUE)
      comparison_table$Species2.p.val <- format(comparison_table$Species2.p.val, digits = 3, scientific = TRUE)
      comparison_table$Max.p.val <- format(comparison_table$Max.p.val, digits = 3, scientific = TRUE)
      
      save(comparison_table, file = comparison_RData)
    }
    
  } else{
    print("File already exists")
    
  }
  
  
}















