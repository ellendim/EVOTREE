
# IDENTIFYING THE BEST UNIQUE CLIQUES PER GENE SET
# This script identifies all unique cliques (i.e. no two cliques overlap in genes) with the lowest FDR sum.

library(tidyverse)


# Choose clique type
clique_types <- c("conserved_genes_COMPLETE") # Change this depending on which gene set!

data <- readRDS(paste0("DATA/cliqueOutput/", clique_types, "_wider.RDS"))

if(clique_types %in% c("conserved_genes_COMPLETE","conserved_genes_PARTIAL" ,"conifer_specific_genes", "differentiated_genes" )){
  data <- data %>%  mutate(Scots = paste0(Scots, "_S"))
}

if(clique_types == "dicot_specific_genes"){
  species <- c("Asp", "Birch","Cher")
}
if(clique_types == "conifer_specific_genes"){
  species <- c("Nor", "Scots","Lodge")
}

if(clique_types %in% c("conserved_genes_COMPLETE", "conserved_genes_PARTIAL", "differentiated_genes" )){
  species <- c("Asp", "Birch","Cher","Nor", "Scots","Lodge")
}


best_unique_cliques <- c()
for(OG in unique(data$OrthoGroup)){
  
  print(OG)
  
  df <- data %>% 
    filter(OrthoGroup == OG) %>% 
    arrange(MeanCliqueSum)
  
  nrows <- nrow(df)
  
  # Best scored clique
  clique <- df[1,]  
  ID <- unname(unlist(clique[, 2]))
  score <- unname(unlist(clique[, 3]))
  vector <-unname(unlist(clique[1, species]))
  
  best_unique_cliques <- append(best_unique_cliques,ID) # add cliqueID to list
  
  #idx <- idx_start + 1
  
  genes <- vector
  
  # Always take the best, then go through the rest
  for(i in 2:nrows){
    
    c_clique <- df[i,]
    c_ID <- unname(unlist(c_clique[, 2]))
    c_score <- unname(unlist(c_clique[, 3]))
    c_vector <-unname(unlist(c_clique[1, species]))
    
    # Overlap between genes in this clique and genes in all the other cliques we have selected so far
    Overlap <- sum(c_vector %in% genes)
    
    if(Overlap == 0){
      print(paste0(" ", c_ID))
      
      best_unique_cliques <- append(best_unique_cliques,c_ID)
      
      genes <- append(genes,c_vector)
    }
  }
}

if(clique_types %in% c("conserved_genes_COMPLETE" ,"conserved_genes_PARTIAL" , "conifer_specific_genes", "differentiated_genes" )){
  data <- data %>% mutate(Scots = gsub("_S", " ", Scots))
}

# Save output
