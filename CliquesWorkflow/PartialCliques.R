
# --- Identifying PARTIAL cliques: CONSERVED ACROSS ALL SPECIES (partial w/ non-significant edges)



# Partial_A and partial_B cliques will be combined to one file.

library(tidyverse)
library(igraph)

# Load necessary files
load("DATA/all_expressed_genes.RData")
load("DATA/ones_and_zeros_all.RData") 


# Table for identifying which orthogroups have expressed genes for the various species.
df_with_sums <- ones_and_zeros_all %>%
  mutate(Angiosperms = rowSums(. [1:3])) %>%
  mutate(Cross = rowSums(. [4:12]))%>%
  mutate(Gymnosperms = rowSums(. [13:15])) %>%
  mutate(Conserved = rowSums(. [1:15]))%>%
  filter(Conserved >= 3 ) 

conserved_complete_OG <- readRDS("DATA/cliqueOutput/conserved_complete_OG.RDS")

pval <- 0.9
expressologs <- expr_genes %>%
  filter(Max.p.Val < pval) %>%
  # filter(!(OrthoGroup %in% discard_list)) %>% 
  filter(!(OrthoGroup %in% conserved_complete_OG)) %>% 
  group_by(OrthoGroup) %>% 
  mutate(OG_size = n()) %>% 
  ungroup() %>% 
  arrange(desc(OG_size)) %>% 
  filter(OG_size >= 3)

OGs_to_use <- unique(expressologs$OrthoGroup)


# Clique algorithm 
expressologs_from_max_cliques <- list() 
for (i in 1:length(OGs_to_use)) {
  # i <- 874
  print(paste0(i, "/", length(OGs_to_use)))
  
  g <- OGs_to_use[i]
  
  expressologs_g <- expressologs %>% 
    filter(OrthoGroup == g) 
  
  expressologs_g <- expressologs_g %>% 
    mutate(UGeneSpecies1 = paste0(expressologs_g$Species1, "-", expressologs_g$GeneSpecies1),
           UGeneSpecies2 = paste0(expressologs_g$Species2, "-", expressologs_g$GeneSpecies2))
  
  nodes <- data.frame(name = unique(c(expressologs_g$UGeneSpecies1, expressologs_g$UGeneSpecies2)))
  edges <- data.frame(from = expressologs_g$UGeneSpecies1, 
                      to = expressologs_g$UGeneSpecies2)
  
  net <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  largest_clique <- max_cliques(net, min = 3) # list of all cliques from size 3 to 6
  
  if (length(largest_clique) >0){
    largest_clique_ok <- largest_clique
  }
  
  for(c in 1:length(largest_clique_ok)){
    
    if (c %% 1000 == 0) {
      cat(c, "\n")
    }
    
    clique <- largest_clique_ok[[c]]
    clique_name <- attr(clique, "names")
    
    clique_species <- str_split_fixed(clique_name, "-", n = 2)[,1]
    clique_genes <-  str_split_fixed(clique_name, "-", n = 2)[,2]
    
    expressologs_in_clique <- expressologs_g %>%
      filter(UGeneSpecies1 %in% clique_name & UGeneSpecies2 %in% clique_name) %>%
      select( OrthoGroup, Species1, GeneSpecies1, Species2, GeneSpecies2, Max.p.Val) %>%
      mutate(cliqueID = paste0(g,"_" ,c))
    
    
    expressologs_from_max_cliques <- append(expressologs_from_max_cliques, list(expressologs_in_clique))
    
  }
  
  
}

clique_genes_unlisted <- plyr::ldply(expressologs_from_max_cliques)

partial_cliques_filterable <- clique_genes_unlisted %>%
  mutate(species_pair = paste0(Species1, Species2)) %>% 
  group_by(cliqueID)%>% 
  mutate(OriginalCliqueSize = n())%>%  
  ungroup() 

saveRDS(partial_cliques_filterable, file = "DATA/cliqueOutput/clique_genes_filterable_PARTIAL_0.9.RDS")


# ------ PARTIAL A --------

partial_A <- partial_cliques_filterable %>% 
  filter(Max.p.Val < 0.1) %>% 
  group_by(cliqueID) %>% 
  mutate(FiltCliqueSize = n()) %>% 
  ungroup() %>% 
  filter(FiltCliqueSize >= 11) # Require at least 11 significant edges (OriginalCLiqueSize will be 15)

# Collect cliqueIDs from partial A
partial_A_cliqueID <- partial_A$cliqueID 

partial <- partial_cliques_filterable %>% 
  filter(cliqueID %in% partial_A_cliqueID)
  
spc1 <- partial %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val, OriginalCliqueSize) %>% 
  mutate(UGeneSpecies1 = paste0(partial$Species1, "-", partial$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- partial %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val, OriginalCliqueSize) %>% 
  mutate(UGeneSpecies2 = paste0(partial$Species2, "-", partial$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)

df <- rbind(spc1,spc2)  

rm(spc1)
rm(spc2)


partial_A_wider <- df %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Asp, Birch, Cher, Nor, Scots, Lodge, OriginalCliqueSize) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup() %>% 
  mutate(conservedType = "partiallySign")

length(unique(partial_A_wider$OrthoGroup))
length(unique(partial_A_wider$cliqueID))
nrow(partial_A_wider)

# ----- PARTIAL B ------

# Read in file with the complete cliques
all_complete_cliques <- readRDS("DATA/cliqueOutput/clique_genes_filterable_COMPLETE.RDS")

partial_B <- all_complete_cliques %>% 
  filter(!(OrthoGroup %in% conserved_complete_OG)) %>% # remove OGs with conserved complete and partial cliques
  filter(!(OrthoGroup%in%partial_A_wider$OrthoGroup)) %>%
  filter(OriginalCliqueSize == 10) # find 5M cliques

spc1 <- partial_B %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val, OriginalCliqueSize) %>% 
  mutate(UGeneSpecies1 = paste0(partial_B$Species1, "-", partial_B$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- partial_B %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val, OriginalCliqueSize) %>% 
  mutate(UGeneSpecies2 = paste0(partial_B$Species2, "-", partial_B$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)

df <- rbind(spc1,spc2)  

rm(spc1)
rm(spc2)

partial_B_wider <- df %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Asp, Birch, Cher, Nor, Scots, Lodge, OriginalCliqueSize) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup() %>% 
  mutate(conservedType = "partiallyPres" ) 

length(unique(partial_B_wider$OrthoGroup))
length(unique(partial_B_wider$cliqueID))
nrow(partial_B_wider)

# ---- COMBINING PARTIAL CONSERVED GENES ACROSS ALL SPECIES -----

p_caas_wider <- rbind(partial_A_wider, partial_B_wider)
nrow(p_caas_wider)

length(unique(p_caas_wider$OrthoGroup))
length(unique(p_caas_wider$cliqueID))

saveRDS(p_caas_wider, file = "DATA/cliqueOutput/conserved_genes_PARTIAL_wider.RDS")





