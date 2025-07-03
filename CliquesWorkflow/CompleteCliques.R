
# --- Identifying COMPLETE cliques: CONSERVED ACROSS ALL SPECIES (complete and partial w/ missing) and LINEAGE-SPECIFIC ---

library(tidyverse)
library(igraph)

setwd("Cliques")

# Load necessary files
load("DATA/all_expressed_genes.RData")
load("DATA/ones_and_zeros_all.RData")  

# Table for identifying which orthogroups have expressed genes for the various species.
df_with_sums <- ones_and_zeros_all %>%
  mutate(Angiosperms = rowSums(. [1:3])) %>%
  mutate(Cross = rowSums(. [4:12]))%>%
  mutate(Gymnosperms = rowSums(. [13:15])) %>%
  mutate(Conserved = rowSums(. [1:15]))%>%
  filter(Conserved >= 3) 

# Filtering expressed genes for co-expressologs
pval <- 0.1
expressologs <- expr_genes %>%
  filter(Max.p.Val < pval) %>%
  group_by(OrthoGroup) %>% 
  mutate(OG_size = n()) %>% 
  ungroup() %>% 
  arrange(desc(OG_size)) %>% 
  filter(OG_size >= 3)

# cat("\n","Largest orthogroup: ", max(expressologs$OG_size), " gene pairs.", "\n", "Average orthogroup: ", 
#     round(mean(expressologs$OG_size), 0), " gene pairs.", "\n", "Smallest orthogroup: ", min(expressologs$OG_size), " gene pairs.")  

OGs_to_use <- unique(expressologs$OrthoGroup)

# Clique algorithm 
expressologs_from_max_cliques <- list() 
for (i in 1:length(OGs_to_use)) {

  print(paste0(i, "/", length(OGs_to_use)))
  
  # For each OG a network is created using genes as nodes and 
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
  
  if (length(largest_clique) > 0){
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

complete_cliques_filterable <- clique_genes_unlisted %>%
  mutate(species_pair = paste0(Species1, Species2)) %>% 
  group_by(cliqueID)%>% 
  mutate(OriginalCliqueSize = n())%>%  
  ungroup()

saveRDS(complete_cliques_filterable, file = "DATA/cliqueOutput/clique_genes_filterable_COMPLETE.RDS")



# --------------- CONSERVED ACROSS ALL SPECIES -----------------------------
# Identify 6M cliques
caas <-readRDS("DATA/cliqueOutput/clique_genes_filterable_COMPLETE.RDS") %>% 
  filter(OriginalCliqueSize == 15) %>% 
  select(Species1, GeneSpecies1, Species2, GeneSpecies2, OrthoGroup, cliqueID, Max.p.Val)

spc1 <- caas %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies1 = paste0(caas$Species1, "-", caas$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- caas %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies2 = paste0(caas$Species2, "-", caas$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)

caas_wider <- rbind(spc1,spc2)  

rm(spc1)
rm(spc2)


caas_wider <- caas_wider %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Asp, Birch, Cher, Nor, Scots, Lodge) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup() %>% 
  mutate(conservedType = "complete")

length(unique(caas_wider$OrthoGroup))
length(unique(caas_wider$cliqueID))

saveRDS(caas_wider, file = "DATA/cliqueOutput/conserved_genes_COMPLETE_wider.RDS")

conserved_complete_OG <- caas_wider %>% distinct(OrthoGroup) %>% pull(OrthoGroup)
saveRDS(conserved_complete_OG, file = "DATA/conserved_complete_OG.RDS") # Save OGs in order to remove conserved genes from other sets.


# ------------------DICOT-SPECIC GENES----------------------

angio_specific_OG <- df_with_sums %>%
  filter(Angiosperms == 3 & Gymnosperms == 0 & Cross == 0) %>% 
  rownames_to_column(var = "OG") %>% 
  pull(OG)

angio <-readRDS("DATA/cliqueOutput/clique_genes_filterable_COMPLETE.RDS") %>% 
  filter(OrthoGroup %in% angio_specific_OG) 

spc1 <- angio %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies1 = paste0(angio$Species1, "-", angio$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- angio %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies2 = paste0(angio$Species2, "-", angio$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)

df <- rbind(spc1,spc2)  

rm(spc1)
rm(spc2)

angio_wider <- df %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Asp, Birch, Cher) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup()

rm(df)
rm(angio)

length(unique(angio_wider$OrthoGroup))
length(unique(angio_wider$cliqueID))

saveRDS(angio_wider, file = "DATA/cliqueOutput/dicto_specific_genes_wider.RDS")

# ------------------GYMNO-SPECIC GENES----------------------
gymno_specific_OG <- df_with_sums %>% 
  filter(Angiosperms == 0 & Gymnosperms == 3 & Cross == 0)%>% 
  rownames_to_column(var = "OG") %>% 
  pull(OG)

gymno <-readRDS("DATA/cliqueOutput/clique_genes_filterable_COMPLETE.RDS") %>% 
  filter(OrthoGroup %in% gymno_specific_OG) 

spc1 <- gymno %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies1 = paste0(gymno$Species1, "-", gymno$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- gymno %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies2 = paste0(gymno$Species2, "-", gymno$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)

df <- rbind(spc1,spc2)  

rm(spc1)
rm(spc2)

gymno_wider <- df %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Nor, Scots, Lodge) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup()

rm(df)
rm(gymno)

length(unique(gymno_wider$OrthoGroup))
length(unique(gymno_wider$cliqueID))

saveRDS(gymno_wider, file = "DATA/cliqueOutput/conifer_specific_genes_wider.RDS")




