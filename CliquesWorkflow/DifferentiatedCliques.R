library(tidyverse)
library(igraph)



caas <- readRDS("DATA/cliqueOutput/conserved_genes_COMPLETE_wider.RDS")
p_caas <- readRDS("DATA/cliqueOutput/conserved_genes_PARTIAL_wider.RDS")
dicot <- readRDS("DATA/cliqueOutput/dicot_specific_genes_wider.RDS")
conifer <- readRDS("DATA/cliqueOutput/conifer_specific_genes_wider.RDS")

# --- Identifying DIFFERENTIATED cliques ---

load("DATA/ones_and_zeros_all.RData")
load("DATA/all_expressed_genes.RData")

# Remove OGs containing conserved or lineage-specific genes
df_with_sums <- ones_and_zeros_all %>%
  mutate(Angiosperms = rowSums(. [1:3])) %>%
  mutate(Cross = rowSums(. [4:12]))%>%
  mutate(Gymnosperms = rowSums(. [13:15])) %>%
  mutate(Conserved = rowSums(. [1:15]))%>%
  filter(Conserved == 15 ) 

none_of_these <- c(unique(caas$OrthoGroup), unique(p_caas$OrthoGroup))

leftoverOG <- df_with_sums %>% 
  filter(!(rownames(df_with_sums)%in%none_of_these)) %>% 
  filter(Angiosperms == 3 & Gymnosperms == 3)

expressologs <- expr_genes %>%
  filter(OrthoGroup %in% rownames(leftoverOG)) %>% 
  group_by(OrthoGroup) %>% 
  mutate(OG_size = n()) %>% 
  ungroup() %>% 
  arrange(desc(OG_size)) %>% 
  filter(OG_size >=15)

OGs_to_use <- unique(expressologs$OrthoGroup)

expressologs_from_max_cliques <- list() 
for (i in 1:length(OGs_to_use)) {
  # i <- 40
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
  
  
  for(c in 1:length(largest_clique)){
    
    if (c %% 1000 == 0) {
      cat(c, "\n")
    }
    
    clique <- largest_clique[[c]]
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

length(unique(complete_cliques_filterable$OrthoGroup))
length(unique(complete_cliques_filterable$cliqueID))
unique(complete_cliques_filterable$OriginalCliqueSize)

# Setting requirements for cliques
speciesPairClade <- tibble(pair = unique(complete_cliques_filterable$species_pair), 
                           pairClade = c("Cross", "Gymno", "Gymno", "Cross", "Cross", "Cross", "Cross", "Angio", "Angio","Gymno", "Cross", "Cross", 
                                         "Cross", "Cross", "Angio"))

dbl_with_clades <- left_join(complete_cliques_filterable, speciesPairClade, by = join_by(species_pair == pair))
dbl_with_clades <- dbl_with_clades %>% pivot_wider(values_from = species_pair, names_from = pairClade)

counting <- ifelse(is.na(dbl_with_clades[9:11]), "0", "1")
counting <- counting %>% as.data.frame() %>% mutate_if(is.character, as.numeric)

dbl_filterable <- cbind(dbl_with_clades[1:8], counting)
dbl_filterable <- dbl_filterable %>% 
  group_by(cliqueID) %>% 
  mutate(AngioSum = sum(Angio))%>% 
  mutate(GymnoSum = sum(Gymno))%>% 
  mutate(CrossSum = sum(Cross)) %>% 
  ungroup() %>% 
  select(-c(9:11))%>% 
  filter(Max.p.Val < 0.1) %>% 
  group_by(cliqueID) %>% 
  filter(OriginalCliqueSize == 15) %>% 
  filter(AngioSum == 3 & GymnoSum == 3)

length(unique(dbl_filterable$OrthoGroup))
length(unique(dbl_filterable$cliqueID))
unique(dbl_filterable$OriginalCliqueSize)

dbl <- complete_cliques_filterable %>% filter(cliqueID %in% dbl_filterable$cliqueID)
  
# --- Formatting ---
spc1 <- dbl %>% 
  select( Species1,GeneSpecies1,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies1 = paste0(dbl$Species1, "-", dbl$GeneSpecies1)) %>% 
  rename("Genes" = GeneSpecies1) %>% 
  rename("UGenes" = UGeneSpecies1)%>% 
  rename("Species" = Species1)

spc2 <- dbl %>% 
  select(Species2,GeneSpecies2,OrthoGroup, cliqueID, Max.p.Val) %>% 
  mutate(UGeneSpecies2 = paste0(dbl$Species2, "-", dbl$GeneSpecies2)) %>% 
  rename("Genes" = GeneSpecies2) %>% 
  rename("UGenes" = UGeneSpecies2) %>% 
  rename("Species" = Species2)


dbl_wider <- rbind(spc1, spc2) %>%
  group_by(cliqueID) %>%
  distinct(UGenes, .keep_all = T) %>%
  mutate(MeanCliqueSum = mean(Max.p.Val)) %>%
  ungroup() %>%
  select(-c(UGenes, Max.p.Val)) %>% 
  pivot_wider(names_from = Species, values_from = Genes) %>% 
  select(OrthoGroup, cliqueID, MeanCliqueSum, Asp, Birch, Cher, Nor, Scots, Lodge) %>% 
  group_by(OrthoGroup) %>% 
  arrange(OrthoGroup,MeanCliqueSum) %>% 
  ungroup()


saveRDS(dbl_wider, file = "DATA/cliqueOutput/differentiated_genes_wider.RDS")



