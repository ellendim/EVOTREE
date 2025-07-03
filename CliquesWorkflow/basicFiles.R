
# BASIC FILES
# Compiling files to use for further analysis

library(tidyverse)

VER <- "V5"
ORTHOLOG_GROUP_FILE <- "DATA/OrthoGroups/Orthogroups_20240823_clean_N1.tsv.gz" 

# List of all comparison files
species_list <- c("Lodge", "Asp", "Nor","Scots","Birch","Cher")
species_list_full_names <- c("Lodgepole pine", "Aspen", "Norway spruce","Scots pine","Birch","Cherry")
combo <- data.frame(t(combn(species_list, 2)))
col_order <- c(  "Asp-Birch",  "Asp-Cher", "Birch-Cher",  "Lodge-Asp", "Lodge-Birch", "Lodge-Cher",  "Asp-Nor","Nor-Birch", "Nor-Cher", "Scots-Birch", "Scots-Cher",  "Asp-Scots" ,  "Lodge-Scots", "Lodge-Nor", "Nor-Scots")


file_list_1 <- c() 
for (i in 1:nrow(combo)){
  
  s1 <- combo[i, "X1"]
  s2 <- combo[i, "X2"]
  
  file_name <- paste0("DATA/comparisonFiles/comparison-",s1,"_",s2,"-pearsonMR0.03no-table-",VER,".RData")
  file_list_1 <- append(file_list_1, file_name, after = length(file_list_1))
  
}

##  --- all_expressed_genes ---
# Combining all expressed gene pairs from  running ComPlEx into one file.

expr_genes <- c()
for (x in file_list_1){
  if(file.exists(x)){
    
    load(x)
    
    key_word_s1 <- sapply(strsplit(x, "_"), "[",1) 
    key_word_s1 <- sapply(strsplit(key_word_s1, "-"), "[",2) 
    key_word_s2 <- sapply(strsplit(x, "_"), "[",2)
    key_word_s2 <- sapply(strsplit(key_word_s2, "-"), "[",1) 
    
    
    genes <- comparison_table %>%
      as.data.frame() %>% 
      mutate_at("Max.p.val", as.numeric) 
    
    expr_genes <- rbind(expr_genes, data.frame(
      OrthoGroup = genes$OrthoGroup,
      Species1 = rep(key_word_s1, nrow(genes)),
      GeneSpecies1 = genes$Species1,
      Species2 = rep(key_word_s2, nrow(genes)),
      GeneSpecies2 = genes$Species2,
      Species1pVal = genes$Species1.p.val, 
      Species2pVal = genes$Species2.p.val,
      Max.p.Val =  genes$Max.p.val))
  }
}

expr_genes <- expr_genes %>% 
  arrange(OrthoGroup)

save(expr_genes, file = "DATA/all_expressed_genes.RData")

## --- ones_and_zeros_all ---
# Orthogroups with binary input for presence (1) or absence (0) of at least one expressed gene pair (per species pair).

ortholog_group_file <- read.delim(ORTHOLOG_GROUP_FILE, header = TRUE, sep = "\t")

ones_and_zeros_all <- ortholog_group_file %>%
  mutate(Pinus_sylvestris_cp = Pinus_sylvestris) %>% 
  rename(
    Aspen = Populus_tremula,
    Birch = Betula_pendula,
    `Norway spruce` = Picea_abies,
    `Scots pine` = Pinus_sylvestris,
    Cherry = Prunus_avium,
    `Lodgepole pine` = Pinus_sylvestris_cp) %>% 
  select(OrthoGroup, all_of(species_list_full_names))

for (x in file_list_1){
  if (file.exists(x)){ 
    load(x) 
    
    key_word_s1 <- sapply(strsplit(x, "_"), "[",1) 
    key_word_s1 <- sapply(strsplit(key_word_s1, "-"), "[",2) 
    key_word_s2 <- sapply(strsplit(x, "_"), "[",2)
    key_word_s2 <- sapply(strsplit(key_word_s2, "-"), "[",1)
    
    
    df <- comparison_table%>%
      select(Species1, Species2, OrthoGroup, Max.p.val) %>%
      mutate_at("Max.p.val", as.numeric) %>%
      group_by(OrthoGroup) %>%
      slice(1) 
    
    new_col_name <-as.character(paste0(key_word_s1, "-", key_word_s2))
    
    colnames(df)[1] = new_col_name
    
    df <- df %>%
      select(-c(Species2, Max.p.val))
    
    new_column <- ones_and_zeros_all %>%
      left_join(df, join_by(OrthoGroup == OrthoGroup)) %>%
      select(all_of(new_col_name))
    
    ones_and_zeros_all <- cbind(ones_and_zeros_all, new_column)
    
    print(new_col_name)
    
  }}

ones_and_zeros_all <- ones_and_zeros_all %>% 
  column_to_rownames(var = "OrthoGroup") %>% 
  select(-c(1:6))

ones_and_zeros_all[is.na(ones_and_zeros_all)] <- 0
ones_and_zeros_all[ones_and_zeros_all != 0] <-1

ones_and_zeros_all <- ones_and_zeros_all %>%
  mutate_if(is.character, as.integer) %>% 
  select(col_order)

save(ones_and_zeros_all, file = "DATA/ones_and_zeros_all.RData")


# --- Sample clusters ---
# Generates sample clusters (for clique HMs)

for(species in species_list){
  file <- paste0(species,"Wood/", species,"Wood_transcriptomics.txt.gz")
  
  expr.wide <- read.delim(file, sep = "\t") %>% 
    column_to_rownames("Genes")
  
  n <- ncol(expr.wide)
  cdist <- as.dist(1-cor(expr.wide, method="pearson"))
  chc <- hclust(cdist,method="ward.D2")
  cdendro <- as.dendrogram(chc)
  
  order <- rep(0, n)
  h <- colnames(expr.wide)
  z <- strsplit(as.character(h), "[.]")
  h <- unlist(lapply(z,"[",2))
  h <- as.integer(h)
  k <- 1
  for(i in 1:max(ncol(expr.wide))) {
    for(j in 1:n) {
      if (h[j] == i) {
        order[j] = k
        k <- k+1
      }
    }
  }
  cdendro <- reorder(cdendro, order, agglo.FUN = mean)
  chc <- as.hclust(cdendro)
  
  cnclust = 4
  cct <- cutree(chc, k=cnclust)
  colsamples <- c("#A3CC51","#51A3CC","#CC5151","#B26F2C")
  ccol <- colsamples[cct]
  
  saveRDS(cct,paste0("DATA/sampleClusters/cct_", species, ".RDS"))
  
  
}










