
source("Stories/plotfuns.R")

# Read OG
OGs <- read_tsv("DATA/OrthoGroups/Orthogroups_20240823_clean.tsv.gz") %>% 
  mutate(Pinus_contorta = Pinus_sylvestris) %>% 
  separate_rows(Arabidopsis_thaliana, sep = ", ")

HOGs <- read_tsv("DATA/OrthoGroups/Orthogroups_20240823_clean_N1.tsv.gz") %>% 
  mutate(Pinus_contorta = Pinus_sylvestris) %>% 
  separate_rows(Arabidopsis_thaliana, sep = ", ")

#######

# Marker genes - Fig 1B

Cliques <- read_tsv("Stories/marker_genes.tsv")
palette_1 <- c("#A3CC51", "mediumpurple3","#51A3CC", "#CC5151",  "#B26F2C")

p <- plot_cliques(Cliques, Cliques$cliqueID, Cliques$cliqueID, XWood, palette_1)
p

pdf("Stories/marker_genes.pdf", width = 20, height = 5)
p
dev.off()

#######

# POK2 - Fig. 4B

Cliques <- read_excel("Supplement/TabS4A_Clique_genes.xlsx", sheet = "Differentiated genes")

hog <- HOGs %>% 
  filter(Arabidopsis_thaliana == "AT3G19050") %>% 
  pull(OrthoGroup)

cliques <- Cliques %>% 
  filter(OrthoGroup == hog) %>% 
  pull(cliqueID)

names <- c("POK2")
palette_1 <- c("#FF8F63", "#FF8F63", "#FF8F63", "#66C2A5", "#66C2A5", "#66C2A5") 

p <- plot_cliques(Cliques, cliques, names, XWood, palette_1)
p

pdf("Stories/POK2.pdf", width = 20, height = 5)
p
dev.off()

#######

# NSTs - Fig. 6B

NSTs <- read_tsv("DATA/Annotations/gene_aliases_20140331.txt", show_col_types = FALSE) %>% 
  filter(str_detect(symbol, "^NST\\d")) %>% 
  dplyr::rename(GeneID = "locus_name", Name = symbol) %>%
  left_join(HOGs %>% select(OrthoGroup, Arabidopsis_thaliana, Populus_tremula, Betula_pendula,
                             Prunus_avium, Picea_abies, Pinus_sylvestris),
            by = c("GeneID"="Arabidopsis_thaliana"))

hogs <- HOGs %>% 
  filter(Arabidopsis_thaliana %in% NSTs$GeneID) %>% 
  distinct(OrthoGroup) %>% 
  pull(OrthoGroup)

names <- "NST"

p <- plot_orthogroups(HOGs, hogs, names, XWood)
p

pdf("Stories/NST.pdf", width = 20, height = 5)
p
dev.off()

