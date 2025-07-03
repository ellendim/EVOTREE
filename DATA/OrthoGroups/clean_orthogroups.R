
library(tidyverse)


#ortho <- read_tsv("Orthogroups_20240823/Orthogroups/Orthogroups.tsv", show_col_types = FALSE)
ortho <- read_tsv("Orthogroups_20240823/Phylogenetic_Hierarchical_Orthogroups/N1.tsv", show_col_types = FALSE)

arab <- read_tsv("../Annotations/Athaliana_167_annotation_info.txt.gz", 
                 show_col_types = FALSE, col_names = FALSE) %>% 
  select(X1, X2) %>% 
  dplyr::rename(Arabidopsis_thaliana = X1, At = X2) %>% 
  mutate(Arabidopsis_thaliana = as.character(Arabidopsis_thaliana))

arab <- ortho %>% 
  select(HOG, Arabidopsis_thaliana) %>% 
  separate_rows(Arabidopsis_thaliana, sep = ", ", convert = FALSE) %>% 
  mutate(Arabidopsis_thaliana = gsub("PAC_", "", Arabidopsis_thaliana)) %>% 
  left_join(arab) %>% 
  select(-Arabidopsis_thaliana) %>% 
  group_by(HOG) %>% 
  summarise(At = paste(At, collapse = ", "))

ortho_new <- ortho %>% 
  left_join(arab) %>% 
  mutate(Arabidopsis_thaliana = At) %>% 
  select(-At) %>% 
  mutate(Populus_tremula = gsub("\\.\\d+", "", Populus_tremula)) %>% 
  mutate(Prunus_avium = Prunus_avium_tieton) %>% 
  select(-Prunus_avium_tieton) %>% 
  mutate(Picea_abies = gsub("\\.mRNA\\.\\d+", "", Picea_abies)) %>% 
  mutate(Pinus_sylvestris = gsub("\\.mRNA\\.\\d+", "", Pinus_sylvestris)) %>%
  dplyr::rename(OrthoGroup = HOG) %>%
  write_tsv("Orthogroups_20240823_clean_N1.tsv.gz")


            