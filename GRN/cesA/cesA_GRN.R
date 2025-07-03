
source("GRN/GRN.R")

# cesA
species <- "aspen"
modules <- list()
modules[["cesA"]] <- c("Potra2n2c4078", "Potra2n6c13749", "Potra2n18c32336", "Potra2n11c23199", "Potra2n4c8952")
genes_ids <- c("PtCesA4", "PtCesA7-A", "PtCesA7-B", "PtCesA8-A", "PtCesA8-B")

regulators_cesA_aspen <- get_regulators(species, modules)

# Filter 
regulators_cesA_aspen <- regulators_cesA_aspen %>% 
  filter(Type != "TFTF" & 
           (MotifPvalue < 0.5 & 
              ExpressionPvalue < 0.05 &
              Correlation > 0.7)) %>% 
  rename(Node = Node_name) 

write_tsv(regulators_cesA_aspen, file = "GRN/cesA/regulators_cesA_aspen.tsv")


species <- "spruce"
modules <- list()
modules[["cesA"]] <- c("PA_chr03_G001205", "PA_chr10_G000911", "PA_chr11_G000197")

regulators_cesA_spruce <- get_regulators(species, modules)

# Filter 
regulators_cesA_spruce <- regulators_cesA_spruce %>% 
  filter(Type != "TFTF" & 
           (MotifPvalue < 0.5 & 
              ExpressionPvalue < 0.05 &
              Correlation > 0.7)) %>% 
  rename(Node = Node_name) 

write_tsv(regulators_cesA_spruce, file = "GRN/cesA/regulators_cesA_spruce.tsv")

