
source("GRN/GRN.R")
species <- "spruce"

load("GRN/modules/RData/comparison-Asp-Nor-pearsonMR0.03no-modules2.RData")

regulators_spruce <- get_regulators(species, module.net.genes)

# Filter 
regulators_spruce <- regulators_spruce %>% 
  filter(Type == "TFTF" | 
           (MotifPvalue < 0.05 & 
              ExpressionPvalue < 0.05 &
              Correlation > 0.7)) %>% 
  rename(Node = Node_name) 

write_tsv(regulators_spruce, file = "GRN/spruce/regulators_spruce.tsv")

modules <- regulators_spruce %>% 
  filter(Type == "TFModule") %>% 
  pull(Module) %>% unique()
length(modules)

regulators_spruce %>% 
  filter(Type == "TFModule") %>% 
  pull(Node) %>% unique() %>% length()

n_genes <- 0
for (m in modules) {
  n_genes <- n_genes + length(module.net.genes[[m]])
}
n_genes

# Filter for figure

modules <- c("Module1", "Module2", "Module3", "Module7")

network <- regulators_spruce %>% 
  filter(Type == "TFTF" | 
           (Module %in% modules & Importance >= 3))

write_tsv(network, file = "GRN/spruce/spruce_network.tsv")

TFs <- network %>% 
  select(Node, locus_name, full_name, symbol_all, Type) %>% 
  group_by(Node) %>% 
  slice(1) %>% 
  ungroup() %>%  
  mutate(Type = c("TF"))

mods <- network %>% 
  select(Module, locus_name, full_name, symbol_all, Type) %>% 
  filter(str_detect(Module, "Module")) %>% 
  rename(Node = Module) %>% 
  group_by(Node) %>% 
  slice(1) %>% 
  ungroup() %>%  
  mutate(locus_name = c(NA), full_name = c(NA),
         symbol_all = c(NA),
         Type = c("Module")) 

write_tsv(rbind(TFs, mods), file = "GRN/spruce/spruce_nodes.tsv")



# Plot modules

expr <- read_tsv("NorWood/NorWood_transcriptomics.txt.gz", show_col_types = FALSE) %>% 
  column_to_rownames(var = "Genes")

modules <- c("Module1", "Module2", "Module3", "Module7")
mcol <- c("#B0C4DE", "#A3CC51", "#CC5151", "#51A3CC") # lightsteelblue3
tree <- "N1"
d <- tibble(xintercept = c(8, 14, 22))

for (i in 1:length(modules)) {
  
  central <- module.net %>% 
    filter(Module == modules[i]) %>% 
    pull(CentralGene)
  
  p <- expr %>% 
    rownames_to_column(var = "Genes") %>% 
    filter(Genes %in% central) %>% 
    pivot_longer(-Genes, names_to = "Sample", values_to = "Expression") %>% 
    separate(Sample, into = c("Tree", "Samples"), sep = "\\.") %>% 
    filter(Tree == tree) %>% 
    mutate(Samples = as.numeric(Samples)) %>%
    ggplot(aes(x = Samples, y = Expression, col = Genes)) +
    scale_color_manual(values = mcol[i]) +
    #geom_line(linewidth = 2) +
    geom_smooth(method = "loess", linewidth = 20, se = FALSE) +
    geom_vline(data = d, mapping = aes(xintercept = xintercept), linetype="dashed", 
               color = "gray68", linewidth=2) +
    theme_void() +
    theme(legend.position = "none")
  
  png(paste0("GRN/spruce/", modules[i], ".png"))
  print(p)
  dev.off()
  
}
