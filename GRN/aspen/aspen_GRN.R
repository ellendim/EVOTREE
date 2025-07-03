
source("GRN/GRN.R")
species <- "aspen"

load("GRN/modules/RData/comparison-Asp-Nor-pearsonMR0.03no-modules1.RData")

regulators_aspen <- get_regulators(species, module.net.genes)

# Filter 
regulators_aspen <- regulators_aspen %>% 
  filter(Type == "TFTF" | 
           (MotifPvalue < 0.05 & 
           ExpressionPvalue < 0.05 &
           Correlation > 0.7)) %>% 
  rename(Node = Node_name) 

write_tsv(regulators_aspen, file = "GRN/aspen/regulators_aspen.tsv")

modules <- regulators_aspen %>% 
  filter(Type == "TFModule") %>% 
  pull(Module) %>% unique()
length(modules)

regulators_aspen %>% 
  filter(Type == "TFModule") %>% 
  pull(Node) %>% unique() %>% length()

n_genes <- 0
for (m in modules) {
  n_genes <- n_genes + length(module.net.genes[[m]])
}
n_genes

# Filter for figure

modules <- c("Module1", "Module5", "Module6", "Module9", "Module11")

network <- regulators_aspen %>% 
  filter(Type == "TFTF" | 
           (Module %in% modules & Importance >= 3))

write_tsv(network, file = "GRN/aspen/aspen_network.tsv")

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

write_tsv(rbind(TFs, mods), file = "GRN/aspen/aspen_nodes.tsv")



# Plot modules

expr <- read_tsv("AspWood/AspWood_transcriptomics.txt.gz", show_col_types = FALSE) %>% 
  column_to_rownames(var = "Genes")

modules <- c("Module1", "Module5", "Module6", "Module9", "Module11")
mcol <- c("#A3CC51", "#B0C4DE", "#B26F2C", "#51A3CC", "#CC5151") # lightsteelblue3
tree <- "A1"
d <- tibble(xintercept = c(6, 11, 20))

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
  
  png(paste0("GRN/aspen/", modules[i], ".png"))
  print(p)
  dev.off()
  
}
