library(tidyverse)
library(readxl)
library(cowplot)
library(scales)
library("grid") 
library("gridExtra")

# READ IN EXPRESSION DATA

XWood <- c("AspWood", "BirchWood", "CherWood", "NorWood", "ScotsWood", "LodgeWood")
Tree <- c("A1", "B1", "C1", "N1", "S1", "L1")

XWood_expr <- list()
XWood_stages <- list()

for (i in 1:length(XWood)) {
  
  species <- XWood[i]
  tree <- Tree[i]
  
  XWood_expr[[species]] <- read_tsv(paste0(species, "/", species, "_transcriptomics.txt.gz"), show_col_types = FALSE) %>% 
    select(Genes, contains(tree))
  
  load(paste0(species, "/gene.cluster.RData"))
  XWood_stages[[species]] <- d
}
remove(d, gene.clusters, i, species, tree)

XWood_stages$AspWood$xintercept <- c(6, 11, 20)
XWood_stages$BirchWood$xintercept <- c(5, 9, 22)
XWood_stages$CherWood$xintercept <- c(4, 8, 19)
XWood_stages$NorWood$xintercept <- c(8, 14, 22)
XWood_stages$ScotsWood$xintercept <- c(6, 13, 21)
XWood_stages$LodgeWood$xintercept <- c(4, 12, 23)

# FUCTIONS

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

scaleFUN <- function(x) sprintf("%2.0f", x)

plot_cliques <- function(Cliques_all, cliques, clique_names, XWood, 
                         color_palette = gg_color_hue(length(cliques))) {
  
  species_name <- c("P. tremula", "B. pendula", "P. avium", "P. abies", "P. sylvestris", "P. contorta")
  
  plots <- list()
  for (j in 1:length(XWood)) {
    
    species1 <- XWood[j]
    species2 <- gsub("Wood", "", species1)
    
    ids <- Cliques_all %>% 
      filter(cliqueID %in% cliques) %>% 
      mutate(cliqueID = factor(cliqueID, levels = cliques)) %>% 
      arrange(cliqueID) %>% 
      pull(species2)
    
    d <- XWood_stages[[species1]]
    
    plots[[length(plots)+1]] <- XWood_expr[[species1]] %>% 
      filter(Genes %in% ids) %>% 
      pivot_longer(-Genes, names_to = "Sample", values_to = "Expression") %>% 
      separate(Sample, into = c("Tree", "Samples"), sep = "\\.") %>% 
      mutate(Samples = as.numeric(Samples)) %>%
      mutate(Genes = factor(Genes, levels = ids)) %>% 
      ggplot(aes(x = Samples, y = Expression, col = Genes)) +
      geom_line(linewidth = 2) +
      #geom_smooth(method = "loess", linewidth = 2, se = FALSE) +
      scale_x_continuous(breaks = c(0, 10, 20)) +
      scale_y_continuous(labels=scaleFUN) +
      #ylim(0,14) +
      scale_color_manual(values = color_palette) +
      ggtitle(species_name[j]) + 
      geom_vline(data = d, mapping = aes(xintercept = xintercept), linetype="dashed", color = "gray68", linewidth=0.5) +
      theme_bw() +
      theme(legend.position = "none", plot.title = element_text(face = "italic"),
            text=element_text(size=20),
            axis.title.x=element_blank(), axis.title.y=element_blank())
      #theme(legend.position = "bottom", legend.title=element_blank())
    
  }
  
  clique_names <- unique(clique_names)
  color_palette <- unique(color_palette)
  
  # Create a dummy plot just for the legend
  p_dummy <- ggplot(data.frame(x = 1:length(clique_names),
                               y = 1:length(clique_names), 
                               names = factor(clique_names, levels = clique_names)), 
                    aes(x, y, color = names, group = 1)) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = color_palette) +
    theme_bw() +
    theme(text=element_text(size=18)) +
    guides(color = guide_legend(direction = "horizontal", title = NULL))
  
  legend <- cowplot::get_legend(p_dummy)
  
  title2 <- ggdraw() + draw_label("Expression (VST)", fontface='plain', size = 20, angle = 90)
  title3 <- ggdraw() + draw_label("Samples", fontface='plain', size = 20)
  
  p <- plot_grid(plotlist = plots, nrow = 1) 
  p <- plot_grid(legend, p, ncol = 1, rel_heights=c(0.1, 0.9))
  p <- plot_grid(p, title3, ncol = 1, rel_heights=c(0.92, 0.08))
  
  plot_grid(title2, p, ncol=2, rel_widths = c(0.025,0.975))
  
}


plot_orthogroups <- function(ortho, OrthoGroups, names, XWood) {
  
  species_name <- c("Populus_tremula", "Betula_pendula", "Prunus_avium", "Picea_abies", "Pinus_sylvestris", "Pinus_contorta")
  species_name2 <- c("P. tremula", "B. pendula", "P. avium", "P. abies", "P. sylvestris", "P. contorta")
  
  all_plots <- list()
  for (i in 1:length(OrthoGroups)) {
    
    og <- OrthoGroups[i]
    
    descr <- og
    
    plots <- list()
    for (j in 1:length(XWood)) {
      
      species1 <- XWood[j]
      species2 <- gsub("Wood", "", species1)
      
      ids <- ortho %>% 
        filter(OrthoGroup %in% og) %>% 
        separate_rows(sym(species_name[j]), sep = ", ") %>% 
        pull(species_name[j]) %>% 
        unique()
      
      d <- XWood_stages[[species1]]
      
      plots[[length(plots)+1]] <- XWood_expr[[species1]] %>% 
        filter(Genes %in% ids) %>% 
        pivot_longer(-Genes, names_to = "Sample", values_to = "Expression") %>% 
        separate(Sample, into = c("Tree", "Samples"), sep = "\\.") %>% 
        mutate(Samples = as.numeric(Samples)) %>% 
        ggplot(aes(x = Samples, y = Expression, col = Genes)) +
        geom_line(linewidth = 2) + # xlab("Samples") + ylab("Expression (VST)") +
        scale_x_continuous(breaks = c(0, 10, 20)) +
        scale_y_continuous(labels=scaleFUN) +
        ggtitle(species_name2[j],) +
        geom_vline(data = d, mapping = aes(xintercept = xintercept), linetype="dashed", 
                   color = "gray68", linewidth=0.5) +
        theme_bw() +
        theme(#legend.position = "none", 
              legend.direction = "vertical", legend.position = "bottom",  legend.title=element_blank(),
              plot.title = element_text(face = "italic"), 
              text=element_text(size=20),
              axis.title.x=element_blank(), axis.title.y=element_blank())
      
    }
    
    title1 <- ggdraw() + draw_label(names[i], fontface='bold', size = 20)
    title2 <- ggdraw() + draw_label("Expression (VST)", fontface='plain', size = 20, angle = 90)
    title3 <- ggdraw() + draw_label("Samples", fontface='plain', size = 20)
    
    p <- plot_grid(plotlist = plots, nrow = 1) 
    p <- plot_grid(title1, p, ncol = 1, rel_heights=c(0.1, 0.9))
    p <- plot_grid(p, title3, ncol = 1, rel_heights=c(0.92, 0.08))
    
    all_plots[[i]] <- plot_grid(title2, p, ncol=2, rel_widths = c(0.025,0.975))
    
  }
  all_plots
}

plot_cliques_ind <- function(cliques, clique_names, XWood) {
  
  species_name <- c("Populus_tremula", "Betula_pendula", "Prunus_avium", "Picea_abies", "Pinus_sylvestris", "Pinus_sylvestris")
  
  all_plots <- list()
  for (i in 1:length(cliques)) {
    
    clique <- cliques[i]
    
    descr <- clique
    
    plots <- list()
    for (j in 1:length(XWood)) {
      
      species1 <- XWood[j]
      species2 <- gsub("Wood", "", species1)
      
      ids <- Cliques %>% 
        filter(cliqueID %in% clique) %>% 
        arrange(cliqueID) %>% 
        pull(species2)
      
      d <- XWood_stages[[species1]]
      
      plots[[length(plots)+1]] <- XWood_expr[[species1]] %>% 
        filter(Genes %in% ids) %>% 
        pivot_longer(-Genes, names_to = "Sample", values_to = "Expression") %>% 
        separate(Sample, into = c("Tree", "Samples"), sep = "\\.") %>% 
        mutate(Samples = as.numeric(Samples)) %>% 
        ggplot(aes(x = Samples, y = Expression, col = Genes)) +
        geom_line(linewidth = 2) +
        ggtitle(species2) + ylab("Expression (VST)") +
        geom_vline(data = d, mapping = aes(xintercept = xintercept), linetype="dashed", 
                   color = "gray68", linewidth=0.5) +
        theme_bw() +
        #theme(legend.position = "none", plot.title = element_text(face = "italic"))
        theme(legend.position = "bottom", legend.title=element_blank(), 
              plot.title = element_text(face = "italic"))
      
    }
    
    title1 <- ggdraw() + draw_label(clique, fontface='bold')
    
    title2 <- ggdraw() + draw_label(clique_names[i], fontface='plain', size = 8)
    
    p <- plot_grid(plotlist = plots, nrow = 1) 
    
    all_plots[[i]] <- plot_grid(title1, title2, p, ncol = 1, rel_heights=c(0.06, 0.04, 1))
    
  }
  all_plots
}
