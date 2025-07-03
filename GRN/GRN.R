library(tidyverse)

path <- ""

#############

# Converting correlation to p-value
p_from_r <- function(r, n) {
  t <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t), df = n - 2)
}

# Link motifs to Arabidopsis genes

get_motif_at <- function() {

  # Get Arabidopsis IDs
  
  sql1 <- read_tsv(paste0(path, "DATA/GRN/JASPAR2024.sql.gz"), col_names = FALSE, show_col_types = FALSE) %>% 
    filter(str_detect(X1, "INSERT INTO MATRIX VALUES")) %>% 
    mutate(X1 = gsub("INSERT INTO MATRIX VALUES\\(", "", X1)) %>% 
    mutate(X1 = gsub("\\);", "", X1)) %>% 
    separate(X1, into = c("ID", NA, NA, NA, "Motif"), sep = ",") %>% 
    mutate(Motif = gsub("\\'", "", Motif))
  
  sql2 <- read_tsv(paste0(path, "DATA/GRN/JASPAR2024.sql.gz"), col_names = FALSE, show_col_types = FALSE) %>% 
    filter(str_detect(X1, "INSERT INTO MATRIX_PROTEIN VALUES")) %>% 
    mutate(X1 = gsub("INSERT INTO MATRIX_PROTEIN VALUES\\(", "", X1)) %>% 
    mutate(X1 = gsub("\\);", "", X1)) %>% 
    separate(X1, into = c("ID", "UNIPROT"), sep = ",") %>% 
    mutate(UNIPROT = gsub("\\'", "", UNIPROT)) %>% 
    filter(UNIPROT != "")
  
  # BLAST
  if (FALSE) {
    inner_join(sql1, sql2) %>% 
      distinct(Motif, UNIPROT) %>% 
      filter(UNIPROT != "") %>% 
      select(UNIPROT) %>% 
      write_tsv("BLAST/uniprot.txt")
  }
  
  UNIPROT <- read_tsv(paste0(path, "DATA/GRN/Uniprot2AGI.txt"), col_names = FALSE, show_col_types = FALSE) %>% 
    separate_rows(X2, sep = "\\;") %>% 
    separate(X2, into = ("X2"), sep = "\\.", extra = "drop") %>% 
    distinct(X2, .keep_all = TRUE) %>% 
    mutate(X2 = toupper(X2)) %>% 
    rename(AtUNIPROT = X2, UNIPROT = X1)
  
  BLAST <- read_tsv(paste0(path, "DATA/GRN/BLAST.tsv"), col_names = FALSE, show_col_types = FALSE) %>% 
    separate(X1, into = c(NA,"X1",NA), sep = "\\|") %>% 
    separate(X2, into = c("X2",NA), sep = "\\.") %>% 
    group_by(X1) %>% 
    arrange(desc(X12)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(X2 = toupper(X2)) %>% 
    select(X1, X2, X12, X11, X3) %>% 
    rename(AtBLAST = X2, UNIPROT = X1, B = X12, E = X11, ID = X3)
  
  sql <- left_join(sql1, sql2) %>% 
    select(-ID) %>% 
    distinct(Motif, UNIPROT) %>% 
    left_join(BLAST) %>%
    left_join(UNIPROT) %>%
    select(-UNIPROT) %>% 
    left_join(read_tsv(paste0(path, "DATA/Annotations/gene_aliases_20140331.txt"), show_col_types = FALSE) %>% 
                filter(str_detect(locus_name, "^AT\\dG|^At\\dg")),
              by = c("Motif"="symbol"), relationship = "many-to-many") %>% 
    mutate(At = if_else(str_detect(Motif, "^AT\\dG|^At\\dg"), toupper(Motif), NA)) %>% 
    filter(!(is.na(AtBLAST) & is.na(AtUNIPROT) & is.na(locus_name) & is.na(At))) %>% 
    distinct(Motif, AtBLAST, AtUNIPROT, locus_name, At, .keep_all = TRUE) %>% 
    rowwise() %>% 
    mutate(At2 = paste(unique(na.omit(c(AtBLAST, AtUNIPROT, locus_name))), collapse = "-"),
           Consistent = length(unique(na.omit(c(AtBLAST, AtUNIPROT, locus_name)))) == 1) %>% 
    mutate(Consistent = if_else(is.na(At), Consistent, TRUE)) %>% 
    group_by(Motif) %>% 
    mutate(n = n()) %>% 
    filter(!(!Consistent & n > 1)) %>% 
    filter(!(!Consistent & E > 1E-100)) %>% 
    mutate(At2 = ifelse(!Consistent, AtBLAST, At2)) %>% 
    mutate(Consistent = if_else(E > 1E-100 & is.na(AtUNIPROT) & is.na(locus_name), FALSE, Consistent)) %>% 
    mutate(At = if_else(is.na(At), At2, At)) %>% 
    group_by(Motif) %>% 
    mutate(Consistent_down = if_else(length(unique(At)) > 1, FALSE, TRUE)) %>% 
    filter(Consistent_down) %>% 
    ungroup()
  
  #write_tsv(sql, "motif_at.tsv")
  
  motif_at <- sql %>% 
    select(Motif, At) %>% 
    left_join(read_tsv(paste0(path, "DATA/Annotations/gene_aliases_20140331.txt"), show_col_types = FALSE) %>% 
              group_by(locus_name) %>% 
              summarise(full_name = full_name[1],
                        symbol_all = paste0(symbol, collapse = "-")),
            by = c("At"="locus_name")) %>% 
    rename(locus_name = At)
      
  
  # Manual
  #unique(fimo$motif_alt_id)[(!unique(fimo$motif_alt_id) %in% motif_at$Motif)]
  motif_at <- rbind(motif_at,
                    tibble(Motif = c("NAC066","NAC030","NAC012"),
                           locus_name = c("AT3G61910", "AT1G71930", "AT1G32770"),
                           full_name = c(NA),
                           symbol_all = c(NA)))

  motif_at <- motif_at %>% 
    left_join(read_tsv(paste0(path, "DATA/GRN/nodes.txt"), show_col_types = FALSE) %>% 
                select(Node, Node_name, Type, Importance),
              by = c("Motif"="Node"))

  idx <- which(motif_at$Motif == "NAC030")
  motif_at$Node_name[idx] <- "VND7"
  motif_at$Importance[idx] <- 4
  
  idx <- which(motif_at$Motif == "NAC012")
  motif_at$Node_name[idx] <- "SND1"
  motif_at$Importance[idx] <- 4
    
  # Filter 
  motif_at <- motif_at %>% 
    drop_na(locus_name)
  
  # Fill in
  motif_at <- motif_at %>%
    mutate(Importance = if_else(is.na(Importance), 1, Importance)) %>% 
    mutate(Node_name = if_else(is.na(Node_name), Motif, Node_name)) %>% 
    mutate(Type = "TF")
  
  motif_at
}

get_regulators <- function(species, modules) {
  
  motif_at <- get_motif_at()
  
  ortho <- read_tsv(paste0(path, "DATA/OrthoGroups/Orthogroups_20240823_clean_N1.tsv.gz"), show_col_types = FALSE) %>% 
    mutate(Pinus_contorta = Pinus_sylvestris) %>% 
    separate_rows(Arabidopsis_thaliana, sep = ", ")
  
  # Enriched motifs in modules
  if (species == "aspen") {
    
    expr <- read_tsv(paste0(path, "AspWood/AspWood_transcriptomics.txt.gz"), show_col_types = FALSE) %>% 
      column_to_rownames(var = "Genes")
    
    fimo <- read_tsv(paste0(path, "DATA/GRN/fimo_aspen_all.tsv.gz"), show_col_types = FALSE) %>%
      separate(sequence_name, into = c("sequence_name"), sep = "_", extra = "drop") %>% 
      filter(sequence_name %in% rownames(expr))
  } else {
    
    expr <- read_tsv(paste0(path, "NorWood/NorWood_transcriptomics.txt.gz"), show_col_types = FALSE) %>% 
      column_to_rownames(var = "Genes")
    
    fimo <- read_tsv(paste0(path, "DATA/GRN/fimo_spruce_all.tsv.gz"), show_col_types = FALSE) %>%
      mutate(sequence_name = paste0("PA_", sequence_name)) %>% 
      filter(sequence_name %in% rownames(expr))
  }
  
  motifs <- unique(fimo$motif_alt_id)
  
  N <- nrow(expr)
  
  module_motif <- tibble()
  for (motif in sort(unique(fimo$motif_alt_id))) {
    
    #cat(motif, "\n")
    
    genes_motifs <- fimo %>% 
      filter(motif_alt_id == motif) %>% 
      pull(sequence_name) %>% 
      unique()
    
    k <- length(genes_motifs)
    
    for (module in sort(names(modules))) {
      
      genes_module <- modules[[module]] %>% unique()
      m <- length(genes_module)
      n <- N-m
      
      x <- length(intersect(genes_module, genes_motifs))
      
      if (m > 5) {
      
        p_val <- 1
        effect_size <- NA
        if (x > 1) {
          p_val <- phyper(x-1, m, n, k, lower.tail = FALSE)
          
          effect_size <- (x/m) / (k/N)
        }
      } else {
        p_val <- 1-x/m
        if (p_val == 0) {
          p_val = 0.1
        }
        effect_size <- x/m
      }
      
      #if (p_val < 0.05) {
        module_motif <- rbind(module_motif,
                              tibble(Module = c(module), Motif = c(motif), 
                                     MotifPvalue = c(p_val), MotifEffectSize = c(effect_size)))
        
      #}
    }
  }
  
  ########
  # Correlation regulator - modules
  
  module_motif <- module_motif %>% 
    left_join(motif_at, relationship = "many-to-many") %>% 
    drop_na(locus_name)
  
  module_motif$Gene <- c(NA)
  module_motif$Correlation <- c(NA)
  module_motif$ExpressionPvalue <- c(NA)
  for(i in 1: nrow(module_motif)) {
    
    module <- module_motif$Module[i]
    
    genes_module <- modules[[module]] %>% unique()
    center <- colSums(expr[genes_module,])/length(genes_module)
    
    if (!is.na(module_motif$locus_name[i])) {
      
      if (species == "aspen") {
        genes <- ortho %>% 
          filter(Arabidopsis_thaliana == module_motif$locus_name[i]) %>% 
          separate_rows(Populus_tremula, sep = ", ") %>%
          filter(Populus_tremula %in% rownames(expr)) %>% # Populus_tremula
          pull(Populus_tremula) 
      } else {
        genes <- ortho %>% 
          filter(Arabidopsis_thaliana == module_motif$locus_name[i]) %>% 
          separate_rows(Picea_abies, sep = ", ") %>%
          filter(Picea_abies %in% rownames(expr)) %>% # Populus_tremula
          pull(Picea_abies) 
      }
      
      if (length(genes) > 0) {
        
        R <- cor(center, t(expr[genes,]))
        idx <- which.max(R)
        gene <- colnames(R)[idx]
        r <- R[idx]
        
        p <- p_from_r(r, 25)
        
        module_motif$Gene[i] <- gene
        module_motif$Correlation[i] <- r
        module_motif$ExpressionPvalue[i] <- p
      }
    }
  }
  
  # Filter unmapped motifs
  module_motif <- module_motif %>% 
    drop_na(Correlation)
  
  # Correct p-values
  if (!(length(modules) == 1 & length(modules[[names(modules)[1]]]) <= 5)) {
    module_motif <- module_motif %>% 
      group_by(Module) %>% 
      mutate(MotifPvalue = p.adjust(MotifPvalue, method = "fdr")) %>% 
      mutate(ExpressionPvalue = p.adjust(ExpressionPvalue, method = "fdr")) %>% 
      ungroup()
  }
  
  # Compute score
  module_motif <- module_motif %>% 
    mutate(Score = (-log10(ExpressionPvalue)/max(-log10(ExpressionPvalue))) * (-log10(MotifPvalue)/max(-log10(MotifPvalue))))
  
  # Filter based on score!
  module_motif <- module_motif %>% 
    arrange(desc(Score)) %>% 
    filter(Score > 0)
  
  module_motif <- module_motif %>% mutate(Type = c("TFModule"))
  
  ######
  # TF-TF

  TFTF <- fimo %>%
    filter(motif_alt_id %in% module_motif$Motif & sequence_name %in% module_motif$Gene) %>%
    group_by(motif_alt_id, sequence_name) %>%
    arrange(`p-value`) %>%
    slice(1) %>%
    ungroup() %>% 
    left_join(module_motif, by = c("sequence_name"="Gene"), relationship = "many-to-many")

  layer1 <- module_motif %>%
    filter(str_detect(Node_name, "VND|NST|BRN|SND1")) %>%
    pull(Motif) %>% 
    unique()

  layer2 <- c("MYB46", "MYB83")

  layer3 <- module_motif %>%
    filter(Importance >= 3) %>% 
    filter(str_detect(Node_name, "^MYB|SND2|SND3")) %>%
    filter(!(Node_name %in% c(layer2, "MYB107"))) %>%
    filter(!str_detect(Node_name, "R")) %>%
    pull(Motif) %>% 
    unique()

  TFTFselect <- TFTF %>%
    filter((motif_alt_id %in% layer1 & Motif %in% layer2) |
             (motif_alt_id %in% layer2 & Motif %in% layer3)) %>%
    select(motif_alt_id, Motif) %>%
    rename(Module = Motif) %>%
    rename(Motif = motif_alt_id) %>%
    select(Module, Motif) %>%
    group_by(Module, Motif) %>%
    slice(1) %>%
    ungroup() %>% 
    left_join(module_motif %>% 
                select(Motif, Node_name) %>% 
                group_by(Motif) %>% 
                slice(1) %>%
                ungroup()) %>% 
    mutate(MotifPvalue = c(1.0), MotifEffectSize = c(1.0), locus_name = c("NA"),
           full_name = c("NA"), symbol_all = c("NA"), 
           Type = c("TFTF"),
           Importance = c(4), Gene = c("NA"), 
           Correlation = c(0), ExpressionPvalue = c(1), 
           Score = c(0))

  TFTFselect$Module <- gsub("NAC010", "SND3", TFTFselect$Module) # SND3
  TFTFselect$Module <- gsub("NAC073", "SND2", TFTFselect$Module) # SND2

  # Add TF-TF links
  module_motif <- rbind(module_motif, TFTFselect)
  
  # Return
  module_motif
  
}
