# EVOTREE

Git repository for the Wood Regulomics article.

To run the scripts, first download DATA.tar.gz (link to figshare!!) containing the input data. 
Scripts in this repo assumes there is a folder named *DATA* in the project root.
All scripts (except .Rmd) should be run with project root as working directory such that data can be loaded with the relative path `DATA/`.



## Folders and files: 


* **XWood (X in {Asp, Birch, Cher, Nor, Scots, Lodge}):** There is a folder per species containing an 
.Rmd file for generating the VST-normalized expression data used in all downstream analysis 
e.g. ScotsWood/ScotsWood.Rmd. These scripts also generate heatmaps, PCAs and example plots.

* **CliquesWorkflow:** Contains all code for identifying cliques (see **Identifying cliques**).
  Additional .Rmd files generate upset plots, sample comparison heatmap, as well as GO enrichment and heatmaps of conserved cliques. The bestUniqueCliques.R script was only used for generating supplementary .xlsx files.

* **GRN:** Scripts for generating the gene regulatory networks (GRNs) for aspen and spruce. 
Run RComPlEx.Rmd and then ModuleNetwork.Rmd in the modules-folder to create modules.
Then run either aspen/aspen_GRN.R or spruce/spruce_GRN.R to create the GRNs.

* **Stories:** Script for generating expression plots.

* **Supplement:** Tables from the supplementary material used to make expression plots.


## Identifying cliques:

**Terminology:**

* *Size of clique* - refers to the number of species included in clique. In this work, size can be six-member (6M), five-member (5M), or three-member (3M).

* *Complete clique* - cliques where all edges have FDR < 0.1 (6M, 5M, 3M).

* *Partial clique* - cliques were we permit some edges to have FDR > 0.1 (6M).

* *Conserved genes* - complete 6M cliques (conserved across all species) OR 3M complete cliques (conserved within lineage).

* *Partially conserved genes* - partial 6M cliques (11-14 edges have FDR < 0.1, the remaining have FDR < 0.9) + complete 5M cliques (gene from one species is missing).

* *Differentiated genes* - partial 6M cliques where all edges between lineage-specific species but no more than four inter-lineage connections have FDR < 0.1. 

**Workflow**

Note: both standard orthogroups (OGs) and Hierarchical Orthogroups (HOGs) with higher resolution were used. Comparison-files under `DATA/comparisonFiles` generated using HOGs are denoted "V5"  while files generated using OGs are denoted "V4". The workflow is the same regardless of type of resolution, however, note that the size of OGs increases memory requirements and runtime significantly. 

1. Run COMPLEX_script.R for all 15 pairs of species. Produces comparison files.
2. Run basicFiles.R, a script for compiling files for downstream analysis.
3. Run the clique algorithms in order explained below. To avoid overlap in orthogroups (OGs) between the different sets conserved/differentiated genes (e.g. a gene should not be classified as both partially conserved, and as differentiated), the clique
   scripts must be run in the following order:
   
   i. **CompleteCliques.R** - identifies genes conserved across all species (6M) and lineage-specific genes (3M). The latter gene sets (i.e. conifer- and dicot-specific) are identified from isolating the OGs with genes only present in the respective lineages.
   
   ii. **PartialCliques.R** - identifies partially conserved genes in two batches: the code starts with removing OGs with complete 6M cliques in order to identity partial 6M cliques (partially significant), then identifies 5M complete cliques (partially present) from
   CompleteCliquesAlgo.R

   iii. **DifferentiatedCliques.R** - removes OGs with conserved, partially conserved, and lineage-specific genes. Identifies cliques with co-expressologs between species within the lineages and between no more than four of the inter-lineage pairs.


