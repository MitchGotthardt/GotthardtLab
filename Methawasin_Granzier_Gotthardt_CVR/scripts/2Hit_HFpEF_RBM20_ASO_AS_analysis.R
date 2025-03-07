################################################################################

# Manuscript: Rbm20 antisense oligonucleotides alleviate diastolic dysfunction in a mouse model of cardiometabolic HFpEF.
# Authors: Mei Methawasin*, Stefan Meinke, Michael H Radke, Gerrie P Farman, Zaynab Hourani, John E. Smith III, Wei Guo, Henk Granzier*, Michael Gotthardt*
# *Corresponding author
# Bioinformatics and data analysis: Stefan Meinke

################################################################################

################################################################################
####################### Alternative Splicing Analysis ##########################
################################################################################

required_libraries <- c("BiocManager",
                        "readr", 
                        "tidyverse", 
                        "magrittr", 
                        "ggplot2", 
                        "dplyr", 
                        "readxl",
                        "eulerr",
                        "ggtranscript",
                        "patchwork",
                        "writexl",
                        "openxlsx",
                        "stringr",
                        "purrr",
                        "UpSetR",
                        "rtracklayer",
                        "org.Mm.eg.db",
                        "GenomicFeatures")

# Check and install/load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


for (package in required_libraries) {
  if (package == "ggtranscript") {
    # Install ggtranscript from GitHub
    devtools::install_github("dzhang32/ggtranscript")
  } else {
    # Check if the package is already installed
    if (!requireNamespace(package, quietly = TRUE)) {
      # Install Bioconductor packages
      if (package %in% BiocManager::available()) {
        BiocManager::install(package)
      } else {
        # Install CRAN packages
        install.packages(package, dependencies = TRUE)
      }
    }
  }
  # Load the package
  library(package, character.only = TRUE)
}



################# load and prepare the alternative Splicing (AS) dataframe ###################

dir <- "path/to/*.MATS.JC.txt/files" 

### Load rMATS output files (JC)
# Get a list of subdirectories in the main directory
subdirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store data frames
dfs_list <- list()

# Iterate through each subdirectory
for (subdir in subdirs) {
  # Get a list of files in the current subdirectory
  files_list <- list.files(path = subdir, pattern = "*.MATS.JC.txt", full.names = TRUE, recursive = TRUE)
  
  # Get the subdirectory name
  subdirectory_name <- basename(subdir)
  
  # Import data frames from files in the current subdirectory
  dfs <- purrr::map(files_list, ~ read_tsv(.x, col_names = TRUE) %>% 
                      mutate(Comparison = subdirectory_name) %>%
                      mutate(Type = gsub("\\.MATS\\.JC\\.txt$", "", basename(.x))))
  
  # Bind the 5 data frames into one within each subdirectory
  dfs_combined <- bind_rows(dfs)
  
  # Append the combined data frame to dfs_list with subdirectory name as the entry
  dfs_list[[subdirectory_name]] <- dfs_combined
}

# Combine all data frames into a single df
AS_all <- bind_rows(dfs_list)

AS_all_filtered <- AS_all %>%
  filter(FDR < 0.05 & abs(IncLevelDifference) > 0.1) %>% 
  dplyr::select(-c(ID...1, ID...12, ID...34))


gene_ID <- AS_all %>% 
  dplyr::select(GeneID, Gene = geneSymbol) %>% 
  unique





################################ Upset Plot ####################################

#################
### Figure 5A ###
#################
Comparisons <- c("2Hit ASO_vs_2Hit PBS", "2Hit ASO_vs_Ctrl PBS",
                 "2Hit PBS_vs_Ctrl PBS", "Ctrl ASO_vs_Ctrl PBS")

upset_list <- lapply(Comparisons, function(comp) {
  AS_all_filtered %>%
    filter(Comparison == comp,
           Type == "SE") %>%
    summarize(geneSymbol = geneSymbol %>% unique()) %>%
    pull(geneSymbol)
})

names(upset_list) <- Comparisons

Figure_5A <- upset(fromList(upset_list),
      order.by = "freq",
      sets.x.label = "number of genes",
      text.scale = 1.5,
      point.size = 2,
      mb.ratio=c(0.7,0.3))





all_genes <- unique(unlist(upset_list))
all_genes <- all_genes[!is.na(all_genes)]

upset_data <- as.data.frame(
  sapply(upset_list, function(x) all_genes %in% x)
)
row.names(upset_data) <- all_genes

# Extract intersections using a custom function
extract_intersections <- function(upset_data) {
  comb_matrix <- as.matrix(upset_data)
  set_names <- colnames(comb_matrix)
  intersect_list <- list()
  
  # Iterate through each combination
  for (i in 1:nrow(comb_matrix)) {
    included_sets <- set_names[comb_matrix[i, ] == TRUE]
    intersection_name <- paste(included_sets, collapse = " & ")
    intersect_list[[intersection_name]] <- c(intersect_list[[intersection_name]], rownames(comb_matrix)[i])
  }
  
  return(intersect_list)
}

# Get intersections
intersect_genes <- extract_intersections(upset_data)


Datafile_S2 <- data.frame(
  Intersection = names(intersect_genes),
  Genes = sapply(intersect_genes, function(x) paste(x, collapse = ", ")),
  Count = sapply(intersect_genes, length)
)



Datafile_S1 <- data.frame(
  `Grouped Comparison` = Comparisons,
  Genes = sapply(Comparisons, function(x) AS_all_filtered %>% filter(Comparison == x, Type == "SE") %>% pull(geneSymbol) %>% unique() %>% paste(collapse = ", ")),
  Count = sapply(Comparisons, function(x) AS_all_filtered %>% filter(Comparison == x, Type == "SE") %>% pull(geneSymbol) %>% unique() %>% length())
)








################################# PSI plot #####################################

files_list_psi <- list.files(path = "path/to/*.psi/files", pattern = "*.psi", full.names = TRUE, recursive = FALSE)

# make a table containing the sample information
sampleTable <- data.frame(
  sample_name = c("3598dRP_S7", "8579LP_S11", "8580dLP_S15", "9859LP_S16",
                  "3598dLP_S6", "8579RP_S12", "8580LP_S13", "8580RP_S14",
                  "8578RP_S9", "3597LP_S3", "3597RP_S4", "7482NP_S8",
                  "8578dLP_S10", "3596LP_S1", "3596dLP_S2", "3597dRP_S5",
                  "P5346dRP_S19", "T63591_S24", "P5496RP_S21", "P5496dRP_S23",
                  "P5346LP_S17", "P5346RP_S18", "P5496dLP_S22", "P5496NP_S20"),
  sample = c("Ctrl PBS 1", "Ctrl PBS 2", "Ctrl PBS 3", "Ctrl PBS 4", 
             "Ctrl ASO 1", "Ctrl ASO 2", "Ctrl ASO 3", "Ctrl ASO 4",
             "2Hit PBS 1", "2Hit PBS 2", "2Hit PBS 3", "2Hit PBS 4",
             "2Hit ASO 1", "2Hit ASO 2", "2Hit ASO 3", "2Hit ASO 4",
             "RBM20 Ctrl PBS 1", "RBM20 Ctrl PBS 2", "RBM20 Ctrl PBS 3", "RBM20 Ctrl PBS 4",
             "RBM20 Ctrl ASO 1", "RBM20 Ctrl ASO 2", "RBM20 Ctrl ASO 3", "RBM20 Ctrl ASO 4"),
  condition = c(rep("Ctrl PBS", 4),
                rep("Ctrl ASO", 4),
                rep("2Hit PBS", 4),
                rep("2Hit ASO", 4),
                rep("RBM20 Ctrl PBS", 4),
                rep("RBM20 Ctrl ASO", 4))
)

# Create a lookup table for files
sample_basenames <- sub(".*_(S[0-9]+)$", "\\1", sampleTable$sample_name)
file_basenames <- sub(".*(S[0-9]+)_PSI\\.psi$", "\\1", basename(files_list_psi))

file_lookup <- data.frame(
  file_path = files_list_psi,
  basename = file_basenames,
  stringsAsFactors = FALSE
)

# Match the basenames and add the file path to the sample table
sampleTable$file_path <- file_lookup$file_path[match(sample_basenames, file_lookup$basename)]

sampleTable %<>%
  filter(condition %in% c("Ctrl PBS", "Ctrl ASO", "2Hit PBS", "2Hit ASO"))

psi_dfs <- lapply(sampleTable$file_path, function(x) read_tsv(x, col_names = c("exonID", "PSI", "low_inclusion_filter")))
names(psi_dfs) <- sampleTable$sample
combined_list_psi <- list()
combined_list_psi <- imap(psi_dfs, ~ mutate(.x, sample = .y))
PSI_all <- bind_rows(combined_list_psi)

PSI_all %<>% 
  left_join(sampleTable, by = "sample")


PSI_all_adj <- PSI_all %>% 
  separate(exonID, into = c("geneID", "exon"), sep = ":") %>% 
  mutate(exon = exon %>% as.numeric) %>% 
  mutate(condition = factor(condition, levels = c("Ctrl PBS", "Ctrl ASO", "2Hit PBS", "2Hit ASO")))



gtf <- import("/path/to/Mus_musculus.GRCm39.110.canonical.gtf")


psi_plot <- function(psi_table, gene){
  
  ID <- gtf %>% 
    as.data.frame() %>% 
    filter(gene_name == gene) %>% 
    pull(gene_id) %>% 
    unique
  
  strand <- unique(strand(subset(gtf, gene_id == ID)))
  
  psi_df <- psi_table %>%
    filter(str_detect(geneID, ID)) %>% 
    group_by(exon, condition) %>%
    dplyr::summarize(mean_PSI = mean(PSI, na.rm = TRUE))
  
  max_exon <- max(psi_df$exon)
  min_exon <- min(psi_df$exon)
  
  if (max_exon > 100) {
    limits <- seq(min_exon, max_exon, 5)
    angle <- 90
    hjust <- 1
  } else {
    limits <- seq(min_exon, max_exon, 1)
    angle <- 0
    hjust <- 0.5
  }
  
  if(strand == "-"){
    psi_df$exon <- rev(psi_df$exon)
  } else {
    limits <- seq(min_exon, max_exon)
  }
  
  p <- ggplot(psi_df, aes(exon, mean_PSI, color = condition)) +
    geom_line(aes(y = mean_PSI, group = condition), size = .6) +
    scale_color_manual(values = c("grey60", "grey10", "red", "darkred")) +
    ggtitle(gene) +
    labs(y = "mean PSI",
         color = "") +
    scale_x_discrete(limits = limits) +
    scale_y_continuous(limits = c(0, 1, .2)) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", size = .2),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.ticks = element_line(color = "grey80", size = .2),
      axis.text.x = element_text(color = "black", angle = angle, vjust = 0.5, hjust = hjust),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.position = "bottom")
  
  return(p)
}


Figure_5B <- psi_plot(PSI_all_adj, "Ttn")
Figure_5C_1 <- psi_plot(PSI_all_adj, "Camk2d")
Figure_5C_2 <- psi_plot(PSI_all_adj, "Ldb3")
Figure_S3A <- psi_plot(PSI_all_adj, "Ank3")
Figure_S3B <- psi_plot(PSI_all_adj, "Ryr2")





# ------------------------------- #
# matrix-associated gene splicing #
# ------------------------------- #

matrix_genes <- read_xlsx("path/to/matrix/genes.xlsx", sheet = "consolidated") %>% # matrix-related genes derived from GO:0030198
  dplyr::rename(Gene = Symbol)

unique(matrix_genes$Gene)

AS_all_matrix <- AS_all %>% 
  filter(geneSymbol %in% matrix_genes$Gene) %>% 
  filter(Comparison %in% Comparisons) %>% 
  filter(Type == "SE")

VennDiag <- list(`matrix genes` = unique(matrix_genes$Gene), 
                 `alternatively spliced` = AS_all_filtered %>% filter(Type == "SE" & Comparison %in% Comparisons) %>% pull(geneSymbol) %>% unique)

VennDiag <- lapply(VennDiag, function(list){
  duplicates <- duplicated(list)
  list <- list[!duplicates]
})

# Create the eulerr object
eul <- euler(VennDiag, shape = "ellipse", proportional = TRUE)

# Generate the Venn diagram
Figure_S4C <- plot(eul, fill = c("#D6604D", "grey70"),
                  quantities = list(type = "counts"),
                  cex = 1)


venn_parts_genes <- list(
  "matrix genes" = setdiff(VennDiag$`matrix genes`, VennDiag$`alternatively spliced`),
  "alternatively spliced" = setdiff(VennDiag$`alternatively spliced`, VennDiag$`matrix genes`),
  "spliced matrix genes" = intersect(VennDiag$`matrix genes`, VennDiag$`alternatively spliced`)
)


matrix <- data.frame(
  Gene = venn_parts_genes$`matrix genes`,
  part = "matrix_genes"
)

alt_spliced <- data.frame(
  Gene = venn_parts_genes$`alternatively spliced`,
  part = "alternatively spliced"
)

overlap <- data.frame(
  Gene = venn_parts_genes$`spliced matrix genes`,
  part = "spliced matrix genes"
)


Datafile_S6 <- bind_rows(matrix, alt_spliced, overlap)



################### retrieve ensembl IDs for each gene name ####################

# Specify the Ensembl dataset and attributes you want to retrieve
ensembl_dataset <- "mmusculus_gene_ensembl" 
ensembl_attributes <- c("ensembl_gene_id", "external_gene_name") 

# Create a biomart object to access the Ensembl database
ensembl <- useMart("ensembl", dataset = ensembl_dataset)#, host = "https://useast.ensembl.org/")

ensembl_gene_names <- Datafile_S6$Gene

# Retrieve gene names using the Ensembl gene IDs
gene_ID <- getBM(attributes = ensembl_attributes, 
                 filters = "external_gene_name", 
                 values = ensembl_gene_names, 
                 mart = ensembl) %>% 
  dplyr::rename(geneID = ensembl_gene_id,
                Gene = external_gene_name)


Datafile_S6 %<>%
  left_join(gene_ID, by = "Gene") %>% 
  relocate(part)


