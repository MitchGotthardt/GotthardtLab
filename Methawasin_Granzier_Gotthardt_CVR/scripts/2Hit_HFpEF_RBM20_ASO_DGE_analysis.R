################################################################################

# Manuscript: Rbm20 antisense oligonucleotides alleviate diastolic dysfunction in a mouse model of cardiometabolic HFpEF.
# Authors: Mei Methawasin*, Stefan Meinke, Michael H Radke, Gerrie P Farman, Zaynab Hourani, John E. Smith III, Wei Guo, Henk Granzier*, Michael Gotthardt*
# *Corresponding author
# Bioinformatics and data analysis: Stefan Meinke

################################################################################

################################################################################
############### Differential Gene Expression Analysis ##########################
################################################################################

required_libraries <- c("BiocManager",
                        "readr", 
                        "tidyverse", 
                        "magrittr", 
                        "ggplot2", 
                        "ggforce",
                        "VennDiagram",
                        "ggvenn",
                        "ggVennDiagram",
                        "dplyr", 
                        "readxl",
                        "writexl",
                        "openxlsx",
                        "ggbeeswarm",
                        "stringr",
                        "scales",
                        "pheatmap",
                        "purrr",
                        "eulerr",
                        "UpSetR",
                        "clusterProfiler",
                        "DESeq2",
                        "edgeR",
                        "rtracklayer",
                        "biomaRt",
                        "org.Mm.eg.db",
                        "GenomicFeatures")

# Check and install/load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


for (package in required_libraries) {
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
  # Load the package
  library(package, character.only = TRUE)
}


################# load and prepare the read counts dataframe ###################

dir <- "/path/to/*ReadsPerGene.out.tab/files" 

# Load read count files
files_list_GE <- list.files(path = dir, 
                            pattern = "*ReadsPerGene.out.tab$", 
                            full.names = TRUE, 
                            recursive = TRUE)

dfs_GE <- lapply(files_list_GE, function(x) read_tsv(x, 
                                                     skip=4, 
                                                     col_names = c("gene", "counts_unstranded", "counts_firstStrand", "counts_secondStrand")))

# name each dataframe based on sample name
dfs_GE <- setNames(dfs_GE, gsub("\\_ReadsPerGene.out.tab$", "", basename(files_list_GE)))

# join all GE dataframes, select only the "gene" and "counts_secondStrand" column of all GE dataframes
dfs_GE <- lapply(dfs_GE, function(x) {
  first_col_name <- colnames(x)[1]
  last_col_name <- colnames(x)[ncol(x)]
  
  # Select the first and last columns using the column names
  selected_cols <- dplyr::select(x, 
                                 geneID = !!first_col_name, 
                                 Sample = !!last_col_name)
})

df_names <- names(dfs_GE)

# rename all "counts_secondStrand" columns according to the dataframe name
dfs_GE <- map2(dfs_GE, df_names, ~setNames(.x, c("geneID", .y)))

# Perform left joins using reduce()
GE_raw <- purrr::reduce(dfs_GE, left_join, by = "geneID")





################### retrieve gene names for each ensembl ID ####################

# Specify the Ensembl dataset and attributes you want to retrieve
ensembl_dataset <- "mmusculus_gene_ensembl" 
ensembl_attributes <- c("ensembl_gene_id", "external_gene_name", "description") 

# Create a biomart object to access the Ensembl database
ensembl <- useMart("ensembl", dataset = ensembl_dataset)#, host = "https://useast.ensembl.org/")

ensembl_gene_ids <- GE_raw$geneID

# Retrieve gene names using the Ensembl gene IDs
gene_ID <- getBM(attributes = ensembl_attributes, 
                 filters = "ensembl_gene_id", 
                 values = ensembl_gene_ids, 
                 mart = ensembl) %>% 
  dplyr::rename(geneID = ensembl_gene_id) %>% 
  mutate(description = str_extract(description, ".*(?=\\[)")) %>% 
  dplyr::rename(Gene = external_gene_name)


GE_raw %<>%
  left_join(gene_ID, by = "geneID") %>% 
  relocate(Gene)







################### Differential Gene Expression using DESeq2 ##################

# prepare the sample/condition Table
sampleTable <- data.frame(
  sample_name = c("3598dRP_S7", "8579LP_S11", "8580dLP_S15", "9859LP_S16",
                  "3598dLP_S6", "8579RP_S12", "8580LP_S13", "8580RP_S14",
                  "8578RP_S9", "3597LP_S3", "3597RP_S4", "7482NP_S8",
                  "8578dLP_S10", "3596LP_S1", "3596dLP_S2", "3597dRP_S5"),
  sample = c("Ctrl PBS 1", "Ctrl PBS 2", "Ctrl PBS 3", "Ctrl PBS 4", 
             "Ctrl ASO 1", "Ctrl ASO 2", "Ctrl ASO 3", "Ctrl ASO 4",
             "2Hit PBS 1", "2Hit PBS 2", "2Hit PBS 3", "2Hit PBS 4",
             "2Hit ASO 1", "2Hit ASO 2", "2Hit ASO 3", "2Hit ASO 4"),
  condition = c(rep("Ctrl PBS", 4),
                rep("Ctrl ASO", 4),
                rep("2Hit PBS", 4),
                rep("2Hit ASO", 4))
)

sampleTable$condition <- factor(sampleTable$condition, levels=c("Ctrl PBS", "Ctrl ASO", "2Hit PBS", "2Hit ASO"))
sampleTable %<>% 
  column_to_rownames("sample_name")

GE_raw_ordered <- GE_raw[, match(rownames(sampleTable), colnames(GE_raw))]
GE_raw_ordered$geneID <- GE_raw$geneID

GE_raw_ordered %<>% 
  column_to_rownames("geneID") %>% 
  drop_na() %>%
  round

#Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = GE_raw_ordered,
                              colData = sampleTable,
                              design = ~ condition)

# Normalize data
dds <- DESeq(dds)

# filter for genes with read counts >=5 in at least 3 samples
dds <- dds[rowSums(counts(dds) >= 5) >= 3, ]

# Test for differential expression of all possible condition combinations, Specify the contrast of interest if there are more than two conditions.
combinations <- combn(dds$condition %>% unique, 2)
rev_combinations <- combinations[c(2, 1), ]

# Create a list to store the results
DESeq_results <- list()
# Iterate over combinations using map2 function
map2(rev_combinations[1, ], rev_combinations[2, ], ~ {
  condition1 <- .x %>% as.character
  condition2 <- .y %>% as.character
  
  # Create the contrast argument dynamically
  contrast <- c("condition", condition1, condition2)
  DESeq_res <- results(dds, contrast = contrast) %>% 
    as.data.frame %>% 
    rownames_to_column("geneID") %>% 
    mutate(Combination = paste(condition1, condition2, sep = " vs "))
  
  DESeq_res
}) %>%
  # Combine the results into a single data frame
  bind_rows(.id = NULL) -> DESeq_results

# View the resulting data frame
DESeq_results %<>% 
  drop_na %>% 
  left_join(gene_ID, by = "geneID") %>%
  relocate(Gene)


DESeq_results_sig_padj <- DESeq_results %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5)





######################## Principal Component analysis ###########################

dds_pca = estimateSizeFactors(dds_pca)
dds_pca = estimateDispersions(object = dds_pca, 
                              fitType = "parametric", 
                              quiet = TRUE)

vsd = varianceStabilizingTransformation(object = dds_pca, 
                                        blind = TRUE,           
                                        fitType = "parametric")

pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE, ntop = 500)
percentVar <- round(100*attr(pcaData, "percentVar"))

pcaData$condition <- factor(pcaData$condition, levels = c("Ctrl PBS", "Ctrl ASO", "2Hit PBS", "2Hit ASO"))

#################
### Figure 6A ###
#################

ggplot(pcaData, aes(PC1, PC2)) +
  geom_point(size = 3, shape = 21, aes(fill = condition)) +
  scale_fill_manual(values = c("grey90", "grey50", "red", "darkred")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  labs(fill = "") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "right")




############################ Calculate TpM values ##############################

# load gtf file
gtf_file <- "//mdc-berlin.net/fs/fast/AG_Gotthardt/smeinke/Annotation_files/Mus_Musculus.GRCm39.110.gtf"

# create the txdb object
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")

# collect the exons per gene id
exons_per_gene <- exonsBy(txdb,by="gene")

# reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic_gene_sizes <- sum(width(reduce(exons_per_gene)))
exonic_gene_sizes <- as.data.frame(exonic_gene_sizes)

# create a data frame for gene lengths.
counts <- GE_raw_ordered

length_df <- exonic_gene_sizes %>% 
  dplyr::rename(length = exonic_gene_sizes)

length_df <- merge(counts, length_df, by = 0, all = FALSE)
length_df %<>% column_to_rownames("Row.names")

#set up the data using the DGEList function, filtering by expression levels, and calculating normalization factors.
dge <- DGEList(counts=counts, genes=data.frame(Length=length_df), group = sampleTable$condition)
dge <- calcNormFactors(dge)

# use edgeR's built-in functions to calculate rpkm, FPKM (Fragments Per Kilobase Million) values and subsequently convert them to TPM.
fpkm <- rpkm(dge, gene.length = dge$genes$Length.length, normalized.lib.sizes = TRUE)
tpm <- fpkm / colSums(fpkm) * 1e6


tpm <- tpm %>% 
  as.data.frame() %>% 
  rownames_to_column("geneID") %>% 
  left_join(gene_ID, by = "geneID") %>% 
  relocate(Gene)

tpm_long <- tpm %>% 
  pivot_longer(cols=-c(Gene,geneID,description),names_to = "Sample", values_to = "TpM") %>% 
  left_join(sampleTable %>% rownames_to_column("Sample"), by = "Sample") %>% 
  dplyr::select(Gene, geneID, sample, TpM, condition)

tpm_long$condition <- factor(tpm_long$condition, levels=c("Ctrl PBS", "Ctrl ASO", "2Hit PBS", "2Hit ASO"))







########################## Fibrosis-associated Genes ###########################
plot_gene <- function(gene_name) {
  stdev <- tpm_long %>%
    filter(Gene == gene_name) %>% 
    group_by(condition) %>% 
    dplyr::summarise(mean = mean(TpM),
                     sd = sd(TpM))
  
  gene_df <- tpm_long %>% 
    filter(Gene == gene_name) %>% 
    merge(., stdev, by = "condition")
  
  ggplot(gene_df %>% dplyr::select(-c(sample, TpM)) %>% unique, aes(condition, mean)) +
    geom_bar(stat="identity", aes(fill=condition), color = "black", alpha=0.7, position = "identity") +
    geom_errorbar(aes(x=condition, ymin=mean-sd, ymax=mean+sd), width = 0.3, alpha=0.9) +
    geom_point(stat = "identity", shape = 21, color = "black", data = gene_df, aes(y = TpM, fill = condition)) +
    scale_fill_manual(values = c("grey90", "grey50", "red", "darkred")) +
    scale_y_continuous(limits = c(0, max((gene_df$mean)+(0.5*gene_df$mean)))) + #expand = c(0,0), 
    labs(x="", y="TpM") +
    ggtitle(gene_name) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"))
}


# Function to generate and arrange multiple plots
plot_all_genes <- function(gene_names) {
  plots <- lapply(gene_names, function(gene) plot_gene(gene))
  Reduce(`+`, plots)  # Arrange plots in a grid
}

# List of genes
genes <- c("Timp1", "Mmp9", "Col3a1", "Col1a1", "Tgfb1", "Loxl2", "Rbm20", "Camk2d", "Ryr2", "Ldb3", "Ank3")

# Generate and arrange plots
all_plots <- plot_all_genes(genes)

# Print or save the combined plot
all_plots






################# z-score normalized fibrosis-associated gene ##################

#################################
### Figure 6F, S3C-D, and S4B ###
#################################

mean_sd <- tpm_long %>%
  filter(Gene %in% genes) %>% 
  group_by(Gene, condition) %>% 
  dplyr::summarise(mean = mean(TpM),
                   sd = sd(TpM))

# Extract the reference mean and sd for "Ctrl PBS"
ctrl_pbs <- mean_sd %>%
  filter(condition == "Ctrl PBS") %>%
  dplyr::select(Gene, mean, sd) %>%
  dplyr::rename(mean_ctrl_pbs = mean, sd_ctrl_pbs = sd)


# Join the reference values back to the original dataframe
mean_sd <- mean_sd %>%
  left_join(ctrl_pbs, by = "Gene")


tpm_mean_df <- tpm_long %>% 
  filter(Gene %in% genes) %>%
  inner_join(mean_sd, by = c("Gene", "condition"))

tpm_mean_df <- tpm_mean_df %>% 
  mutate(normalized_value = TpM / mean_ctrl_pbs,
         normalized_sd = sd / mean_ctrl_pbs)


mean_norm <- tpm_mean_df %>% 
  group_by(Gene, condition) %>% 
  dplyr::summarize(mean_norm = mean(normalized_value))

tpm_mean_df <- tpm_mean_df %>% 
  left_join(mean_norm, by = c("Gene","condition"))


plot_gene <- function(gene_name) {
  
  gene_df <- tpm_mean_df %>% 
    filter(Gene == gene_name)
  
  ggplot(gene_df, aes(condition, mean_norm)) +
    geom_bar(stat="identity", aes(fill=condition), color = "black", alpha=0.7, position = "identity") +
    geom_errorbar(aes(x=condition, ymin=mean_norm-normalized_sd, ymax=mean_norm+normalized_sd), width = 0.3, alpha=0.9) +
    geom_point(stat = "identity", shape = 21, color = "black", data = gene_df, aes(y = normalized_value, fill = condition)) +
    # geom_beeswarm(shape = 21, color = "black", data = gene_df, aes(y = normalized_value, fill = condition), cex = 3)+ 
    scale_fill_manual(values = c("grey90", "grey50", "red", "darkred")) +
    scale_y_continuous(limits = c(0, max((gene_df$mean_norm)+(0.5*gene_df$mean_norm)))) + 
    labs(x="", y="normalized expression") +
    ggtitle(gene_name) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 10),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 10),
          axis.text = element_text(color = "black", size = 10),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"))
}

# Generate and arrange plots
all_plots <- plot_all_genes(genes)



### statistics
genes <- c("Ank3", "Camk2d", "Ryr2", "Ldb3", "Col3a1", "Col1a1", "Tgfb1", "Timp1", "Mmp9", "Loxl2")

signif_results <- list()
for(gene in genes){
  df <- tpm_long %>% 
    filter(Gene == gene) %>% 
    filter(condition %in% c("2Hit ASO", "2Hit PBS", "Ctrl ASO", "Ctrl PBS")) %>% 
    separate(col = "condition", into = c("condition1", "condition2"), sep = " ")
  
  # Two-way ANOVA
  two_way_anova_result <- aov(TpM ~ condition1 * condition2, data = df)
  
  tukey_result <- TukeyHSD(two_way_anova_result)
  
  signif_results[[gene]] <- tukey_result
}








################################# Upset plot ###################################

#################
### Figure 6B ###
#################
Comparisons <- c("2Hit PBS vs Ctrl PBS",
                 "Ctrl ASO vs Ctrl PBS",
                 "2Hit ASO vs 2Hit PBS",
                 "2Hit ASO vs Ctrl PBS")

# create the upset_list
upset_list <- map(Comparisons, ~ DESeq_results_sig_padj %>%
                    filter(Combination == .x) %>%
                    pull(geneID) %>%
                    unique)

# Name the elements of the list with the comparison names
names(upset_list) <- Comparisons

Figure_6B <- upset(fromList(upset_list),
      order.by = "freq",
      sets.x.label = "number of genes",
      text.scale = 1.5,
      point.size = 2,
      mb.ratio=c(0.7,0.3))




################
### Table S7 ###
################
# Create a named vector for each set
all_genes <- unique(unlist(upset_list))
all_genes <- all_genes[!is.na(all_genes)]
upset_data <- as.data.frame(
  sapply(upset_list, function(x) all_genes %in% x)
)
row.names(upset_data) <- all_genes

# Extract intersections
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

# Create vector for mapping geneIDs with gene names
gene_id_to_name <- setNames(gene_ID$Gene, gene_ID$geneID)

# Function to map gene IDs to gene names
map_gene_names <- function(gene_ids) {
  gene_names <- gene_id_to_name[gene_ids]
  return(gene_names)
}

# Apply the function to the intersect_genes list
intersect_genes_with_names <- lapply(intersect_genes, map_gene_names)


Table_S8 <- data.frame(
  Intersection = names(intersect_genes_with_names), 
  Genes = sapply(intersect_genes_with_names, function(x) paste(unname(x), collapse = ", ")), 
  Count = sapply(intersect_genes_with_names, function(x) length(unname(x))) 
)







Comparisons <- unique(DESeq_results_sig_padj %>% 
                        filter(Combination %in% c("2Hit PBS vs Ctrl PBS", "Ctrl ASO vs Ctrl PBS", "2Hit ASO vs 2Hit PBS", "2Hit ASO vs Ctrl PBS")) %>% 
                        pull(Combination))

Table_S7 <- data.frame(
  `Grouped Comparison` = Comparisons,
  Genes = sapply(Comparisons, function(x) DESeq_results_sig_padj %>% filter(Combination == x) %>% pull(Gene) %>% unique() %>% paste(collapse = ", ")),
  Count = sapply(Comparisons, function(x) DESeq_results_sig_padj %>% filter(Combination == x) %>% pull(Gene) %>% unique() %>% length())
)


######################## Correlation Matrix ####################################

#################
### Figure 6C ###
#################

HitPBS_CtrlASO <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "2Hit PBS vs Ctrl ASO") %>% pull(geneID)))
HitASO_HitPBS <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "2Hit ASO vs 2Hit PBS") %>% pull(geneID)))
CtrlASO_CtrlPBS <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "Ctrl ASO vs Ctrl PBS") %>% pull(geneID)))  
HitASO_CtrlPBS <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "2Hit ASO vs Ctrl PBS") %>% pull(geneID)))
HitPBS_CtrlPBS <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "2Hit PBS vs Ctrl PBS") %>% pull(geneID)))
HitASO_CtrlASO <- length(unique(DESeq_results_sig_padj %>% filter(Combination == "2Hit ASO vs Ctrl ASO") %>% pull(geneID)))

corr_df <- data.frame(row.names = c("Ctrl ASO", "2Hit PBS", "2Hit ASO"),
                      `Ctrl PBS` = c(CtrlASO_CtrlPBS, HitPBS_CtrlPBS, HitASO_CtrlPBS),
                      `Ctrl ASO` = c(0, HitPBS_CtrlASO, HitASO_CtrlASO),
                      `2Hit PBS` = c(0,0,HitASO_HitPBS))

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

get_lower_tri(corr_df)

corr_melt <- get_lower_tri(corr_df) %>% 
  rownames_to_column("Var1") %>% 
  pivot_longer(cols = -Var1, names_to = "Var2", values_to = "n") %>% 
  mutate(Var2 = case_when(
    Var2 == "Ctrl.PBS" ~ "Ctrl PBS",
    Var2 == "X2Hit.PBS" ~ "2Hit PBS",
    Var2 == "Ctrl.ASO" ~ "Ctrl ASO"
  ))


corr_melt$Var1 <- factor(corr_melt$Var1, levels = c("Ctrl ASO", "2Hit PBS", "2Hit ASO"))
corr_melt$Var2 <- factor(corr_melt$Var2, levels = c("Ctrl PBS", "Ctrl ASO", "2Hit PBS"))


Figure_6C <- ggplot(data = corr_melt, aes(x=Var1, y=Var2, fill=n)) + 
  geom_tile(color = "black") +
  geom_text(aes(x = Var1, y = Var2, label = n), size = 4) +
  labs(x="", y="") +
  scale_fill_gradient(high="grey50",low = "white", na.value = "white", limits=c(200, 1000)) +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 90, vjust = 2, hjust = .5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color = "black"),
        legend.position = c(0,1),
        legend.justification = c(-0.4,1.15))




################
### Table S9 ###
################
HitASO_CtrlPBS <- DESeq_results_sig_padj %>% 
  filter(Combination == "2Hit ASO vs Ctrl PBS") %>% 
  dplyr::select(Gene, geneID)

HitPBS_CtrlPBS <- DESeq_results_sig_padj %>% 
  filter(Combination == "2Hit PBS vs Ctrl PBS") %>% 
  dplyr::select(Gene, geneID)

rev <- HitPBS_CtrlPBS %>% 
  filter(!(geneID %in% HitASO_CtrlPBS$geneID)) %>% 
  mutate(Regulation = "reverted")

mis <- HitASO_CtrlPBS %>% 
  filter(!(geneID %in% HitPBS_CtrlPBS$geneID)) %>% 
  mutate(Regulation = "misregulated")

co <- HitASO_CtrlPBS %>% 
  filter(geneID %in% HitPBS_CtrlPBS$geneID) %>% 
  mutate(Regulation = "coregulated")

nrow(rev)
nrow(mis)
nrow(co)
length(unique(rev$Gene))
length(unique(mis$Gene))


Table_S9 <- bind_rows(rev, mis, co)




############################# Venn Diagram #####################################

#################
### Figure 6E ###
#################

VennDiag <- list(`2Hit PBS vs Ctrl PBS` = DESeq_results_sig_padj %>% filter(Combination == "2Hit PBS vs Ctrl PBS") %>% pull(geneID),
                 `2Hit ASO vs Ctrl PBS` = DESeq_results_sig_padj %>% filter(Combination == "2Hit ASO vs Ctrl PBS") %>% pull(geneID))

VennDiag <- lapply(VennDiag, function(list){
  duplicates <- duplicated(list)
  list <- list[!duplicates]
})

# Create the eulerr object
eul <- euler(VennDiag, shape = "ellipse", proportional = TRUE)

# Generate the Venn diagram
Figure_6E <- plot(eul, fill = c("#D6604D", "grey70"),
                  quantities = list(type = "counts"),
                  cex = 1)




# Extract the sets for each part of the Venn diagram
venn_parts <- eul$original.values

# To view the specific genes for each part of the diagram:
venn_parts_genes <- list(
  "reverted" = setdiff(VennDiag$`2Hit PBS vs Ctrl PBS`, VennDiag$`2Hit ASO vs Ctrl PBS`),
  "misregulated" = setdiff(VennDiag$`2Hit ASO vs Ctrl PBS`, VennDiag$`2Hit PBS vs Ctrl PBS`),
  "coregulated" = intersect(VennDiag$`2Hit PBS vs Ctrl PBS`, VennDiag$`2Hit ASO vs Ctrl PBS`)
)


reverted <- data.frame(
  geneID = venn_parts_genes$reverted,
  Regulation = "reverted"
)

coregulated <- data.frame(
  geneID = venn_parts_genes$coregulated,
  Regulation = "coregulated"
)

misregulated <- data.frame(
  geneID = venn_parts_genes$misregulated,
  Regulation = "misregulated"
)


Table_S9 <- bind_rows(reverted, misregulated, coregulated) %>% 
  left_join(gene_ID %>% dplyr::select(-description), by = "geneID") %>% 
  relocate(Gene)







####################### Enrichment Analysis ####################################

#################
### Figure 6E ###
#################

### KEGG enrichment ###
KEGG_plot <- function(genes,background){
  
  kegg <- enrichKEGG(gene = genes,
                     organism = "mmu",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = background)
  
  
  kegg_df <- as.data.frame(kegg) %>%
    mutate(Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description)) %>% 
    arrange(p.adjust) %>% 
    mutate(Description = factor(Description, levels = rev(unique(Description))),
           padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
           `-log10(padj)` = -log10(p.adjust),
           enrichment = "KEGG")
  
  kegg_plot <- ggplot(head(kegg_df, 25), aes(x = Description, y =`-log10(padj)`)) +
    geom_bar(stat = "identity", fill = "grey90", color = "black") +
    geom_text(aes(label = Count), vjust = -1, hjust=1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, max(kegg_df$`-log10(padj)`+.5))) +
    labs(y = "-log10(padj)",
         x = "Enriched Term",
         size = "#Genes",
         fill = "") +
    ggtitle("KEGG enrichment") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = c(1, 0), legend.justification = c(1, 0))
  
  return(kegg_plot)
  
  return(kegg_df)
}


# for Wikipathway analysis
download.file("https://data.wikipathways.org/current/gmt/wikipathways-20240910-gmt-Mus_musculus.gmt", destfile = "Mus_musculus.gmt") # for original analysis used https://data.wikipathways.org/current/gmt/wikipathways-20240710-gmt-Mus_musculus.gmt

# Load the downloaded GMT file
gmt_data <- read.gmt("Mus_musculus.gmt")

wiki_plot <- function(df){
  ggplot(head(df, 25), aes(x = Description, y = `-log10(padj)`)) +
    geom_bar(stat = "identity", fill = "grey90", color = "black") +
    geom_text(aes(label = Count), vjust = -1, hjust=1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, max(df$`-log10(padj)`+.5))) +
    labs(y = "-log10(padj)",
         x = "Enriched Term",
         size = "#Genes",
         fill = "") +
    ggtitle("Wikipathways") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = c(1, 0), legend.justification = c(1, 0))
}
  

# Generate the enrichment plots
Regulation <- unique(Table_S9$Regulation)
background_genes <- unique(DESeq_results$geneID)

KEGG_plots <- list()
KEGG_dfs <- list()
Wiki_plots <- list()
Wiki_dfs <- list()

for (i in Regulation) {
  
  # Filter the geneIDs based on Regulation
  geneIDs <- Table_S9 %>% 
    filter(Regulation == i) %>% 
    pull(geneID)
  
  # Convert geneIDs to ENTREZ IDs for analysis
  entrez_ids <- bitr(geneIDs, 
                     fromType = "ENSEMBL", 
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
  
  # Convert background genes to ENTREZ IDs for background reference
  background_ids <- bitr(background_genes,
                         fromType = "ENSEMBL",
                         toType = "ENTREZID",
                         OrgDb = org.Mm.eg.db)
  
  ### Perform KEGG enrichment once
  kegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                     organism = "mmu",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     universe = background_ids$ENTREZID)
  
  # Prepare KEGG data frame
  kegg_df <- as.data.frame(kegg) %>%
    mutate(Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description)) %>%
    arrange(p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))),
           padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
           `-log10(padj)` = -log10(p.adjust),
           enrichment = "KEGG")
  
  # Save the KEGG data frame
  KEGG_dfs[[i]] <- kegg_df
  
  # Generate KEGG plot from the already prepared kegg_df
  kegg_plot <- ggplot(head(kegg_df, 25), aes(x = Description, y = `-log10(padj)`)) +
    geom_bar(stat = "identity", fill = "grey90", color = "black") +
    geom_text(aes(label = Count), vjust = -1, hjust = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, max(kegg_df$`-log10(padj)` + .5))) +
    labs(y = "-log10(padj)", x = "Enriched Term", size = "#Genes", fill = "") +
    ggtitle("KEGG enrichment") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = c(1, 0), legend.justification = c(1, 0))
  
  # Save the KEGG plot
  KEGG_plots[[i]] <- kegg_plot
  
  
  ### Perform WikiPathways enrichment
  wiki <- enricher(gene = entrez_ids$ENTREZID, TERM2GENE = gmt_data)
  
  # Prepare WikiPathways data frame
  wiki_df <- as.data.frame(wiki) %>%
    arrange(p.adjust) %>%
    mutate(padj = ifelse(`p.adjust` < 0.05, "padj < 0.05", "padj > 0.05"),
           `-log10(padj)` = -log10(p.adjust),
           enrichment = "Wikipathway")
  
  # Clean up the description column
  wiki_df$Description <- sub("%.*", "", wiki_df$Description)
  wiki_df$Description <- factor(wiki_df$Description, levels = unique(wiki_df$Description))
  
  # Save the WikiPathways data frame
  Wiki_dfs[[i]] <- wiki_df
  
  # Generate WikiPathways plot from the already prepared wiki_df
  wikipathway_plot <- ggplot(head(wiki_df, 25), aes(x = Description, y = `-log10(padj)`)) +
    geom_bar(stat = "identity", fill = "grey90", color = "black") +
    geom_text(aes(label = Count), vjust = -1, hjust = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_y_continuous(limits = c(0, max(wiki_df$`-log10(padj)` + .5))) +
    labs(y = "-log10(padj)", x = "Enriched Term", size = "#Genes", fill = "") +
    ggtitle("Wikipathways") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = c(1, 0), legend.justification = c(1, 0))
  
  # Save the WikiPathways plot
  Wiki_plots[[i]] <- wikipathway_plot
}






dfs_list <- list(
  Wiki_misregulated = Wiki_dfs$misregulated %>% dplyr::select(Description, p.adjust, Count, padj, `-log10(padj)`, enrichment),
  KEGG_coregulated = KEGG_dfs$coregulated %>% dplyr::select(Description, p.adjust, Count, padj, `-log10(padj)`, enrichment),
  Wiki_reverted = Wiki_dfs$reverted %>% dplyr::select(Description, p.adjust, Count, padj, `-log10(padj)`, enrichment)
)


enrichment_df <- bind_rows(dfs_list)
enrichment_df$Description <- factor(enrichment_df$Description, levels = unique(enrichment_df$Description))



Figure_6E <- ggplot(enrichment_df, aes(x = Description, y = `-log10(padj)`)) +
  geom_bar(stat = "identity", fill = "grey90", color = "black") +
  geom_text(aes(label = Count), vjust = -1, hjust=.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_y_continuous(limits = c(0, max(enrichment_df$`-log10(padj)`+.5))) +
  labs(y = "-log10(padj)",
       x = "") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = c(1, 0), legend.justification = c(1, 0))







Go_enrich <- function(genes,background,ont){
  go <- enrichGO(gene = genes,
                 universe = background,
                 OrgDb = org.Mm.eg.db,
                 ont = ont,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
  
  enriched_terms <- go$Description
  p_values <- go$pvalue
  padj <- go$p.adjust
  gene_counts <- go$Count
  
  result_df <- data.frame(
    Enriched_Term = enriched_terms,
    P_Value = p_values,
    Adjusted_P_Value = padj,
    Gene_Count = gene_counts)
  
  result_df %<>% 
    arrange(P_Value)
  
  top10 <- head(result_df, 10)
  
  # Reorder the Enriched_Term factor variable based on the P_Value column
  top10$Enriched_Term <- factor(top10$Enriched_Term, levels = top10$Enriched_Term)
  
  return(top10)
}


upregulated <- DESeq_results_sig_padj %>% filter(Combination == "Ctrl ASO vs Ctrl PBS" &
                                                   log2FoldChange > 0)

downregulated <- DESeq_results_sig_padj %>% filter(Combination == "Ctrl ASO vs Ctrl PBS" &
                                                     log2FoldChange < 0)

background_genes <- unique(DESeq_results$geneID)

geneIDs_upregulated <- upregulated %>% pull(geneID) %>% unique()
geneIDs_downregulated <- downregulated %>% pull(geneID) %>% unique()
  
entrez_ids_up <- bitr(geneIDs_upregulated, 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

entrez_ids_dn <- bitr(geneIDs_downregulated, 
                      fromType = "ENSEMBL", 
                      toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)

background_ids <- bitr(background_genes,
                       fromType = "ENSEMBL",
                       toType = "ENTREZID",
                       OrgDb = org.Mm.eg.db)
 
up_df <- Go_enrich(entrez_ids_up$ENTREZID, background_ids$ENTREZID, "BP")
dn_df <- Go_enrich(entrez_ids_dn$ENTREZID, background_ids$ENTREZID, "BP")

enrich_df <- bind_rows(
  up_df %>% mutate(Regulation = "upregulated"),
  dn_df %>% mutate(Regulation = "downregulated")
)

enrich_df %<>% 
  mutate(Regulation = factor(Regulation, levels = c("upregulated", "downregulated")))



Figure_S4A <- ggplot(enrich_df, aes(x = Adjusted_P_Value, y = reorder(Enriched_Term, -Adjusted_P_Value))) + 
  geom_bar(stat = "identity",color = "black", aes(fill = Regulation)) +
  geom_text(aes(label = Gene_Count), vjust = 0.5, hjust=-0.5) +
  scale_fill_manual(values = c("#D6604D","#4393C3")) +
  labs(x = "adjusted P-Value", y = "")+ 
  ggtitle("GO Biological Process (Ctrl-ASO/Ctrl-PBS)")+
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(trans = trans_reverser("log10")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 10),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = c(1, 0), legend.justification = c(1, 0))
