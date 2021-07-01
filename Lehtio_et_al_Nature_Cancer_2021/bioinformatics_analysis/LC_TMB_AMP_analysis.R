## LC TMB - APM relationship / Immune gene sets

# Libraries
library(openxlsx) 
library(GSA)
library(GSVA)
library(reshape2)
library(dplyr)


# Data         
load('files/LC_colours')

proteomics_init <- read.table('files/LC_Landscape_Proteome.txt')

# Complete cases
proteomics <- proteomics_init[complete.cases(proteomics_init),]

LC_metadata_init <- read.xlsx('files/Supplementary_table_1.xlsx', rowNames = TRUE)


# Remove sample with NA TMB
idx <- which(is.na(LC_metadata_init$TMB))
LC_metadata <- LC_metadata_init[-idx, ]
proteomics <- proteomics[, -idx]


# Log2 transform
LC_metadata$TMB <- log2(LC_metadata$TMB)


msig <- GSA.read.gmt('files/c2.cp.kegg.v6.2.symbols.gmt.txt')
enrich_list <- msig$genesets
names(enrich_list) <- msig$geneset.names
antigen_list <- enrich_list['KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION']


# Interferon hallmark
msig <- GSA.read.gmt('files/h.all.v6.2.symbols.gmt.txt')
enrich_list <- msig$genesets
names(enrich_list) <- msig$geneset.names
interferon_list <-  enrich_list[c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")]


# Immune-list taken from Charoentong, P. et al. Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade. Cell Rep 18, 248-262, doi:10.1016/j.celrep.2016.12.019 (2017).

immunedt <- read.xlsx('files/Immunophenogram_1-s2.0-S2211124716317090-mmc3.xlsx',startRow = 2)

immune_list <- lapply(unique(immunedt$Cell.type), function(i) {
  immunedt$Metagene[immunedt$Cell.type == i]
})

names(immune_list) <- unique(immunedt$Cell.type)

all_lists <- c(antigen_list, interferon_list, immune_list)

# overlap_genes <- lapply(all_lists, intersect, rownames(proteomics))

names(all_lists)[match(c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "Effector memeory CD4 T cell", "Effector memeory CD8 T cell", "Gamma delta T cell"), names(all_lists))] <- c("IFN-alpha", "IFN-gamma", "Effector memory CD4 T cell", "Effector memory CD8 T cell", "Gamma-delta T cell")

# GSEA
scores <- gsva(as.matrix(proteomics), gset.idx.list = all_lists, method = "ssgsea")

# Immune pathway enrichment
toplot2 <- as.data.frame(t(scores[-1,]))
toplot2$Subtype <- LC_metadata$Proteome.Subtype
per_subtype_mean <- as.data.frame(toplot2 %>% group_by(Subtype) %>% summarise_all(mean))
per_subtype_mean <- melt(per_subtype_mean, id.vars = 'Subtype')

per_subtype_mean$variable <- factor(per_subtype_mean$variable, 
                                    levels =  rev(c('Activated CD8 T cell', 'Central memory CD8 T cell', 'Effector memory CD8 T cell', 'Activated CD4 T cell', 'Central memory CD4 T cell', 'Effector memory CD4 T cell', 'T follicular helper cell', 'Gamma-delta T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Type 2 T helper cell', 'Regulatory T cell', 'Activated B cell', 'Immature  B cell', 'Memory B cell','Natural killer cell', 'CD56bright natural killer cell', 'CD56dim natural killer cell', 'MDSC', 'Natural killer T cell', 'Activated dendritic cell', 'Plasmacytoid dendritic cell', 'Immature dendritic cell', 'Macrophage', 'Eosinophil', 'Mast cell', 'Monocyte', 'Neutrophil', "IFN-alpha", "IFN-gamma")))

per_subtype_mean$group <- ifelse(per_subtype_mean$variable %in% c('Activated CD8 T cell', 'Central memory CD8 T cell', 'Effector memory CD8 T cell', 'Activated CD4 T cell', 'Central memory CD4 T cell', 'Effector memory CD4 T cell', 'T follicular helper cell', 'Gamma-delta T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Type 2 T helper cell', 'Regulatory T cell'), 'T-cell subsets', ifelse(per_subtype_mean$variable %in% c('Activated B cell', 'Immature  B cell', 'Memory B cell'), 'B-cell subsets', ifelse(per_subtype_mean$variable %in% c('Natural killer cell', 'CD56bright natural killer cell', 'CD56dim natural killer cell'), 'NK-cell subsets', ifelse(per_subtype_mean$variable %in% c('MDSC', 'Natural killer T cell', 'Activated dendritic cell', 'Plasmacytoid dendritic cell', 'Immature dendritic cell', 'Macrophage', 'Eosinophil', 'Mast cell', 'Monocyte', 'Neutrophil'), 'Other', 'Hallmark gene sets'))))

per_subtype_mean$group <- factor(per_subtype_mean$group, levels = c('T-cell subsets', 'B-cell subsets', 'NK-cell subsets', 'Other', 'Hallmark gene sets'))



## K-means TMB
tmb_kmean <- kmeans(LC_metadata$TMB,centers = c(mean(sort(LC_metadata$TMB)[1:5]), mean(sort(LC_metadata$TMB, decreasing = TRUE)[1:5])))
tmb_cluster <- tmb_kmean$cluster


# K-means APM
apm <- scores['KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION', ]
apm_kmean <- kmeans(apm,centers = c(mean(sort(apm)[1:5]), mean(sort(apm, decreasing = TRUE)[1:5])))
apm_cluster <- apm_kmean$cluster


# Plot
toplot <- data.frame('LU' = rownames(LC_metadata), 
                     'apm' = apm, 
                     'tmb' = LC_metadata$TMB, 
                     'apm_group' = apm_cluster, 
                     'tmb_group' = tmb_cluster, 
                     'Subtype' = as.factor(LC_metadata$Proteome.Subtype))

toplot$Cluster <- ifelse(toplot$apm_group == 1 & toplot$tmb_group == 1, 'TMB-L/APM-L', ifelse(toplot$apm_group == 1 & toplot$tmb_group == 2, 'TMB-H/APM-L', ifelse(toplot$apm_group == 2 & toplot$tmb_group == 1, 'TMB-L/APM-H', 'TMB-H/APM-H')))

# Finding minimum
apm_line <- as.data.frame(toplot %>% group_by(apm_group) %>% summarize(min(apm)))
tmb_line <- as.data.frame(toplot %>% group_by(tmb_group) %>% summarize(min(tmb)))



# Enrichment test for each quadrant
fisherdt <- do.call(rbind, lapply(unique(toplot$Cluster), function(i) {
  fisher_res <- sapply(unique(toplot$Subtype), function(j) {
    n1 <- nrow(subset(toplot, Cluster == i & Subtype == j))
    n2 <- nrow(subset(toplot, Cluster == i & Subtype != j))
    n3 <- nrow(subset(toplot, Cluster != i & Subtype == j))
    n4 <- nrow(subset(toplot, Cluster != i & Subtype != j))
    tab <- matrix(c(n1, n2, n3, n4), nrow = 2, byrow = TRUE)
    res <- c(fisher.test(tab)$p.value, log2(( n1/(n1+n2) ) /( (n1+n3)/(n1+n2+n3+n4) )))
  })
  fisher_res[1, ] <- p.adjust(fisher_res[1, ], method = 'BH')
  fisher_res <- as.data.frame(t(fisher_res))
  colnames(fisher_res) <- c('P-value', 'log2FE')
  tb <- cbind.data.frame(fisher_res, 'Cluster' = i, 'Subtype' = 1:6)
}))

###