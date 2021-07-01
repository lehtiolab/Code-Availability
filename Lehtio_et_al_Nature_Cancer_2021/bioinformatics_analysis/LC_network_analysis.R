### LC Network analysis

## Libraries 
library(openxlsx)
library(Seurat)
library(PCAtools)
library(ggplot2)
library(GSA)
library(reshape2)
library(clusterProfiler)
library(dplyr)


## Data
load('files/LC_colours')
proteomics <- read.table('files/LC_Landscape_Proteome.txt')
LC_metadata <- read.xlsx('files/Supplementary_table_1.xlsx', rowNames = TRUE)


## Proteins differentially expressed 
resDeqMS_fc_lk <- read.table('files/DEqMS_LC_log2FC_all.txt', header = TRUE, row.names = 1)
resDeqMS_pval_lk <- read.table('files/DEqMS_LC_Pvalues_all.txt', header = TRUE, row.names = 1)

fc_value <- 0.5
p_value <-  0.01

fc_thres <- which(abs(resDeqMS_fc_lk) > fc_value, arr.ind = TRUE)
p_thres <- which(resDeqMS_pval_lk < p_value, arr.ind = TRUE)

distinct_proteins <- intersect(apply(fc_thres, 1, paste, collapse = '_'),
                               apply(p_thres, 1, paste, collapse = '_'))
distinct_proteins  <- do.call(rbind, strsplit(distinct_proteins, '_'))
distinct_proteins <- unique(distinct_proteins[, 1])
distinct_proteins <- rownames(resDeqMS_fc_lk)[as.numeric(distinct_proteins)]

feature_count <- apply(proteomics, 1, function(x){ sum(!is.na(x)) })

sample_perc <- 0.7
feature_count <- names(feature_count)[feature_count >= sample_perc * ncol(proteomics)]

distinct_proteins <- intersect(feature_count, distinct_proteins)

coreg <- cor(t(proteomics[distinct_proteins, ]), use = 'pairwise.complete.obs', method = 'pearson')


## Select Principal components
coreg_seurat <- CreateSeuratObject(as.matrix(coreg))
all_proteins <-  colnames(coreg)

## Scale correlation values
coreg_seurat <- ScaleData(coreg_seurat, features = all_proteins)

## Run pca
coreg_seurat <- RunPCA(coreg_seurat, features = colnames(coreg))

# ElbowPlot(coreg_seurat, ndims = 50)
pca_std <- Stdev(object = coreg_seurat, reduction = 'pca')
chosen <- findElbowPoint(pca_std)
##########################################################################################################
# # Dimensionality reduction
# # UMAP
# lapply(c(seq(5, 30, 5)), function(k) {
#   lapply(c(0.001, seq(0.005, 0.5, 0.05), 0.5), function(l) {
#       coreg_seurat <- RunUMAP(coreg_seurat, dims = 1:8, check_duplicates = FALSE, n.neighbors = k, min.dist = l)
# 
#   pdf(paste0('files/LC_UMAP_coreg_',k, '_', l, '.pdf'), width = 6, height = 6)
#   p <- DimPlot(coreg_seurat, reduction = "umap", cols = 'black')
#   
#   p <- p + labs(title = paste0('k = ', k, ', ', 'l = ', l))  + guides(col = FALSE) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
# 
#   plot(p)
#   dev.off()
# 
#   })
# })

## Final choice
k = 20
l = 0.2
coreg_seurat <- RunUMAP(coreg_seurat, dims = 1:8, n.neighbors = k, min.dist = l)
###################################################################################################
# Clustering assessed by significance for enrichment of MSigDB Hallmark gene sets (Liberzon, A. et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst 1, 417-425, doi:10.1016/j.cels.2015.12.004 (2015))
input_list <- GSA.read.gmt('files/h.all.v6.2.symbols.gmt.txt')
enrich_list <- input_list$genesets
names(enrich_list) <- input_list$geneset.names
 enrich_list <- melt(enrich_list)
colnames(enrich_list) = c('gene', 'name')
enrich_list <- enrich_list[, c('name', 'gene')]

coreg_seurat <- FindNeighbors(coreg_seurat, dims = 1:8, features = NULL,
                              k.param = 20,
                              compute.SNN = TRUE,
                              prune.SNN = 1/15,
                              nn.method = "rann",
                              annoy.metric = "euclidean",
                              nn.eps = 0,
                              verbose = TRUE,
                              force.recalc = FALSE)


# resolution_seq <- seq(0.2, 4, 0.2)
# lapply(resolution_seq, function(i) {
# 
#   coreg_seurat <- FindClusters(coreg_seurat, resolution = i)
#   tosubset <- paste0('RNA_snn_res.', i)
# 
#   modules <- coreg_seurat[[tosubset]]
#   idx_sorted <- sort(unique(modules[, 1]))
# 
#   module_enrichment <- do.call(rbind, lapply(idx_sorted, function(i) {
#     # print(i)
#     res <- enricher(gene = rownames(modules)[modules[,1] %in%  i],
#                     pvalueCutoff = 0.05,
#                     pAdjustMethod = "BH",
#                     universe = rownames(modules),
#                     minGSSize = 10,
#                     maxGSSize = 500,
#                     qvalueCutoff = 0.05,
#                     TERM2GENE = enrich_list,
#                     TERM2NAME = NA )
#     res <- as.data.frame(res)
# 
#     if(nrow(res) == 0) {
#       NA} else {res$Module <- i}
# 
#     res
# 
#   }))
# 
# 
# 
#   p_res <- as.data.frame(table(modules))
#   p_res$Enriched <- ifelse(p_res$modules %in% module_enrichment$Module, '*', '')
# 
# 
#   p <- ggplot(p_res ,aes(x = modules, y = Freq)) +
#     geom_bar(stat = 'identity') +
#     geom_text(aes(label = Enriched), vjust = -0.2, size = 8) +
#     scale_y_continuous(expand = c(0,0, 0.1, 0), ) +
#     scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
#     labs(x = 'Clusters', y ='# Proteins', title = paste0('Resolution = ', i)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5))
# 
#   pdf(paste0('files/LC_Louvain_',i,'.pdf'), width = 7, height = 4)
#   plot(p)
#   dev.off()
# 
#   })

## Plot network
n_resolution <- 0.6
coreg_seurat <- FindClusters(coreg_seurat, resolution = n_resolution)


modules <- coreg_seurat[[paste0('RNA_snn_res.', n_resolution)]]
freq_modules <- table(modules[,1])
idx_sorted <- sort(unique(modules[,1]))


module_enrichment <- do.call(rbind, lapply(idx_sorted, function(i) {
  # print(i)
  res <- enricher(gene = rownames(modules)[modules[,1] %in%  i], 
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", 
                  universe = rownames(modules),
                  minGSSize = 10, 
                  maxGSSize = 500, 
                  qvalueCutoff = 0.05, 
                  TERM2GENE = enrich_list,
                  TERM2NAME = NA)
  res <- as.data.frame(res)
  if(nrow(res) == 0) {res[1, ] = NA} else { 
    res$Module <- i
    return(res)}
}))


module_enrichment <- module_enrichment[order(module_enrichment$Count, decreasing = TRUE), ]
##################################################################################################
## Cell type annotation (Travaglini, K. J. et al. A molecular cell atlas of the human lung from single-cell RNA sequencing. Nature 587, 619-625, doi:10.1038/s41586-020-2922-4 (2020))
cells <- read.table('files/LC_project_cell_types.txt', header = TRUE)
cells <- cells[, c("cell_type", "cell_genes")]
colnames(cells) <- c('name', 'gene')

cell_enrichment <- do.call(rbind, lapply(idx_sorted, function(i) {
  # print(i)
  res <- enricher(gene = rownames(modules)[modules[,1] %in%  i], 
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH", 
                  universe = rownames(modules),
                  minGSSize = 10, 
                  maxGSSize = 500, 
                  qvalueCutoff = 0.05, 
                  TERM2GENE = cells,
                  TERM2NAME = NA)
  res <- as.data.frame(res)
  if(nrow(res) == 0) {res[1, ] = NA} else { 
    res$Module <- unique(module_enrichment$Module[match(i,module_enrichment$old_module)])
    return(res)}
}))

cell_enrichment <- cell_enrichment[order(cell_enrichment$Count, decreasing = TRUE), ]
############################################################################################################
## Group-wise median abundances
module_enrichment_unique <-  module_enrichment[!duplicated(module_enrichment$Module), ]

dat <- coreg_seurat
plotUMAP <- as.data.frame(Embeddings(object = dat[['umap']]))
plotUMAP$Module <- module_enrichment_unique$Module[match(dat$seurat_clusters, module_enrichment_unique$Module)]


plotUMAP[, paste0('Subtype', sort(unique(LC_metadata$Proteome.Subtype)))] <- sapply(sort(unique(LC_metadata$Proteome.Subtype)), function(i) {
  prot_sub <- proteomics[rownames(plotUMAP), LC_metadata$Proteome.Subtype == i]
  apply(prot_sub, 1, median, na.rm = TRUE)
})

plotUMAP$proteins <- rownames(plotUMAP)
plotUMAP <- reshape2::melt(plotUMAP, id.vars =  c('proteins', "Module", "UMAP_1", "UMAP_2"))

# Average ratios
avg_ratios <- do.call(cbind, lapply(unique(plotUMAP$Module), function(i) {
  res <- subset(plotUMAP, plotUMAP$Module == i)
  fc_res <- res %>% group_by(variable) %>% summarize('avg' = mean(value))
  fc_res$avg
}))

###