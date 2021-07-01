### LC consensus clustering

# Libraries
library(ConsensusClusterPlus)
library(reshape2)
library(openxlsx)
library(NMF)


# Data 
load('files/LC_colours')
proteomics <- read.table('files/LC_Landscape_Proteome.txt')
LC_metadata <- read.xlsx('files/Supplementary_table_1.xlsx', rowNames = TRUE)


# Remove NAs
proteomics <- proteomics[complete.cases(proteomics), ]

# Consensus clustering
results <- ConsensusClusterPlus(d = as.matrix(proteomics),
                                maxK = 11,  
                                pItem=0.8, 
                                pFeature = 1, 
                                reps = 1000, 
                                clusterAlg = 'hc',
                                innerLinkage = 'ward.D2', 
                                finalLinkage = 'ward.D2', 
                                distance = "spearman",
                                title= paste0('files/LC_Consensus_Cluster'), plot="pdf",  
                                seed = 123)


## Choose K
no_cluster <- 6
clust =  results[[no_cluster]]$consensusClass
names(clust) = colnames(proteomics)
  
  

## Choose K
## Consensus index
icl <- calcICL(results, plot = FALSE)

icl_samples <- as.data.frame(icl[[2]])
icl_samples <- subset(icl_samples, k == no_cluster)

icl_samples <- icl_samples[order(icl_samples$cluster), ]
samples <- as.character(unique(icl_samples$item))


# Normalize Consensus index to 1
icl_samples <- as.data.frame(t(sapply(samples, function(i) {
  dat <- subset(icl_samples, item == i)
  dat$itemConsensus / sum(dat$itemConsensus)
})))

colnames(icl_samples) <- paste0(1:no_cluster)
icl_samples$samples <- samples
icl_samples$cluster <- clust[icl_samples$samples]
#############################################################################################################
## This is code taken from ConsensusCluster R package (Wilkerson, M.D., Hayes, D.N. (2010). ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking. Bioinformatics, 2010 Jun 15;26(12):1572-3) to reproduce CDF and Delta area core figures
triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  
  
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  
  fm = t(nm)+nm
  diag(fm) = diag(m)
  
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
  
}

ml <- lapply(2:11, function(i) results[[i]][['ml']])
breaks=100
k=length(ml)
areaK = c()
hist_dat <- list()
for (i in 1:length(ml)){
  
  v=triangle(ml[[i]],mode=1)
  
  #empirical CDF distribution. default number of breaks is 100    
  h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
  h$counts = cumsum(h$counts)/sum(h$counts)
  
  hist_dat[[i]] <- data.frame(k=i+1, 'consensus_index' = h$mids, cumulative_sum = h$counts)
  #calculate area under CDF curve, by histogram method.
  thisArea=0
  for (bi in 1:(length(h$breaks)-1)){
    thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
    bi = bi + 1
  }
  areaK = c(areaK,thisArea)
  
}

# Plot area under CDF change.
deltaK=areaK[1] #initial auc at k=2
for(i in 2:(length(areaK))){
  #proportional increase relative to prior K.
  deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
}
############################################################################################################
### LC NMF clustering
# Split into two datasets
proteomics_plus <- proteomics
proteomics_plus[proteomics_plus < 0] <- 0

proteomics_negative <- proteomics
proteomics_negative[proteomics_negative > 0] <- 0
proteomics_negative <- abs(proteomics_negative)

proteomics_integrated <- rbind.data.frame(proteomics_plus, proteomics_negative)


## Perform runs for each value of r in range 2:11
estim.r <- nmf(x = proteomics_integrated, 
               rank = 2:11,
               method = 'brunet', 
               nrun=100, 
               seed=1234, 
               .opt = 'v')

# plot(estim.r)

## Cophenetic correlation
K = 6
nmf_clusters <- predict(estim.r$fit[[K-1]], what = 'consensus',  prob = TRUE)
LC_metadata$NMF <- nmf_clusters[match(LC_metadata$TMT.set.and.label, names(nmf_clusters))]

# Change labels
new_annotation <- c('1' = '6', '2' = '5', '3' = '4', '4' = '1', '5' = '2', '6' ='3')
LC_metadata$NMF <- new_annotation[match(LC_metadata$NMF, names(new_annotation))]
cross_tab <- as.data.frame(table(LC_metadata$Proteome.Subtype, LC_metadata$NMF))


# Cophonetic correlation
coph_cor <- sapply(1:10,function(i) {
  cophcor(estim.r$consensus[[i]])
})


## Membership index
# H matrix 
h <- as.data.frame(t(coef(estim.r$fit[[K-1]])))
h_membership <- as.data.frame(t(apply(h, 1, function(i) {
  i/sum(i)
})))

new_membership <- c('V1' = 'Subtype_6', 
                    'V2' = 'Subtype_4',
                    'V3' = 'Subtype_1', 
                    'V4' = 'Subtype_5', 
                    'V5' = 'Subtype_3',
                    'V6' = 'Subtype_2' )

## Choose K
colnames(h_membership) <- new_membership[colnames(h_membership)]
h_membership$samples <- rownames(h_membership)

h_membership$subtype <- colnames(h_membership)[apply(h, 1, which.max)]
h_membership$CC_subtype <- paste('Subtype',LC_metadata[h_membership$samples, 'Proteome.Subtype'],sep = '_')

###