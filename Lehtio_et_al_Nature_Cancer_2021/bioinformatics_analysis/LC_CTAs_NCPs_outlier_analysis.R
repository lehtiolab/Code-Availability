## LC Cancer-testis antigens / Non-canonical peptides outlier analysis


# Libraries
library(openxlsx) 
library(reshape2)
library(dplyr)


# Data
load('files/LC_colours')
peptide_per_protein <- read.delim('files/LC_symbols_set1_16_Peptides.txt')
proteomics <- read.table('files/LC_Landscape_Proteome.txt')
LC_metadata <- read.xlsx('files/Supplementary_table_1.xlsx', rowNames = TRUE)


# Cancer-testis antigen outlier detection
# CTdatabase (Almeida, L. G. et al. CTdatabase: a knowledge-base of high-throughput and curated data on cancer-testis antigens. Nucleic Acids Res 37, D816-819, doi:10.1093/nar/gkn673 (2009))
CT_dat <- read.delim('files/CT_database.txt', sep ='\t')
CT_dat_genes <- unique(CT_dat$Family.member)
CT_dat_genes[CT_dat_genes == 'cyclin A1'] <-  'CCNA1'
CT_dat_genes[CT_dat_genes == 'HIWI, MIWI, PIWI'] <-  'PIWIL1'

CT_dat_genes <- sapply(strsplit(CT_dat_genes, '/'), function(i) i[[1]])
CT_dat_genes[CT_dat_genes == "GOLGAGL2 FA"]  <- "GOLGAGL2"
CT_dat_genes[CT_dat_genes == "EDAG, NDR"] <- "HEMGN"
CT_dat_genes[CT_dat_genes == "LAGE-1b"] <- "CTAG2"
CT_dat_genes[CT_dat_genes == "IMP-3"] <- "IMP3"

CT_dat_genes <- toupper(CT_dat_genes)
CT_dat_genes <- unique(CT_dat_genes)


# Human protein atlas (www.proteinatlas.org) testis-enriched gene symbols
HPA_testis <- read.delim('files/HPA_Testis_Tissue.tsv', sep ='\t')


# Threshold
outlier_thres <- log2(3)

dat_3_a <- data.frame('Gene_symbol' = rownames(proteomics))
dat_3_a$CT_database <- ifelse(dat_3_a$Gene_symbol %in% CT_dat_genes, 'Yes', NA)
dat_3_a$Human_Protein_Atlas_Testis_Tissue <- ifelse(dat_3_a$Gene_symbol %in% HPA_testis$Gene, 'Yes', NA)
dat_3_a$Maximum_number_of_peptides <- peptide_per_protein$MaxPeps[match(dat_3_a$Gene_symbol, rownames(peptide_per_protein))]

outlyingSamples <- apply(proteomics, 1, function(i) any(i > outlier_thres, na.rm = TRUE))
outlyingSamples <- names(outlyingSamples)[outlyingSamples == 'TRUE']
dat_3_a$Outlier_abundance <- ifelse(dat_3_a$Gene_symbol %in% outlyingSamples, 'Yes', NA)

dat_3_a$Testis_specificity <- HPA_testis$RNA.tissue.distribution[match(dat_3_a$Gene_symbol, HPA_testis$Gene)]

dat_3_a$Filtered <- ifelse((dat_3_a$CT_database == 'Yes' | dat_3_a$Human_Protein_Atlas_Testis_Tissue == 'Yes') & dat_3_a$Maximum_number_of_peptides > 1 & dat_3_a$Outlier_abundance == 'Yes' & (dat_3_a$Testis_specificity %in% c("Detected in single", "Detected in some") | is.na(dat_3_a$Testis_specificity) ), 'Yes' ,NA)

dat_3_a <- dat_3_a[order(dat_3_a$Filtered), ]


cta_prot <- proteomics[dat_3_a$Gene_symbol[dat_3_a$Filtered == 'Yes' & !is.na(dat_3_a$Filtered)], ]


outlier_freq <- apply(cta_prot, 2,  function(i) sum(i > outlier_thres, na.rm = TRUE))

CTAs <- data.frame('cta_outliers' = outlier_freq, 
                     'Proteome_subtype' = LC_metadata$Proteome.Subtype[match(names(outlier_freq), LC_metadata$TMT.set.and.label)])
####################################################################################################
# Non-canonical peptides outlier detection
novelpeptides <- read.delim('files/LC_6ft_rerun_specaifix_peptides.txt', sep = '\t', row.names = 1)
novelpeptidesquant <- novelpeptides[, LC_metadata$TMT.set.and.label]


# Exclude non-quantified across all sets (non-informative)
idx <- apply(novelpeptidesquant, 1, function(i) sum(is.na(i)) == 141)
novelpeptidesquant <- novelpeptidesquant[!idx, ]
novelpeptides <- novelpeptides[!idx, ]

## Exclude mapped to known variants (SNPdb) and SpectrumAI fails by majority
spectrumAI_res <- sapply(1:nrow(novelpeptides), function(i) {
  idx <- novelpeptides$blastp_category[i] == 'map to known protein with 1 aa mismatch'
  if(!idx) {
    return('Not applicable')} else {
      dt_specAI <- novelpeptides[i, grep('SpectrumAI_result', colnames(novelpeptides))]
      
      # Remove set17
      dt_specAI <- dt_specAI[!grepl('set17', names(dt_specAI))]
      dt_specAI <- dt_specAI[order(names(dt_specAI))]
      dt_specAI <- as.character(unlist(dt_specAI))
      
      dt_specAI[is.na(dt_specAI)] = 'NA'
      dt_quant <- novelpeptides[i, colnames(novelpeptides) %in% LC_metadata$TMT]
      dt_quant <- suppressMessages(melt(dt_quant))
      dt_quant$set <- gsub('_.*', '', dt_quant$variable)
      dt_quant <- dt_quant %>% group_by(set) %>% summarise(n = mean(value, na.rm = TRUE))
      dt_quant <- dt_quant$n
      
      res <- table(dt_specAI[!is.na(dt_quant)])
      names(res)[which.max(res)]
    }
})

toexclude <- novelpeptides$blastp_category == 'match to known protein' |
    spectrumAI_res == 'FAIL'| novelpeptides$from.SNPdb != 'No match in SNP-DB'
  
novelpeptides <- novelpeptides[!toexclude, ]
novelpeptidesquant <- novelpeptidesquant[!toexclude, ]
  

# Normalize per set median
set_names <- sort(unique(gsub('_.*', '', colnames(novelpeptidesquant))))
  
novelpeptidesquant <- do.call(cbind, lapply(set_names, function(i) {
  res <- novelpeptidesquant[, grep(i, colnames(novelpeptidesquant))]
  res <- sweep(res, 1, apply(res, 1, median, na.rm = TRUE), '/')
  }))
  
  
common_samples <- intersect(LC_metadata$TMT, colnames(novelpeptidesquant))
novelpeptidesquant <- as.data.frame(novelpeptidesquant[, common_samples])
  
# Log2 transformation
novelpeptidesquant <- log2(novelpeptidesquant)

# Outliers
ncp_outlier_freq <- apply(novelpeptidesquant, 2,  function(i) sum(i > outlier_thres, na.rm = TRUE))
NCPs <- data.frame('ncp_outliers' = ncp_outlier_freq, 
                    'Proteome_subtype' = LC_metadata$Proteome.Subtype)

############################################################################################
## Annotate by annovar using genomic coordinates of peptides (hg19)
## Create input file
# annovar_input <- novelpeptides[, c("chr",   "start", "end" )]
# annovar_input[, c('ref', 'alt')] <- 0 

# write.table(annovar_input, 'files/LC_novpep_annovar.avinput', sep = ' ', row.names = FALSE, col.names = FALSE, quote = FALSE)


## Annovar
# Refseq, ucsc, ensembl, gencode hg19 
# annotate_variation.pl -out annovar_output_refseq -build hg19 LC_novpep_annovar.avinput --geneanno -dbtype refGene humandb/
# annotate_variation.pl -out annovar_output_UCSC -build hg19 LC_novpep_annovar.avinput --geneanno -dbtype knownGene humandb/
# annotate_variation.pl -out annovar_output_ensembl -build hg19 LC_novpep_annovar.avinput --geneanno -dbtype ensGene humandb/
# annotate_variation.pl -out annovar_output_gencode -build hg19 LC_novpep_annovar.avinput --geneanno -dbtype wgEncodeGencodeBasicV19 humandb/


## LncGencode
# annotate_variation.pl -regionanno -dbtype gff3 -gff3dbfile lncipedia_5_2_hc_hg19.gff.txt LC_novpep_annovar.avinput humandb/ -out  annovar_output_lincpediaLnc --build hg19
# annotate_variation.pl -regionanno -dbtype gff3 -gff3dbfile gencode.v34.long_noncoding_RNAsLifted.gff3 LC_novpep_annovar.avinput humandb/ -out  annovar_output_gencodeLnc --build hg19

## Pseudogenes
# annotate_variation.pl -regionanno -dbtype gff3 -gff3dbfile gencode.v34.2wayconspseudosLifted.gff3 LC_novpep_annovar.avinput humandb/ -out  annovar_output_gencodePseudo --build hg19


## Refseq, ucsc gencode ensembl hg19, 
annotation_names <- c('refseq', 'ucsc', 'gencode', 'ensembl')
annotation_files <- lapply(annotation_names, function(i) {
  
  res <- read.delim(paste0('files/annovar_output_', i,'.variant_function'), sep = '\t', header = FALSE, strip.white=TRUE)
  res$V1 <- as.character(res$V1)
  res
  
})

# Merge output
annotation_merged <- Reduce(cbind.data.frame, annotation_files)

# Remove superfluous columns
annotation_merged <- annotation_merged[, -seq(6,length(annotation_merged), 3)]


colnames(annotation_merged) <- c('refseq_annot', 'refseq_gene', 'pos', 'ucsc_annot', 'ucsc_gene', 'gencode_annot', 'gencode_gene', 'ensembl_annot', 'ensembl_gene')


# Lnc/Pseudo genes
gencodeLnc1 <- read.delim('files/annovar_output_gencodeLnc.hg19_gff3', sep = '\t', header = FALSE)
gencodeLnc2 <- read.delim('files/annovar_output_lincpediaLnc.hg19_gff3', sep = '\t', header = FALSE)
gencodeLnc <- rbind.data.frame(gencodeLnc1, gencodeLnc2)
gencodeLnc <- gencodeLnc[!duplicated(gencodeLnc$V3), ]

gencodePseudo <- read.delim('files/annovar_output_gencodePseudo.hg19_gff3', sep = '\t', header = FALSE)

annotation_merged$Lncgene_annot <- ifelse(annotation_merged$pos %in% gencodeLnc$V3, "lncrna", NA)
annotation_merged$Pseudo_annot <- ifelse(annotation_merged$pos %in% gencodePseudo$V3, "pseudogene", NA)


# AltORf prediction
altorf <- read.delim('files/LC_novpep_thres6ft_annotated_colpeps_protinfo.tsv', sep = '\t')
altorf <- altorf[!is.na(altorf$ProteinInfo), ]

# Re-annotate
# term_freq <- table(gsub('_.*', '', unlist(strsplit(as.character(altorf$ProteinInfo), ','))))c

altorf$meta <- unlist(sapply(altorf$ProteinInfo, function(i) {
  
  res <- as.character(i)
  res_freq <- sapply(strsplit(res, ','), function(i) gsub('_.*', '', i))
  ifelse(all(res_freq == 'altorf'), 'AltOrf', 'Rest')
}))

altorf <- altorf[altorf$meta == 'AltOrf', ]

annotation_merged$pep <- rownames(novelpeptides)

annotation_merged$AltOrf_annot[gsub('[^A-Za-z]', '', annotation_merged$pep) %in% altorf$Peptide] <-  'AltOrf'


##  Retroelements using matching protein names in blastp from uniprot 
uniprot <- read.delim('files/uniprot-filtered-organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+review--.tab', sep = '\t')
uniprot_retro <- uniprot[grep(paste(c('Endogenous retrovirus', 'Retrotransposon', 'retrotransposable', 'transposon'), collapse = '|'), uniprot$Protein.names), ]


annotation_merged$ERV_annot <- NA
annotation_merged$ERV_annot[grep(paste0(uniprot_retro$Entry.name, collapse = '|'),
                            novelpeptides[annotation_merged$pep, 'blastp_match'])] <- 'ERV'

# Precedence rules
final_rules <- unique(unlist(annotation_merged[, grep('annot', colnames(annotation_merged))]))
final_rules <- final_rules[!is.na(final_rules)]

# Re-order
final_rules <- c('AltOrf', "ERV", "pseudogene", "exonic", "splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_intronic", "lncrna", "UTR5", "UTR3", "UTR5;UTR3", "intronic", "upstream", "downstream", "upstream;downstream", "intergenic")

annotation_merged$meta_annot <- apply(annotation_merged[, grep('annot', colnames(annotation_merged))], 1, function(i) {
  final_rules[min(which(final_rules %in% unlist(i)))]}
)

annotation_merged$no_outliers <- sapply(annotation_merged$pep, function(i) {
  sum(novelpeptidesquant[i, ] > outlier_thres, na.rm = TRUE) 
  })

####################################################################################################
###  NCP correlation with metadata
LC_metadata$MKI67 <-  as.numeric(proteomics['MKI67', ])
LC_metadata$TP53 <- factor(LC_metadata$TP53_mut, levels = c(0,1), labels = c('wt', 'Mut'))
# LC_metadata$Proteomics <- as.factor(LC_metadata$Proteomics)
LC_metadata$TMB <- log2(LC_metadata$TMB)
LC_metadata$NCPs <- log2(NCPs$ncp_outliers)

# Linear models
lm1 <- lm('NCPs ~ TMB + Purity + TP53 + MKI67' , data = LC_metadata)
res_lm <- summary(lm1)
#############################################################################################################
## Neoantigen burden calculation
# Min-max normalization
min.max.norm<- function(x)
{
  return((x- min(x, na.rm = TRUE)) /(max(x, na.rm = TRUE)-min(x,  na.rm = TRUE)))
}

TMB_CTA_NCP <- data.frame('Proteome_subytpe' = LC_metadata$Proteome.Subtype,
                          'TMB' = LC_metadata$TMB, 
                          'CTA' = CTAs$cta_outliers,
                          'NCP' = NCPs$ncp_outliers)

TMB_CTA_NCP[, paste0('score_', c('TMB', 'CTA','NCP'))] <- as.data.frame(apply(TMB_CTA_NCP[, c('TMB', "CTA",'NCP')],2, min.max.norm))


# Sum and min-max normalization
TMB_CTA_NCP$Neoantigen_burden <- min.max.norm(apply(TMB_CTA_NCP[, c("score_TMB", "score_CTA", "score_NCP")], 1, sum))


# Median per subtype
neoantigen_burden_median <- as.data.frame(TMB_CTA_NCP[,c('Proteome_subytpe',"score_TMB", "score_CTA", "score_NCP", 'Neoantigen_burden')] %>%
                                     group_by(Proteome_subytpe) %>% 
                                     summarise_all(median, na.rm = TRUE))

###