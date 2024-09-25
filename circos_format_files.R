# Functions and code that will run to format the input files to be ready for the 
# make_circos.CMVD.rmd

rm(list = ls())

################################################################################
#                                 functions                                    #
################################################################################

#' Formats data for gene input
#' @param data data.frame containing biofilter output

#data = genes_gr_c
format_genes <- function(data){
  
  data <- as.data.frame(data)
  data = data[order(data$P),]
  data2 <- data[,c(2,3,3,16,24)] # label = search gene (not gene)
  data2$POS.1 <- data2$POS.1 + 1
  names(data2)[1:5] <- c("chr", "start", "end", "label","group") 
  data2$chr <- paste0("chr", data2$chr)
  data2$start <- as.numeric(data2$start)
  data3 = unique(data2)
  data3 = data3[!duplicated(data3$label), ]
  
  data.bed <- data3[data3$chr != "chrNA", ]
  
  return(data.bed)
}


format_manhattan <- function (data, annotations, 
                             p.threshold,
                             col = "black",
                             anno.col = "#009E73"){
  
  
  data <- as.data.frame(data)
  annotations <- as.data.frame(annotations)
  # data$BP <- data$position
  # data$CHR <- data$chromosome
  data$end <- data$BP + 1
  data$value1 <- -log10(data$pvalue)
  data$value2 <- col
  
  data.fmt <- data[data$value1 > -log10(p.threshold) | data$value2 == 1,
                   c("CHR", "BP", "end", "value1", "value2")]
  
  names(data.fmt)[1:2] <- c("chr", "start")
  
  # data.fmt$chr <- paste0("chr", data.fmt$chr)
  
  return(data.fmt)  
}

################################################################################
#                                 file paths                                   #
################################################################################

# inputs 
# Four input files here
gwas.path <- file.path("")
twas.path <- file.path("")
pwas.path <- file.path("")
anno.path <- file.path("")

# 4 outputs
genes.out <- file.path("", "gene_track_data_hfm")
manhattan.out <- file.path("", "manhattan_data_hfm")
twas.out <- file.path("", "twas_data_hfm")
pwas.out <- file.path("", "pwas_data_hfm")

################################################################################
#                               loading zone                                   #
################################################################################

# Load the data for the GWAS annotation
gwas.data <- data.table::fread(gwas.path)
colnames(gwas.data) <- gsub("#", "", colnames(gwas.data))
#gene.data <- data.table::fread(gene.path)

anno.data2 <- data.table::fread(anno.path)
# colnames(anno.data2) <- gsub("#", "", colnames(anno.data2))
#anno.data <- readxl::read_excel(anno.path)

twas.data <- data.table::fread(twas.path)
pwas.data <- data.table::fread(pwas.path)
#rsid <- data.table::fread(rsid.path)



################################################################################
#                            format for gene track                             #
################################################################################

genes_tp = anno.data2
# genes_tp$chr <- sub("^", "chr", genes_tp$CHR)
genes_tp$end = genes_tp$start + 1
# genes_tp$label = genes_tp$ANNOTATION
# genes_tp$group = genes_tp$PHENOTYPE
genes.bed <- genes_tp[,c(6,7,8,2,3)]

# Create final formatted dataset
genes.bed$start = as.numeric(genes.bed$start)
genes.bed$end = as.numeric(genes.bed$end)


save(genes.bed, file = genes.out)

################################################################################
#                      format for manhattan plot track                         #
################################################################################

mega.bed_gr <- format_manhattan(data = gwas.data,
                                annotations = anno.data2,
                                p.threshold = 1) # changed from 0.01 to 1 to include all SNPs
save(mega.bed_gr, file = manhattan.out)



################################################################################
#                              format for twas                              #
################################################################################
# p-value threshold
twas_gr = twas.data[twas.data$pvalue<=0.05,]
twas_gr4 = twas_gr[,c("gene_name_x", "Tissue", "pvalue", "gene_chr", "gene_start")]
twas_gr4$log10p <- -log10(twas_gr4$pvalue)
# Flip the orientation for the plot
twas_gr4$end = twas_gr4$gene_start + 1
twas_gr4$tissue = ifelse(twas_gr4$Tissue=="Heart_Left_Ventricle", "pink", 
                         ifelse(twas_gr4$Tissue=="Artery_Aorta", "darkgreen",
                                ifelse(twas_gr4$Tissue=="Heart_Atrial_Appendage", "red",
                                       ifelse(twas_gr4$Tissue=="Artery_Tibial", "darkred",
                                              ifelse(twas_gr4$Tissue=="Whole_Blood", "maroon",
                                                     ifelse(twas_gr4$Tissue=="Adipose_Visceral_Omentum", "white",
                                                            ifelse(twas_gr4$Tissue=="Adipose_Subcutaneous", "magenta",
                                                                   ifelse(twas_gr4$Tissue=="Kidney_Cortex", "orange",
                                                                          ifelse(twas_gr4$Tissue=="Artery_Coronary", "purple",
                                                                                 ifelse(twas_gr4$Tissue=="Liver", "#FFC107", ""))))))))))


twas_gr5 = twas_gr4[,c(4,5,7,6,8)]

colnames(twas_gr5) <- c("chr", "start", "end", "value1", "value2")
twas_gr5$value1 <- twas_gr5$value1 * (-1)


save(twas_gr5, file = twas.out)

################################################################################
#                              format for pwas                              #
################################################################################
# pwas_gr = pwas.data[grep(pheno, pwas.data$phenotype),]
# p-value threshold
pwas_gr = pwas.data[pwas.data$pvalue<=1,]
pwas_gr4 = pwas_gr[,c("gene_name_x", "Model", "pvalue", "gene_chr", "gene_start")]
pwas_gr4$log10p <- -log10(pwas_gr$pvalue)
# Flip the orientation for the plot
pwas_gr4$end_position = pwas_gr4$gene_start + 1
pwas_gr4$tissue = ifelse(pwas_gr4$Model=="EA", "darkblue", 
                         ifelse(pwas_gr4$Model=="AA", "lightblue", ""))


pwas_gr5 = pwas_gr4[,c(4,5,7,6,8)]
colnames(pwas_gr5) <- c("chr", "start", "end", "value1", "value2")
pwas_gr5$value1 <- pwas_gr5$value1 * (-1)


save(pwas_gr5, file = pwas.out)

