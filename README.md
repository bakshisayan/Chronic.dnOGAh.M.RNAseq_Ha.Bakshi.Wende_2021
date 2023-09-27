## Chronic.dnOGAh.M.RNAseq_Ha.Bakshi.Wende_2021 ##

### Author: 
Sayan Bakshi, M.Sc. (Integrated) | PhD Trainee
### First Author of the manuscript: 
Dr. Chae-Myeong Ha, PhD
### Corresponding Author / PI: 
Dr. Adam R. Wende, PhD, Department of Pathology, University of Alabama at Birmingham, Birmingham, Alabama, USA

## Introduction
TBD

## Methods 
TBD 

## Loading Libraries 
```r
## Setting file path ##
set.seed(011)
file.path="~/Box/200420 GlcNAc dnOGAh and iOGT/_Sayan.analysis/5_Manuscript/Final_Analysis/_Scripts/"
setwd(file.path)

## Loading packages ##
basic.packages <- c("dplyr", "tidyr", "openxlsx")
viz.packages <- c("pheatmap", "RColorBrewer", "FactoClass")
RNAseq.packages <- c("biomaRt", "DESeq2", "clusterProfiler", "ReactomePA", 
                     "enrichplot")
genome.packages <- c("org.Mm.eg.db")
packages <- c(basic.packages, viz.packages, RNAseq.packages, genome.packages)

  #install.packages("pacman")
pacman::p_load(packages, character.only = TRUE)
```
## Input Sample Data and Count Data
```r
dir.create("../1_Input")
dir.create("../2_Output")
dir.create("../3_Results")

mypathOP <- "../2_Output/"
mypathRes <- "../3_Results/"

## RNASeq Data Input
#(Code modified from https://www.biostars.org/p/241602/)
raw.input <- list.files(path = "../1_Input/", 
                       pattern = "*ReadsPerGene.out.tab$", full.names = T)
counts.files <- lapply(raw.input, read.table)
counts <- as.data.frame(sapply(counts.files, function(x)x[,4]))
raw.input.2 <- gsub("../1_Input//", "", raw.input)
raw.input.2 <- gsub("[_][A-Z][0-9]+[_][A-Z][0-9]+[_][0-9]+[.]ReadsPerGene[.]out[.]tab", "", raw.input.2)
colnames(counts) <- raw.input.2
row.names(counts) <- counts.files[[1]]$V1
dnOGAh_counts <- counts[-c(1:4),]
dnOGAh_counts <- dnOGAh_counts[,order(colnames(dnOGAh_counts))]

## Sample (Mus musculus) data input
sample.info <- read.xlsx("../1_Input/Sample.Information.xlsx",
                         sheet="Sample.Data_Final")
sample.info <- sample.info[order(sample.info$Name),]
sample.info$Genotype <- factor(sample.info$Genotype,levels = c("Con", "dnOGAh"))
sample.info$Project <- factor(sample.info$Project,levels = c("2wk_ON", "6mo_ON"))
sample.info$Duration <- factor(sample.info$Duration,levels = c("2wk", "24wk"))
sample.info$Group <- paste(sample.info$Duration, sample.info$Genotype, sep = "_")
sample.info$Group <- factor(sample.info$Group, 
                            levels = c("2wk_Con", "2wk_dnOGAh", "24wk_Con", "24wk_dnOGAh"))

#For M 2wk ON: dnOGAh vs Con
ON_2wk.M <- sample.info[sample.info$Project=="2wk_ON",] 
#For M 24wk ON: dnOGAh vs Con
ON_6mo.M <- sample.info[sample.info$Project=="6mo_ON",] 

## Annotation Input ##
## Organism: Mouse (Mus musculus)
mm39 <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id", "ensembl_gene_id_version", "strand",
                         "refseq_mrna", "external_gene_name", "chromosome_name", 
                         "start_position", "end_position", "entrezgene_id", 
                         "entrezgene_description"),  mart=mm39)
```
## Differential Analysis for Changes in dnOGAh mice relative to Control mice 
```r
## Comparison #1: dnOGAh vs. Con | Male | 2 weeks:
phenodata = ON_2wk.M

## DESeq2 Dataset preparation by CountMatrix Input:
phenodata <- data.frame(phenodata)
cts <- as.matrix(dnOGAh_counts[,phenodata$Name])
coldata <- phenodata[,c(2:ncol(phenodata))]
rownames(coldata) <- phenodata$Name

## DESeq Matrix:
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, 
                              design = ~Genotype)

## DESeq Analysis: (Wald)
dds <- DESeq(dds)

## Prefiltering:
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

## Extracting transformed values
rld <- rlog(dds, blind=FALSE)

## Annotation & Data Filtering
res <- results(dds, name="Genotype_dnOGAh_vs_Con")

## Normalized Counts:
Norm_counts <- counts(dds, normalized=T)
res$ensembl_gene_id <- rownames(res)

## New Results Table with Normalized counts: (not trimmed)
res_2 <- merge(res, Norm_counts, by=0)
rownames(res_2) <- res_2$Row.names
res_2$ensembl_gene_id <- gsub("\\..*","",res_2$ensembl_gene_id)
res_tmp.0 <- res_2[,c(2:8)]
res_tmp <- res_2[,c(9:ncol(res_2))]
res_select <- res_tmp[which(colnames(res_tmp) %in% phenodata$Name)]
res_select$ensembl_gene_id <- rownames(res_select)
res_select$ensembl_gene_id <- gsub("\\..*","",res_select$ensembl_gene_id)
res_3 <- merge(res_tmp.0, res_select, by="ensembl_gene_id")

res_annot <- merge(bm, res_3, by=c("ensembl_gene_id")) %>% mutate_all(na_if,"")
res_clean <- res_annot[!duplicated(res_annot$ensembl_gene_id),]

## Filtering
res_q0.1FC1.5 <- dplyr::filter(res_clean, abs(log2FoldChange)>log2(1.5) & padj<0.1)
nrow(res_q0.1FC1.5)
res_q0.1FC1.5.Up <- dplyr::filter(res_q0.1FC1.5, log2FoldChange > 0)
nrow(res_q0.1FC1.5.Up)
res_q0.1FC1.5.Down <- dplyr::filter(res_q0.1FC1.5, log2FoldChange < 0)
nrow(res_q0.1FC1.5.Down)

DEG_2wk_All <- res_clean
DEG_2wk_FDR <- res_q0.1FC1.5
DEG_2wk_FDR.Up <- res_q0.1FC1.5.Up
DEG_2wk_FDR.Down <- res_q0.1FC1.5.Down
  
## Comparison #2: dnOGAh vs. Con | Male | 24 weeks:
phenodata = ON_6mo.M

## DESeq2 Dataset preparation by CountMatrix Input:
phenodata <- data.frame(phenodata)
cts <- as.matrix(dnOGAh_counts[,phenodata$Name])
coldata <- phenodata[,c(2:ncol(phenodata))]
rownames(coldata) <- phenodata$Name

## DESeq Matrix:
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, 
                              design = ~Genotype)

## DESeq Analysis: (Wald)
dds <- DESeq(dds)

## Prefiltering:
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

## Extracting transformed values
rld <- rlog(dds, blind=FALSE)

## Annotation & Data Filtering
res <- results(dds, name="Genotype_dnOGAh_vs_Con")

## Normalized Counts:
Norm_counts <- counts(dds, normalized=T)
res$ensembl_gene_id <- rownames(res)

## New Results Table with Normalized counts: (not trimmed)
res_2 <- merge(res, Norm_counts, by=0)
rownames(res_2) <- res_2$Row.names
res_2$ensembl_gene_id <- gsub("\\..*","",res_2$ensembl_gene_id)
res_tmp.0 <- res_2[,c(2:8)]
res_tmp <- res_2[,c(9:ncol(res_2))]
res_select <- res_tmp[which(colnames(res_tmp) %in% phenodata$Name)]
res_select$ensembl_gene_id <- rownames(res_select)
res_select$ensembl_gene_id <- gsub("\\..*","",res_select$ensembl_gene_id)
res_3 <- merge(res_tmp.0, res_select, by="ensembl_gene_id")

res_annot <- merge(bm, res_3, by=c("ensembl_gene_id")) %>% mutate_all(na_if,"")
res_clean <- res_annot[!duplicated(res_annot$ensembl_gene_id),]

## Filtering
res_q0.1FC1.5 <- dplyr::filter(res_clean, abs(log2FoldChange)>log2(1.5) & padj<0.1)
nrow(res_q0.1FC1.5)
res_q0.1FC1.5.Up <- dplyr::filter(res_q0.1FC1.5, log2FoldChange > 0)
nrow(res_q0.1FC1.5.Up)
res_q0.1FC1.5.Down <- dplyr::filter(res_q0.1FC1.5, log2FoldChange < 0)
nrow(res_q0.1FC1.5.Down)

DEG_24wk_All <- res_clean
DEG_24wk_FDR <- res_q0.1FC1.5
DEG_24wk_FDR.Up <- res_q0.1FC1.5.Up
DEG_24wk_FDR.Down <- res_q0.1FC1.5.Down
```
## Heatmap: All Genes Combined
```r
## Combined Data: dnOGAh vs. Con Analysis for both 2 weeks and 24 weeks
DEG_Geno.2wk24wk_Merge <- merge(DEG_2wk_All, DEG_24wk_All, by="ensembl_gene_id")
names(DEG_Geno.2wk24wk_Merge) <- gsub(names(DEG_Geno.2wk24wk_Merge), 
                                        pattern = "\\.x", replacement = "\\.Geno.2wk") %>%
                                   gsub(names(DEG_Geno.2wk24wk_Merge), 
                                        pattern = "\\.y", replacement = "\\.Geno.24wk")
selection <- c(names(DEG_Geno.2wk24wk_Merge[,2:10])) %>% 
  str_replace(pattern = "\\.Geno.2wk", replacement = "")
names(DEG_Geno.2wk24wk_Merge)[2:10] <- selection
 
## Heatmap
hm <- data.matrix(DEG_Geno.2wk24wk_Merge%>%dplyr::select(starts_with("dnOGAh")))
row.names(hm) <- DEG_Geno.2wk24wk_Merge$ensembl_gene_id
ann_col <- sample.info%>%dplyr::select(Group)
row.names(ann_col) <- sample.info$Name
ann_colors <- list(Group = c(`2wk_Con`="darkgreen", `2wk_dnOGAh`="black",
                             `24wk_Con`="blue", `24wk_dnOGAh`="red"))

callback <-  function(hc, mat){
sv = svd(t(mat))$v[,1]
dend = reorder(as.dendrogram(hc), wts = sv)
as.hclust(dend)
}

hm.Plot <- pheatmap(hm, cluster_rows=T, show_rownames=F,
         cluster_cols=T, annotation_colors = ann_colors, 
         clustering_callback = callback, 
         kmeans_k = 1000, cutree_rows = 2, cutree_cols = 2,
         annotation_col=ann_col, scale="row", fontsize = 7, 
         main = "All Genes: 2wk and 24wk ON.M_dnOGAh.v.Con",
         colorRampPalette(c("blue4", "white", "darkred"))(50))
hm.Plot

```
## Principal Component Analysis (PCA) Plot: All Genes Combined
```r
phenodata <- data.frame(sample.info)
cts <- as.matrix(dnOGAh_counts[,phenodata$Name])
coldata <- phenodata[,c(2:ncol(phenodata))]
rownames(coldata) <- phenodata$Name
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, 
                              design = ~Group)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
rld <- rlog(dds, blind=FALSE)

mat <- assay(rld)
pca <- prcomp(t(mat))
Var <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100*Var, digits = 2)

pca <- as.data.frame(pca$x)
pca$Name <- rownames(pca)
pca.comp <- merge(pca, phenodata, by="Name")
row.names(pca.comp) <- pca.comp$Name
pca.comp$Group <- factor(pca.comp$Group,
                         levels = c("2wk_Con", "2wk_dnOGAh", "24wk_Con", "24wk_dnOGAh"))

shapes <- c(15,16,17,18) 
shapes <- shapes[as.numeric(pca.comp$Group)]
colors <- c("darkgreen", "black", "blue", "red")
colors <- colors[as.numeric(pca.comp$Group)]

s3d <- scatterplot3d(pca.comp[,2:4], grid=F, box=F,
                     main="3D PCA Plot",
                     xlab = paste0("PC1: ",percentVar[1],"% variance"),
                     ylab = paste0("PC2: ",percentVar[2],"% variance"),
                     zlab = paste0("PC3: ",percentVar[3],"% variance"))
addgrids3d(pca.comp[,2:4], grid = c("xy", "xz", "yz"), angle = 40)
s3d$points3d(pca.comp[,2:4], pch = shapes, col = colors, cex=2)
legend("bottom", legend = levels(pca.comp$Group),
       col = c("darkgreen", "black", "blue", "red"), cex = 1,
       pch = c(15,16,17,18), inset = -0.25, xpd = T, horiz = T)

```
## CNET Plot from Pathway Enrichment Analysis against Reactome Database
```r
## Combining Q < 0.1 and |FC| > 1.5 -significant DEGs:
DEG_Merge_FDR <- merge(DEG_2wk_FDR, DEG_24wk_FDR, by="ensembl_gene_id")
names(DEG_Merge_FDR) <- gsub(names(DEG_Merge_FDR), 
                                        pattern = "\\.x", replacement = "\\.Geno.2wk") %>%
                                   gsub(names(DEG_Merge_FDR), 
                                        pattern = "\\.y", replacement = "\\.Geno.24wk")
selection <- c(names(DEG_Merge_FDR[,2:10])) %>% 
  str_replace(pattern = "\\.Geno.2wk", replacement = "")
names(DEG_Merge_FDR)[2:10] <- selection

DEG_Merge_FDR.Up <- DEG_Merge_FDR %>% 
  filter(log2FoldChange.Geno.2wk>0 & log2FoldChange.Geno.24wk>0)
DEG_Merge_FDR.Inv <- DEG_Merge_FDR %>% 
  filter((log2FoldChange.Geno.2wk>0 & log2FoldChange.Geno.24wk<0)|
           (log2FoldChange.Geno.2wk<0 & log2FoldChange.Geno.24wk>0))
DEG_Merge_FDR.Down <- DEG_Merge_FDR %>% 
  filter(log2FoldChange.Geno.2wk<0 & log2FoldChange.Geno.24wk<0)

DEG_Geno.2wk.only <- anti_join(DEG_2wk_FDR, DEG_24wk_FDR, by="ensembl_gene_id")
DEG_Geno.2wk.only.Up <- DEG_Geno.2wk.only %>% filter(log2FoldChange>0)
DEG_Geno.2wk.only.Down <- DEG_Geno.2wk.only %>% filter(log2FoldChange<0)

DEG_Geno.24wk.only <- anti_join(DEG_24wk_FDR, DEG_2wk_FDR, by="ensembl_gene_id")
DEG_Geno.24wk.only.Up <- DEG_Geno.24wk.only %>% filter(log2FoldChange>0)
DEG_Geno.24wk.only.Down <- DEG_Geno.24wk.only %>% filter(log2FoldChange<0)

## Pathway Enrichment Analysis against Reactome Database:
## CNET Plot: DEG - Changes in dnOGAh mice relative to Con ONLY at 24 weeks
res_q0.1FC1.5 <- DEG_Geno.24wk.only 
reactome.data <- res_q0.1FC1.5%>%dplyr::select(entrezgene_id, log2FoldChange)
reactome.data <- na.omit(reactome.data)
reactome.list <- reactome.data$log2FoldChange
names(reactome.list) <- reactome.data$entrezgene_id

DEG.enrich <- enrichPathway(gene = names(reactome.list), 
                               organism = "mouse", pvalueCutoff=0.05,
                               qvalueCutoff=0.1, minGSSize =1, readable = T)
min = min(res_q0.1FC1.5$log2FoldChange)
max = max(res_q0.1FC1.5$log2FoldChange)

cnetplot(DEG.enrich, categorySize="pvalue", foldChange = reactome.list,
         colorEdge=T, circular=F, node_label="all", cex_category= 1,
         cex_gene=1, cex_label_category=1.2, cex_label_gene=0.8, layout="kk") + 
  ggplot2::scale_color_gradient2(midpoint=0, low="darkblue", mid="khaki1",
                                   high="firebrick4", limits=c(min,max)) +
  labs(color='log2FoldChange')

## Pathway Enrichment Analysis against Reactome Database:
## CNET Plot: DEG - Downregulated Changes in dnOGAh mice relative to Con ONLY at 24 weeks
res_q0.1FC1.5 <- DEG_Geno.24wk.only.Down
reactome.data <- res_q0.1FC1.5%>%dplyr::select(entrezgene_id, log2FoldChange)
reactome.data <- na.omit(reactome.data)
reactome.list <- reactome.data$log2FoldChange
names(reactome.list) <- reactome.data$entrezgene_id

DEG.enrich <- enrichPathway(gene = names(reactome.list), 
                               organism = "mouse", pvalueCutoff=0.05,
                               qvalueCutoff=0.1, minGSSize =1, readable = T)
min = min(res_q0.1FC1.5$log2FoldChange)
max = max(res_q0.1FC1.5$log2FoldChange)

cnetplot(DEG.enrich, categorySize="pvalue", foldChange = reactome.list,
         colorEdge=T, circular=F, node_label="all", cex_category= 1,
         cex_gene=1, cex_label_category=1.2, cex_label_gene=0.8, layout="kk") + 
  ggplot2::scale_color_gradient2(midpoint=0, low="darkblue", mid="khaki1",
                                   high="firebrick4", limits=c(min,max)) +
  labs(color='log2FoldChange')

```
## Session Information
```r
sessioninfo::session_info()
```
| ─ Session info ───────────────────────────────────────────────────────────────────────────────── 	|
|--------------------------------------------------------------------------------------------------	|
| R version 4.0.3 (2020-10-10)                                                                     	|
| Platform: x86_64-apple-darwin17.0 (64-bit)                                                       	|
| Running under: macOS Big Sur 10.16                                                               	|
|                                                                                                  	|
|                                                                                                  	|
| Matrix products: default                                                                         	|
| LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib              	|
|                                                                                                  	|
|                                                                                                  	|
| locale:                                                                                          	|
| [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8                                	|
|                                                                                                  	|
|                                                                                                  	|
| ─ Packages ───────────────────────────────────────────────────────────────────────────────────── 	|
| package              * version     date       lib source                                         	|
| abind                  1.4-5       2016-07-21 [1] CRAN (R 4.0.2)                                 	|
| ade4                 * 1.7-18      2021-09-16 [1] CRAN (R 4.0.2)                                 	|
| affy                   1.68.0      2020-10-27 [1] Bioconductor                                   	|
| affyio                 1.60.0      2020-10-27 [1] Bioconductor                                   	|
| annotate               1.68.0      2020-10-27 [1] Bioconductor                                   	|
| AnnotationDbi        * 1.52.0      2020-10-27 [1] Bioconductor                                   	|
| askpass                1.1         2019-01-13 [1] CRAN (R 4.0.2)                                 	|
| assertthat             0.2.1       2019-03-21 [1] CRAN (R 4.0.2)                                 	|
| backports              1.2.1       2020-12-09 [1] CRAN (R 4.0.2)                                 	|
| beeswarm               0.3.1       2021-03-07 [1] CRAN (R 4.0.2)                                 	|
| Biobase              * 2.50.0      2020-10-27 [1] Bioconductor                                   	|
| BiocFileCache          1.14.0      2020-10-27 [1] Bioconductor                                   	|
| BiocGenerics         * 0.36.1      2021-04-16 [1] Bioconductor                                   	|
| BiocManager            1.30.15     2021-05-11 [1] CRAN (R 4.0.2)                                 	|
| BiocParallel           1.24.1      2020-11-06 [1] Bioconductor                                   	|
| biomaRt              * 2.46.3      2021-02-11 [1] Bioconductor                                   	|
| Biostrings             2.58.0      2020-10-27 [1] Bioconductor                                   	|
| bit                    4.0.4       2020-08-04 [1] CRAN (R 4.0.2)                                 	|
| bit64                  4.0.5       2020-08-30 [1] CRAN (R 4.0.2)                                 	|
| bitops                 1.0-7       2021-04-24 [1] CRAN (R 4.0.2)                                 	|
| blob                   1.2.1       2020-01-20 [1] CRAN (R 4.0.2)                                 	|
| broom                  0.7.6       2021-04-05 [1] CRAN (R 4.0.2)                                 	|
| cachem                 1.0.5       2021-05-15 [1] CRAN (R 4.0.2)                                 	|
| car                    3.0-10      2020-09-29 [1] CRAN (R 4.0.2)                                 	|
| carData                3.0-4       2020-05-22 [1] CRAN (R 4.0.2)                                 	|
| cellranger             1.1.0       2016-07-27 [1] CRAN (R 4.0.2)                                 	|
| checkmate              2.0.0       2020-02-06 [1] CRAN (R 4.0.2)                                 	|
| circlize             * 0.4.12      2021-01-08 [1] CRAN (R 4.0.3)                                 	|
| cli                    2.5.0       2021-04-26 [1] CRAN (R 4.0.2)                                 	|
| clusterProfiler      * 3.18.1      2021-02-11 [1] Bioconductor                                   	|
| colorspace             2.0-1       2021-05-04 [1] CRAN (R 4.0.2)                                 	|
| cowplot                1.1.1       2020-12-30 [1] CRAN (R 4.0.2)                                 	|
| crayon                 1.4.1       2021-02-08 [1] CRAN (R 4.0.2)                                 	|
| curl                   4.3.1       2021-04-30 [1] CRAN (R 4.0.2)                                 	|
| data.table             1.14.0      2021-02-21 [1] CRAN (R 4.0.2)                                 	|
| DBI                    1.1.1       2021-01-15 [1] CRAN (R 4.0.2)                                 	|
| dbplyr                 2.1.1       2021-04-06 [1] CRAN (R 4.0.2)                                 	|
| DelayedArray           0.16.3      2021-03-24 [1] Bioconductor                                   	|
| DESeq2               * 1.30.1      2021-02-19 [1] Bioconductor                                   	|
| digest                 0.6.27      2020-10-24 [1] CRAN (R 4.0.2)                                 	|
| DO.db                  2.9         2020-12-18 [1] Bioconductor                                   	|
| DOSE                   3.16.0      2020-10-27 [1] Bioconductor                                   	|
| downloader             0.4         2015-07-09 [1] CRAN (R 4.0.2)                                 	|
| dplyr                * 1.0.6       2021-05-05 [1] CRAN (R 4.0.2)                                 	|
| ellipsis               0.3.2       2021-04-29 [1] CRAN (R 4.0.2)                                 	|
| enrichplot           * 1.10.2      2021-01-28 [1] Bioconductor                                   	|
| enrichR              * 3.0         2021-02-02 [1] CRAN (R 4.0.2)                                 	|
| eulerr               * 6.1.0       2020-03-09 [1] CRAN (R 4.0.2)                                 	|
| evaluate               0.14        2019-05-28 [1] CRAN (R 4.0.1)                                 	|
| FactoClass           * 1.2.7       2018-10-01 [1] CRAN (R 4.0.2)                                 	|
| fansi                  0.4.2       2021-01-15 [1] CRAN (R 4.0.2)                                 	|
| farver                 2.1.0       2021-02-28 [1] CRAN (R 4.0.2)                                 	|
| fastmap                1.1.0       2021-01-25 [1] CRAN (R 4.0.2)                                 	|
| fastmatch              1.1-0       2017-01-28 [1] CRAN (R 4.0.2)                                 	|
| fgsea                * 1.16.0      2020-10-27 [1] Bioconductor                                   	|
| forcats                0.5.1       2021-01-27 [1] CRAN (R 4.0.2)                                 	|
| foreign                0.8-80      2020-05-24 [2] CRAN (R 4.0.3)                                 	|
| genefilter             1.72.1      2021-01-21 [1] Bioconductor                                   	|
| geneplotter            1.68.0      2020-10-27 [1] Bioconductor                                   	|
| generics               0.1.0       2020-10-31 [1] CRAN (R 4.0.2)                                 	|
| GenomeInfoDb         * 1.26.7      2021-04-08 [1] Bioconductor                                   	|
| GenomeInfoDbData       1.2.4       2020-12-18 [1] Bioconductor                                   	|
| GenomicRanges        * 1.42.0      2020-10-27 [1] Bioconductor                                   	|
| ggbeeswarm           * 0.6.0       2017-08-07 [1] CRAN (R 4.0.2)                                 	|
| ggforce                0.3.3       2021-03-05 [1] CRAN (R 4.0.2)                                 	|
| ggfortify            * 0.4.11      2020-10-02 [1] CRAN (R 4.0.2)                                 	|
| ggnewscale           * 0.4.5       2021-01-11 [1] CRAN (R 4.0.2)                                 	|
| ggplot2              * 3.3.3       2020-12-30 [1] CRAN (R 4.0.2)                                 	|
| ggpubr               * 0.4.0       2020-06-27 [1] CRAN (R 4.0.2)                                 	|
| ggraph                 2.0.5       2021-02-23 [1] CRAN (R 4.0.2)                                 	|
| ggrepel              * 0.9.1       2021-01-15 [1] CRAN (R 4.0.2)                                 	|
| ggsignif               0.6.1       2021-02-23 [1] CRAN (R 4.0.2)                                 	|
| GlobalOptions          0.1.2       2020-06-10 [1] CRAN (R 4.0.2)                                 	|
| glue                   1.4.2       2020-08-27 [1] CRAN (R 4.0.2)                                 	|
| GO.db                  3.12.1      2020-12-18 [1] Bioconductor                                   	|
| GOSemSim               2.16.1      2020-10-29 [1] Bioconductor                                   	|
| graph                  1.68.0      2020-10-27 [1] Bioconductor                                   	|
| graphite               1.36.0      2020-10-27 [1] Bioconductor                                   	|
| graphlayouts           0.7.1       2020-10-26 [1] CRAN (R 4.0.2)                                 	|
| gridExtra            * 2.3         2017-09-09 [1] CRAN (R 4.0.2)                                 	|
| gtable                 0.3.0       2019-03-25 [1] CRAN (R 4.0.2)                                 	|
| haven                  2.4.1       2021-04-23 [1] CRAN (R 4.0.2)                                 	|
| hms                    1.1.0       2021-05-17 [1] CRAN (R 4.0.2)                                 	|
| htmltools              0.5.1.1     2021-01-22 [1] CRAN (R 4.0.2)                                 	|
| httr                   1.4.2       2020-07-20 [1] CRAN (R 4.0.2)                                 	|
| igraph                 1.2.6       2020-10-06 [1] CRAN (R 4.0.2)                                 	|
| IRanges              * 2.24.1      2020-12-12 [1] Bioconductor                                   	|
| KEGGgraph              1.50.0      2020-10-27 [1] Bioconductor                                   	|
| KEGGREST               1.30.1      2020-11-23 [1] Bioconductor                                   	|
| KernSmooth             2.23-17     2020-04-26 [2] CRAN (R 4.0.3)                                 	|
| knitr                  1.33        2021-04-24 [1] CRAN (R 4.0.2)                                 	|
| labeling               0.4.2       2020-10-20 [1] CRAN (R 4.0.2)                                 	|
| lattice                0.20-41     2020-04-02 [2] CRAN (R 4.0.3)                                 	|
| lifecycle              1.0.0       2021-02-15 [1] CRAN (R 4.0.2)                                 	|
| limma                  3.46.0      2020-10-27 [1] Bioconductor                                   	|
| locfit                 1.5-9.4     2020-03-25 [1] CRAN (R 4.0.2)                                 	|
| magrittr               2.0.1       2020-11-17 [1] CRAN (R 4.0.2)                                 	|
| MASS                   7.3-53      2020-09-09 [2] CRAN (R 4.0.3)                                 	|
| Matrix                 1.2-18      2019-11-27 [2] CRAN (R 4.0.3)                                 	|
| MatrixGenerics       * 1.2.1       2021-01-30 [1] Bioconductor                                   	|
| matrixStats          * 0.58.0      2021-01-29 [1] CRAN (R 4.0.2)                                 	|
| memoise                2.0.0       2021-01-26 [1] CRAN (R 4.0.2)                                 	|
| munsell                0.5.0       2018-06-12 [1] CRAN (R 4.0.2)                                 	|
| openssl                1.4.4       2021-04-30 [1] CRAN (R 4.0.2)                                 	|
| openxlsx             * 4.2.3       2020-10-27 [1] CRAN (R 4.0.2)                                 	|
| org.Hs.eg.db           3.12.0      2020-12-18 [1] Bioconductor                                   	|
| org.Mm.eg.db         * 3.12.0      2021-04-19 [1] Bioconductor                                   	|
| pacman                 0.5.1       2019-03-11 [1] CRAN (R 4.0.2)                                 	|
| pathview             * 1.30.1      2020-12-10 [1] Bioconductor                                   	|
| pheatmap             * 1.0.12      2019-01-04 [1] CRAN (R 4.0.2)                                 	|
| pillar                 1.6.1       2021-05-16 [1] CRAN (R 4.0.2)                                 	|
| pkgconfig              2.0.3       2019-09-22 [1] CRAN (R 4.0.2)                                 	|
| plyr                   1.8.6       2020-03-03 [1] CRAN (R 4.0.2)                                 	|
| png                    0.1-7       2013-12-03 [1] CRAN (R 4.0.2)                                 	|
| polyclip               1.10-0      2019-03-14 [1] CRAN (R 4.0.2)                                 	|
| preprocessCore         1.52.1      2021-01-08 [1] Bioconductor                                   	|
| prettyunits            1.1.1       2020-01-24 [1] CRAN (R 4.0.2)                                 	|
| progress               1.2.2       2019-05-16 [1] CRAN (R 4.0.2)                                 	|
| purrr                  0.3.4       2020-04-17 [1] CRAN (R 4.0.2)                                 	|
| qvalue                 2.22.0      2020-10-27 [1] Bioconductor                                   	|
| R6                     2.5.0       2020-10-28 [1] CRAN (R 4.0.2)                                 	|
| rappdirs               0.3.3       2021-01-31 [1] CRAN (R 4.0.2)                                 	|
| RColorBrewer         * 1.1-2       2014-12-07 [1] CRAN (R 4.0.2)                                 	|
| Rcpp                   1.0.6       2021-01-15 [1] CRAN (R 4.0.2)                                 	|
| RCurl                  1.98-1.3    2021-03-16 [1] CRAN (R 4.0.2)                                 	|
| reactome.db            1.74.0      2020-12-18 [1] Bioconductor                                   	|
| ReactomePA           * 1.34.0      2020-10-27 [1] Bioconductor                                   	|
| readxl                 1.3.1       2019-03-13 [1] CRAN (R 4.0.2)                                 	|
| reshape2               1.4.4       2020-04-09 [1] CRAN (R 4.0.2)                                 	|
| Rgraphviz              2.34.0      2020-10-27 [1] Bioconductor                                   	|
| rio                    0.5.26      2021-03-01 [1] CRAN (R 4.0.2)                                 	|
| rjson                  0.2.20      2018-06-08 [1] CRAN (R 4.0.2)                                 	|
| rlang                  0.4.11      2021-04-30 [1] CRAN (R 4.0.2)                                 	|
| rmarkdown              2.8         2021-05-07 [1] CRAN (R 4.0.2)                                 	|
| RSQLite                2.2.7       2021-04-22 [1] CRAN (R 4.0.2)                                 	|
| rstatix                0.7.0       2021-02-13 [1] CRAN (R 4.0.2)                                 	|
| rstudioapi             0.13        2020-11-12 [1] CRAN (R 4.0.2)                                 	|
| rvcheck                0.1.8       2020-03-01 [1] CRAN (R 4.0.2)                                 	|
| S4Vectors            * 0.28.1      2020-12-09 [1] Bioconductor                                   	|
| scales                 1.1.1       2020-05-11 [1] CRAN (R 4.0.2)                                 	|
| scatterpie             0.1.6       2021-04-23 [1] CRAN (R 4.0.2)                                 	|
| scatterplot3d        * 0.3-41      2018-03-14 [1] CRAN (R 4.0.2)                                 	|
| sessioninfo            1.1.1       2018-11-05 [1] CRAN (R 4.0.2)                                 	|
| shadowtext             0.0.8       2021-04-23 [1] CRAN (R 4.0.2)                                 	|
| shape                  1.4.6       2021-05-19 [1] CRAN (R 4.0.3)                                 	|
| sp                     1.4-5       2021-01-10 [1] CRAN (R 4.0.2)                                 	|
| stringi                1.6.2       2021-05-17 [1] CRAN (R 4.0.2)                                 	|
| stringr              * 1.4.0       2019-02-10 [1] CRAN (R 4.0.2)                                 	|
| SummarizedExperiment * 1.20.0      2020-10-27 [1] Bioconductor                                   	|
| survival               3.2-7       2020-09-28 [2] CRAN (R 4.0.3)                                 	|
| tibble                 3.1.2       2021-05-16 [1] CRAN (R 4.0.2)                                 	|
| tidygraph              1.2.0       2020-05-12 [1] CRAN (R 4.0.2)                                 	|
| tidyr                * 1.1.3       2021-03-03 [1] CRAN (R 4.0.2)                                 	|
| tidyselect             1.1.1       2021-04-30 [1] CRAN (R 4.0.2)                                 	|
| tinytex                0.31        2021-03-30 [1] CRAN (R 4.0.2)                                 	|
| tweenr                 1.0.2       2021-03-23 [1] CRAN (R 4.0.2)                                 	|
| utf8                   1.2.1       2021-03-12 [1] CRAN (R 4.0.2)                                 	|
| vctrs                  0.3.8       2021-04-29 [1] CRAN (R 4.0.2)                                 	|
| venndir              * 0.0.15.9000 2021-06-17 [1] Github (jmw86069/venndir@38b3865)              	|
| vipor                  0.4.5       2017-03-22 [1] CRAN (R 4.0.2)                                 	|
| viridis                0.6.1       2021-05-11 [1] CRAN (R 4.0.2)                                 	|
| viridisLite            0.4.0       2021-04-13 [1] CRAN (R 4.0.2)                                 	|
| vsn                  * 3.58.0      2020-10-28 [1] Bioconductor                                   	|
| withr                  2.4.2       2021-04-18 [1] CRAN (R 4.0.2)                                 	|
| xfun                   0.23        2021-05-15 [1] CRAN (R 4.0.2)                                 	|
| XML                    3.99-0.6    2021-03-16 [1] CRAN (R 4.0.2)                                 	|
| xml2                   1.3.2       2020-04-23 [1] CRAN (R 4.0.2)                                 	|
| xtable               * 1.8-4       2019-04-21 [1] CRAN (R 4.0.2)                                 	|
| XVector                0.30.0      2020-10-28 [1] Bioconductor                                   	|
| yaml                   2.2.1       2020-02-01 [1] CRAN (R 4.0.2)                                 	|
| zip                    2.1.1       2020-08-27 [1] CRAN (R 4.0.2)                                 	|
| zlibbioc               1.36.0      2020-10-28 [1] Bioconductor                                   	|
|                                                                                                  	|
|                                                                                                  	|
| [1] /Users/sayanbakshi/Library/R/4.0/library                                                     	|
| [2] /Library/Frameworks/R.framework/Versions/4.0/Resources/library                               	|
