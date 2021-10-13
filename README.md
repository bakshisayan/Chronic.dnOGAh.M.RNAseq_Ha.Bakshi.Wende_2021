## Chronic.dnOGAh.M.RNAseq_Ha.Bakshi.Wende_2021 ##

### Author: 
Sayan Bakshi, M.Sc. (Integrated) | PhD Trainee | sbakshi@uab.edu
### First Author of the manuscript: 
Dr. Chae-Myeong Ha, PhD | cha@uabmc.edu
### Corresponding Author / PI: 
Dr. Adam R. Wende, PhD, Department of Pathology, University of Alabama at Birmingham, Birmingham, Alabama, USA | adamwende@uabmc.edu

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
sessioninfo()
```
| R version 4.0.3 (2020-10-10)                                                                  |
|-----------------------------------------------------------------------------------------------|
| Platform: x86_64-apple-darwin17.0 (64-bit)                                                    |
| Running under: macOS Big Sur 10.16                                                            |
|                                                                                               |
|                                                                                               |
| Matrix products: default                                                                      |
| LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib           |
|                                                                                               |
|                                                                                               |
| locale:                                                                                       |
| [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8                             |
|                                                                                               |
|                                                                                               |
| attached base packages:                                                                       |
| [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base      |
|                                                                                               |
|                                                                                               |
| other attached packages:                                                                      |
| [1] FactoClass_1.2.7      |      scatterplot3d_0.3-41  |     xtable_1.8-4                      |
| [4] ade4_1.7-18           |      eulerr_6.1.0          |      venndir_0.0.15.9000               |
| [7] org.Mm.eg.db_3.12.0   |      AnnotationDbi_1.52.0  |      enrichplot_1.10.2                 |
| [10] circlize_0.4.12      |       ReactomePA_1.34.0    |       pathview_1.30.1                  |
| [13] fgsea_1.16.0         |       enrichR_3.0          |       clusterProfiler_3.18.1           |
| [16] vsn_3.58.0           |       DESeq2_1.30.1        |       SummarizedExperiment_1.20.0      |
| [19] Biobase_2.50.0       |       MatrixGenerics_1.2.1 |       matrixStats_0.58.0               |
| [22] GenomicRanges_1.42.0 |       GenomeInfoDb_1.26.7  |       IRanges_2.24.1                   |
| [25] S4Vectors_0.28.1     |       BiocGenerics_0.36.1  |       biomaRt_2.46.3                   |
| [28] ggnewscale_0.4.5     |       ggbeeswarm_0.6.0     |       gridExtra_2.3                    |
| [31] ggfortify_0.4.11     |       ggpubr_0.4.0         |       ggrepel_0.9.1                    |
| [34] RColorBrewer_1.1-2   |       pheatmap_1.0.12      |       ggplot2_3.3.3                    |
| [37] stringr_1.4.0        |       openxlsx_4.2.3       |       tidyr_1.1.3                      |
| [40] dplyr_1.0.6          |                                                                    |
|                                                                                               |
|                                                                                               |
| loaded via a namespace (and not attached):                                                    |
| [1] pacman_0.5.1         |  utf8_1.2.1          |   tidyselect_1.1.1   |    RSQLite_2.2.7        |
| [5] grid_4.0.3           |  BiocParallel_1.24.1 |   scatterpie_0.1.6   |    munsell_0.5.0        |
| [9] preprocessCore_1.52.1 | withr_2.4.2         |   colorspace_2.0-1   |    GOSemSim_2.16.1      |
| [13] knitr_1.33          |   rstudioapi_0.13    |    ggsignif_0.6.1    |     DOSE_3.16.0         |
| [17] labeling_0.4.2      |   KEGGgraph_1.50.0   |    GenomeInfoDbData_1.2.4 | polyclip_1.10-0     |
| [21] bit64_4.0.5         |   farver_2.1.0       |    downloader_0.4       |  vctrs_0.3.8         |
| [25] generics_0.1.0      |   xfun_0.23          |    BiocFileCache_1.14.0 |  R6_2.5.0            |
| [29] graphlayouts_0.7.1  |   locfit_1.5-9.4     |    bitops_1.0-7         |  cachem_1.0.5        |
| [33] DelayedArray_0.16.3 |   assertthat_0.2.1   |    scales_1.1.1         |  ggraph_2.0.5        |
| [37] beeswarm_0.3.1      |   gtable_0.3.0       |    affy_1.68.0          |  tidygraph_1.2.0     |
| [41] rlang_0.4.11        |   genefilter_1.72.1  |    GlobalOptions_0.1.2  |  splines_4.0.3       |
| [45] rstatix_0.7.0       |   broom_0.7.6        |    checkmate_2.0.0      |  BiocManager_1.30.15 |
| [49] yaml_2.2.1          |   reshape2_1.4.4     |    abind_1.4-5          |  backports_1.2.1     |
| [53] qvalue_2.22.0       |   tools_4.0.3        |    affyio_1.60.0        |  ellipsis_0.3.2      |
| [57] Rcpp_1.0.6          |   plyr_1.8.6         |    progress_1.2.2       |  zlibbioc_1.36.0     |
| [61] purrr_0.3.4         |   RCurl_1.98-1.3     |    prettyunits_1.1.1    |  openssl_1.4.4       |
| [65] viridis_0.6.1       |   cowplot_1.1.1      |    haven_2.4.1          |  tinytex_0.31        |
| [69] magrittr_2.0.1      |   data.table_1.14.0  |    DO.db_2.9            |  reactome.db_1.74.0  |
| [73] hms_1.1.0           |   evaluate_0.14      |    XML_3.99-0.6         |  rio_0.5.26          |
| [77] readxl_1.3.1        |   shape_1.4.6        |    compiler_4.0.3       |  tibble_3.1.2        |
| [81] KernSmooth_2.23-17  |   crayon_1.4.1       |    shadowtext_0.0.8     |  htmltools_0.5.1.1   |
| [85] geneplotter_1.68.0  |   DBI_1.1.1          |    tweenr_1.0.2         |  dbplyr_2.1.1        |
| [89] MASS_7.3-53         |   rappdirs_0.3.3     |    Matrix_1.2-18        |  car_3.0-10          |
| [93] cli_2.5.0           |   igraph_1.2.6       |    forcats_0.5.1        |  pkgconfig_2.0.3     |
| [97] rvcheck_0.1.8       |   sp_1.4-5           |    foreign_0.8-80       |  xml2_1.3.2          |
| [101] annotate_1.68.0    |    vipor_0.4.5       |     XVector_0.30.0      |   digest_0.6.27      |
| [105] graph_1.68.0       |    Biostrings_2.58.0 |     rmarkdown_2.8       |   cellranger_1.1.0   |
| [109] fastmatch_1.1-0    |    curl_4.3.1        |     graphite_1.36.0     |   rjson_0.2.20       |
| [113] lifecycle_1.0.0    |    carData_3.0-4     |     viridisLite_0.4.0   |   askpass_1.1        |
| [117] limma_3.46.0       |    fansi_0.4.2       |     pillar_1.6.1        |   lattice_0.20-41    |
| [121] KEGGREST_1.30.1    |    fastmap_1.1.0     |     httr_1.4.2          |   survival_3.2-7     |
| [125] GO.db_3.12.1       |    glue_1.4.2        |     zip_2.1.1           |   png_0.1-7          |
| [129] bit_4.0.4          |    Rgraphviz_2.34.0  |     ggforce_0.3.3       |   stringi_1.6.2      |
| [133] blob_1.2.1         |    org.Hs.eg.db_3.12.0 |    memoise_2.0.0      |                      |
