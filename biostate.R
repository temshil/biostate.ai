library(tximport)
library(edgeR)
library(sva)
library(limma)
library(DESeq2)
library(RUVSeq)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(rtracklayer)

gtf0 <- import("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.primary_assembly.annotation.gtf.gz")
gtf0 <- gtf0[gtf0$type=="exon"]
gtf1 <- gtf0[!duplicated(gtf0$transcript_id)]
gtf1b <- gtf1[!is.na(gtf1$transcript_id)]
all_tx_info <- mcols(gtf1b)[,c("transcript_id", "gene_id", "gene_name")]
save(all_tx_info, file="all_tx_info.RData")
#load("all_tx_info.RData")

samples <- read.csv("samples.csv")
samples$genotype <- factor(samples$genotype)
samples$timepoint <- factor(samples$timepoint)
samples$batch <- factor(samples$batch)
samples$condition <- factor(samples$condition)
samples$sample_id <- factor(samples$sample_id)

files <- file.path("data", samples$folder, "quant.sf")
names(files) <- samples$sample_id
txi <- tximport(files, type = "salmon", txOut = TRUE) 
txi <- summarizeToGene(txi, all_tx_info[,1:2], countsFromAbundance = "lengthScaledTPM")

genes <- data.frame(gene_id=rownames(txi$counts))
gene_names <- unique(all_tx_info[,c(2,3)])
genes_merge <- merge(genes, gene_names, by.x="gene_id", by.y="gene_id")

y <- DGEList(txi$counts, samples = samples, genes=genes_merge)

plotMDS(y)

genotype <- relevel(factor(samples$genotype), ref = "wt")
timepoint <-  factor(samples$timepoint)
batch <-factor(samples$batch)

isexpr <- filterByExpr(y, group = interaction(genotype, timepoint))
table(isexpr)

hasannot <- rowSums(is.na(y$genes))==0
y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]
dim(y)

y <- normLibSizes(y)

design <- model.matrix(~batch + genotype * timepoint)

y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

colnames(design) <- make.names(colnames(design))
colnames(design)

cont.matrix <- makeContrasts(
  DelvsWT_0  = genotypedel,
  DelvsWT_7  = genotypedel + genotypedel.timepoint7,
  DelvsWT_14 = genotypedel + genotypedel.timepoint14,
  DelvsWT_23 = genotypedel + genotypedel.timepoint23,
  levels = design
)

for (i in 1:4) {
  qlf <- glmQLFTest(fit,contrast=cont.matrix[, i])
  res <- as.data.frame(topTags(qlf, n=Inf))
  sig <- res[res$FDR<0.05,]
  print(paste0(i, " ",nrow(sig)))
}

# [1] "1 42"
# [1] "2 6"
# [1] "3 9"
# [1] "4 0"

#sva
design0 <- model.matrix(~1, data=samples)

svobj <- sva(y$counts, design, design0)

design.sva <- cbind(design, svobj$sv)
colnames(design.sva) <- c(colnames(design), paste0("SV", 1:ncol(svobj$sv)))
colnames(design.sva) <- make.names(colnames(design.sva))
colnames(design.sva)

y.sva <- estimateDisp(y, design.sva, robust=TRUE)
plotBCV(y.sva)

fit.sva <- glmQLFit(y.sva, design.sva, robust=TRUE)
plotQLDisp(fit.sva)

nSV <- ncol(svobj$sv)

cont.matrix.sva <- rbind(cont.matrix, matrix(0, nrow=nSV, ncol=ncol(cont.matrix)))

for (i in 1:4) {
  qlf.sva <- glmQLFTest(fit.sva,contrast=cont.matrix.sva[, i])
  res.sva <- as.data.frame(topTags(qlf.sva, n=Inf))
  print(paste0(i, " ",res.sva[res.sva$gene_name=="Mov10",]))
  sig.sva <- res.sva[res.sva$FDR<0.05,]
  print(paste0(i, " ",nrow(sig.sva)))
}

# [1] "1 32"
# [1] "2 3"
# [1] "3 12"
# [1] "4 1"

#sva ComBat
counts <- txi$counts      
batch <- samples$batch     
condition <-  factor(samples$condition)

counts_combat <- ComBat_seq(counts = counts,
                            batch = batch,
                            group = condition)

ycombat <- DGEList(counts_combat, samples = samples, genes=genes_merge)

plotMDS(ycombat)

isexprcombat <- filterByExpr(ycombat, group=condition)
table(isexprcombat)

hasannotcombat <- rowSums(is.na(ycombat$genes))==0
ycombat <- ycombat[isexprcombat & hasannotcombat, , keep.lib.sizes=FALSE]
dim(ycombat)

ycombat <- normLibSizes(ycombat)

designcombat <- model.matrix(~condition)

colnames(designcombat) <- levels(condition)
ycombat <- estimateDisp(ycombat, designcombat, robust=TRUE)
plotBCV(ycombat)

fitcombat <- glmQLFit(ycombat, designcombat, robust=TRUE)
plotQLDisp(fitcombat)

cont.matrixcombat <- makeContrasts(DelvsWT.0 = D0-W0,
                             DelvsWT.7 = D7-W7,
                             DelvsWT.14 = D14-W14,
                             DelvsWT.23 = D23-W23,
                             levels = designcombat)

for (i in 1:4) {
  qlf <- glmQLFTest(fitcombat,contrast=cont.matrixcombat[, i])
  res <- as.data.frame(topTags(qlf, n=Inf))
  print(paste0(i, " ",res[res$gene_name=="Mov10",]))
  sig <- res[res$FDR<0.05,]
  print(paste0(i, " ",nrow(sig)))
}

# [1] "1 18616"
# [1] "2 41"
# [1] "3 13"
# [1] "4 10"

#DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ batch + genotype + timepoint + genotype:timepoint)

dds$genotype <- relevel(dds$genotype, ref = "wt")      

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup=c("genotype","timepoint"))

# resultsNames(dds)

res0 <- na.omit(results(dds, name = "genotype_del_vs_wt"))
sig0 <- res0[res0$padj<0.05,]
nrow(sig0)
#[1] 40

res7 <- na.omit(results(dds, contrast=list("genotype_del_vs_wt","genotypedel.timepoint7")))
sig7 <- res7[res7$padj<0.05,]
nrow(sig7)
#[1] 38

res14 <- na.omit(results(dds, contrast=list("genotype_del_vs_wt","genotypedel.timepoint14")))
sig14 <- res14[res14$padj<0.05,]
nrow(sig14)
#[1] 31

res23 <- na.omit(results(dds, contrast=list("genotype_del_vs_wt","genotypedel.timepoint23")))
sig23 <- res23[res23$padj<0.05,]
nrow(sig23)
#[1] 36

#RUVseq
set <- newSeqExpressionSet(counts(dds))
res <- as.data.frame(results(dds))
diff <- makeGroups(samples$condition)
empirical <- rownames(tail(arrange(res, padj),n = 5000))

set <- RUVs(set, empirical, k = 2, diff)
pData(set)

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds$sample_id, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

samples$W_1 <- set$W_1
samples$W_2 <- set$W_2

ddsruv <- DESeqDataSetFromTximport(txi, colData = samples, design = ~W_1 + genotype + timepoint + genotype:timepoint)
ddsruv$genotype <- relevel(dds$genotype, ref = "wt")   

smallestGroupSize <- 3
keepruv <- rowSums(counts(ddsruv) >= 10) >= smallestGroupSize
ddsruv <- ddsruv[keepruv,]

ddsruv <- DESeq(ddsruv)
vsdruv <- vst(ddsruv, blind = TRUE)
plotPCA(vsdruv, intgroup=c("genotype","timepoint"))

res0ruv <- na.omit(results(ddsruv, name = "genotype_del_vs_wt"))
sig0ruv <- res0ruv[res0ruv$padj<0.05,]
nrow(sig0ruv)
#[1] 31
nrow(sig0)

res7ruv <- na.omit(results(ddsruv, contrast=list("genotype_del_vs_wt","genotypedel.timepoint7")))
sig7ruv <- res7ruv[res7ruv$padj<0.05,]
nrow(sig7ruv)
#[1] 42
nrow(sig7)

res14ruv <- na.omit(results(ddsruv, contrast=list("genotype_del_vs_wt","genotypedel.timepoint14")))
sig14ruv <- res14ruv[res14ruv$padj<0.05,]
nrow(sig14ruv)
#[1] 29
nrow(sig14)

res23ruv <- na.omit(results(ddsruv, contrast=list("genotype_del_vs_wt","genotypedel.timepoint23")))
sig23ruv <- res23ruv[res23ruv$padj<0.05,]
nrow(sig23ruv)
#[1] 43
nrow(sig23)

samples_del <- samples[c(25,26,27,28,29,30),]

files_del <- file.path("data", samples_del$folder, "quant.sf")
txi_del <- tximport(files_del, type = "salmon", txOut = TRUE) 
txi_del <- summarizeToGene(txi_del, all_tx_info[,1:2], countsFromAbundance = "lengthScaledTPM")
dds_del <- DESeqDataSetFromTximport(txi_del, colData = samples_del, design = ~ condition)
smallestGroupSize <- 3
keep_del <- rowSums(counts(dds_del) >= 10) >= smallestGroupSize
dds_del <- dds_del[keep_del,]
dds_del <- DESeq(dds_del)
res_del <- results(dds_del,contrast = c("condition","D0","W0"))

up_del <- res_del[!is.na(res_del$padj)&res_del$padj<0.05&res_del$log2FoldChange>0,]
down_del <- res_del[!is.na(res_del$padj)&res_del$padj<0.05&res_del$log2FoldChange<0,]
all_tx_info[all_tx_info$gene_name=="Numa1",]
'ENSMUSG00000066306.14'%in%rownames(down_del)

up0 <- sig0[sig0$log2FoldChange>0,]
down0 <- sig0[sig0$log2FoldChange<0,]

up0ruv <- sig0ruv[sig0ruv$log2FoldChange>0,]
down0ruv <- sig0ruv[sig0ruv$log2FoldChange<0,]

up_del_orig <- read.csv("up.csv", header = TRUE)
down_del_orig <- read.csv("down.csv", header = TRUE)

up_all <- data.frame(gene=unique(c(rownames(up_del),rownames(up0),rownames(up0ruv),up_del_orig$gene_id)))
up_all$up_del_orig <- up_all$gene%in%up_del_orig$gene_id
up_all$up_del <- up_all$gene%in%rownames(up_del)
up_all$up0 <- up_all$gene%in%rownames(up0)
up_all$up0ruv <-  up_all$gene%in%rownames(up0ruv)
vennDiagram(up_all[,c(2,3,4,5)])

down_all <- data.frame(gene=unique(c(rownames(down_del),rownames(down0),rownames(down0ruv),down_del_orig$gene_id)))
down_all$down_del_orig <- down_all$gene%in%down_del_orig$gene_id
down_all$down_del <- down_all$gene%in%rownames(down_del)
down_all$down0 <- down_all$gene%in%rownames(down0)
down_all$down0ruv <-  down_all$gene%in%rownames(down0ruv)
vennDiagram(down_all[,c(2,3,4,5)])

#Exporting DEGs
gene_map <- as.data.frame(all_tx_info[, c("gene_id", "gene_name")])
gene_map_unique <- gene_map[!duplicated(gene_map$gene_id), ]

gene_map_unique$gene_id <- as.character(gene_map_unique$gene_id)
gene_map_unique$gene_name <- as.character(gene_map_unique$gene_name)

add_gene_name <- function(res, gene_map) {
  res$gene_id <- rownames(res)
  res <- merge(res, gene_map, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
  rownames(res) <- res$gene_id
  return(res)
}

sig0 <- add_gene_name(as.data.frame(sig0), gene_map_unique)
sig7 <- add_gene_name(as.data.frame(sig7), gene_map_unique)
sig14 <- add_gene_name(as.data.frame(sig14), gene_map_unique)
sig23 <- add_gene_name(as.data.frame(sig23), gene_map_unique)

sig0[sig0$gene_name=="Mov10",]
sig7[sig7$gene_name=="Mov10",]
sig14[sig14$gene_name=="Mov10",]
sig23[sig23$gene_name=="Mov10",]

res0_new <- add_gene_name(as.data.frame(res0), gene_map_unique)
res7_new <- add_gene_name(as.data.frame(res7), gene_map_unique)
res14_new <- add_gene_name(as.data.frame(res14), gene_map_unique)
res23_new <- add_gene_name(as.data.frame(res23), gene_map_unique)

res0_new[res0_new$gene_name=="Mov10",]
res7_new[res7_new$gene_name=="Mov10",]
res14_new[res14_new$gene_name=="Mov10",]
res23_new[res23_new$gene_name=="Mov10",]

all_degs <- rbind(sig0,sig7,sig14,sig23)

all_degs <- all_degs[!duplicated(all_degs$gene_id), ]
deg_genes <- all_degs$gene_id

all_degs$sig0 <- rownames(all_degs)%in%rownames(sig0)
all_degs$sig7 <- rownames(all_degs)%in%rownames(sig7)
all_degs$sig14 <- rownames(all_degs)%in%rownames(sig14)
all_degs$sig23 <- rownames(all_degs)%in%rownames(sig23)
vennDiagram(all_degs[,c(9,10,11,12)])

all_degs[all_degs$sig0==TRUE&all_degs$sig7==TRUE&all_degs$sig14==TRUE&all_degs$sig23==TRUE,]

write.csv(sig0, file="sig0.csv")
write.csv(sig7, file="sig7.csv")
write.csv(sig14, file="sig14.csv")
write.csv(sig23, file="sig23.csv")

panther_sig0 <- read.table("pantherGeneList_sig0.txt", sep="\t", header = FALSE)
colnames(panther_sig0) <- c('Gene_ID','gene_name','GeneName','PANTHER_Family_Subfamily','PANTHER_Protein_Class')
sig0 <- merge(sig0, panther_sig0, by = "gene_name", all = TRUE)

panther_sig7 <- read.table("pantherGeneList_sig7.txt", sep="\t", header = FALSE)
colnames(panther_sig7) <- c('Gene_ID','gene_name','GeneName','PANTHER_Family_Subfamily','PANTHER_Protein_Class')
sig7 <- merge(sig7, panther_sig7, by = "gene_name", all = TRUE)

panther_sig14 <- read.table("pantherGeneList_sig14.txt", sep="\t", header = FALSE)
colnames(panther_sig14) <- c('Gene_ID','gene_name','GeneName','PANTHER_Family_Subfamily','PANTHER_Protein_Class')
sig14 <- merge(sig14, panther_sig14, by = "gene_name", all = TRUE)

panther_sig23 <- read.table("pantherGeneList_sig23.txt", sep="\t", header = FALSE)
colnames(panther_sig23) <- c('Gene_ID','gene_name','GeneName','PANTHER_Family_Subfamily','PANTHER_Protein_Class')
sig23 <- merge(sig23, panther_sig23, by = "gene_name", all = TRUE)

write.csv(sig0, file="sig0.csv")
write.csv(sig7, file="sig7.csv")
write.csv(sig14, file="sig14.csv")
write.csv(sig23, file="sig23.csv")

#Clustering
expr_mat <- assay(vsd)

colnames(expr_mat) <- samples$sample_id

expr_mat_scaled <- t(scale(t(expr_mat)))

expr_mat_scaled_clean <- expr_mat_scaled[complete.cases(expr_mat_scaled), ]
#deg_genes <- unique(rownames(all_degs))
#expr_mat_scaled_clean <- expr_mat_scaled_clean[deg_genes,]

dist_mat <- dist(expr_mat_scaled_clean)
clust_res <- hclust(dist_mat, method = "ward.D2")
plot(clust_res)
clusters <- cutree(clust_res, k = 5)

#save(dist_mat, clust_res, clusters, file="clusters.RData")

ordered_genes <- names(sort(clusters))
expr_mat_ordered <- expr_mat_scaled_clean[ordered_genes,c(25,26,27,1,2,3,4,5,6,7,8,9,10,11,12,28,29,30,13,14,15,16,17,18,19,20,21,22,23,24)]

annotation_col <- data.frame(
  timepoint = factor(as.character(samples$timepoint))
)
rownames(annotation_col) <- samples$sample_id 

annotation_row <- data.frame(
  cluster = as.factor(clusters[rownames(expr_mat_ordered)])
)
rownames(annotation_row) <- rownames(expr_mat_ordered)

ann_colors <- list(
  timepoint = c("0" = "#a6cee3", "7" = "#1f78b4", "14" = "#b2df8a", "23" = "#33a02c")
)


breaks <- seq(-3, 3, length.out = 100)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaks) - 1)

img <- pheatmap(expr_mat_ordered,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_col = annotation_col,
                annotation_row = annotation_row,
                annotation_colors = ann_colors,
                show_rownames = FALSE,
                color = colors,
                breaks = breaks,
                main = "All genes")

save_pheatmap_png <- function(x, filename, width=8, height=6.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height, res=600, units="in")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(img, "biostate_all_genes.png")

genes_by_cluster <- split(names(clusters), clusters)
expr_cluster1 <- data.frame(gene_id=rownames(expr_mat_scaled_clean[genes_by_cluster$`1`, ]))
expr_cluster2 <- data.frame(gene_id=rownames(expr_mat_scaled_clean[genes_by_cluster$`2`, ]))
expr_cluster3 <- data.frame(gene_id=rownames(expr_mat_scaled_clean[genes_by_cluster$`3`, ]))
expr_cluster4 <- data.frame(gene_id=rownames(expr_mat_scaled_clean[genes_by_cluster$`4`, ]))
expr_cluster5 <- data.frame(gene_id=rownames(expr_mat_scaled_clean[genes_by_cluster$`5`, ]))

write.csv(expr_cluster1, file="expr_cluster1.csv", row.names = FALSE, quote = FALSE)
write.csv(expr_cluster2, file="expr_cluster2.csv", row.names = FALSE, quote = FALSE)
write.csv(expr_cluster3, file="expr_cluster3.csv", row.names = FALSE, quote = FALSE)
write.csv(expr_cluster4, file="expr_cluster4.csv", row.names = FALSE, quote = FALSE)
write.csv(expr_cluster5, file="expr_cluster5.csv", row.names = FALSE, quote = FALSE)

for (i in seq(1,5)){
  name = paste0("expr_cluster",i)
  data <- read.table(paste0(name,".txt"), sep="\t", skip = 11, header = TRUE)
  colnames(data)[6] <- "fold"
  data <- data[order(data$fold),]
  data <- tail(data[,c(1,6,8)], n = 10)
  colnames(data) <- c("mf","fold", "fdr")
  data$log10 <- -log10(data$fdr)
  data$mf <- substr(data$mf,1,nchar(data$mf)-13)
  data$cont <- c(1:nrow(data))

  
  
  ggplot(data, aes(y=log10, x=cont))+
    geom_bar(stat="identity", color="black", position=position_dodge())+
    geom_hline(yintercept=1.3, color = "red")+
    ggtitle(name)+
    ylab("Negative Log10 FDR")+
    xlab("Fold Enrichment")+
    xlim(0,10)+
    coord_flip()+
    scale_x_continuous(breaks = data$cont,
                       labels = data$fold,
                       sec.axis = dup_axis(breaks = data$cont,
                                           labels = data$mf,
                                           name=""),
                       expand = c(0,0))
  
    ggsave(paste0(name,".png"), 
         width = 300,
         height = 150,
         dpi = 300,
         units = c("mm"))
}

sink("session_info.txt")
sessionInfo()
sink()