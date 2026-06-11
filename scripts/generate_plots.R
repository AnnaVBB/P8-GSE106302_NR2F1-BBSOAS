library(DESeq2)
library(ggplot2)
library(pheatmap)

# 1. Carregar objetos
deg <- read.csv("results/deseq2/DEG_results.csv", row.names = 1)
dds <- readRDS("results/deseq2/dds_object.rds")
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# 2. Gráfico PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="genotype", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
  geom_point(size=3) + xlab(paste0("PC1: ",percentVar,"% variance")) +
  ylab(paste0("PC2: ",percentVar,"% variance")) + theme_minimal()
ggsave("results/plots/PCA_plot.png", width=6, height=5)

# 3. Volcano Plot
png("results/plots/volcano_plot.png", width=800, height=600)
plot(deg$log2FoldChange, -log10(deg$pvalue), pch=20, col="grey", main="Volcano Plot Nr2f1")
points(deg$log2FoldChange[deg$padj<0.05 & abs(deg$log2FoldChange)>1], 
       -log10(deg$pvalue)[deg$padj<0.05 & abs(deg$log2FoldChange)>1], pch=20, col="red")
dev.off()

# 4. Heatmap de DEGs
sig_genes <- rownames(deg[which(deg$padj < 0.05)[1:50], ]) # Top 50 genes
mat <- assay(vsd)[sig_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=as.data.frame(colData(dds)["genotype"]), 
         filename="results/plots/heatmap_DEG.png", width=6, height=8)

# 5. MA Plot
png("results/plots/MA_plot.png", width=800, height=600)
plotMA(results(dds), ylim=c(-5,5), main="MA Plot Nr2f1")
dev.off()

# 6. Criar Relatório HTML de controle básico
writeLines("<html><body><h1>Relatório de RNA-seq Nr2f1</h1><p>Análise concluída. Verifique os plots salvos na pasta results/plots/.</p></body></html>", "results/plots/RELATORIO_FINAL.html")

