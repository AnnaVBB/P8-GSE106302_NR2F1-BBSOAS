library(clusterProfiler)
library(org.Mm.eg.db)


# 1. Ler a tabela de resultados do DESeq2
deg <- read.csv("results/deseq2/DEG_results.csv", row.names = 1)

# 2. Filtrar apenas os genes significativamente alterados (padj < 0.05)
sig_genes <- rownames(deg[which(deg$padj < 0.05), ])

# 3. Criar a pasta de resultados se não existir
dir.create("results/enrichment", showWarnings = FALSE, recursive = TRUE)

# 4. Enriquecimento de Ontologia Gênica (GO) - Processos Biológicos
go_res <- enrichGO(
    gene          = sig_genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
)
write.csv(as.data.frame(go_res), "results/enrichment/GO_results.csv", row.names = FALSE)

# 5. Enriquecimento de Vias KEGG
# O KEGG exige conversão de ENSEMBL para códigos ENTREZ ID
gene_conversion <- bitr(sig_genes, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)

kegg_res <- enrichKEGG(
    gene         = gene_conversion$ENTREZID,
    organism     = "mmu", # mmu significa Mus musculus
    pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_res), "results/enrichment/KEGG_results.csv", row.names = FALSE)

message("Enriquecimento GO e KEGG concluído com sucesso!")
