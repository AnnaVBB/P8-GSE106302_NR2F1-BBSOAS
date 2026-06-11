library(DESeq2)

# 1. Ler os metadados e estruturar os fatores biológicos
samples <- read.table("samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples$sample_id

# Garantir que a coluna genotype seja interpretada como Fator (Factor)
# E definir "WT" como o grupo de referência (controle)
samples$genotype <- factor(samples$genotype, levels = c("WT", "HET"))

# 2. Carregar o objeto estruturado pelo tximport
txi <- readRDS("results/deseq2/txi.rds")

# 3. Construir o objeto de análise do DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ genotype)

# 4. Rodar a análise estatística de expressão diferencial
dds <- DESeq(dds)

# 5. Extrair os resultados da comparação (HET vs WT)
res <- results(dds, contrast=c("genotype", "HET", "WT"))

# Ordenar os resultados pelos menores p-valores ajustados (mais significativos)
resOrdered <- res[order(res$padj), ]

# Converter em dataframe e salvar em formato CSV
resData <- as.data.frame(resOrdered)
write.csv(resData, "results/deseq2/DEG_results.csv", row.names = TRUE)
saveRDS(dds, "results/deseq2/dds_object.rds")

message("DESeq2 concluído com sucesso! Tabela DEG_results.csv gerada.")
