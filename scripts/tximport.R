library(tximport)
library(readr)

# 1. Ler o arquivo de metadados das amostras
samples <- read.table("samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 2. Mapear os caminhos dos arquivos quant.sf do Salmon
files <- file.path("results/salmon", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# 3. Ler a tabela de conversão transcrito-para-gene
tx2gene <- read.csv("reference/tx2gene.csv")

# 4. Executar a consolidação dos dados ao nível de gene
txi <- tximport(
    files, 
    type = "salmon", 
    tx2gene = tx2gene, 
    ignoreTxVersion = TRUE
)

# 5. Criar a pasta de destino caso não exista e salvar o objeto R estruturado
dir.create("results/deseq2", showWarnings = FALSE, recursive = TRUE)
saveRDS(txi, "results/deseq2/txi.rds")

# 6. Gravar o arquivo de texto final exigido pelo Snakemake como prova de sucesso
writeLines("Tximport concluído com sucesso. Dados consolidados ao nível de gene.", "results/deseq2/tximport_complete.txt")
