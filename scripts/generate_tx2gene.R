library(GenomicFeatures)

# Ajustado para a nova versão exigida pelo Bioconductor
txdb <- txdbmaker::makeTxDbFromGFF(
    "reference/annotation.gtf",
    format="gtf"
)

k <- keys(txdb, keytype="TXNAME")

tx2gene <- select(
    txdb,
    keys=k,
    columns="GENEID",
    keytype="TXNAME"
)

write.csv(
    tx2gene,
    "reference/tx2gene.csv",
    row.names=FALSE
)
