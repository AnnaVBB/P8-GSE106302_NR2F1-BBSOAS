library(snakemake)

#A função makeTxDbFromGFF pega o caminho GTF, extrai as informações (genes, transcritos, éxons,etc) e organiza em um objeto TxDB (transcript database)

txdb <- txdbmaker::makeTxDbFromGFF(
    snakemake@input[[1]],
    format="gtf"
)

#extrai todos os IDs de transcritos no banco. Para ver as outras chaves, use keytypes(txdb)
k <- keys(txdb, keytype="TXNAME")

#seleciona o gene correspondente a cada transcript_ID
tx2gene <- select(
    txdb,
    keys=k,
    columns="GENEID",
    keytype="TXNAME"
)

write.csv(
    tx2gene,
    snakemake@output[[1]],
    row.names=FALSE
)
