#########################################
# CONFIGURAÇÃO (URLs do genoma)
#########################################
config = {
    "fasta": "https://ftp.ensembl.org/pub/current_fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
    "gtf": "https://ensembl.org"
}

#########################################
# P8-GSE106302 - RNA-seq Nr2f1/BBSOAS
#########################################

rule all:
    input:
        "results/enrichment/GO_results.csv",
        "results/enrichment/KEGG_results.csv"

#########################################
# REFERÊNCIA
#########################################

rule download_transcriptome:
    output:
        "reference/transcriptome.fa"

rule download_gtf:
    output:
        "reference/annotation.gtf"

rule generate_tx2gene:
    input:
        "reference/annotation.gtf"
    output:
        "reference/tx2gene.csv"

rule salmon_index:
    input:
        "reference/transcriptome.fa"
    output:
        directory("reference/salmon_index")

#########################################
# DADOS BRUTOS
#########################################

rule download_fastq:
    output:
        "data/raw/download_complete.txt"

#########################################
# QC
#########################################

rule fastqc_raw:
    input:
        "data/raw/download_complete.txt"
    output:
        "results/fastqc/raw_complete.txt"

rule fastp:
    input:
        "results/fastqc/raw_complete.txt"
    output:
        "results/fastqc/trim_complete.txt"

rule fastqc_trimmed:
    input:
        "results/fastqc/trim_complete.txt"
    output:
        "results/fastqc/post_trim_complete.txt"

rule multiqc:
    input:
        "results/fastqc/post_trim_complete.txt"
    output:
        "results/fastqc/multiqc_complete.txt"

#########################################
# QUANTIFICAÇÃO
#########################################

rule salmon_quant:
    input:
        "results/fastqc/multiqc_complete.txt",
        "reference/salmon_index"
    output:
        "results/salmon/quant_complete.txt"

#########################################
# ANÁLISE ESTATÍSTICA
#########################################

rule tximport:
    input:
        "results/salmon/quant_complete.txt"
    output:
        "results/deseq2/tximport_complete.txt"

rule deseq2:
    input:
        "results/deseq2/tximport_complete.txt"
    output:
        "results/deseq2/DEG_results.csv"

#########################################
# ENRIQUECIMENTO
#########################################

rule enrichment:
    input:
        "results/deseq2/DEG_results.csv"
    output:
        "results/enrichment/GO_results.csv",
        "results/enrichment/KEGG_results.csv"
