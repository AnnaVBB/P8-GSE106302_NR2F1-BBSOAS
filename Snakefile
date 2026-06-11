#########################################
# CONFIGURAÇÃO 
#########################################
configfile: "config.yaml"

import pandas as pd

samples = pd.read_csv(
    "samples.tsv",
    sep="\t"
)

SAMPLES = samples["sample_id"].tolist()
#########################################
# P8-GSE106302 - RNA-seq Nr2f1/BBSOAS
#########################################
rule all:
    input:
        expand(
            "results/fastqc/raw/{sample}_1_fastqc.html",
            sample=SAMPLES
        ),
        expand(
            "results/fastqc/raw/{sample}_2_fastqc.html",
            sample=SAMPLES
        ),
        "results/enrichment/GO_results.csv",
        "results/enrichment/KEGG_results.csv"

#########################################
# REFERÊNCIA
#########################################
rule salmon_index:
    input:
        "reference/transcriptome.fa"
    output:
        directory("reference/salmon_index")
    conda:
        "envs/salmon.yaml"
    threads:
        config["resources"]["salmon"]["threads"]
    shell:
        """
        salmon index \
            -t {input} \
            -i {output} \
            -p {threads}
        """

#########################################
# DADOS BRUTOS
#########################################
rule download_fastq:
    output:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"

    conda:
        "envs/sra.yaml"

    shell:
        """
        fasterq-dump {wildcards.sample} \
            -O data/raw \
            --split-files

        gzip data/raw/{wildcards.sample}_1.fastq
        gzip data/raw/{wildcards.sample}_2.fastq
        """
#########################################
# QC
#########################################
rule fastqc_raw:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"

    output:
        html1="results/fastqc/raw/{sample}_1_fastqc.html",
        zip1="results/fastqc/raw/{sample}_1_fastqc.zip",
        html2="results/fastqc/raw/{sample}_2_fastqc.html",
        zip2="results/fastqc/raw/{sample}_2_fastqc.zip"

    conda:
        "envs/fastqc.yaml"

    shell:
        """
        fastqc \
            {input.r1} \
            {input.r2} \
            -o results/fastqc/raw
        """

rule fastp:
    output:
        "results/fastqc/trim_complete.txt"

rule fastqc_trimmed:
    output:
        "results/fastqc/post_trim_complete.txt"

rule multiqc:
    output:
        "results/fastqc/multiqc_complete.txt"

#########################################
# QUANTIFICAÇÃO
#########################################
rule salmon_quant:
    output:
        "results/salmon/quant_complete.txt"

#########################################
# ANÁLISE ESTATÍSTICA
#########################################
rule tximport:
    output:
        "results/deseq2/tximport_complete.txt"

rule deseq2:
    output:
        "results/deseq2/DEG_results.csv"

#########################################
# ENRIQUECIMENTO
#########################################
rule enrichment:
    output:
        "results/enrichment/GO_results.csv",
        "results/enrichment/KEGG_results.csv"
