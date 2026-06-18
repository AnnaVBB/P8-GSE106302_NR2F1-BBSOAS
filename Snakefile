#########################################
# CONFIGURAÇÃO
#########################################
import pandas as pd
configfile: "config.yaml"

samples = pd.read_csv(
    "samples.tsv",
    sep="\t"
)

SAMPLES = samples["sample_id"].tolist()

FASTA_URL = config["reference"]["fasta"]
GTF_URL = config["reference"]["gtf"]

TRANSCRIPTOME = "reference/transcriptome.fa"
ANNOTATION = "reference/annotation.gtf"
TX2GENE = "reference/tx2gene.csv"
#########################################
# P8-GSE106302 - RNA-seq Nr2f1/BBSOAS
#########################################
rule all:
    input:
        "reference/salmon_index", 
        # Execução do QC inicial e final
        "results/fastqc/multiqc_report.html",
        "results/fastqc/trimmed/multiqc_trimmed_report.html",
        # Quantificação
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        # Análise estatística e tabelas de expressão
        "results/deseq2/tximport_complete.txt",
        "results/deseq2/DEG_results.csv",
        # Enriquecimento funcional
        "results/enrichment/GO_results.csv",
        "results/enrichment/KEGG_results.csv",
        # Relatório final
        "results/plots/RELATORIO_FINAL.html"

########################################
# REFERÊNCIA
#########################################
# DOWNLOAD DAS REFERÊNCIAS
# Baixar os arquivos de referência do transcriptoma (.fa)
rule download_transcriptome:
    output:
        "reference/transcriptome.fa.gz"

    shell:
        """
        mkdir -p reference

        wget \
            -O {output} \
            {FASTA_URL}
        """

rule decompress_transcriptome:
    input:
        "reference/transcriptome.fa.gz"

    output:
        TRANSCRIPTOME

    shell:
        """
        gunzip -c {input} > {output}
        """

# Baixar anotação dos genes (.gtf)
rule download_gtf:
    output:
        "reference/annotation.gtf.gz"

    shell:
        """
        wget \
            -O {output} \
            {GTF_URL}
        """


rule decompress_gtf:
    input:
        "reference/annotation.gtf.gz"

    output:
        ANNOTATION

    shell:
        """
        gunzip -c {input} > {output}
        """
########################################
# TX2GENE
rule generate_tx2gene:
    input:
        ANNOTATION

    output:
        TX2GENE

    conda:
        "envs/rnaseq.yaml"

    script:
        "scripts/generate_tx2gene.R"
#####################################
rule salmon_index:
    input:
        fasta=TRANSCRIPTOME
    output:
        directory("reference/salmon_index")
    conda:
        "envs/salmon.yaml"
    threads:
        config["resources"]["salmon"]["threads"]
    shell:
        """
        salmon index \
            -t {input.fasta} \
            -i {output} \
            -p {threads}
        """

#########################################
# DADOS BRUTOS
#########################################
rule download_fastq:
    output:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"

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
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"
    output:
        html1 = "results/fastqc/{sample}_1_fastqc.html",
        html2 = "results/fastqc/{sample}_2_fastqc.html",
        zip1 = "results/fastqc/{sample}_1_fastqc.zip",
        zip2 = "results/fastqc/{sample}_2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    threads: 2
    shell:
        "fastqc -t {threads} -o results/fastqc/ {input.r1} {input.r2}"

###################################
rule fastp:
    input:
        r1 = "data/raw/{sample}_1.fastq.gz",
        r2 = "data/raw/{sample}_2.fastq.gz"
    output:
        r1 = "results/fastqc/trimmed/{sample}_1_trimmed.fastq.gz",
        r2 = "results/fastqc/trimmed/{sample}_2_trimmed.fastq.gz",
        json = "results/fastqc/trimmed/{sample}_fastp.json",
        html = "results/fastqc/trimmed/{sample}_fastp.html"
    conda:
        "envs/fastp.yaml"
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -j {output.json} -h {output.html} \
              --thread {threads} --detect_adapter_for_pe
        """

##########################################
rule fastqc_trimmed:
    input:
        r1 = "results/fastqc/trimmed/{sample}_1_trimmed.fastq.gz",
        r2 = "results/fastqc/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        html1 = "results/fastqc/trimmed/{sample}_1_trimmed_fastqc.html",
        html2 = "results/fastqc/trimmed/{sample}_2_trimmed_fastqc.html",
        zip1 = "results/fastqc/trimmed/{sample}_1_trimmed_fastqc.zip",
        zip2 = "results/fastqc/trimmed/{sample}_2_trimmed_fastqc.zip"
    conda:
        "envs/fastqc.yaml"  # Reutiliza o mesmo ambiente que você já criou para o FastQC!
    threads: 2
    shell:
        "fastqc -t {threads} -o results/fastqc/trimmed/ {input.r1} {input.r2}"

################################################
rule multiqc:
    input:
        expand("results/fastqc/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("results/fastqc/{sample}_2_fastqc.html", sample=SAMPLES)
    output:
        "results/fastqc/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc results/fastqc/ -o results/fastqc/ -n multiqc_report.html"
##################################################################
rule multiqc_trimmed:
    input:
        expand(
            "results/fastqc/trimmed/{sample}_1_trimmed_fastqc.html", sample=SAMPLES),
        expand(
            "results/fastqc/trimmed/{sample}_2_trimmed_fastqc.html", sample=SAMPLES)
    output:
        "results/fastqc/trimmed/multiqc_trimmed_report.html"
    conda:
        "envs/multiqc.yaml"  # Reutiliza o ambiente do MultiQC!
    shell:
        "multiqc results/fastqc/trimmed/ -o results/fastqc/trimmed/ -n multiqc_trimmed_report.html"


#########################################
# QUANTIFICAÇÃO
#########################################
rule salmon_quant:
    input:
        index = "reference/salmon_index",
        r1 = "results/fastqc/trimmed/{sample}_1_trimmed.fastq.gz",
        r2 = "results/fastqc/trimmed/{sample}_2_trimmed.fastq.gz"
    output:
        sf = "results/salmon/{sample}/quant.sf",
        dir = directory("results/salmon/{sample}")
    conda:
        "envs/salmon.yaml"
    threads: 6
    shell:
        """
        salmon quant -i {input.index} -l A \
            -1 {input.r1} -2 {input.r2} \
            -p {threads} --validateMappings \
            -o {output.dir}
        """

#########################################
# ANÁLISE ESTATÍSTICA
#########################################
rule tximport:
    input:
        sf = expand("results/salmon/{sample}/quant.sf", sample=SAMPLES),
        tx2gene = TX2GENE,
        samples = "samples.tsv"
    output:
        txt = "results/deseq2/tximport_complete.txt",
        rds = "results/deseq2/txi.rds"
    conda:
        "envs/rnaseq.yaml"
    script:
        "scripts/tximport.R"
###############################################
rule deseq2:
    input:
        rds = "results/deseq2/txi.rds",
        samples = "samples.tsv"
    output:
        # Tabela final contendo os genes, p-valores e log2FoldChange
        csv = "results/deseq2/DEG_results.csv"
    conda:
        "envs/rnaseq.yaml"
    script:
        "scripts/deseq2.R"


#########################################
# ENRIQUECIMENTO
#########################################
rule enrichment:
    input:
        csv = "results/deseq2/DEG_results.csv"
    output:
        go = "results/enrichment/GO_results.csv",
        kegg = "results/enrichment/KEGG_results.csv"
    conda:
        "envs/rnaseq.yaml"
    script:
        "scripts/enrichment.R"

#########################################
# RELATÓRIO CIENTÍFICO FINAL
#########################################
rule generate_report:
    input:
        csv = "results/deseq2/DEG_results.csv",
        dds = "results/deseq2/dds_object.rds",
        go = "results/enrichment/GO_results.csv",
        kegg = "results/enrichment/KEGG_results.csv"
    output:
        html = "results/plots/RELATORIO_FINAL.html"
    conda:
        "envs/rnaseq.yaml"
    script:
        "scripts/render_report.R"

