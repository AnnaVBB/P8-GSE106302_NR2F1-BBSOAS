# 💻 Pipeline de RNA-seq para Análise da Síndrome BBSOAS (GSE106302) 🧠

Este repositório contém um pipeline de bioinformática automatizado e reprodutível para análise de dados de RNA-seq do hipocampo de camundongos (*Mus musculus*). O objetivo é investigar as alterações transcriptômicas decorrentes da haploinsuficiência do fator de transcrição **Nr2f1**, modelo biológico associado à Síndrome de Bosch-Boonstra-Schaaf (BBSOAS).

## 🧬 Contexto Biológico 🧠

A Síndrome BBSOAS é uma doença neurodesenvolvimental rara causada por mutações autossômicas dominantes ou deleções no gene *NR2F1*. Pacientes apresentam atraso global no desenvolvimento, deficiência intelectual e defeitos visuais corticais. Este pipeline reanalisa os dados públicos do dataset **GSE106302** (amostras HET vs WT) para identificar os alvos moleculares downstream afetados.

## 🟣 Arquitetura do Pipeline e Ferramentas

O fluxo de trabalho foi construído utilizando a ferramenta de gerência de fluxos **Snakemake** combinada com ambientes isolados **Conda** para garantir a máxima reprodutibilidade científica.

*   **Download de Dados:** `fasterq-dump` (SRA Toolkit) para obtenção automática das amostras.
*   **Controle de Qualidade (QC):** `FastQC` e `MultiQC` para avaliação das leituras brutas e processadas.
*   **Tratamento de Leituras (Trimming):** `fastp` para remoção de adaptadores e filtros de qualidade por janela deslizante.
*   **Quantificação:** `Salmon` (indexação e pseudo-alinhamento) utilizando o transcriptoma de referência do camundongo GRCm38.
*   **Expressão Diferencial:** `DESeq2` (R/Bioconductor) para normalização, cálculo de variância e identificação de DEGs.
*   **Enriquecimento Funcional:** `clusterProfiler` e `enrichplot` para análises de Ontologia Gênica (GO), Vias KEGG e GSEA.
*   **Relatório Consolidado:** `R Markdown` compilado em um documento único HTML interativo.

## 📁 Estrutura do Repositório

```text
├── config.yaml          # Configurações de recursos e threads do pipeline
├── samples.tsv          # Metadados e identificadores das amostras SRA
├── Snakefile            # Regras e automação do pipeline Snakemake
├── relatorio.Rmd        # Documento fonte do relatório científico final
├── envs/                # Ambientes virtuais Conda (.yaml) para cada ferramenta
└── scripts/             # Scripts auxiliares em linguagem R (DESeq2, Enriquecimento)
```
## 🟣 Como visualizar os resultados?
Para analisar a discussão científica, as tabelas de impacto clínico e as figuras geradas (PCA, Volcano Plot, Heatmap e GSEA), você não precisa instalar nada:
*   Acesse [clique aqui para acessar o RELATORIO_FINAL.html](https://annavbb.github.io/P8-GSE106302_NR2F1-BBSOAS/index.html)
   
## 🟣 Como Reproduzir este Projeto

### 1. Pré-requisitos
Certifique-se de ter o Conda (ou Miniconda) e o Snakemake instalados em seu ambiente Linux.

### 2. Clonar o Repositório
```bash
git clone https://github.com
cd P8-GSE106302_NR2F1-BBSOAS
```

### 3. Preparar as Referências
Insira o arquivo fasta do transcriptoma em `reference/transcriptome.fa` e a tabela de conversão em `reference/tx2gene.csv`.
#### 3.1 Obtenha os metadados 
Através do estudo GEO, crie um arquivo .tsv com os metadados do estudo. 
NOTA: para outras amostras, altere o arquivo "samples.tsv", que pode ser gerado com ajuda do script "gerar_metadados.py", baseando-se nos arquivos (Summary e RunInfo) do estudo GEO desejado.  

#### 3.2 Baixe as referencias 
Obtenha as referencias em https://ftp.ensembl.org/pub/ 
```bash
wget https://ftp.ensembl.org/pub/release-116/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz
gunzip Mus_musculus.GRCm39.cdna.all.fa.gz
mv Mus_musculus.GRCm39.cdna.all.fa transcriptome.fa
(resultado: reference/transcriptome.fa)
```
```bash
wget https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz 
gunzip Mus_musculus.GRCm39.*.gtf.gz
mv Mus_musculus.GRCm39.*.gtf annotation.gtf
(resultado: reference/annotation.gtf)
```
### 3.3 Gere o tx2gene.csv 
Com o script "generate_tx2gene.R" e o "annotation.gtf" gere o arquivo .csv
```bash
# Criar um ambiente separado para o R:
conda deactivate
conda create -n rnaseq -c conda-forge -c bioconda r-base bioconductor-genomicfeatures bioconductor-txdbmaker -y
conda activate rnaseq
Rscript scripts/generate_tx2gene.R

#para voltar para o outro ambiente
conda deactivate
conda env list
```

### 4. Executar o Pipeline Completo
O Snakemake resolverá todas as dependências, criará os ambientes virtuais isolados e processará os dados em lote:
```bash
snakemake --use-conda -j 4
```

## 📊 Principais Descobertas e Métricas
*   **Eficiência de Alinhamento:** O mapeamento pelo Salmon apresentou altas taxas de eficiência (>80%), confirmando a adequação do transcriptoma de referência escolhido.
*   **Perfil Transcritômico:** Adotando critérios estatísticos rigorosos ($p\text{-ajustado} < 0.05$ e $|\log_2\text{FC}| > 1$), foram identificados **221 Genes Diferencialmente Expressos (DEGs)**, revelando uma marcante assimetria biológica voltada para a regulação positiva (215 genes *up-regulated* e apenas 6 *down-regulated*).
*   **Assinatura Funcional (GSEA/KEGG):** Os resultados demonstraram o enriquecimento de vias estruturais de matriz extracelular (interação ECM-receptor), citoesqueleto e montagem ciliar. Contudo, observou-se a supressão severa da via de transporte anterógrado intraciliar. Essa desconexão molecular sugere a formação de cílios primários neuronais disfuncionais no hipocampo mutante, fornecendo uma justificativa mecânica para o deficit cognitivo da síndrome.

## 📂 Resultados Gerados
Após a execução, o pipeline estruturará os resultados na pasta `results/`:
*   `results/fastqc/`: Relatórios de qualidade brutos e pós-trimming.
*   `results/salmon/`: Arquivos de quantificação de expressão gênica.
*   `results/deseq2/`: Tabela de DEGs (`DEG_results.csv`) e objeto binário (`dds_object.rds`).
*   `results/enrichment/`: Tabelas de termos funcionais GO e KEGG.
*   **`results/plots/RELATORIO_FINAL.html`**: Relatório científico definitivo contendo todos os textos, tabelas de impacto clínico e gráficos embutidos (PCA, Volcano Plot com símbolos oficiais, Heatmap e GSEA).
