# Carrega explicitamente os pacotes gráficos para evitar falhas de inicialização no ambiente do terminal
library(rmarkdown)
library(EnhancedVolcano)
library(ggplot2)

message("Iniciando a renderizacao automatica do relatorio.Rmd...")

# Executa o processamento do documento e direciona a saida para a pasta correta
rmarkdown::render(
    input = "relatorio.Rmd",
    output_file = "RELATORIO_FINAL.html",
    output_dir = "results/plots",
    clean = TRUE
)

message("Relatorio gerado com sucesso em results/plots/RELATORIO_FINAL.html!")
