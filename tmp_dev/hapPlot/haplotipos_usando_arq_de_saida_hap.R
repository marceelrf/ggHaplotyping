#Script para tentar criar um haplotipo de KIR

#pacotes
library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(readr)
library(reshape2) #melt
library(gggenes)
library(plotly)

# Criar um vetor de cores para cada gene
gene_colors <- c("KIR2DL1" = "#CDB5CD", "KIR2DL2" = "#ADFF2F", 
                 "KIR2DL3" = "#436EEE", "KIR2DL4" = "#8B8B7A", 
                 "KIR2DL5AB" = "#A52A2A", "KIR2DP1" = "#C1FFC1", 
                 "KIR2DS1" = "#00F5FF", "KIR2DS2" = "#EE6A50", 
                 "KIR2DS3" = "#CDC673", "KIR2DS4" = "orange", 
                 "KIR2DS5" = "#EEA2AD", "KIR3DL1" = "#FF83FA", 
                 "KIR3DL2" = "#AB82FF", "KIR3DL3" = "#B0C4DE", 
                 "KIR3DP1" = "#FFFF00", "KIR3DS1" = "#008B00")

# Criar um mapeamento de cores para 0 e 1 para cada gene
color_mapping <- unlist(sapply(names(gene_colors), function(gene) {
  c(`1` = gene_colors[[gene]], `0` = "white")
}))

#----ARQUIVO DE ENTRADA
table <- read.table("haplotype_all_merged_two_chromosomes_perl_line.db.txt", h = T, sep = "\t")

#filtrando uma unica amostra para ver se funciona o codigo
#sample <-  c("446303051", "1016611301", "210800083", "1012617701", "926224051", "423320091")
# Filtrando as amostras específicas
#data_filtered <- table %>% 
 # filter(SAMPLE %in% sample)


#Filtrando as colunas necessarias para o app -H1 e -H2 de cada amostra
table_filtered <- table %>% 
                  select(Sample, matches("^KIR.*h[12]")) %>%       # Seleciona "Sample" e colunas que começam com "KIR" e contêm "h1" E "h2"
                  rename_with(~ str_replace_all(., "\\.+", "_") %>% str_replace("_$", ""), matches("^KIR.*h[12]")) %>%  # Substitui "." por "_" e remove "_" final
                  mutate(across(starts_with("KIR"), 
                    ~ ifelse(str_detect(., "\\*null$"), 0, 1),   # Verifica se o valor contém "*null" e substitui por 0, senão 1
                    .names = "presence_{.col}"))   # Cria as colunas de presença


# Transformando para o formato longo

#-------------------------PRESENCA/AUSENCIA
presence_long <- table_filtered %>%
                 select(Sample, matches("presence_KIR.*_h[12]")) %>%  # Seleciona as colunas de presença/ausência
                 pivot_longer(cols = matches("presence_KIR.*_h[12]"),  # Transforma as colunas de presença/ausência
                             names_to = "Name",  # Coloca o nome completo na coluna 'Name'
                             values_to = "Presence") %>%  # Coloca os valores de presença na coluna 'Presence'
                 mutate(Name = str_replace(Name, "presence_", "")) %>%  # Remove o prefixo "presence_"
                 separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
                 select(Sample, Gene, Haplo, Presence)  # Seleciona apenas as colunas desejadas

#--------------------------ALELOS
alleles_long <- table_filtered %>%
                select(Sample, matches("^KIR.*_h[12]$")) %>%  # Seleciona as colunas de alelos (h1 e h2)
                pivot_longer(cols = matches("^KIR.*_h[12]$"),  # Transforma as colunas de alelos
                names_to = "Name",  # Coloca o nome completo na coluna 'Name'
                values_to = "Alelo") %>%  # Coloca os valores dos alelos na coluna 'Alelo'
                separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
                select(Sample, Gene, Haplo, Alelo) %>%  # Seleciona apenas as colunas desejadas
                mutate(Alelo = str_remove(Alelo, "\\(version \\d+\\)")) %>%  #remove o version de todos os alelos
                mutate(Alelo = str_remove(Alelo, ".*(?=\\*)")) %>% 
                mutate(Alelo = str_replace(Alelo, "(\\w+)\\.(\\w+)\\.(.+)", "\\1*\\2\\3")) %>% #antes era KIR3DL2.new.30 por exemplo, então retirou o ".", substitui e fica assim KIR3DL2*new30
                mutate(Alelo = str_remove(Alelo, "^\\w+")) #Remove o nome do gene no início, deixando apenas o que vem após o *

MANTEM_os_alelos <-  table_filtered %>%
  select(Sample, matches("^KIR.*_h[12]$")) %>%  # Seleciona as colunas de alelos (h1 e h2)
  pivot_longer(cols = matches("^KIR.*_h[12]$"),  # Transforma as colunas de alelos
               names_to = "Name",  # Coloca o nome completo na coluna 'Name'
               values_to = "Alelo") %>%  # Coloca os valores dos alelos na coluna 'Alelo'
  separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
  select(Sample, Gene, Haplo, Alelo) %>%  # Seleciona apenas as colunas desejadas
  mutate(Alelo = str_remove(Alelo, "\\(version \\d+\\)")) %>%  #remove o version de todos os alelos
  mutate(Alelo = str_replace(Alelo, "(\\w+)\\.(\\w+)\\.(.+)", "\\1*\\2\\3"))  #antes era KIR3DL2.new.30 por exemplo, então retirou o ".", substitui e fica assim KIR3DL2*new30




# Juntando os dados de Presence e Alelo
combined_data <- left_join(presence_long, MANTEM_os_alelos, by = c("Sample", "Gene", "Haplo"))
combined_data <- left_join(presence_long, alleles_long, by = c("Sample", "Gene", "Haplo"))

#---------GRAFICOOOOOOOO
# Como plotar a legenda com a cor dos genes KIR na ordem
plot_data_5_samples <- combined_data %>%
                       filter(Sample %in% sample(unique(combined_data$Sample), 5)) %>%  # Seleciona 25 amostras aleatórias
                       mutate(Gene_Haplo = paste(Gene, Haplo, sep = "_")) %>%
  mutate(
   xpos = case_when(
     Gene == "KIR3DL3" ~ 1,
     Gene == "KIR2DS2" ~ 2.7,  # Aumenta a distância para acomodar caixas alargadas
     Gene == "KIR2DL2" ~ 4.2,
     Gene == "KIR2DL3" ~ 4.2,  # Mesma posição para KIR2DL2 e KIR2DL3
     Gene == "KIR2DL5AB" ~ 5.7,
     Gene == "KIR2DS3" ~ 7.2,  # Mesma posição para KIR2DS3 e KIR2DS5
     Gene == "KIR2DS5" ~ 7.2,
     Gene == "KIR2DP1" ~ 8.7,
     Gene == "KIR2DL1" ~ 10.2,
     Gene == "KIR3DP1" ~ 11.7,
     Gene == "KIR2DL4" ~ 13.2,
     Gene == "KIR3DL1" ~ 14.7,  # Mesma posição para KIR3DL1 e KIR3DS1
     Gene == "KIR3DS1" ~ 14.7,
     Gene == "KIR2DS4" ~ 16.2,  # Mesma posição para KIR2DS4 e KIR2DS1
     Gene == "KIR2DS1" ~ 16.2,
     Gene == "KIR3DL2" ~ 17.7,
     TRUE ~ NA_real_
),
    start = xpos,
    end = xpos + 1
  ) %>% 
  mutate(Gene = factor(Gene, levels = unique(Gene)))


plot_data_5_samples %>%
  filter(Presence == 1) %>%  # Filtra os dados para incluir apenas onde Presence == 1
  mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                        "H1", "H2"),
         Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%  # Cria uma nova coluna combinando Sample e Haplo
  mutate(Group = rep(1:ceiling(n()/5), each = 7, length.out = n())) %>%  # Cria grupos de 5 amostras
  ggplot(aes(xmin = start,
             xmax = end, 
             y = Sample, 
             fill = Gene)) +
  geom_segment(aes(x = 0.9, xend = 18.5),  # Linha horizontal 
               color = "black", size = 0.5) +
  geom_rect(aes(xmin = start - 0.2,  # Alarga 0.6 à esquerda
                xmax = end + 0.2,    # Alarga 0.6 à direita
                ymin = as.numeric(Sample) - 0.01,  # Altura da caixa
                ymax = as.numeric(Sample) + 0.01),  # Altura da caixa
            color = "black", 
            linewidth = 0.5) +
 geom_gene_label(aes(label = Alelo),
                 size = 20) +  # Ajusta o tamanho do texto no rótulo
  facet_wrap(Sample ~ Haplo, 
             scales = "free_y", 
             ncol = 1) +  # Facetar por amostra e haplótipo (H1 e H2) em linhas
  scale_fill_manual(values = gene_colors) +
  theme_minimal() +
  guides(fill = guide_legend(nrow = 2), ncol = 2) +  # Deixa os itens da legenda na horizontal
  labs(x = "Genes", 
       y = "Sample", 
       title = "KIR haplotypes (H1 and H2) for each sample") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        strip.text = element_blank()) +
  geom_label(aes(x = -1.5, 
                 label = Sample_Haplo),  # Usa a nova coluna Sample_Haplo para mostrar Sample e Haplo
             size = 3, 
             fill = NA,  # Remove a cor de fundo da caixa do nome da amostra
             show.legend = FALSE)  # Adiciona rótulo com o nome da amostra e haplótipo

plot_data_5_samples %>%
  filter(Presence == 1) %>%  # Filtra os dados para incluir apenas onde Presence == 1
  mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                        "H1", "H2"),
         Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%  # Cria uma nova coluna combinando Sample e Haplo
  mutate(Group = rep(1:ceiling(n()/5), each = 5, length.out = n())) %>%  # Cria grupos de 5 amostras
  ggplot(aes(xmin = start,
             xmax = end, 
             y = Sample, 
             fill = Gene)) +
  geom_segment(aes(x = 0.9, xend = 18.5),  # Linha horizontal 
               color = "black", size = 0.5) +
  geom_rect(aes(xmin = start - 0.2,  # Alarga 0.6 à esquerda
                xmax = end + 0.2,    # Alarga 0.6 à direita
                ymin = as.numeric(Sample) - 0.01,  # Altura da caixa
                ymax = as.numeric(Sample) + 0.01),  # Altura da caixa
            color = "black", 
            linewidth = 0.5) +
  geom_gene_label(aes(label = Alelo),
                  size = 30) +  # Ajusta o tamanho do texto no rótulo
  facet_wrap(Sample ~ Haplo, 
             scales = "free_y", 
             ncol = 1) +  # Facetar por amostra e haplótipo (H1 e H2) em linhas
  scale_fill_manual(values = gene_colors) +
  theme_minimal() +
  guides(fill = guide_legend(nrow = 2), ncol = 2) +  # Deixa os itens da legenda na horizontal
  labs(x = "Genes", 
       y = "Sample", 
       title = "KIR haplotypes (H1 and H2) for each sample") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_blank()) +
  geom_label(aes(x = -1.5, 
                 label = Sample_Haplo),  # Usa a nova coluna Sample_Haplo para mostrar Sample e Haplo
             size = 3, 
             fill = NA,  # Remove a cor de fundo da caixa do nome da amostra
             show.legend = FALSE)  # Adiciona rótulo com o nome da amostra e haplótipo

#-geom_rect()#-------------------------------------------------------------------------------------------------------------------------------------
#filtrando uma unica amostra para ver se funciona o codigo
#sample <-  c("446303051", "1016611301", "210800083", "1012617701", "926224051", "423320091")
# Filtrando as amostras específicas
#data_filtered <- table %>% 
# filter(SAMPLE %in% sample)