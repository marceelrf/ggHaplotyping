#Script para tentar criar um haplotipo de KIR
install.packages("reshape2")
install.packages("tidyverse")
#pacotes
library(tidyverse)
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(readr)
library(reshape2) #melt

table <- read.table("presence.table.txt", h = T, sep = "\t")

#filtrando uma unica amostra para ver se funciona o codigo
sample <-  c("446303051", "1016611301", "210800083", "1012617701", "926224051", "423320091")
# Filtrando as amostras específicas
data_filtered <- table %>% 
  filter(SAMPLE %in% sample)


# Transformar os dados para o formato longo
data_long <- melt(data_filtered, id.vars = "SAMPLE") %>% 
             filter(variable != "X5UPKIR", variable != "HLA.E", variable != "HLA.G") %>% 
             filter(!grepl("^LIL", variable)) 

# Definir SAMPLE como um fator para garantir que seja tratado como discreto
data_long$SAMPLE <- factor(data_long$SAMPLE, levels = unique(data_long$SAMPLE))


# Criar um vetor de cores para cada gene
gene_colors <- c("KIR2DL1" = "#CDB5CD", "KIR2DL2" = "#ADFF2F", 
                 "KIR2DL3" = "#436EEE", "KIR2DL4" = "#8B8B7A", 
                 "KIR2DL5AB" = "#A52A2A", "KIR2DP1" = "black", 
                 "KIR2DS1" = "#00F5FF", "KIR2DS2" = "#EE6A50", 
                 "KIR2DS3" = "#CDC673", "KIR2DS4" = "orange", 
                 "KIR2DS5" = "#EEA2AD", "KIR3DL1" = "#FF83FA", 
                 "KIR3DL2" = "#AB82FF", "KIR3DL3" = "#B0C4DE", 
                 "KIR3DP1" = "#FFFF00", "KIR3DS1" = "#008B00")

# Criar um mapeamento de cores para 0 e 1 para cada gene
color_mapping <- unlist(sapply(names(gene_colors), function(gene) {
  c(`1` = gene_colors[[gene]], `0` = "white")
}))

# Criar o gráfico
data_long %>% 
  ggplot(aes(x = variable, y = SAMPLE)) +
  geom_tile(aes(fill = interaction(value, variable)), color = "white") +  # Interação de value e variable
  scale_fill_manual(values = color_mapping) +
  theme_minimal() +
  labs(x = "Genes", y = "Amostra", fill = "Presença") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_discrete(limits = unique(data_long$SAMPLE))  # Define as amostras únicas


# Criar o gráfico
table %>% 
  select(SAMPLE, starts_with("KIR")) %>% 
  slice_sample(n = 25) %>% 
  mutate(ypos = row_number(.)) %>% 
  pivot_longer(starts_with("KIR"), names_to = "Gene",values_to = "Presence") %>% 
  filter(Presence == 1) %>% 
  mutate(Gene = factor(Gene)) %>% 
  mutate(xpos = case_when(
    Gene == "KIR3DL3" ~ 1,
    Gene == "KIR2DS2" ~ 2,
    Gene == "KIR2DL2" ~ 3,
    Gene == "KIR2DL3" ~ 3,
    Gene == "KIR2DL5AB" ~ 4,
    Gene == "KIR2DS3" ~ 5,
    Gene == "KIR2DS5" ~ 5,
    Gene == "KIR2DP1" ~ 6,
    Gene == "KIR2DL1" ~ 7,
    Gene == "KIR3DP1" ~ 8,
    Gene == "KIR2DL4" ~ 9,
    Gene == "KIR3DL1" ~ 10,
    Gene == "KIR3DS1" ~ 10,
    Gene == "KIR2DL5AB" ~ 11,
    Gene == "KIR2DS3" ~ 12,
    Gene == "KIR2DS5" ~ 12,
    Gene == "KIR2DS4" ~ 13,
    Gene == "KIR2DS1" ~ 13,
    Gene == "KIR3DL2" ~ 14,
    TRUE~NA
  )) %>% 
  mutate(SAMPLE = factor(SAMPLE)) %>% 
  ggplot(aes(y=ypos,fill = Gene)) +
  geom_hline(aes(yintercept = ypos)) +
  annotate(geom = "rect",xmin = .5, xmax = 8.4,
           ymin=.5,ymax = 25+.5, alpha = .25,
           fill = "grey",
           color = "grey") +
  annotate(geom = "rect",xmin = 8.6, xmax = 14.5,
           ymin=.5,ymax = 25+.5, alpha = .25,
           fill = "black",
           color = "black") +
  geom_rect(aes(xmin = xpos -.25, xmax = xpos+.25, ymin=ypos-.25,ymax=ypos+.25),
            alpha = 1,
            color = "black") +
  scale_fill_manual(values = gene_colors) +
  theme_minimal() +
  labs(x = "Genes", y = "Amostra", fill = "Presença",
       title = "KIR haplotypes - 25 samples") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()) +
  geom_label(aes(x = -1.5,y = ypos, label = SAMPLE, fill = NULL),
             show.legend = F)

ggsave(filename = "kir_25_samples.png",bg = "white")
