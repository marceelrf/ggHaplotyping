if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
library(pacman)


# List of packages to load
packages <- c("shiny",
              "tidyverse",
              "ggplot2", 
              "plotly",
              "gggenes",
              "viridis",
              "reshape2",
              "shinyalert")

# Function to install and load packages if missing
p_load(packages, character.only = TRUE)

# Função para processar o arquivo e gerar o gráfico
processar_gerar_grafico <- function(input_file) {
  
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
  
  # Ler o arquivo de entrada
  table <- read.table(input_file, h = T, sep = "\t")
  
  # Filtrando as colunas necessárias para o app -H1 e -H2 de cada amostra
  table_filtered <- table %>% 
    select(Sample, matches("^KIR.*h[12]")) %>%  # Seleciona "Sample" e colunas que começam com "KIR" e contêm "h1" e "h2"
    rename_with(~ str_replace_all(., "\\.+", "_") %>% str_replace("_$", ""), matches("^KIR.*h[12]")) %>%  # Substitui "." por "_"
    mutate(across(starts_with("KIR"), 
                  ~ ifelse(str_detect(., "\\*null$"), 0, 1),   # Verifica se o valor contém "*null" e substitui por 0
                  .names = "presence_{.col}"))  # Cria as colunas de presença
  
  # Transformando para o formato longo
  presence_long <- table_filtered %>%
    select(Sample, matches("presence_KIR.*_h[12]")) %>%  # Seleciona as colunas de presença/ausência
    pivot_longer(cols = matches("presence_KIR.*_h[12]"),  # Transforma as colunas de presença/ausência
                 names_to = "Name",  # Coloca o nome completo na coluna 'Name'
                 values_to = "Presence") %>%  # Coloca os valores de presença na coluna 'Presence'
    mutate(Name = str_replace(Name, "presence_", "")) %>%  # Remove o prefixo "presence_"
    separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
    select(Sample, Gene, Haplo, Presence)  # Seleciona apenas as colunas desejadas
  
  alleles_long <- table_filtered %>%
    select(Sample, matches("^KIR.*_h[12]$")) %>%  # Seleciona as colunas de alelos (h1 e h2)
    pivot_longer(cols = matches("^KIR.*_h[12]$"),  # Transforma as colunas de alelos
                 names_to = "Name",  # Coloca o nome completo na coluna 'Name'
                 values_to = "Alelo") %>%  # Coloca os valores dos alelos na coluna 'Alelo'
    separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
    select(Sample, Gene, Haplo, Alelo) %>%  # Seleciona apenas as colunas desejadas
    mutate(Alelo = str_remove(Alelo, "\\(version \\d+\\)")) %>%  # Remove a versão dos alelos
    mutate(Alelo = str_remove(Alelo, ".*(?=\\*)")) %>% 
    mutate(Alelo = str_replace(Alelo, "(\\w+)\\.(\\w+)\\.(.+)", "\\1*\\2\\3")) %>% # Ajusta formato dos alelos
    mutate(Alelo = str_remove(Alelo, "^\\w+")) # Remove o nome do gene no início
  
  combined_data <- left_join(presence_long, alleles_long, by = c("Sample", "Gene", "Haplo"))
  
  plot_data_5_samples <- combined_data %>%
    filter(Sample %in% head(unique(combined_data$Sample), 5)) %>%
    mutate(Gene_Haplo = paste(Gene, Haplo, sep = "_")) %>%
    mutate(
      xpos = case_when(
        Gene == "KIR3DL3" ~ 1,
        Gene == "KIR2DS2" ~ 2.7, 
        Gene == "KIR2DL2" ~ 4.2,
        Gene == "KIR2DL3" ~ 4.2,  
        Gene == "KIR2DL5AB" ~ 5.7,
        Gene == "KIR2DS3" ~ 7.2,  
        Gene == "KIR2DS5" ~ 7.2,
        Gene == "KIR2DP1" ~ 8.7,
        Gene == "KIR2DL1" ~ 10.2,
        Gene == "KIR3DP1" ~ 11.7,
        Gene == "KIR2DL4" ~ 13.2,
        Gene == "KIR3DL1" ~ 14.7,  
        Gene == "KIR3DS1" ~ 14.7,
        Gene == "KIR2DS4" ~ 16.2,  
        Gene == "KIR2DS1" ~ 16.2,
        Gene == "KIR3DL2" ~ 17.7,
        TRUE ~ NA_real_
      ),
      start = xpos,
      end = xpos + 1
    ) %>% 
    mutate(Gene = factor(Gene, levels = unique(Gene)))
  
  plot_presence <- plot_data_5_samples %>%
    filter(Presence == 1) %>%  # Filtra os dados para incluir apenas onde Presence == 1
    mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                          "H1", "H2"),
           Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%  # Cria uma nova coluna combinando Sample e Haplo
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
    # geom_gene_label(aes(label = Alelo),
    #                size = 20) +  # Ajusta o tamanho do texto no rótulo
    facet_wrap(Sample ~ Haplo, 
               scales = "free_y", 
               ncol = 1) +  # Facetar por amostra e haplótipo (H1 e H2) em linhas
    scale_fill_manual(values = gene_colors) +
    theme_minimal() +
    guides(fill = guide_legend(nrow = 2), ncol = 2) +  # Deixa os itens da legenda na horizontal
    labs(x = "Genes", 
         y = "Sample", 
         title = "") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          strip.text = element_blank(),
          legend.text = element_text(size = 18),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 20)) +
    geom_label(aes(x = -1.5, 
                   label = Sample_Haplo),  # Usa a nova coluna Sample_Haplo para mostrar Sample e Haplo
               size = 3.8,  # Reduz o tamanho do texto
               fill = NA,  # Remove a cor de fundo da caixa do nome da amostra
               show.legend = FALSE,
               nudge_x = 1.3)  # Ajusta a posição do texto para a direita, se necessário
  
  
  plot_alleles <- plot_data_5_samples %>%
    filter(Presence == 1) %>%  # Filtra os dados para incluir apenas onde Presence == 1
    mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                          "H1", "H2"),
           Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%  # Cria uma nova coluna combinando Sample e Haplo
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
         title = "") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.position = "top",
          strip.text = element_blank(),
          legend.text = element_text(size = 18),
          plot.margin = margin(t = 5, r = 10, b = 5, l = 20)) +
    geom_label(aes(x = -1.5, 
                   label = Sample_Haplo),  # Usa a nova coluna Sample_Haplo para mostrar Sample e Haplo
               size = 3.8,  # Reduz o tamanho do texto
               fill = NA,  # Remove a cor de fundo da caixa do nome da amostra
               show.legend = FALSE,
               nudge_x = 1.3)  # Ajusta a posição do texto para a direita, se necessário
  
  return(plots = list(plot_presence = plot_presence,plot_alleles = plot_alleles))
}


ui <- fluidPage(
  tags$div(
    style = "display: flex; align-items: center; margin-bottom: 10px;",
    tags$h1("KIR Haplotypes App", style = "margin-right: 15px;"),
    tags$img(src = "https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSFDRAnTfp1qWUVX1BaXrlLze1pDJImh3nr7A&s", 
             height = "75px")
  ),
  sidebarLayout( 
    sidebarPanel( #criar o painel do lado esquerdo
      width = 3,  # Ajusta o tamanho do painel lateral (3 é a largura em colunas do layout Bootstrap)
      fileInput("file1", "Upload a data file", accept = c(".txt"), width = "100%"),  # Largura do painel
      actionButton("goButton", "Generate Graph", style = "width: 100%;"),
      tags$hr(), #add uma linha horizontal
      
      radioButtons("select", "Select which information you want:",  #montar a caixinha para selecionar
                   choices = c("Only Presence", "With the alleles"),
                   selected = "Only Presence"), # "Plot Presence" é a opção selecionada por padrão
      tags$hr(), #add uma linha horizontal
      
      textInput(inputId = "txt", label = "Define the title of the haplotype:", value = "KIR haplotypes (h1 and h2) for each sample"), #caixinha para o usr escrever o titulo
      
      tags$hr(),
      
      downloadButton("savePlot", "Save Current Plot", style = "width: 100%;")
      ),
    
    mainPanel( #painel principal
      h3(textOutput("txt", container = span)), #escreve o titulo
      plotOutput("grafico", height = "800px", width = "100%")  # Tamanho do gráfico
    )
  )
)


server <- function(input, output)
{
  # Armazenar a lista de gráficos reativamente
  stored_plots <- reactiveVal(NULL)  # Para armazenar os gráficos
  
  observeEvent(input$goButton, {
    
    req(input$file1)  # Verifica se o arquivo foi carregado
    
    # Iniciar a barra de progresso
    withProgress(message = 'Checking columns...', value = 0, {
      
      # Atualiza o progresso para 50% (metade do processo)
      incProgress(0.5, detail = "Processing ...")  
      
      # Ler o arquivo
      file_data <- read.table(input$file1$datapath, header = TRUE, sep = "\t", check.names = TRUE)
      
      # Atualiza o progresso para 30% (arquivo lido)
      incProgress(0.2, detail = "Reading the file...")
      
      # Identificar as colunas que começam com "KIR" e seguem o padrão esperado
      # Ajuste da expressão regular para pegar as colunas no formato 'KIR3DL3..h1.' ou similar
      kir_columns <- grep("^KIR\\d+[A-Z]*\\.\\.h[12]\\.$", names(file_data), value = TRUE)
      
      # Definir as colunas obrigatórias
      required_columns <- c("Sample", kir_columns)
      
      # Checar se todas as colunas obrigatórias estão no arquivo
      missing_columns <- setdiff(required_columns, names(file_data))
      
      # Se as colunas estiverem faltando, exibe o alerta
      if (length(missing_columns) > 0) 
      {
        # Cria uma lista de colunas faltantes para exibir no warning
        missing_columns_text <- paste(missing_columns, collapse = ", ")
        
        shinyalert(
          title = "Error: Columns Missing",
          text = paste("The following columns are missing from the file:", missing_columns_text),
          type = "warning",
          confirmButtonCol = "#DD6B55",
          confirmButtonText = "OK"
        )
        return() # Interrompe a execução caso as colunas estejam faltando
      }
      
      # Atualiza o progresso para 50% (colunas verificadas)
      incProgress(0.2, detail = "Completed checking the columns...")    
      
      # Gerar gráficos após a verificação das colunas
      incProgress(0.5, detail = "Creating graphics...")
      plots <- processar_gerar_grafico(input$file1$datapath)
      
      # Armazenar os gráficos
      stored_plots(plots)
      
      # Finaliza a barra de progresso
      incProgress(0.2, detail = "Complete processing")
      
    })  
    
    output$grafico <- renderPlot({
      if (input$select == "With the alleles") {
        return(stored_plots()$plot_alleles)  # Exibe o gráfico de alelos
      } else {
        return(stored_plots()$plot_presence)  # Exibe o gráfico de presença
      }
    })
    
    
  })
  
  content = function(file) {
    plots <- stored_plots()
    if (!is.null(plots)) {
      if (input$select == "With the alleles") {
        ggsave(file, plot = plots$plot_alleles, width = 11.5, height = 7, bg = "white")
      } else {
        ggsave(file, plot = plots$plot_presence, width = 11.5, height = 7, bg = "white")
      }
    }}
  
  #----------------titulo do grafico
  output$txt <- renderText({
    req(input$txt)  # Só exibe se o campo de entrada não estiver vazio
    if (input$txt != "") {
      return(input$txt)  # Exibe o título digitado
    } else {
      return(NULL)  # Caso contrário, não exibe nada
    }
  })
}

shinyApp(ui = ui, server = server)


# 