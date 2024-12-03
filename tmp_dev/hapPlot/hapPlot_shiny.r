if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
library(pacman)


# List of packages to load
packages <- c("shiny",
              "tidyverse",
              "ggplot2",
              "gggenes",
              "viridis",
              "reshape2",
              "shinyalert")

# Function to install and load packages if missing
p_load(packages, character.only = TRUE)

# Criar um vetor de cores para cada gene
gene_colors <- c("KIR2DL1" = "#CDB5CD", "KIR2DL2" = "#ADFF2F",
                 "KIR2DL3" = "#436EEE", "KIR2DL4" = "#8B8B7A",
                 "KIR2DL5A" = "#A52A2A", "KIR2DP1" = "#C1FFC1",
                 "KIR2DS1" = "#00F5FF", "KIR2DS2" = "#EE6A50",
                 "KIR2DS3" = "#CDC673", "KIR2DS4" = "orange",
                 "KIR2DS5" = "#EEA2AD", "KIR3DL1" = "#FF83FA",
                 "KIR3DL2" = "#AB82FF", "KIR3DL3" = "#B0C4DE",
                 "KIR3DP1" = "#FFFF00", "KIR3DS1" = "#008B00",
                 "KIR2DL5B" = "#CD8C95")

color_mapping <- unlist(sapply(names(gene_colors), function(gene) {
  c(`1` = gene_colors[[gene]], `0` = "white")
}))

prepared_data <- function(input_file)
{
  table <- read.table(input_file, h = T, sep = "\t")

  # Filtrando as colunas necessárias para o app -H1 e -H2 de cada amostra

  table_filtered <- table %>%
    select(Sample, matches("^KIR.*h[12]")) %>%  # Seleciona "Sample" e colunas que começam com "KIR" e contêm "h1" e "h2"
    rename_with(~ str_replace_all(., "\\.+", "_") %>% str_replace("_$", ""), matches("^KIR.*h[12]")) %>%  # Substitui "." por "_"
    mutate(
      KIR2DL5A_h1 = ifelse(grepl("KIR2DL5A", KIR2DL5AB_h1),
                           sub(".*(KIR2DL5A\\*\\w+).*", "\\1", KIR2DL5AB_h1), #Verifica se o Allele KIR2DL5A está presente na coluna KIR2DL5AB_h1.
                           "KIR2DL5AB*null"),                                 #Extrai apenas o Allele KIR2DL5A (e.g., KIR2DL5A*001), caso contrário, KIR2DL5AB*null
      KIR2DL5B_h1 = ifelse(grepl("KIR2DL5B", KIR2DL5AB_h1),
                           sub(".*(KIR2DL5B\\*\\w+).*", "\\1", KIR2DL5AB_h1),
                           "KIR2DL5AB*null"),
      KIR2DL5A_h2 = ifelse(grepl("KIR2DL5A", KIR2DL5AB_h2),
                           sub(".*(KIR2DL5A\\*\\w+).*", "\\1", KIR2DL5AB_h2),
                           "KIR2DL5AB*null"),
      KIR2DL5B_h2 = ifelse(grepl("KIR2DL5B", KIR2DL5AB_h2),
                           sub(".*(KIR2DL5B\\*\\w+).*", "\\1", KIR2DL5AB_h2),
                           "KIR2DL5AB*null")
    ) %>%
    select(-KIR2DL5AB_h1) %>% # Remove a coluna original, opcional -eu retirei pq agr são coisas distintas
    select(-KIR2DL5AB_h2) %>%  # Remove a coluna original, opcional
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


  Alleles_long <- table_filtered %>%
    select(Sample, matches("^KIR.*_h[12]$")) %>%  # Seleciona as colunas de Alleles (h1 e h2)
    pivot_longer(cols = matches("^KIR.*_h[12]$"),  # Transforma as colunas de Alleles
                 names_to = "Name",  # Coloca o nome completo na coluna 'Name'
                 values_to = "Allele") %>%  # Coloca os valores dos Alleles na coluna 'Allele'
    separate(Name, into = c("Gene", "Haplo"), sep = "_", extra = "merge", fill = "right") %>%  # Divide o nome em Gene e Haplo
    select(Sample, Gene, Haplo, Allele) %>%  # Seleciona apenas as colunas desejadas
    mutate(Allele = str_remove(Allele, "\\(version \\d+\\)")) %>%  # Remove a versão dos Alleles
    mutate(Allele = str_remove(Allele, ".*(?=\\*)")) %>%
    mutate(Allele = str_replace(Allele, "(\\w+)\\.(\\w+)\\.(.+)", "\\1*\\2\\3")) %>% # Ajusta formato dos Alleles
    mutate(Allele = str_remove(Allele, "^\\w+")) # Remove o nome do gene no início

  combined_data <- left_join(presence_long, Alleles_long, by = c("Sample", "Gene", "Haplo"))

  return(combined_data)

}


gerar_grafico <- function(combined_data, grafico_tipo, selected_samples)

{
  plot_data_samples <- combined_data %>%
    filter(Sample %in% selected_samples) %>%
    mutate(Gene_Haplo = paste(Gene, Haplo, sep = "_")) %>%
    mutate(
      xpos = case_when(
        Gene == "KIR3DL3" ~ 1,
        Gene == "KIR2DS2" ~ 3.0,
        Gene == "KIR2DL2" ~ 5.0,
        Gene == "KIR2DL3" ~ 5.0,  # Mesmo valor que KIR2DL2
        Gene == "KIR2DL5B" ~ 7.0,
        Gene == "KIR2DS3" ~ 9.0,
        Gene == "KIR2DS5" ~ 9.0,  # Mesmo valor que KIR2DS3
        Gene == "KIR2DP1" ~ 11.0,
        Gene == "KIR2DL1" ~ 13.0,
        Gene == "KIR3DP1" ~ 15.0,
        Gene == "KIR2DL4" ~ 17.0,
        Gene == "KIR3DL1" ~ 19.0,
        Gene == "KIR3DS1" ~ 21.0,
        Gene == "KIR2DL5A" ~ 23.0,
        Gene == "KIR2DS4" ~ 25.0,
        Gene == "KIR2DS1" ~ 25.0,  # Mesmo valor que KIR2DS4
        Gene == "KIR3DL2" ~ 27.0,
        TRUE ~ NA_real_
      ),
      start = xpos,
      end = xpos + 1
    ) %>%
    mutate(Gene = factor(Gene, levels = unique(Gene)))

  #---------GRAFICOS
  if (grafico_tipo == "Only Presence") {
    plot <- plot_data_samples %>%
      filter(Presence == 1) %>%
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 27.5), #linha para conectar as caixinhas
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      facet_wrap(Sample ~ Haplo,
                 scales = "free_y",
                 ncol = 1) +
      scale_fill_manual(values = gene_colors) +
      theme_minimal() +
      guides(fill = guide_legend(nrow = 2), ncol = 3) +   #ajustar a legenda e controla a ordem da legenda - fill = guide_legend(nrow = 1), ncol = 0 ncol é a quantidade de colunas dividindo os itens da legenda, nrow são as linhas
      labs(x = "Genes",
           y = "Sample",
           title = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "top", #legenda acima do grafico
            strip.text = element_blank(),   #não aparece o nome das facetas
            legend.text = element_text(size = 18),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 20)) +
      geom_label(aes(x = -1.5,
                     label = Sample_Haplo),
                 size = 3.8,
                 fill = NA,
                 show.legend = FALSE,
                 nudge_x = 1.3)
  } else if (grafico_tipo == "With the Alleles") {
    plot <- plot_data_samples %>%
      filter(Presence == 1) %>%
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 27.5),
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      geom_text(aes(x = (start + end) / 2,  # Posiciona o texto no centro do retângulo
                    y = as.numeric(Sample),
                    label = Allele),
                size = 4,  # Ajusta o tamanho do texto
                color = "black") +  # Define a cor do texto
      facet_wrap(Sample ~ Haplo,
                 scales = "free_y",
                 ncol = 1) +
      scale_fill_manual(values = gene_colors) +
      theme_minimal() +
      guides(fill = guide_legend(nrow = 2), ncol = 2) +
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
                     label = Sample_Haplo),
                 size = 3.8,
                 fill = NA,
                 show.legend = FALSE,
                 nudge_x = 1.3)
  }

  return(plot)

}

#-----------------Telomeric e Centromeric
type_hap <- function(combined_data, haplotipo_tipo, selected_samples)
{
  plot_data_samples <- combined_data %>%
    filter(Sample %in% selected_samples) %>%
    mutate(Gene_Haplo = paste(Gene, Haplo, sep = "_")) %>%
    mutate(
      xpos = case_when(
        Gene == "KIR3DL3" ~ 1,
        Gene == "KIR2DS2" ~ 3.0,
        Gene == "KIR2DL2" ~ 5.0,
        Gene == "KIR2DL3" ~ 5.0,  # Mesmo valor que KIR2DL2
        Gene == "KIR2DL5AB" ~ 7.0,
        Gene == "KIR2DS3" ~ 9.0,
        Gene == "KIR2DS5" ~ 9.0,  # Mesmo valor que KIR2DS3
        Gene == "KIR2DP1" ~ 11.0,
        Gene == "KIR2DL1" ~ 13.0,
        Gene == "KIR3DP1" ~ 15.0,
        Gene == "KIR2DL4" ~ 17.0,
        Gene == "KIR3DL1" ~ 19.0,
        Gene == "KIR3DS1" ~ 21.0,
        Gene == "KIR2DL5AB" ~ 23.0,  # KIR2DL5AB repetido após KIR3DS1
        Gene == "KIR2DS4" ~ 25.0,
        Gene == "KIR2DS1" ~ 25.0,  # Mesmo valor que KIR2DS4
        Gene == "KIR3DL2" ~ 27.0,
        TRUE ~ NA_real_
      ),
      start = xpos,
      end = xpos + 1
    ) %>%
    mutate(Gene = factor(Gene, levels = unique(Gene)))

  #GRAFICO
  if (haplotipo_tipo == "All the Genes") {
    hap_plot <- plot_data_samples %>%
      filter(Presence == 1) %>%
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 27.5), #linha para conectar as caixinhas
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      facet_wrap(Sample ~ Haplo,
                 scales = "free_y",
                 ncol = 1) +
      scale_fill_manual(values = gene_colors) +
      theme_minimal() +
      guides(fill = guide_legend(nrow = 2), ncol = 3) +   #ajustar a legenda e controla a ordem da legenda - fill = guide_legend(nrow = 1), ncol = 0 ncol é a quantidade de colunas dividindo os itens da legenda, nrow são as linhas
      labs(x = "Genes",
           y = "Sample",
           title = "") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            legend.position = "top", #legenda acima do grafico
            strip.text = element_blank(),   #não aparece o nome das facetas
            legend.text = element_text(size = 18),
            plot.margin = margin(t = 5, r = 10, b = 5, l = 20)) +
      geom_label(aes(x = -1.5,
                     label = Sample_Haplo),
                 size = 3.8,
                 fill = NA,
                 show.legend = FALSE,
                 nudge_x = 1.3)
  } else if (haplotipo_tipo == "Only Centromeric KIR Genes") {
    hap_plot <- plot_data_samples %>%
      filter(Presence == 1,
             Gene %in% c("KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5AB",
                         "KIR2DS3", "KIR2DS5", "KIR2DP1", "KIR2DL1", "KIR3DP1")) %>% #filtra só os genes centromericos
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 27.5), #linha para conectar as caixinhas
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      geom_text(aes(x = (start + end) / 2,  # Posiciona o texto no centro do retângulo
                    y = as.numeric(Sample),
                    label = Allele),
                size = 4,  # Ajusta o tamanho do texto
                color = "black") +  # Define a cor do texto
      facet_wrap(Sample ~ Haplo,
                 scales = "free_y",
                 ncol = 1) +
      scale_fill_manual(values = gene_colors) +
      theme_minimal() +
      guides(fill = guide_legend(nrow = 2), ncol = 2) +
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
                     label = Sample_Haplo),
                 size = 3.8,
                 fill = NA,
                 show.legend = FALSE,
                 nudge_x = 1.3)
  } else if (haplotipo_tipo == "Only Telomeric KIR Genes") {
    hap_plot <- plot_data_samples %>%
      filter(Presence == 1,
             Gene %in% c("KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5AB", "KIR2DS4",
                         "KIR2DS1", "KIR3DL2")) %>%  #filtra só os genes centromericos
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 27.5), #linha para conectar as caixinhas
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      geom_text(aes(x = (start + end) / 2,  # Posiciona o texto no centro do retângulo
                    y = as.numeric(Sample),
                    label = Allele),
                size = 4,  # Ajusta o tamanho do texto
                color = "black") +  # Define a cor do texto
      facet_wrap(Sample ~ Haplo,
                 scales = "free_y",
                 ncol = 1) +
      scale_fill_manual(values = gene_colors) +
      theme_minimal() +
      guides(fill = guide_legend(nrow = 2), ncol = 2) +
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
                     label = Sample_Haplo),
                 size = 3.8,
                 fill = NA,
                 show.legend = FALSE,
                 nudge_x = 1.3)
  }
  return(hap_plot)
}

#-----------------------------app
ui <- fluidPage(
  tags$div(
    style = "display: flex; align-items: center; margin-bottom: 10px;",
    tags$h1("KIR Haplotypes App", style = "margin-right: 15px;"),
    tags$img(src = "https://www.fmb.unesp.br/Home/Pesquisa/unidadedepesquisaexperimental/dedicados/lgmb-no-site-771x460.png",
             height = "75px")
  ),
  sidebarLayout(
    sidebarPanel( #criar o painel do lado esquerdo
      width = 3,  # Ajusta o tamanho do painel lateral (3 é a largura em colunas do layout Bootstrap)
      fileInput("file1", "Upload a data file", accept = c(".txt", ".csv"), width = "100%"),  # Largura do painel

      tags$hr(), #add uma linha horizontal

      selectizeInput("amostras_input", "Select the Samples", choices = NULL, multiple = TRUE, options = list(placeholder = "The first 5 Samples of the input")),

      radioButtons("select", "Select which information you want:",  #montar a caixinha para selecionar
                   choices = c("Only Presence", "With the Alleles"),
                   selected = "Only Presence"), # "Plot Presence" é a opção selecionada por padrão
      tags$hr(), #add uma linha horizontal

      # radioButtons("type", "Type of Haplotype",  #montar a caixinha para selecionar
      #              choices = c("All the Genes", "Only Centromeric KIR Genes", "Only Telomeric KIR Genes"),
      #              selected = "All the Genes"), # "Plot Presence" é a opção selecionada por padrão

      tags$hr(), #add uma linha horizontal

      actionButton("generate", "Generate Graph", style = "width: 100%;"),

      tags$hr(), #add uma linha horizontal

      textInput(inputId = "txt", label = "Define the title of the haplotype:", value = "KIR haplotypes (h1 and h2) for each sample"), #caixinha para o usr escrever o titulo

      tags$hr(),

      downloadButton("savePlot", "Save Current Plot", style = "width: 100%;")
    ),

    mainPanel( #painel principal
      tabsetPanel(
        tabPanel("Graph", # Aba do gráfico
                 h3(textOutput("txt", container = span)),  # Título dinâmico)
                 plotOutput("grafico", height = "800px", width = "100%")),  # Gráfico
        tabPanel("About",  # Aba para "About"
                 h3("About This Application"),
                 p("This tab provides information about the application.")))

    )
  ))

server <- function(input, output, session) {

  # Armazenar os gráficos gerados
  plots <- reactiveValues(plot_Alleles = NULL, plot_presence = NULL)
  # hap_plot <- reactiveValues(plot_all = NULL, plot_cen = NULL, plot_tel = NULL)

  # Função auxiliar para verificar colunas e exibir alertas
  verificar_colunas <- function(file_data) {
    kir_columns <- grep("^KIR\\d+[A-Z]*\\.\\.h[12]\\.$", names(file_data), value = TRUE)
    required_columns <- c("Sample", kir_columns)
    missing_columns <- setdiff(required_columns, names(file_data))

    if (length(missing_columns) > 0) {
      shinyalert(
        title = "Error: Columns Missing",
        text = paste("The following columns are missing from the file:", paste(missing_columns, collapse = ", ")),
        type = "warning",
        confirmButtonCol = "#DD6B55",
        confirmButtonText = "OK"
      )
      return(NULL)
    }
    return(required_columns)
  }

  # Atualizar opções de amostras após carregar o arquivo
  observeEvent(input$file1, {
    req(input$file1)
    file_data <- read.table(input$file1$datapath, header = TRUE, sep = "\t", check.names = TRUE)

    required_columns <- verificar_colunas(file_data)
    if (is.null(required_columns)) return() # Interrompe se faltarem colunas

    # Atualizar as opções de seleção de amostras
    updateSelectizeInput(session, "amostras_input", choices = unique(file_data$Sample), server = TRUE)
  })

  # Gerar gráficos após clicar no botão "generate"
  observeEvent(input$generate, {
    req(input$file1)

    withProgress(message = 'Processing...', value = 0, {
      incProgress(0.2, detail = "Reading the file...")
      file_data <- read.table(input$file1$datapath, header = TRUE, sep = "\t", check.names = TRUE)

      required_columns <- verificar_colunas(file_data)
      if (is.null(required_columns)) return() # Interrompe se faltarem colunas

      incProgress(0.2, detail = "Completed checking the columns...")

      # Selecionar amostras ou pegar as 5 primeiras por padrão
      selected_samples <- if (length(input$amostras_input) > 0) {
        input$amostras_input
      } else {
        head(unique(file_data$Sample), 5)
      }

      incProgress(0.4, detail = "Creating graphics...")
      data <- prepared_data(input$file1$datapath)

      # Gerar gráficos
      plots$plot_Alleles <- gerar_grafico(data, "With the Alleles", selected_samples)
      plots$plot_presence <- gerar_grafico(data, "Only Presence", selected_samples)
#
#       hap_plot$plot_all <- type_hap(data, "All the Genes", selected_samples)
#       hap_plot$plot_cen <- type_hap(data, "Only Centromeric KIR Genes", selected_samples)
#       hap_plot$plot_tel <- type_hap(data, "Only Telomeric KIR Genes", selected_samples)

      incProgress(0.2, detail = "Done!")

      # Renderizar gráfico selecionado
      output$grafico <- renderPlot({
        if (input$select == "With the Alleles") {
          plots$plot_Alleles
        } else {
          plots$plot_presence
        }
      })

})
  })

  # Botão para salvar o gráfico exibido atualmente
  output$savePlot <- downloadHandler(
    filename = function() {
      if (input$select == "With the Alleles") {
        return("Alleles_plot.png")
      } else {
        return("presence_plot.png")
      }
    },
    content = function(file) {
      selected_plot <- if (input$select == "With the Alleles") {
        plots$plot_Alleles
      } else {
        plots$plot_presence
      }
      ggsave(file, plot = selected_plot, width = 11.5, height = 7, bg = "white")
    }
  )

  # Renderizar o título do gráfico
  output$txt <- renderText({
    if (!is.null(input$txt) && input$txt != "") {
      return(input$txt)
    } else {
      return(NULL)
    }
  })
}

shinyApp(ui = ui, server = server)

