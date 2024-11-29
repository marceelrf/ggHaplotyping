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

# Criar um vetor de cores para cada gene
gene_colors <- c("KIR2DL1" = "#CDB5CD", "KIR2DL2" = "#ADFF2F",
                 "KIR2DL3" = "#436EEE", "KIR2DL4" = "#8B8B7A",
                 "KIR2DL5AB" = "#A52A2A", "KIR2DP1" = "#C1FFC1",
                 "KIR2DS1" = "#00F5FF", "KIR2DS2" = "#EE6A50",
                 "KIR2DS3" = "#CDC673", "KIR2DS4" = "orange",
                 "KIR2DS5" = "#EEA2AD", "KIR3DL1" = "#FF83FA",
                 "KIR3DL2" = "#AB82FF", "KIR3DL3" = "#B0C4DE",
                 "KIR3DP1" = "#FFFF00", "KIR3DS1" = "#008B00")

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
      geom_segment(aes(x = 0.9, xend = 18.5),
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
  } else if (grafico_tipo == "With the alleles") {
    plot <- plot_data_samples %>%
      filter(Presence == 1) %>%
      mutate(Haplo = ifelse(grepl("h1", Gene_Haplo),
                            "H1", "H2"),
             Sample_Haplo = paste(Sample, Haplo, sep = " - ")) %>%
      ggplot(aes(xmin = start,
                 xmax = end,
                 y = Sample,
                 fill = Gene)) +
      geom_segment(aes(x = 0.9, xend = 18.5),
                   color = "black", size = 0.5) +
      geom_rect(aes(xmin = start - 0.2,
                    xmax = end + 0.2,
                    ymin = as.numeric(Sample) - 0.01,
                    ymax = as.numeric(Sample) + 0.01),
                color = "black",
                linewidth = 0.5) +
      geom_text(aes(x = (start + end) / 2,  # Posiciona o texto no centro do retângulo
                    y = as.numeric(Sample),
                    label = Alelo),
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
                   choices = c("Only Presence", "With the alleles"),
                   selected = "Only Presence"), # "Plot Presence" é a opção selecionada por padrão
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
                 p("This tab provides information about the application."))))


    )
  )


server <- function(input, output, session) {


  # Armazenar os gráficos gerados
  plots <- reactiveValues(plot_alleles = NULL, plot_presence = NULL)


  # Ao carregar o arquivo, processe os dados e atualize as amostras disponíveis para seleção
  observeEvent(input$file1, {
    req(input$file1)
    data <- prepared_data(input$file1$datapath)
    updateSelectizeInput(session, "amostras_input", choices = unique(data$Sample), server = TRUE)
  })

  observeEvent(input$generate, {

      req(input$file1)

    # Barra de progresso
      withProgress(message = 'Processing...', value = 0,
     {
      incProgress(0.2, detail = "Reading the file...")

      # Ler e processar o arquivo
      file_data <- read.table(input$file1$datapath, header = TRUE, sep = "\t", check.names = TRUE)

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

      # Obter as amostras selecionadas
      selected_samples <- input$amostras_input
      if (length(selected_samples) == 0) {
        data <- prepared_data(input$file1$datapath)
        selected_samples <- head(unique(data$Sample), 5) # Pega as 5 primeiras amostras, caso nenhuma seja selecionada
      }

      incProgress(0.6, detail = "Creating graphics...")

      data <- prepared_data(input$file1$datapath)

      incProgress(0.2, detail = "Done!")

      # Gerar o gráfico com as amostras selecionadas
      output$grafico <- renderPlot({
                         gerar_grafico(data, input$select, selected_samples)})


      # Gerar os gráficos e armazená-los
      plots$plot_alleles <- gerar_grafico(data, "With the alleles", selected_samples)
      plots$plot_presence <- gerar_grafico(data, "Only Presence", selected_samples)

    })


      # Botão para salvar o gráfico exibido atualmente
      output$savePlot <- downloadHandler(
        filename = function() {
          if (input$select == "With the alleles") {
            return("alleles_plot.png")
          } else {
            "presence_plot.png"
          }
        },
        content = function(file) {
          if (input$select == "With the alleles") {
              ggsave(file, plot = plots$plot_alleles, width = 11.5, height = 7, bg = "white")
            } else {
              ggsave(file, plot = plots$plot_presence, width = 11.5, height = 7, bg = "white")
            }})


  # Renderizar o título do gráfico
  output$txt <- renderText({
    req(input$txt)  # Exibe apenas se o campo de entrada não estiver vazio
    if (input$txt != "") {
      return(input$txt)
    } else {
      return(NULL)}})

})
}



shinyApp(ui = ui, server = server)

