#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinysky)
library(shinydashboard)
library(tidyverse)
library(ggtranscript)
library(ggpubr)
library(scales)
library(plotly)
library(DT)

# Create the header object
header <- dashboardHeader(
  title = "Ebbert Lab"
)

# Create the sidebar with different tabs
sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "tabs",
    menuItem("Transcript Expression", tabName = "TEx", icon = icon("stats", lib = "glyphicon")),
    menuItem("New Gene Bodies", tabName = "ngb", icon = icon("plus", lib = "glyphicon")),
    menuItem("New Transcripts", tabName = "ntx", icon = icon("plus", lib = "glyphicon")),
    menuItem("Ebbert Lab Website", href = "https://ebbertlab.com/", icon = icon("new-window", lib = "glyphicon"))
  )
)

# load data for the app
# transcript data, formatted for plot
transcript_data <- read_tsv("../../../../../mlpa241/Downloads/cDNA_files_r_shiny/annotation_r_shiny.tsv") %>% 
  rename(seqnames = chr) %>%
  mutate(gene_name = coalesce(gene_name, gene_id))

seqnames_levels = c("1", "2", "3", "4", "5", 
                    "6", "7", "8", "9", "10", 
                    "11", "12", "13", "14", "15", 
                    "16", "17", "18", "19", "20", 
                    "21", "22", "MT", "X", "Y", 
                    "ERCC-00002", "ERCC-00003", "ERCC-00004", 
                    "ERCC-00009", "ERCC-00012", "ERCC-00013", "ERCC-00014", "ERCC-00016", 
                    "ERCC-00017", "ERCC-00019", "ERCC-00022", "ERCC-00024", "ERCC-00025", 
                    "ERCC-00028", "ERCC-00031", "ERCC-00033", "ERCC-00034", "ERCC-00035", 
                    "ERCC-00039", "ERCC-00040", "ERCC-00041", "ERCC-00042", "ERCC-00043",
                    "ERCC-00044", "ERCC-00046", "ERCC-00048", "ERCC-00051", "ERCC-00053", 
                    "ERCC-00054", "ERCC-00057", "ERCC-00058", "ERCC-00059", "ERCC-00060",
                    "ERCC-00061", "ERCC-00062", "ERCC-00067", "ERCC-00069", "ERCC-00071", 
                    "ERCC-00073", "ERCC-00074", "ERCC-00075", "ERCC-00076", "ERCC-00077", 
                    "ERCC-00078", "ERCC-00079", "ERCC-00081", "ERCC-00083", "ERCC-00084", 
                    "ERCC-00085", "ERCC-00086", "ERCC-00092", "ERCC-00095", "ERCC-00096", 
                    "ERCC-00097", "ERCC-00098", "ERCC-00099", "ERCC-00104", "ERCC-00108", 
                    "ERCC-00109", "ERCC-00111", "ERCC-00112", "ERCC-00113", "ERCC-00116",
                    "ERCC-00117", "ERCC-00120", "ERCC-00123", "ERCC-00126", "ERCC-00130", 
                    "ERCC-00131", "ERCC-00134", "ERCC-00136", "ERCC-00137", "ERCC-00138", 
                    "ERCC-00142", "ERCC-00143", "ERCC-00144", "ERCC-00145", "ERCC-00147", 
                    "ERCC-00148", "ERCC-00150", "ERCC-00154", "ERCC-00156", "ERCC-00157", 
                    "ERCC-00158", "ERCC-00160", "ERCC-00162", "ERCC-00163", "ERCC-00164", 
                    "ERCC-00165", "ERCC-00168", "ERCC-00170", "ERCC-00171", "GL000009.2", 
                    "GL000194.1", "GL000195.1", "GL000205.2", "GL000213.1", "GL000216.2", 
                    "GL000218.1", "GL000219.1", "GL000220.1", "GL000224.1", "GL000225.1", 
                    "KI270442.1", "KI270711.1", "KI270713.1", "KI270721.1", "KI270726.1", 
                    "KI270727.1", "KI270728.1", "KI270731.1", "KI270733.1", "KI270734.1", 
                    "KI270744.1", "KI270750.1")

transcript_data$seqnames <- fct_rev(factor(transcript_data$seqnames, levels = seqnames_levels))

new_genes <- transcript_data %>%
  filter(type == "transcript", startsWith(gene_id, 'gene.')) %>%
  arrange(seqnames)

new_transcripts <- transcript_data %>%
  filter(type == "transcript", startsWith(gene_id, 'E'), startsWith(transcript_id, 'tx.')) %>%
  arrange(seqnames)

# expression data, formatted for plot
expression_data <- read_tsv("../../../../../mlpa241/Downloads/cDNA_files_r_shiny/expression_matrix_r_shiny.tsv") %>%
  left_join(transcript_data %>% select(transcript_id, annotation_status, discovery_category, transcript_biotype) %>% distinct(), by = "transcript_id") %>%
  pivot_longer(cols=3:10, names_to = "sample_expression_type", values_to = "number") %>%
  extract(sample_expression_type, into = c("sample_name", "expression_type"), "(.*)_([^_]+$)") %>%
  pivot_wider(names_from = expression_type, values_from = number)

# Density plot data
density_data <- read_tsv("../../../../../mlpa241/Downloads/cDNA_files_r_shiny/expression_matrix_r_shiny_GENE.tsv")


anti_sense_data <- read_tsv("../../../../../mlpa241/Downloads/cDNA_files_r_shiny/AllNovelTranscripts_overlap_All_HG38_v107_Transcripts.tsv") %>% 
  filter(StrandType == "OppositeStrand", startsWith(NewGeneID, "gene"), startsWith(Hg38v107_GeneID, "E")) %>% 
  select(!Hg38v107_transcriptID) %>% 
  distinct()

# contains just gene_id and gene_names
just_genes <- transcript_data %>% select(gene_id, gene_name) %>% distinct()

# table of new genes that could potentially be anti-sense
is_antisense <- left_join(anti_sense_data %>% 
                      select(NewGeneID, Hg38v107_GeneID), 
                      just_genes, 
                      by = c("NewGeneID" = "gene_id")) %>%
       rename(gene_id = NewGeneID) %>%
       rename(antisense_id = Hg38v107_GeneID) #%>%
       # add_column(is_antisense = "true") %>%
       # add_column(has_antisense = "false")

# table of know genes that we potentially found an anti-sense transcript for
has_antisense <- left_join(anti_sense_data %>%
                      select(Hg38v107_GeneID, NewGeneID),
                      just_genes, 
                      by = c("Hg38v107_GeneID" = "gene_id")) %>%
       rename(gene_id = Hg38v107_GeneID) %>%
       rename(antisense_id = NewGeneID) #%>%
       # add_column(has_antisense = "true") %>%
       # add_column(is_antisense = "false")

# combine into single table for lookup
just_genes_antisense <- bind_rows(is_antisense, has_antisense)

# gene search space
lookup <- as.data.frame(
             transcript_data %>%
               select(gene_id, gene_name) %>%
               distinct(gene_id, gene_name)
         )

template <- HTML("<strong><em>{{gene_name}}</em></strong> <p> <font size=\"2\"><em>{{gene_id}}</em></font> ")

# Create the body of the app
body <- dashboardBody(
  tabItems(
    # the first sidebar tab, displaying the transcript expression
    tabItem(
      tabName = "TEx",
      fluidRow(
        h1("Page Title", align='center'),
        p("title of publication and reference", align='center'),
        column(
          width = 12,
          box(
            status = "primary",
            uiOutput("plot.ui"),
            width = NULL,
            div(
              align='right',
              actionButton('linkButton', "Checkout new gene bodies we discovered"),
              downloadButton('downloadFig', 'Download', width=NULL)
            )
          )
          
        ),
      ),
      fluidRow(
        # search bar for genes
        column(
          width = 6,
          box(
            width = NULL,
            status = "primary",
            helpText("Search for gene of interest:"),
            textInput.typeahead(
              id="geneSearch",
              placeholder="HUGO Ensembl ID",
              local=lookup,
              valueKey = "gene_id",
              tokens=c(1,2),
              template= HTML("")
            )
          ),
          box(
            width = NULL,
            status = "warning",
            helpText(h4("Use the below options to customize the figures",align='center'))
          ),
          column(
            width = 6,
            box(
              # option for number of transcripts to display
              width = NULL,
              solidHeader = TRUE,
              radioButtons("numTransRadio",
                           label = "Transcripts to Display",
                           choices = c("Top 5" = "top_5", "Top 10" = "top_10", "Top 50%" = "top_half", "All" = "all"),
                           selected = "all"),
            ),
            box(
              # option for how to color the figure
              width = NULL,
              solidHeader = TRUE,
              radioButtons("colorRadio", label = "Color by", choices = c("Discovery" = "annotation_status", "Type of new transcript" = "discovery_category", "Transcript Biotype" = "transcript_biotype"))
            ),
          ),
          column(
            width = 6,
            box(
              # option for how to display expression (CPM vs Total counts)
              width = NULL,
              solidHeader = TRUE,
              radioButtons("expressionRadio", label = "Expression display", choices = c("Counts Per Million (CPM)" = "CPM", "Total counts" = "counts")),
              
            )
          ),
        ),
        #   # plot with transcripts and expression data
        #   # below page looks like it will help with the dynamic height
        #   # https://stackoverflow.com/questions/31882463/dynamic-plot-height-in-shiny
        column(
          width = 6,
          box(
            width = NULL,
            status = "primary",
            plotOutput("densityPlot"),
            div(align='right', downloadButton('downloadDensityFig', 'Download'))
          )
        ),
      ),
      fluidRow(
        # slider to determine width of download plot
        box(
          width = 6,
          solidHeader = TRUE,
          helpText(h5("Download dimensions for Transcripts and Expression plot",align='center')),
          sliderInput("widthSlider", "Width:", 5, 49, 12),
          # slider to determine height of download plot
          sliderInput("heightSlider", "Height:", 3, 49, 6)
        ),
        box(
          width = 6,
          solidHeader = TRUE,
          helpText(h5("Download dimensions for Density plot",align='center')),
          sliderInput("widthDensitySlider", "Width:", 3, 15, 7),
          # slider to determine height of download plot
          sliderInput("heightDensitySlider", "Height:", 5, 15, 8)
        ),
      ),
      fluidRow(
        box(
          width = 12,
          title = "References",
          status = "info",
          p("Chang W, Borges Ribeiro B (2021). _shinydashboard: Create Dashboards with 'Shiny'_. R package version 0.7.2, <https://CRAN.R-project.org/package=shinydashboard>."),
          p("Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2022). _shiny: Web Application Framework for R_. R package version 1.7.3, <https://CRAN.R-project.org/package=shiny>."),
          p(" ZJ D (2019). _shinysky: A Set of Shiny Components and Widgets_. R package version 0.1.3, <https://github.com/AnalytixWare/ShinySky>."),
          p("Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>."),
          p("Gustavsson EK, Zhang D, Reynolds RH, Garcia-Ruiz S, Ryten M (2022). “ggtranscript: an R package for the visualization and interpretation of transcript isoforms using ggplot2.” _Bioinformatics_. doi:10.1093/bioinformatics/btac409 <https://doi.org/10.1093/bioinformatics/btac409>, <https://academic.oup.com/bioinformatics/article/38/15/3844/6617821>."),
          p("Kassambara A (2022). _ggpubr: 'ggplot2' Based Publication Ready Plots_. R package version 0.5.0, <https://CRAN.R-project.org/package=ggpubr>.")
        )
      ),
    ),
    
    
    
    tabItem(tabName = "ngb",
      tags$head(tags$style('
        #my_tooltip {
        position: absolute;
        width: 300px;
        z-index: 100;
        padding: 0;
        }
      ')),
            
      tags$script('
        $(document).ready(function() {
          // id of the plot
          $("#testing").mousemove(function(e) { 
          
            // ID of uiOutput
            $("#my_tooltip").show();         
            $("#my_tooltip").css({             
              top: (e.pageY + 5) + "px",             
              left: (e.pageX + 5) + "px"         
            });     
          });     
        });
      '),  
            
      fluidRow(
        h1("Page Title", align='center'),
        p("title of publication and reference", align='center'),
        box(
          width = NULL,
          plotOutput("testing", 
             hover = hoverOpts(id = 'plot_hover', delay = 50), 
             brush = brushOpts(id = "plot_brush", resetOnNew = TRUE)),
          div(
            align='right',
            actionButton('resetButton', "Reset Graph"),
          )
        ),
        box(
          width = NULL,
          DT::dataTableOutput("selected")
        )
      ),
      fluidRow(
        box(
          width = 12,
          title = "References",
          status = "info",
          p("Chang W, Borges Ribeiro B (2021). _shinydashboard: Create Dashboards with 'Shiny'_. R package version 0.7.2, <https://CRAN.R-project.org/package=shinydashboard>."),
          p("Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2022). _shiny: Web Application Framework for R_. R package version 1.7.3, <https://CRAN.R-project.org/package=shiny>."),
          p(" ZJ D (2019). _shinysky: A Set of Shiny Components and Widgets_. R package version 0.1.3, <https://github.com/AnalytixWare/ShinySky>."),
          p("Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>."),
          p("Gustavsson EK, Zhang D, Reynolds RH, Garcia-Ruiz S, Ryten M (2022). “ggtranscript: an R package for the visualization and interpretation of transcript isoforms using ggplot2.” _Bioinformatics_. doi:10.1093/bioinformatics/btac409 <https://doi.org/10.1093/bioinformatics/btac409>, <https://academic.oup.com/bioinformatics/article/38/15/3844/6617821>."),
          p("Kassambara A (2022). _ggpubr: 'ggplot2' Based Publication Ready Plots_. R package version 0.5.0, <https://CRAN.R-project.org/package=ggpubr>.")
        )
      ),
      uiOutput("my_tooltip")
    ),
    
    
    tabItem(tabName = "ntx",
            tags$head(tags$style('
        #my_tooltip_tx {
        position: absolute;
        width: 300px;
        z-index: 100;
        padding: 0;
        }
      ')),
            
            tags$script('
        $(document).ready(function() {
          // id of the plot
          $("#new_transcripts_plot").mousemove(function(e) { 
          
            // ID of uiOutput
            $("#my_tooltip_tx").show();         
            $("#my_tooltip_tx").css({             
              top: (e.pageY + 5) + "px",             
              left: (e.pageX + 5) + "px"         
            });     
          });     
        });
      '),  
            
            fluidRow(
              h1("Page Title", align='center'),
              p("title of publication and reference", align='center'),
              box(
                width = NULL,
                plotOutput("new_transcripts_plot", 
                           hover = hoverOpts(id = 'plot_hover_tx', delay = 50), 
                           brush = brushOpts(id = "plot_brush_tx", resetOnNew = TRUE)),
                div(
                  align='right',
                  actionButton('resetButton_tx', "Reset Graph"),
                )
              ),
              box(
                width = NULL,
                DT::dataTableOutput("selected_transcripts")
              )
            ),
            fluidRow(
              box(
                width = 12,
                title = "References",
                status = "info",
                p("Chang W, Borges Ribeiro B (2021). _shinydashboard: Create Dashboards with 'Shiny'_. R package version 0.7.2, <https://CRAN.R-project.org/package=shinydashboard>."),
                p("Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2022). _shiny: Web Application Framework for R_. R package version 1.7.3, <https://CRAN.R-project.org/package=shiny>."),
                p(" ZJ D (2019). _shinysky: A Set of Shiny Components and Widgets_. R package version 0.1.3, <https://github.com/AnalytixWare/ShinySky>."),
                p("Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>."),
                p("Gustavsson EK, Zhang D, Reynolds RH, Garcia-Ruiz S, Ryten M (2022). “ggtranscript: an R package for the visualization and interpretation of transcript isoforms using ggplot2.” _Bioinformatics_. doi:10.1093/bioinformatics/btac409 <https://doi.org/10.1093/bioinformatics/btac409>, <https://academic.oup.com/bioinformatics/article/38/15/3844/6617821>."),
                p("Kassambara A (2022). _ggpubr: 'ggplot2' Based Publication Ready Plots_. R package version 0.5.0, <https://CRAN.R-project.org/package=ggpubr>.")
              )
            ),
            uiOutput("my_tooltip_tx")
    )
  ),
)


# Pass everything to the UI
ui <- dashboardPage(
  header,
  sidebar,
  body
)



# Define server logic 
server <- function(input, output, session) {
  # print(citation("shinydashboard"))
  # print(citation("shiny"))
  # print(citation("shinysky"))
  # print(citation("tidyverse"))
  # print(citation("ggtranscript"))
  # print(citation("ggpubr"))
  
  # if 0 (default), the no antisense transcript is involved with current gene,
  # if 1, then the gene has a potential antisense transcript,
  # if 2, the new gene IS a potential antisense transcript to know gene
  # if 3, reset to no anti-sense stuff
  # switch_for_antisense <<- 0
  # display_antisense <- reactiveVal(0)
  # anti_sense_gene_id <- reactiveVal("")
  current_ids <- reactiveValues(gene_id = "ENSG00000166295", # ANAPC16
                                 anti_sense_id = "", 
                                 showing_antisense = FALSE)
  
  #id <- "ENSG00000130203" # APOE
  #id <- "ENSG00000166295" # ANAPC16
  #id <- "ENSG00000168209" # DDIT4
  #id <- "ENSG00000219545" # UMAD1
  
  # the base data to be displayed
  plot_data <- reactive({
    id <- current_ids$gene_id
    antisense_id <- current_ids$anti_sense_id

    # will need to arrange in order
    gene_transcripts <- transcript_data %>% 
      filter(gene_id == id)
    gene_expression <- expression_data %>% 
      filter(gene_id == id)
    
    # #grab the gene_name to print out later
    selected_gene_name <- unique(gene_transcripts$gene_name)
    
    #grab the gene_name to print out later
    antisense_strand <- ""
    antisense_transcripts <- c()

    if (current_ids$showing_antisense == TRUE){
      gene_antisense_transcripts <- transcript_data %>%
        filter(gene_id == antisense_id)
      gene_antisense_expression <- expression_data %>%
        filter(gene_id == antisense_id)
      antisense_strand <- unique(gene_antisense_transcripts$strand)
      antisense_transcripts <- unique(gene_antisense_transcripts$transcript_id)
      
      gene_antisense_transcripts <- gene_antisense_transcripts %>%
            mutate(strand = unique(gene_transcripts$strand))
      
      gene_transcripts <- bind_rows(gene_transcripts, gene_antisense_transcripts)
      gene_expression <- bind_rows(gene_expression, gene_antisense_expression)
    }
    
    return(
      list(
        txs = gene_transcripts,
        expr = gene_expression,
        antisense_strand = antisense_strand,
        antisense_txs = antisense_transcripts,
        gene_name = selected_gene_name
      )
    )
  })
  
  
  # Change the data in the plot data to reflect if the user wants to see antisense or not
  observeEvent(input$linkButton, {
    print(current_ids$gene_id)
    if(current_ids$showing_antisense == FALSE & current_ids$anti_sense_id != "") {
      current_ids$showing_antisense = TRUE
      updateLinkButton('Reset')
    } else if (current_ids$anti_sense_id != ""){
      current_ids$showing_antisense = FALSE
      if (startsWith(current_ids$gene_id, "E")) {
        updateLinkButton('PS. Checkout the potential anti-sense transcript(s) we found for this gene')
      } else {
        updateLinkButton('PS. This gene body has potential anti-sense transcript(s) for a known gene')
      }
    } else {
      updateTabItems(session, "tabs", "ngb")
    }
  })
  
  
  #function to update the linkButton text
  updateLinkButton <- function(btn_text){
    updateActionButton(
        session,
        "linkButton",
        label = btn_text
      )
  }
  
  
  #Update linkButton text on gene change
  observeEvent(input$geneSearch, {
    req(input$geneSearch)
    current_ids$gene_id <- input$geneSearch
    
    if (nrow(tib <- just_genes_antisense %>% filter(gene_id == input$geneSearch)) > 0) {
      current_ids$anti_sense_id <- tib$antisense_id
      if (startsWith(current_ids$gene_id, "E")) {
        updateLinkButton('PS. Checkout the potential anti-sense transcript(s) we found for this gene')
        # if it is a new gene body
      } else {
        updateLinkButton('PS. This gene body has potential anti-sense transcript(s) for a known gene')
      }
    } else {
      current_ids$anti_sense_id <- ""
      updateLinkButton('Checkout new gene bodies we discovered')
    }
    
    current_ids$showing_antisense = FALSE
  })
  
  
  selected_new_genes <- reactiveValues(data = new_genes)
  selected_new_transcripts <- reactiveValues(data = new_transcripts)

  
  # Used for the gene search bar
  updateTextInput.typeahead(session, "geneSearch", lookup, "gene_id", c(lookup$gene_name, lookup$gene_id), template, 
                            placeholder = "HUGO Ensembl ID")
  
  
  # Use the brush to zoom in on the new gene bodies graph\
  # TODO: need to fix the hover tool tip to allow hovering on more than just the exact y line and X start
  observeEvent(input$plot_brush, {
    brushed <- brushedPoints(new_genes, input$plot_brush, xvar = "start")
    if (nrow(brushed) > 0){
      selected_new_genes$data <- brushed
    }
  })
  
  
  # Reset the new gene bodies graph on click
  observeEvent(input$resetButton, {
    selected_new_genes$data <- new_genes
  })
  
  
  # render the new gene bodies plot
  output$testing <- renderPlot({
      plot_new_genes <- selected_new_genes$data
      #plot transcripts
      ggplot(plot_new_genes, aes(
        xstart = start,
        xend = end,
        y = seqnames
      )) +
      geom_range(
        aes(fill = !! sym(input$colorRadio)),
      ) +
      scale_x_continuous(labels = comma) + 
      theme(legend.position="none")
      #   theme (
      #     axis.title.y = element_blank()
      #   )
  })
  
  
  output$my_tooltip <- renderUI({
    hover <- input$plot_hover 
    y <- nearPoints(selected_new_genes$data, input$plot_hover, xvar = "start", maxpoints = 1)
    req(nrow(y) != 0)
    verbatimTextOutput("vals")
  })
  
  
  output$vals <- renderPrint({
    hover <- input$plot_hover
    y <- nearPoints(selected_new_genes$data, input$plot_hover, xvar = "start", maxpoints = 1)
    req(nrow(y) != 0)
    print(paste0( "Location: ", y$seqnames, ":", y$start, "-", y$end))
  })  
  
  
  output$selected <- DT::renderDataTable({
    res <- brushedPoints(selected_new_genes$data, input$plot_brush, xvar = "start")
    
    datatable(
      selected_new_genes$data %>% 
        arrange(desc(seqnames)) %>% 
        select(seqnames, start, end, strand, gene_id, transcript_id)
      )
  })
  
  
  # Use the brush to zoom in on the new gene bodies graph\
  # TODO: need to fix the hover tool tip to allow hovering on more than just the exact y line and X start
  observeEvent(input$plot_brush_tx, {
    brushed <- brushedPoints(new_transcripts, input$plot_brush_tx, xvar = "start")
    if (nrow(brushed) > 0){
      selected_new_transcripts$data <- brushed
    }
  })
  
  
  # Reset the new transcripts graph on click
  observeEvent(input$resetButton_tx, {
    selected_new_transcripts$data <- new_transcripts
  })
  
  
  # render the new transcripts plot
  output$new_transcripts_plot <- renderPlot({
    plot_new_transcripts <- selected_new_transcripts$data
    #plot transcripts
    ggplot(plot_new_transcripts, aes(
      xstart = start,
      xend = end,
      y = seqnames
    )) +
      geom_range(
        aes(fill = discovery_category),
      ) +
      scale_x_continuous(labels = comma) + 
      theme(legend.position="bottom")
    #   theme (
    #     axis.title.y = element_blank()
    #   )
  })
  
  
  output$my_tooltip_tx <- renderUI({
    hover <- input$plot_hover_tx 
    y <- nearPoints(transcript_data, input$plot_hover_tx, xvar = "start", maxpoints = 1)
    req(nrow(y) != 0)
    verbatimTextOutput("vals_tx")
  })
  
  
  output$vals_tx <- renderPrint({
    hover <- input$plot_hover_tx
    y <- nearPoints(transcript_data, input$plot_hover_tx, xvar = "start", maxpoints = 1)
    req(nrow(y) != 0)
    print(paste0( "Gene_ID: ", y$gene_id))
  })  
  
  
  output$selected_transcripts <- DT::renderDataTable({
    res <- brushedPoints(selected_new_transcripts$data, input$plot_brush_tx, xvar = "start")
    
    datatable(
      selected_new_transcripts$data %>% 
        arrange(desc(seqnames)) %>% 
        select(seqnames, start, end, strand, gene_id, transcript_id)
    )
  })
  

  # create the transcript/expression plot
  # this will be used to render the plot on the page and also to download
  comb_plot <- reactive({
    id <- current_ids$gene_id
    gene_transcripts <- plot_data()$txs
    gene_expression <- plot_data()$expr
    selected_gene_name <- plot_data()$gene_name
    antisense_strand <- plot_data()$antisense_strand
    antisense_transcripts <- plot_data()$antisense_txs

    factor_order <- gene_expression %>% 
      select(transcript_id, gene_id, annotation_status, discovery_category, sample_name, input$expressionRadio) %>%
      group_by(transcript_id) %>%
      mutate(exp_avg = mean(!! sym(input$expressionRadio))) %>% 
      select(transcript_id, exp_avg) %>%
      distinct() %>%
      arrange(exp_avg) %>%
      select(transcript_id)
    
    factor_order <- factor_order$transcript_id
    
    gene_expression$transcript_id <- factor(gene_expression$transcript_id, levels = factor_order)
    gene_transcripts$transcript_id <- factor(gene_transcripts$transcript_id, levels = factor_order)
    
    gene_transcripts <- gene_transcripts %>%
      arrange(transcript_id)
    
    # filter number of transcripts
    num_transcripts = length(factor_order)
    transcripts_to_display = factor_order
    
    if (input$numTransRadio == 'top_5'){
      if (num_transcripts >= 5) {
        transcripts_to_display <- tail(factor_order, 5)
      }
    } else if (input$numTransRadio == 'top_10') {
      if (num_transcripts >= 10) {
        transcripts_to_display <- tail(factor_order, 10)
      }
    } else if (input$numTransRadio == 'top_half') {
      transcripts_to_display <- tail(factor_order, ceiling(num_transcripts*.5))
    }
    
    plotHeight <- 200 + 50*length(transcripts_to_display)
    
    gene_transcripts <- filter(gene_transcripts, transcript_id %in% transcripts_to_display)
    gene_expression <- filter(gene_expression, transcript_id %in% transcripts_to_display)

    # format transcripts for ggplot
    gene_exons <- gene_transcripts %>% filter(type=="exon")
    gene_cds_regions <- gene_transcripts %>% filter(type == "CDS")
    
    # find the difference between the start and ending positions of the exons and cds.
    # This will be used to calculate the re-scaled cds regions later
    cds_exon_diff <- left_join(
      gene_cds_regions %>%
        rename(c_start = start) %>%
        rename(c_end = end) %>%
        select(!type),
      gene_exons %>%
        rename(e_start = start) %>%
        rename(e_end = end) %>%
        select(!type)
    ) %>%
      mutate(d_start = abs(e_start - c_start)) %>%
      mutate(d_end = abs(e_end - c_end))
    
    # re-scale the exons and introns
    gene_rescaled <- shorten_gaps(
      gene_exons, 
      to_intron(gene_exons, "transcript_id"), 
      group_var = "transcript_id"
    )
    
    if(current_ids$showing_antisense == TRUE){
      # fix the stand of the antisense data
      gene_rescaled <- gene_rescaled %>%
        mutate(strand = replace(strand, transcript_id %in% antisense_transcripts, antisense_strand))
    }
    
    gene_rescaled_exons <- gene_rescaled %>% filter(type == "exon")
    gene_rescaled_introns <- gene_rescaled %>% filter(type == 'intron')
    gene_rescaled_cds <- left_join(
      cds_exon_diff %>%
        mutate(type = "CDS") %>%
        select(!c(c_start, c_end, e_start, e_end)),
      gene_rescaled_exons %>%
        rename(e_start = start)%>%
        rename(e_end = end) %>%
        select(!type)
    ) %>%
      mutate(start = e_start+d_start) %>%
      mutate(end = e_end - d_end) %>%
      select(!c(e_start, e_end, d_start, d_end))
    
    # 
    if(input$colorRadio =="annotation_status") {
      fill_legend_name <- "Annotation Status"
    } else if (input$colorRadio == "transcript_biotype") {
      fill_legend_name <- "Transcript Biotype"
    } else {
      fill_legend_name <- "Discovery Category"
    }
    
    displayCategories <- transcript_data %>%
            select(!! sym(input$colorRadio)) %>%
            distinct(!! sym(input$colorRadio)) %>%
            pull(!! sym(input$colorRadio))
    
    colorLevels <- setNames(hue_pal()(length(displayCategories)), levels(as.factor(displayCategories)))
    
    #plot transcripts
    transcriptPlt <- ggplot(gene_rescaled_exons, aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
      geom_range(
        aes(fill = !! sym(input$colorRadio)),
        height = 0.25
      ) +
      geom_range(
        aes(fill = !! sym(input$colorRadio)),
        data = gene_rescaled_cds
      ) +
      geom_intron(
        data = gene_rescaled_introns,
        arrow.min.intron.length = 75*length(transcripts_to_display), # variable based on number of transcripts we are displaying
        aes(strand = strand)
      ) + 
      theme (
        axis.title.y = element_blank()
      ) + 
      scale_fill_manual(values = colorLevels) +
      guides(fill=guide_legend(title=fill_legend_name))

    expressionPlt <- ggplot(gene_expression, aes(transcript_id, !! sym(input$expressionRadio))) + 
      geom_boxplot(
        aes(
          fill = !! sym(input$colorRadio)
          ), 
        outlier.shape = NA
        ) + 
      geom_jitter() + 
      coord_flip() + 
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank()
      ) +
      scale_fill_manual(values = colorLevels)
    
    if(input$expressionRadio =="CPM") {
      exp_type <- "CPM"
    } else {
      exp_type <- "Total Counts"
    }
    
    # return plot, gene_name, and gene_id as a list so that we can access them all later for rendering the plot
    # and for downloading it
    return(
      list(
        plot = annotate_figure(
          ggarrange(transcriptPlt, expressionPlt, ncol=2, common.legend = TRUE, legend="bottom"),
            top = text_grob(paste("\n", selected_gene_name, " (", id,"): ","Transcripts and Expression (", exp_type,")\n", sep = ""), face = "bold", size = 14)
        ),
        gene_name = selected_gene_name,
        gene_id = id,
        plotHeight = plotHeight,
        num_transcripts = length(transcripts_to_display)
      )
    )
  })
  
  exprDensityPlot <- reactive({
    
    density_data <- density_data %>%
      rename(CPM = average_CPM) %>%
      rename(counts = total_counts)%>%
      filter(!! sym(input$expressionRadio) > 0)%>%
      mutate(log_comb_exp = log10(!! sym(input$expressionRadio)))

    plt <- ggplot(density_data, aes(x=log_comb_exp)) +
      geom_density() +
      geom_vline(xintercept = density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id], colour = "#F8766D") + 
      ggtitle(paste("Gene Expression: ", comb_plot()$gene_name)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = 'bold')) +
      ylab("Density")
      
    y_ranges <- layer_scales(plt)$y$range$range
    x_ranges <- layer_scales(plt)$x$range$range
    
    if (input$expressionRadio == "CPM"){
      if (! comb_plot()$gene_id %in% density_data$gene_id) {
        plt <- plt + geom_text(aes(x_ranges[2], y_ranges[2], label = paste(comb_plot()$gene_name, 'is not \nexpressed')), nudge_x = -2, nudge_y = -0.05, colour = "#F8766D")
      } else if (density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id] > 4) {
        #need to flip which side the name is on
        plt <- plt + geom_text(aes(density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id], y_ranges[2], label = comb_plot()$gene_name), nudge_x = -1, colour = "#F8766D")
      } else {
        plt <- plt + geom_text(aes(density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id], y_ranges[2], label = comb_plot()$gene_name), nudge_x = 1, colour = "#F8766D")
      }
      plt <- plt + xlab("log10 (average CPM)")
        
    } else {
      if (! comb_plot()$gene_id %in% density_data$gene_id) {
        plt <- plt + geom_text(aes(x_ranges[2], y_ranges[2], label = paste(comb_plot()$gene_name, 'is not \nexpressed')), nudge_x = -2, nudge_y = -0.05, colour = "#F8766D")
      } else if (density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id] > 6){
        plt <- plt + geom_text(aes(density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id], y_ranges[2], label = comb_plot()$gene_name), nudge_x = -1, colour = "#F8766D")
      } else {
        plt <- plt + geom_text(aes(density_data$log_comb_exp[density_data$gene_id == comb_plot()$gene_id], y_ranges[2], label = comb_plot()$gene_name), nudge_x = 1, colour = "#F8766D")
      }
      plt <- plt + xlab("log10 (combined total counts)")
    }
    
    print(plt)
  })
  
  # Set the recommended size for downloading the plot
  observe({
    heightValue <- 2 + ceiling(comb_plot()$num_transcripts * 0.5)
    widthValue <- 10 + ceiling(comb_plot()$num_transcripts * 0.1)
    
    if (heightValue > 49){
      heightValue <- 49
    }
    
    if (widthValue > 49){
      widthValue <- 49
    }

    updateSliderInput(session, inputId = "widthSlider", value = widthValue)
    updateSliderInput(session, inputId = "heightSlider", value = heightValue)
  })
  
  output$densityPlot <- renderPlot({
    exprDensityPlot()
  })

  output$transcriptPlot <- renderPlot({
    comb_plot()$plot
  })
  
  output$plot.ui <- renderUI({
    plotOutput("transcriptPlot", height = comb_plot()$plotHeight)
  })

  output$downloadFig <- downloadHandler(
    filename = function() {
      paste('Transcripts-and-expression-', comb_plot()$gene_name, "-", comb_plot()$gene_id  ,'.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, comb_plot()$plot, width = input$widthSlider, height = input$heightSlider, dpi = 300)
    }
  )
  
  output$downloadDensityFig <- downloadHandler(
    filename = function() {
      paste('Density-plot-highlight-', comb_plot()$gene_name, "-", comb_plot()$gene_id  ,'.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, exprDensityPlot(), width = input$widthDensitySlider, height = input$heightDensitySlider, dpi = 300)

    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
