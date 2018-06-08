library(shiny)
library(shinyjs)
library(visNetwork)
library(shinycssloaders)

source("grandforest-web-common/enrichment.R")
source("grandforest-web-common/targets.R")

tooltip_label <- function(text, tooltip) {
  HTML(sprintf('<span data-toggle="tooltip" data-placement="top" title="%s">%s <i class="fa fa-question-circle"></i></span>', tooltip, text))
}

shinyUI(tagList(
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="style.css"),
    tags$link(rel="stylesheet", type="text/css", href="loader.css"),
    tags$script(type="text/javascript", '/* load tooltips */ $(function () { $(\'[data-toggle="tooltip"]\').tooltip() })')
  ),
  useShinyjs(),
  div(id="loading-content", h2("Loading..."), div(class="loader", "Loading")),
  navbarPage("Grand Forest • Unsupervised", inverse=TRUE,
    footer=column(width=12, hr(), p(paste0("Grand Forest • Unsupervised workflow • Version ", APP_VERSION))),
    tabPanel(HTML("Analysis</a></li><li><a href=\"https://grandforest.compbio.sdu.dk/guide/unsupervised\" target=_blank>User guide</a></li><li><a href=\"https://grandforest.compbio.sdu.dk/#cite\" target=_blank>Cite"),
      sidebarLayout(
        sidebarPanel(width=3,
          tags$h3("Upload data", class="sidebar-top-heading"),
          checkboxInput("useExampleData", "Use example data"),
          conditionalPanel("input.useExampleData == false",
            fileInput("file", tooltip_label("Expression data", "See `User guide` for a description of supported file formats.")),
            checkboxInput("hasSurvival", "Include survival information", value=FALSE),
            conditionalPanel("input.hasSurvival == true",
              fluidRow(
                column(width=6, textInput("timeVar", tooltip_label("Time variable name", "Name of column containing survival times. Values must be numeric."))),
                column(width=6, textInput("statusVar", tooltip_label("Status variable name", "Name of column containing status/event values. 1=event, 0=right censored.")))
              )
            ),
            textInput("clusterVar", tooltip_label("Known cluster variable name (optional)", "Name of column containing known clusters."))
          ),
          selectInput("species", "Species", list(
            "Homo sapiens" = "human",
            "Mus musculus" = "mouse"
          )),
          uiOutput("graphSelect"),
          conditionalPanel("input.graph == 'custom'",
            fileInput("graphFile", tooltip_label("Network file", "See `User guide` for a description of supported file formats."))
          ),
          actionButton("uploadButton", "Submit", class="btn-primary"),
          conditionalPanel("output.hasModel == true",
            h3("Split parameters"),
            uiOutput("summary"),
            numericInput("ntrees", "Number of decision trees", DEFAULT_NUM_TREES, min = MIN_NUM_TREES, max = MAX_NUM_TREES),
            sliderInput("nfeatures", "Number of features to split on", min=MIN_NUM_FEATURES, max=MAX_NUM_FEATURES, value=DEFAULT_NUM_FEATURES, step=1),
            sliderInput("nclusters", "Number of clusters to split into", min=MIN_NUM_CLUSTERS, max=MAX_NUM_CLUSTERS, value=DEFAULT_NUM_CLUSTERS, step=1),
            actionButton("splitButton", "Split selected node", class="btn-primary"),
            actionButton("collapseButton", "Collapse selected node", class="btn-primary")
          )
        ),
        mainPanel(width=9, conditionalPanel("output.hasModel == true",
          fluidRow(
            column(width = 6,
              h3("Split tree"),
              p("Click on a node to select it, then click \"Split selected node\" from the sidebar to split it."),
              wellPanel(
                withSpinner(visNetworkOutput("splitTree", height=500)),
                downloadButton("dlSplitTree", "Download tree", class="btn-sm")
              )
            ),
            column(width = 6,
              h3("Heatmap"),
              p("Heatmap of sample subgroups clustered on computed feature subnetwork."),
              wellPanel(
                withSpinner(plotOutput("featureHeatmap", height=500)),
                conditionalPanel("output.hasSplitSelected == true",
                  fluidRow(
                    column(width=4, downloadButton("dlFeatureHeatmap", "Download heatmap", class="btn-sm")),
                    column(width=4, checkboxInput("featureHeatmapGeneSymbols", "Show gene symbol", value=TRUE))
                  )
                )
              )
            )
          ),
          h3("Feature subnetwork"),
          wellPanel(
            withSpinner(visNetworkOutput("featureGraph")),
            conditionalPanel("output.hasSplitSelected == true",
              fluidRow(
                column(width=4, downloadButton("dlFeatureGraph", "Download network", class="btn-sm")),
                column(width=4, checkboxInput("featureGraphGeneSymbols", "Show gene symbols", value=TRUE))
              )
            )
          ),
          h3("Genes"),
          conditionalPanel("output.hasSplitSelected != true",
            div(class="alert alert-info", p("Please select a split node to see selected genes."))
          ),
          conditionalPanel("output.hasSplitSelected == true",
            div(class="body-tabs", tabsetPanel(
              tabPanel("Feature table",
                dataTableOutput("featureTable"),
                downloadButton("dlFeatureTable", "Download table", class="btn-sm")
              ),
              tabPanel("Gene set enrichment",
                fluidRow(
                  column(width=4, uiOutput("enrichmentTypeSelect")),
                  column(width=4, numericInput("enrichmentPvalueCutoff", "p-value cutoff", value=0.05, min=0, max=1, step=0.01)),
                  column(width=4, numericInput("enrichmentQvalueCutoff", "q-value cutoff", value=0.2, min=0, max=1, step=0.01))
                ),
                actionButton("enrichmentButton", "Run enrichment analysis", class="btn-primary"),
                conditionalPanel("output.hasEnrichmentTable == true",
                  hr(),
                  tabsetPanel(
                    tabPanel("Table",
                      dataTableOutput("enrichmentTable"),
                      downloadButton("dlEnrichmentTable", "Download table", class="btn-sm")
                    ),
                    tabPanel("Dot plot",
                      withSpinner(plotOutput("enrichmentPlot"))
                    )
                  )
                )
              ),
              tabPanel("Drug/miRNA targets",
                conditionalPanel("output.currentSpecies == 'human'",
                  selectInput("targetsType", "Database", gene_target_sources()),
                  actionButton("targetsButton", "Get gene targets", class="btn-primary"),
                  conditionalPanel("output.hasTargetsTable == true",
                    hr(),
                    tabsetPanel(type="tabs",
                      tabPanel("Table",
                        withSpinner(dataTableOutput("targetsTable")),
                        downloadButton("dlTargetsTable", "Download table", class="btn-sm")
                      ),
                      tabPanel("Network",
                        withSpinner(visNetworkOutput("targetsNetwork")),
                        checkboxInput("targetsNetworkSymbols", "Show gene symbols", value=TRUE)
                      )
                    )
                  )
                ),
                conditionalPanel("output.currentSpecies != 'human'",
                  div(class="alert alert-info", "Drug/miRNA target search only supported for Homo sapiens.")
                )
              )
            ))
          ),
          conditionalPanel("output.hasSurvivalData == true",
            h3("Survival curves"),
            wellPanel(
              radioButtons("survivalPlotType", "Plot curves for:", inline=TRUE, choices = list(
                "Selected split" = "selected",
                "All leaf nodes" = "leaves"
              )),
              checkboxInput("survivalPlotShowKnown", "Plot known clusters", value=FALSE),
              withSpinner(plotOutput("survivalPlot"))
            )
          )
        )
      ))
    )
  )
))
