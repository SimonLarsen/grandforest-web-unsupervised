library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinysky)
library(visNetwork)

source("grandforest-web-common/enrichment.R")

shinyUI(navbarPage("Grand Forest • Unsupervised", theme=shinytheme("cosmo"),
  footer=column(width=12, hr(), p("Grand Forest • Unsupervised workflow • Version 0.1")),
  tabPanel("Analysis",
    tags$head(
      tags$link(rel="stylesheet", type="text/css", href="style.css")
    ),
    sidebarLayout(
      sidebarPanel(width = 3,
        h3("Upload data"),
        fileInput("file", "Expression table (.csv file)", accept = "text/csv"),
        checkboxInput("hasSurvival", "Include survival information", value=FALSE),
        conditionalPanel("input.hasSurvival == true",
          textInput("timeVar", "Time variable name"),
          textInput("statusVar", "Status variable name")
        ),
        textInput("clusterVar", "Known cluster variable name (optional)"),
        selectInput("graph", "Genetic interaction network",
          list(
            "IID, Human, Experimental only" = "iidexp",
            "IID, Human, Full" = "iidall",
            "RegNetwork" = "regnetwork",
            "BioGRID" = "biogrid",
            "HTRIdb" = "htri"
          )
        ),
        actionButton("uploadButton", "Submit", styleclass = "primary"),
        conditionalPanel("output.hasModel == true",
          h3("Split parameters"),
          numericInput("ntrees", "Number of decision trees", DEFAULT_NUM_TREES, min = MIN_NUM_TREES, max = MAX_NUM_TREES),
          sliderInput("nfeatures", "Number of features to split on", min=MIN_NUM_FEATURES, max=MAX_NUM_FEATURES, value=DEFAULT_NUM_FEATURES, step=1),
          sliderInput("nclusters", "Number of clusters to split into", min=MIN_NUM_CLUSTERS, max=MAX_NUM_CLUSTERS, value=DEFAULT_NUM_CLUSTERS, step=1),
          actionButton("splitButton", "Split selected node", styleclass = "primary"),
          actionButton("collapseButton", "Collapse selected node", styleclass="primary")
        )
      ),
      mainPanel(conditionalPanel("output.hasModel == true",
        fluidRow(
          column(width = 6,
            h3("Split tree"),
            p("Click on a node to select it, then click \"Split selected node\" from the sidebar to split it."),
            wellPanel(
              visNetworkOutput("splitTree", height=500),
              downloadButton("dlSplitTree", "Download tree", class="btn-sm")
            )
          ),
          column(width = 6,
            h3("Heatmap"),
            p("Heatmap of sample subgroups clustered on computed feature subgraph."),
            wellPanel(
              plotOutput("featureHeatmap", height=500),
              conditionalPanel("output.hasSplitSelected == true",
                downloadButton("dlFeatureHeatmap", "Download heatmap", class="btn-sm")
              )
            )
          )
        ),
        fluidRow(
          column(width = 6,
            h3("Feature subgraph"),
            wellPanel(
              visNetworkOutput("featureGraph"),
              fluidRow(
                conditionalPanel("output.hasSplitSelected == true",
                  column(width = 4, downloadButton("dlFeatureGraph", "Download network", class="btn-sm")),
                  column(width = 4, checkboxInput("featureGraphGeneSymbols", "Show gene symbols", value=TRUE))
                )
              )
            )
          ),
          column(width = 6,
            h3("Top features"),
            wellPanel(
              dataTableOutput("featureTable"),
              conditionalPanel("output.hasSplitSelected == true",
                downloadButton("dlFeatureTable", "Download table", class="btn-sm")
              )
            )
          )
        ),
        conditionalPanel("output.hasSurvivalData == true",
          h3("Survival curves"),
          conditionalPanel("output.hasSplitSelected != true",
            tags$div(class="alert alert-info",
              p("Please select a split node to show survival curves.")
            )
          ),
          conditionalPanel("output.hasSplitSelected == true",
            plotOutput("survivalPlot")
          )
        ),
        h3("Gene set enrichment"),
        conditionalPanel("output.hasSplitSelected != true",
          tags$div(class="alert alert-info",
            p("Please select a split node to perform gene set enrichment")
          )
        ),
        conditionalPanel("output.hasSplitSelected == true",
          wellPanel(
            fluidRow(
              column(width=4, selectInput("enrichmentType", "Enrichment type", gene_set_enrichment_types())),
              column(width=4, numericInput("enrichmentPvalueCutoff", "p-value cutoff", value=0.05, min=0, max=1, step=0.01)),
              column(width=4, numericInput("enrichmentQvalueCutoff", "q-value cutoff", value=0.2, min=0, max=1, step=0.01))
            ),
            actionButton("enrichmentButton", "Run enrichment analysis", styleclass="primary"),
            conditionalPanel("output.hasEnrichmentTable == true",
              hr(),
              tabsetPanel(
                tabPanel("Table",
                  dataTableOutput("enrichmentTable"),
                  downloadButton("dlEnrichmentTable", "Download table", class="btn-sm")
                ),
                tabPanel("Plot",
                  plotOutput("enrichmentPlot")
                )
              )
            )
          )
        )
      )
    ))
  ),
  tabPanel("User guide"),
  tabPanel("Cite")
))
