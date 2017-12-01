library(data.table)
library(magrittr)
library(grandforest)
library(igraph)
library(visNetwork)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(ggplot2)
library(survival)
library(survminer)
library(gridExtra)

source("grandforest-web-common/get_network.R")
source("grandforest-web-common/enrichment.R")
source("grandforest-web-common/feature_graph.R")

topn <- function(x, n) {
  x[order(x, decreasing=TRUE)][1:n]
}

split_importance <- function(D, edges, ntrees) {
  bg <- apply(D, 2, function(x) sample(x, length(x), TRUE))
  all <- rbind(D, bg)
  all[[GROUP_FEATURE_NAME]] <- as.factor(c(rep("fg", nrow(D)), rep("bg", nrow(bg))))
  
  fit <- grandforest(data=all, dependent.variable.name=GROUP_FEATURE_NAME, graph_data=edges, num.trees=ntrees, importance="impurity")
  importance(fit)
}

split_node <- function(tree, i, D, edges, ntrees, nfeatures, nclusters) {
  rows <- tree[[i]]$rows
  if(is.null(tree[[i]]$importance) || tree[[i]]$ntrees != ntrees) {
    tree[[i]]$importance <- split_importance(D[rows,], edges, ntrees)
  }
  
  features <- names(topn(tree[[i]]$importance, nfeatures))
  cl <- kmeans(D[rows,features,with=FALSE], centers=nclusters)
  cluster <- cl$cluster
  
  tlength <- length(tree)
  
  tree[[i]]$ntrees <- ntrees
  tree[[i]]$nfeatures <- nfeatures
  tree[[i]]$cluster <- cluster
  tree[[i]]$children <- seq(tlength+1, tlength+nclusters)
  
  for(j in 1:nclusters) {
    tree[[tlength+j]] <- list(rows = rows[cluster == j])
  }
  
  return(tree)
}

split_tree_reachable_nodes <- function(tree) {
  visit <- c(1)
  reachable <- c()
  i <- 1
  while(i <= length(visit)) {
    node_id <- visit[i]
    if(!is.null(tree[[node_id]]$children)) {
      visit <- c(visit, tree[[node_id]]$children)
    }
    reachable <- c(reachable, node_id)
    i <- i+1
  }
  return(reachable)
}

split_tree_network <- function(tree) {
  if(length(tree) == 1) {
    node_ids <- 1
    edges <- data.frame(from=c(), to=c())
  } else {
    node_ids <- split_tree_reachable_nodes(tree)
    edges <- do.call(rbind, lapply(node_ids, function(i) {
      if(!is.null(tree[[i]]$children)) {
        data.frame(from = i, to = tree[[i]]$children)
      }
    }))
  }
  return(list(nodes=node_ids, edges=edges))
}

shinyServer(function(input, output, session) {
  currentData <- reactiveVal()
  currentSurvivalData <- reactiveVal()
  currentKnownClusters <- reactiveVal()
  currentEdges <- reactiveVal()
  currentTree <- reactiveVal(list())
  lastSplitID <- reactiveVal()
  currentEnrichmentTable <- reactiveVal()
  
  output$hasModel <- reactive({
    req(currentTree())
    return(length(currentTree()) > 0)
  })
  outputOptions(output, "hasModel", suspendWhenHidden=FALSE)
  
  output$hasSplitSelected <- reactive({
    req(input$splitTree_selected)
    node_id <- as.numeric(input$splitTree_selected)
    tree <- currentTree()
    req(tree[[node_id]]$importance)
    return(TRUE)
  })
  outputOptions(output, "hasSplitSelected", suspendWhenHidden=FALSE)
  
  output$hasEnrichmentTable <- reactive({
    req(currentEnrichmentTable())
    return(TRUE)
  })
  outputOptions(output, "hasEnrichmentTable", suspendWhenHidden=FALSE)
  
  output$hasSurvivalData <- reactive({
    req(currentSurvivalData())
    return(TRUE)
  })
  outputOptions(output, "hasSurvivalData", suspendWhenHidden=FALSE)
  
  observeEvent(input$uploadButton, {
    if(!isTruthy(input$file)) {
      showNotification("Please select an expression data file and wait for it to upload before submitting.", type="error")
      return()
    }
    
    withProgress(value=0, message = "Parsing data", {
      setProgress(value=0, detail="Parsing expression data")
      D <- fread(input$file$datapath, header=TRUE, sep=",")
      
      if(GROUP_FEATURE_NAME %in% colnames(D)) {
        showNotification(paste0("Column name \"", GROUP_FEATURE_NAME, "\" is not allowed in expression table."), type="error")
        return()
      }
      
      clusters <- NULL
      if(input$clusterVar != "") {
        if(!(input$clusterVar %in% colnames(D))) {
          showNotification("Known cluster variable name does not match any column name in data file.", type="error")
          return()
        }
        
        clusters <- as.factor(D[[input$clusterVar]])
        cluster_col <- which(colnames(D) == input$clusterVar)
        D <- D[,-cluster_col,with=FALSE]
      }
      
      survival <- NULL
      if(input$hasSurvival) {
        if(!(input$timeVar %in% colnames(D))) {
          showNotification("Survival time variable name does not match any column name in data file.", type="error")
          return()
        }
        if(!(input$statusVar %in% colnames(D))) {
          showNotification("Survival status variable name does not match any column name in data file.", type="error")
          return()
        }
        
        survival <- list()
        survival$time <- as.numeric(D[[input$timeVar]])
        survival$status <- as.numeric(D[[input$statusVar]])
        time_col <- which(colnames(D) == input$timeVar)
        status_col <- which(colnames(D) == input$statusVar)
        D <- D[,-c(time_col,status_col),with=FALSE]
      }
      
      # scale and mean center data
      setProgress(value=0.7, detail="Normalizing data")
      D <- tryCatch({
        data.table(scale(D, center=TRUE, scale=TRUE))
      }, error = function(e) {
        showNotification("Normalization failed. Not all columns are numeric.", type="error")
        req(FALSE)
      })
      
      setProgress(value=0.8, detail="Preparing network")
      graph.path <- get_network_file(input$graph)
      edges <- fread(graph.path, header=FALSE, sep="\t", colClasses=rep("character", 2))
      colnames(edges) <- c("from","to")
      
      setProgress(value=0.9, detail="Finishing up")
      currentData(D)
      currentKnownClusters(clusters)
      currentSurvivalData(survival)
      currentEdges(edges)
      currentTree(list(list(rows=1:nrow(D))))
      lastSplitID(1)
    })
  })
  
  observeEvent(input$splitButton, {
    if(!isTruthy(input$splitTree_selected)) {
      showNotification("Please select a node to split.", type="error")
      return()
    }
    if(input$ntrees > MAX_NUM_TREES || input$ntrees < MIN_NUM_TREES) {
      showNotification(paste0("Number of trees must be >= ", MIN_NUM_FEATURES, " and <= ", MAX_NUM_FEATURES), type="error")
      return()
    }
    if(input$nfeatures > MAX_NUM_FEATURES || input$nfeatures < MIN_NUM_FEATURES) {
      showNotification(paste0("Number of features to split on must be >= ", MIN_NUM_FEATURES, " and <= ", MAX_NUM_FEATURES, "."), type="error")
      return()
    }
    
    withProgress(value=0, message="Splitting node", {
      setProgress(value=0.0, detail = "Preparing data")
      tree <- currentTree()
      node_id <- as.numeric(input$splitTree_selected)
      
      setProgress(value=0.1, detail = "Training model")
      tree <- split_node(tree, node_id, currentData(), currentEdges(), input$ntrees, input$nfeatures, input$nclusters)
      currentTree(tree)
      lastSplitID(node_id)
    })
  })
  
  observeEvent(input$collapseButton, {
    if(!isTruthy(input$splitTree_selected)) {
      showNotification("Please select a node to collapse.", type="error")
      return()
    }
    
    tree <- currentTree()
    node_id <- as.numeric(input$splitTree_selected)
    if(is.null(tree[[node_id]]$children)) {
      showNotification("You can only collapse nodes that have been split.", type="error")
      return()
    }
    
    tree[[node_id]] <- list(rows = tree[[node_id]]$rows)
    
    currentTree(tree)
    lastSplitID(node_id)
  })
  
  currentFeatures <- reactive({
    req(input$splitTree_selected)
    node_id <- as.numeric(input$splitTree_selected)
    tree <- currentTree()
    req(tree[[node_id]]$importance)
    topn(tree[[node_id]]$importance, tree[[node_id]]$nfeatures)
  })
  
  featureTable <- reactive({
    importance <- currentFeatures()
    names <- mapIds(org.Hs.eg.db, names(importance), "SYMBOL", "ENTREZID")
    data.table(gene=names(importance), name=names, importance=importance)
  })
  
  output$splitTree <- renderVisNetwork({
    tree <- currentTree()
    if(length(tree) == 0) return(NULL)
    
    network <- split_tree_network(tree)
    node_ids <- network$nodes
    edges <- network$edges
    
    labels <- paste0(node_ids, ", n=", sapply(node_ids, function(i) length(tree[[i]]$rows)))
    
    if(isTruthy(currentKnownClusters())) {
      known_clusters <- currentKnownClusters()
      labels <- sapply(1:length(node_ids), function(i) {
        node_id <- node_ids[i]
        rows <- tree[[node_id]]$rows
        clusters <- known_clusters[rows]
        most_common <- names(which.max(table(clusters)))
        most_common_pct <- sum(clusters == most_common) / length(clusters) * 100
        paste0(labels[i], sprintf("\n%s (%.1f%%)", most_common, most_common_pct))
      })
    }
    nodes <- data.frame(id=node_ids, label=labels)
    
    visNetwork(nodes, edges) %>%
      visNodes(shape="circle") %>%
      visEdges(arrows="to") %>%
      visInteraction(dragNodes=FALSE) %>%
      visHierarchicalLayout(sortMethod="directed") %>%
      visOptions(nodesIdSelection=list(enabled=TRUE, useLabels=FALSE, selected=isolate(lastSplitID())))
  })
  
  output$featureGraph <- renderVisNetwork({
    req(input$splitTree_selected)
    node_id <- as.numeric(input$splitTree_selected)
    features <- names(currentFeatures())
    
    edges <- currentEdges()
    tree <- currentTree()
    D <- currentData()[tree[[node_id]]$rows,]
    
    labels <- features
    if(input$featureGraphGeneSymbols) {
      labels <- mapIds(org.Hs.eg.db, features, "SYMBOL", "ENTREZID")
    }
    
    group_names <- tree[[node_id]]$children
    groups <- group_names[tree[[node_id]]$cluster]
    feature_graph(D, edges, features, labels, groups)
  })
  
  output$featureTable <- renderDataTable({
    featureTable()
  }, options=list(
    scrollX = TRUE,
    searching = FALSE,
    pageLength = 10
  ))

  featureHeatmapPlot <- reactive({
    req(input$splitTree_selected)
    
    node_id <- as.numeric(input$splitTree_selected)
    tree <- currentTree()
    
    rows <- tree[[node_id]]$rows
    features <- names(currentFeatures())
    gene_names <- mapIds(org.Hs.eg.db, features, "SYMBOL", "ENTREZID")
    group_names <- tree[[node_id]]$children
    groups <- group_names[tree[[node_id]]$cluster]
    
    D <- currentData()
    D <- D[rows,features,with=FALSE]
    colnames(D) <- paste0(features, " (", gene_names, ")")
    
    col.ramp <- colorRamp2(c(-2, 0, 2), c("magenta", "black", "green"))
    
    hm <- Heatmap(D, name="expression", split=groups, col=col.ramp)
    
    # add known cluster annotation if available
    if(isTruthy(currentKnownClusters())) {
      clusters <- currentKnownClusters()
      cluster_levels <- levels(clusters)
      colors <- rainbow(length(cluster_levels))
      anno <- data.frame(cluster=clusters)
      cluster_col <- setNames(colors, cluster_levels)
      hm <- hm + rowAnnotation(anno, col=list(cluster=cluster_col))
    }
    return(hm)
  })
  
  output$featureHeatmap <- renderPlot({
    featureHeatmapPlot()
  })
  
  output$survivalPlot <- renderPlot({
    survival <- req(currentSurvivalData())
    tree <- currentTree()
    
    if(input$survivalPlotType == "selected") {
      req(input$splitTree_selected)
      node_id <- as.numeric(input$splitTree_selected)
      req(!is.null(tree[[node_id]]$importance))
      
      rows <- tree[[node_id]]$rows
      group_names <- tree[[node_id]]$children
      cluster <- as.factor(group_names[tree[[node_id]]$cluster])
    } else if(input$survivalPlotType == "leaves") {
      node_ids <- split_tree_reachable_nodes(tree)
      is_leaf <- sapply(node_ids, function(i) is.null(tree[[i]]$importance))
      node_ids <- node_ids[is_leaf]
      
      rows <- do.call(c, lapply(node_ids, function(i) tree[[i]]$rows))
      cluster <- as.factor(do.call(c, lapply(node_ids, function(i) rep(i, length(tree[[i]]$rows)))))
    } else {
      return()
    }
    
    D <- data.frame(time = survival$time[rows], status = survival$status[rows], cluster = cluster)
    fit <- survfit(Surv(time, status)~cluster, data=D)
    pval <- paste0("p = ", surv_pvalue(fit, data=D)$pval)
    p <- ggsurvplot(fit, data=D, pval=pval)
    
    if(input$survivalPlotShowKnown && isTruthy(currentKnownClusters())) {
      D$cluster <- as.factor(currentKnownClusters()[rows])
      fit <- survfit(Surv(time, status)~cluster, data=D)
      pval <- paste0("p = ", surv_pvalue(fit, data=D)$pval)
      p2 <- ggsurvplot(fit, data=D, pval=pval)
      p <- grid.arrange(p$plot, p2$plot, ncol=2)
    }
    
    return(p)
  })
  
  output$survivalPlotKnown <- renderPlot({
    D <- survivalPlotData()
    D$cluster <- as.factor(req(currentKnownClusters()))
    fit <- survfit(Surv(time, status)~cluster, data=D)
    pval <- paste0("p = ", surv_pvalue(fit, data=D)$pval)
    ggsurvplot(fit, data=D, pval=pval)
  })
  
  observeEvent(input$enrichmentButton, {
    withProgress(message="Performing gene set enrichment", {
      setProgress(value=0.1, detail="Preparing data")
      genes <- names(currentFeatures())
      D <- currentData()
      universe <- colnames(D)
      
      setProgress(value=0.1, detail="Computing enrichment")
      out <- gene_set_enrichment(genes, universe, input$enrichmentType, input$enrichmentPvalueCutoff, input$enrichmentQvalueCutoff)
      
      setProgress(value=0.1, detail="Finishing up")
      currentEnrichmentTable(out)
    })
  })
  
  output$enrichmentTable <- renderDataTable({
    D <- req(as.data.frame(currentEnrichmentTable()))
    D$ID <- gene_set_enrichment_get_links(D$ID, isolate(input$enrichmentType))
    return(D)
  }, options = list(
    scrollX = TRUE,
    pageLength = 10
  ), escape = FALSE)
  
  output$enrichmentPlot <- renderPlot({
    D <- req(currentEnrichmentTable())
    DOSE::dotplot(D, showCategory=20)
  })
  
  output$dlSplitTree <- downloadHandler(
    filename = "split_tree.csv",
    content = function(file) {
      tree <- currentTree()
      node_ids <- split_tree_reachable_nodes(tree)
      D <- do.call(rbind, lapply(node_ids, function(i) {
        data.frame(id=i, children=paste(tree[[i]]$children, collapse=";"), rows=paste(tree[[i]]$rows, collapse=";"))
      }))
      write.csv(D, file, row.names=FALSE)
    }
  )
  
  output$dlFeatureGraph <- downloadHandler(
    filename = function() {
      paste0("network_top", length(currentFeatures()), ".sif")
    },
    content = function(file) {
      features <- names(currentFeatures())

      edges <- currentEdges()
      edges <- subset(edges, from %in% features & to %in% features)
      D <- data.frame(from=edges$from, type=".", to=edges$to)
      if(input$featureGraphGeneSymbols) {
        D$from <- mapIds(org.Hs.eg.db, as.character(D$from), "SYMBOL", "ENTREZID")
        D$to <- mapIds(org.Hs.eg.db, as.character(D$to), "SYMBOL", "ENTREZID")
      }

      write.table(D, file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
    }
  )
  
  output$dlFeatureHeatmap <- downloadHandler(
    filename = function() {
      paste0("heatmap_", input$splitTree_selected, ".pdf")
    },
    content = function(file) {
      pdf(file=file, width=10, height=10)
      print(featureHeatmapPlot())
      dev.off()
    }
  )

  output$dlFeatureTable <- downloadHandler(
    filename = function() {
      paste0("genes_cluster", input$splitTree_selected, ".csv")
    },
    content = function(file) {
      write.csv(featureTable(), file, row.names=FALSE)
    }
  )
  
  output$dlEnrichmentTable <- downloadHandler(
    filename = "enrichment.csv",
    content = function(file) {
      write.csv(as.data.frame(currentEnrichmentTable()), file, row.names=FALSE)
    }
  )
})
