GetMTcluster <- function(data, ClustNum, ColName){
  MT_cluster <- kmeans(data$PropPos, centers = ClustNum)
  MT_cluster$Cluster <- MT_cluster$cluster
  data[[ColName]] <- "NonMT_PD"
  if(ClustNum == 2){
    data[[ColName]][MT_cluster$Cluster == order(MT_cluster$centers[,1])[1]] <- "MT_PD"
  } else if(ClustNum == 3){
    data[[ColName]][MT_cluster$Cluster == order(MT_cluster$centers[,1])[2]] <- "MT_PDa"
    data[[ColName]][MT_cluster$Cluster == order(MT_cluster$centers[,1])[1]] <- "MT_PDb"
  }
  return(data)
}

PlotClusters <- function(Data = Metadata, GroupCol = "GroupPD", 
                         ClusterColumn,colors, title, showLegend = T){
  ggplot(Data,aes_string(GroupCol,
                         "PercentPositive",
                         color = ClusterColumn)) +
    theme_classic() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x="", y = "Positive neurons (%)", title = title) + 
    geom_quasirandom(show.legend = showLegend) +
    scale_color_manual(values = colors) +
    facet_wrap(~Cohort2, scales = "free_x", nrow = 1)
}

Count2CPM <- function(countData){
  apply(countData, 2, function(smp){
    TotalCount = sum(smp)
    (10^6)*smp/TotalCount
  })
}

PCAgo <- function(data, sampleRegEx = "SL", GOterm){
  data %<>% filter(GeneSymbol %in% GOterm)
  PCA <- prcomp(t(data %>%
                    select(matches(sampleRegEx))), scale = T)
  #Remove genes with different sign
  if(sum(PCA$rotation[,1]) > 0){
    GenesRM <- names(PCA$rotation[,1])[PCA$rotation[,1] < 0]
  } else {
    GenesRM <- names(PCA$rotation[,1])[PCA$rotation[,1] > 0]
  }
  
  data %<>% filter(!Probe %in% GenesRM)
  PCA <- prcomp(t(data %>% select(matches(sampleRegEx))), scale = T)
  #Fix direction 
  if(sum(PCA$rotation[,1]) < 0){
    PCA$x[,1]  = -1*PCA$x[,1]
    PCA$rotation[,1]  = -1*PCA$rotation[,1]
  }
  
  PCA_EffectSize <- t(data %>% select(matches(sampleRegEx))) %*% rescale(PCA$rotation[,1],
                                                                         c(0, 1/nrow(data)))
  
  return(list(GenesRM = GeneSymbolAll %>%
                filter(Probe %in% GenesRM) %>% .$GeneSymbol,
              GeneIn = data$GeneSymbol,
              PCA = PCA,
              EffectSize = PCA_EffectSize))
}

AddGOmgpTometa <- function(Meta, MGP, IDcol = SampleIDcol){
  temp <- MGP[match(Meta[[IDcol]], names(MGP))]
  rescale(temp, c(0,1))
}

DESeq2RUN <- function(data, Meta, model, ContGenes = NULL){
  Meta %<>% droplevels()
  data <- data[,Meta[[SampleIDcol]]]
  if(is.null(ContGenes)){
    ContGenes = 1:nrow(data)
  }
  ModelMatrix <- model.matrix(model, Meta)
  DESeqDS <- DESeqDataSetFromMatrix(countData = data, colData = Meta, design = ModelMatrix)
  DESeqDS <- estimateSizeFactors(DESeqDS, controlGenes = ContGenes)
  DESeqDS <- estimateDispersions(DESeqDS, fitType = "local")
  DESeqOut <- nbinomWaldTest(DESeqDS, modelMatrix = ModelMatrix)
  return(DESeqOut)
}

GetDESeq2Results <- function(DESeqOut, coef, alpha = 0.05, indepFilter = TRUE){
  DEresults <- results(DESeqOut, name = coef, alpha = alpha,
                       format = "DataFrame", independentFiltering = indepFilter) 
  print(summary(DEresults))
  DEresults %<>% data.frame()
  DEresults$EnsemblID <- rownames(DEresults)
  DEresults <- merge(DEresults, geneNames %>% filter(!duplicated(gene_id)) %>%
                       select(gene_id2, gene_name, gene_type), by.x = "EnsemblID", by.y = "gene_id2") %>% arrange(padj)
  return(DEresults)
}

GetOneSidedPval <- function(ResultsObj, adjust = "BH", logFCcol = "log2FoldChange", GeneCol = "GeneSymbol", pvalCol = "pvalue"){
  DESeqResultsDF <- data.frame(ResultsObj)
  DESeqResultsDF$DownPval <- apply(DESeqResultsDF %>% select_(.dots = c(logFCcol, pvalCol)), 1, function(x){
    if(x[1] < 0){
      x[2]/2
    } else {
      1-x[2]/2
    }
  })
  DESeqResultsDF$DownPvalAdj <- p.adjust(DESeqResultsDF$DownPval, adjust)
  DESeqResultsDF$UpPval <- apply(DESeqResultsDF %>% select_(.dots = c(logFCcol, pvalCol)), 1, function(x){
    if(x[1] > 0){
      x[2]/2
    } else {
      1-x[2]/2
    }
  })
  DESeqResultsDF$UpPvalAdj <- p.adjust(DESeqResultsDF$UpPval, adjust)
  rownames(DESeqResultsDF) <- as.character(DESeqResultsDF[[GeneCol]])
  return(DESeqResultsDF)
}


RunEnrich <- function(DEobj, GeneCol = "gene_name", ScoreCol = "pvalue",  minsize = 15, maxsize = 500, method = "both", ErmineJmod = "gsr", Pathways = PathwaysCombined){
  DEresultFiltered <- DEobj %>% filter(!is.na(padj)) %>% filter(!duplicated(gene_name))
  rownames(DEresultFiltered) <- DEresultFiltered[[GeneCol]]
  EnrichList <- list()
  if(method != "ErmineJ"){
    Ranks = deframe(DEresultFiltered %>% select("gene_name", ScoreCol))
    Pathways2 <- Pathways[lapply(Pathways, function(Ptwy){
      length(Ptwy) >= minsize & length(Ptwy) <= maxsize
    }) %>% unlist]
    EnrichList$fgsea <- fgseaMultilevel(pathways=Pathways2, stats=Ranks, nPermSimple = 1000)
    EnrichList$fgsea$GeneSet <- sapply(EnrichList$fgsea$pathway, function(x) strsplit(x, "__")[[1]][1])
    EnrichList$fgsea$pathway <- sapply(EnrichList$fgsea$pathway, function(x) strsplit(x, "__")[[1]][2])
  }
  if(method != "FGSEA"){
    if(ErmineJmod == "gsr"){
      if(ScoreCol == "stat"){
        ScoreColUp = "stat"
        ScoreColDown = "stat"
        bigIsBetterUp = T
        bigIsBetterDown = F
        logTrans = F
      } else if (ScoreCol == "pvalue"){
        DEresultFiltered <- GetOneSidedPval(DEresultFiltered, adjust = "BH", GeneCol = GeneCol)
        ScoreColUp = "UpPval"
        ScoreColDown = "DownPval"
        bigIsBetterUp = F
        bigIsBetterDown = F
        logTrans = T
      }
      EnrichUp =  gsr(scores = DEresultFiltered, scoreColumn = ScoreColUp, customGeneSets = Pathways, bigIsBetter = bigIsBetterUp, logTrans = logTrans, minClassSize = minsize, maxClassSize = maxsize)
      EnrichUp$results %<>% select(-Name, -"Same as", -NumProbes) %>% mutate(Direction = "Up")
      EnrichUp$results$GeneSet <- sapply(EnrichUp$results$ID, function(x) strsplit(x, "__")[[1]][1])
      EnrichUp$results$ID <- sapply(EnrichUp$results$ID, function(x) strsplit(x, "__")[[1]][2])
      EnrichDown = gsr(scores = DEresultFiltered, scoreColumn = ScoreColDown, customGeneSets = Pathways, bigIsBetter = bigIsBetterDown, logTrans = logTrans,  minClassSize = minsize, maxClassSize = maxsize)
      EnrichDown$results %<>% select(-Name, -"Same as", -NumProbes) %>% mutate(Direction = "Down")
      EnrichDown$results$GeneSet <- sapply(EnrichDown$results$ID, function(x) strsplit(x, "__")[[1]][1])
      EnrichDown$results$ID <- sapply(EnrichDown$results$ID, function(x) strsplit(x, "__")[[1]][2])
      EnrichList$ErmineJ <- list(EnrichUp = EnrichUp, EnrichDown = EnrichDown)
    }
  }
  return(EnrichList)
}


GetCountMatrix = function(){
  readRDS("../BigRNAproject/Data/txi.Rds") %>%
    .$counts
}

MetaFilterFun <- function(Metadata){
  Metadata %>% filter(Group2 %in% c("Control", "PD"), Group != "LRRK2_PD",
                      Cohort2 %in% c("Norway", "Barcelona Brain Bank"),
                      Age > 40) %>%
    droplevels()
}

GetCellularChanges <- function(BasicModel = "~MitoType + Sex + Age  +  PMI_binned + Cohort2", Meta = studyFinal$Metadata, Ylim = c(-0.5, 0.5)){
  GroupChanges <- list()
  GroupChangesDF <- list()
  GroupChanges[[1]] <- sapply(names(Meta)[grepl("_Genes", names(Meta))],
                              function(celltype){
                                lm(as.formula(paste0(celltype, BasicModel)),
                                   data = Meta)
                              }, simplify = F)
  
  GroupChangesDF[[1]] <- sapply(names(GroupChanges[[1]]), function(CellType){
    lm.data <- GroupChanges[[1]][[CellType]]
    data <- lm.data  %>% summary %>% .$coef %>% data.frame()
    names(data)[4] <- "pValue"
    DF <- data[2:3,c(1,4)] %>% data.frame() 
    ConfInt <- confint.lm(lm.data)[2:3,]
    DF <- cbind(DF, ConfInt)
    DF %<>% mutate(CellType = CellType,
                   Group = row.names(.))
    DF$Group <- sapply(DF$Group, function(x){
      gsub("MitoType", "", x)
    })
    DF$CellTypeName <- sapply(DF$CellType, function(x){
      gsub("_Genes", "", x)
    })
    names(DF)[3:4] <- c("Low", "High")
    DF %>% mutate(across(where(is.numeric), signif, digits = 2)) %>%
      filter(!grepl("Synaps|OxP", .$CellType))
  }, simplify = F) %>% rbindlist() %>% data.frame() %>% arrange(desc(Estimate))
  
  GroupChangesDF[[1]]$Color = "black"
  
  CellTypeNameLevels <- unique(GroupChangesDF[[1]]$CellTypeName)
  GroupChangesDF[[1]]$CellTypeName <- factor(GroupChangesDF[[1]]$CellTypeName, levels = CellTypeNameLevels) 
  
  names(GroupChangesDF)[1] <- BasicModel
  
  TopCell = GroupChangesDF[[1]] %>% arrange(pValue) %>% head(1) %>% .$CellType
  i = 2
  
  while(nrow(GroupChangesDF[[i-1]] %>% filter(pValue < 0.05, !CellType %in% TopCell[-length(TopCell)], Group == "CI-PD")) > 0){
    Model = paste0(BasicModel, " + ", paste0(TopCell, collapse = " + "))
    GroupChanges[[i]] <- sapply(names(Meta)[grepl("_Genes", names(Meta))],
                                function(celltype){
                                  lm(as.formula(paste0(celltype, Model)),
                                     data = studyFinal$Metadata)
                                }, simplify = F)
    
    GroupChangesDF[[i]] <- sapply(names(GroupChanges[[i]]), function(CellType){
      lm.data <- GroupChanges[[i]][[CellType]]
      data <- lm.data  %>% summary %>% .$coef %>% data.frame()
      names(data)[4] <- "pValue"
      DF <- data[2:3,c(1,4)] %>% data.frame() 
      ConfInt <- confint.lm(lm.data)[2:3,]
      DF <- cbind(DF, ConfInt)
      DF %<>% mutate(CellType = CellType,
                     Group = row.names(.))
      DF$Group <- sapply(DF$Group, function(x){
        gsub("MitoType", "", x)
      })
      DF$CellTypeName <- sapply(DF$CellType, function(x){
        gsub("_Genes", "", x)
      })
      names(DF)[3:4] <- c("Low", "High")
      DF %>% mutate(across(where(is.numeric), signif, digits = 2)) %>%
        filter(!grepl("Synaps|OxP", .$CellType))
    }, simplify = F) %>% rbindlist() %>% data.frame() %>% arrange(pValue)
    GroupChangesDF[[i]]$CellTypeName <- factor(GroupChangesDF[[i]]$CellTypeName, levels = CellTypeNameLevels)
    GroupChangesDF[[i]]$Color = "black"
    GroupChangesDF[[i]]$Color[GroupChangesDF[[i]]$CellType %in% TopCell] <- NA
    TopCell = c(TopCell, GroupChangesDF[[i]] %>%
                  filter(!CellType %in% TopCell, Group == "CI-PD") %>% arrange(pValue) %>% head(1) %>% .$CellType )
    names(GroupChangesDF)[i] <- Model
    i = i+1
  }
  
  PlotList <- sapply(names(GroupChangesDF), function(ModelName){
    StatDF <- GroupChangesDF[[ModelName]]
    ggplot(StatDF, aes(CellTypeName, Estimate)) +
      theme_bw() +
      ylim(Ylim) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(x = "", title = ModelName) +
      geom_hline(yintercept = 0, color = "red", lty = "dashed") +
      geom_point(color = StatDF$Color) +
      geom_errorbar(aes(ymin = Low, ymax = High), color = StatDF$Color) +
      facet_wrap(~Group)
  }, simplify = F)
  FinalPlot <- ggarrange(plotlist = PlotList, nrow = 2)
  return(list(ChangesList = GroupChangesDF, Plot = FinalPlot))
}

PlotCountsMeta <- function(Gene, dds, Grp = "MT_GroupAll2", CompareVar = "DV200", colorBy = NULL,
                           colorBy2 = "PDgroup",
                           colors = c("grey", "orange", MoviePalettes$BugsLife[c(3,2)])){
  
  if(is.null(colorBy)){
    colorBy = Grp
  }
  geneID = geneNames %>% filter(gene_name == Gene, !duplicated(gene_name)) %>% .$gene_id2
  temp <- plotCounts(dds, intgroup = Grp, gene = geneID, normalized = T, returnData = T)
  temp <- merge(dds@colData %>% data.frame %>%
                  select(-Grp), temp, 
                by.y = "row.names", by.x = "RNAseq_id_ParkOme2", sort = F)
  temp %<>% mutate(logCount = log(count + 1)) %>% droplevels()
  Plot1 <- ggplot(temp, aes_string(CompareVar, "logCount", color = colorBy)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(title.position = "top")) +
    geom_point() +
    scale_color_manual(values = colors, name = "Group")
  
  Plot2 <- ggplot(temp, aes_string(Grp, "logCount")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "") +
    guides(colour = guide_legend(title.position = "top")) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, aes_string(color = colorBy2))
  
  Plot3 <- annotate_figure(ggarrange(Plot1, Plot2), fig.lab = Gene, fig.lab.size = 14)
  return(list(Data = temp, Plot1 = Plot1, Plot2 = Plot2, Plot3 = Plot3))
}

PlotGeneByGroup <- function(Gene){
  Data <-  PlotCountsMeta(Gene = Gene, CompareVar = "PercentPositive",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "Cohort2")$Data
  ggplot(Data, aes(MitoType_plot, logCount)) +
    theme_classic() +
    labs(x = "", y = "Expression (log2)", title = paste0(Gene," gene expression")) +
    geom_boxplot(outlier.shape = NA, aes_string(fill = ClusterCol), show.legend = F) +
    #geom_jitter(width = 0.2, height = 0) +
    scale_fill_manual(values = c("grey90", "orange",
                                 MoviePalettes$BugsLife[c(3,2)]), name = "")
  ggsave(paste0(ResultsPath, Gene, "ByGroup.pdf"), device = "pdf",
         width = 4, height = 2.5, useDingbats = F)
}

Enrichheatmap <- function(GeneSet, Ptwy, annoCol = AnnoDF , SampleNaming = "Paper_ID",
                          OrderCols = c("Cluster", "OrdPCA", "Groups"),
                          cutree_rows = 3, cutree_cols = 3,
                          show_rownames = T, show_colnames = T,
                          method = "complete"){
  OrderCols = match.arg(OrderCols)
  temp <- studyFinal$ExpHigh %>% filter(GeneSymbol %in% PathwaysList[[GeneSet]][[Ptwy]]) %>%
    select(matches("SL"))
  rownames(temp) <- studyFinal$ExpHigh %>%
    filter(GeneSymbol %in% PathwaysList[[GeneSet]][[Ptwy]]) %>% .$GeneSymbol
  PCA <- prcomp(t(temp), scale = T)
  if(sum(PCA$rotation[,1]) < 0){
    PCA$x[,1] <- -1*PCA$x[,1]
  }
  if(OrderCols == "Cluster"){
    cluster_cols = T
  } else {
    cluster_cols = F
    if(OrderCols == "OrdPCA"){
      Order = PCA$x[,1] %>% sort %>% names
    } else if(OrderCols == "Groups"){
      Order = studyFinal$Metadata %>% arrange(MitoType, OxPhosES) %>% .$RNAseq_id_ParkOme2
    }
    temp %<>% select(Order)
  }
  
  tempScaled <- apply(temp, 1, scale) %>% t
  colnames(tempScaled) <- studyFinal$Metadata[[SampleNaming]][match(colnames(temp), studyFinal$Metadata$RNAseq_id_ParkOme2)] 
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  
  pheatmap(tempScaled, angle_col = 90, na_col = "white",
           border_color = NA, cluster_cols = cluster_cols, breaks = seq(-4, 4, by = 0.1),
           cutree_rows = cutree_rows, cutree_cols = cutree_cols,
           show_rownames = show_rownames,  show_colnames = show_colnames,
           clustering_method = method, silent = T,
           fontsize_row = 6, clustering_callback = callback, 
           color = colorRampPalette(c("darkblue", "grey10", "gold2"))(length( seq(-4, 4, by = 0.1))+1),
           annotation_col = annoCol, main = Ptwy, 
           annotation_colors = annoColors)
}
