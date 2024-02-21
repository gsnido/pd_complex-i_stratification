packageF("ape")
packageF("phytools")

GetJaccardIndx <- function(Vec1, Vec2){
  Intersection <- intersect(Vec1, Vec2)
  Union <- union(Vec1, Vec2)
  length(Intersection)/length(Union)
}

GetEnrichSimilarity <- function(EnrichResult, FilterCol = "CorrectedPvalue", threshold = 0.05){
  EnrichResult <- EnrichResult %>% data.frame()
  EnrichResult$FilterColumn <- EnrichResult[[FilterCol]]
  
  SignifTerms <- EnrichResult %>% filter(FilterColumn < threshold) %>% arrange(FilterColumn) %>%
    select(ID, NumGenes, CorrectedPvalue, GeneMembers, FilterColumn, GeneSet)
  Duplicates <- SignifTerms %>% filter(duplicated(ID)) %>% .$ID %>% unique
  
  SignifTerms$ID[SignifTerms$ID %in% Duplicates] <- sapply(Duplicates, function(Term){
    SignifTerms %>% filter(ID == Term) %>% mutate(ID = paste0(ID, "_", GeneSet)) %>% .$ID
  }, simplify = F) %>% unlist

  SimilarityMatrix <- sapply(SignifTerms$ID, function(Term){
    Genes1 <- sapply(SignifTerms %>% filter(ID == Term) %>% .$GeneMembers, function(x){
      strsplit(x, "\\|")[[1]]
    }) %>% as.character()
    sapply(SignifTerms$ID, function(Term2){
      Genes2 <- sapply(SignifTerms %>% filter(ID == Term2) %>% .$GeneMembers, function(x){
        strsplit(x, "\\|")[[1]]
      }) %>% as.character()
      GetJaccardIndx(Genes1, Genes2)
    })
  }, simplify = F) %>% do.call(cbind, .) %>% data.frame()
  
  #This is to fix the conversion of "-" to "."
  colnames(SimilarityMatrix) <- rownames(SimilarityMatrix)
  
  return(list(Terms = SignifTerms %>% select(-GeneMembers), SimilarityMatrix = SimilarityMatrix))
} 

PlotPhyloTree <- function(ObjName, colors = MoviePalettes$SpiritedAway, edge.width = NULL,
                          h = 1.46, cex = 0.4, ShowLabel = T, Save = F){
  
  SimilarityMatrix <- get(ObjName)$SimilarityMatrix
  mat_cluster_cols <- hclust(as.dist(1-SimilarityMatrix),method = "ward.D2")
  
  #packageF("dendsort")
  #sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  #mat_cluster_cols <- sort_hclust(mat_cluster_cols)
  
  
  Clusters = cutree(mat_cluster_cols, h = h)
  
  
  OrderedLabels <- match(mat_cluster_cols$labels[mat_cluster_cols$order], names(Clusters[!duplicated(Clusters)]))
  ClustersOrder <-  Clusters[!duplicated(Clusters)][OrderedLabels[!is.na(OrderedLabels)]]
  Clusters = factor(Clusters, levels = ClustersOrder)
  
  temp <- as.phylo(mat_cluster_cols)
  temp2 <- temp
  temp2$edge <- cbind(temp$edge, rep(0, nrow(temp$edge)))
  temp2$edge <- cbind(temp2$edge, rep("black", nrow(temp$edge)))
  MedJC <- sapply(levels(Clusters), function(Clust){
    subData <- SimilarityMatrix[names(Clusters[Clusters == Clust]),
                                names(Clusters[Clusters == Clust])]
    median(subData %>% unlist, na.rm = T)
  })
  
  Low <- which(MedJC == min(MedJC))
  colors2 <- c(colors,  MoviePalettes$MadMaxDesert)
  colors2[Low] <- "grey"
  
  for(i in 1:length(levels(Clusters))){
    temp2$edge[which.edge(temp, group = names(Clusters[Clusters == ClustersOrder[i]])),3] <- 3*MedJC[i] + 1
  }
  
  if(is.null(edge.width)){
    edge.width =  ceiling(as.numeric(temp2$edge[,3]))
    edge.width[temp2$edge[,3] == "0"] <- 1
  }
  
  for(i in 1:length(levels(Clusters))){
    temp2$edge[which.edge(temp, group = names(Clusters[Clusters == ClustersOrder[i]])),4] <- colors2[i]
  }
  
  if(Save){
   pdf(paste0(ResultsPath, ObjName, ".pdf"), width = 8, height = 8) 
  }

  plot.phylo(temp, type = "fan", cex = cex, show.tip.label = ShowLabel, label.offset = 0.1, no.margin	 = T, 
             tip.color = colors2[Clusters], plot = T, 
             edge.width = edge.width, edge.color = temp2$edge[,4])
  
  add.simmap.legend(colors=unique(colors2), x=0.9*par()$usr[1],leg = signif(MedJC, digits = 2), main = "Median Jaccard Similarity",
                    y=0.9*par()$usr[4], prompt=F,fsize=0.9)

  text(x = 0.9*par()$usr[1], y=0.95*par()$usr[4], "Median Jaccard Similarity", adj = 0)
  
  if(Save){
    closeDev()
  }
}

EnrichSimilarityNonMT2Down <- GetEnrichSimilarity(EnrichResult = Enrichment_NonMT2$ErmineJ$EnrichDown$results %>%
                                                    filter(GeneSet != "MitoCarta"))
diag(EnrichSimilarityNonMT2Down$SimilarityMatrix) <- NA

PlotPhyloTree(ObjName = "EnrichSimilarityNonMT2Down", edge.width = 2, cex = 0.5, Save = T, h = 1.39)

EnrichSimilarityNonMT2Up <- GetEnrichSimilarity(EnrichResult = Enrichment_NonMT2$ErmineJ$EnrichUp$results %>%
                                                  filter(GeneSet != "MitoCarta"))
diag(EnrichSimilarityNonMT2Up$SimilarityMatrix) <- NA


PlotPhyloTree(ObjName = "EnrichSimilarityNonMT2Up", edge.width = 2,
              h = 1.05, colors = MoviePalettes$BugsLife[c(2:6,9)], cex = 0.8, Save = T)


#Repeat for BBB only
EnrichSimilarityNonMT3Down_BBB <- GetEnrichSimilarity(EnrichResult = Enrichment_BBB_nonMT$ErmineJ$EnrichDown$results %>%
                                                    filter(GeneSet != "MitoCarta"))
diag(EnrichSimilarityNonMT3Down_BBB$SimilarityMatrix) <- NA

PlotPhyloTree(ObjName = "EnrichSimilarityNonMT3Down_BBB", edge.width = 2, cex = 0.4, Save = T, h = 1.75)


EnrichSimilarityNonMT3Up_BBB <- GetEnrichSimilarity(EnrichResult = Enrichment_BBB_nonMT$ErmineJ$EnrichUp$results %>%
                                                      filter(GeneSet != "MitoCarta"))
diag(EnrichSimilarityNonMT3Up_BBB$SimilarityMatrix) <- NA


PlotPhyloTree(ObjName = "EnrichSimilarityNonMT3Up_BBB", edge.width = 2,
              h = 1.6, colors = MoviePalettes$BugsLife, cex = 0.5, Save = T)



#Corrletaion plots

CreateCompareDF <- function(MTgroup, Direction){
  Combined <- get(paste0("Enrichment_", MTgroup, "3")) %>% .$ErmineJ %>% .[[paste0("Enrich", Direction)]] %>% .$results
  NTB <- get(paste0("Enrichment_BBB_", MTgroup)) %>% .$ErmineJ %>% .[[paste0("Enrich", Direction)]] %>% .$results
  SignifCombined <- Combined %>% filter(CorrectedPvalue < 0.05) %>% .$ID
  SignifNTB <- NTB %>% filter(CorrectedPvalue < 0.05) %>% .$ID
  Combined_Sub <- Combined %>% filter(ID %in% c(SignifCombined, SignifNTB)) %>% data.frame() %>%
    select(ID, Pval, CorrectedPvalue) %>% mutate(logP = -1*log10(Pval), logAdjP = -1*log10(CorrectedPvalue))
  NTB_Sub <- NTB %>% filter(ID %in% c(SignifCombined, SignifNTB)) %>% data.frame() %>%
    select(ID, Pval, CorrectedPvalue) %>% mutate(logP = -1*log10(Pval), logAdjP = -1*log10(CorrectedPvalue))
  merge(Combined_Sub, NTB_Sub, by = "ID", suffixes = c("_combined", "_NTB")) %>%
    arrange(Pval_combined) %>% filter(!duplicated(ID))
  }


ComareDF_MTPD_up <- CreateCompareDF(MTgroup = "MTPD", Direction = "Up") 
ComareDF_MTPD_down <- CreateCompareDF(MTgroup = "MTPD", Direction = "Down")

ComareDF_NonMT_up <- CreateCompareDF(MTgroup = "NonMT", Direction = "Up")
ComareDF_NonMT_down <- CreateCompareDF(MTgroup = "NonMT", Direction = "Down")


ggplot(ComareDF_NonMT_up, aes(logAdjP_combined, logAdjP_NTB)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = -log10(0.05), color = "red")


ggplot(ComareDF_NonMT_down, aes(logP_combined, logP_NTB)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = -log10(0.05), color = "red")


ggplot(ComareDF_MTPD_up, aes(logAdjP_combined, logAdjP_NTB)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = -log10(0.05), color = "red")


ggplot(ComareDF_MTPD_down, aes(logP_combined, logP_NTB)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = -log10(0.05), color = "red")

#Just exploring
#Repeat for mitoPD
# EnrichSimilarityMTPD3Up <- GetEnrichSimilarity(EnrichResult = Enrichment_MTPD3$ErmineJ$EnrichUp$results %>%
#                                                   filter(GeneSet != "MitoCarta"))
# diag(EnrichSimilarityMTPD3Up$SimilarityMatrix) <- NA
# 
# PlotPhyloTree(ObjName = "EnrichSimilarityMTPD3Up", edge.width = 2,
#               h = 1.55, colors = MoviePalettes$BugsLife, cex = 0.4, Save = T)
# 
# #Due to the large number of significant terms for the downreagulated genes,
# #running separately for each gene set type
# 
# EnrichSimilarityMTPD3DownGObp <- GetEnrichSimilarity(EnrichResult = Enrichment_MTPD3$ErmineJ$EnrichDown$results %>%
#                                                     filter(GeneSet == "GObp"))
# diag(EnrichSimilarityMTPD3DownGObp$SimilarityMatrix) <- NA
# 
# PlotPhyloTree(ObjName = "EnrichSimilarityMTPD3DownGObp", edge.width = 2,
#               h = 2.28, cex = 0.2, Save = T)
# 
# 
# EnrichSimilarityMTPD3DownKEGG <- GetEnrichSimilarity(EnrichResult = Enrichment_MTPD3$ErmineJ$EnrichDown$results %>%
#                                                        filter(GeneSet == "KEGG"))
# diag(EnrichSimilarityMTPD3DownKEGG$SimilarityMatrix) <- NA
# 
# PlotPhyloTree(ObjName = "EnrichSimilarityMTPD3DownKEGG", edge.width = 2,
#               h = 1.05, cex = 0.6, Save = F)
# 
# 
# packageF("wordcloud")
# packageF("wordcloud2")
# 
# PlotCloud <- function(EnrichObject, ...){
#   temp <-  sapply(EnrichObject %>% .$ID, function(x){
#     strsplit(x, "_")[[1]]
#   }, simplify = F) %>% unlist %>% as.character() %>% toupper()
#   
#   temp[grepl("Mitochon", temp, ignore.case = T)] <- "MITOCHONDRIA"
#   temp[grepl("synaps|synapt", temp, ignore.case = T)] <- "SYNAPSE"
#   temp[grepl("lysosom", temp, ignore.case = T)] <- "LYSOSOME"
#   temp[grepl("axon", temp, ignore.case = T)] <- "AXON"
#   temp[grepl("metabol", temp, ignore.case = T)] <- "METABOLISM"
#   
#   
#   temp %<>% table
#   
#   tempDF <- data.frame(word = names(temp), freq = as.numeric(temp))
#   tempDF %<>% filter(!word %in% c("OF", "TO", "FROM", "VIA", "IN", "AND", "THE", "REGULATION", "PROCESS", "POSITIVE", "NEGATIVE")) %>% arrange(desc(freq))
#   rownames(tempDF) <- tempDF$word
#   
#   wordcloud(words = tempDF$word, freq = tempDF$freq, min.freq = 1,
#             max.words=300, random.order=FALSE, rot.per=0.35,
#             colors=brewer.pal(8, "Dark2"))
#   
# }
#   
# PlotCloud(EnrichObject = Enrichment_Controls$ErmineJ$EnrichDown$results %>% filter(CorrectedPvalue < 0.05))
# 
# 
# EnrichSimilarityMTPD3_GOcc <- GetEnrichSimilarity(EnrichResult = Enrichment_MTPD3$ErmineJ$EnrichDown$results %>%
#                                                filter(GeneSet == "GOcc"))
# 
# diag(EnrichSimilarityMTPD3_GOcc$SimilarityMatrix) <- NA
# mat_cluster_cols <- hclust(dist(EnrichSimilarityMTPD3_GOcc$SimilarityMatrix),method = "ward.D2")
# clus10 = cutree(mat_cluster_cols, 10)
# plot(as.phylo(mat_cluster_cols), type = "fan", cex = 0.4, label.offset = 0.2, tip.color = colors[clus10], main = "GOcc")
# 
# EnrichSimilarityMTPD3_GObp <- GetEnrichSimilarity(EnrichResult = Enrichment_MTPD3$ErmineJ$EnrichDown$results %>%
#                                                     filter(GeneSet == "GObp"))
# 
# diag(EnrichSimilarityMTPD3_GObp$SimilarityMatrix) <- NA
# mat_cluster_cols <- hclust(dist(EnrichSimilarityMTPD3_GObp$SimilarityMatrix),method = "ward.D2")
# clus15 = cutree(mat_cluster_cols, 15)
# EnrichSimilarityMTPD3_GObp$Terms$Cluster <- clus15[match(EnrichSimilarityMTPD3_GObp$Terms$ID, names(clus15))]
# plot(as.phylo(mat_cluster_cols), type = "fan", cex = 0.2, label.offset = 0.2, tip.color = colors[clus10], main = "GObp")
