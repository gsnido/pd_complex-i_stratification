ResultsPath = paste0("Results_", Study, "_", strsplit(AssemblyFilename, "\\.")[[1]][2], "_RNAseq")

if(!ResultsPath %in% list.dirs(full.names = F, recursive = T)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

plotMA = DESeq2::plotMA

studyFinal$Metadata$MitoType <- studyFinal$Metadata[[ClusterCol]]

SubjNumber <- studyFinal$Metadata %>% group_by(MitoType) %>%
  summarise(n = n()) %>% data.frame() %>%
  mutate(MTgroupPlot = paste0(MitoType, "\n(n=", n, ")"))

studyFinal$Metadata$MitoType_plot <- SubjNumber$MTgroupPlot[match(studyFinal$Metadata$MitoType,
                                                                  SubjNumber$MitoType)] 

studyFinal$Metadata$MitoType_plot <- factor(studyFinal$Metadata$MitoType_plot, levels = SubjNumber$MTgroupPlot)


CovarCor = cor(studyFinal$Metadata %>% select(PC1, DV200, RIN, PMI,  matches("Genes")))


##### Not part of the manuscript anymore
# ggplot(studyFinal$Metadata, aes(PercentPositive, OxPhosES)) +
#   theme_minimal() +
#   geom_point(aes(color = MitoType)) +
#   scale_color_manual(values = c("grey70", "orange",
#                                 MoviePalettes$BugsLife[c(2)]), name = "") +
#   geom_smooth(method = "lm")
# 
# ggplot(studyFinal$Metadata, aes(MitoType_plot, OxPhosES,fill = MitoType, color = MitoType)) +
#   theme_minimal() +
#   labs(x = "", y = "Expression (log2)") +
#   geom_boxplot(outlier.shape = NA, width = 0.12, alpha = 0.5, color = "black") +
#   geom_jitter(width = 0.05, height = 0, size = 0.5, color = "black") +
#   stat_halfeye(adjust = 0.5, justification = -0.15, .width = 0, width = 0.9, point_colour = NA) +
#   scale_fill_manual(values = c("grey70", "orange",
#                                MoviePalettes$BugsLife[c(2)]), name = "") +
#   scale_color_manual(values = c("grey70", "orange",
#                                MoviePalettes$BugsLife[c(2)]), name = "")
# 
# ggsave(paste0(ResultsPath, "OxPhosByGroup.pdf"), device = "pdf",
#        width = 8, height = 4, useDingbats = F)


SampleCor2 <- cor(studyFinal$ExpHigh %>% select(matches("SL")))
diag(SampleCor2) <- NA

colnames(SampleCor2) <- studyFinal$Metadata$Paper_ID[match(colnames(SampleCor2), studyFinal$Metadata[[SampleIDcol]])]
rownames(SampleCor2) <- studyFinal$Metadata$Paper_ID[match(rownames(SampleCor2), studyFinal$Metadata[[SampleIDcol]])]

annoCol = data.frame(MT_Type = studyFinal$Metadata$kmeans_Final,
                     DV200 = studyFinal$Metadata$DV200,
                     Age = studyFinal$Metadata$Age,
                     Sex = studyFinal$Metadata$Sex,
                     row.names = studyFinal$Metadata$Paper_ID)


annoColors$Age = c("grey90", MoviePalettes$BugsLife[8])

pheatmap(SampleCor2, angle_col = 90, na_col = "white",border_color = NA, clustering_method = "ward.D2",cutree_cols = 3, 
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         show_rownames = T, show_colnames = T,
         annotation_col = annoCol,fontsize_col = 7, fontsize_row = 7,
         annotation_colors = annoColors,
         filename = paste0(ResultsPath, "SampleCorPaperB.pdf"),
         width = 12, height = 10)
closeDev()

#Get cellular composition changes
CellularChanges <- GetCellularChanges()
ggsave(paste0(ResultsPath, "CellularChanges.pdf"), plot = CellularChanges$Plot, device = "pdf", width = 8, height = 10, useDingbats = F)

CellularChangesDF <- sapply(names(CellularChanges$ChangesList), function(Mod){
  DF <- CellularChanges$ChangesList[[Mod]] %>% arrange(Group, pValue) %>%
    mutate(Model = Mod) %>% select(-Color, -CellType)
  }, simplify = F) %>% rbindlist() %>% data.frame() %>% arrange(Model, Group, CellTypeName)

CellularChangesDF$OutPut <- apply(CellularChangesDF %>% select(Estimate, Low, High),1, function(x){
  paste0(round(x[1], digits = 2), " (",
         round(x[2], digits = 2), ", ",
         round(x[3], digits = 2),")")
})

CellularChangesDF %<>% mutate(pValueOut = as.character(pValue))


CellularChangesDFlong <- pivot_longer(CellularChangesDF %>%
                                        select(Model, Group, CellTypeName, Estimate, pValue,
                                               OutPut, pValueOut), cols = c(6:7)) %>% data.frame()

MTtypes = unique(CellularChangesDFlong$Group)
CellularChangesDFlong <- rbind(CellularChangesDFlong, data.frame(Model = unique(CellularChangesDFlong$Model),
                                                                 CellTypeName = "Measure",
                                                                 Estimate = 0, pValue = min(CellularChangesDFlong$pValue) - 0.0021,
                                                                 Group = c(rep(MTtypes[1], 6),
                                                                           rep(MTtypes[2],6)),
                                                                 name = rep(c(rep("OutPut", 3),
                                                                              rep("pValueOut", 3)),2),
                                                                 value = rep(c(rep("Group Coeficient (95%C.I.)", 3),
                                                                            rep("pValue", 3)),2)))

CellularChangesDFplot <- sapply(unique(CellularChangesDFlong$Model), function(Mod){
  DF <- CellularChangesDFlong %>% filter(Model == Mod) %>% arrange(Group, desc(pValue)) 
  
  DF$CellTypeName <- factor( DF$CellTypeName, levels = unique( DF$CellTypeName))
  Cov = strsplit(Mod, "\\+")[[1]]
  Cov <- sapply(Cov, function(x){
    gsub(" |_Genes", "", x)
  })
  DF %<>% filter(!CellTypeName %in% Cov)
  
  DF$Estimate[DF$CellTypeName == "Measure"] <- NA
  Title = gsub("_binned|2", "", Mod)
  Title = gsub("Genes", "MGP", Title)
  Title = gsub("MitoType", "MT_Type", Title)
  
  Plot <- ggplot(DF, aes(name, CellTypeName, fill = Estimate)) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "", y = "", title = Title) +
    geom_tile() +
    geom_text(aes(label = value)) +
    scale_fill_gradient2(low = MoviePalettes$BugsLife[6], mid = "white",
                         high =  MoviePalettes$BugsLife[4], midpoint = 0, na.value = "white", limits = c(-0.35, 0.35)) +
    geom_vline(xintercept = 1.7) +
    geom_hline(yintercept = length(unique(DF$CellTypeName)) - 0.5) +
    facet_wrap(~Group, ncol = 2)
  
  Plot <- ChangeFacetLabels(Plot, FillCol = c("#ADEF8B", "#BB9BEA")) %>% ggarrange()
  return(Plot)
}, simplify = F)

ggarrange(plotlist = CellularChangesDFplot, nrow = 3, common.legend = T, legend = "right")
ggsave(paste0(ResultsPath, "CellularChangesDF.pdf"), device = "pdf", height = 9, width = 10)


## Run DE analysis
PathwaysList2 <- sapply(names(PathwaysList), function(GeneSetType){
  GeneSet <- PathwaysList[[GeneSetType]]
  names(GeneSet) <- sapply(names(GeneSet), function(Ptwy){
    gsub("^GO_|^KEGG_|^HALL.*_", "", Ptwy) %>% gsub(" ", "_", .)
  })
  names(GeneSet) <- paste0(GeneSetType, "__", names(GeneSet))
  GeneSet
}, simplify = F)

PathwaysCombined <- c(PathwaysList2$GObp,
                      PathwaysList2$GOmf,
                      PathwaysList2$GOcc,
                      PathwaysList2$MitoCarta,
                      PathwaysList2$KEGG)

Model = as.formula(~MitoType + Sex + Age + Cohort2 + GabaVIPReln_Genes + Oligo_Genes + PMI_binned)

# DESeqOut0 <- DESeq2RUN(data =  studyFinal$countMatrix,
#                       Meta = studyFinal$Metadata, model = Model)
# 
# DESeqResults_nonMT0 <- GetDESeq2Results(DESeqOut0, coef = "MitoTypeNonMT_PD")
# DESeqResults_MTPD0 <- GetDESeq2Results(DESeqOut0, coef = "MitoTypeMT_PD")
# 
# 
# Enrichment_NonMT0 <- RunEnrich(DESeqResults_nonMT0, method = "ErmineJ")
# Enrichment_MTPD0 <- RunEnrich(DESeqResults_MTPD0,  method = "ErmineJ")

#Normalizing to protein coding genes

DESeqOut <- DESeq2RUN(data =  studyFinal$countMatrix,
                      Meta = studyFinal$Metadata, model = Model,
                      ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
                                                                            filter(gene_type == "protein_coding") %>%
                                                                            .$gene_id2))

DESeqResults_nonMT <- GetDESeq2Results(DESeqOut, coef = "MitoTypenCI.PD")
DESeqResults_MTPD <- GetDESeq2Results(DESeqOut, coef = "MitoTypeCI.PD")


Enrichment_NonMT <- RunEnrich(DESeqResults_nonMT, method = "ErmineJ")
Enrichment_MTPD <- RunEnrich(DESeqResults_MTPD,  method = "ErmineJ")



#Correcting for oligodendrocytes only
Model2 = as.formula(~MitoType + Sex + Age + Cohort2 +  Oligo_Genes + PMI_binned)

DESeqOut2 <- DESeq2RUN(data =  studyFinal$countMatrix,
                       Meta = studyFinal$Metadata, model = Model2,
                       ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
                                                                             filter(gene_type == "protein_coding") %>%
                                                                             .$gene_id2)) 

DESeqResults_nonMT2 <- GetDESeq2Results(DESeqOut2, coef = "MitoTypenCI.PD")
DESeqResults_MTPD2 <- GetDESeq2Results(DESeqOut2, coef = "MitoTypeCI.PD")

Enrichment_NonMT2 <- RunEnrich(DESeqResults_nonMT2,  method = "ErmineJ")
Enrichment_MTPD2 <- RunEnrich(DESeqResults_MTPD2, method = "ErmineJ")


PlotRanks <- function(DESeqResults, Pathway = PathwaysList$KEGG$KEGG_OXIDATIVE_PHOSPHORYLATION, PathwayName = "KEGG_OXPHOS"){
  temp <- DESeqResults %>% select(gene_name, log2FoldChange, padj) %>% mutate(pPadj = -log10(padj))
  temp$Direction <- sapply(temp$log2FoldChange, function(x){
    if(x < 0){
      "Downregulated"
    } else {
      "Upregulated"
    }
  })
  
  temp$Pathway <- sapply(temp$gene_name, function(Gene){
    if(Gene %in% Pathway){
      PathwayName
    } else {
      "Other"
    }
  })
  
  nGenes <- temp %>% group_by(Direction, Pathway) %>% summarise(n = n()) %>% data.frame() %>%
    mutate(Direction_Pathway = paste(Direction, Pathway, sep = "_"))
  
  temp %<>% mutate(Direction_Pathway = paste(Direction, Pathway, sep = "_"))
  temp$Pathway2 <- nGenes$n[match(temp$Direction_Pathway, nGenes$Direction_Pathway)]
  temp %<>% mutate(Pathway2 = paste0(Pathway, "\nn = ", Pathway2))
  
  ggplot(temp, aes(Pathway2, pPadj, color = Pathway, fill = Pathway)) +
    theme_minimal() +
    labs(x = "", y = "-log10(padj)") +
    geom_boxplot(width = 0.12, alpha = 0.5, color = "black", position= position_nudge(x=-.1)) +
    stat_halfeye(adjust = 0.5, justification = -0.15, .width = 0, width = 0.8, position= position_nudge(x= -.1), point_colour = NA) +
    scale_fill_manual(values = c("cornflowerblue", "grey80")) +
    scale_color_manual(values = c("cornflowerblue", "grey80")) +
    geom_hline(yintercept = -log10(0.05), lty = "dashed") +
    facet_wrap(~Direction, scales = "free_x")
      
}

OxPhosPlot <- lapply(list(MT2 = DESeqResults_MTPD2, MT = DESeqResults_MTPD,
                          nonMT2 = DESeqResults_nonMT2,nonMT = DESeqResults_nonMT), function(DESeqResults){
                            PlotRanks(DESeqResults)
                          })
ggarrange(plotlist = OxPhosPlot, common.legend = T, labels = names(OxPhosPlot))
ggsave(paste0(ResultsPath, "OxPhosRanks.pdf"), device = "pdf", width = 8, height = 4, useDingbats = F)

#This part is commented since it is not used for the paper
# #Repeating the analysis with the two mtPD groups separated
# Model4 = as.formula(~MT_GroupPD3 + Sex + Age + Cohort2 + Astrocyte_Genes + GabaVIPReln_Genes + Oligo_Genes + PMI_binned)
# 
# DESeqOut4 <- DESeq2RUN(data =  studyFinal$countMatrix,
#                        Meta = studyFinal$Metadata, model = Model4,
#                        ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
#                                                                              filter(gene_type == "protein_coding") %>%
#                                                                              .$gene_id2)) 
# 
# DESeqResults_nonMT4 <- GetDESeq2Results(DESeqOut4, coef = "MT_GroupPD3NonMT_PD")
# DESeqResults_MTa <- GetDESeq2Results(DESeqOut4, coef = "MT_GroupPD3MT_PDa")
# DESeqResults_MTb <- GetDESeq2Results(DESeqOut4, coef = "MT_GroupPD3MT_PDb")
# 
# Enrichment_MTa <- RunEnrich(DESeqResults_MTa,  method = "ErmineJ")
# Enrichment_MTb <- RunEnrich(DESeqResults_MTb,  method = "ErmineJ")
# 
# #Without correcting for astrocytes
# Model5 = as.formula(~MT_GroupPD3 + Sex + Age + Cohort2 + Oligo_Genes + PMI_binned)
# 
# DESeqOut5 <- DESeq2RUN(data =  studyFinal$countMatrix,
#                        Meta = studyFinal$Metadata, model = Model5,
#                        ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
#                                                                              filter(gene_type == "protein_coding") %>%
#                                                                              .$gene_id2)) 
# 
# DESeqResults_nonMT5 <- GetDESeq2Results(DESeqOut5, coef = "MT_GroupPD3NonMT_PD")
# DESeqResults_MTa2 <- GetDESeq2Results(DESeqOut5, coef = "MT_GroupPD3MT_PDa")
# DESeqResults_MTb2 <- GetDESeq2Results(DESeqOut5, coef = "MT_GroupPD3MT_PDb")
# 
# Enrichment_MTa2 <- RunEnrich(DESeqResults_MTa2,  method = "ErmineJ")
# Enrichment_MTb2 <- RunEnrich(DESeqResults_MTb2,  method = "ErmineJ")


#Combine the outputs of all models in both groups
CombinedResultsDE <-  merge(DESeqResults_nonMT %>% select(-(lfcSE)), DESeqResults_MTPD %>%
                              select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                            by = "EnsemblID", suffixes = c("_NonMTmodel1", "_MTPDmodel1"))


CombinedResultsDE <-  merge(CombinedResultsDE, merge(DESeqResults_nonMT2 %>%
                                                       select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                     DESeqResults_MTPD2 %>%
                                                       select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                     by = "EnsemblID", suffixes = c("_NonMTmodel2", "_MTPDmodel2")), by = "EnsemblID") %>%
  select(EnsemblID, gene_name, gene_type, baseMean,
         matches("stat"), matches("log2"), matches("pval"), matches("padj"))



sapply(ls(pat = "DESeqResults_.*MT.?[0-9]?$|CombinedResults"), function(x){
  write.table(get(x), paste0(ResultsPath, x, ".tsv"), row.names = F, col.names = T, sep = "\t")
})

#Save enrichmet results
sapply(ls(pat = "Enrichment_"), function(EnrGrp){
  Data = get(EnrGrp)
  #browser()
  DataDF <- lapply(Data$ErmineJ, function(List){
    List$results %>% data.frame() %>% select(-GeneMembers) %>% filter(CorrectedPvalue < 0.05)
  }) %>% rbindlist() %>% data.frame()
  write.table(DataDF, paste(ResultsPath, EnrGrp, ".tsv"), col.names = T, row.names = F, sep = "\t")
})


#Look at the enrichment for Lysosomal and OxPhos pathways

EnrichPathWays = sapply(ls(pat = "Enrichment_.*[0-3]?$"), function(Mod){
  Data = get(Mod)
  DataUp <- Data$ErmineJ$EnrichUp$results %>% filter(ID %in% c("LYSOSOME", "OXPHOS", "OXIDATIVE_PHOSPHORYLATION"), GeneSet %in% c("KEGG", "MitoCarta"))
  DataDown <- Data$ErmineJ$EnrichDown$results %>% filter(ID %in% c("LYSOSOME", "OXPHOS", "OXIDATIVE_PHOSPHORYLATION"), GeneSet %in% c("KEGG", "MitoCarta"))
  rbind(DataUp, DataDown) %>%
    select(-GeneMembers) %>% data.frame() %>%
    arrange(Pval) %>% filter(!duplicated(ID)) %>%
    mutate(GroupMod = gsub("Enrichment_", "", Mod),
           pLog10p = -1*log10(Pval))
}, simplify = F) %>% rbindlist %>% data.frame() %>% .[!grepl("a|b|0", .$GroupMod),]

EnrichPathWays$Group <- sapply(EnrichPathWays$GroupMod, function(x){
  if(grepl("Non", x)){
    "NonMT_PD"
  } else {
    "MT_PD"
  }
})

EnrichPathWays$Mod <- sapply(EnrichPathWays$GroupMod, function(x){
  if(grepl("2", x)){
    "Model_2"
  } else if(grepl("3", x)) {
    "Model_3"
  } else {
    "Model_1"
  }
})

EnrichPathWays$pLog10p <- apply(EnrichPathWays %>% select(Direction, pLog10p), 1, function(x){
  if(x[1] == "Down"){
    -1*as.numeric(x[2])
  } else {
    as.numeric(x[2])
  }
})

EnrichPathWays$Signif <- sapply(EnrichPathWays$CorrectedPvalue, function(x){
  if(x < 0.05){
    "Yes"
  } else {
    "No"
  }
})


EnrichPathWays$ID <- sapply(EnrichPathWays$ID, function(x){
  if(x == "LYSOSOME"){
    paste0(x, " (KEGG)")
  } else if(x == "OXPHOS"){
    paste0(x, " (MitoCarta)")
  } else {
    "OXPHOS (KEGG)"
  }
})

PatwayPlot <- ggplot(EnrichPathWays, aes(Mod, pLog10p, fill = Direction, alpha = Signif)) +
  theme_bw(base_size = 14, base_family = "serif") +
  labs(y = "-log10(pValue)", x = "") +
  geom_bar(stat = "identity") +
  facet_grid(Group~ID, scales = "free") +
  scale_fill_manual(values = MoviePalettes$BugsLife[c(6, 4)]) +
  scale_alpha_manual(values = c(0.5, 1), name = "FDR significant", guide = "none") +
  #geom_hline(yintercept = c(log10(0.05),  -log10(0.05)), lty = "dashed") +
  coord_flip() +
  scale_y_continuous(breaks = pretty(EnrichPathWays$pLog10, n = 6),
                     labels = abs(pretty(EnrichPathWays$pLog10p, n = 6)))

ChangeFacetLabels(PatwayPlot, FillCol = c("#ADEF8B", "#BB9BEA", rep("grey90", 3))) %>% ggarrange()
ggsave(paste0(ResultsPath, "LysoOxphosEnrich.pdf"), device = "pdf", width = 10, height = 4)



# PlotCountsMeta(Gene = "MALAT1", CompareVar = "RIN",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "GroupPD")$Plot3
# PlotCountsMeta(Gene = "SIRT1", CompareVar = "OxPhosES",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "GroupPD")$Plot3
# 
# GM2A <- PlotCountsMeta(Gene = "GM2A", CompareVar = "Oligo_Genes",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "Cohort2")
# LAMP1 <- PlotCountsMeta(Gene = "LAMP1", CompareVar = "PercentPositive",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "Cohort2")
# TPP1 <- PlotCountsMeta(Gene = "TPP1", CompareVar = "Oligo_Genes",  dds = DESeqOut2, Grp = ClusterCol, colorBy2 = "Cohort2")
# 
# 
# PlotGeneByGroup("LAMP1")
# PlotGeneByGroup ("GM2A")
# PlotGeneByGroup ("TPP1")

LysoGenes <- CombinedResultsDE %>% filter(gene_name %in% c("LAMP1", "TPP1", "GM2A")) %>% select(matches("gene_name|padj|log2")) %>%
  pivot_longer(-gene_name, names_pattern = "(.*)_(.*)", names_to = c("measure", "Group")) %>%
  pivot_wider(names_from = "measure") %>% data.frame()

LysoGenes$Model <- sapply(LysoGenes$Group, function(x){
  gsub("NonMT|MTPD", "", x)
})

LysoGenes$Group <- sapply(LysoGenes$Group, function(x){
  gsub("model.", "", x)
})

LysoGenes %<>% mutate(log2FoldChange = round(log2FoldChange, digits = 2),
                      padj = signif(padj, digits = 2)) %>% arrange(gene_name, Group, Model)

write.table(LysoGenes, paste0(ResultsPath, "LysoGenes.tsv"), sep = "\t", row.names = F, col.names = T)

save.image(paste0(ResultsPath,"RNAseq.RData"))
saveRDS(studyFinal, paste0(ResultsPath, "studyFinal.Rds"))
saveRDS(PCAresults, paste0(ResultsPath, "PCAresults.Rds"))

AnnoDF = data.frame(MT_GroupFinal = studyFinal$Metadata$kmeans_Final,
                    PercentPositive = studyFinal$Metadata$PercentPositive,
                    DV200 = studyFinal$Metadata$DV200,
                    OxPhos = studyFinal$Metadata$OxPhosES,
                    Age = studyFinal$Metadata$Age,
                    Cohort = studyFinal$Metadata$Cohort2,
                    row.names = studyFinal$Metadata$Paper_ID)
       



PlotList <- sapply(grep("CI |CII |CIII |CIV |CV ", names(PathwaysList$MitoCarta), value = T), function(ptwy){
  PlotObj <- Enrichheatmap(GeneSet = "MitoCarta", Ptwy = ptwy,
                           show_rownames = F, show_colnames = F,
                           cutree_rows = 1, OrderCols = "Groups")
  geneNum = length(PlotObj$tree_row$labels)
  Plot <- cowplot::gtable_remove_grobs(PlotObj$gtable, names = c("col_annotation",
                                                                "col_annotation_names",
                                                                "main"))
  list(PlotObj = PlotObj, geneNum = geneNum, Plot = Plot)
}, simplify = F)


GeneOrder <- sapply(names(PlotList), function(ptwy){
  Data = PlotList[[ptwy]]$PlotObj$tree_row
  Complex = strsplit(ptwy, " ")[[1]][1] %>% gsub("C", "Complex ", .)
  Type = strsplit(ptwy, " ")[[1]][2]
  data.frame(GeneNames = Data$labels[Data$order], Complex = Complex, Type = Type)
}, simplify = F) 

for(i in c(1, 2, 6:8)){
  GeneOrder[[i]] <- GeneOrder[[i]][nrow(GeneOrder[[i]]):1,]
}

GeneOrderDF <- rbindlist(GeneOrder) %>% data.frame()
GeneOrderDF[GeneOrderDF$GeneNames %in% (GeneOrderDF %>%
                                          filter(duplicated(GeneNames)) %>%
                                          .$GeneNames),]$Complex <- "Mixed"
GeneOrderDF %<>% filter(!duplicated(GeneNames))
rownames(GeneOrderDF) <- GeneOrderDF$GeneNames
GeneOrderDF %<>% arrange(desc(Type), Complex)


temp <- studyFinal$ExpHigh %>% filter(GeneSymbol %in% GeneOrderDF$GeneNames) %>%
  select(matches("SL"))

rownames(temp) <- studyFinal$ExpHigh %>%
  filter(GeneSymbol %in%  GeneOrderDF$GeneNames) %>% .$GeneSymbol


ColOrder <- studyFinal$Metadata %>% arrange(Cohort2, MitoType, OxPhosES) %>% .$RNAseq_id_ParkOme2
temp %<>% select(all_of(ColOrder))
colnames(temp) <- studyFinal$Metadata$Paper_ID[match(colnames(temp),
                                                     studyFinal$Metadata$RNAseq_id_ParkOme2)] 
temp <- temp[GeneOrderDF$GeneNames,]


tempNum <- GeneOrderDF %>% filter(Type == "subunits") %>%  .$Complex %>% table()
annoColors$Complex <- c("Complex I" =  MoviePalettes$MoonRiseKingdomColors[2],
                        "Complex II" =  MoviePalettes$MoonRiseKingdomColors[10],
                        "Complex III" =  MoviePalettes$MoonRiseKingdomColors[4],
                        "Complex IV" =  MoviePalettes$MoonRiseKingdomColors[5],
                        "Complex V" =  MoviePalettes$MoonRiseKingdomColors[7],
                        "Mixed" = "black")
annoColors$Type <- c("subunits" = MoviePalettes$AmericanBeauty[5],
                     "assembly" = MoviePalettes$AmericanBeauty[6])

pheatmap(temp, scale = "row", angle_col = 90, na_col = "white",
         border_color = NA, cluster_cols = F, cluster_rows = F,
         gaps_row = c(tempNum[1], sum(tempNum[1:2]),
                      sum(tempNum[1:3]), sum(tempNum[1:4]),
                      sum(tempNum[1:5]), sum(tempNum[1:5])),
         gaps_col = AnnoDF %>% filter(Cohort != "Norway") %>% nrow, 
         cutree_rows = 3, cutree_cols = 3,
         show_colnames = F, show_rownames = F,
         clustering_method = "ward.D2", clustering_callback = callback, 
         color = colorRampPalette(c("darkblue", "grey10", "gold2"))(999),
         annotation_col = AnnoDF,
         annotation_row = GeneOrderDF %>% select(-GeneNames),
         annotation_colors = annoColors,
         main = "OxPhos gene expression",
         filename = paste0(ResultsPath, "AllOxphosHeatMapFinal.pdf"), width = 8, height = 10)

#With gene and sample names
pheatmap(temp %>% filter(rownames(.) %in% (GeneOrderDF %>% filter(Type == "assembly") %>% .$GeneNames)),
         scale = "row", angle_col = 90, na_col = "white",
         border_color = NA, cluster_cols = F, cluster_rows = F,
         gaps_col = AnnoDF %>% filter(Cohort != "Norway") %>% nrow, 
         cutree_rows = 3, cutree_cols = 3,
         show_colnames = T, show_rownames = T,
         clustering_method = "complete", clustering_callback = callback, 
         color = colorRampPalette(c("darkblue", "grey10", "gold2"))(999),
         annotation_col = AnnoDF,
         annotation_row = GeneOrderDF %>% select(-GeneNames),
         annotation_colors = annoColors, fontsize_row = 8, fontsize_col = 8,
         main = "OxPhos gene expression",
         filename = paste0(ResultsPath, "AssemplyOxphosHeatMapIDs.pdf"), width = 14, height = 9)

pheatmap(temp %>% filter(!rownames(.) %in% (GeneOrderDF %>% filter(Type == "assembly") %>% .$GeneNames)),
         scale = "row", angle_col = 90, na_col = "white",
         border_color = NA, cluster_cols = F, cluster_rows = F,
         gaps_row = c(tempNum[1], sum(tempNum[1:2]),
                      sum(tempNum[1:3]), sum(tempNum[1:4]),
                      sum(tempNum[1:5])),
         gaps_col = AnnoDF %>% filter(Cohort != "Norway") %>% nrow, 
         cutree_rows = 3, cutree_cols = 3,
         show_colnames = T, show_rownames = T,
         clustering_method = "complete", clustering_callback = callback, 
         color = colorRampPalette(c("darkblue", "grey10", "gold2"))(999),
         annotation_col = AnnoDF,
         annotation_row = GeneOrderDF %>% select(-GeneNames),
         annotation_colors = annoColors, fontsize_row = 8, fontsize_col = 8,
         main = "OxPhos gene expression",
         filename = paste0(ResultsPath, "SubunitsOxphosHeatMapIDs.pdf"), width = 14, height = 11)
