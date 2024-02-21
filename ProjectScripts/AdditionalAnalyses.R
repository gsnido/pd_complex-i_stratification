#Look at the sample clustering in each cohort separetely
MetaBBB <- studyFinal$Metadata %>% filter(Cohort2 == "Barcelona")

SampleCorBBB <- cor(studyFinal$ExpHigh %>% select(MetaBBB$RNAseq_id_ParkOme2))
diag(SampleCorBBB) <- NA

colnames(SampleCorBBB) <- MetaBBB$Paper_ID[match(colnames(SampleCorBBB), MetaBBB[[SampleIDcol]])]
rownames(SampleCorBBB) <- MetaBBB$Paper_ID[match(rownames(SampleCorBBB), MetaBBB[[SampleIDcol]])]

annoColBBB = data.frame(MT_GroupFinal = MetaBBB$kmeans_Final,
                        DV200 = MetaBBB$DV200,
                        OxPhos_Genes = MetaBBB$OxPhosES,
                        Age = MetaBBB$Age,
                        Sex = MetaBBB$Sex,
                        row.names = MetaBBB$Paper_ID)

                    


pheatmap(SampleCorBBB, angle_col = 90, na_col = "white",border_color = NA, clustering_method = "ward.D2",
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         show_rownames = T, show_colnames = T,
         annotation_col = annoColBBB,fontsize_col = 7, fontsize_row = 7,
         annotation_colors = annoColors, main = "ESP",
         filename = paste0(ResultsPath, "SampleCorBBB.pdf"),
         width = 12, height = 10)



MetaPW <- studyFinal$Metadata %>% filter(Cohort2 == "Norway")

SampleCorPW <- cor(studyFinal$ExpHigh %>% select(MetaPW$RNAseq_id_ParkOme2))
diag(SampleCorPW) <- NA

colnames(SampleCorPW) <- MetaPW$Paper_ID[match(colnames(SampleCorPW), MetaPW[[SampleIDcol]])]
rownames(SampleCorPW) <- MetaPW$Paper_ID[match(rownames(SampleCorPW), MetaPW[[SampleIDcol]])]

annoColPW = data.frame(MT_GroupFinal = MetaPW$kmeans_Final,
                        DV200 = MetaPW$DV200,
                        OxPhos_Genes = MetaPW$OxPhosES,
                        Age = MetaPW$Age,
                        Sex = MetaPW$Sex,
                        row.names = MetaPW$Paper_ID)



pheatmap(SampleCorPW, angle_col = 90, na_col = "white",border_color = NA, clustering_method = "ward.D2", 
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         show_rownames = T, show_colnames = T,
         annotation_col = annoColPW,fontsize_col = 7, fontsize_row = 7,
         annotation_colors = annoColors, main = "NOR",
         filename = paste0(ResultsPath, "SampleCorPW.pdf"),
         width = 12, height = 10)




#Cellular changes
CellularChanges <- GetCellularChanges(BasicModel = "~MitoType + Sex + Age  +  PMI_binned", Meta = MetaBBB)
ggsave(paste0(ResultsPath, "CellularChangesBBB.pdf"), plot = CellularChanges$Plot, device = "pdf", width = 8, height = 10, useDingbats = F)

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
  
  Plot <- ChangeFacetLabels(Plot, FillCol = c(MoviePalettes$BugsLife[2], "orange")) %>% ggarrange()
  return(Plot)
}, simplify = F)

ggarrange(plotlist = CellularChangesDFplot, nrow = 3, common.legend = T, legend = "right")
ggsave(paste0(ResultsPath, "CellularChangesDF_BBB.pdf"), device = "pdf", height = 9, width = 10)

#Repeat in PW
CellularChanges <- GetCellularChanges(BasicModel = "~MitoType + Sex + Age  +  PMI_binned", Meta = MetaPW, Ylim = c(-0.7, 0.7))
ggsave(paste0(ResultsPath, "CellularChangesPW.pdf"), plot = CellularChanges$Plot, device = "pdf", width = 8, height = 10, useDingbats = F)

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
                         high =  MoviePalettes$BugsLife[4], midpoint = 0, na.value = "white", limits = c(-0.5, 0.5)) +
    geom_vline(xintercept = 1.7) +
    geom_hline(yintercept = length(unique(DF$CellTypeName)) - 0.5) +
    facet_wrap(~Group, ncol = 2)
  
  Plot <- ChangeFacetLabels(Plot, FillCol = c(MoviePalettes$BugsLife[2], "orange")) %>% ggarrange()
  return(Plot)
}, simplify = F)

ggarrange(plotlist = CellularChangesDFplot, nrow = 3, common.legend = T, legend = "right")
ggsave(paste0(ResultsPath, "CellularChangesDF_PW.pdf"), device = "pdf", height = 9, width = 10)


#Repeating DE analysis for for BBB samples only 
Model_BBB = as.formula(~MitoType + Sex + Age + Astrocyte_Genes + GabaVIPReln_Genes + Oligo_Genes + PMI_binned)
DESeqOut_BBB <- DESeq2RUN(data =  studyFinal$countMatrix,
                          Meta = MetaBBB, model = Model_BBB,
                          ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
                                                                                filter(gene_type == "protein_coding") %>%
                                                                                .$gene_id2)) 

DESeqResults_BBB_nonMT <- GetDESeq2Results(DESeqOut_BBB, coef = "MitoTypenCI.PD")
DESeqResults_BBB_MTPD <- GetDESeq2Results(DESeqOut_BBB, coef = "MitoTypeCI.PD")

Enrichment_BBB_nonMT <- RunEnrich(DESeqResults_BBB_nonMT,  method = "ErmineJ")
Enrichment_BBB_MTPD <- RunEnrich(DESeqResults_BBB_MTPD,  method = "ErmineJ")


#Repeat without correcting for GabaVIPReln MGP
Model_BBB2 = ~MitoType + Sex + Age +  Oligo_Genes + Astrocyte_Genes +  PMI_binned
DESeqOut_BBB2 <- DESeq2RUN(data =  studyFinal$countMatrix, 
                            Meta = MetaBBB, model = Model_BBB2,
                            ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
                                                                                  filter(gene_type == "protein_coding") %>%
                                                                                  .$gene_id2)) 

DESeqResults_BBB_nonMT2 <- GetDESeq2Results(DESeqOut_BBB2, coef = "MitoTypenCI.PD")
DESeqResults_BBB_MTPD2 <- GetDESeq2Results(DESeqOut_BBB2, coef = "MitoTypeCI.PD")

Enrichment_BBB_nonMT2 <- RunEnrich(DESeqResults_BBB_nonMT2,  method = "ErmineJ")
Enrichment_BBB_MTPD2 <- RunEnrich(DESeqResults_BBB_MTPD2,  method = "ErmineJ")

#Correcting for oligodendrocytes only
Model_BBB3 = ~MitoType + Sex + Age +  Oligo_Genes + PMI_binned

DESeqOut_BBB3 <- DESeq2RUN(data =  studyFinal$countMatrix, 
                           Meta = MetaBBB, model = Model_BBB3,
                           ContGenes = rownames(studyFinal$countMatrix) %in% c(geneNames %>%
                                                                                 filter(gene_type == "protein_coding") %>%
                                                                                 .$gene_id2)) 

DESeqResults_BBB_nonMT3 <- GetDESeq2Results(DESeqOut_BBB3, coef = "MitoTypenCI.PD")
DESeqResults_BBB_MTPD3 <- GetDESeq2Results(DESeqOut_BBB3, coef = "MitoTypeCI.PD")

Enrichment_BBB_nonMT3 <- RunEnrich(DESeqResults_BBB_nonMT3,  method = "ErmineJ")
Enrichment_BBB_MTPD3 <- RunEnrich(DESeqResults_BBB_MTPD3,  method = "ErmineJ")


CombinedResultsDE_BBB <-  merge(DESeqResults_BBB_nonMT2 %>% select(-(lfcSE)), DESeqResults_BBB_MTPD %>%
                                  select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                by = "EnsemblID", suffixes = c("_NonMTmodel1", "_MTPDmodel1"))


CombinedResultsDE_BBB <-  merge(CombinedResultsDE_BBB, merge(DESeqResults_BBB_nonMT2 %>%
                                                               select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                             DESeqResults_BBB_MTPD2 %>%
                                                               select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                             by = "EnsemblID", suffixes = c("_NonMTmodel2", "_MTPDmodel2")), by = "EnsemblID")

CombinedResultsDE_BBB <- merge(CombinedResultsDE_BBB, merge(DESeqResults_BBB_nonMT3 %>%
                                                              select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                            DESeqResults_BBB_MTPD3 %>%
                                                              select(EnsemblID, log2FoldChange, stat, pvalue, padj),
                                                            by = "EnsemblID", suffixes = c("_NonMTmodel3", "_MTPDmodel3")), by = "EnsemblID") %>%
  select(EnsemblID, gene_name, gene_type, baseMean,
         matches("stat"), matches("log2"), matches("pval"), matches("padj"))




write.table(CombinedResultsDE_BBB, paste0(ResultsPath, "CombinedResults_BBB.tsv"), row.names = F, col.names = T, sep = "\t")


#Save enrichment results
sapply(ls(pat = "Enrichment_BBB"), function(EnrGrp){
  Data = get(EnrGrp)
  #browser()
  DataDF <- lapply(Data$ErmineJ, function(List){
    List$results %>% data.frame() %>% select(-GeneMembers) %>% filter(CorrectedPvalue < 0.05)
  }) %>% rbindlist() %>% data.frame()
  write.table(DataDF, paste(ResultsPath, EnrGrp, ".tsv"), col.names = T, row.names = F, sep = "\t")
})
