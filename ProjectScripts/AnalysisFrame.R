#BiocManager::install("devtools")
library(devtools)
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/generalFunc.R?raw=T")

packageF("org.Hs.eg.db")
packageF("GenomicFeatures")
packageF("AnnotationDbi")
packageF("pheatmap")
packageF("DESeq2")
packageF("cluster")
packageF("fgsea")
packageF("tibble")
packageF("ggbeeswarm")
packageF("ggdist")
packageF("grid")
packageF("gridExtra")
packageF("ggnewscale")
#install_github("https://github.com/PavlidisLab/ermineR")
library(ermineR)
packageF("factoextra")
packageF("NbClust")
packageF("betareg")

SampleIDcol = "RNAseq_id_ParkOme2"

AnnoLoc = "Annotations"
AssemblyFilename = "gencode.v35.annotation.gff3.gz"

#Load custom functions
source("ProjectScripts/ProjectFunctions.R")


#Get genome annotations
source("ProjectScripts/Annotations.R")

if(!file.exists(geneNameFile)){
  geneNames <- rtracklayer::import.gff3(GTF_file) %>%
    data.frame() %>% filter(type == "transcript")  %>%
    select(seqnames, start, end, ID,
           gene_id, gene_name, gene_type, transcript_type)
  saveRDS(geneNames, geneNameFile)
} else {
  geneNames <- readRDS(geneNameFile)
}

geneNames$gene_id2 <- sapply(geneNames$gene_id, function(x){
  strsplit(x, "\\.")[[1]][1]
})


ResultsPath = paste0("GeneralResults_", strsplit(AssemblyFilename, "\\.")[[1]][2])
if(!ResultsPath %in% list.dirs(full.names = F, recursive = T)){
  dir.create(ResultsPath)
}
ResultsPath = paste0(ResultsPath, "/")

Study = "GeneralAnalysis"

#Get metadata
Metadata <- readRDS("Data/Metadata.Rds")

Metadata %<>% mutate(NeuExpRegion = OrgRegion,
                     Series_sample_id = Filename,
                     Study = Study)


IHCall <- read.table("Data/IHCdataPaper.txt", header = T, sep = "\t") %>% mutate(kmeans_Final_combined = kmeans_Final)
IHCall$kmeans_Final_combined <- sapply(IHCall$kmeans_Final, function(x){
  strsplit(x, " ")[[1]][1]
}) 

IHCall$kmeans_Final <- factor(IHCall$kmeans_Final, levels = c("Control", "nCI-PD", "CI-PD mild", "CI-PD severe") )
IHCall$kmeans_Final_combined <- factor(IHCall$kmeans_Final_combined, levels = c("Control", "nCI-PD", "CI-PD"))


#Add MT cluster information to the metadata
Metadata <- merge(Metadata, IHCall, by.x = "Biobank_ID", by.y = "BiobankID", all.x = T)
# Metadata <- merge(Metadata, IHCpd %>% select(BiobankID, matches("MT_Group")),
#                   by.x = "Biobank_ID", by.y = "BiobankID", all.x = T)

#Metadata %<>% mutate(PercentPositive = 100*PropPos)
Metadata %<>% mutate(PercentPositive = PFC_median_Irene) %>% filter(!is.na(PercentPositive))




#Plot the relevant PD clusters
sapply(grep("MT_Group|kmeans", names(Metadata), value = T), function(ClusterCol){
  Data = Metadata %>% mutate(MitoType = Metadata[[ClusterCol]])
  SubjNumber <- Data %>% group_by(MitoType) %>%
    summarise(n = n()) %>% data.frame() %>%
    mutate(MTgroupPlot = paste0(MitoType, "\n(n=", n, ")"))
  
  Data$MitoType_plot <- SubjNumber$MTgroupPlot[match(Data$MitoType, SubjNumber$MitoType)] 
  
  Data$MitoType_plot <- factor(Data$MitoType_plot, levels = SubjNumber$MTgroupPlot)
  
  PlotCohortCombined <- ggplot(Data, aes(GroupPD, PercentPositive, color = MitoType_plot)) +
    theme_classic() +
    theme(legend.text = element_text(size=6)) +
    labs(x="", y = "", title = "Combined data") + 
    geom_quasirandom() +
    scale_color_manual(values = c("grey", "orange", MoviePalettes$BugsLife[c(3,2)]),
                       name = "Cluster")
  
  ggarrange(PlotClusters(Data = Data, GroupCol = "GroupPD",  ClusterColumn = ClusterCol, title = NULL,
                         colors = c("grey", "orange", MoviePalettes$BugsLife[c(3,2)]),
                         showLegend = F),
            PlotCohortCombined, ncol = 2, widths = c(1, 1.1))
  ggsave(paste0(ResultsPath, "IHCplotCombined_", ClusterCol, ".pdf"), device = "pdf", width = 6,
         height = 2.5, dpi = 300, useDingbats = F, scale = 1.3)
  
})



#Remove low quality RNA sample  and samples without RNAseq data
Metadata %<>% filter(DV200 > 70)

# Get counts
countMatrix <- readRDS("Data/CountMatrix.Rds")

#Keep only samples that have RNAseq data
Metadata <- Metadata[Metadata[[SampleIDcol]] %in% colnames(countMatrix),]

# Match count matrix names to metadata names
countMatrix <- countMatrix[,Metadata[[SampleIDcol]]] 

#Remove genes with maximal count < 2 and mitochondrial genes
Max <- apply(countMatrix, 1, max)
mitoGenes <- geneNames %>% filter(seqnames == "chrM")

countMatrix <- countMatrix[Max > 10,]

countMito <- countMatrix[rownames(countMatrix) %in% mitoGenes$gene_id2,]
countMitoSum <- countMito %>% apply(2, sum)
TotLibSize <- countMatrix  %>% apply(2, sum) 
ShortGeneSum <-  countMatrix %>% data.frame %>% filter(rownames(.) %in% (geneNames %>% mutate(size = end - start) %>% arrange(size) %>% filter(size < 300) %>% .$gene_id2)) %>% apply(2, sum)


MitoCountFiltered <- countMatrix[!rownames(countMatrix) %in% mitoGenes$gene_id2,]
MitoFiltCountSum = apply(MitoCountFiltered, 2, sum)

Metadata$MitoCount <- countMitoSum[match(Metadata[[SampleIDcol]], names(countMitoSum))]
Metadata$TotLibSize <- TotLibSize[match(Metadata[[SampleIDcol]], names(TotLibSize))]
Metadata$ShortGeneSum <- ShortGeneSum[match(Metadata[[SampleIDcol]], names(ShortGeneSum))]


#Look at sample correlation after mitochndrial gene removal 
SampleCor <- cor(MitoCountFiltered)

annoCol = data.frame(Age = Metadata$Age,
                     Batch = Metadata$Batch,
                     Sex = Metadata$Sex,
                     Cohort = Metadata$Cohort2,
                     MT_GroupFinal = Metadata$kmeans_Final,
                     MitoCount = Metadata$MitoCount,
                     LibSize = Metadata$TotLibSize,
                     DV200 = Metadata$DV200,
                     row.names = Metadata[[SampleIDcol]])

annoColors <- list(MT_GroupFinal = c("CI-PD severe" = "#8AD662",
                                     "CI-PD mild" = "#006633",
                                     "nCI-PD" = "#BB9BEA",
                                     "Control" = "grey90"),
                   MT_Type = c("CI-PD" = "#006633",
                               "NonMT_PD" = "#BB9BEA",
                               "Control" = "grey90"),  
                   Sex = c(F = "#FFEED9", M = "#E0B6BE"),
                   Age = c("darkseagreen1", "darkorchid4"),
                   DV200 = c("black", "orange"),
                   Cohort = c("Norway" = MoviePalettes$BugsLife[4],
                              "Barcelona" = MoviePalettes$BugsLife[6]))


pheatmap(SampleCor, angle_col = 90, na_col = "white",border_color = NA,clustering_method = "ward.D2",
         color = colorRampPalette(c("darkblue", "gold2"))(999),
         annotation_col = annoCol,
         annotation_colors = annoColors)


#Look at the most highly expressed genes
TopFiveProportion <- sapply(names(MitoFiltCountSum), function(sbj){
  SubMatrix = data.frame(genes = rownames(MitoCountFiltered), Counts = MitoCountFiltered[,sbj])
  TopFive = SubMatrix %>% arrange(desc(Counts)) %>% head(5)
  TopFive %<>%  mutate(Proportion = Counts/MitoFiltCountSum[sbj])
  Genes <- geneNames[match(TopFive$genes, geneNames$gene_id2),]  %>% select(gene_name, gene_type)
  Genes$Filename = sbj
  temp <- cbind(Genes, TopFive)
  names(temp)[names(temp) == "genes"] <- "ensemblID"
  temp
}, simplify = FALSE) %>% rbindlist()


TopFiveSum <- TopFiveProportion %>% group_by(Filename) %>%
  summarise(TotProp = sum(Proportion)) %>%
  data.frame %>% arrange(TotProp)


TopFiveGeneFreq <- TopFiveProportion %>% group_by(ensemblID) %>%
  summarise(n = n()) %>%
  data.frame

TopFiveGeneFreq <- merge(TopFiveGeneFreq, geneNames,
                         by.x = "ensemblID", by.y = "gene_id2", all.x = T, all.y = F, sort = F) 

TopFiveGeneFreq$ensemblID2 <- sapply(TopFiveGeneFreq$ensemblID,  function(x){
  strsplit(x, "\\.")[[1]][1]
})

TopFiveGeneFreq %<>% mutate(ensemblID2 = paste0(ensemblID2,
                                                " (", gene_name, ", ", n, ")"))
TopFiveGeneFreq %<>% arrange(desc(n))

TopFiveProportion$ensemblID2 <- TopFiveGeneFreq$ensemblID2[match(TopFiveProportion$ensemblID,
                                                                 TopFiveGeneFreq$ensemblID)]

#Order subjects based on library size after filtering of mtDNA genes
TopFiveProportion <- merge(TopFiveProportion, Metadata %>%
                             select(Filename, Batch, Cohort,
                                    DV200, RIN, MitoCount),
                           by = "Filename", sort = F)

TopFiveProportion$Filename <- factor(TopFiveProportion$Filename,
                                     levels = TopFiveSum$Filename)
TopFiveProportion$ensemblID2 <- factor(TopFiveProportion$ensemblID2,
                                       levels = unique(as.character(TopFiveGeneFreq$ensemblID2)))



TopFivePlot <- ggplot(TopFiveProportion %>% filter(Proportion > 0.006), aes(Filename, Proportion, fill = ensemblID2)) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, size = 8)) +
  labs(y = "Proportion of reads", x = "Sample", title = "Top five genes with the highest read count") + 
  scale_fill_manual(values = c(MoviePalettes$MoonRiseKingdomColors[3:9], gray.colors(6)),
                    name = "GeneID (n)") +
  geom_bar(stat = "identity")

ggarrange(TopFivePlot,
          ggplot(TopFiveProportion %>% group_by(Filename) %>%
                   summarise(TopFiveProp = sum(Proportion),
                             DV200 = mean(DV200),
                             mtDNAcount = mean(MitoCount)) %>% data.frame(),
                 aes(DV200, TopFiveProp, color = mtDNAcount)) +
            geom_point(),
          nrow = 2, heights = c(1.5, 1))
ggsave(paste0(ResultsPath, "TopFiveGenesProp.pdf"),
       device = "pdf", width = 10, height = 8, dpi = 300, useDingbats = F)


#Get the common genes with the highest count in majority of the samples and remove them from the count matrix
CommonTopGenes <- TopFiveGeneFreq[TopFiveGeneFreq$n > 0.5*ncol(MitoCountFiltered),] %>% filter(!duplicated(ensemblID))

CommonTopGenesSum <- apply(MitoCountFiltered[rownames(MitoCountFiltered )%in% CommonTopGenes$ensemblID,], 2, sum)

countMatrixFiltered <- MitoCountFiltered[!rownames(MitoCountFiltered) %in% as.character(CommonTopGenes$ensemblID),]

# Remove gene with less than 5 counts in  > 80% of the samples
ZeroCount <- apply(countMatrixFiltered, 1, function(x){
  sum(x < 5)
})

countMatrixFiltered <- countMatrixFiltered[ZeroCount < 0.8*ncol(countMatrixFiltered),]

#Create log2 CPM matrix after removal of mitochondria-encoded genes
cpmMatrixFiltered <- Count2CPM(countMatrixFiltered) %>% data.frame()
cpmMatrixFiltered <- apply(cpmMatrixFiltered, c(1,2), function(x) log2(x+1)) %>% data.frame()
cpmMatrixFiltered <- cbind(rownames(countMatrixFiltered), cpmMatrixFiltered)
colnames(cpmMatrixFiltered)[1] <- "genes"


#Add gene symbols
GeneSymbolAll <- data.frame(GeneSymbol = geneNames$gene_name[match(rownames(cpmMatrixFiltered), geneNames$gene_id2)],
                            Probe = rownames(cpmMatrixFiltered),
                            ensemblID = rownames(cpmMatrixFiltered))


ExpDataCPM <- cbind(GeneSymbolAll, cpmMatrixFiltered[-1])



studyFinal <- PreProccessRNAseq(Metadata = Metadata, expData = ExpDataCPM,
                                SexCol = "Sex", Combat = FALSE, resultsPath = ResultsPath)

SummaryTable <- studyFinal$Metadata %>% group_by(Cohort2) %>% summarise(RINmean = mean(RIN),
                                                                        RINmin = min(RIN),
                                                                        RINmax = max(RIN),
                                                                        DV200mean = mean(DV200),
                                                                        DV200min = min(DV200),
                                                                        DV200max = max(DV200)) %>% data.frame()
write.table(SummaryTable, paste0(ResultsPath, "GeneralRNAseQCsumary.tsv"), sep = "\t", row.names = F, col.names = T)

#Look at sample correlation 
SampleCor <- cor(studyFinal$ExpHigh %>% select(matches("SL")))
diag(SampleCor) <- NA

colnames(SampleCor) <- studyFinal$Metadata$Paper_ID[match(colnames(SampleCor), studyFinal$Metadata[[SampleIDcol]])]
rownames(SampleCor) <- studyFinal$Metadata$Paper_ID[match(rownames(SampleCor), studyFinal$Metadata[[SampleIDcol]])]


annoCol = data.frame(MT_GroupFinal = Metadata$kmeans_Final,
                     Age = Metadata$Age,
                     Sex = Metadata$Sex,
                     DV200 = Metadata$DV200,
                     Cohort = Metadata$Cohort2,
                     row.names = Metadata$Paper_ID)



pheatmap::pheatmap(SampleCor, angle_col = 90, na_col = "white", border_color = NA,
                   clustering_method = "ward.D2",cutree_cols = 3, 
                   color = colorRampPalette(c("darkblue", "gold2"))(999),
                   show_rownames = T, show_colnames = T,
                   annotation_col = annoCol,fontsize_col = 7, fontsize_row = 7,
                   annotation_colors = annoColors, filename = paste0(ResultsPath, "SampleCor.pdf"),
                   width = 12, height = 10)
closeDev()


#Plotting the PMI/Percent positive/DV200 group difference in each cohort
Colors <- c("#B3B3B3", "#BB9BEA", "#006633", "#7AD151")

PlotPMI <- ggplot(studyFinal$Metadata, aes(kmeans_Final, PMI)) +
  theme_bw() +
  labs(x = "", title = "PMI") +
  geom_boxplot(outlier.shape = NA, aes(fill = kmeans_Final), alpha = 0.6) +
  scale_fill_manual(values = Colors, name = "") +
  geom_jitter(width = 0.2) +
  facet_wrap(~Cohort2)

PlotPercentPositive <- ggplot(studyFinal$Metadata, aes(kmeans_Final, PercentPositive, fill = kmeans_Final)) +
  theme_bw() +
  labs(x = "", title = "IHC - %Positive") +
  geom_boxplot(outlier.shape = NA, aes(fill = kmeans_Final), alpha = 0.6) +
  scale_fill_manual(values = Colors, name = "") +
  geom_jitter(width = 0.2) +
  facet_wrap(~Cohort2)

PlotDV200 <- ggplot(studyFinal$Metadata, aes(kmeans_Final, DV200, fill = kmeans_Final)) +
  theme_bw() +
  labs(x = "", title = "DV200") +
  geom_boxplot(outlier.shape = NA, aes(fill = kmeans_Final), alpha = 0.6) +
  scale_fill_manual(values = Colors, name = "") +
  geom_jitter(width = 0.2) +
  facet_wrap(~Cohort2)

ggarrange(PlotPMI, PlotDV200, PlotPercentPositive, nrow = 3, common.legend = T, legend = "right")
ggsave(paste0(ResultsPath, "2KmeanCluster_PMI_DV200_IHC.pdf"), device = "pdf", width = 9, height = 8, useDingbats = F)



#Get cell type MGPS
source_url("https://github.com/ltoker/GeneralRscripts/blob/main/Cell_type_PCA.R?raw=T")

region = studyFinal$Metadata$NeuExpRegion %>% unique
CellType_genes <- GetMarkers(region)

#Exclude GabaPV genes which are not neuron specific in human (Darmanis) data
CellType_genes$GabaPV_Genes <- CellType_genes$GabaPV_Genes[!CellType_genes$GabaPV_Genes %in% c("WIF1", "TMEM132C", "BTN2A2")]

#Bootstrap with replacement the samples (90% of the samples)
SampleNames <- as.character(studyFinal$Metadata$Filename)

PCAresults <- sapply(paste0("Boot_", 1:100), function(boot){
  BootSamples <- sample(SampleNames, 0.9*length(SampleNames), replace = F)
  dataSub <- studyFinal$ExpHigh %>% select(c("GeneSymbol", BootSamples))
  PCA_genes_All_based(dataset_id=Study,
                      dataset=dataSub,
                      CellType_genes=CellType_genes,
                      contName = "SL",SampleReg = "SL",
                      NoiseThershold = studyFinal$NoiseThreshold)
}, simplify = F)

PCA_resultsMean <- sapply(names(PCAresults[[1]]$modified), function(celltype){
  temp <- data.frame(CommonName = names(PCAresults[[1]]$modified[[celltype]]$x[,1]),
                     Rot = PCAresults[[1]]$modified[[celltype]]$x[,1])
  for(i in 2:length(PCAresults)){
    temp <- merge(temp, data.frame(CommonName = names(PCAresults[[i]]$modified[[celltype]]$x[,1]),
                                   Rot = PCAresults[[i]]$modified[[celltype]]$x[,1]), by = "CommonName", all = TRUE, sort = F)
  }
  names(temp)[2:ncol(temp)] <- paste0("Rot", c(1:c(ncol(temp)-1)))
  temp$MeanRot <- rowMeans(temp[-1], na.rm = T)
  temp
}, simplify=FALSE)


#Add estimation to Metadata 
AllEstimates <- lapply(PCA_resultsMean, function(x){
  x$MeanRot
}) %>% do.call(cbind, .) %>% data.frame()
AllEstimates$CommonName <- PCA_resultsMean[[1]]$CommonName

#Remove the cell types for which MGPs we not calculated
AllEstimates <- AllEstimates[,apply(AllEstimates, 2, function(x) sum(is.na(x)) < nrow(AllEstimates))]

studyFinal$Metadata <- merge(studyFinal$Metadata, AllEstimates, by = "CommonName", sort = F)
studyFinal$Metadata <- studyFinal$Metadata[,!apply(studyFinal$Metadata, 2, function(x) sum(is.na(x))) == nrow(studyFinal$Metadata)]


studyFinal$countMatrix <- countMatrixFiltered[rownames(countMatrixFiltered) %in% studyFinal$ExpHigh$ensemblID,]

#match the sample order (yes, again.. check it later)
studyFinal$countMatrix <- studyFinal$countMatrix[,studyFinal$Metadata$Filename]

#Convert to integers
studyFinal$countMatrix <- apply(studyFinal$countMatrix, c(1,2), as.integer)
studyFinal$Metadata <- studyFinal$Metadata[,!apply(studyFinal$Metadata, 2, function(x) sum(is.na(x))) == nrow(studyFinal$Metadata)]


PCAall <- prcomp(t(studyFinal$ExpHigh %>% select(matches("SL"))), scale = T)
studyFinal$Metadata <- merge(studyFinal$Metadata, PCAall$x[,c(1:5)],
                             by.x = SampleIDcol, by.y = "row.names", sort = F)


#Get synaptic and OxPhos MGPs
PathwaysMitocartaDF <- read.table("GOannotations/Human.MitoCarta3.0.txt", header = T, sep = "\t")

pathwaysMitocarta <- sapply(PathwaysMitocartaDF$MitoPathway, function(Path){
  temp <- PathwaysMitocartaDF %>% filter(MitoPathway == Path) %>% .$Genes
  strsplit(temp, ",")[[1]] %>% sapply(function(x) gsub(" ", "", x)) %>% as.character()
}, simplify  = F)


PathwaysList <- list(GObp = gmtPathways("GOannotations/c5.go.bp.v7.2.symbols.xls"),
                     GOmf = gmtPathways("GOannotations/c5.go.mf.v7.2.symbols.xls"),
                     GOcc = gmtPathways("GOannotations/c5.go.cc.v7.2.symbols.xls"),
                     KEGG = gmtPathways("GOannotations/c2.cp.kegg.v7.2.symbols.xls"),
                     MitoCarta = pathwaysMitocarta)

##### This part is not directly included in the manuscript, but interesting to look at
PCApresynapse <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                       GOterm = PathwaysList$GOcc$GO_PRESYNAPSE)

PCApostsynapse <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                        GOterm = PathwaysList$GOcc$GO_POSTSYNAPSE)

PrePostIntersect <- intersect(PCApresynapse$GeneIn,
                              PCApostsynapse$GeneIn)

PreSynUnique <- PCApresynapse$GeneIn[!PCApresynapse$GeneIn %in% PrePostIntersect]

PCApresynapseUnique <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                             GOterm = PreSynUnique)

PostSynUnique <- PCApostsynapse$GeneIn[!PCApostsynapse$GeneIn %in% PrePostIntersect]

PCApostynapseUnique <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                             GOterm = PostSynUnique)

PCAoxphos <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                   GOterm = PathwaysList$KEGG$KEGG_OXIDATIVE_PHOSPHORYLATION)

PCAlysosome <- PCAgo(studyFinal$ExpHigh, sampleRegEx = "SL",
                     GOterm = PathwaysList$KEGG$KEGG_LYSOSOME)


studyFinal$Metadata$PostSynapse_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                        MGP = PCApostsynapse$PCA$x[,1])
studyFinal$Metadata$PreSynapse_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                       MGP = PCApresynapse$PCA$x[,1])
studyFinal$Metadata$PostSynapseUnique_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                              MGP = PCApostynapseUnique$PCA$x[,1])
studyFinal$Metadata$PreSynapseUnique_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                             MGP = PCApresynapseUnique$PCA$x[,1])

studyFinal$Metadata$PreSynapseUniqueES <- PCApresynapseUnique$EffectSize[match(studyFinal$Metadata[[SampleIDcol]],
                                                                               rownames(PCApresynapseUnique$EffectSize))]
studyFinal$Metadata$OxPhos_Genes <- AddGOmgpTometa(studyFinal$Metadata,
                                                   MGP = PCAoxphos$PCA$x[,1])
studyFinal$Metadata$OxPhosES <- PCAoxphos$EffectSize[match(studyFinal$Metadata[[SampleIDcol]],
                                                           rownames(PCAoxphos$EffectSize))]

studyFinal$Metadata$LysosomeES <- PCAlysosome$EffectSize[match(studyFinal$Metadata[[SampleIDcol]],
                                                               rownames(PCAlysosome$EffectSize))]

studyFinal$Metadata$PMI[is.na(studyFinal$Metadata$PMI)] <- median(studyFinal$Metadata$PMI,
                                                                  na.rm = T)

MitoGenes <- read.table("GOannotations/MRC_subunits_Haris_18.04.20.txt", header = T, sep = "\t") %>%
  filter(!is.na(Entrez.ID), !is.na(MRC.complex)) %>% droplevels()
AssemblyFactors <- MitoGenes %>%
  filter(Localization.function.in.complex == "assembly factor") %>%
  .$EnsemblGeneID %>% as.character()

MitoGenesList <- sapply(unique(MitoGenes$MRC.complex), function(x){
  genes <- MitoGenes %>% filter(MRC.complex == x, !EnsemblGeneID %in% (AssemblyFactors)) %>%
    .$EnsemblGeneID %>% as.character
  temp <- studyFinal$ExpHigh %>% filter(ensemblID %in% genes)
  temp2 <- temp %>% select(matches("SL"))
  rownames(temp2) <- as.character(temp$GeneSymbol)
  temp2
}, simplify = F)
names(MitoGenesList) <- sapply(names(MitoGenesList), function(x) gsub(" ", "_", x))

PCAcomplex <- lapply(MitoGenesList, function(Complex){
  temp <- prcomp(t(Complex), scale = T)
  if(sum(temp$rotation[,1]) < 0){
    temp$rotation[,1] <- -1*temp$rotation[,1]
    temp$x[,1] <- -1*temp$x[,1]
  }
  temp
})

PCAcomplexEffectSize <- sapply(names(MitoGenesList), function(Complex){
  Matrix <- MitoGenesList[[Complex]]
  t(Matrix) %*% rescale(PCAcomplex[[Complex]]$rotation[,1], c(0, 1/nrow(Matrix)))
}, simplify = F)

for(name in names(PCAcomplexEffectSize)){
  studyFinal$Metadata[[paste0(name, "_logRNA")]] <- PCAcomplexEffectSize[[name]][match(studyFinal$Metadata$RNAseq_id_ParkOme2,
                                                                                       rownames(PCAcomplexEffectSize[[name]])),1]
}

ProtCodGene <- geneNames %>% filter(gene_type == "protein_coding") %>%
  filter(!duplicated(gene_id)) %>% .$gene_id2 %>% as.character

lncRNA <- geneNames %>% filter(gene_type == "lncRNA") %>%
  filter(!duplicated(gene_id)) %>% .$gene_id2 %>% as.character

PseudoGene <- geneNames[geneNames$gene_type == "processed_pseudogene",] %>%
  filter(!duplicated(gene_id)) %>% .$gene_id2 %>% as.character

snoRNAGene <- geneNames[geneNames$gene_type == "snoRNA",] %>%
  filter(!duplicated(gene_id)) %>% .$gene_id2 %>% as.character

studyFinal$Metadata$ProtCodingSum <- apply(studyFinal$countMatrix[rownames(studyFinal$countMatrix) %in% ProtCodGene,], 2, sum)
studyFinal$Metadata$lncRNASum <- apply(studyFinal$countMatrix[rownames(studyFinal$countMatrix) %in% lncRNA,], 2, sum)
studyFinal$Metadata$PseudoGeneSum <- apply(studyFinal$countMatrix[rownames(studyFinal$countMatrix) %in% PseudoGene,], 2, sum)
studyFinal$Metadata$snoRNAGeneSum <- apply(studyFinal$countMatrix[rownames(studyFinal$countMatrix) %in% snoRNAGene,], 2, sum)

studyFinal$Metadata$FinalLibSize <- apply(studyFinal$countMatrix, 2, sum)

studyFinal$Metadata$PMI_binned <- cut(studyFinal$Metadata$PMI,
                                      breaks=c(seq(from=0, to=48, by=6),
                                               max(studyFinal$Metadata$PMI)),
                                      include.lowest=TRUE, labels=1:9) %>% as.integer


rm(countMatrix, ExpDataCPM, GeneSymbolAll, cpmMatrixFiltered,
   mitoGenes, countMito, countMatrixFiltered, MitoCountFiltered)


# Main analysis
Study = "MT_PD_3KnmeanClusters_combined"
ClusterCol = "kmeans_Final_combined"

annoColors$MT_Type = annoColors$MT_GroupFinal

source("ProjectScripts/RNAseq.R")


#Repeat the analysis in BBB (NTB) cohort separately , and look at the cohorts separately

#source("ProjectScripts/AdditionalAnalyses.R")

sessionInfo() %>% saveRDS("SessionInfo.Rds")