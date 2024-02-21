require("tidyverse")
require("MAST")
require("lme4")
options(readr.show_col_types = FALSE)

set.seed(1)

main_cell_types <- c('oligodendrocytes', 'ex', 'in', 'astrocytes')


# Read count data and subset CI-PD ("mitoPD") and nCI-PD ("PD")
if (!file.exists("./RData/sc_sce.Rds")) {
    ## Retrieve counts from Synapse.org repository
    ##  install.packages("synapser", repos = c("http://ran.synapse.org/", "http://cran.fhcrc.org/"))
    ##  https://r-docs.synapse.org/articles/manageSynapseCredentials.html
    require("synapser") 
    synLogin(authToken = "") # Create a token in Synapse.org
    syn53645160 <- synGet(entity = 'syn53645160', downloadLocation = "./RData/") 
}
sc_sce <- readRDS("./RData/sc_sce.Rds")

# Subset to PD cases only
sc_sce <- sc_sce[,colData(sc_sce)$Diagnosis %in% c("PD", "mitoPD")]
colData(sc_sce)$Diagnosis <- relevel(droplevels(colData(sc_sce)$Diagnosis), "PD")
gc()

# Metadata <- colData(obj) %>% as_tibble() %>% select(sample_id:Diagnosis) %>% distinct()

libsizes <- colSums(assays(sc_sce)$counts)/1e4

allClusters <- unique(sc_sce$cell_type2)

zlmCT_subclusters <- lapply(1:length(allClusters), function(clusterNum){
    lrtString <- "DiagnosismitoPD"
    ct <- allClusters[clusterNum]
    message(paste0(ct, " (", clusterNum, "/", length(allClusters), ")"))
    outputFile <- paste0("./Tables/DE_MASTRE_mitoPD_vs_PD_", lrtString, "_", ct, ".tsv")
    if (file.exists(outputFile)) {
        message("Output file found, SKIPPING")
        return(0)
    }
    minCells <- 0.25
    sce.sub <- sc_sce[,colData(sc_sce)$cell_type2 == ct]
    sce.sub <- sce.sub[rowSums(assays(sce.sub)$counts > 0) > (ncol(sce.sub)*minCells),]
    libsizes.sub <- libsizes[colData(sc_sce)$cell_type2 == ct]
    colData(sce.sub)$ngeneson <- scale(colSums(assay(sce.sub) > 0))
    assays(sce.sub)$et <- log2(sweep(assays(sce.sub)$counts, 2, libsizes.sub, `/`)+1)
    gc()
    sce.sub <- SceToSingleCellAssay(sce.sub, check_sanity = TRUE)
    gc()
    message(paste(c("GENES:", "BARCODES:"), dim(sce.sub), collapse = "\n"))

    zlm_model <- zlm(formula = as.formula("~ ngeneson + Sex + Age + PMI + Diagnosis + (1 | sample_id)"),
                     sca = sce.sub, method = 'glmer', strictConvergence = FALSE, ebayes = FALSE, silent = TRUE, parallel = TRUE)
    
    summaryDt <- MAST::summary(zlm_model, doLRT = FALSE, parallel = FALSE)$datatable
    lrtCond <- lrTest(zlm_model, Hypothesis(lrtString))
    fcHurdle <- summaryDt %>% as_tibble %>% filter(contrast == lrtString, component == 'logFC') %>%
        select(primerid, coef, ci.hi, ci.lo, z) %>%
        left_join(lrtCond[,,3] %>% as_tibble(rownames = 'primerid') %>% select(primerid, pval = hurdle), by = 'primerid') %>%
        na.omit
    Result <- fcHurdle %>% mutate(p_adj = p.adjust(pval)) %>%
        mutate(contrast = lrtString)
    message(paste0(" WRITING TO: ", outputFile))
    write_tsv(Result, file = outputFile)
    gc()
    return(1)
})

