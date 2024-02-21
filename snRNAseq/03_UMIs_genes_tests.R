require('tidyverse')
require("ggpmisc")
require("SingleCellExperiment")
require("betareg")
require("lme4")
require("lmerTest")

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
combine <- dplyr::combine

options(readr.show_col_types = FALSE)

source("./functions.R")


# Read data
if (!file.exists("./RData/sc_sce.Rds")) {
    ## Retrieve counts from Synapse.org repository
    ##  install.packages("synapser", repos = c("http://ran.synapse.org/", "http://cran.fhcrc.org/"))
    ##  https://r-docs.synapse.org/articles/manageSynapseCredentials.html
    require("synapser") 
    synLogin(authToken = "") # Create a token in Synapse.org
    syn53645160 <- synGet(entity = 'syn53645160', downloadLocation = "./RData/") 
}
sc_sce <- readRDS("./RData/sc_sce.Rds")


colData(sc_sce)$Diagnosis <- relevel(colData(sc_sce)$Diagnosis, "Control")
assays(sc_sce)$logcounts <- NULL
gc()

Metadata <- colData(sc_sce) %>% as_tibble() %>% distinct(sample_id, Biobank_ID, Cohort, Age, Sex, PMI, Diagnosis)
Metadata %>% print(n = Inf)



ordered.diseases <- c("mitoPD", "PD")
ordered.diagnosis <- c("Control", ordered.diseases)
ordered.major.cts <- c("ex", "in", "oligodendrocytes", "astrocytes", "microglia", "OPC", "endothelial")
ordered.individuals <- Metadata %>% arrange(sample_id) %>% pull(sample_id)
ct.order2 <- c("ex_01", "ex_02", "ex_06", "ex_09", "ex_11", "ex_17", "ex_19", "ex_21", "ex_22", "ex_23", "ex_24", "ex_30",
               "in_08", "in_10", "in_12", "in_15", "in_20", "in_26", "in_27", "in_28",
               "oligodendrocytes_00", "oligodendrocytes_03", "oligodendrocytes_04",
               "astrocytes_05", "astrocytes_18", "microglia_13", "endothelial_25", "OPC_07")


# Test for differences in nb of nuclei captured
nbNuc <- colData(sc_sce) %>% as_tibble() %>%
    group_by(sample_id, Diagnosis) %>%
    tally()
aov(n ~ Diagnosis, data = nbNuc) %>% summary()






# Test differences in the cell type-specific numbers of nuclei

txData <- colData(sc_sce) %>% as_tibble(rownames = "cell") %>%
    select(cell, sample_id, Biobank_ID, Diagnosis, nCount_RNA, nFeature_RNA, Sex, Age, PMI, cell_type, cell_type2) %>%
    mutate(sample_id = as.factor(sample_id)) %>%
    mutate(cell_type2 = factor(cell_type2, levels = ct.order2)) %>%
    mutate(cell_type = factor(cell_type, levels = ordered.major.cts)) %>%
    group_by(sample_id) %>%
    mutate(nCount_norm = log2(nCount_RNA / median(nCount_RNA))) %>%
    mutate(nFeature_norm = log2(nFeature_RNA)) %>%
    ungroup()



# Test association of disease status with transcriptional output in each cell
# type cluster

tx_out_all <- lapply(ct.order2, function(CT) {
    message(CT)
    tmpData <- txData %>% mutate(Age = scale2(Age), PMI = scale2(PMI), Sex = factor(Sex, levels = c("M", "F"))) %>% filter(cell_type2 == CT)
    fit1 <- lmer(formula = nFeature_norm ~ PMI + Sex + Age + Diagnosis + (1|sample_id),
                data = tmpData)
    confint1 <- confint(fit1)
    colnames(confint1) <- c("ci.lo", "ci.hi")
    res1 <- summary(fit1)$coefficients %>% as_tibble(rownames = "term") %>%
        left_join(confint1 %>% as_tibble(rownames = "term") %>% filter(!str_starts(term, "\\.")), by = "term") %>%
        mutate(cell_type = CT, dep_var = "nFeature")
    fit2 <- lmer(formula = nCount_norm ~ PMI + Sex + Age + Diagnosis + (1|sample_id),
                data = tmpData)
    confint2 <- confint(fit2)
    colnames(confint2) <- c("ci.lo", "ci.hi")
    res2 <- summary(fit2)$coefficients %>% as_tibble(rownames = "term") %>%
        left_join(confint2 %>% as_tibble(rownames = "term") %>% filter(!str_starts(term, "\\.")), by = "term") %>%
        mutate(cell_type = CT, dep_var = "nCount")
    bind_rows(res1, res2)
}) %>% Reduce("bind_rows", .)

colnames(tx_out_all) <- c("term", "coef", "stderr", "df", "stat", "pval", "ci.lo", "ci.hi", "cell_type", "dep_var")


# > svg("./Figs/tx_output_assoc_disease.svg", height = 7, width = 9)
# > png("./Figs/tx_output_assoc_disease.png", height = 1200, width = 1300, res = 150)
tx_out_all %>% filter(term %in% c("DiagnosisPD", "DiagnosismitoPD")) %>%
    filter(dep_var == "nFeature") %>%
    mutate(cell_type = factor(cell_type, levels = rev(ct.order2))) %>%
    mutate(term = factor(as.character(term), levels = c("DiagnosisPD", "DiagnosismitoPD"))) %>%
    mutate(padj = p.adjust(pval, method = 'fdr')) %>%
    mutate(stars = stars.pvalue(padj)) %>%
    mutate(stars = ifelse(stars == '.', '', stars)) %>%
    mutate(stars.left = ifelse(coef < 0, stars, '')) %>%
    mutate(stars.right = ifelse(coef > 0, stars, '')) %>%
    mutate(pos.star = ifelse(coef < 0, ci.lo - 0.05, ci.hi + 0.05)) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    geom_pointrange(aes(x = cell_type, y = coef, ymin = ci.lo, ymax = ci.hi), size = .2, linewidth = .4) +
    geom_text(aes(x = cell_type, y = pos.star, label = stars.left), size = 5, hjust = 1, vjust = 0.75) +
    geom_text(aes(x = cell_type, y = pos.star, label = stars.right), size = 5, hjust = 0, vjust = 0.75) +
    facet_grid(~term) +
    #scale_y_continuous(name = "coefficient", limits = symmetric_limits) +
    scale_y_continuous(name = "coefficient", limits = c(-1.6, 1.6)) +
    coord_flip() +
    xlab("")
dev.off()
# > dev.off()

tx_out_all %>% filter(term %in% c("DiagnosisPD", "DiagnosismitoPD")) %>%
    filter(dep_var == "nCount") %>%
    mutate(cell_type = factor(cell_type, levels = rev(ct.order2))) %>%
    mutate(term = factor(as.character(term), levels = c("DiagnosisPD", "DiagnosismitoPD"))) %>%
    mutate(padj = p.adjust(pval, method = 'fdr')) %>%
    #mutate(stars = stars.pvalue(pval)) %>%
    mutate(stars = stars.pvalue(padj)) %>%
    mutate(stars = ifelse(stars == '.', '', stars)) %>%
    mutate(stars.left = ifelse(coef < 0, stars, '')) %>%
    mutate(stars.right = ifelse(coef > 0, stars, '')) %>%
    mutate(pos.star = ifelse(coef < 0, ci.lo - 0.05, ci.hi + 0.05)) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    geom_pointrange(aes(x = cell_type, y = coef, ymin = ci.lo, ymax = ci.hi), size = .2, linewidth = .4) +
    geom_text(aes(x = cell_type, y = pos.star, label = stars.left), size = 5, hjust = 1, vjust = 0.75) +
    geom_text(aes(x = cell_type, y = pos.star, label = stars.right), size = 5, hjust = 0, vjust = 0.75) +
    facet_grid(~term) +
    #scale_y_continuous(name = "coefficient", limits = symmetric_limits) +
    scale_y_continuous(name = "coefficient", limits = c(-2, 2)) +
    coord_flip() +
    xlab("")
dev.off()

