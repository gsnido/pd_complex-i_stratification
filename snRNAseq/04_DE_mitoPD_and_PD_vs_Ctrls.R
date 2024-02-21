require("tidyverse")


options(readr.show_col_types = FALSE)
source("./functions.R")



geneNames.filename <- "./RData/geneNames.Rds"
geneNames <- readRDS(geneNames.filename)
prot_coding_genes <- geneNames %>% filter(!str_starts(gene_name, "MT-")) %>%
    filter(gene_type == "protein_coding") %>%
    pull(gene_name)


ordered.diagnosis <- c("mitoPD", "PD")
ordered.diseases <- c("mitoPD")



DE_files <- Sys.glob("./Tables/DE_MASTRE_Diagnosis*.tsv")
names(DE_files) <- gsub(pattern = "\\.tsv$", replacement = "", gsub(pattern = "^DE_", replacement = "", basename(DE_files)))
#all(file.exists(DE_files))
#DE_files


# Read all DE files
dex <- lapply(names(DE_files), function(defile){
    DE <- read_tsv(DE_files[defile])
    regex <- paste0("^MASTRE_",DE$contrast[1], "_(.+)$")
    Diagnosis <- str_replace(DE$contrast[1], "Diagnosis", "")
    Cluster <- str_extract(defile, regex, group = 1)
    DE %>% mutate(cluster = Cluster, diagnosis = Diagnosis) %>% select(-contrast)
}) %>% Reduce("bind_rows",.) %>%
    filter(primerid %in% filter(geneNames, gene_type == "protein_coding")$gene_name)
colnames(dex) <- c("gene_id", "coef", "ci.hi", "ci.lo", "z", "pval", "p_adj", "cluster", "diagnosis")

autosomes <- paste0("chr", 1:22)
genes.in.autosomes <- geneNames %>% filter(seqnames %in% autosomes) %>% pull(gene_name) %>% unique()

# nb of significant genes in at least one cell types
dex %>% filter(gene_id %in% genes.in.autosomes) %>%
    filter(p_adj < 0.05) %>% distinct(gene_id, diagnosis) %>%
    group_by(diagnosis) %>%
    tally()

# nb of DE genes common to the two groups
dex %>% filter(gene_id %in% genes.in.autosomes) %>%
    filter(p_adj < 0.05) %>% distinct(gene_id, diagnosis) %>%
    group_by(gene_id) %>%
    tally() %>%
    filter(n == 2) %>%
    nrow()

# DEGs in each cell type cluster
ct.order <- c("ex_01", "ex_02", "ex_06", "ex_09", "ex_11", "ex_17", "ex_19", "ex_21", "ex_22", "ex_23", "ex_24", "ex_30",
              "in_08", "in_10", "in_12", "in_15", "in_20", "in_26", "in_27", "in_28",
              "oligodendrocytes_00", "oligodendrocytes_03", "oligodendrocytes_04",
              "astrocytes_05", "astrocytes_18", "microglia_13", "endothelial_25", "OPC_07")

common_degs <- dex %>% filter(p_adj < 0.05) %>%
    filter(gene_id %in% genes.in.autosomes) %>%
    group_by(gene_id, cluster) %>%
    tally() %>% ungroup() %>%
    filter(n == 2) %>%
    group_by(cluster) %>%
    tally() %>% ungroup() %>%
    mutate(cluster = factor(cluster, levels = ct.order)) %>%
    arrange(cluster) %>%
    select(cluster, Common = n)

per_diagnosis_degs <- dex %>% filter(p_adj < 0.05) %>%
    filter(gene_id %in% genes.in.autosomes) %>%
    group_by(diagnosis, cluster) %>%
    tally() %>% ungroup() %>%
    pivot_wider(names_from = "diagnosis", values_from = "n") %>%
    mutate(cluster = factor(cluster, ct.order)) %>%
    arrange(cluster)

deg_count <- full_join(per_diagnosis_degs, common_degs, by = join_by(cluster)) %>%
    replace_na(replace = list(PD = 0, mitoPD = 0, Common = 0)) %>%
    complete(cluster, fill = list(PD = 0, mitoPD = 0, Common = 0)) %>%
    select(cluster, "CI-PD" = mitoPD, "nCI-PD" = PD, Common)

# > pdf("./Figs/F7c_barplot.pdf", height = 5.5, width = 8.5)
# > png("./Figs/F7c_barplot.png", height = 1100, width = 1700, res = 200)
# > svg("./Figs/F7c_barplot.svg", height = 5.5, width = 8.5)
# > tiff("./Figs/F7c_barplot.tiff", height = 1100, width = 1700, res = 200)
deg_count %>% 
    mutate(cluster = fct_rev(cluster)) %>%
    pivot_longer(-cluster, names_to = "Names", values_to = "counts") %>%
    mutate(Names = factor(Names, levels = c("CI-PD", "nCI-PD", "Common"))) %>%
    ggplot() +
    geom_col(aes(counts, cluster, fill = Names), width = 0.9) +
    scale_fill_manual(values = c('darkgreen', 'purple', 'grey50')) +
    facet_grid(~Names) +
    theme_bw() +
    theme(axis.title.y = element_blank(), legend.position = "none") +
    labs(x = "number of differentially expressed genes (P adj. < 0.05)")
# > dev.off()



mitocarta <- openxlsx::read.xlsx("./Data/Human.MitoCarta3.0.xlsx", sheet = 2) %>% as_tibble() %>%
    select(gene_name = Symbol, pathway = MitoCarta3.0_MitoPathways, chr = hg19_Chromosome) %>%
    filter(!is.na(chr), pathway != "0") %>%
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X", "M")))) %>%
    mutate(oxphos_subunits = str_extract(pathway, 'OXPHOS > (Complex [IV]+) > (.+) subunits', group = 1)) %>%
    mutate(oxphos_subunits = factor(oxphos_subunits, levels = c("Complex I", "Complex II", "Complex III", "Complex IV", "Complex V", "Cytochrome C"))) %>%
    mutate(metabolism = str_trim(str_extract(pathway, 'Metabolism > ([^>\\|]+)', group = 1))) %>%
    mutate(central_dogma = str_trim(str_extract(pathway, 'Mitochondrial central dogma > ([^>\\|]+)', group = 1))) %>%
    mutate(oxphos_factors = str_trim(str_extract(pathway, 'OXPHOS > OXPHOS (assembly factors)', group = 1))) %>%
    mutate(prot_things = str_trim(str_extract(pathway, 'Protein import, sorting and homeostasis > ([^>\\|]+)', group = 1))) %>%
    mutate(signaling = str_trim(str_extract(pathway, 'Signaling > ([^>\\|]+)', group = 1))) %>%
    mutate(mt_dynamics = str_trim(str_extract(pathway, 'Mitochondrial dynamics and surveillance > ([^>\\|]+)', group = 1))) %>%
    mutate(small_mol_trans = str_trim(str_extract(pathway, 'Small molecule transport > ([^>\\|]+)', group = 1))) 
mitocarta[mitocarta$pathway == "Signaling",]$signaling <- "Signaling"
mitocarta[mitocarta$pathway == "Small molecule transport",]$small_mol_trans <- "Small molecule transport"
mitocarta[mitocarta$pathway == "Mitochondrial dynamics and surveillance",]$mt_dynamics <- "Mitochondrial dynamics and surveillance"
mitocarta[mitocarta$pathway == "Small molecule transport | Signaling",]$small_mol_trans <- "Small molecule transport"
mitocarta[mitocarta$pathway == "Small molecule transport | Signaling",]$signaling <- "Signaling"


# HEATMAPS CI-PD vs Ctrl / nCI-PD vs Ctrl

# [S10] All nuclear prot-coding genes
Plot1a <- dex %>% filter(!str_starts(gene_id, "MT-")) %>%
    plotHeatmap(Diagnosis = "mitoPD", min_coef = log2(1.25), min_p = 0.05, Cap = 2, title = "CI-PD")
# > png("./Figs/S10_A-mitoPD_vs_Ctrl.png", height = 1950, width = 1600, res = 200)
# > draw(Plot1a$plot, annotation_legend_list = list(Plot1a$legend))
# > dev.off()
# > tiff("./Figs/S10_A-mitoPD_vs_Ctrl.tiff", height = 1950, width = 1600, res = 200)
# > draw(Plot1a$plot, annotation_legend_list = list(Plot1a$legend))
# > dev.off()
# > svg("./Figs/S10_A-mitoPD_vs_Ctrl.svg", height = 9.5, width = 8)
# > draw(Plot1a$plot, annotation_legend_list = list(Plot1a$legend))
# > dev.off()
# > pdf("./Figs/S10_A-mitoPD_vs_Ctrl.pdf", height = 9.5, width = 8)
# > draw(Plot1a$plot, annotation_legend_list = list(Plot1a$legend))
# > dev.off()
draw(Plot1a$plot, annotation_legend_list = list(Plot1a$legend))

## dex %>% filter(!str_starts(gene_id, "MT-")) %>%
##     plotBarplot(Diagnosis = "mitoPD", min_coef = 0.322, min_p = 0.05, which.p = "fdr")

Plot1b <- dex %>% filter(!str_starts(gene_id, "MT-")) %>%
    plotHeatmap(Diagnosis = "PD", min_coef = log2(1.25), min_p = 0.05, Cap = 2, title = "nCI-PD")
# > png("./Figs/S10_A-non-mitoPD_vs_Ctrl.png", height = 3900, width = 1600, res = 200)
# > draw(Plot1b$plot, annotation_legend_list = list(Plot1b$legend))
# > dev.off()
# > tiff("./Figs/S10_A-non-mitoPD_vs_Ctrl.tiff", height = 3900, width = 1600, res = 200)
# > draw(Plot1b$plot, annotation_legend_list = list(Plot1b$legend))
# > dev.off()
# > svg("./Figs/S10_A-non-mitoPD_vs_Ctrl.svg", height = 19.5, width = 8)
# > draw(Plot1b$plot, annotation_legend_list = list(Plot1b$legend))
# > dev.off()
# > pdf("./Figs/S10_A-non-mitoPD_vs_Ctrl.pdf", height = 19.5, width = 8)
# > draw(Plot1b$plot, annotation_legend_list = list(Plot1b$legend))
# > dev.off()
draw(Plot1b$plot, annotation_legend_list = list(Plot1b$legend))





# OXPHOS mt

# [S12, S13] OXPHOS nuc
Plot4a <- plotHeatmap(dexData = dex, Diagnosis = "mitoPD", min_coef = 0, min_p = 1,
                     Genes = mitocarta %>% filter(!is.na(oxphos_subunits), chr != "chrM") %>% arrange(oxphos_subunits) %>% pull(gene_name))
# > png("./Figs_mitoPD/DE_mitoPD_vs_Control_OXPHOS_nuc.png", height = 1660, width = 1600, res = 200)
# > draw(Plot4a$plot, annotation_legend_list = list(Plot4a$legend))
# > dev.off()
# > tiff("./Figs_mitoPD/DE_mitoPD_vs_Control_OXPHOS_nuc.tiff", height = 1660, width = 1600, res = 200)
# > draw(Plot4a$plot, annotation_legend_list = list(Plot4a$legend))
# > dev.off()
# > svg("./Figs_mitoPD/DE_mitoPD_vs_Control_OXPHOS_nuc.svg", height = 8, width = 8)
# > draw(Plot4a$plot, annotation_legend_list = list(Plot4a$legend))
# > dev.off()
# > pdf("./Figs_mitoPD/DE_mitoPD_vs_Control_OXPHOS_nuc.pdf", height = 8, width = 8)
# > draw(Plot4a$plot, annotation_legend_list = list(Plot4a$legend))
# > dev.off()
draw(Plot4a$plot, annotation_legend_list = list(Plot4a$legend))

Plot4b <- plotHeatmap(dexData = dex, Diagnosis = "PD", min_coef = 0, min_p = 1,
                     Genes = mitocarta %>% filter(!is.na(oxphos_subunits), chr != "chrM") %>% arrange(oxphos_subunits) %>% pull(gene_name))
# > png("./Figs_mitoPD/DE_PD_vs_Control_OXPHOS_nuc.png", height = 1660, width = 1600, res = 200)
# > draw(Plot4b$plot, annotation_legend_list = list(Plot4b$legend))
# > dev.off()
# > tiff("./Figs_mitoPD/DE_PD_vs_Control_OXPHOS_nuc.tiff", height = 1660, width = 1600, res = 200)
# > draw(Plot4b$plot, annotation_legend_list = list(Plot4b$legend))
# > dev.off()
# > svg("./Figs_mitoPD/DE_PD_vs_Control_OXPHOS_nuc.svg", height = 8, width = 8)
# > draw(Plot4b$plot, annotation_legend_list = list(Plot4b$legend))
# > dev.off()
# > pdf("./Figs_mitoPD/DE_PD_vs_Control_OXPHOS_nuc.pdf", height = 8, width = 8)
# > draw(Plot4b$plot, annotation_legend_list = list(Plot4b$legend))
# > dev.off()
draw(Plot4b$plot, annotation_legend_list = list(Plot4b$legend))

