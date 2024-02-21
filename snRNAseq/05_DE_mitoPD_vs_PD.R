require("tidyverse")
require('openxlsx')

options(readr.show_col_types = FALSE)
source("./functions.R")



geneNames.filename <- "./RData/geneNames.Rds"
geneNames <- readRDS(geneNames.filename)
prot_coding_genes <- geneNames %>% filter(!str_starts(gene_name, "MT-")) %>%
    filter(gene_type == "protein_coding") %>%
    pull(gene_name)


ordered.diagnosis <- c("mitoPD", "PD")
ordered.diseases <- c("mitoPD")



DE_files <- Sys.glob("./Tables/DE_MASTRE_mitoPD_vs_PD_DiagnosismitoPD*.tsv")
names(DE_files) <- gsub(pattern = "\\.tsv$", replacement = "", gsub(pattern = "^DE_MASTRE_mitoPD_vs_PD_", replacement = "", basename(DE_files)))
#all(file.exists(DE_files))
#DE_files


# Read all DE files

dex <- lapply(names(DE_files), function(defile){
    message(defile)
    DE <- read_tsv(DE_files[defile])
    regex <- paste0("^DiagnosismitoPD_(.+)$")
    cluster <- str_extract(defile, regex, group = 1)
    DE %>% mutate(cluster = cluster, diagnosis = "mitoPD") %>% select(-contrast) 
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
dex %>% filter(p_adj < 0.05) %>% distinct(gene_id, diagnosis) %>%
    group_by(gene_id) %>%
    tally() %>%
    filter(n == 2)





# Expression heatmaps and DEG count barplots
# For heatmaps, plotting only genes that exhibit an absolute log2 FC above 0.5,
# regardless of significance level. For barplots, count genes below 5% FDR.


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



# All prot-coding genes
Plot1 <- plotHeatmap(dexData = dex, Diagnosis = "mitoPD", min_coef = 0.5, min_p = 0.05, Cap = 2)
# > png("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD.png", height = 4500, width = 1700, res = 200)
# > svg("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD.svg", height = 21, width = 8)
# > dev.off()
draw(Plot1$plot, annotation_legend_list = list(Plot1$legend))
#plotBarplot(dexData = dex, Diagnosis = "mitoPD", min_coef = 0.322, min_p = 0.05, which.p = "fdr")



# [S11] OXPHOS nuc
Plot4 <- plotHeatmap(dexData = dex, Diagnosis = "mitoPD", min_coef = 0, min_p = 1,
                     Genes = mitocarta %>% filter(!is.na(oxphos_subunits), chr != "chrM") %>% arrange(oxphos_subunits) %>% pull(gene_name))
# > png("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD_OXPHOS_nuc.png", height = 1660, width = 1600, res = 200)
# > draw(Plot4$plot, annotation_legend_list = list(Plot4$legend))
# > dev.off()
# > tiff("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD_OXPHOS_nuc.tiff", height = 1660, width = 1600, res = 200)
# > draw(Plot4$plot, annotation_legend_list = list(Plot4$legend))
# > dev.off()
# > svg("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD_OXPHOS_nuc.svg", height = 9.5, width = 8)
# > draw(Plot4$plot, annotation_legend_list = list(Plot4$legend))
# > dev.off()
# > pdf("./Figs_mitoPD/DE_mitoPD_vs_non-mitoPD_OXPHOS_nuc.pdf", height = 9.5, width = 8)
# > draw(Plot4$plot, annotation_legend_list = list(Plot4$legend))
# > dev.off()
# > dev.off()
draw(Plot4$plot, annotation_legend_list = list(Plot4$legend))


# Number of dysregulated genes per cell type
dex %>% filter(gene_id %in% genes.in.autosomes) %>%
    mutate(significant = ifelse(p_adj < 0.05, "SIGN", "NO")) %>%
    group_by(cluster, significant) %>%
    tally() %>%
    filter(significant == "SIGN") %>%
    print(n = Inf)

