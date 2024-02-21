require("tidyverse")
require("fgsea")


options(readr.show_col_types=FALSE)


set.seed(1)



# Pathways
gmt.files <- c(KEGG = "./Data/c2.cp.kegg.v7.4.symbols.gmt")
all(file.exists(gmt.files))


DE_files <- Sys.glob("./Tables/DE_MASTRE_mitoPD_vs_PD_*.tsv")
names(DE_files) <- gsub(pattern = "\\.tsv$", replacement = "", gsub(pattern = "^DE_MASTRE_mitoPD_vs_PD_", replacement = "", basename(DE_files)))
#all(file.exists(DE_files))
#DE_files

runORA <- function(res_obj, gmtFiles, mode = c("nondirectional", "up", "down"), pval_thr = 0.05, minSize = 5, maxSize = 200,
                    pval_colname = "pval", coef_colname = "coef") {
    mode <- match.arg(mode)
    message(paste0(" >> Running ORA ", mode, "-regulation"))
    if (mode == "nondirectional") {
        #topGenes <- filter(res_obj, pval < pval_thr)$gene_id %>% unique
        topGenes <- filter(res_obj, (!!as.symbol(pval_colname)) < pval_thr)$gene_id %>% unique
        Universe <- unique(res_obj$gene_id)
    } else if (mode == "up"){
        topGenes <- filter(res_obj, pval < pval_thr, coef > 0)$gene_id %>% unique
        topGenes <- filter(res_obj, (!!as.symbol(pval_colname)) < pval_thr, (!!as.symbol(coef_colname)) > 0)$gene_id %>% unique
        Universe <- unique(res_obj$gene_id)
    } else if (mode == "down") {
        topGenes <- filter(res_obj, pval < pval_thr, coef < 0)$gene_id %>% unique
        topGenes <- filter(res_obj, (!!as.symbol(pval_colname)) < pval_thr, (!!as.symbol(coef_colname)) < 0)$gene_id %>% unique
        Universe <- unique(res_obj$gene_id)
    } else {
        stop("NO")
    }
    message(paste0("    ", length(topGenes), "/", length(Universe), " (", round(100*length(topGenes)/(length(topGenes)+length(Universe)), digits = 2), "%) top genes"))
    if (length(topGenes)<10) {
        warning(" Less than 10 top genes... returning NA")
        return(tibble(pathway = NA, pval = NA, padj = NA, overlap = NA, size = NA, overlapGenes = NA, pathway_db = NA, type = mode))
    }
    # if (length(topGenes)/(length(topGenes)+length(Universe)) > 0.25) {
    #     warning(" Top genes are more than 25% of the Universe, consider lowering p-value")
    # }
    lapply(names(gmtFiles), function(Pw) {
        pathways <- gmtPathways(gmt.files[Pw])
        #message(paste0("  - ", Pw, " ORA..."))
        gseaRes <- fgsea::fora(pathways, genes = topGenes, universe = Universe, minSize = minSize, maxSize = maxSize) %>%
            as_tibble %>%
            mutate(pathway_db = Pw, type = mode)
    }) %>% do.call("bind_rows",.)
}



# Create folders for GSEA and ORA
for (db in names(gmt.files)) {
    dir.create(paste0("./Tables/ORA_mitoPD_vs_PD/nondirectional/", db), showWarnings = FALSE, recursive = TRUE)
}

ora_tests <- lapply(1:length(DE_files), function(nContrast){
    contrast <- names(DE_files)[nContrast]
    message(paste0(nContrast, "/", length(DE_files), " - ", contrast))
    dex <- read_tsv(DE_files[contrast], show_col_types = FALSE)
    colnames(dex)[1] <- "gene_id"
    ct <- str_match(contrast, "_(oligo.+|ex.*|in.*|astro.+|endo.+|OPC.*|micro.+)$")[,2]
    test_column <- 'pval'
    # nondir
    a_gsea <- runORA(res_obj = dex, gmtFiles = gmt.files, minSize = 5, maxSize = 1000, pval_thr = 0.05, pval_colname = test_column, mode = "nondirectional") %>%
        mutate(contrast = contrast, ct = ct, ORA_on = test_column) %>%
        rowwise() %>% mutate(overlapGenes = str_flatten_comma(overlapGenes)) %>% ungroup()
    for (db_name in unique(a_gsea$pathway_db)) {
        if (is.na(db_name)) next;
        a_gsea %>% filter(pathway_db == db_name) %>%
            mutate(padj = p.adjust(padj, method = "fdr")) %>% 
            arrange(padj, pval) %>%
            write_tsv(file = paste0("./Tables/ORA_mitoPD_vs_PD/nondirectional/", db_name, "/ORA_nondirectional_", contrast, "_", db_name, ".tsv"))
    }
    return(list(nondirectional = a_gsea))
})
names(ora_tests) <- names(DE_files)



# Combine ex, combine in, combine glia

# ex
appCon <- names(DE_files)[str_detect(names(DE_files), "DiagnosismitoPD_ex_")]
dex_list <- lapply(appCon, function(Contrast) {
    dex <- read_tsv(DE_files[Contrast], show_col_types = FALSE)
    ct <- str_match(Contrast, "_(oligo.+|ex.*|in.*|astro.+|endo.+|OPC.*|micro.+)$")[,2]
    dex <- dex %>% mutate(ct = ct)
    colnames(dex)[1] <- "gene_id"
    return(dex)
}) %>% Reduce("bind_rows", .)
L <- nrow(dex_list)

dex_table <- dex_list %>% select(gene_id, ct, pval) %>% pivot_wider(names_from = "ct", values_from = "pval")
dex_table$p.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
})
dex_table$signif.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harm <- harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
    ifelse(harm < 0.05*((1/L)*length(ps)), 1, 0)
})

Res <- dex_list %>% group_by(gene_id) %>%
    summarise(fdr.intersect = sum(p_adj < 0.05)) %>%
    left_join(dex_table, by = "gene_id") %>%
    select(gene_id, fdr.intersect, p.hmp, signif.hmp)

Res %>% mutate(fdr.intersect = ifelse(fdr.intersect > 0, 1, 0)) %>%
    group_by(fdr.intersect, signif.hmp) %>% tally()

ex_ora1 <- lapply(names(gmt.files), function(Pw) {
    pathways <- gmtPathways(gmt.files[Pw])
    oraRes <- fgsea::fora(pathways, genes = Res[Res$fdr.intersect > 0,]$gene_id,
                          universe = Res$gene_id, minSize = 5, maxSize = 1000) %>%
    as_tibble %>%
    mutate(pathway_db = Pw, type = "nondirectional_intersection")
}) %>% do.call("bind_rows",.) %>%
    mutate(contrast = "DiagnosismitoPDvsPD", ct = "all_ex", ora_on = "intersection_fdr") %>%
    rowwise() %>% mutate(overlapGenes = str_flatten_comma(overlapGenes)) %>% ungroup()

for (db_name in unique(ex_ora1$pathway_db)) {
    if (is.na(db_name)) next;
    ex_ora1 %>% filter(pathway_db == db_name) %>%
        mutate(padj = p.adjust(padj, method = "fdr")) %>% 
        arrange(padj, pval) %>%
        write_tsv(file = paste0("./Tables/ORA_mitoPD_vs_PD/nondirectional/", db_name, "/ORA_nondirectional_mitoPD_vs_PD_ex_fdr_intersection_", db_name, ".tsv"))
}




# in
appCon <- names(DE_files)[str_detect(names(DE_files), "DiagnosismitoPD_in_")]
dex_list <- lapply(appCon, function(Contrast) {
    dex <- read_tsv(DE_files[Contrast], show_col_types = FALSE)
    ct <- str_match(Contrast, "_(oligo.+|ex.*|in.*|astro.+|endo.+|OPC.*|micro.+)$")[,2]
    dex <- dex %>% mutate(ct = ct)
    colnames(dex)[1] <- "gene_id"
    return(dex)
}) %>% Reduce("bind_rows", .)
L <- nrow(dex_list)

dex_table <- dex_list %>% select(gene_id, ct, pval) %>% pivot_wider(names_from = "ct", values_from = "pval")
dex_table$p.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
})
dex_table$signif.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harm <- harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
    ifelse(harm < 0.05*((1/L)*length(ps)), 1, 0)
})

Res <- dex_list %>% group_by(gene_id) %>%
    summarise(fdr.intersect = sum(p_adj < 0.05)) %>%
    left_join(dex_table, by = "gene_id") %>%
    select(gene_id, fdr.intersect, p.hmp, signif.hmp)

Res %>% mutate(fdr.intersect = ifelse(fdr.intersect > 0, 1, 0)) %>%
    group_by(fdr.intersect, signif.hmp) %>% tally()

in_ora1 <- lapply(names(gmt.files), function(Pw) {
    pathways <- gmtPathways(gmt.files[Pw])
    oraRes <- fgsea::fora(pathways, genes = Res[Res$fdr.intersect > 0,]$gene_id,
                          universe = Res$gene_id, minSize = 5, maxSize = 1000) %>%
    as_tibble %>%
    mutate(pathway_db = Pw, type = "nondirectional_intersection")
}) %>% do.call("bind_rows",.) %>%
    mutate(contrast = "DiagnosismitoPDvsPD", ct = "all_in", ora_on = "intersection_fdr") %>%
    rowwise() %>% mutate(overlapGenes = str_flatten_comma(overlapGenes)) %>% ungroup()

for (db_name in unique(in_ora1$pathway_db)) {
    if (is.na(db_name)) next;
    in_ora1 %>% filter(pathway_db == db_name) %>%
        mutate(padj = p.adjust(padj, method = "fdr")) %>% 
        arrange(padj, pval) %>%
        write_tsv(file = paste0("./Tables/ORA_mitoPD_vs_PD/nondirectional/", db_name, "/ORA_nondirectional_mitoPD_vs_PD_in_fdr_intersection_", db_name, ".tsv"))
}





# glia
appCon <- names(DE_files)[str_detect(names(DE_files), "DiagnosismitoPD_(astro.+|endo.+|OPC.*|micro.+)_")]
dex_list <- lapply(appCon, function(Contrast) {
    dex <- read_tsv(DE_files[Contrast], show_col_types = FALSE)
    ct <- str_match(Contrast, "_(oligo.+|ex.*|in.*|astro.+|endo.+|OPC.*|micro.+)$")[,2]
    dex <- dex %>% mutate(ct = ct)
    colnames(dex)[1] <- "gene_id"
    return(dex)
}) %>% Reduce("bind_rows", .)
L <- nrow(dex_list)

dex_table <- dex_list %>% select(gene_id, ct, pval) %>% pivot_wider(names_from = "ct", values_from = "pval")
dex_table$p.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
})
dex_table$signif.hmp <- apply(dex_table, 1, function(Ps) {
    ps <- Ps[2:length(Ps)]
    ps <- ps[!is.na(ps)]
    harm <- harmonicmeanp::p.hmp(ps, rep(1/L, length(ps)), L = L)
    ifelse(harm < 0.05*((1/L)*length(ps)), 1, 0)
})

Res <- dex_list %>% group_by(gene_id) %>%
    summarise(fdr.intersect = sum(p_adj < 0.05)) %>%
    left_join(dex_table, by = "gene_id") %>%
    select(gene_id, fdr.intersect, p.hmp, signif.hmp)

Res %>% mutate(fdr.intersect = ifelse(fdr.intersect > 0, 1, 0)) %>%
    group_by(fdr.intersect, signif.hmp) %>% tally()

glia_ora1 <- lapply(names(gmt.files), function(Pw) {
    pathways <- gmtPathways(gmt.files[Pw])
    oraRes <- fgsea::fora(pathways, genes = Res[Res$fdr.intersect > 0,]$gene_id,
                          universe = Res$gene_id, minSize = 5, maxSize = 1000) %>%
    as_tibble %>%
    mutate(pathway_db = Pw, type = "nondirectional_intersection")
}) %>% do.call("bind_rows",.) %>%
    mutate(contrast = "DiagnosismitoPDvsPD", ct = "all_glia", ora_on = "intersection_fdr") %>%
    rowwise() %>% mutate(overlapGenes = str_flatten_comma(overlapGenes)) %>% ungroup()

for (db_name in unique(glia_ora1$pathway_db)) {
    if (is.na(db_name)) next;
    glia_ora1 %>% filter(pathway_db == db_name) %>%
        mutate(padj = p.adjust(padj, method = "fdr")) %>% 
        arrange(padj, pval) %>%
        write_tsv(file = paste0("./Tables/ORA_mitoPD_vs_PD/nondirectional/", db_name, "/ORA_nondirectional_mitoPD_vs_PD_glia_fdr_intersection_", db_name, ".tsv"))
}

