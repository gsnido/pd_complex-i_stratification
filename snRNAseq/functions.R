require("ComplexHeatmap")
require("RColorBrewer")
require("circlize")

scale2 <- function(x, na.rm=FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

stars.pvalue <- function(pval) {
    stars <- rep('', length(pval))
    stars[pval<=0.1] <- "."
    stars[pval<=0.05] <- "*"
    stars[pval<=0.01] <- "**"
    stars[pval<=0.001] <- "***"
    return(stars)
}

# Colorblind-friendly colours
cb_red <- "#DC3220"
cb_blue <- "#005AB5"

make_matrix <- function(df){
    my_matrix <- as.matrix(df[-1])
    rownames(my_matrix) <- c(df[[1]])
    my_matrix
}


# Function to plot a barplot with the number of significant genes

plotHeatmap <- function(dexData, Diagnosis, min_coef, min_p, Cap = NA, Genes = NULL, Reorder = TRUE, genenames.side = "right", title = NA) {
    row.annot.tmp <- dexData %>% filter(diagnosis == Diagnosis, p_adj < min_p, abs(coef) > min_coef, cluster %in% rownames(ct.order)) %>%
        select(gene_id, cluster, p_adj) %>% 
        group_by(gene_id) %>% 
        arrange(p_adj) %>% 
        slice_head(n = 1) %>% ungroup()
    if (!is.null(Genes)) {
        which.sig.tmp <- Genes
    } else {
        which.sig.tmp <- pull(row.annot.tmp, gene_id)
    }


    row.annot.df.tmp <- as.data.frame(select(row.annot.tmp, cluster) %>% left_join(as_tibble(ct.order, rownames = "cluster"), by = "cluster"))
    #all(row.annot.tmp$cluster == row.annot.df.tmp$cluster)
    rownames(row.annot.df.tmp) <- row.annot.tmp$gene_id

    cb_red <- "#DC3220"
    cb_blue <- "#005AB5"

    coefs <- lapply(c(rownames(ct.order)), function(Cluster) {
        res.tmp <- dexData %>% filter(gene_id %in% which.sig.tmp, cluster == Cluster, diagnosis == Diagnosis) %>%
            select(gene_id, coef)
        colnames(res.tmp) <- c("gene_id", Cluster)
        return(res.tmp)
    }) %>% Reduce(function(x, y) full_join(x, y, by = "gene_id"), .) %>% make_matrix
    
    zs <- lapply(c(rownames(ct.order)), function(Cluster) {
        res.tmp <- dexData %>% filter(gene_id %in% which.sig.tmp, cluster == Cluster, diagnosis == Diagnosis) %>%
            select(gene_id, z)
        colnames(res.tmp) <- c("gene_id", Cluster)
        return(res.tmp)
    }) %>% Reduce(function(x, y) full_join(x, y, by = "gene_id"), .) %>% make_matrix
    
    pvals <- lapply(c(rownames(ct.order)), function(Cluster) {
        res.tmp <- dexData %>% filter(gene_id %in% which.sig.tmp, cluster == Cluster, diagnosis == Diagnosis) %>%
            select(gene_id, pval)
        colnames(res.tmp) <- c("gene_id", Cluster)
        return(res.tmp)
    }) %>% Reduce(function(x, y) full_join(x, y, by = "gene_id"), .) %>% make_matrix
    
    padjs <- lapply(c(rownames(ct.order)), function(Cluster) {
        res.tmp <- dexData %>% filter(gene_id %in% which.sig.tmp, cluster == Cluster, diagnosis == Diagnosis) %>%
            select(gene_id, p_adj)
        colnames(res.tmp) <- c("gene_id", Cluster)
        return(res.tmp)
    }) %>% Reduce(function(x, y) full_join(x, y, by = "gene_id"), .) %>% make_matrix
    
    if (Reorder) {
        if (is.null(Genes)) {
            coefs2 <- coefs[match(row.annot.tmp$gene_id, rownames(coefs)),]
            row.order <- lapply(unique(row.annot.df.tmp$ct), function(CT) {
                xsub <- coefs2[row.annot.df.tmp$ct == CT,, drop = FALSE]
                xsub <- xsub[,colnames(xsub) %in% row.annot.df.tmp[row.annot.df.tmp$ct == CT,]$cluster, drop = FALSE]
                xsub <- sign(replace(xsub, is.na(xsub), 0))
                if (nrow(xsub) == 0) {
                    message(paste0("WARNING: No genes found in ", CT))
                    return(NULL)
                } else if (nrow(xsub) <= 2) {
                    message(paste0("WARNING: Few genes (1 or 2) found in ", CT, " (not clustering)"))
                    return(rownames(xsub))
                } else {
                    hc <- hclust(dist(replace(xsub, is.na(xsub), 0)))
                    sv <- svd(t(replace(xsub, is.na(xsub), 0)))$v[,1]
                    hc2 <- as.hclust(reorder(as.dendrogram(hc), wts = sv))
                    return(rownames(xsub)[hc2$order])
                }
            })
            names(row.order) <- unique(row.annot.df.tmp$ct)
            row.order <- do.call("c", row.order[major.ct])
        } else {
            xsub <- sign(replace(coefs, is.na(coefs), 0))
            if(nrow(xsub) == 0) {
                message(paste0("WARNING: No genes found, returning NULL"))
                return(NULL)
            } else if (nrow(xsub) <= 2) {
                message(paste0("WARNING: Few genes (1 or 2) found (not clustering)"))
                row.order <- rownames(xsub)
            } else {
                hc <- hclust(dist(replace(xsub, is.na(xsub), 0)))
                sv <- svd(t(xsub))$v[,1]
                hc2 <- as.hclust(reorder(as.dendrogram(hc), wts = sv))
                row.order <- rownames(xsub)[hc2$order]
            }
            names(row.order) <- rep("ex_00", length(row.order))
        }
    } else {
        if (is.null(Genes)) {
            row.order <- unique(dexData$gene_id)[unique(dexData$gene_id) %in% rownames(coefs)]
            names(row.order) <- rep("ex_00", length(row.order))
        } else {
            row.order <- unique(Genes)
        }
    }
    
    coefs <- coefs[row.order,]
    zs <- zs[row.order,]
    pvals <- pvals[row.order,]
    padjs <- padjs[row.order,]
    
    #zs.sig <- zs[keep.sig.rows,]
    #coefs.sig <- coefs[keep.sig.rows,]
    main_cts <- ct.order[!rownames(ct.order) %in% c("ex", "in", "oligodendrocytes", "astrocytes"),, drop = FALSE]
    ha_c2 <- HeatmapAnnotation(cell_type = factor(c(main_cts$ct), levels = c(major.ct)), show_annotation_name = FALSE,
                          col = list(cell_type = setNames(brewer.pal(length(c(major.ct)), "Pastel1"), c(major.ct))))
    lgd2 <- Legend(labels = c("FDR < 5%", "p < 0.05"), title = "", 
                   graphics = list(
                                   function(x, y, w, h) grid.points(x, y, gp = gpar(col = "black", fill = "black"), size = unit(2, "mm"), pch = 21),
                                   function(x, y, w, h) grid.points(x, y, gp = gpar(col = "black"), size = unit(1.6, "mm"), pch = 4)))
    Rowsplit <- NULL
    if (Reorder) {Rowsplit <- factor(str_replace(names(row.order), "\\d+$", ""), levels = major.ct)}
    Layerfun <- NULL
    if (Reorder) {
        Layerfun <- function(j, i, x, y, width, height, fill) {
            if(str_replace(names(row.order), "\\d+$", "")[i][1] == c(main_cts$ct)[j][1]) {
                grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
                }
        }
    }
    Cap = ifelse(is.na(Cap), signif(max(abs(range(coefs[!is.na(coefs)]))), 2), Cap)
    message(paste0("CAP = ", Cap))
    if (is.infinite(Cap)) {
        warning("Capping at 2 because infinite values")
        Cap <- 2
    }


    ht2 <- Heatmap(coefs,
            column_title = ifelse(is.na(title), paste0(Diagnosis), title), 
            name = paste0("logFC (capped at +/- ", Cap, ")"),
            col = colorRamp2(c(-Cap, 0, Cap), c(cb_blue, "white", cb_red)), # For capping logFC
            na_col = "grey80",
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            rect_gp = gpar(col = "grey40", lwd = 1),
            row_split = Rowsplit,
            row_gap = unit(1, "mm"),
            row_names_side = genenames.side,
            row_names_gp = gpar(fontsize = 8, col = "black"),
            row_title = NULL,
            column_split = factor(c(main_cts$ct), levels = c(major.ct)),
            column_names_rot = 45,
            cluster_column_slices = FALSE,
            column_names_gp = gpar(fontsize = 8),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if (!(is.na(padjs[i,j])) & (padjs[i, j] < 0.05))
                    grid.points(x, y, pch = 21, size = unit(2, "mm"), gp = gpar(col = "black", fill = "black"))
                else if (!(is.na(pvals[i,j])) & (pvals[i, j] < 0.05))
                    grid.points(x, y, pch = 4, size = unit(1.6, "mm"))
            },
            layer_fun = Layerfun,
            top_annotation = ha_c2
    )
    return(list(plot = ht2, legend = lgd2, n_genes = nrow(coefs)))
}

