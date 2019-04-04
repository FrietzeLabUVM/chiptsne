#' Title
#'
#' @param bfc
#' @param rname
#' @param FUN
#' @param force_overwrite
#' @param check_only
#'
#' @return
#' @export
#'
#' @examples
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE, check_only = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        message("results do exist.")
        if(check_only){
            return(TRUE)
        }
        message("loading results.")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        if(!file.exists(cache_path)){
            message("results do not exists.")
        }else{
            message("results being overwritten.")
        }
        if(check_only){
            return(FALSE)
        }
        message("executing function...")
        res = FUN()
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

#' Title
#'
#' @param gr
#' @param ref_gr
#' @param gr_size
#'
#' @return
#' @export
#'
#' @examples
my_annotate = function(gr,
                       ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz",
                                                        feature.type = "transcript", format = "gtf"),
                       gr_size = 2000){
    max_dist = 5e3
    dists = distanceToNearest(ref_gr, resize(gr, gr_size, fix = "center"))
    dists = subset(dists, distance <= max_dist)
    anno_dt = cbind(as.data.table(ref_gr[queryHits(dists)])[, .(gene_name, gene_id, transcript_id)],
                    as.data.table(gr[subjectHits(dists)]), distance = mcols(dists)[[1]])
    anno_dt
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
symbol2uniprot = function(x){
    bres = bitr(x, fromType = "SYMBOL",
                toType = c("UNIPROT"),
                OrgDb = org.Hs.eg.db)
    bres[!duplicated(bres$SYMBOL),]$UNIPROT
}

#' Title
#'
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @examples
sampleCap = function(x, n = 500){
    n = min(n, length(unique(x)))
    out = sample(unique(x), n)
    if(is.factor(out)) out = as.character(out)
    out
}

#' Title
#'
#' @param x
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
norm1 = function(x, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    (x - min(xrng)) / (max(xrng) - min(xrng))
}

#' Title
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
mybin = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(norm1(x, xrng) * (n_points-.00001))+1
}

#' Title
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
mybin_centers = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    # xrng = range(x)
    xspc = diff(xrng)/n_points/2
    xs = seq(min(xrng)+xspc, max(xrng)-xspc, diff(xrng)/(n_points))
    xs
}

movingAverage = function(x, n = 1, centered = TRUE) {

    if (centered) {
        before <- floor((n - 1)/2)
        after <- ceiling((n - 1)/2)
    } else {
        before <- n - 1
        after <- 0
    }

    # Track the sum and count of number of non-NA items
    s <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data
    new <- x
    # Add to count list wherever there isn't a
    count <- count + (!is.na(new))
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new <- c(rep(NA, i), x[seq_len(length(x) - i)])

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new <- c(x[(i + 1):length(x)], rep(NA, i))

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # return sum divided by count
    s/count
}

#' Title
#'
#' @param p
#' @param rects
#'
#' @return
#' @export
#'
#' @examples
annotate_rects = function(p, rects){
    anns = matrix(unlist(rects), ncol = 4, byrow = TRUE)
    p_alpha = .1

    for(i in seq_len(nrow(anns))){
        p = p +
            annotate("rect",
                     xmin = anns[i, 1],
                     xmax = anns[i, 2],
                     ymin = anns[i, 3],
                     ymax = anns[i, 4], color = "black", fill = NA) +
            annotate("label",
                     # x = (anns[i, 1] + anns[i, 2])/2,
                     # y = (anns[i, 3] + anns[i, 4])/2,
                     x = anns[i, 1],
                     y = anns[i, 4],
                     hjust = 1,
                     vjust = 1,
                     label = i)


    }
    p
}


#' Title
#'
#' @param tsne_res
#' @param cell_a
#' @param cell_b
#' @param n_points
#'
#' @return
#' @export
#'
#' @examples
calc_delta = function(tsne_res, cell_a, cell_b, n_points){
    v_dt = dcast(tsne_res[cell %in% c(cell_a, cell_b)], "id~cell", value.var = c("tx", "ty"))
    colnames(v_dt) = sub(cell_a, "cell_a", colnames(v_dt))
    colnames(v_dt) = sub(cell_b, "cell_b", colnames(v_dt))
    v_dt$bx_cell_a = mybin(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    xs = mybin_centers(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    v_dt$btx_cell_a = xs[v_dt$bx_cell_a]

    v_dt$by_cell_a = mybin(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    ys = mybin_centers(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    v_dt$bty_cell_a = ys[v_dt$by_cell_a]

    av_dt = v_dt[, list(tx_cell_b = mean(tx_cell_b), ty_cell_b = mean(ty_cell_b), N = .N), list(bx_cell_a, by_cell_a)]
    av_dt$tx_cell_a = xs[av_dt$bx_cell_a]
    av_dt$ty_cell_a = ys[av_dt$by_cell_a]
    return(list(velocity_dt = v_dt, agg_velocity_dt = av_dt))
}


#' Title
#'
#' @param qcell
#' @param as_facet
#'
#' @return
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @export
#'
#' @examples
make_tss_plot = function(qcell, as_facet = TRUE){
    xrng = range(tsne_res$tx)
    yrng = range(tsne_res$ty)

    pr_img_res = make_tsne_img(
        profiles_dt = pr_dt,
        position_dt = tsne_res[cell %in% qcell],
        apply_norm = FALSE,
        #force_rewrite = TRUE,
        xrng = xrng,
        yrng = yrng,
        n_points = n_points, line_colors = c("signal" = "black")
    )
    if(as_facet){
        p_pr_density = plot_tsne_img_byCell(pr_img_res$images_dt,
                                            pr_img_res$tsne_dt[cell %in% qcell],
                                            n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
        p_h7_density = plot_tsne_img_byCell(img_res$images_dt,
                                            img_res$tsne_dt[cell %in% qcell],
                                            n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }else{
        p_pr_density = plot_tsne_img(pr_img_res$images_dt,
                                     # pr_img_res$tsne_dt[cell %in% qcell],
                                     n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
        p_h7_density = plot_tsne_img(img_res$images_dt,
                                     # img_res$tsne_dt[cell %in% qcell],
                                     n_points = n_points, N_ceiling = NULL)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }


    pg = cowplot::plot_grid(p_h7_density + labs(title = paste(paste(qcell, collapse = ", "), ": k4me3+k27me3")),
                            p_pr_density + labs(title = paste(paste(qcell, collapse = ", "), ": tss frequency")))
    # ggsave(paste0("tmp_", qcell, "_tss.pdf"), plot = pg, width = 8, height = 4)
    pg
    # head(pr_img_res$tsne_dt[cell == cell, .N, by = list(bx, by)][order(N, decreasing = TRUE)])
    # head(img_res$tsne_dt[cell == cell, .N, by = list(bx, by)][order(N, decreasing = TRUE)])
}

#' Title
#'
#' @param data_dt
#' @param qgr
#' @param qcells
#' @param tss_ids
#'
#' @return
#' @export
#'
#' @examples
plot_profiles_selected = function(data_dt, qgr, qcells, tss_ids,
                                  color_mapping = c("H3K27me3" = "firebrick",
                                                    "H3K4me3" = "forestgreen")){
    # qcells = c("H7", "CD34", "Kasumi1", "mm1s", "Nalm6")
    # tss_ids = subset(qgr, gene_name == "RUNX1")$id[1]
    stopifnot(qcells %in% unique(data_dt$cell))
    if(!tss_ids %in% data_dt$id){
        tmp = unlist(strsplit(tss_ids, " "))
        if(length(tmp) > 1){
            tss_ids = tmp[1]
            tmp = as.numeric(tmp[-1])
            tss_ids = subset(qgr, gene_name %in% tss_ids)$id[tmp]
        }else{
            tss_ids = subset(qgr, gene_name %in% tss_ids)$id
        }

    }
    names(qgr) = qgr$id
    message(paste(as.character(qgr[tss_ids]), collapse = "\n"))
    plot_dt = data_dt[id %in% tss_ids & cell %in% qcells]
    plot_dt$cell = factor(plot_dt$cell, levels = qcells)
    plot_dt$id = factor(plot_dt$id, levels = tss_ids)

    p = ggplot(plot_dt, aes(x = x, ymin = 0, ymax = y, y = y, color = mark, fill = mark)) +
        facet_grid("cell~id", switch = "y") +
        geom_ribbon(alpha = .5) +
        geom_path(show.legend = FALSE) +
        theme_classic() +
        theme(strip.background = element_blank(), strip.placement = "outside",
              strip.text.y = element_text(angle = 180)) +
        labs(x = "", y = "") +
        scale_color_manual(values = color_mapping) +
        scale_fill_manual(values = color_mapping)
    p
}
