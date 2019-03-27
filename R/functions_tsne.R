library(Rtsne)
library(ggimage)

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
#' @param qdt
#' @param qgr
#' @param qwin
#' @param qmet
#' @param cap_value
#' @param high_on_right
#' @param bfc
#' @param n_cores
#' @param rname
#' @param force_overwrite
#'
#' @return
#' @export
#'
#' @examples
fetch_tsne_mat = function(qdt,
                          qgr,
                          qwin = 50,
                          qmet = "summary",
                          cap_value = 20,
                          high_on_right = TRUE,
                          bfc = BiocFileCache::BiocFileCache(),
                          n_cores = getOption("mc.cores", 1),
                          rname = digest::digest(list(qgr, qdt[, 1:3], qwin, qmet, cap_value, high_on_right)),
                          force_overwrite = FALSE){
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(is.null(qdt$norm_factor)){
        qdt$norm_factor = 1
    }
    if(!all(file.exists(qdt[[1]]))){
        stop(paste(sep = "\n", "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt[[1]][!file.exists(qdt[[1]])])))

    }

    bw_fetch = function(){
        ssvFetchBigwig(qdt[, 1:3], qgr,
                       return_data.table = TRUE,
                       win_method = qmet,
                       win_size = qwin, n_cores = n_cores)
    }

    bw_dt = bfcif(bfc, rname, bw_fetch, force_overwrite = force_overwrite)
    bw_dt$sample = NULL
    bw_dt = bw_dt[, list(y = mean(y)), list(cell, id, mark, x)]
    if(!all(qdt$norm_factor == 1)){
        norm_dt = unique(qdt[, list(cell, mark, norm_factor)])
        bw_dt = merge(bw_dt, norm_dt)
        bw_dt[, y := y * norm_factor]
        bw_dt$norm_factor = NULL
    }

    bw_dt[y > cap_value, y := cap_value]



    if(high_on_right){
        balance_dt = bw_dt[, list(right_sum = sum(y[x > 0]), left_sum = sum(y[x < 0])), by = list(cell, id)]
        balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum, cell, id)]
        most_flipped = balance_dt[, list(fraction_flipped = sum(needs_flip) / .N), by = list(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        strand(qgr) = "+"
        strand(qgr)[most_flipped$flip_strand] = "-"
        bw_dt = merge(bw_dt, balance_dt, by = c("id", "cell"))
        remove(balance_dt)
        bw_dt[needs_flip == TRUE, x := -x]
        bw_dt$needs_flip = NULL
    }

    bw_dt$x = round(bw_dt$x, 3)

    for(m in unique(qdt$mark)){
        if(!exists("...tsne_mat")){
            dt = dcast( bw_dt[mark == m], id+cell~x, value.var = "y")
            ...tsne_mat = as.matrix(dt[, -1:-2])
            rn = paste(dt$id, dt$cell)
        }else{
            dt = dcast( bw_dt[mark == m], id+cell~x, value.var = "y")
            stopifnot(all(paste(dt$id, dt$cell) == rn))
            ...tsne_mat = cbind(...tsne_mat, as.matrix(dt[, -1:-2]))


        }
    }
    rownames(...tsne_mat) = rn
    return(list(bw_dt = bw_dt, tsne_mat = ...tsne_mat, query_gr = qgr))
}

#' Title
#'
#' @param tsne_mat
#' @param perplexity
#' @param n_cores
#' @param high_topright
#' @param norm1
#' @param bfc
#' @param rname
#' @param force_overwrite
#'
#' @return
#' @export
#' @import Rtsne
#' @examples
run_tsne = function(tsne_mat, perplexity = 100,
                    n_cores = getOption("mc.cores", 1),
                    high_topright = TRUE,
                    norm1 = TRUE,
                    bfc = BiocFileCache::BiocFileCache(),
                    rname = digest::digest(list(
                        tsne_mat[sample(1:nrow(tsne_mat), 20),
                                 sample(1:ncol(tsne_mat), 20)],
                        perplexity
                    )),
                    force_overwrite = FALSE){
    set.seed(0)
    res_tsne = bfcif(bfc, rname, force_overwrite = force_overwrite,
                     FUN = function(){
                         Rtsne(tsne_mat, num_threads = n_cores, perplexity = perplexity, check_duplicates = FALSE)
                     })

    tdt = as.data.table(res_tsne$Y)
    colnames(tdt) = c("tx", "ty")
    tdt$rn = rownames(tsne_mat)
    tdt[, c("id", "cell") := tstrsplit(rn, " ", keep = 1:2)]

    if(norm1){
        tdt$tx = norm1(tdt$tx)-.5
        tdt$ty = norm1(tdt$ty)-.5
    }


    if(high_topright){
        rs = rowSums(tsne_mat)
        tdt$rs = rs[tdt$rn]
        x_cutoff = mean(range(tdt$tx))
        x_flip = sum(tdt[tx > x_cutoff]$rs) < sum(tdt[tx < x_cutoff]$rs)
        if(x_flip){
            tdt[, tx := max(tx) - tx + min(tx)]
        }
        y_cutoff = mean(range(tdt$ty))
        y_flip = sum(tdt[ty > y_cutoff]$rs) < sum(tdt[ty < y_cutoff]$rs)
        if(y_flip){
            tdt[, ty := max(ty) - ty + min(ty)]
        }
        tdt$rs = NULL
    }

    tdt
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
            tss_ids = subset(qgr, gene_name == tss_ids)$id[tmp]
        }else{
            tss_ids = subset(qgr, gene_name == tss_ids)$id
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






