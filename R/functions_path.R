

#' Title
#'
#' @param tsne_dt
#' @param qgr
#' @param qcells
#' @param tss_ids
#' @param grp_var
#' @param line_type
#' @param label_type
#'
#' @return
#' @export
#'
#' @examples
plot_path = function(tsne_dt,
                     qgr,
                     qcells,
                     tss_ids,
                     grp_var = c("id", "grp")[1],
                     line_type = c("curve", "spline", "straight")[2],
                     label_type = c("text", "label", "none")[2],
                     bg_points = 5000,
                     arrow_FUN = NULL) {
    stopifnot(qcells %in% unique(tsne_dt$cell))
    names(qgr) = qgr$id
    stopifnot(tss_ids %in% tsne_dt$id)

    lines_dt = tsne_dt[cell %in% qcells & id %in% tss_ids]

    lines_dt$cell = factor(lines_dt$cell, levels = qcells)
    lines_dt = lines_dt[order(cell)][order(id)][]
    lines_dt[, cell_o := seq(.N), by = list(id)]
    # lines_dt
    lines_dt$id = as.character(qgr[lines_dt$id])
    p =     ggplot() +
        geom_point(data = tsne_dt[sampleCap(seq(nrow(tsne_dt)), bg_points), ],
                   aes(x = tx, y = ty), color = "gray") +
        labs(title = paste(qcells, collapse = ", ")) +
        theme_classic() +
        scale_color_brewer(palette = "Dark2")
    switch(line_type,
           curve = {
               plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)], list(tx, ty, id, cell_o)],
                               lines_dt[seq_along(qcells)[-1], list(tx_end = tx,
                                                                    ty_end = ty,
                                                                    id,
                                                                    cell_o = cell_o - 1)])
               switch(grp_var,
                      grp = {
                          p = p +
                              geom_curve(
                                  data = plot_dt,
                                  aes(
                                      x = tx,
                                      y = ty,
                                      xend = tx_end,
                                      yend = ty_end,
                                      color = id
                                  ),
                                  size = 1,
                                  arrow = arrow_FUN
                              )
                      },
                      id = {
                          p = p +
                              geom_curve(
                                  data = plot_dt[cell_o < max(cell_o)],
                                  aes(
                                      x = tx,
                                      y = ty,
                                      xend = tx_end,
                                      yend = ty_end,
                                      color = id
                                  ),
                                  size = 1
                              ) +
                              geom_curve(
                                  data = plot_dt[cell_o == max(cell_o)],
                                  aes(
                                      x = tx,
                                      y = ty,
                                      xend = tx_end,
                                      yend = ty_end,
                                      color = id
                                  ),
                                  size = 1,
                                  arrow = arrow_FUN
                              )
                      })
           },
           spline = {
               n = 20
               sp_y = lines_dt[, spline(x = cell_o,
                                        y = ty,
                                        n = n * (length(qcells) - 1)), by = id][, list(pid = seq(.N), ty = y), by = list(id)]
               sp_x = lines_dt[, spline(x = cell_o,
                                        y = tx,
                                        n = n * (length(qcells) - 1)), by = id][, list(pid = seq(.N), tx = y), by = list(id)]
               sp_dt = merge(sp_x, sp_y, by = c("id", "pid"))
               ceiling(sp_dt$pid / n)

               sp_dt[, grp := ceiling(pid / n)]
               sp_dt[, grp_o := seq(.N), by = list(grp, id)]
               start_dt = merge(lines_dt[cell_o < length(qcells), list(tx, ty, grp = cell_o, id)],
                                unique(sp_dt[, list(id, grp)]))[, grp_o := 0]
               end_dt = merge(lines_dt[cell_o > 1 &
                                           cell_o < length(qcells), list(tx, ty, grp = cell_o - 1, id)],
                              unique(sp_dt[, list(id, grp = grp)]))[, grp_o := n +
                                                                        1]
               plot_dt = rbind(sp_dt[, list(grp, id, tx, ty, grp_o)],
                               start_dt,
                               end_dt)[order(grp_o)][order(id)][order(grp)]
               switch(grp_var,
                      grp = {
                          p = p +
                              geom_path(
                                  data = plot_dt,
                                  aes(
                                      x = tx,
                                      y = ty,
                                      color = id,
                                      group = paste(grp, id)
                                  ),
                                  arrow = arrow_FUN,
                                  size = 1.2,
                                  alpha = 1,
                                  show.legend = FALSE
                              )
                      },
                      id = {
                          p = p +
                              geom_path(
                                  data = plot_dt,
                                  aes(
                                      x = tx,
                                      y = ty,
                                      color = id,
                                      group = id
                                  ),
                                  arrow = arrow_FUN,
                                  size = 1.2,
                                  alpha = 1,
                                  show.legend = FALSE
                              )
                      })

           },
           straight = {
               switch(grp_var,
                      grp = {
                          plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)], list(tx, ty, id, cell_o)],
                                          lines_dt[seq_along(qcells)[-1], list(tx_end = tx,
                                                                               ty_end = ty,
                                                                               id,
                                                                               cell_o = cell_o - 1)])
                          p = p +
                              geom_segment(
                                  data = plot_dt,
                                  aes(
                                      x = tx,
                                      y = ty,
                                      xend = tx_end,
                                      yend = ty_end,
                                      color = id
                                  ),
                                  size = 1,
                                  arrow = arrow_FUN
                              )
                      },
                      id = {
                          plot_dt = lines_dt
                          p = p + geom_path(data = plot_dt,
                                            aes(x = tx, y = ty),
                                            arrow = arrow_FUN)
                      })

           })
    p = p + geom_point(
        data = lines_dt,
        aes(x = tx, y = ty, color = id),
        size = 3,
        shape = 21,
        fill = "white"
    )
    switch(label_type,
           text = {
               p = p + ggrepel::geom_text_repel(
                   data = lines_dt,
                   aes(
                       x = tx,
                       y = ty,
                       color = id,
                       label = cell
                   ),
                   show.legend = FALSE
               )
           },
           label = {
               p = p + ggrepel::geom_label_repel(
                   data = lines_dt,
                   aes(
                       x = tx,
                       y = ty,
                       color = id,
                       label = cell
                   ),
                   fill = "white",
                   show.legend = FALSE
               )
           },
           none = {
               p = p
           })
    p
}

#' Title
#'
#' @param tsne_dt
#' @param qgr
#' @param qcells
#' @param tss_ids
#' @param grp_var
#' @param xrng
#' @param yrng
#' @param line_type
#' @param label_type
#' @param bg_points
#' @param arrow_FUN
#'
#' @return
#' @export
#'
#' @examples
plot_outline = function(tsne_dt,
                        qgr,
                        qcells,
                        p = NULL,
                        id_to_plot = NULL,
                        xrng = c(-.5, .5),
                        yrng = c(-.5, .5),
                        bg_color = "gray",
                        line_color_mapping = "black",
                        fill_color_mapping = "gray",
                        label_type = c("text", "label", "none")[3],
                        bg_points = 5000,
                        arrow_FUN = NULL) {
    stopifnot(qcells %in% unique(tsne_dt$cell))
    names(qgr) = qgr$id
    if(is.numeric(label_type)){
        label_type = c("text", "label", "none")[label_type]
    }
    if(is.null(id_to_plot)){
        id_to_plot = unique(tsne_dt$id)
    }
    stopifnot(id_to_plot %in% tsne_dt$id)

    lines_dt = tsne_dt[cell %in% qcells & id %in% id_to_plot]

    lines_dt$cell = factor(lines_dt$cell, levels = qcells)
    lines_dt = lines_dt[order(cell)][order(id)]

    lo = (seq(id_to_plot) %% length(line_color_mapping))+1
    line_color_mapping = line_color_mapping[lo]
    names(line_color_mapping) = id_to_plot

    fo = (seq(id_to_plot) %% length(fill_color_mapping))+1
    fill_color_mapping = fill_color_mapping[fo]
    names(fill_color_mapping) = id_to_plot

    # lines_dt
    if(bg_points < 0){
        id_tp = character()
    }else if(bg_points == 0){
        id_tp = id_to_plot
    }else{
        id_tp = sampleCap(unique(tsne_dt$id), bg_points)
        id_tp = union(id_tp, id_to_plot)
    }

    if(is.null(p)){
        p = ggplot() +
            labs(title = paste(qcells, collapse = ", ")) +
            theme_classic() +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }
    p = p +
        annotate("point",
                 x = tsne_dt[id %in% id_tp,]$tx,
                 y = tsne_dt[id %in% id_tp,]$ty,
                 color = bg_color)

    ch_dt = lines_dt[, .(ch_i  = chull(tx, ty)), .(id)]

    lines_dt[, ch_i := seq(.N), by = .(id)]
    ch_res = lines_dt[, .(ch_i  = chull(tx, ty)), by = .(id)]
    ch_res$o = seq(nrow(ch_res))
    poly_dt = merge(lines_dt, ch_res)
    poly_dt = poly_dt[order(o)]

    for(tid in unique(poly_dt$id)){
        p = p +
            annotate("polygon",
                     x = poly_dt[id == tid]$tx,
                     y = poly_dt[id == tid]$ty,
                     color = line_color_mapping[tid],
                     fill = fill_color_mapping[tid])
    }
    lab_dt = lines_dt[, .(tx = mean(tx), ty = mean(ty)), by = .(id)]
    switch(label_type,
           text = {
               p = p + ggrepel::geom_text_repel(
                   data = lab_dt,
                   aes(
                       x = tx,
                       y = ty,
                       label = id
                   ),
                   color = "black",
                   show.legend = FALSE
               )
           },
           label = {
               p = p + ggrepel::geom_label_repel(
                   data = lab_dt,
                   aes(
                       x = tx,
                       y = ty,
                       label = id
                   ),
                   color = "black",
                   fill = "white",
                   show.legend = FALSE
               )
           },
           none = {
               p = p
           })
    p
}
