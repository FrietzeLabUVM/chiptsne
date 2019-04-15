#' plot_path
#'
#' Traces a path through the t-sne space in cell line order defined by qcells
#' for ids in id_to_plot
#'
#' Good for looking at a small number of ids in a modest number of cells.
#'
#' @param tsne_dt data.table with tsne info (tx, ty, id, and cell)
#' @param qcells character vector of items in tsne_dt$cell
#' @param id_to_plot character vector of ids in tsne_dt$id
#' @param p exiting ggplot to add a layer onto.  Default of NULL
#' creates a new ggplot.
#' @param xrng numeric of length 2. passed to coord_cartesian xlim.  Not used if p is specified.  Default is c(-.5, .5).
#' @param yrng numeric of length 2. passed to coord_cartesian ylim.  Not used if p is specified.  Default is c(-.5, .5).
#' @param arrowhead_position character, must be one of "each" or "end".
#'   Determines if arrowheads are drawn for each segment or only on the final
#'   segment.
#' @param line_type character vector describing type of line to connect qcells.
#'   One of : curve, spline, or straight
#' @param label_type character vector describing labelling method for points
#'   along lines.  One of : text, label, or none.
#' @param bg_points number of background id points to plot.
#' @param arrow_FUN result of grid::arrow().  Default of NULL does not draw arrowheads.
#'
#' @return ggplot showing how individual ids behave across qcells.
#' @export
#' @importFrom stats spline
#' @importFrom ggrepel geom_text_repel geom_label_repel
#'
#' @examples
#' data(tsne_dt)
#' plot_path(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3])
#' plot_path(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3],
#'     arrowhead_position = "each", label_type = "none")
#' plot_path(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3],
#'     arrowhead_position = "end", label_type = "none", line_type = "spline",
#'     arrow_FUN = arrow())
plot_path = function(tsne_dt,
                     qcells,
                     id_to_plot,
                     p = NULL,
                     xrng = c(-.5, .5),
                     yrng = c(-.5, .5),
                     arrowhead_position = c("end", "each")[1],
                     line_type = c("curve", "spline", "straight")[2],
                     label_type = c("text", "label", "none")[2],
                     bg_points = 5000,
                     arrow_FUN = NULL) {
    stopifnot(qcells %in% unique(tsne_dt$cell))
    stopifnot(arrowhead_position %in% c("end", "each"))
    stopifnot(line_type %in% c("curve", "spline", "straight"))
    stopifnot(label_type %in% c("text", "label", "none"))
    stopifnot(id_to_plot %in% tsne_dt$id)

    lines_dt = tsne_dt[cell %in% qcells & id %in% id_to_plot]

    lines_dt$cell = factor(lines_dt$cell, levels = qcells)
    lines_dt = lines_dt[order(cell)][order(id)][]
    lines_dt[, cell_o := seq(.N), by = list(id)]
    # lines_dt
    if(is.null(p)){
        p = ggplot() +
            geom_point(data = tsne_dt[sampleCap(seq(nrow(tsne_dt)), bg_points), ],
                       aes(x = tx, y = ty), color = "gray") +
            labs(title = paste(qcells, collapse = ", ")) +
            theme_classic() +
            scale_color_brewer(palette = "Dark2") +
            coord_cartesian(xlim = xrng, ylim = yrng)
    }

    switch(line_type,
           curve = {
               # plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)], list(tx, ty, id, cell_o)],
               #                 lines_dt[seq_along(qcells)[-1], list(tx_end = tx,
               #                                                      ty_end = ty,
               #                                                      id,
               #                                                      cell_o = cell_o - 1)])
               plot_dt = merge(lines_dt[cell_o != length(qcells), list(tx, ty, id, cell_o)],
                               lines_dt[cell_o != 1, list(tx_end = tx,
                                                          ty_end = ty,
                                                          id,
                                                          cell_o = cell_o - 1)])
               switch(arrowhead_position,
                      each = {
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
                      end = {
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
               sp_y = lines_dt[, stats::spline(x = cell_o,
                                        y = ty,
                                        n = n * (length(qcells) - 1)), by = id][, list(pid = seq(.N), ty = y), by = list(id)]
               sp_x = lines_dt[, stats::spline(x = cell_o,
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
               switch(arrowhead_position,
                      each = {
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
                      end = {
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
               switch(arrowhead_position,
                      each = {
                          plot_dt = merge(lines_dt[cell_o != length(qcells), list(tx, ty, id, cell_o)],
                                          lines_dt[cell_o != 1, list(tx_end = tx,
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
                                      color = id,
                                      group = id
                                  ),
                                  size = 1,
                                  arrow = arrow_FUN
                              )
                      },
                      end = {
                          plot_dt = lines_dt
                          p = p + geom_path(data = plot_dt,
                                            aes(x = tx, y = ty, color = id, group = id),
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





#' plot_outline
#'
#' a ggplot where the position of id in every cell specified by qcells is
#' connected in a polygon.  Allows the identification of both regions where ids
#' are stable/dynamic and individual ids that are particularly dynamic.
#'
#' Good for looking at large numbers of ids with a modest number of cells.
#'
#' @param tsne_dt data.table with tsne info (tx, ty, id, and cell)
#' @param qcells character vector of items in tsne_dt$cell
#' @param id_to_plot character vector of ids in tsne_dt$id
#' @param p exiting ggplot to add a layer onto.  Default of NULL creates a new
#'   ggplot.
#' @param xrng numeric of length 2. passed to coord_cartesian xlim.  Not used if
#'   p is specified.  Default is c(-.5, .5).
#' @param yrng numeric of length 2. passed to coord_cartesian ylim.  Not used if
#'   p is specified.  Default is c(-.5, .5).
#' @param bg_color character. color to use for background points. Default is
#'   "gray"
#' @param line_color_mapping character that is valid color. If less than length
#'   of id_to_plot, recycled across specified id_to_plot.  Can be named vector
#'   to completely specify id_to_plot.
#' @param fill_color_mapping character that is valid color. If less than length
#'   of id_to_plot, recycled across specified id_to_plot.  Can be named vector
#'   to completely specify id_to_plot.
#' @param label_type  character.  one of c("text", "label", "none").  controls
#'   how, if at all, plot objects are labelled.
#' @param bg_points number of points to plot in background.  if 0, only points
#'   corresponding to id_to_plot are drawn.  if -1, no points at all are drawn.
#' @param arrow_FUN result of grid::arrow().  Default of NULL does not draw arrowheads.
#'
#' @return a ggplot
#' @export
#' @importFrom grDevices chull
#'
#' @examples
#' data(tsne_dt)
#' plot_outline(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3])
#' plot_outline(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3],
#'     label_type = "none")
#' plot_outline(tsne_dt, unique(tsne_dt$cell), unique(tsne_dt$id)[1:3],
#'     label_type = "label")
plot_outline = function(tsne_dt,
                        qcells,
                        id_to_plot = NULL,
                        p = NULL,
                        xrng = c(-.5, .5),
                        yrng = c(-.5, .5),
                        bg_color = "gray",
                        line_color_mapping = "black",
                        fill_color_mapping = "gray",
                        label_type = c("text", "label", "none")[3],
                        bg_points = 5000,
                        arrow_FUN = NULL) {
    stopifnot(qcells %in% unique(tsne_dt$cell))
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

    ch_dt = lines_dt[, .(ch_i  = grDevices::chull(tx, ty)), .(id)]

    lines_dt[, ch_i := seq(.N), by = .(id)]
    ch_res = lines_dt[, .(ch_i  = grDevices::chull(tx, ty)), by = .(id)]
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
