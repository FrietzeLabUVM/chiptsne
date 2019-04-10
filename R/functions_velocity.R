# enforce domain of [0,360)
ang_cap = function(angle){
    angle%%360
}

# return angle between (x1, y1) and (x2, y2)
xy2deg = function(x1, y1, x2, y2){
    x = x2 - x1
    y = y2 - y1
    deg = atan(y/x) * 180 / pi + 180
    deg[x < 0] = deg[x < 0] + 180
    # deg[deg > 360] = deg[deg > 360] - 360
    deg = 360 - (deg + 90)
    ang_cap(deg)
}

# return euclidean distance between (x1, y1) and (x2, y2)
xy2dist = function(x1, y1, x2, y2){
    x = x2 - x1
    y = y2 - y1
    (x^2 + y^2)^.5
}

#' Title
#'
#' @param tsne_dt
#' @param cell_a
#' @param cell_b
#' @param n_points
#' @param p
#' @param min_N
#' @param id_to_plot
#'
#' @return
#' @export
#'
#' @examples
plot_regional_velocity = function(tsne_dt,
                                  cell_a, cell_b,
                                  n_points,
                                  points_as_background = FALSE,
                                  p = NULL,
                                  min_N = 0,
                                  id_to_plot = NULL,
                                  strategy = c("by_direction", "by_destination", "individual_recentered")[2]
){
    if(is.numeric(strategy)){
        strategy = c("by_destination", "by_direction", "individual_recentered")[strategy]
    }
    if(is.null(id_to_plot)){
        tsne_dt.tp = copy(tsne_dt)
    }else{
        tsne_dt.tp = tsne_dt[id %in% id_to_plot]
    }

    av_dt = calc_delta(tsne_dt.tp, cell_a, cell_b, n_points, strategy = strategy)
    # if(is.null(av_dt$N)) av_dt$N = 1
    if (is.null(p))
        p = ggplot()
    if (points_as_background) {
        p = p + geom_point(data = tsne_dt.tp[cell %in% c(cell_a, cell_b)],
                           aes(x = tx, y = ty, color = cell))
    }
    if (is.null(av_dt$N)) {
        p = p + geom_segment(
            data = av_dt,
            aes(
                x = tx_cell_a,
                xend = tx_cell_b,
                y = ty_cell_a,
                yend = ty_cell_b
            ),
            arrow = arrow(length = unit(x = .02, units = "npc"))
        ) +
            coord_cartesian(xlim = range(tsne_dt$tx), ylim = range(tsne_dt$ty)) +
            scale_fill_gradient2(low = "gray", high = "black") +
            labs(x = "x", y = "y", fill = "density") +
            theme_classic()
    }else{
        p = p + geom_segment(
            data = av_dt[N >= min_N],
            aes(
                x = tx_cell_a,
                xend = tx_cell_b,
                y = ty_cell_a,
                yend = ty_cell_b,
                size = N
            ),
            arrow = arrow(length = unit(x = .02, units = "npc"))
        ) +
            coord_cartesian(xlim = range(tsne_dt$tx), ylim = range(tsne_dt$ty)) +
            scale_fill_gradient2(low = "gray", high = "black") +
            labs(x = "x", y = "y", fill = "density") +
            scale_size_continuous(range = c(.5, 2),
                                  breaks = range(av_dt$N)) +
            theme_classic()
    }

    # p_velocity = ggplot() +
    # geom_density2d(data = tsne_dt, aes(x = tsne_dt$tx, y = tsne_dt$ty, color = "lightgray") +
    # geom_density2d(data = tsne_dt, aes(x = tx, y = ty), color = "lightgray") +
    # stat_density_2d(data = tsne_dt, aes(x = tx, y = ty, fill = stat(level)), geom = "polygon", bins = 7) +
    # geom_point(data = tsne_dt.tp[cell %in% c(cell_a, cell_b)][id %in% sampleCap(id, 500)], aes(x = tx, y = ty, color = cell)) +
    p
}

#' prep_velocity
#'
#' @param tsne_dt
#' @param cell_a
#' @param cell_b
#' @param id_to_plot
#' @param max_plotted
#' @param return_data
#' @param delta.min
#' @param delta.max
#' @param angle.min
#' @param angle.max
#' @param drop_backgroud
#'
#' @return
#' @export
#'
#' @examples
prep_velocity = function(tsne_dt,
                         cell_a,
                         cell_b,
                         id_to_plot = NULL,
                         max_plotted = 500,
                         return_data = FALSE,
                         delta.min = 0,
                         delta.max = Inf,
                         angle.min = 0,
                         angle.max = 360,
                         drop_backgroud = FALSE) {

    v_dt = calc_delta(tsne_dt, cell_a, cell_b, 10)

    # manual selection
    if(is.null(id_to_plot)){
        v_dt.tp = copy(v_dt)
    }else{
        v_dt.tp = v_dt[id %in% id_to_plot]
    }
    # random selection
    if(max_plotted < Inf){
        v_dt.tp = v_dt.tp[id %in% sampleCap(v_dt.tp$id, max_plotted)]
    }
    # angle selection
    v_dt.tp[, angle := xy2deg(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
    if(angle.min > angle.max){
        v_dt.tp[, foreground := angle <= angle.min & angle >= angle.max]
    }else{
        v_dt.tp[, foreground := angle >= angle.min & angle <= angle.max]
    }
    # distance selection
    v_dt.tp[, distance := xy2dist(x1 = tx_cell_a, x2 = tx_cell_b, y1 = ty_cell_a, y2 = ty_cell_b)]
    v_dt.tp = v_dt.tp[distance < delta.min | distance > delta.max, foreground := FALSE]
    if(drop_backgroud){
        v_dt.tp = v_dt.tp[foreground == TRUE]
    }
    v_dt.tp[]
}

#' Title
#'
#' @param velocity_dt
#'
#' @return
#' @export
#'
#' @examples
plot_velocity_centered = function(velocity_dt){
    ggplot(velocity_dt, aes(x = angle, xend = angle, y = 0, yend = distance, color = angle)) +
        geom_segment(stat = "identity", alpha = .5, size = 1) +
        coord_polar() +
        scale_color_gradientn(colours = c("orange", "red", "purple", "blue",
                                          "green", "orange"), limits = c(0, 360), breaks = 0:4*90) +
        scale_x_continuous(limits = c(0, 360), breaks = 0:4*90)
}

#' Title
#'
#' @param velocity_dt
#'
#' @return
#' @export
#'
#' @examples
plot_velocity_bins = function(velocity_dt, bins = 36, bin_FUN = list(length, sum, mean)[[1]]){#use_distance = FALSE){
    velocity_dt[, angle_bin := ceiling((angle)/(360/bins))]
    # if(use_distance){
    #     b_dt = velocity_dt[, .(N = sum(distance)), angle_bin]
    # }else{
    #     b_dt = velocity_dt[, .N, angle_bin]
    # }
    if(is.null(bin_FUN)){
        bin_FUN = length
    }
    b_dt = velocity_dt[, .(bin_value = bin_FUN(distance)), angle_bin]

    p_key = ggplot(b_dt, aes(x = angle_bin, y = bin_value, fill = angle_bin)) +
        geom_bar(width = 1, stat = "identity") + coord_polar() +
        scale_fill_gradientn(colours = c("orange", "red", "purple", "blue",
                                         "green", "orange"), limits = c(0, 360)/(360/bins),
                             breaks = 0:4*90/(360/bins), labels = function(x)x*360/bins) +
        scale_x_continuous(labels = function(x)x*360/bins, breaks = 0:4*90/(360/bins), limits = c(0, bins))
    p_key
}

#' Title
#'
#' @param tsne_dt
#' @param cell_a
#' @param cell_b
#' @param p
#' @param id_to_plot
#' @param max_plotted
#' @param delta.min
#' @param delta.max
#' @param angle.min
#' @param angle.max
#'
#' @return
#' @export
#'
#' @examples
plot_velocity_arrows = function(velocity_dt,
                                p = NULL,
                                id_to_plot = NULL,
                                max_plotted = 500,
                                delta.min = 0,
                                delta.max = Inf,
                                angle.min = 0,
                                angle.max = 360){
    bg = velocity_dt[foreground == FALSE & distance >= delta.min]
    if(is.null(p)) p = ggplot()
    p_arrows = p +
        annotate("segment",
                 x = bg$tx_cell_a, xend = bg$tx_cell_b,
                 y = bg$ty_cell_a, yend = bg$ty_cell_b,
                 color = "lightgray") +
        geom_segment(data = velocity_dt[foreground == TRUE],
                     aes(x = tx_cell_a, xend = tx_cell_b,
                         y = ty_cell_a, yend = ty_cell_b,
                         color = angle),
                     arrow = arrow(length = unit(0.1,"cm"))) +

        scale_color_gradientn(colours = c("orange", "red", "purple", "blue",
                                          "green", "orange"), limits = c(0, 360), breaks = 0:4*90)
    p_arrows
}

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
plot_velocity_arrows_selected = function(tsne_dt, qgr, qcells, tss_ids,
                                         grp_var = c("id", "grp")[1],
                                         line_type = c("curve", "spline", "straight")[2],
                                         label_type = c("text", "label", "none")[2]){
    # qcells = c("H7", "CD34", "Kasumi1", "mm1s", "Nalm6")
    # tss_ids = subset(qgr, gene_name == "RUNX1")$id[1]
    stopifnot(qcells %in% unique(tsne_dt$cell))
    if(!tss_ids %in% tsne_dt$id){
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
    stopifnot(tss_ids %in% tsne_dt$id)

    lines_dt = tsne_dt[cell %in% qcells & id %in% tss_ids]

    lines_dt$cell = factor(lines_dt$cell, levels = qcells)
    lines_dt = lines_dt[order(cell)][order(id)][]
    lines_dt[, cell_o := seq(.N), by = list(id)]
    # lines_dt
    lines_dt$id = as.character(qgr[lines_dt$id])


    # lines_dt$tx = scales::rescale(c(1,1, 2,2,1,1:5), to = c(-.5,.5))
    # lines_dt$ty = scales::rescale(c(0,1, 1.5,2.5,3,1:5), to = c(-.5,.5))
    #
    #

    # lines_dt[, list(tx = spline(x = pid, y = tx, n = n*(length(qcells)-1)),
    # ty = spline(x = pid, y = ty, n = n*(length(qcells)-1))), by = id]
    p =     ggplot() +
        geom_point(data = tsne_dt[sample(seq(nrow(tsne_dt)), 5000),],
                   aes(x = tx, y = ty), color = "gray") +
        labs(title = paste(qcells, collapse = ", ")) +
        theme_classic() +
        scale_color_brewer(palette = "Dark2")
    switch(line_type,
           curve = {
               plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)],list(tx, ty, id, cell_o)],
                               lines_dt[seq_along(qcells)[-1], list(tx_end = tx, ty_end = ty, id, cell_o = cell_o -1)])
               switch(grp_var,
                      grp = {
                          p = p +
                              geom_curve(data = plot_dt,
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1, arrow = arrow())
                      },
                      id = {
                          p = p +
                              geom_curve(data = plot_dt[cell_o < max(cell_o)],
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1) +
                              geom_curve(data = plot_dt[cell_o == max(cell_o)],
                                         aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                         size = 1, arrow = arrow())
                      })
           },
           spline = {
               n = 20
               sp_y = lines_dt[, spline(x = cell_o, y = ty,
                                        n = n*(length(qcells)-1)), by = id][
                                            , list(pid = seq(.N), ty = y), by = list(id)]
               sp_x = lines_dt[, spline(x = cell_o, y = tx,
                                        n = n*(length(qcells)-1)), by = id][
                                            , list(pid = seq(.N), tx = y), by = list(id)]
               sp_dt = merge(sp_x, sp_y, by = c("id", "pid"))
               ceiling(sp_dt$pid/n)

               sp_dt[, grp := ceiling(pid / n)]
               sp_dt[, grp_o := seq(.N), by = list(grp, id)]
               start_dt = merge(lines_dt[cell_o < length(qcells), list(tx, ty, grp = cell_o, id)],
                                unique(sp_dt[, list(id, grp)]))[, grp_o := 0]
               end_dt = merge(lines_dt[cell_o > 1 & cell_o < length(qcells), list(tx, ty, grp = cell_o-1, id)],
                              unique(sp_dt[, list(id, grp = grp)]))[, grp_o := n+1]
               plot_dt = rbind(
                   sp_dt[, list(grp, id, tx, ty, grp_o)],
                   start_dt,
                   end_dt)[order(grp_o)][order(id)][order(grp)]
               switch(grp_var,
                      grp = {
                          p = p +
                              geom_path(data = plot_dt,
                                        aes(x = tx, y = ty, color = id, group = paste(grp,id)),
                                        arrow = arrow(),
                                        size = 1.2, alpha = 1,
                                        show.legend = FALSE)
                      },
                      id = {
                          p = p +
                              geom_path(data = plot_dt,
                                        aes(x = tx, y = ty, color = id, group = id),
                                        arrow = arrow(),
                                        size = 1.2, alpha = 1,
                                        show.legend = FALSE)
                      })

           },
           straight = {
               switch(grp_var,
                      grp = {
                          plot_dt = merge(lines_dt[seq_along(qcells)[-length(qcells)],list(tx, ty, id, cell_o)],
                                          lines_dt[seq_along(qcells)[-1], list(tx_end = tx, ty_end = ty, id, cell_o = cell_o -1)])
                          p = p +
                              geom_segment(data = plot_dt,
                                           aes(x = tx, y = ty, xend = tx_end, yend = ty_end, color = id),
                                           size = 1, arrow = arrow())
                      },
                      id = {
                          plot_dt = lines_dt
                          p = p + geom_path(data = plot_dt, aes(x = tx, y = ty), arrow = arrow())
                      })

           })
    p = p + geom_point(data = lines_dt,
                       aes(x = tx, y = ty, color = id),
                       size = 3, shape = 21, fill = "white")
    switch(label_type,
           text = {
               p = p + ggrepel::geom_text_repel(data = lines_dt,
                                                aes(x = tx, y = ty, color = id, label = cell),
                                                show.legend = FALSE)
           },
           label = {
               p = p + ggrepel::geom_label_repel(data = lines_dt,
                                                 aes(x = tx, y = ty, color = id, label = cell),
                                                 fill = "white", show.legend = FALSE)
           },
           none = {
               p = p
           })
    p
}
