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
plot_recentered_velocity = function(tsne_dt,
                                    cell_a, cell_b,
                                    n_points,
                                    points_as_background = FALSE,
                                    p = NULL,
                                    min_N = 0,
                                    id_to_plot = NULL,
                                    arrow_FUN = NULL,
                                    angle_as_color = FALSE
){
    if (is.null(id_to_plot)) {
        tsne_dt.tp = copy(tsne_dt)
    } else{
        tsne_dt.tp = tsne_dt[id %in% id_to_plot]
    }
    strategy = "individual_recentered"
    av_dt = calc_delta(tsne_dt.tp, cell_a, cell_b, n_points, strategy = strategy)

    # if(is.null(av_dt$N)) av_dt$N = 1
    if (is.null(p))
        p = ggplot()
    if (points_as_background) {
        p = p + geom_point(data = tsne_dt.tp[cell %in% c(cell_a, cell_b)],
                           aes(x = tx, y = ty, color = cell))
    }
    if (is.null(av_dt$N)) {
        if (angle_as_color) {
            av_dt[, angle := xy2deg(tx_cell_a, ty_cell_a, tx_cell_b, ty_cell_b)]
            p = p + geom_segment(
                data = av_dt,
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    color = angle
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                scale_color_gradientn(
                    colours = c("orange", "red", "purple", "blue",
                                "green", "orange"),
                    limits = c(0, 360),
                    breaks = 0:4 * 90
                ) +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                theme_classic()
        } else{
            p = p + geom_segment(
                data = av_dt,
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                theme_classic()
        }
    } else{
        if (angle_as_color) {
            av_dt[, angle := xy2deg(tx_cell_a, ty_cell_a, tx_cell_b, ty_cell_b)]
            p = p + geom_segment(
                data = av_dt[N >= min_N],
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    size = N,
                    color = angle
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                scale_size_continuous(range = c(.5, 2),
                                      breaks = range(av_dt$N)) +
                scale_color_gradientn(
                    colours = c("orange", "red", "purple", "blue",
                                "green", "orange"),
                    limits = c(0, 360),
                    breaks = 0:4 * 90
                ) +
                theme_classic()
        } else{
            p = p + geom_segment(
                data = av_dt[N >= min_N],
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    size = N
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                scale_size_continuous(range = c(.5, 2),
                                      breaks = range(av_dt$N)) +
                theme_classic()
        }
    }
    p
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
                                  strategy = c("by_direction", "by_destination")[2],
                                  arrow_FUN = NULL,
                                  angle_as_color = FALSE
){
    if(is.numeric(strategy)) {
        strategy = c("by_destination",
                     "by_direction")[strategy]
    }
    if (is.null(id_to_plot)) {
        tsne_dt.tp = copy(tsne_dt)
    } else{
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
        if (angle_as_color) {
            av_dt[, angle := xy2deg(tx_cell_a, ty_cell_a, tx_cell_b, ty_cell_b)]
            p = p + geom_segment(
                data = av_dt,
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    color = angle
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                scale_color_gradientn(
                    colours = c("orange", "red", "purple", "blue",
                                "green", "orange"),
                    limits = c(0, 360),
                    breaks = 0:4 * 90
                ) +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                theme_classic()
        } else{
            p = p + geom_segment(
                data = av_dt,
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                theme_classic()
        }
    } else{
        if (angle_as_color) {
            av_dt[, angle := xy2deg(tx_cell_a, ty_cell_a, tx_cell_b, ty_cell_b)]
            p = p + geom_segment(
                data = av_dt[N >= min_N],
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    size = N,
                    color = angle
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                scale_size_continuous(range = c(.5, 2),
                                      breaks = range(av_dt$N)) +
                scale_color_gradientn(
                    colours = c("orange", "red", "purple", "blue",
                                "green", "orange"),
                    limits = c(0, 360),
                    breaks = 0:4 * 90
                ) +
                theme_classic()
        } else{
            p = p + geom_segment(
                data = av_dt[N >= min_N],
                aes(
                    x = tx_cell_a,
                    xend = tx_cell_b,
                    y = ty_cell_a,
                    yend = ty_cell_b,
                    size = N
                ),
                arrow = arrow_FUN
            ) +
                coord_cartesian(xlim = range(tsne_dt$tx),
                                ylim = range(tsne_dt$ty)) +
                scale_fill_gradient2(low = "gray", high = "black") +
                labs(x = "x",
                     y = "y",
                     fill = "density") +
                scale_size_continuous(range = c(.5, 2),
                                      breaks = range(av_dt$N)) +
                theme_classic()
        }
    }
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
plot_velocity_bins = function(velocity_dt,
                              bins = 36,
                              bin_FUN = list(length, sum, mean)[[1]]) {
    #use_distance = FALSE){
    velocity_dt[, angle_bin := floor((ang_cap(angle)) / (360 / bins))]
    # if(use_distance){
    #     b_dt = velocity_dt[, .(N = sum(distance)), angle_bin]
    # }else{
    #     b_dt = velocity_dt[, .N, angle_bin]
    # }
    if(is.null(bin_FUN)){
        bin_FUN = length
    }
    b_dt = velocity_dt[, .(bin_value = bin_FUN(distance)), angle_bin]
    ang_per_bin = 360/bins
    b_dt[, xmin := angle_bin*ang_per_bin]
    b_dt[, xmax := (angle_bin + 1)*ang_per_bin]

    # p_key = ggplot(b_dt, aes(x = angle_bin, y = bin_value, fill = angle_bin)) +
    #     geom_bar(width = 1, stat = "identity") + coord_polar() +
    #     scale_fill_gradientn(colours = c("orange", "red", "purple", "blue",
    #                                      "green", "orange"), limits = c(0, 360)/(360/bins),
    #                          breaks = 0:4*90/(360/bins), labels = function(x)x*360/bins) +
    #     scale_x_continuous(labels = function(x)x*360/bins, breaks = 0:4*90/(360/bins), limits = c(0, bins))
    #
    p_key = ggplot(b_dt, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = bin_value, fill = angle_bin)) +
        geom_rect() +
        coord_polar() +
        scale_fill_gradientn(colours = c("orange", "red", "purple", "blue",
                                         "green", "orange"), limits = c(0, 360)/(360/bins),
                             breaks = 0:4*90/(360/bins), labels = function(x)x*360/bins) +
        scale_x_continuous(breaks = 0:4*90, limits = c(0, 360))
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
                                angle.max = 360,
                                angle_as_color = TRUE,
                                background_color = "lightgray"){
    bg = velocity_dt[foreground == FALSE & distance >= delta.min]
    if(is.null(p)) p = ggplot()
    if(angle_as_color){
        p_arrows = p +
            annotate("segment",
                     x = bg$tx_cell_a, xend = bg$tx_cell_b,
                     y = bg$ty_cell_a, yend = bg$ty_cell_b,
                     color = background_color) +
            geom_segment(data = velocity_dt[foreground == TRUE],
                         aes(x = tx_cell_a, xend = tx_cell_b,
                             y = ty_cell_a, yend = ty_cell_b,
                             color = angle),
                         arrow = arrow(length = unit(0.1,"cm"))) +

            scale_color_gradientn(colours = c("orange", "red", "purple", "blue",
                                              "green", "orange"), limits = c(0, 360), breaks = 0:4*90)
    }else{
        p_arrows = p +
            annotate("segment",
                     x = bg$tx_cell_a, xend = bg$tx_cell_b,
                     y = bg$ty_cell_a, yend = bg$ty_cell_b,
                     color = background_color) +
            geom_segment(data = velocity_dt[foreground == TRUE],
                         aes(x = tx_cell_a, xend = tx_cell_b,
                             y = ty_cell_a, yend = ty_cell_b),
                         arrow = arrow(length = unit(0.1,"cm")))
    }
    p_arrows
}

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
