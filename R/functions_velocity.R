#' plot_recentered_velocity
#'
#' @param tsne_dt tidy data.table of tsne results.  variables names must
#'   include: tx, ty, id, and tall_var
#' @param tall_var_a an item in tsne_dt$tall_var. the origin for calculating
#'   angle/distance.
#' @param tall_var_b an item in tsne_dt$tall_var. the destination for calculating
#'   angle/distance.
#' @param x_points number of equally spaced bins in the x-dimension. Required.
#' @param y_points number of equally spaced bins in the y-dimension. Default is
#'   x_points.
#' @param points_as_background if TRUE, show points in background.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param id_to_plot character.  ids in tsne_dt$id. Default of NULL causes all
#'   ids to be used.
#' @param arrow_FUN result of grid::arrow().  Default of NULL does not draw
#'   arrowheads.
#' @param angle_as_color if TRUE, a rainbow like scale is applied to angle.
#'
#' @return plot where all velocities in a bin transposed to the same origin.
#'
#' @examples
#' data("tsne_dt")
#' plot_recentered_velocity(tsne_dt,
#'     unique(tsne_dt$tall_var)[1], unique(tsne_dt$tall_var)[2], 4)
#' plot_recentered_velocity(tsne_dt,
#'     unique(tsne_dt$tall_var)[1], unique(tsne_dt$tall_var)[2], 1,
#'     angle_as_color = TRUE)
plot_recentered_velocity = function(tsne_dt,
                                    tall_var_a, tall_var_b,
                                    x_points,
                                    y_points = x_points,
                                    points_as_background = FALSE,
                                    p = NULL,
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
    av_dt = calc_delta(tsne_dt.tp,
                       tall_var_a, tall_var_b,
                       x_points, y_points,
                       strategy = strategy)

    # if(is.null(av_dt$N)) av_dt$N = 1
    if (is.null(p))
        p = ggplot()
    if (points_as_background) {
        p = p + geom_point(data = tsne_dt.tp[tall_var %in% c(tall_var_a, tall_var_b)],
                           aes(x = tx, y = ty, color = tall_var))
    }
    if (angle_as_color) {
        av_dt[, angle := xy2deg(tx_tall_var_a, ty_tall_var_a, tx_tall_var_b, ty_tall_var_b)]
        p = p + geom_segment(
            data = av_dt,
            aes(
                x = tx_tall_var_a,
                xend = tx_tall_var_b,
                y = ty_tall_var_a,
                yend = ty_tall_var_b,
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
                x = tx_tall_var_a,
                xend = tx_tall_var_b,
                y = ty_tall_var_a,
                yend = ty_tall_var_b
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
    p
}

#' plot_regional_velocity
#'
#' @param tsne_dt tidy data.table of tsne results.  variables names must
#'   include: tx, ty, id, and tall_var
#' @param tall_var_a an item in tsne_dt$tall_var. the origin for calculating
#'   angle/distance.
#' @param tall_var_b an item in tsne_dt$tall_var. the destination for calculating
#'   angle/distance.
#' @param x_points number of equally spaced bins in the x-dimension. Required.
#' @param y_points number of equally spaced bins in the y-dimension. Default is
#'   x_points.
#' @param points_as_background if TRUE, show points in background.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param min_N regional bins with fewer than min_N counts will not be shown.
#' @param id_to_plot character.  ids in tsne_dt$id. Default of NULL causes all
#'   ids to be used.
#' @param strategy character.  must be one of "by_direction" or "by_destination"
#' @param arrow_FUN result of grid::arrow().  Default of NULL does not draw
#'   arrowheads.
#' @param angle_as_color if TRUE, a rainbow like scale is applied to angle.
#'
#' @return ggplot where regionally binned velocities are summarized.
#'
#' @examples
#' data("tsne_dt")
#' plot_regional_velocity(tsne_dt, unique(tsne_dt$tall_var)[1],
#'     unique(tsne_dt$tall_var)[2], 4, min_N = 4, arrow_FUN = arrow())
plot_regional_velocity = function(tsne_dt,
                                  tall_var_a,
                                  tall_var_b,
                                  x_points,
                                  y_points = x_points,
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

    av_dt = calc_delta(tsne_dt.tp, tall_var_a, tall_var_b, x_points, y_points, strategy = strategy)

    # if(is.null(av_dt$N)) av_dt$N = 1
    if (is.null(p))
        p = ggplot()
    if (points_as_background) {
        p = p + geom_point(data = tsne_dt.tp[tall_var %in% c(tall_var_a, tall_var_b)],
                           aes(x = tx, y = ty, color = tall_var))
    }
    if (is.null(av_dt$N)) {
        if (angle_as_color) {
            av_dt[, angle := xy2deg(tx_tall_var_a, ty_tall_var_a, tx_tall_var_b, ty_tall_var_b)]
            p = p + geom_segment(
                data = av_dt,
                aes(
                    x = tx_tall_var_a,
                    xend = tx_tall_var_b,
                    y = ty_tall_var_a,
                    yend = ty_tall_var_b,
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
                    x = tx_tall_var_a,
                    xend = tx_tall_var_b,
                    y = ty_tall_var_a,
                    yend = ty_tall_var_b
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
            av_dt[, angle := xy2deg(tx_tall_var_a, ty_tall_var_a, tx_tall_var_b, ty_tall_var_b)]
            p = p + geom_segment(
                data = av_dt[N >= min_N],
                aes(
                    x = tx_tall_var_a,
                    xend = tx_tall_var_b,
                    y = ty_tall_var_a,
                    yend = ty_tall_var_b,
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
                    x = tx_tall_var_a,
                    xend = tx_tall_var_b,
                    y = ty_tall_var_a,
                    yend = ty_tall_var_b,
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

#' plot_velocity_centered
#'
#' @param velocity_dt data.table with angle and distance info.
#' @param line_alpha numeric [0,1] alpha value of lines
#' @param line_size numeric >0 size value of lines
#'
#' @return plot of
#'
#' @examples
#' data("tsne_dt")
#' vel_dt = prep_velocity(tsne_dt, unique(tsne_dt$tall_var)[1], unique(tsne_dt$tall_var)[2])
#' plot_velocity_centered(vel_dt)
plot_velocity_centered = function(velocity_dt, line_alpha = .5, line_size = 1){
    ggplot(velocity_dt, aes(x = angle, xend = angle, y = 0, yend = distance, color = angle)) +
        geom_segment(stat = "identity", alpha = line_alpha, size = line_size) +
        coord_polar() +
        scale_color_gradientn(colours = c("orange", "red", "purple", "blue",
                                          "green", "orange"), limits = c(0, 360), breaks = 0:4*90) +
        scale_x_continuous(limits = c(0, 360), breaks = 0:4*90)
}

#' Title
#'
#' @param velocity_dt data.table with angle and distance info.
#' @param bins number of bins to use when summarizing changes.
#' @param bin_FUN function to apply per bin to summarize. will be
#' passed a numeric vector and must return single value.  Suggested methods
#' are length, sum, mean, and median.
#'
#' @return a ggplot with a polar coordinate barplot.
#'
#' @examples
#' data("tsne_dt")
#' vel_dt = prep_velocity(
#'     tsne_dt,
#'     unique(tsne_dt$tall_var)[1],
#'     unique(tsne_dt$tall_var)[2]
#'     )
#' plot_velocity_bins(vel_dt, bins = 18, bin_FUN = sum)
#' plot_velocity_bins(vel_dt, bins = 8, bin_FUN = median)
plot_velocity_bins = function(velocity_dt,
                              bins = 36,
                              bin_FUN = sum) {
    bin_value = NULL #data.table bindin
    velocity_dt[, angle_bin := floor((ang_cap(angle)) / (360 / bins))]
    if(is.null(bin_FUN)){
        bin_FUN = length
    }
    b_dt = velocity_dt[, .(bin_value = bin_FUN(distance)), angle_bin]
    ang_per_bin = 360/bins
    b_dt[, xmin := angle_bin*ang_per_bin]
    b_dt[, xmax := (angle_bin + 1)*ang_per_bin]

    p_key = ggplot(b_dt,
                   aes(
                       xmin = xmin,
                       xmax = xmax,
                       ymin = 0,
                       ymax = bin_value,
                       fill = angle_bin
                   )) +
        geom_rect() +
        coord_polar() +
        scale_fill_gradientn(
            colours = c("orange", "red", "purple", "blue",
                        "green", "orange"),
            limits = c(0, 360) / (360 / bins),
            breaks = 0:4 * 90 / (360 / bins),
            labels = function(x)
                x * 360 / bins
        ) +
        scale_x_continuous(breaks = 0:4 * 90, limits = c(0, 360))
    p_key
}

#' plot_velocity_arrows
#'
#' @param velocity_dt data.table with angle and distance info.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param id_to_plot character.  ids in tsne_dt$id. Default of NULL causes all
#'   ids to be used.
#' @param angle_as_color if TRUE, a rainbow like scale is applied to angle.
#' @param background_color character. color to use when plotting background
#'   points.
#'
#' @return ggplot of arrows indicating position in tall_var line a to tall_var line b.
#'
#' @examples
#' data("tsne_dt")
#' vel_dt = prep_velocity(tsne_dt,
#'     unique(tsne_dt$tall_var)[2], unique(tsne_dt$tall_var)[3])
#' plot_velocity_arrows(vel_dt)
plot_velocity_arrows = function(velocity_dt,
                                p = NULL,
                                id_to_plot = NULL,
                                angle_as_color = TRUE,
                                background_color = "lightgray"){
    # bg = velocity_dt[foreground == FALSE & distance >= delta.min]
    bg = velocity_dt[foreground == FALSE]
    if(is.null(p)) p = ggplot()
    if(angle_as_color){
        p_arrows = p +
            annotate("segment",
                     x = bg$tx_tall_var_a, xend = bg$tx_tall_var_b,
                     y = bg$ty_tall_var_a, yend = bg$ty_tall_var_b,
                     color = background_color) +
            geom_segment(data = velocity_dt[foreground == TRUE],
                         aes(x = tx_tall_var_a, xend = tx_tall_var_b,
                             y = ty_tall_var_a, yend = ty_tall_var_b,
                             color = angle),
                         arrow = arrow(length = unit(0.1,"cm"))) +

            scale_color_gradientn(colours = c("orange", "red", "purple", "blue",
                                              "green", "orange"), limits = c(0, 360), breaks = 0:4*90)
    }else{
        p_arrows = p +
            annotate("segment",
                     x = bg$tx_tall_var_a, xend = bg$tx_tall_var_b,
                     y = bg$ty_tall_var_a, yend = bg$ty_tall_var_b,
                     color = background_color) +
            geom_segment(data = velocity_dt[foreground == TRUE],
                         aes(x = tx_tall_var_a, xend = tx_tall_var_b,
                             y = ty_tall_var_a, yend = ty_tall_var_b),
                         arrow = arrow(length = unit(0.1,"cm")))
    }
    p_arrows
}

#' prep_velocity
#'
#' @param tsne_dt tidy data.table of tsne results.  variables names must
#'   include: tx, ty, id, and tall_var
#' @param tall_var_a an item in tsne_dt$tall_var. the origin for calculating
#'   angle/distance.
#' @param tall_var_b an item in tsne_dt$tall_var. the destination for calculating
#'   angle/distance.
#' @param id_to_plot character.  ids in tsne_dt$id. Default of NULL causes all
#'   ids to be used.
#' @param max_plotted the maximum number of ids to plot.
#' @param delta.min numeric. ids below this delta marked as background.
#' @param delta.max numeric. ids above this delta marked as background.
#' @param angle.min numeric. ids with angle below this marked as background.
#' @param angle.max numeric. ids with angle above this marked as background.
#' @param drop_backgroud if TRUE, all background velocities are dropped. Default
#'   is FALSE.
#'
#' @return a tidy data.table of velocity measurements from tall_var_a to tall_var_b
#'
#' @examples
#' data("tsne_dt")
#' vel_dt = prep_velocity(tsne_dt,
#'     unique(tsne_dt$tall_var)[2], unique(tsne_dt$tall_var)[3])
#' vel_dt
prep_velocity = function(tsne_dt,
                         tall_var_a,
                         tall_var_b,
                         id_to_plot = NULL,
                         max_plotted = 500,
                         delta.min = 0,
                         delta.max = Inf,
                         angle.min = 0,
                         angle.max = 360,
                         drop_backgroud = FALSE) {

    v_dt = calc_delta(tsne_dt, tall_var_a, tall_var_b, 10)

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
    v_dt.tp[, angle := xy2deg(x1 = tx_tall_var_a, x2 = tx_tall_var_b, y1 = ty_tall_var_a, y2 = ty_tall_var_b)]
    if(angle.min > angle.max){
        v_dt.tp[, foreground := angle <= angle.min & angle >= angle.max]
    }else{
        v_dt.tp[, foreground := angle >= angle.min & angle <= angle.max]
    }
    # distance selection
    v_dt.tp[, distance := xy2dist(x1 = tx_tall_var_a, x2 = tx_tall_var_b, y1 = ty_tall_var_a, y2 = ty_tall_var_b)]
    v_dt.tp = v_dt.tp[distance < delta.min | distance > delta.max, foreground := FALSE]
    if(drop_backgroud){
        v_dt.tp = v_dt.tp[foreground == TRUE]
    }
    v_dt.tp[]
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
