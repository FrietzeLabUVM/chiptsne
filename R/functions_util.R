#' wraps scales::rescale to enforce range of \code{to} on output
#'
#' @param x continuous vector of values to manipulate.
#' @param to output range (numeric vector of length two)
#' @param from input range (vector of length two). If not given, is calculated from the range of x
#'
#' @return x rescaled from \code{from} domain to \code{to} domain, within bounds of \code{to}
#' @importFrom scales rescale
#' @examples
#' #behaves identically to scales::rescale when x is within 'from' domain
#' chiptsne:::rescale_capped(0:10, to = c(0, 1), c(0, 10))
#' scales::rescale(0:10, to = c(0, 1), c(0, 10))
#' #when x exceeds 'from' domain, results are still within 'to' domain
#' chiptsne:::rescale_capped(0:10, to = c(0,1), c(0,5))
#' #not true for scales::rescale
#' scales::rescale(0:10, to = c(0,1), c(0,5))
rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

#' bin_values
#'
#' @param x Values to assign to bins
#' @param n_bins Number of bins to assign values to
#' @param xrng Optional numeric of length 2. Defines the domain to define bins
#'   for. Defaults to range(x).
#'
#' @return bin assignments parall to x.
#' @examples
#' bin_values(0:10, 3)
bin_values = function(x, n_bins, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_bins-.00001))+1
}

#' bin_values_centers
#'
#' @param n_bins Number of bins to assign values to
#' @param rng Numeric to define bins
#'   for. Defaults to range(x).
#'
#' @return numeric of length n_bins containing centers of bins.
#' @examples
#' bin_values_centers(3, 0:10)
#' bin_values_centers(3, range(0:10))
bin_values_centers = function(n_bins, rng){
    if(length(rng) != 2)
        rng = range(rng)
    stopifnot(length(rng) == 2)
    xspc = diff(rng)/n_bins/2
    xs = seq(min(rng)+xspc, max(rng)-xspc, diff(rng)/(n_bins))
    xs
}

# x: the vector n: the number of samples centered: if FALSE, then average
# current sample and previous (n-1) samples if TRUE, then average
# symmetrically in past and
# future. (If n is even, use one more sample from future.)
# from http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
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

# adds a seriees of rectangles to ggplot p
# rects is matrix, 1 rectangle per row, columns of xmin, xmax, ymin, ymax
annotate_rects = function(p,
                          rects,
                          rect_fill = NA,
                          rect_color = "black",
                          text_color = "black",
                          text_position = c("topleft", "center", "none")
){
    anns = matrix(unlist(rects), ncol = 4, byrow = TRUE)
    p_alpha = .1

    for(i in seq_len(nrow(anns))){
        p = p +
            annotate("rect",
                     xmin = anns[i, 1],
                     xmax = anns[i, 2],
                     ymin = anns[i, 3],
                     ymax = anns[i, 4],
                     color = rect_color,
                     fill = rect_fill)
        if(!is.null(text_color)){
            switch(text_position,
                   topleft = {
                       p = p +
                           annotate("label",
                                    # x = (anns[i, 1] + anns[i, 2])/2,
                                    # y = (anns[i, 3] + anns[i, 4])/2,
                                    x = anns[i, 1],
                                    y = anns[i, 4],
                                    hjust = 1,
                                    vjust = 1,
                                    label = i,
                                    color = text_color)
                   },
                   center = {
                       p = p +
                           annotate("label",
                                    # x = (anns[i, 1] + anns[i, 2])/2,
                                    # y = (anns[i, 3] + anns[i, 4])/2,
                                    x = (anns[i, 1]+anns[i, 2])/2,
                                    y = (anns[i, 3]+anns[i, 4])/2,
                                    hjust = .5,
                                    vjust = .5,
                                    label = i,
                                    color = text_color)
                   }
            )
        }
    }
    p
}


#' calc_delta
#'
#' calculates delta (distance/difference) between points in tall_var_a and tall_var_b
#' according to strategy.
#'
#' @param tsne_dt result of stsRunTsne()
#' @param tall_var_a character that matches single item in tsne_dt$tall_var.  The origin
#'   of deltas.
#' @param tall_var_b character that matches single item in tsne_dt$tall_var  The
#'   destination of deltas.
#' @param x_points number of equally spaced bins in the x-dimension. Required.
#' @param y_points number of equally spaced bins in the y-dimension. Default is
#'   x_points.
#' @param strategy character.  One of c("by_destination", "by_direction",
#'   "individual_recentered", "normal")
#'
#' @return data.table with directional and angle information added.
#'
#' @examples
#' data(tsne_dt)
#' calc_delta(tsne_dt, unique(tsne_dt$tall_var)[1], unique(tsne_dt$tall_var)[2], x_points= 4)
calc_delta = function(tsne_dt, tall_var_a, tall_var_b, x_points, y_points = x_points, strategy = "normal"){
    stopifnot(tall_var_a %in% unique(tsne_dt$tall_var))
    stopifnot(tall_var_b %in% unique(tsne_dt$tall_var))
    v_dt = dcast(tsne_dt[tall_var %in% c(tall_var_a, tall_var_b)], "id~tall_var", value.var = c("tx", "ty"))
    colnames(v_dt) = sub(tall_var_a, "tall_var_a", colnames(v_dt))
    colnames(v_dt) = sub(tall_var_b, "tall_var_b", colnames(v_dt))

    v_dt$bx_tall_var_a = bin_values(v_dt$tx_tall_var_a, n_bins = x_points, xrng = range(tsne_dt$tx))
    xs = bin_values_centers(n_bins = x_points, rng = range(tsne_dt$tx))
    v_dt$btx_tall_var_a = xs[v_dt$bx_tall_var_a]

    v_dt$by_tall_var_a = bin_values(v_dt$ty_tall_var_a, n_bins = y_points, xrng = range(tsne_dt$ty))
    ys = bin_values_centers(n_bins = y_points, rng = range(tsne_dt$ty))
    v_dt$bty_tall_var_a = ys[v_dt$by_tall_var_a]
    # strategy = "by_destination"
    # strategy = "by_direction"

    v_dt = switch(strategy,
                  "by_destination" = {
                      v_dt = v_dt[, list(tx_tall_var_b = mean(tx_tall_var_b), ty_tall_var_b = mean(ty_tall_var_b), N = .N), list(bx_tall_var_a, by_tall_var_a)]
                      v_dt$tx_tall_var_a = xs[v_dt$bx_tall_var_a]
                      v_dt$ty_tall_var_a = ys[v_dt$by_tall_var_a]
                      v_dt
                  },
                  "by_direction" = {
                      v_dt = v_dt[, list(tx_tall_var_b = mean(tx_tall_var_b - tx_tall_var_a), ty_tall_var_b = mean(ty_tall_var_b - ty_tall_var_a), N = .N), list(bx_tall_var_a, by_tall_var_a)]
                      v_dt$tx_tall_var_b = xs[v_dt$bx_tall_var_a] + v_dt$tx_tall_var_b
                      v_dt$ty_tall_var_b = ys[v_dt$by_tall_var_a] + v_dt$ty_tall_var_b
                      v_dt$tx_tall_var_a = xs[v_dt$bx_tall_var_a]
                      v_dt$ty_tall_var_a = ys[v_dt$by_tall_var_a]
                      v_dt
                  },
                  "individual_recentered" = {
                      v_dt = copy(v_dt)
                      v_dt = v_dt[, tx_tall_var_b := tx_tall_var_b - tx_tall_var_a]
                      v_dt = v_dt[, ty_tall_var_b := ty_tall_var_b - ty_tall_var_a]
                      v_dt$tx_tall_var_b = xs[v_dt$bx_tall_var_a] + v_dt$tx_tall_var_b
                      v_dt$ty_tall_var_b = ys[v_dt$by_tall_var_a] + v_dt$ty_tall_var_b
                      v_dt$tx_tall_var_a = xs[v_dt$bx_tall_var_a]
                      v_dt$ty_tall_var_a = ys[v_dt$by_tall_var_a]
                      v_dt
                  },
                  "normal" = {
                      v_dt
                  },
                  {
                      stop("unrecognized delta strategy ", strategy)
                  })


    return(v_dt[])
}

#' plot_profiles_selected
#'
#' @param profile_dt a tidy data.table containing profile data.  variable names
#' are id, tall_var, wide_var, x, y.
#' @param qtall_vars character vector of tall_vars to show data for.
#' @param id_to_plot character vector of ids to show data for.
#' @param color_mapping color_mapping for wide_vars.
#'
#' @return a ggplot
#'
#' @examples
#' data(profile_dt)
#' plot_profiles_selected(profile_dt,
#'     qtall_vars = unique(profile_dt$tall_var)[1:2],
#'     id_to_plot = unique(profile_dt$id)[1:3]) +
#'     theme(panel.spacing = unit(.04, "npc"))
plot_profiles_selected = function(profile_dt,
                                  qtall_vars,
                                  id_to_plot,
                                  color_mapping = NULL) {
    if(is.null(color_mapping)){
        color_mapping = seqsetvis::safeBrew(length(unique(profile_dt$wide_var)))
        names(color_mapping) = unique(profile_dt$wide_var)
    }
    stopifnot(qtall_vars %in% unique(profile_dt$tall_var))
    stopifnot(id_to_plot %in% profile_dt$id)
    plot_dt = profile_dt[id %in% id_to_plot & tall_var %in% qtall_vars]
    plot_dt$tall_var = factor(plot_dt$tall_var, levels = qtall_vars)
    plot_dt$id = factor(plot_dt$id, levels = id_to_plot)

    p = ggplot(plot_dt, aes(x = x, ymin = 0, ymax = y, y = y, color = wide_var, fill = wide_var)) +
        facet_grid("tall_var~id", switch = "y") +
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


#' plot_flat
#'
#' plot to a temporary png and redraw as a flat image
#'
#' useful when svg would be too expensive due to high number of items
#'
#' @param p a ggplot object to render
#' @param dpi dots per inch for temp png, default is 150
#' @param width width in inches for temp png, default it 8
#' @param height height in inches for temp png, default is 8
#'
#' @return a ggplot object containing flat image of original
#'
#' @examples
#' chiptsne:::plot_flat(ggplot())
plot_flat = function(p, dpi = 150, width = dev.size()[1], height = dev.size()[2], panel_only = FALSE){
    if(panel_only){
        pan = p +
            theme(
                title = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.margin = margin(0,0,0,0),
                legend.position="none"

            )
        ggpan = ggplot_build(pan)
        xlim = ggpan$layout$panel_scales_x[[1]]$range$range
        ylim = ggpan$layout$panel_scales_y[[1]]$range$range

        exp_f = ifelse(ggpan$layout$coord$expand, .05, 0)
        exlim = scales::expand_range(xlim, exp_f)
        eylim = scales::expand_range(ylim, exp_f)
        # pan + annotate("point", x= exlim, y = eylim, color = "red") +
        #     coord_cartesian(xlim = xlim, ylim = ylim, expand = ggpan$layout$coord$expand)
        ggsave("tmp.png", plot = pan, width = width, height = height, units = "in", dpi = dpi)
        df = data.frame(image = "tmp.png", xmin = min(exlim), ymin = min(eylim), xmax = max(exlim), ymax = max(eylim))
        pf = ggplot(df,
                    aes(image = image, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
            geom_image.rect() +
            coord_cartesian(expand = FALSE)
        pf$labels = p$labels
        # pf$scales = p$scales
        # pf2 = pf + geom_point(aes(x = xmin, y = ymin, color = xmin)) + scale_color_continuous()

    }else{
        ggsave("tmp.png", plot = p, width = width, height = height, units = "in", dpi = dpi)
        df = data.frame(image = "tmp.png", xmin = 0, ymin = 0, xmax = 1, ymax = 1)
        pf = ggplot(df,
                    aes(image = image, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
            geom_image.rect() +
            coord_cartesian(expand = FALSE) +
            theme_void()

    }
    print(pf)
}
