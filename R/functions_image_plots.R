#' stsPlotSummaryProfiles
#'
#' @param profile_dt a tidy data.table for profile data as retrieved by
#'   stsFetchTsneInput.  Expected variable names are id, tall_var, wide_var, x, and y.
#' @param position_dt a tidy data.table containing t-sne embedding.  Expected
#'   variable names are tx, ty, id, and tall_var.
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#'   Defaults to same value as x_points.
#' @param q_tall_vars character vector of tall_vars to plot. Default of NULL plots all.
#' @param q_wide_vars character vector of wide_vars to plot.  Default of NULL plots all.
#' @param xrng view domain in x dimension, default is range of position_dt$tx.
#' @param yrng view domain in y dimension, default is range of position_dt$ty.
#' @param plot_type character, must be one of "glyph" or "raster".  raster
#' uses ggimage::geom_image to embed images where glyph uses GGally::glyphs
#' @param rname prefix for image files.  existing image files are used if
#'   present.
#' @param odir output directory for image files.
#' @param force_rewrite if TRUE, images are overwritten even if they exist.
#' @param n_cores number of cores to use when writing images.  Default is value
#'   of mc.cores option if set or 1.
#' @param apply_norm if TRUE, y values are trimmed to 95th percentile and
#'   transformed ot domain of [0,1]. Default is TRUE.
#' @param ylim y-limits of regional summary plots.  Default of c(0, 1) is
#'   compatible with apply_norm = TRUE.
#' @param ma_size moving average size when smoothing profiles.
#' @param n_splines number of points to interpolate with splines.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param facet_byCell boolean. If TRUE, plots are facetted by tall_var.
#' @param line_color_mapping named vector of line color.  Names correspond to
#'   values of profile_dt 'wide_var' variable and values are colors.
#' @param vertical_facet_mapping named vector of groups.  vertical_facet_mapping
#'   names are wide_vars and order of first occurence determines vertical position
#'   top to bottom.
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#' @param return_data if TRUE, data.table that would have been used to create
#'   ggplot is returned instead.
#'
#' @return returns a ggplot containing images summarizing t-sne space at
#'   resolution determined by x_points and y_points with the space density
#'   mapped to size.
#' @export
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' stsPlotSummaryProfiles(profile_dt, tsne_dt, x_points = 4)
stsPlotSummaryProfiles = function(## basic inputs
    profile_dt,
    position_dt,
    ## resolution
    x_points,
    y_points = x_points,
    ## selection
    q_tall_vars = NULL,
    q_wide_vars = NULL,
    ## view domain
    xrng = range(position_dt$tx),
    yrng = range(position_dt$ty),
    ## plotting strategy
    plot_type = c("glyph", "raster")[1],
    ## image file path and output
    rname = NULL,
    odir = file.path(tempdir(), rname),
    force_rewrite = FALSE,
    n_cores = getOption("mc.cores", 1),
    ## preplot transformation
    apply_norm = TRUE,
    ylim = c(0, 1),
    ma_size = 2,
    n_splines = 10,
    ## plot organization
    p = NULL,
    facet_byCell = FALSE,
    line_color_mapping = NULL,
    vertical_facet_mapping = NULL,
    ## sizing and level of detail
    N_floor = 0,
    N_ceiling = NULL,
    min_size = .3,
    ## output
    return_data = FALSE) {
    if(is.null(q_tall_vars)){
        if(is.factor(profile_dt$tall_var)){
            q_tall_vars = levels(profile_dt$tall_var)
        }else{
            q_tall_vars = sort(unique(profile_dt$tall_var))
        }
    }
    if(is.null(q_wide_vars)){
        if(is.factor(profile_dt$wide_var)){
            q_wide_vars = levels(profile_dt$wide_var)
        }else{
            q_wide_vars = sort(unique(profile_dt$wide_var))
        }
    }
    if(is.null(line_color_mapping)){
        line_color_mapping = seqsetvis::safeBrew(length(unique(profile_dt$wide_var)))
        names(line_color_mapping) = unique(profile_dt$wide_var)
    }
    line_color_mapping = line_color_mapping[names(line_color_mapping) %in% q_wide_vars]
    if(is.factor(profile_dt$tall_var)){
        stopifnot(q_tall_vars %in% levels(profile_dt$tall_var))
    }else{
        stopifnot(q_tall_vars %in% unique(profile_dt$tall_var))
    }
    if(is.factor(profile_dt$wide_var)){
        stopifnot(q_wide_vars %in% levels(profile_dt$wide_var))
    }else{
        stopifnot(q_wide_vars %in% unique(profile_dt$wide_var))
    }
    if(!plot_type %in% c("glyph", "raster")){
        stop("plot_type (\"", plot_type, "\") must be one of \"glyph\" or \"raster\".")
    }

    prof_dt = copy(profile_dt[tall_var %in% q_tall_vars & wide_var %in% q_wide_vars])
    pos_dt = copy(position_dt[tall_var %in% q_tall_vars])

    if(is.null(rname)){
        rname = digest::digest(
            list(
                prof_dt,
                pos_dt,
                x_points,
                y_points,
                xrng,
                yrng,
                apply_norm,
                ylim,
                line_color_mapping,
                n_splines,
                ma_size,
                facet_byCell
            )
        )
    }

    prof_dt$tall_var = factor(prof_dt$tall_var, levels = q_tall_vars)
    prof_dt$wide_var = factor(prof_dt$wide_var, levels = q_wide_vars)
    pos_dt$tall_var = factor(pos_dt$tall_var, levels = q_tall_vars)
    if(plot_type == "raster"){
        #prepare images
        if (!facet_byCell) {
            summary_dt = prep_summary(prof_dt,
                                      pos_dt,
                                      x_points,
                                      y_points,
                                      xrng,
                                      yrng,
                                      NULL)
            img_res = prep_images(
                summary_dt = summary_dt,
                xrng = xrng,
                yrng = yrng,
                x_points = x_points,
                y_points = y_points,
                rname = rname,
                odir = odir,
                force_rewrite = force_rewrite,
                apply_norm = apply_norm,
                ylim = ylim,
                ma_size = ma_size,
                n_splines = n_splines,
                n_cores = n_cores,
                line_color_mapping = line_color_mapping,
                vertical_facet_mapping = vertical_facet_mapping
            )
            plot_summary_raster(
                image_dt = img_res$image_dt,
                xrng = xrng,
                yrng = yrng,
                x_points = x_points,
                y_points = y_points,
                p = p,
                line_color_mapping = img_res$line_color_mapping,
                N_floor = N_floor,
                N_ceiling = N_ceiling,
                min_size = min_size,
                return_data = return_data
            )
        } else{#do facet by tall_var
            summary_dt = prep_summary(prof_dt,
                                      pos_dt,
                                      x_points,
                                      y_points,
                                      xrng,
                                      yrng,
                                      "tall_var")
            img_res = prep_images(
                summary_dt = summary_dt,
                xrng = xrng,
                yrng = yrng,
                x_points = x_points,
                y_points = y_points,
                rname = rname,
                odir = odir,
                force_rewrite = force_rewrite,
                apply_norm = apply_norm,
                ylim = ylim,
                facet_by = "tall_var",
                ma_size = ma_size,
                n_splines = n_splines,
                n_cores = n_cores,
                line_color_mapping = line_color_mapping,
                vertical_facet_mapping = vertical_facet_mapping
            )
            plot_summary_raster_byCell(
                image_dt = img_res$image_dt,
                xrng = xrng,
                yrng = yrng,
                x_points = x_points,
                y_points = y_points,
                p = p,
                line_color_mapping = img_res$line_color_mapping,
                N_floor = N_floor,
                N_ceiling = N_ceiling,
                min_size = min_size,
                return_data = return_data
            )
        }
    }else if(plot_type == "glyph"){
        if (!facet_byCell) {
            summary_dt = prep_summary(prof_dt,
                                      pos_dt,
                                      x_points,
                                      y_points,
                                      xrng,
                                      yrng,
                                      NULL)
            plot_summary_glyph(
                summary_dt,
                x_points = x_points,
                y_points = y_points,
                xrng = xrng,
                yrng = yrng,
                ylim = ylim,
                N_floor = N_floor,
                N_ceiling = N_ceiling,
                min_size = min_size,
                color_mapping = line_color_mapping,
                return_data = return_data
            )
        }else{
            summary_dt_l = lapply(q_tall_vars, function(cl){
                prep_summary(prof_dt[tall_var == cl],
                             pos_dt,
                             x_points,
                             y_points,
                             xrng,
                             yrng,
                             NULL)
            })
            names(summary_dt_l) = q_tall_vars
            summary_dt = rbindlist(summary_dt_l, use.names = TRUE, idcol = "tall_var")
            if(return_data){
                plot_summary_glyph(
                    summary_dt,
                    x_points = x_points,
                    y_points = y_points,
                    xrng = xrng,
                    yrng = yrng,
                    ylim = ylim,
                    N_floor = N_floor,
                    N_ceiling = N_ceiling,
                    min_size = min_size,
                    color_mapping = line_color_mapping,
                    return_data = return_data
                )
            }else{
                plot_summary_glyph(
                    summary_dt,
                    x_points = x_points,
                    y_points = y_points,
                    xrng = xrng,
                    yrng = yrng,
                    ylim = ylim,
                    N_floor = N_floor,
                    N_ceiling = N_ceiling,
                    min_size = min_size,
                    color_mapping = line_color_mapping
                ) +
                    facet_wrap("tall_var")
            }
        }

    }

}

#' prep_summary
#'
#' @param profile_dt a tidy data.table for profile data as retrieved by
#'   stsFetchTsneInput.  Expected variable names are id, tall_var, wide_var, x, and y.
#' @param position_dt a tidy data.table containing t-sne embedding.  Expected
#'   variable names are tx, ty, id, and tall_var.
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#'   Defaults to same value as x_points.
#' @param xrng view domain in x dimension, default is range of position_dt$tx.
#' @param yrng view domain in y dimension, default is range of position_dt$ty.
#' @param facet_by character. variable name to facet profile_dt by when
#'   constructing images. The only valid non-null value with chiptsne functions
#'   is "tall_var".
#'
#' @return summary of profiles binned across tsne space according to x_points,
#'   y_points, and within xrng and yrng
#' @export
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' summary_dt = prep_summary(profile_dt, tsne_dt, 4)
#' img_res = prep_images(summary_dt, 4)
#' #zoom on top-right quadrant
#' summary_dt.zoom = prep_summary(profile_dt, tsne_dt, 4,
#'     xrng = c(0, .5), yrng = c(0, .5))
#' img_res.zoom = prep_images(summary_dt.zoom, 4,
#'     xrng = c(0, .5), yrng = c(0, .5))
prep_summary = function(profile_dt,
                        position_dt,
                        x_points,
                        y_points = x_points,
                        xrng = range(position_dt$tx),
                        yrng = range(position_dt$ty),
                        facet_by = NULL){
    position_dt = copy(position_dt[tx >= min(xrng) &
                                       tx <= max(xrng) &
                                       ty >= min(yrng) &
                                       ty <= max(yrng)])
    #use positional info from position_dt to bin points
    position_dt = position_dt[id %in% unique(profile_dt$id)]
    if(is.null(position_dt$bx))
        position_dt[, bx := bin_values(tx, x_points, xrng = xrng)]
    if(is.null(position_dt$by))
        position_dt[, by := bin_values(ty, y_points, xrng = yrng)]
    #merge binning info to profiles
    summary_dt = merge(
        profile_dt,
        position_dt[, list(bx, by, tall_var, id)],
        allow.cartesian = TRUE,
        by = intersect(colnames(profile_dt), c("tall_var", "id"))
    )#, by = c("tall_var", "id"))
    if (is.null(summary_dt$wide_var))
        summary_dt$wide_var = "signal"

    if (is.null(facet_by)) {
        summary_dt = summary_dt[, list(y = mean(y)), list(bx, by, x, wide_var)]
    } else{
        summary_dt = summary_dt[, list(y = mean(y)), list(bx, by, x, wide_var, get(facet_by))]
        colnames(summary_dt)[colnames(summary_dt) == "get"] = facet_by
    }
    N_dt = position_dt[, .(.N), by = .(bx, by)]
    summary_dt = merge(summary_dt, N_dt, by = c("bx", "by"))
    #each combination of bx and by is a unique plot_id
    summary_dt[, plot_id := paste(bx, by, sep = "_")]
    summary_dt[]
}

#' prep_images
#'
#' Prepares images summarizing profiles in t-sne regions and returns ggimage
#' compatible data.frame for plotting.
#'
#' @param summary_dt a tidy data.table from prep_summary()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#'   Defaults to same value as x_points.
#' @param xrng view domain in x dimension, default is range of position_dt$tx.
#' @param yrng view domain in y dimension, default is range of position_dt$ty.
#' @param rname prefix for image files.  existing image files are used if
#'   present.
#' @param odir output directory for image files.
#' @param force_rewrite if TRUE, images are overwritten even if they exist.
#' @param apply_norm if TRUE, y values are trimmed to 95th percentile and
#'   transformed ot domain of [0,1]. Default is TRUE.
#' @param ylim y-limits of regional summary plots.  Default of c(0, 1) is
#'   compatible with apply_norm = TRUE.
#' @param facet_by character. varaible name to facet profile_dt by when
#'   constructing images. The only valid non-null value with chiptsne functions
#'   is "tall_var".
#' @param n_splines number of points to interpolate with splines.
#' @param ma_size moving average size when smoothing profiles.
#' @param n_cores number of cores to use when writing images.  Default is value
#'   of mc.cores option if set or 1.
#' @param line_color_mapping named vector of line color.  Names correspond to
#'   values of profile_dt 'wide_var' variable and values are colors.
#' @param vertical_facet_mapping named vector of vertical facet for data
#'
#' @return data.table with variables
#' @export
#' @importFrom seqsetvis applySpline
#' @importFrom stats quantile
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' summary_dt = prep_summary(profile_dt, tsne_dt, 4)
#' img_res = prep_images(summary_dt, 4)
#' #zoom on top-right quadrant
#' summary_dt.zoom = prep_summary(profile_dt, tsne_dt, 4,
#'     xrng = c(0, .5), yrng = c(0, .5))
#' img_res.zoom = prep_images(summary_dt.zoom, 4,
#'     xrng = c(0, .5), yrng = c(0, .5))
#' #use results with plot_summary_raster() to make plots
prep_images = function(summary_dt,
                       x_points,
                       y_points = x_points,
                       xrng = c(-.5, .5),
                       yrng = c(-.5, .5),
                       rname = NULL,
                       odir = file.path(tempdir(), rname),
                       force_rewrite = FALSE,
                       apply_norm = TRUE,
                       ylim = c(0, 1),
                       facet_by = NULL,
                       # view_rect = list(),
                       ma_size = 2,
                       n_splines = 10,
                       n_cores = getOption("mc.cores", 1),
                       line_color_mapping = NULL,
                       vertical_facet_mapping = NULL) {
    if(is.null(line_color_mapping)){
        line_color_mapping = seqsetvis::safeBrew(length(unique(summary_dt$wide_var)))
        names(line_color_mapping) = unique(summary_dt$wide_var)
    }
    if (!all(unique(summary_dt$wide_var) %in% names(line_color_mapping))) {
        missing_colors = setdiff(unique(summary_dt$wide_var), names(line_color_mapping))
        stop(
            "line_color_mapping is missing assignments for: ",
            paste(missing_colors, collapse = ", ")
        )
    }
    if(is.null(rname)){
        rname = digest::digest(
            list(
                summary_dt,
                x_points,
                y_points,
                apply_norm,
                ylim,
                line_color_mapping,
                n_splines,
                ma_size,
                facet_by
            )
        )
    }

    # position_dt = copy(position_dt[tx >= min(xrng) & tx <= max(xrng) &
    #                                    ty >= min(yrng) & ty <= max(yrng)])
    # position_dt[, bx := bin_values(tx, x_points, xrng = xrng)]
    # position_dt[, by := bin_values(ty, y_points, xrng = yrng)]

    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    if (is.null(facet_by)) {
        img_dt = unique(summary_dt[, list(bx, by, plot_id)])
        img_dt[, png_file := file.path(odir, paste0(plot_id, ".png"))]
    } else{
        img_dt = unique(summary_dt[, list(bx, by, plot_id, get(facet_by))])
        colnames(img_dt)[4] = facet_by
        img_dt[, png_file := file.path(odir, paste0(get(facet_by), "_", plot_id, ".png"))]
    }

    xs = bin_values_centers(x_points, rng = xrng)
    ys = bin_values_centers(y_points, rng = yrng)

    # plot(expand.grid(xs, ys), xlim = xrng, ylim = yrng)
    # rect(min(xrng), min(yrng), max(xrng), max(yrng), col = rgb(0,0,1,.1))

    img_dt[, tx := xs[bx]]
    img_dt[, ty := ys[by]]


    if (apply_norm) {
        summary_dt[, ynorm := y / stats::quantile(y, .95), by = list(wide_var)]
        summary_dt[ynorm > 1, ynorm := 1]
    } else{
        summary_dt[, ynorm := y]
    }


    if (force_rewrite) {
        file.remove(img_dt$png_file[file.exists(img_dt$png_file)])
    }

    if (is.null(facet_by)) {
        count_dt = unique(summary_dt[, .(bx, by, N)])
        img_dt = merge(img_dt, count_dt, by = c("bx", "by"))
    } else{
        count_dt = unique(summary_dt[, .(bx, by, get(facet_by), N)])
        colnames(count_dt)[colnames(count_dt) == "V3"] = facet_by
        img_dt = merge(img_dt, count_dt, by = c("bx", "by", facet_by))
    }

    if (is.null(vertical_facet_mapping)) {
        summary_dt$group = 1
    } else{
        if (!all(unique(summary_dt$wide_var) %in% names(vertical_facet_mapping))) {
            missing_groups = setdiff(unique(summary_dt$wide_var), names(vertical_facet_mapping))
            stop(
                "vertical_facet_mapping is missing assignments for: ",
                paste(missing_groups, collapse = ", ")
            )
        }
        summary_dt$group = factor(vertical_facet_mapping[summary_dt$wide_var], levels = unique(vertical_facet_mapping))
    }


    if (any(!file.exists(img_dt$png_file))) {
        plot_info = lapply(which(!file.exists(img_dt$png_file)), function(i) {
            # hidden = parallel::mclapply(which(!file.exists(img_dt$png_file)), function(i){}, mc.cores = n_cores)
            fpath = img_dt$png_file[i]
            p_id = img_dt$plot_id[i]
            if (is.null(facet_by)) {
                pdt = summary_dt[plot_id == p_id]
            } else{
                pdt = summary_dt[plot_id == p_id & get(facet_by) == img_dt[[facet_by]][i]]
            }

            # pdt[, ysm := seqsetvis:::movingAverage(ynorm, n = ma_size), by = list(wide_var)]
            pdt[, ysm := movingAverage(ynorm, n = ma_size), by = list(wide_var)]

            pdt = seqsetvis::applySpline(pdt,
                                         n = n_splines,
                                         by_ = "wide_var",
                                         y_ = "ysm")
            list(pdt, fpath)
        })
        # hidden = lapply(plot_info, function(x){
        # figure out how not to copy global env
        hidden = parallel::mclapply(plot_info, function(x) {
            fpath = x[[2]]
            pdt = x[[1]]
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = list(wide_var)]

            p = ggplot(pdt,
                       aes(
                           x = x,
                           y = ysm,
                           ymin = 0,
                           ymax = ysm,
                           color = wide_var,
                           fill = wide_var
                       )) +
                geom_ribbon(alpha = .3) +
                geom_path(size = .6, alpha = 1) +
                scale_color_manual(values = line_color_mapping) +
                scale_fill_manual(values = line_color_mapping) +
                theme_void() + guides(color = "none", fill =
                                          'none') +
                facet_wrap("group", ncol = 1) +
                theme(strip.background = element_blank(),
                      strip.text = element_blank()) +
                coord_cartesian(
                    ylim = ylim,
                    xlim = c(-.5, .5),
                    expand = FALSE
                )
            ggsave(
                fpath,
                p,
                width = 2,
                height = 2,
                units = "cm"
            )
            # p
            NULL
            # })
        }, mc.cores = n_cores)
    }

    img_dt[, png_file := normalizePath(png_file)]

    # if (is.factor(position_dt$tall_var)) {
    #     if (!is.null(img_dt$tall_var)) {
    #         img_dt$tall_var = factor(img_dt$tall_var, levels = levels(position_dt$tall_var))
    #     }
    #     if (!is.null(summary_dt$tall_var)) {
    #         summary_dt$tall_var = factor(summary_dt$tall_var, levels = levels(position_dt$tall_var))
    #     }
    # }

    return(
        list(
            image_dt = img_dt[],
            summary_profile_dt = summary_dt[],
            x_points = x_points,
            y_points = y_points,
            xrng = xrng,
            yrng = yrng,
            line_color_mapping = line_color_mapping
        )
    )
}


set_size = function(dt, N_floor, N_ceiling, size.name = "img_size"){
    tmp_var = NULL #binding for data.table
    stopifnot("N" %in% colnames(dt))
    if (is.null(N_ceiling)) {
        N_ceiling = max(dt$N)
    }
    dt[, tmp_var := N]
    dt[tmp_var > N_ceiling, tmp_var := N_ceiling]
    dt[tmp_var < N_floor, tmp_var := N_floor]

    dt[, tmp_var := tmp_var - N_floor]
    dt[, tmp_var := tmp_var / N_ceiling]
    dt[[size.name]] = dt$tmp_var
    dt$tmp_var = NULL
    dt
}

#' set_image_rects
#'
#' configures result of prep_images by setting rectangle parameters
#'
#' @param image_dt $image_dt of result from prep_images()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#'
#' @return image_dt with rect aesthetics (xmin, xmax, ymin, and ymax) added.
#' @export
#'
#' @examples
#' library(ggplot2)
#' data("profile_dt")
#' data("tsne_dt")
#' summary_dt = prep_summary(profile_dt,
#'                       tsne_dt,
#'                       x_points = 4)
#' img_res = prep_images(summary_dt, 4)
#' img_rect = set_image_rects(img_res$image_dt,
#'                            x_points = img_res$x_points,
#'                            y_points = img_res$y_points,
#'                            xrng = img_res$xrng,
#'                            yrng = img_res$yrng)
#' ggplot(img_rect, aes(xmin = xmin, xmax = xmax,
#'     ymin = ymin, ymax = ymax)) + geom_rect()
#' ggplot(img_rect, aes(xmin = xmin, xmax = xmax,
#'     ymin = ymin, ymax = ymax, image = png_file)) + geom_image.rect()
set_image_rects = function(image_dt,
                           x_points,
                           y_points,
                           xrng,
                           yrng,
                           N_floor = 0,
                           N_ceiling = NULL,
                           min_size = .3) {
    image_dt = set_size(image_dt, N_floor, N_ceiling, size.name = "img_size")
    # image_dt[, img_size := N]
    # image_dt[img_size > N_ceiling, img_size := N_ceiling]
    # image_dt[img_size < N_floor, img_size := N_floor]
    #
    # image_dt[, img_size := img_size - N_floor]
    # image_dt[, img_size := img_size / N_ceiling]
    image_dt = image_dt[img_size >= min_size]

    xspc = diff(xrng) / x_points / 2
    yspc = diff(yrng) / y_points / 2

    image_dt[, xmin := tx - xspc * img_size]
    image_dt[, xmax := tx + xspc * img_size]
    image_dt[, ymin := ty - yspc * img_size]
    image_dt[, ymax := ty + yspc * img_size]
    image_dt[]
}

#' plot_summary_raster
#'
#' @param image_dt $image_dt of result from prep_images()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param line_color_mapping named character vector for scale_color_manual, recommend $line_color_mapping from prep_images()
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#' @param return_data if TRUE, data.table that would have been used to create
#'   ggplot is returned instead.
#'
#' @return ggplot containing images summarizing t-sne regions.  if return_data
#'   is TRUE, instead the data.table containing plot info is returned.
#' @export
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' summary_dt = prep_summary(profile_dt, tsne_dt, 4)
#' img_res = prep_images(summary_dt, 4)
#' #zoom on top-right quadrant
#'
#' summary_dt.zoom = prep_summary(profile_dt, tsne_dt, 8, xrng = c(0, .5), yrng = c(0, .5))
#' img_res.zoom = prep_images(summary_dt.zoom, 8, xrng = c(0, .5), yrng = c(0, .5))
#' plot_summary_raster(img_res$image_dt,
#'               x_points = img_res$x_points)
#' plot_summary_raster(img_res.zoom$image_dt,
#'               x_points = img_res.zoom$x_points,
#'               xrng = img_res.zoom$xrng,
#'               yrng = img_res.zoom$yrng, min_size = 0)
plot_summary_raster = function(image_dt,
                               x_points,
                               y_points = x_points,
                               xrng = c(-.5, .5),
                               yrng = c(-.5, .5),
                               p = NULL,
                               line_color_mapping = NULL,
                               N_floor = 0,
                               N_ceiling = NULL,
                               min_size = .3,
                               return_data = FALSE) {
    image_dt = copy(image_dt)
    image_dt = set_image_rects(
        image_dt,
        x_points = x_points,
        y_points = y_points,
        xrng = xrng,
        yrng = yrng,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_size
    )
    if (return_data) {
        return(image_dt)
    }
    if (is.null(p))
        p = ggplot()
    if(!is.null(line_color_mapping)){
        col_dt = image_dt[, .(tx, ty, wide_var = rep(names(line_color_mapping), each = nrow(image_dt)))]
        p = p + geom_point(data = col_dt,
                           aes(x = tx, y = ty, color = wide_var)) +
            scale_color_manual(values = line_color_mapping)
    }
    p = p +
        geom_image.rect(data = image_dt,
                        aes(
                            xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax,
                            image = png_file
                        )) + #, color = rgb(0,0,1,.2)) +
        geom_rect(
            data = image_dt,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax
            ),
            fill = NA,
            color = "black"
        ) +
        coord_cartesian(xlim = xrng, ylim = yrng)
    p
}

#' plot_summary_raster
#'
#' @param image_dt $image_dt of result from prep_images()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param line_color_mapping named character vector for scale_color_manual, recommend $line_color_mapping from prep_images()
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#' @param return_data if TRUE, data.table that would have been used to create
#'   ggplot is returned instead.
#'
#' @return ggplot containing images summarizing t-sne regions.  if return_data
#'   is TRUE, instead the data.table containing plot info is returned.
#' @export
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' summary_dt = prep_summary(profile_dt, tsne_dt, 4, facet_by = "tall_var")
#' img_res = prep_images(summary_dt, 4, facet_by = "tall_var")
#' #zoom on top-right quadrant
#' summary_dt.zoom = prep_summary(profile_dt, tsne_dt, 4, facet_by = "tall_var",
#'     xrng = c(0, .5), yrng = c(0, .5))
#' img_res.zoom = prep_images(summary_dt.zoom, 4, facet_by = "tall_var",
#'     xrng = c(0, .5), yrng = c(0, .5))
#' plot_summary_raster_byCell(img_res$image_dt,
#'               x_points = img_res$x_points)
#' plot_summary_raster_byCell(img_res.zoom$image_dt,
#'               x_points = img_res.zoom$x_points,
#'               xrng = img_res.zoom$xrng,
#'               yrng = img_res.zoom$yrng)
plot_summary_raster_byCell = function(image_dt,
                                      x_points,
                                      y_points = x_points,
                                      p = NULL,
                                      line_color_mapping = NULL,
                                      xrng = c(-.5, .5),
                                      yrng = c(-.5, .5),
                                      N_floor = 0,
                                      N_ceiling = NULL,
                                      min_size = .3,
                                      return_data = FALSE) {
    tall_var = bx = by = NULL

    # image_dt = merge(image_dt[, list(bx, by, plot_id, png_file, tx, ty)],
    #                  postion[, list(N = .N), list(tall_var, bx, by)], by = c("bx", "by"), allow.cartesian = TRUE)
    image_dt = set_image_rects(
        image_dt,
        x_points = x_points,
        y_points = y_points,
        xrng = xrng,
        yrng = yrng,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_size
    )
    if (return_data)
        return(image_dt)
    if (is.null(p))
        p = ggplot()
    if(!is.null(line_color_mapping)){
        col_dt = image_dt[, .(tx, ty, wide_var = rep(names(line_color_mapping), each = nrow(image_dt)))]
        p = p + geom_point(data = col_dt,
                           aes(x = tx, y = ty, color = wide_var)) +
            scale_color_manual(values = line_color_mapping)
    }
    p = p +
        geom_image.rect(data = image_dt,
                        aes(
                            xmin = xmin,
                            xmax = xmax,
                            ymin = ymin,
                            ymax = ymax,
                            image = png_file
                        )) + #, color = rgb(0,0,1,.2)) +
        geom_rect(
            data = image_dt,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax
            ),
            fill = NA,
            color = "black"
        ) +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        facet_wrap("tall_var", drop = FALSE)
    p
}

#' plot_summary_glyph
#'
#' @param summary_dt results from prep_summary()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
#' @param ylim ylimits per glyph
#' @param p an existing ggplot to overlay images onto.  Default of NULL starts a
#'   new plot.
#' @param N_floor The value of N to consider 0.  bins with N values <= N_floor
#'   will be ignored.
#' @param N_ceiling The value of N to consider 1.  bins with N values >=
#'   N_ceiling will have images drawn at full size.
#' @param min_size Numeric (0, 1]. The minimum size images to draw.  The default
#'   of .3 draws images for all bins with N values >= 30% of the way from
#'   N_floor to N_ceiling.
#' @param return_data if TRUE, data.table that would have been used to create
#'   ggplot is returned instead.
#' @param color_mapping mapping for scale_color_manual
#'
#' @return a ggplot containing glyphs of local profile summaries arranged in
#'   t-sne space.
#' @export
#' @importFrom GGally glyphs
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' n_points = 12
#' summary_dt = prep_summary(profile_dt = profile_dt,
#'     position_dt = tsne_dt, x_points = n_points)
#' plot_summary_glyph(summary_dt,
#'     x_points = n_points)
plot_summary_glyph = function(
    summary_dt,
    x_points,
    y_points = x_points,
    xrng = c(-.5, .5),
    yrng = c(-.5, .5),
    ylim = NULL,
    p = NULL,
    N_floor = 0,
    N_ceiling = NULL,
    min_size = .3,
    return_data = FALSE,
    color_mapping = NULL){
    group_size = gx = gy = gid = NULL #bindings for data.table
    summary_dt = set_size(summary_dt, N_floor, N_ceiling, size.name = "group_size")
    summary_dt = summary_dt[group_size >= min_size]

    if(is.null(ylim)){
        ylim = range(summary_dt$y)
    }
    ylim = range(ylim)
    #apply ylimits
    summary_dt[y < min(ylim), y := min(ylim)]
    summary_dt[y > max(ylim), y := max(ylim)]

    #apply sizing
    summary_dt[, x := x * group_size]
    summary_dt[, y := y * group_size]

    xs = bin_values_centers(x_points, rng = xrng)
    ys = bin_values_centers(y_points, rng = yrng)
    summary_dt[, tx := xs[bx]]
    summary_dt[, ty := ys[by]]

    down_scale = max(summary_dt$group_size)

    glyph_dt = as.data.table(GGally::glyphs(summary_dt,
                                            x_major = "tx",
                                            x_minor = "x",
                                            y_major = "ty",
                                            y_minor = "y",
                                            width = diff(xrng)/x_points*.95*down_scale,
                                            height = diff(yrng)/y_points*.95*down_scale))
    if(return_data){
        return(glyph_dt)
    }

    if(is.null(color_mapping)){
        if(is.factor(summary_dt$wide_var)){
            uwide_vars = levels(summary_dt$wide_var)
        }else if(is.character(summary_dt$wide_var)){
            uwide_vars = unique(summary_dt$wide_var)
        }

        color_mapping = seqsetvis::safeBrew(length(uwide_vars))
    }
    if(is.null(p)){
        p = ggplot()
    }
    p +
        geom_path(data = glyph_dt,
                  aes(gx, gy, group = paste(gid, wide_var), color = wide_var)) +
        labs(x = "tx", y = "ty") +
        scale_color_manual(values =  color_mapping) +
        coord_cartesian(xrng, yrng)
}
