#' stsPlotSummaryProfiles
#'
#' @param profile_dt a tidy data.table for profile data as retrieved by
#'   stsFetchTsneInput.  Expected variable names are id, cell, mark, x, and y.
#' @param position_dt a tidy data.table containing t-sne embedding.  Expected
#'   variable names are tx, ty, id, and cell.
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#'   Defaults to same value as x_points.
#' @param q_cells character vector of cells to plot. Default of NULL plots all.
#' @param q_marks character vector of marks to plot.  Default of NULL plots all.
#' @param xrng view domain in x dimension, default is range of position_dt$tx.
#' @param yrng view domain in y dimension, default is range of position_dt$ty.
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
#' @param facet_by character. varaible name to facet profile_dt by when
#'   constructing images. The only valid non-null value with seqtsne functions
#'   is "cell".
#' @param line_color_mapping named vector of line color.  Names correspond to values of
#'   profile_dt 'mark' variable and values are colors.
#' @param vertical_facet_mapping named vector of groups.  vertical_facet_mapping names are marks and order
#'   of first occurence determines vertical position top to bottom.
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
    q_cells = NULL,
    q_marks = NULL,
    ## view domain
    xrng = range(position_dt$tx),
    yrng = range(position_dt$ty),
    ## image file path and output
    rname = digest::digest(
        list(
            profile_dt,
            position_dt,
            x_points,
            y_points,
            apply_norm,
            ylim,
            line_color_mapping,
            n_splines,
            ma_size,
            facet_by
        )
    ),
    odir = file.path("tsne_images/", rname),
    force_rewrite = FALSE,
    n_cores = getOption("mc.cores", 1),
    ## preplot transformation
    apply_norm = TRUE,
    ylim = c(0, 1),
    ma_size = 2,
    n_splines = 10,
    ## plot organization
    p = NULL,
    facet_by = NULL,
    line_color_mapping = c("H3K4me3" = "forestgreen",
                    "H3K27me3" = "firebrick1"),
    vertical_facet_mapping = NULL,
    ## sizing and level of detail
    N_floor = 0,
    N_ceiling = NULL,
    min_size = .3,
    ## output
    return_data = FALSE) {
    if(is.null(q_cells)){
        if(is.factor(profile_dt$cell)){
            q_cells = levels(profile_dt$cell)
        }else{
            q_cells = sort(unique(profile_dt$cell))
        }
    }
    if(is.null(q_marks)){
        if(is.factor(profile_dt$mark)){
            q_marks = levels(profile_dt$mark)
        }else{
            q_marks = sort(unique(profile_dt$mark))
        }
    }
    stopifnot(q_cells %in% profile_dt$cell)
    stopifnot(q_marks %in% profile_dt$mark)
    prof_dt = copy(profile_dt[cell %in% q_cells & mark %in% q_marks])
    pos_dt = copy(position_dt[cell %in% q_cells])
    prof_dt$cell = factor(prof_dt$cell, levels = q_cells)
    prof_dt$mark = factor(prof_dt$mark, levels = q_marks)
    pos_dt$cell = factor(pos_dt$cell, levels = q_cells)
    #prepare images
    img_res = prep_images(
        profile_dt = prof_dt,
        position_dt = pos_dt,
        xrng = xrng,
        yrng = yrng,
        x_points = x_points,
        y_points = y_points,
        rname = rname,
        odir = odir,
        force_rewrite = force_rewrite,
        apply_norm = apply_norm,
        ylim = ylim,
        facet_by = facet_by,
        ma_size = ma_size,
        n_splines = n_splines,
        n_cores = n_cores,
        line_color_mapping = line_color_mapping,
        vertical_facet_mapping = vertical_facet_mapping
    )
    if (is.null(facet_by)) {
        plot_tsne_img(
            image_dt = img_res$image_dt,
            xrng = xrng,
            yrng = yrng,
            x_points = x_points,
            y_points = y_points,
            p = p,
            N_floor = N_floor,
            N_ceiling = N_ceiling,
            min_size = min_size,
            return_data = return_data
        )
    } else{
        plot_tsne_img_byCell(
            image_dt = img_res$image_dt,
            tsne_dt = img_res$tsne_dt,
            xrng = xrng,
            yrng = yrng,
            x_points = x_points,
            y_points = y_points,
            p = p,
            N_floor = N_floor,
            N_ceiling = N_ceiling,
            min_size = min_size,
            return_data = return_data
        )
    }


}

#' prep_images
#'
#' Prepares images summarizing profiles in t-sne regions and returns ggimage
#' compatible data.frame for plotting.
#'
#' @param profile_dt a tidy data.table for profile data as retrieved by
#'   stsFetchTsneInput.  Expected variable names are id, cell, mark, x, and y.
#' @param position_dt a tidy data.table containing t-sne embedding.  Expected
#'   variable names are tx, ty, id, and cell.
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
#'   constructing images. The only valid non-null value with seqtsne functions
#'   is "cell".
#' @param n_splines number of points to interpolate with splines.
#' @param ma_size moving average size when smoothing profiles.
#' @param n_cores number of cores to use when writing images.  Default is
#' value of mc.cores option if set or 1.
#' @param line_color_mapping named vector of line color.  Names correspond to values
#' of profile_dt 'mark' variable and values are colors.
#' @param vertical_facet_mapping named vector of vertical facet for data
#'
#' @return data.table with variables
#' @export
#' @importFrom seqsetvis applySpline
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' img_res = prep_images(profile_dt, tsne_dt, 4)
#' #zoom on top-right quadrant
#' img_res.zoom = prep_images(profile_dt, tsne_dt, 4, xrng = c(0, .5), yrng = c(0, .5))
#' #use results with plot_tsne_img() to make plots
prep_images = function(profile_dt,
                       position_dt,
                       x_points,
                       y_points = x_points,
                       xrng = range(position_dt$tx),
                       yrng = range(position_dt$ty),
                       rname = digest::digest(
                           list(
                               profile_dt,
                               position_dt,
                               x_points,
                               y_points,
                               apply_norm,
                               ylim,
                               line_color_mapping,
                               n_splines,
                               ma_size,
                               facet_by
                           )
                       ),
                       odir = file.path("tsne_images/", rname),
                       force_rewrite = FALSE,
                       apply_norm = TRUE,
                       ylim = c(0, 1),
                       facet_by = NULL,
                       # view_rect = list(),
                       ma_size = 2,
                       n_splines = 10,
                       n_cores = getOption("mc.cores", 1),
                       line_color_mapping = c("H3K4me3" = "forestgreen",
                                       "H3K27me3" = "firebrick1"),
                       vertical_facet_mapping = NULL) {
    if (!all(unique(profile_dt$mark) %in% names(line_color_mapping))) {
        missing_colors = setdiff(unique(profile_dt$mark), names(line_color_mapping))
        stop(
            "line_color_mapping is missing assignments for: ",
            paste(missing_colors, collapse = ", ")
        )
    }
    stopifnot()
    # stopifnot(is.list(view_rect))
    # if(is.null(view_rect$xmin))
    position_dt = copy(position_dt)
    position_dt = position_dt[tx >= min(xrng) &
                                  tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng)]
    #use positional info from position_dt to bin points
    position_dt[, bx := bin_values(tx, x_points, xrng = xrng)]
    position_dt[, by := bin_values(ty, y_points, xrng = yrng)]
    #merge binning info to profiles
    mdt = merge(
        profile_dt,
        position_dt[, list(bx, by, cell, id)],
        allow.cartesian = TRUE,
        by = intersect(colnames(profile_dt), c("cell", "id"))
    )#, by = c("cell", "id"))
    if (is.null(mdt$mark))
        mdt$mark = "signal"

    if (is.null(facet_by)) {
        mdt = mdt[, list(y = mean(y)), list(bx, by, x, mark)]
    } else{
        mdt = mdt[, list(y = mean(y)), list(bx, by, x, mark, get(facet_by))]
        colnames(mdt)[colnames(mdt) == "get"] = facet_by
    }
    #each combination of bx and by is a unique plot_id
    mdt[, plot_id := paste(bx, by, sep = "_")]
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    if (is.null(facet_by)) {
        img_dt = unique(mdt[, list(bx, by, plot_id)])
        img_dt[, png_file := file.path(odir, paste0(plot_id, ".png"))]
    } else{
        img_dt = unique(mdt[, list(bx, by, plot_id, get(facet_by))])
        colnames(img_dt)[4] = facet_by
        img_dt[, png_file := file.path(odir, paste0(get(facet_by), "_", plot_id, ".png"))]
    }

    xs = bin_values_centers(position_dt$tx, x_points, xrng = xrng)
    ys = bin_values_centers(position_dt$ty, y_points, xrng = yrng)

    # plot(expand.grid(xs, ys), xlim = xrng, ylim = yrng)
    # rect(min(xrng), min(yrng), max(xrng), max(yrng), col = rgb(0,0,1,.1))

    img_dt[, tx := xs[bx]]
    img_dt[, ty := ys[by]]


    if (apply_norm) {
        mdt[, ynorm := y / quantile(y, .95), by = list(mark)]
        mdt[ynorm > 1, ynorm := 1]
    } else{
        mdt[, ynorm := y]
    }


    if (force_rewrite) {
        file.remove(img_dt$png_file[file.exists(img_dt$png_file)])
    }

    if (is.null(facet_by)) {
        img_dt = merge(img_dt, position_dt[, .N, list(bx, by)])
    } else{
        tmp_dt = position_dt[, .N, list(bx, by, get(facet_by))]
        colnames(tmp_dt)[colnames(tmp_dt) == "get"] = facet_by
        img_dt = merge(img_dt, tmp_dt)
    }

    if (is.null(vertical_facet_mapping)) {
        mdt$group = 1
    } else{
        if (!all(unique(mdt$mark) %in% names(vertical_facet_mapping))) {
            missing_groups = setdiff(unique(mdt$mark), names(vertical_facet_mapping))
            stop(
                "vertical_facet_mapping is missing assignments for: ",
                paste(missing_groups, collapse = ", ")
            )
        }
        mdt$group = factor(vertical_facet_mapping[mdt$mark], levels = unique(vertical_facet_mapping))
    }


    if (any(!file.exists(img_dt$png_file))) {
        plot_info = lapply(which(!file.exists(img_dt$png_file)), function(i) {
            # hidden = parallel::mclapply(which(!file.exists(img_dt$png_file)), function(i){}, mc.cores = n_cores)
            fpath = img_dt$png_file[i]
            p_id = img_dt$plot_id[i]
            if (is.null(facet_by)) {
                pdt = mdt[plot_id == p_id]
            } else{
                pdt = mdt[plot_id == p_id & get(facet_by) == img_dt[[facet_by]][i]]
            }

            # pdt[, ysm := seqsetvis:::movingAverage(ynorm, n = ma_size), by = list(mark)]
            pdt[, ysm := movingAverage(ynorm, n = ma_size), by = list(mark)]

            pdt = seqsetvis::applySpline(pdt,
                                         n = n_splines,
                                         by_ = "mark",
                                         y_ = "ysm")
            list(pdt, fpath)
        })
        # hidden = lapply(plot_info, function(x){
        # figure out how not to copy global env
        hidden = parallel::mclapply(plot_info, function(x) {
            fpath = x[[2]]
            pdt = x[[1]]
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = list(mark)]

            p = ggplot(pdt,
                       aes(
                           x = x,
                           y = ysm,
                           ymin = 0,
                           ymax = ysm,
                           color = mark,
                           fill = mark
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

    if (is.factor(position_dt$cell)) {
        if (!is.null(img_dt$cell)) {
            img_dt$cell = factor(img_dt$cell, levels = levels(position_dt$cell))
        }
        if (!is.null(mdt$cell)) {
            mdt$cell = factor(mdt$cell, levels = levels(position_dt$cell))
        }
    }

    return(
        list(
            image_dt = img_dt[],
            summary_profile_dt = mdt[],
            tsne_dt = position_dt[],
            x_points = x_points,
            y_points = y_points,
            xrng = xrng,
            yrng = yrng
        )
    )
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
#'
#' img_res = prep_images(profile_dt,
#'                       tsne_dt,
#'                       x_points = 4)
#' img_rect = set_image_rects(img_res$image_dt,
#'                            x_points = img_res$x_points,
#'                            y_points = img_res$y_points,
#'                            xrng = img_res$xrng,
#'                            yrng = img_res$yrng)
#' ggplot(img_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + geom_rect()
#' ggplot(img_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, image = png_file)) + geom_image.rect()
set_image_rects = function(image_dt,
                           x_points,
                           y_points,
                           xrng,
                           yrng,
                           N_floor = 0,
                           N_ceiling = NULL,
                           min_size = .3) {
    if (is.null(N_ceiling)) {
        N_ceiling = max(image_dt$N)
    }
    image_dt[, img_size := N]
    image_dt[img_size > N_ceiling, img_size := N_ceiling]
    image_dt[img_size < N_floor, img_size := N_floor]

    image_dt[, img_size := img_size - N_floor]
    image_dt[, img_size := img_size / N_ceiling]
    image_dt[, img_size := img_size]
    image_dt = image_dt[img_size >= min_size]

    xspc = diff(xrng) / x_points / 2
    yspc = diff(yrng) / y_points / 2

    image_dt[, xmin := tx - xspc * img_size]
    image_dt[, xmax := tx + xspc * img_size]
    image_dt[, ymin := ty - yspc * img_size]
    image_dt[, ymax := ty + yspc * img_size]
    image_dt[]
}

#' plot_tsne_img
#'
#' @param image_dt $image_dt of result from prep_images()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
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
#'
#' @return ggplot containing images summarizing t-sne regions.  if return_data
#'   is TRUE, instead the data.table containing plot info is returned.
#' @export
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' img_res = prep_images(profile_dt, tsne_dt, 4)
#' #zoom on top-right quadrant
#' img_res.zoom = prep_images(profile_dt, tsne_dt, 4, xrng = c(0, .5), yrng = c(0, .5))
#' plot_tsne_img(img_res$image_dt,
#'               x_points = img_res$x_points)
#' plot_tsne_img(img_res.zoom$image_dt,
#'               x_points = img_res.zoom$x_points,
#'               xrng = img_res.zoom$xrng,
#'               yrng = img_res.zoom$yrng)
plot_tsne_img = function(image_dt,
                         x_points,
                         y_points = x_points,
                         xrng = c(-.5, .5),
                         yrng = c(-.5, .5),
                         p = NULL,

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

#' plot_tsne_img
#'
#' @param image_dt $image_dt of result from prep_images()
#' @param x_points numeric.  number of grid points to use in x dimension.
#' @param y_points numeric.  number of grid points to use in y dimension.
#' @param xrng view domain in x dimension.
#' @param yrng view domain in y dimension.
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
#'
#' @return ggplot containing images summarizing t-sne regions.  if return_data
#'   is TRUE, instead the data.table containing plot info is returned.
#' @export
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' img_res = prep_images(profile_dt, tsne_dt, 4, facet_by = "cell")
#' #zoom on top-right quadrant
#' img_res.zoom = prep_images(profile_dt, tsne_dt, 4, xrng = c(0, .5), yrng = c(0, .5))
#' plot_tsne_img_byCell(img_res$image_dt,
#'               x_points = img_res$x_points)
#' plot_tsne_img_byCell(img_res.zoom$image_dt,
#'               x_points = img_res.zoom$x_points,
#'               xrng = img_res.zoom$xrng,
#'               yrng = img_res.zoom$yrng)
plot_tsne_img_byCell = function(image_dt,
                                tsne_dt,
                                x_points,
                                y_points = x_points,
                                p = NULL,
                                xrng = c(-.5, .5),
                                yrng = c(-.5, .5),
                                N_floor = 0,
                                N_ceiling = NULL,
                                min_size = .3,
                                return_data = FALSE) {
    cell = bx = by = NULL

    # image_dt = merge(image_dt[, list(bx, by, plot_id, png_file, tx, ty)],
    #                  tsne_dt[, list(N = .N), list(cell, bx, by)], by = c("bx", "by"), allow.cartesian = TRUE)
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
        facet_wrap("cell", drop = FALSE)
    p
}
