
#' run_tsne
#'
#' @param profile_dt Tidy data.table of profile information. As returned by seqsetvis::ssvFetchBam.
#' @param perplexity perplexity value for t-SNE. Passed to Rtsne::Rtsne.
#' @param n_cores Number of threads for running t-SNE.
#' @param high_topright If TRUE, flip tx/ty as needed such that most points are in the top-right quadrant.
#' @param norm1 If TRUE, rescale tx and ty to be centered at 0,0 with total size of 1.
#' @param Y_init Optional initial coordinates for tx and ty.
#' @param wide_var Wide variable.  Spreads matrix into columns. Each regions tSNE position will be based on one profile for every value of wide_var.
#' @param tall_var Tall variable. Repeats matrix entries down rows. Each region will appear once per value of tall_var in final tSNE.
#'
#' @return Return data.table mapping each id_var entry to tSNE space: defined by tx and ty
#' @import Rtsne
#'
#' @examples
#' data("profile_dt")
#' setalloccol(profile_dt)
#' chiptsne:::run_tsne(profile_dt, wide_var = "wide_var", tall_var = "tall_var")
run_tsne = function (profile_dt,
                     perplexity = 100,
                     n_cores = getOption("mc.cores", 1),
                     high_topright = TRUE,
                     norm1 = TRUE,
                     Y_init = NULL,
                     x_var = "x",
                     y_var = "y",
                     id_var = "id",
                     wide_var = "name",
                     tall_var = "tall_none",
                     fun.aggregate = mean) {
    valid_vars = c(ifelse(tall_var == "tall_none", character(), tall_var), wide_var, x_var, y_var, id_var)
    valid_vars = valid_vars[!is.na(valid_vars)]
    stopifnot(valid_vars %in% colnames(profile_dt))

    if(is.null(profile_dt[[tall_var]])){
        set(profile_dt, j = tall_var, value = "none")
    }
    tsne_mat = dt2mat(profile_dt,
                      unique(profile_dt[[wide_var]]),
                      x_var = x_var,
                      y_var = y_var,
                      id_var = id_var,
                      wide_var_ = wide_var,
                      tall_var_ = tall_var,
                      fun.aggregate = fun.aggregate)
    bad_col = apply(tsne_mat, 2, stats::var) == 0
    if (any(bad_col)) {
        warning("zero variance columns detected in tsne matrix input.",
                "\n", round(sum(bad_col)/length(bad_col) * 100,
                            2), "% of columns removed.")
        tsne_mat = tsne_mat[, !bad_col, drop = FALSE]
    }
    if (is.data.table(Y_init)) {
        if(tall_var != "tall_none"){
            Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init[[id_var]],
                                                                           Y_init[[tall_var]]))
        }else{
            Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init[[id_var]], "none"))
        }

        Y_init = Y_init[rownames(tsne_mat), ]
        stopifnot(ncol(Y_init) == 2)
    }
    else {
        stopifnot(is.null(Y_init))
    }
    max_perplexity = floor(nrow(tsne_mat)/4)
    if(perplexity > max_perplexity){
        perplexity = max_perplexity
        warning("Reducing perplexity to ", perplexity, " to accommodate data of ", nrow(tsne_mat), " rows.")
    }
    res_tsne = Rtsne::Rtsne(tsne_mat,
                            Y_init = Y_init,
                            num_threads = n_cores,
                            perplexity = perplexity,
                            check_duplicates = FALSE)
    tsne_dt = as.data.table(res_tsne$Y)
    setnames(tsne_dt, c("tx", "ty"))
    tsne_dt$rn = rownames(tsne_mat)

    tsne_dt[, `:=`(c(id_var, tall_var), tstrsplit(rn, " ", keep = seq(2)))]


    if (norm1) {
        tsne_dt$tx = rescale_capped(tsne_dt$tx) - 0.5
        tsne_dt$ty = rescale_capped(tsne_dt$ty) - 0.5
    }
    if (high_topright) {#flip tx/ty if needed so that
        rs = rowSums(tsne_mat)
        tsne_dt$rs = rs[tsne_dt$rn]
        x_cutoff = mean(range(tsne_dt$tx))
        x_flip = sum(tsne_dt[tx > x_cutoff]$rs) < sum(tsne_dt[tx < x_cutoff]$rs)
        if (x_flip) {
            tsne_dt[, `:=`(tx, max(tx) - tx + min(tx))]
        }
        y_cutoff = mean(range(tsne_dt$ty))
        y_flip = sum(tsne_dt[ty > y_cutoff]$rs) < sum(tsne_dt[ty < y_cutoff]$rs)
        if (y_flip) {
            tsne_dt[, `:=`(ty, max(ty) - ty + min(ty))]
        }
        tsne_dt$rs = NULL
    }
    tsne_dt$rn = NULL
    tsne_dt[]
}

rescale_capped = function (x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

dt2mat = function (prof_dt,
                   wide_values,
                   x_var = "x",
                   y_var = "y",
                   id_var = "id",
                   wide_var_ = "name",
                   tall_var_ = "tall_none",
                   fun.aggregate = mean){
    if(is.null(prof_dt[[tall_var_]])){
        if(tall_var_ == "tall_none"){
            prof_dt[[tall_var_]] = "none"
        }
    }
    stopifnot(c(id_var, wide_var_, tall_var_, x_var, y_var) %in% colnames(prof_dt))
    dt = reshape2::dcast(data = prof_dt[get(wide_var_) %in% wide_values],
                    formula = paste0(id_var, "+", tall_var_, "~", x_var, "+", wide_var_),
                    value.var = y_var,
                    fun.aggregate = fun.aggregate)
    wide_mat = as.matrix(dt[, -seq_len(2)])
    rownames(wide_mat) = paste(dt[[id_var]], dt[[tall_var_]])
    if(!all(table(colnames(wide_mat)) == 1)){
        stop("Not all colnames of matrix were unique.")
    }
    if(!all(table(rownames(wide_mat)) == 1)){
        stop("Not all colnames of matrix were unique.")
    }
    wide_mat
}
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
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' stsPlotSummaryProfiles(profile_dt, tsne_dt, x_points = 4)
stsPlotSummaryProfiles = function (profile_dt,
                                   position_dt,
                                   x_points,
                                   y_points = x_points,
                                   x_var = "x",
                                   y_var = "y",
                                   extra_vars = character(),
                                   id_var = "id",
                                   wide_var = "name",
                                   tall_var = "tall_none",
                                   q_tall_values = NULL,
                                   q_wide_values = NULL,
                                   xrng = NULL,
                                   yrng = NULL,
                                   plot_type = c("glyph", "raster")[1],
                                   rname = NULL,
                                   odir = NULL,
                                   force_rewrite = FALSE,
                                   n_cores = getOption("mc.cores", 1),
                                   apply_norm = TRUE,
                                   ylim = c(0, 1),
                                   ma_size = 2,
                                   n_splines = 10,
                                   p = NULL,
                                   facet_byCell = FALSE,
                                   line_color_mapping = NULL,
                                   vertical_facet_mapping = NULL,
                                   N_floor = 0,
                                   N_ceiling = NULL,
                                   min_size = 0.3,
                                   return_data = FALSE) {
    if(is.null(profile_dt[[tall_var]])){
        profile_dt[[tall_var]] = "none"
    }
    if (is.null(q_tall_values)) {
        if (is.factor(profile_dt[[tall_var]])) {
            q_tall_values = levels(profile_dt[[tall_var]])
        }
        else {
            q_tall_values = sort(unique(profile_dt[[tall_var]]))
        }
    }
    if (is.null(q_wide_values)) {
        if (is.factor(profile_dt[[wide_var]])) {
            q_wide_values = levels(profile_dt[[wide_var]])
        }
        else {
            q_wide_values = sort(unique(profile_dt[[wide_var]]))
        }
    }
    if (is.null(line_color_mapping)) {
        line_color_mapping = seqsetvis::safeBrew(length(unique(profile_dt[[wide_var]])))
        names(line_color_mapping) = unique(profile_dt[[wide_var]])
    }
    line_color_mapping = line_color_mapping[names(line_color_mapping) %in%
                                                q_wide_values]
    if (is.factor(profile_dt[[tall_var]])) {
        stopifnot(q_tall_values %in% levels(profile_dt[[tall_var]]))
    }
    else {
        stopifnot(q_tall_values %in% unique(profile_dt[[tall_var]]))
    }
    if (is.factor(profile_dt[[wide_var]])) {
        stopifnot(q_wide_values %in% levels(profile_dt[[wide_var]]))
    }
    else {
        stopifnot(q_wide_values %in% unique(profile_dt[[wide_var]]))
    }
    if (!plot_type %in% c("glyph", "raster")) {
        stop("plot_type (\"", plot_type, "\") must be one of \"glyph\" or \"raster\".")
    }
    k = profile_dt[[tall_var]] %in% q_tall_values &
        profile_dt[[wide_var]] %in% q_wide_values
    prof_dt = copy(profile_dt[k == TRUE])
    pos_dt = copy(position_dt[position_dt[[tall_var]] %in% q_tall_values])
    if(is.null(xrng)){
        xrng = range(pos_dt$tx)
    }
    if(is.null(yrng)){
        yrng = range(pos_dt$ty)
    }
    if (is.null(rname)) {
        rname = digest::digest(list(prof_dt,
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
                                    facet_byCell))
    }
    if(is.null(odir)){
        odir = file.path(tempdir(), rname)
    }
    prof_dt[[tall_var]] = factor(prof_dt[[tall_var]], levels = q_tall_values)
    prof_dt[[wide_var]] = factor(prof_dt[[wide_var]], levels = q_wide_values)
    pos_dt[[tall_var]] = factor(pos_dt[[tall_var]], levels = q_tall_values)
    #### raster ####
    if (plot_type == "raster") {
        #### raster - no facet ####
        if (!facet_byCell) {
            summary_dt = prep_summary(profile_dt = prof_dt,
                                      position_dt = pos_dt,
                                      x_points = x_points,
                                      y_points = y_points,
                                      xrng = xrng,
                                      yrng = yrng,
                                      facet_by = NULL,
                                      x_var = x_var,
                                      y_var = y_var,
                                      id_var = id_var,
                                      wide_var = wide_var,
                                      tall_var = tall_var,
                                      extra_vars = extra_vars)
            img_res = prep_images(summary_dt = summary_dt,
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
                                  vertical_facet_mapping = vertical_facet_mapping,
                                  wide_var = wide_var,
                                  x_var = x_var,
                                  y_var = y_var)
            plot_summary_raster(image_dt = img_res$image_dt,
                                xrng = xrng,
                                yrng = yrng,
                                x_points = x_points,
                                y_points = y_points,
                                p = p,
                                line_color_mapping = img_res$line_color_mapping,
                                N_floor = N_floor,
                                N_ceiling = N_ceiling,
                                min_size = min_size,
                                return_data = return_data,
                                wide_var = wide_var)
            #### raster - facet ####
        } else {
            summary_dt = prep_summary(profile_dt = prof_dt,
                                      position_dt = pos_dt,
                                      x_points = x_points,
                                      y_points = y_points,
                                      xrng = xrng,
                                      yrng = yrng,
                                      facet_by = tall_var,
                                      x_var = x_var,
                                      y_var = y_var,
                                      id_var = id_var,
                                      wide_var = wide_var,
                                      extra_vars = extra_vars)
            img_res = prep_images(summary_dt = summary_dt,
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
                                  vertical_facet_mapping = vertical_facet_mapping,
                                  wide_var = wide_var,
                                  x_var = x_var,
                                  y_var = y_var)
            plot_summary_raster_byCell(image_dt = img_res$image_dt,
                                       xrng = xrng,
                                       yrng = yrng,
                                       x_points = x_points,
                                       y_points = y_points,
                                       p = p,
                                       line_color_mapping = img_res$line_color_mapping,
                                       N_floor = N_floor,
                                       N_ceiling = N_ceiling,
                                       min_size = min_size,
                                       return_data = return_data)
        }
    }
    #### glyph ####
    else if (plot_type == "glyph") {
        #### glyph - no facet ####
        if (!facet_byCell) {
            summary_dt = prep_summary(profile_dt = prof_dt,
                                      position_dt = pos_dt,
                                      x_points = x_points,
                                      y_points = y_points,
                                      xrng = xrng,
                                      yrng = yrng,
                                      facet_by = NULL,
                                      x_var = x_var,
                                      y_var = y_var,
                                      id_var = id_var,
                                      wide_var = wide_var,
                                      tall_var = tall_var,
                                      extra_vars = extra_vars)
            plot_summary_glyph(summary_dt = summary_dt,
                               xrng = xrng,
                               yrng = yrng,
                               x_points = x_points,
                               y_points = y_points,
                               p = p,
                               ylim = ylim,
                               N_floor = N_floor,
                               N_ceiling = N_ceiling,
                               min_size = min_size,
                               color_mapping = line_color_mapping,
                               return_data = return_data,
                               x_var = x_var,
                               y_var = y_var,
                               wide_var = wide_var,
                               extra_vars = extra_vars)
            #### glyph - facet ####
        } else {
            summary_dt_l = lapply(q_tall_values, function(cl) {
                prep_summary(prof_dt[get(tall_var) == cl],
                             position_dt = pos_dt,
                             x_points = x_points,
                             y_points = y_points,
                             xrng = xrng,
                             yrng = yrng,
                             facet_by = NULL,
                             x_var = x_var,
                             y_var = y_var,
                             id_var = id_var,
                             wide_var = wide_var,
                             tall_var = tall_var,
                             extra_vars = extra_vars)
            })
            names(summary_dt_l) = q_tall_values
            summary_dt = rbindlist(summary_dt_l,
                                   use.names = TRUE,
                                   idcol = tall_var)
            if (return_data) {
                plot_summary_glyph(summary_dt = summary_dt,
                                   xrng = xrng,
                                   yrng = yrng,
                                   x_points = x_points,
                                   y_points = y_points,
                                   ylim = ylim,
                                   N_floor = N_floor,
                                   N_ceiling = N_ceiling,
                                   min_size = min_size,
                                   color_mapping = line_color_mapping,
                                   return_data = return_data,
                                   x_var = x_var,
                                   y_var = y_var,
                                   wide_var = wide_var)
            }
            else {
                plot_summary_glyph(summary_dt,
                                   x_points = x_points,
                                   y_points = y_points,
                                   xrng = xrng,
                                   yrng = yrng,
                                   ylim = ylim,
                                   N_floor = N_floor,
                                   N_ceiling = N_ceiling,
                                   min_size = min_size,
                                   color_mapping = line_color_mapping,
                                   x_var = x_var,
                                   y_var = y_var,
                                   wide_var = wide_var,
                                   extra_vars = extra_vars) +
                    facet_wrap(paste0("~", tall_var))
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
prep_summary = function (profile_dt,
                         position_dt,
                         x_points,
                         y_points = x_points,
                         xrng = range(position_dt$tx),
                         yrng = range(position_dt$ty),
                         facet_by = NULL,
                         x_var = "x",
                         y_var = "y",
                         id_var = "id",
                         wide_var = "name",
                         tall_var = "tall_none",
                         extra_vars = character()){
    position_dt = copy(position_dt[tx >= min(xrng) & tx <= max(xrng) &
                                       ty >= min(yrng) & ty <= max(yrng)])
    position_dt = position_dt[get(id_var) %in% unique(profile_dt[[id_var]])]
    if (is.null(position_dt$bx))
        position_dt[, `:=`(bx, bin_values(tx, x_points,
                                          xrng = xrng))]
    if (is.null(position_dt$by))
        position_dt[, `:=`(by, bin_values(ty, y_points,
                                          xrng = yrng))]
    summary_dt = merge(profile_dt, position_dt[, c("bx", "by", tall_var, id_var), with = FALSE],
                       allow.cartesian = TRUE,
                       by = intersect(colnames(profile_dt), c(tall_var, id_var)))
    if (is.null(summary_dt[[wide_var]]))
        summary_dt[[wide_var]] = "signal"
    if (is.null(facet_by)) {
        summary_dt = summary_dt[, list(y_tmp_ = mean(get(y_var))), c(unique(c("bx", "by", x_var, wide_var, extra_vars)))]
    }
    else {
        summary_dt = summary_dt[, list(y_tmp_ = mean(get(y_var))), c(unique(c("bx", "by", x_var, wide_var, facet_by, extra_vars)))]
    }
    setnames(summary_dt, "y_tmp_", y_var)
    N_dt = position_dt[, .(.N), by = .(bx, by)]
    summary_dt = merge(summary_dt, N_dt, by = c("bx", "by"))
    summary_dt[, `:=`(plot_id, paste(bx, by, sep = "_"))]
    summary_dt[]
}

bin_values = function (x, n_bins, xrng = range(x))
{
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_bins - 1e-05)) +
        1
}

bin_values_centers = function (n_bins, rng)
{
    if (length(rng) != 2)
        rng = range(rng)
    stopifnot(length(rng) == 2)
    xspc = diff(rng)/n_bins/2
    xs = seq(min(rng) + xspc, max(rng) - xspc, diff(rng)/(n_bins))
    xs
}
