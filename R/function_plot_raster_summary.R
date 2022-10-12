

#' aggregate_signals
#'
#' aggregate_signals is retrieves a single representative value for a genomic
#' region per sample by default this is the maximum anywhere in the region but
#' this can be overriden using xmin/xmax and agg_FUN
#' @param profile_dt
#' @param agg_FUN
#' @param xmin
#' @param xmax
#' @param by_
#'
#' @return
#'
#' @examples
#' aggregate_signals(profile_dt)
#' aggregate_signals(profile_dt, xmin = -.2, xmax = .2, agg_FUN = mean)
aggregate_signals = function(profile_dt,
                             agg_FUN = max,
                             y_ = "y",
                             yout_ = "ynorm",
                             xmin = -Inf,
                             xmax = Inf,
                             by_ = c("tall_var", "wide_var")){
    agg_dt = profile_dt[x >= xmin & x <= xmax,
                        .(val_ = agg_FUN(get(y_))),
                        by = c("id", by_)]
    agg_dt[[yout_]] = agg_dt$val_
    agg_dt[, c(yout_, "id", by_), with = FALSE]
}

#' plot_binned_aggregates
#'
#' @param agg_dt
#' @param xbins
#' @param ybins
#' @param xrng
#' @param yrng
#' @param val
#' @param bxval
#' @param byval
#' @param facet_
#' @param bin_met
#' @param min_size
#' @param return_data
#'
#' @return
#'
#' @examples
#' agg_dt = aggregate_signals(profile_dt)
#' agg_dt = merge(agg_dt, tsne_dt, by = c("id", "tall_var"))
#' #control resolution in x and y separately
#' plot_raster_summary(agg_dt, facet_ = c("tall_var", "wide_var"),
#'     xbins = 2, ybins = 3)
#' #filter out sparse regions
#' plot_raster_summary(agg_dt, facet_ = c("tall_var", "wide_var"),
#'     xbins = 12, min_size = 0, xrng =  c(-.5, .5))
plot_binned_aggregates = function(agg_dt,
                               xbins = 50,
                               ybins = xbins,
                               xrng = NULL,
                               yrng = NULL,
                               val = "y",
                               bxval = "tx",
                               byval = "ty",
                               facet_ = "wide_var",
                               extra_vars = character(),
                               bin_met = mean,
                               min_size = 1, return_data = FALSE){



    if(is.null(xrng)) xrng = range(agg_dt[[bxval]])
    if(is.null(yrng)) yrng = range(agg_dt[[byval]])
    agg_dt[bxval >= min(xrng) & bxval <= max(xrng) &
               byval >= min(yrng) & byval <= max(yrng)]
    agg_dt[, bx := bin_values(get(bxval), n_bins = xbins, xrng = xrng)]
    agg_dt[, by := bin_values(get(byval), n_bins = ybins, xrng = yrng)]

    bin_dt = agg_dt[, .(y = bin_met(get(val)), N = .N), c(unique(c(facet_, extra_vars, "bx", "by")))]
    # if(min_size > 1){
    #     bin_dt = bin_dt[N >= min_size]
    # }
    bxvc = bin_values_centers(n_bins = xbins, xrng)
    w = diff(bxvc[1:2])
    byvc = bin_values_centers(n_bins = ybins, yrng)
    h = diff(byvc[1:2])
    bin_dt[, tx := bxvc[bx]]
    bin_dt[, ty := byvc[by]]
    if(return_data){
        return(bin_dt)
    }
    ggplot(bin_dt[N >= min_size], aes(x = tx, y = ty, fill = y)) +
        geom_tile(width = w, height = h) + facet_wrap(facet_) +
        scale_fill_viridis_c() +
        coord_cartesian(xlim = xrng, ylim = yrng)
}
