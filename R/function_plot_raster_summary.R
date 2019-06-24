
#' plot_raster_summary of a value
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
#' @export
#'
#' @examples
#' agg_dt = profile_dt[, .(y = mean(y)), .(id, tall_var, wide_var)]
#' agg_dt = merge(agg_dt, tsne_dt, by = c("id", "tall_var"))
#' #control resolution in x and y separately
#' plot_raster_summary(agg_dt, facet_ = c("tall_var", "wide_var"),
#'     xbins = 2, ybins = 3)
#' #filter out sparse regions
#' plot_raster_summary(agg_dt, facet_ = c("tall_var", "wide_var"),
#'     xbins = 2, min_size = 10, xrng =  c(-.5, .5))
plot_raster_summary = function(agg_dt,
                               xbins = 50,
                               ybins = xbins,
                               xrng = NULL,
                               yrng = NULL,
                               val = "y",
                               bxval = "tx",
                               byval = "ty",
                               facet_ = "wide_var",
                               bin_met = mean,
                               min_size = 1, return_data = FALSE){



    if(is.null(xrng)) xrng = range(agg_dt[[bxval]])
    if(is.null(yrng)) yrng = range(agg_dt[[byval]])
    agg_dt[, bx := bin_values(tx, n_bins = xbins, xrng = xrng)]
    agg_dt[, by := bin_values(ty, n_bins = ybins, xrng = yrng)]

    bin_dt = agg_dt[, .(y = bin_met(y), N = .N), c(facet_, "bx", "by")]
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
