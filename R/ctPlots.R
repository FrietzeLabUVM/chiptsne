
#' ctPlotBinAggregates
#'
#' @param sts
#' @param feature_name
#' @param signal_name
#' @param xmin
#' @param xmax
#' @param agg_FUN
#' @param xbins
#' @param ybins
#'
#' @return
#' @export
#'
#' @examples
ctPlotBinAggregates = function(sts,
                               feature_name = NULL,
                               signal_name = NULL,
                               xmin = -Inf,
                               xmax = Inf,
                               agg_FUN = max,
                               xbins = 20,
                               ybins = xbins){
    if(is.null(feature_name)){
        feature_name = names(sts$signal_data)[1]
    }
    if(!feature_name %in% names(sts$signal_data)){
        stop(feature_name, " not a valid feature name. Valid: ",
             paste(names(sts$signal_data), collapse = ", "))
    }
    if(is.null(signal_name)){
        signal_name = names(sts$signal_data[[feature_name]])[1]
    }
    if(!signal_name %in% names(sts$signal_data[[feature_name]])){
        stop(signal_name, " not a valid signal name. Valid: ",
             paste(names(sts$signal_data[[feature_name]]), collapse = ", "))
    }

    prof_dt = sts$signal_data$CTCF_features$CTCF_signal$signal_data
    prof_dt[, tall_var := "none"]
    prof_dt[, wide_var := name]

    tsne_dt = sts$signal_data$CTCF_features$CTCF_signal$xy_data

    #aggregate_signals is retrieves a single representative value for a genomic region per sample
    #by default this is the maximum anywhere in the region but this can be
    #overriden using xmin/xmax and agg_FUN
    agg_dt = chiptsne::aggregate_signals(
        prof_dt,
        yout_ = "y",
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax)
    agg_dt = merge(agg_dt, tsne_dt, by = "id")

    chiptsne::plot_binned_aggregates(
        agg_dt = agg_dt,
        xbins = xbins,
        ybins = ybins)
}
