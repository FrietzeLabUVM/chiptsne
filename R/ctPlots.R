.prepare_plot_inputs = function(sts, feature_name, signal_name){
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

    prof_dt = sts$signal_data[[feature_name]][[signal_name]]$signal_data
    prof_dt[, tall_var := "none"]
    prof_dt[, wide_var := name]

    tsne_dt = sts$signal_data[[feature_name]][[signal_name]]$xy_data
    return(list(
        feature_name = feature_name,
        signal_name = signal_name,
        prof_dt = prof_dt,
        tsne_dt = tsne_dt
    ))
}

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
#' data(sts.test)
#' ctPlotBinAggregates(sts)
ctPlotBinAggregates = function(sts,
                               feature_name = NULL,
                               signal_name = NULL,
                               xmin = -Inf,
                               xmax = Inf,
                               agg_FUN = max,
                               xbins = 20,
                               ybins = xbins){
    args2 = .prepare_plot_inputs(sts, feature_name, signal_name)
    for(var_name in names(args2)){
        assign(var_name, args2[[var_name]])
    }

    #aggregate_signals is retrieves a single representative value for a genomic region per sample
    #by default this is the maximum anywhere in the region but this can be
    #overriden using xmin/xmax and agg_FUN
    agg_dt = aggregate_signals(
        prof_dt,
        yout_ = "y",
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax)
    agg_dt = merge(agg_dt, tsne_dt, by = "id")

    plot_binned_aggregates(
        agg_dt = agg_dt,
        xbins = xbins,
        ybins = ybins)
}


#' ctPlotSummaryProfiles
#'
#' @param sts
#' @param feature_name
#' @param signal_name
#' @param xbins
#' @param ybins
#' @param q_tall_values
#' @param q_wide_values
#' @param xrng
#' @param yrng
#' @param plot_type
#' @param rname
#' @param odir
#' @param force_rewrite
#' @param n_cores
#' @param ylim
#' @param ma_size
#' @param n_splines
#' @param p
#' @param line_color_mapping
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param return_data
#'
#' @return
#' @export
#'
#' @examples
ctPlotSummaryProfiles = function(sts,
                                 feature_name = NULL,
                                 signal_name = NULL,
                                 xbins = 20,
                                 ybins = xbins,
                                 ###
                                 q_tall_values = NULL,
                                 q_wide_values = NULL,
                                 xrng = NULL,
                                 yrng = NULL,
                                 plot_type = c("glyph", "raster")[1],
                                 rname = NULL,
                                 odir = NULL,
                                 force_rewrite = FALSE,
                                 n_cores = getOption("mc.cores", 1),
                                 ylim = c(0, NA),
                                 ma_size = 2,
                                 n_splines = 10,
                                 p = NULL,
                                 line_color_mapping = NULL,
                                 N_floor = 0,
                                 N_ceiling = NULL,
                                 min_size = 0.3,
                                 return_data = FALSE
){
    args2 = chiptsne:::.prepare_plot_inputs(sts,
                                            feature_name = feature_name,
                                            signal_name = signal_name)
    for(var_name in names(args2)){
        assign(var_name, args2[[var_name]])
    }
    chiptsne:::stsPlotSummaryProfiles(
        profile_dt = prof_dt,
        position_dt = tsne_dt,
        x_points = xbins,
        y_points = ybins,
        ###
        q_tall_values = q_tall_values,
        q_wide_values = q_wide_values,
        xrng = xrng,
        yrng = yrng,
        plot_type = plot_type,
        rname = rname,
        odir = odir,
        force_rewrite = force_rewrite,
        n_cores = n_cores,
        apply_norm = FALSE,
        ylim = ylim,
        ma_size = ma_size,
        n_splines = n_splines,
        p = p,
        facet_byCell = FALSE,
        line_color_mapping = line_color_mapping,
        vertical_facet_mapping = NULL,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_size,
        return_data = return_data
    )
}
