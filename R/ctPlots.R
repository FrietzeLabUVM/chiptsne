

#' ctPlotBinAggregates
#'
#' Divide t-SNE space into xbins*ybins and summarize signal in each bean as the
#' mean value results from agg_FUN on xmin:xmax of signal profiles present.
#'
#' @param sts A ssvTSNE object with ssvQC.prepSignal already called.
#' @param feature_name Feature name present in sts.  With default of NULL, first
#'   feature name will be used.
#' @param signal_name Signal name present in sts.  With default of NULL, first
#'   signal name will be used.
#' @param xmin The min range of x-values to apply agg_FUN to.
#' @param xmax The max range of x-values to apply agg_FUN to.
#' @param profile_value Value to use for profiles. Default is the plot_value
#'   defined in the signal config of sts.
#' @param profile_value_label Label to use for profile scale. Default is the
#'   plot_value label defined in the signal config of sts.
#' @param agg_FUN A function appled to all y_ values in xmin to xmax range. Must
#'   accept single numeric vector and return 1 value.
#' @param xbins Number of bins (pixel) in the x direction
#' @param ybins Number of bins (pixel) in the y direction. Default reuses xbins.
#' @param bg_color Color to use for plot background. Default is "gray60".
#' @param min_size Bins must contain at least this many points to appear in
#'   final plot.
#'
#' @return ggplot summarizing signal profiles across t-SNE space in heatmap
#'   style plot.
#' @export
#'
#' @examples
#' data(ex_sts)
#' ctPlotBinAggregates(sts)
ctPlotBinAggregates = function(sts,
                               feature_name = NULL,
                               signal_name = NULL,
                               xmin = -Inf,
                               xmax = Inf,
                               profile_value = ssvQC:::val2var[sts$signal_config$plot_value],
                               profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value],
                               agg_FUN = max,
                               xbins = 50,
                               ybins = xbins,
                               bg_color = "gray60",
                               min_size = 5,
                               extra_vars = character()){
    .prepare_plot_inputs(
        sts,
        feature_name,
        signal_name
    )

    #aggregate_signals is retrieves a single representative value for a genomic region per sample
    #by default this is the maximum anywhere in the region but this can be
    #overriden using xmin/xmax and agg_FUN
    agg_dt = aggregate_signals(
        prof_dt,
        y_ = profile_value,
        yout_ = profile_value,
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax)
    agg_dt = merge(agg_dt, tsne_dt, by = "id")
    m_vars = setdiff(colnames(sts$signal_config$meta_data), colnames(agg_dt))
    agg_dt = as.data.table(merge(sts$signal_config$meta_data[, m_vars], agg_dt, by.y = "wide_var", by.x = "name"))
    facet_str = paste0(sts$signal_config$color_by, "~", sts$signal_config$run_by)
#TODO add extra vars
    #TODO add ssvQC win_size
    extra_vars =  c(
        sts$signal_config$run_by,
        sts$signal_config$color_by,
        "name",
        "name_split"
    )
    extra_vars = extra_vars[extra_vars %in% colnames(agg_dt)]
    plot_binned_aggregates(
        agg_dt = agg_dt,
        xbins = xbins,
        ybins = ybins,
        val = profile_value,
        extra_vars = extra_vars,
        facet_ = "name",
        min_size = min_size
    ) +
        labs(fill = profile_value_label) +
        facet_grid(facet_str) +
        theme(
            panel.background = element_rect(fill = bg_color),
            panel.grid = element_blank()
        )
}


#' ctPlotSummaryProfiles
#'
#'
#' @param sts A ssvTSNE object with ssvQC.prepSignal already called.
#' @param feature_name Feature name present in sts.  With default of NULL, first
#'   feature name will be used.
#' @param signal_name Signal name present in sts.  With default of NULL, first
#'   signal name will be used.
#' @param xbins Number of bins (pixel) in the x direction
#' @param ybins Number of bins (pixel) in the y direction. Default reuses xbins.
#' @param profile_value Value to use for profiles. Default is the plot_value
#'   defined in the signal config of sts.
#' @param profile_value_label Label to use for profile scale. Default is the
#'   plot_value label defined in the signal config of sts.
#' @param xrng numeric vector of length 2 defining range of x-axis. Default is full range.
#' @param yrng numeric vector of length 2 defining range of y-axis. Default is full range.
#' @param ylim Range of y-axis within profile glyphs. Default is 0 to max.
#' @param ma_size Number of profile x-values to use for moving average.
#' @param n_splines Number of splines to use to smooth after moving average.
#' @param p Previous ggplot to plot over.
#' @param N_floor Profiles with N_floor points will have the minimum size.
#' @param N_ceiling Profiles with N_ceiling or more points will have the maximum
#'   size.
#' @param min_fraction Minimum fraction between N_floor and N_ceiling required
#'   to be included. For example 1) if N_floor is 5, and min_fraction is 0,
#'   profiles with fewer than 5 points will be omitted 2) if N_ceiling is 25 and
#'   min_fraction is 1, all profiles will be plotted at max size and profiles
#'   with fewer than 25 points will be dropped 3) N_floor is 5 and N_ceiling is
#'   25 and min_fraction is .5, profiles with 15 points will be half size and
#'   profiles will fewer than 15 points will be omitted.
#' @param return_data
#'
#' @return ggplot where profiles present in bins are summarized by glyphs.
#' @export
#'
#' @examples
#' data(ex_sts)
#' ctPlotSummaryProfiles(sts)
ctPlotSummaryProfiles = function(sts,
                                 feature_name = NULL,
                                 signal_name = NULL,
                                 xbins = 10,
                                 ybins = xbins,
                                 profile_value = ssvQC:::val2var[sts$signal_config$plot_value],
                                 profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value],
                                 ###
                                 xrng = NULL,
                                 yrng = NULL,
                                 ylim = c(0, NA),
                                 ma_size = 2,
                                 n_splines = 10,
                                 p = NULL,
                                 N_floor = 0,
                                 N_ceiling = NULL,
                                 min_fraction = 0.2,
                                 return_data = FALSE,
                                 extra_vars = c("name", "name_split")
){
    q_tall_values = NULL
    q_wide_values = NULL
    #these are only used by raster
    plot_type = "glyph" #raster not supported due to facetting.
    force_rewrite = FALSE
    n_cores = getOption("mc.cores", 1)
    rname = NULL
    odir = NULL

    .prepare_plot_inputs(
        sts,
        feature_name = feature_name,
        signal_name = signal_name
    )

    extra_vars =  union(
        extra_vars,
        c(
            sts$signal_config$run_by,
            sts$signal_config$color_by
        )
    )

    extra_vars = extra_vars[extra_vars %in% colnames(prof_dt)]

    p = stsPlotSummaryProfiles(
        profile_dt = prof_dt,
        position_dt = tsne_dt,
        x_points = xbins,
        y_points = ybins,
        y_var = profile_value,
        extra_vars = extra_vars,
        wide_var = sts$signal_config$color_by,
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
        line_color_mapping = unlist(sts$signal_config$color_mapping),
        vertical_facet_mapping = NULL,
        N_floor = N_floor,
        N_ceiling = N_ceiling,
        min_size = min_fraction,
        return_data = return_data
    )
    facet_str = paste0("~", sts$signal_config$run_by)
    p +  facet_wrap(facet_str)
}

#' ctPlotPoints
#'
#'
#' @param sts A ssvTSNE object with ssvQC.prepSignal already called.
#' @param feature_name Feature name present in sts.  With default of NULL, first
#'   feature name will be used.
#' @param signal_name Signal name present in sts.  With default of NULL, first
#'   signal name will be used.
#' @param xmin The min range of x-values to apply agg_FUN to.
#' @param xmax The max range of x-values to apply agg_FUN to.
#' @param profile_value Value to use for profiles. Default is the plot_value
#'   defined in the signal config of sts.
#' @param profile_value_label Label to use for profile scale. Default is the
#'   plot_value label defined in the signal config of sts.
#' @param bg_color Color to use for plot background. Default is "gray60".
#' @param agg_FUN A function appled to all y_ values in xmin to xmax range. Must
#'   accept single numeric vector and return 1 value.
#' @param point_size Size to plot points at. Default is 1.
#'
#' @return ggplot where each profile is summarized by a single point in t-SNE space.
#' @export
#'
#' @examples
#' data(ex_sts)
#' ctPlotPoints(sts)
#'
ctPlotPoints = function(
        sts,
        feature_name = NULL,
        signal_name = NULL,
        xmin = -Inf,
        xmax = Inf,
        profile_value = ssvQC:::val2var[sts$signal_config$plot_value],
        profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value],
        bg_color = "gray60",
        agg_FUN = max,
        point_size = 1
){
    .prepare_plot_inputs(sts, feature_name = feature_name, signal_name = signal_name)
    agg_dt = aggregate_signals(
        prof_dt,
        y_ = profile_value,
        yout_ = profile_value,
        agg_FUN = agg_FUN,
        xmin = xmin,
        xmax = xmax)
    agg_dt = merge(agg_dt, tsne_dt, by = "id")
    m_vars = setdiff(colnames(sts$signal_config$meta_data), colnames(agg_dt))
    agg_dt = as.data.table(merge(sts$signal_config$meta_data[, m_vars], agg_dt, by.y = "wide_var", by.x = "name"))
    facet_str = paste0(sts$signal_config$color_by, "~", sts$signal_config$run_by)

    ggplot(agg_dt, aes_string(x = "tx", y = "ty", color = profile_value)) +
        geom_point(size = point_size) +
        facet_grid(facet_str) +
        scale_color_viridis_c() +
        theme(
            panel.background = element_rect(fill = bg_color),
            panel.grid = element_blank()) +
        labs(color = profile_value_label)
}

#' ctPlotPointsAnnotation
#'
#' @param sts A ssvTSNE object with ssvQC.prepSignal already called.
#' @param feature_name Feature name present in sts.  With default of NULL, first
#'   feature name will be used.
#' @param signal_name Signal name present in sts.  With default of NULL, first
#'   signal name will be used.
#' @param meta_data A data.frame containing "id" and anno_var. Will use cluster assignment if not set.
#' @param anno_var Variable to extract from meta_data for plotting.
#' @param anno_var_label Label to use for color legend. Default is anno_var.
#' @param bg_color Color to use for plot background. Default is "gray60".
#' @param return_data
#'
#' @return ggplot with color of t-SNE points determined by annotation in meta_data.
#' @export
#'
#' @examples
#' data(ex_sts)
#' ctPlotPointsAnnotation(sts)
#'
#' #example using peak call
#' query_gr = sts$features_config$assessment_features[[1]]
#' meta_data = data.table(id = names(query_gr), chromosome = as.character(seqnames(query_gr)))
#' ctPlotPointsAnnotation(sts, meta_data = meta_data, anno_var = "chromosome")
ctPlotPointsAnnotation = function(
        sts,
        meta_data = NULL,
        feature_name = NULL,
        signal_name = NULL,
        anno_var = NULL,
        anno_var_label = anno_var,
        bg_color = "gray60",
        return_data = FALSE,
        point_size = 1
){
    .prepare_plot_inputs(sts, feature_name = feature_name, signal_name = signal_name)
    if(is.null(meta_data)){
        meta_data = sts$signal_data[[feature_name]][[signal_name]]$assignment_data
        anno_var = "cluster_id"
        message("No meta_data given, annotating points with clustering assignment.")
    }
    if(is.null(anno_var)){
        stop("anno_var must be set.")
    }
    if(!anno_var %in% colnames(meta_data)){
        stop("anno_var must be present in meta_data.")
    }
    # if(anno_var %in% colnames(tsne_dt)){
    #     stop(anno_var, " anno_var must not be present in tsne_dt")
    # }
    # sts$signal_data[[feature_name]][[signal_name]]$xy_data
    # ko_cn = setdiff(colnames(meta_data), colnames(tsne_dt))
    # ko_cn = union(ko_cn, "id")

    ko_meta_data = setdiff(colnames(meta_data), c("tx", "ty"))
    ko_tsne_dt = setdiff(colnames(tsne_dt), c(anno_var))

    tsne_dt = merge(meta_data[, ko_meta_data, with = FALSE], tsne_dt[, ko_tsne_dt, with = FALSE], by = "id")
    # facet_str = paste0(sts$signal_config$color_by, "~", sts$signal_config$run_by)

    if(return_data){
        return(tsne_dt)
    }

    ggplot(tsne_dt, aes_string(x = "tx", y = "ty", color = anno_var)) +
        geom_point(size = point_size) +
        # facet_grid(facet_str) +
        # scale_color_viridis_c() +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank()) +
        labs(color = anno_var_label)
}

.prepare_plot_inputs = function(sts, feature_name, signal_name, env = parent.frame()){
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
    args = list(
        feature_name = feature_name,
        signal_name = signal_name,
        prof_dt = prof_dt,
        tsne_dt = tsne_dt
    )

    for(var_name in names(args)){
        assign(var_name, args[[var_name]], pos = env)
    }

    invisible(list(
        feature_name = feature_name,
        signal_name = signal_name,
        prof_dt = prof_dt,
        tsne_dt = tsne_dt
    ))
}


#' ctClusterPoints
#'
#' @param sts A ssvTSNE object with ssvQC.prepSignal already called.
#' @param n_clust number of nearest neighbor clusters to calculate.  Reached to calling too many clusters and iteratively combining most similar clusters.
#'
#' @return A data.table containing id and cluster_id information. Suitable for ctPlotPointsAnnotation.
#' @export
#'
#' @examples
#' data(ex_sts)
#' clust_dt =ctClusterPoints(sts, n_clust = 3)
#'
#' ctPlotPointsAnnotation(sts, meta_data = clust_dt, anno_var = "cluster_id")
ctClusterPoints = function(sts, feature_name = NULL, signal_name = NULL, n_clust = 6){
    .prepare_plot_inputs(sts, feature_name = feature_name, signal_name = signal_name)

    nn_clust.k = function(tsne_dt, tall_var = "tall_var", n_clust = n_clust){
        chiptsne:::nn_clust(tsne_dt,
                            tall_var = tall_var,
                            nn = ceiling(nrow(tsne_dt)/n_clust))
    }

    # tsne_dt.clust = chiptsne:::nn_clust(tsne_dt, tall_var = "tall_none", nn = 500)
    # tsne_dt.clust$cluster_id
    tsne_dt.clust = nn_clust.k(tsne_dt, tall_var = "tall_none", n_clust = n_clust*2)

    n_found = length(unique(tsne_dt.clust$cluster_id))
    if(n_found > n_clust){
        tsne_dt.clust = chiptsne:::combine_most_similar(tsne_dt.clust, prof_dt, n_times = n_found - n_clust)
    }else if(n_found < n_clust){
        message("Too few clusters found by nn_clust to reach requested n_clust.")
    }
    # ggplot(tsne_dt.clust, aes(x = tx, y = ty, color = cluster_id)) +
    #     geom_point()
    tsne_dt.clust[, .(id, cluster_id)]
}

