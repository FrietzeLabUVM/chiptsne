#' @importClassesFrom ssvQC ssvQC.complete
#' @importMethodsFrom ssvQC ssvQC.runAll ssvQC.prepSignal ssvQC.plotSignal
setClass("ChIPtSNE",
         representation = list(
             perplexity = "numeric",
             n_glyphs_x = "numeric",
             n_glyphs_y = "numeric",
             n_heatmap_pixels_x = "numeric",
             n_heatmap_pixels_y = "numeric"
         ),
         contains = "ssvQC.complete")

#' ChIPtSNE
#'
#' @param features_config Controls features configuration.  May be a:
#'   QcConfigFeatures object, path to a file defining configuration via
#'   QcConfigFeatures.parse, features files to define via
#'   QcConfigFeatures.files, or a data.frame to pass to QcConfigFeatures.
#' @param signal_config Controls signal configuration.  May be a: QcConfigSignal
#'   object, path to a file defining configuration via QcConfigSignal.parse,
#'   features files to define via QcConfigSignal.files, or a data.frame to pass
#'   to QcConfigSignal.
#' @param out_dir NYI
#' @param bfc BiocFileCache object to use for caching. If NULL, default
#'   new_cache() will be used.
#' @param perplexity passed to Rtsne::Rtsne()
#' @param n_glyphs_x number of glyphs across x-axis to summarize profiles in t-SNE space. Default is 8.
#' @param n_glyphs_y number of glyphs across y-axis. Default is the same as n_glyphs_x.
#' @param n_heatmap_pixels_x number of pixels to use in heatmap across x-axis in t-SNE space. Default is 25.
#' @param n_heatmap_pixels_y number of pixels to use in heatmap across y-axis in t-SNE space. Default is same as n_heatmap_pixels_x.
#'
#' @return a valid ChIPtSNE object.
#' @export
#'
#' @examples
#' library(chiptsne)
#' library(ssvQC)
#' query_bed_file = system.file(package = "chiptsne", "extdata/query_gr.bed")
#' query_gr = rtracklayer::import.bed(query_bed_file)
#' bam_files = dir(system.file(package = "chiptsne", "extdata"), pattern = "bam$", full.names = TRUE)
#' qc_config_features = QcConfigFeatures.GRanges(query_gr)
#' qc_config_signal = QcConfigSignal.files(bam_files)
#' qc_config_signal$flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left
#' qc_config_signal$center_signal_at_max = TRUE
#' sts = ChIPtSNE(qc_config_features, qc_config_signal)
#' sts = ChIPtSNE.runTSNE(sts)
#'
#' sts$n_glyphs_x = 7
#' sts$n_glyphs_y = 5
#'
#' sts$n_heatmap_pixels_x = 5
#' sts$n_heatmap_pixels_y = 9
#' sts.replot = ssvQC.plotSignal(sts)
#'
#' sts$plots$TSNE$regional_glyphs$query_features$All_signal
#' sts.replot$plots$TSNE$regional_glyphs$query_features$All_signal
#'
#' sts$plots$TSNE$regional_heatmap$query_features$All_signal
#' sts.replot$plots$TSNE$regional_heatmap$query_features$All_signal
#'
ChIPtSNE = function(features_config = NULL,
                   signal_config = NULL,
                   out_dir = getwd(),
                   bfc = NULL,
                   perplexity = 50,
                   n_glyphs_x = 8,
                   n_glyphs_y = n_glyphs_x,
                   n_heatmap_pixels_x = 25,
                   n_heatmap_pixels_y = n_heatmap_pixels_x){
    matched_only = FALSE
    if(is.null(features_config) & is.null(signal_config)){
        stop("At least one of features_config or signal_config must be specified.")
    }

    features_config = ssvQC:::.prep_features_config(features_config)
    signal_config = ssvQC:::.prep_signal_config(signal_config)

    signal_config@cluster_value = SQC_SIGNAL_VALUES$linearQuantile
    signal_config@sort_value = SQC_SIGNAL_VALUES$linearQuantile
    signal_config@plot_value = SQC_SIGNAL_VALUES$linearQuantile

    if(is.null(bfc)){
        bfc = ssvQC:::new_cache()
    }

    dir.create(out_dir, showWarnings = FALSE)

    new("ChIPtSNE",
        features_config = features_config,
        signal_config = signal_config,
        signal_data = list(),
        other_data = list(),
        out_dir = out_dir,
        bfc = bfc,
        saving_enabled = TRUE,
        perplexity = perplexity,
        n_glyphs_x = n_glyphs_x,
        n_glyphs_y = n_glyphs_y,
        n_heatmap_pixels_x = n_heatmap_pixels_x,
        n_heatmap_pixels_y = n_heatmap_pixels_y,
        matched_only = matched_only)
}

# ssvQC.runAll = ssvQC::ssvQC.runAll
setMethod("ssvQC.runAll", "ChIPtSNE", function(object){
    object = callNextMethod()
    message("NYI")
    object
})

# ssvQC.prepSignal = ssvQC::ssvQC.prepSignal
setMethod("ssvQC.prepSignal", "ChIPtSNE", function(object){
    object = callNextMethod()
    object@signal_data = lapply(object@signal_data, function(signal_clust_objs){
        lapply(signal_clust_objs, function(clust_obj){
            ClusteredSignal_TSNE.from_ClusteredSignal(clust_obj, object)
        })
    })
    object
})


# ssvQC.plotSignal = ssvQC::ssvQC.plotSignal
setMethod("ssvQC.plotSignal", "ChIPtSNE", function(object){
    object = callNextMethod()

    n_glyphs_x = object@n_glyphs_x
    n_glyphs_y = object@n_glyphs_y

    n_heatmap_pixels_x = object@n_heatmap_pixels_x
    n_heatmap_pixels_y = object@n_heatmap_pixels_y

    object@plots$TSNE$regional_glyphs = lapply(object@signal_data, function(signal_data_groups){
        lapply(signal_data_groups, function(signal_data){
            if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
                stop("Stored signal data is not of class ClusteredSignal_TSNE.")
            }
            extra_vars =  c(
                sts$signal_config$run_by,
                sts$signal_config$color_by,
                "name",
                "name_split"
            )
            extra_vars = extra_vars[extra_vars %in% colnames(signal_data@signal_data)]

            p_summary_profiles = stsPlotSummaryProfiles(
                signal_data@signal_data,
                signal_data@xy_data,
                x_points = n_glyphs_x,
                y_points = n_glyphs_y,
                y_var = ssvQC:::val2var[object@signal_config@plot_value],
                extra_vars = extra_vars,
                wide_var = sts$signal_config$color_by,
                line_color_mapping = unlist(sts$signal_config$color_mapping)
            )
            p_summary_profiles
        })
    })

    object@plots$TSNE$cluster_glyphs = lapply(object@signal_data, function(signal_data_groups){
        lapply(signal_data_groups, function(signal_data){
            if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
                stop("Stored signal data is not of class ClusteredSignal_TSNE.")
            }

            use_tsne_clust = TRUE
            if(use_tsne_clust){
                cluster_result = nn_clust(signal_data@xy_data, tall_var = "tall_none")
                p_cluster_profiles = stsPlotClusterProfiles(
                    signal_data@signal_data,
                    cluster_result,
                    wide_var = object$signal_config$color_by)
            }else{
                p_cluster_profiles = stsPlotClusterProfiles(
                    signal_data@signal_data,
                    merge(
                        signal_data@signal_data,
                        signal_data@xy_data, by = "id"),
                    wide_var = object$signal_config$color_by)
            }
            p_cluster_profiles + scale_color_manual(values = unlist(object@signal_config$color_mapping))
        })
    })

    object@plots$TSNE$regional_heatmap = lapply(object@signal_data, function(signal_data_groups){
        lapply(signal_data_groups, function(signal_data){
            if(!"ClusteredSignal_TSNE" %in% class(signal_data)){
                stop("Stored signal data is not of class ClusteredSignal_TSNE.")
            }
            prof_dt = signal_data@signal_data
            xy_dt = signal_data@xy_data

            x_var = "x"
            y_var = "y"
            id_var = "id"
            wide_var = object$signal_config$color_by
            tall_var = "tall_none"
            agg_FUN = max

            y_var = ssvQC:::val2var[object@signal_config@plot_value]
            y_lab = names(y_var)

            xy_dt[, xbin := bin_values(tx, n_heatmap_pixels_x, c(-.5, .5))]
            xy_dt[, ybin := bin_values(ty, n_heatmap_pixels_y, c(-.5, .5))]

            heat_dt = merge(prof_dt, xy_dt[, .(id, xbin, ybin)], by = "id")[, .(y = agg_FUN(get(y_var))), c(id_var, wide_var, "xbin", "ybin")]

            facet_str = paste0("~", wide_var)
            ggplot(heat_dt, aes(x = xbin, y = ybin, fill = y)) +
                geom_tile() +
                facet_wrap(facet_str) +
                scale_fill_viridis_c()
        })
    })


    saveRDS(object, "dev_object.ssvQC.plotSignal.Rds")
    object
})


### $ Accessor
setMethod("names", "ChIPtSNE",
          function(x)
          {
              c("plots", "signal_data", "signal_config", "features_config", "SCC", "FRIP", "correlation",
                "perplexity",
                "n_glyphs_x",
                "n_glyphs_y",
                "n_heatmap_pixels_x",
                "n_heatmap_pixels_y")

          })


setMethod("$", "ChIPtSNE",
          function(x, name)
          {
              switch (name,
                      plots = x@plots,
                      signal_data = x@signal_data,
                      SCC = x@other_data$SCC,
                      FRIP = x@other_data$FRIP,
                      correlation = list(read_count = x@other_data$read_count_correlation, signal_profile = x@other_data$signal_profile_correlation),
                      bfc = x@bfc,
                      features_config = x@features_config,
                      signal_config = x@signal_config,
                      perplexity = x@perplexity,
                      n_glyphs_x = x@n_glyphs_x,
                      n_glyphs_y = x@n_glyphs_y,
                      n_heatmap_pixels_x = x@n_heatmap_pixels_x,
                      n_heatmap_pixels_y = x@n_heatmap_pixels_y

              )
          })

setReplaceMethod("$", "ChIPtSNE",
                 function(x, name, value)
                 {
                     warn_msg = "This assignment is not supported.  No effect."
                     switch (name,
                             features_config = {
                                 x@features_config = value
                             },
                             signal_config = {
                                 x@signal_config = value
                             },
                             perplexity = {
                                 x@perplexity = value
                             },
                             n_glyphs_x = {
                                 x@n_glyphs_x = value
                             },
                             n_glyphs_y = {
                                 x@n_glyphs_y = value
                             },
                             n_heatmap_pixels_x = {
                                 x@n_heatmap_pixels_x = value
                             },
                             n_heatmap_pixels_y = {
                                 x@n_heatmap_pixels_y = value
                             },
                             {warning(warn_msg)}

                     )

                     #TODO, some assignments may be appropriate
                     x
                 })

#' ChIPtSNE.runTSNE
#'
#' @param sts A valid ChIPtSNE object, needs ssvQC.prepSignal to be called.
#'
#' @return A valid ChIPtSNE object with t-SNE called and ready for plots.
#' @export
#'
#' @examples
#' library(chiptsne)
#' library(ssvQC)
#' query_bed_file = system.file(package = "chiptsne", "extdata/query_gr.bed")
#' query_gr = rtracklayer::import.bed(query_bed_file)
#' bam_files = dir(system.file(package = "chiptsne", "extdata"), pattern = "bam$", full.names = TRUE)
#' qc_config_features = QcConfigFeatures.GRanges(query_gr)
#' qc_config_signal = QcConfigSignal.files(bam_files)
#' qc_config_signal$flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left
#' qc_config_signal$center_signal_at_max = TRUE
#' sts = ChIPtSNE(qc_config_features, qc_config_signal)
#' sts = ChIPtSNE.runTSNE(sts)
ChIPtSNE.runTSNE = function(sts){
    sts = ssvQC.plotSignal(sts)
}

#### QcConfigFeatures ####

#' @export
#' @inherit ssvQC::QcConfigFeatures
ConfigFeatures = ssvQC::QcConfigFeatures

#' @export
#' @inherit ssvQC::QcConfigFeatures.GRanges
ConfigFeatures.GRanges = ssvQC::QcConfigFeatures.GRanges

#' @export
#' @inherit ssvQC::QcConfigFeatures.files
ConfigFeatures.files = ssvQC::QcConfigFeatures.files

#' @export
#' @inherit ssvQC::QcConfigFeatures.parse
ConfigFeatures.parse = ssvQC::QcConfigFeatures.parse

#' @export
#' @inherit ssvQC::QcConfigFeatures.save_config
ConfigFeatures.save_config = ssvQC::QcConfigFeatures.save_config

#### QcConfigSignal ####

#' @export
#' @inherit ssvQC::QcConfigSignal
ConfigSignal = ssvQC::QcConfigSignal

#' @export
#' @inherit ssvQC::QcConfigSignal.files
ConfigSignal.files = ssvQC::QcConfigSignal.files

#' @export
#' @inherit ssvQC::sampleCap
sampleCap = ssvQC::sampleCap

