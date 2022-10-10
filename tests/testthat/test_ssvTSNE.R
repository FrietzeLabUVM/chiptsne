testthat::context("ssvTSNE")
library(ssvQC)
library(chiptsne)
library(testthat)
options(mc.cores = 10)
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = TRUE
set.seed(0)
features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
features_config = QcConfigFeatures.parse(features_config_file)

bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
bam_config = QcConfigSignal.parse(bam_config_file)
bam_config$view_size = 600

bam_config$color_by
bam_config$run_by
bam_config$meta_data

# TSNE plot paramters can be controlled at ssvTSNE object creation or set later
# if the y settings aren't specified they are the same as x
sts = ssvTSNE(
    features_config,
    bam_config,
    n_glyphs_x = 3,
    n_heatmap_pixels_x = 5)
#
#
# sts$perplexity
#
# # sts = ssvQC.plotFeatures(sts)
# # sts$perplexity = 10
# sts = ssvQC.prepFetch(sts)
# # debug(ssvQC:::ClusteredSignal_TSNE.from_ClusteredSignal)
# sts = ssvQC.prepSignal(sts)
# sts$signal_data
# sts = ssvQC.referenceUsesSameScale(sts)
# sts = ssvQC.prepSignal(sts)
# sts = ssvQC.plotSignal(sts)

suppressWarnings({
    sts = ssvQC.plotSignal(sts)
})

sts$plots$signal$heatmaps$CTCF_features$CTCF_signal
sts$plots$TSNE$regional_glyphs$CTCF_features$CTCF_signal

sts$signal_data$CTCF_features$CTCF_signal$signal_data
sts$signal_data$CTCF_features$CTCF_signal$xy_data

test_that("sts TSNE plots", {
    expect_is(sts$plots$TSNE, "list")
    expect_is(sts$plots$TSNE$regional_glyphs$CTCF_features$CTCF_signal, "ggplot")
    expect_is(sts$plots$TSNE$cluster_glyphs$CTCF_features$CTCF_signal, "ggplot")
    expect_is(sts$plots$TSNE$regional_heatmap$CTCF_features$CTCF_signal, "ggplot")
})

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

sts

stsPlotBinAggregates(sts, xbins = 5)

prof_dt = sts$signal_data$CTCF_features$CTCF_signal$signal_data
prof_dt[, tall_var := "none"]
prof_dt[, wide_var := name]
tsne_dt = sts$signal_data$CTCF_features$CTCF_signal$xy_data
agg_dt = chiptsne::aggregate_signals(prof_dt)
agg_dt = merge(agg_dt, tsne_dt, by = "id")
# undebug(chiptsne::plot_binned_aggregates)
chiptsne::plot_binned_aggregates(agg_dt = agg_dt, val = "ynorm", xbins = 3)

#save the test ssvTSNE object
if(FALSE){
    save(sts, file = "data/sts.test.rda")
}
