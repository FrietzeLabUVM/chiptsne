testthat::context("ssvTSNE")
library(ssvQC)
library(chiptsne)
library(testthat)
options(mc.cores = 10)
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = TRUE
# SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = FALSE
set.seed(0)
features_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
features_config = QcConfigFeatures.parse(features_config_file)

bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
bam_config = QcConfigSignal.parse(bam_config_file)
bam_config$center_signal_at_max = TRUE
bam_config$flip_signal_mode = "high_on_left"
bam_config$view_size = 600
bam_config$heatmap_limit_values = c(0, 3)

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

debug(chiptsne:::plot_binned_aggregates)

ctPlotBinAggregates(sts)
ctPlotBinAggregates(sts, xbins = 5)
ctPlotBinAggregates(sts, xbins = 5, xmin = -.1, xmax = .1)
ctPlotBinAggregates(sts, xbins = 5, xmin = -Inf, xmax = -.3)

ctPlotBinAggregates(sts, xbins = 5, xmin = .3, xmax = Inf) +
    facet_wrap(~cell+mark)

ctPlotSummaryProfiles(sts)
ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1)
ctPlotSummaryProfiles(sts, plot_type = "raster")
ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1, plot_type = "raster")

ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1) +
    facet_wrap(~cell)

ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1, plot_type = "raster") +
    facet_wrap(~mark)

ctPlotPoints(sts)
ctPlotPoints(sts) +
    facet_wrap(~cell+mark)
ctPlotPoints(sts) +
    facet_wrap(~name_split)

chiptsne:::stsPlotSummaryProfiles(profile_dt = )
chiptsne:::prep_summary
chiptsne:::plot_summary_glyph

#save the test ssvTSNE object
if(FALSE){
    save(sts, file = "data/sts.test.rda")
}
