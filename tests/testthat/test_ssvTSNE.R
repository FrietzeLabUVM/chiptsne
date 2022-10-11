testthat::context("ssvTSNE")
library(ssvQC)
library(chiptsne)
library(testthat)
options(mc.cores = 10)
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = TRUE
SQC_OPTIONS$SQC_FORCE_CACHE_OVERWRITE = FALSE
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


ctPlotBinAggregates(sts)
ctPlotBinAggregates(sts, xbins = 5)
ctPlotBinAggregates(sts, xbins = 5, xmin = -.1, xmax = .1)
ctPlotBinAggregates(sts, xbins = 5, xmin = -Inf, xmax = -.3)
ctPlotBinAggregates(sts, xbins = 5, xmin = .3, xmax = Inf)

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

undebug(ctPlotSummaryProfiles)
debug(chiptsne:::stsPlotSummaryProfiles)

ctPlotSummaryProfiles(sts)
ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1)
ctPlotSummaryProfiles(sts, plot_type = "raster")
ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 1, plot_type = "raster")

chiptsne:::stsPlotSummaryProfiles(profile_dt = )
chiptsne:::prep_summary
chiptsne:::plot_summary_glyph

#save the test ssvTSNE object
if(FALSE){
    save(sts, file = "data/sts.test.rda")
}
