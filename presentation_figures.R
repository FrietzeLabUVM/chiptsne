#### setup ####
library(chiptsne)
library(seqsetvis)
library(ggplot2)
library(magrittr)
library(GenomicRanges)
library(ssvQC)
library(ssvTracks)
options(mc.cores = 20)

base_point_size = .3
my_cluster_value = SQC_SIGNAL_VALUES$linearQuantile
my_plot_value = SQC_SIGNAL_VALUES$linearQuantile
# treat_cols = c("GFP" = "black", "IK" = "chartreuse2")
treat_cols = c("GFP" = "#5aae61", "IK" = "#762a83")

odir = "output_presentation_figures"
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
res_file = function(f){
    file.path(odir, f)
}
my_ggsave = function(filename, plot = last_plot(), width = 6, height = 6){
    filename = sub(".pdf$", "", filename)
    filename = sub(".png$", "", filename)
    ggsave(filename = res_file(paste0(filename, ".pdf")), plot = plot, width = width, height = height)
    ggsave(filename = res_file(paste0(filename, ".png")), plot = plot, width = width, height = height, dpi = "print")
}

# pdf("presentation_figures.pdf")


chip_bam_cfg_dt = IKdata::chipseq.setup_bam_files()
chip_np_cfg_dt = IKdata::chipseq.setup_peak_files()

col_map = IKdata:::mark_colors
col_map["H3K27ac"] = safeBrew(8, "set1")[2]
col_map["H3K4me1"] = safeBrew(8, "set1")[4]
col_map["input"] = "gray80"
    theme_update(panel.background = element_blank(), panel.grid = element_blank())

    get_bed_files = function(d){
        bed_files = d %>%
            dir(., full.names = TRUE) %>%
            dir(., pattern = "csv$", full.names = TRUE)
        if(length(bed_files) < 1){
            message("skipping no csaw results found.")
            next
        }
        comparisons = sub(".csaw_regions.csv", "", basename(bed_files))
        marks = basename(dirname(bed_files))
        names(bed_files) = paste(marks, comparisons, sep = ".")

        bed_files
    }
    bed_files.chipseq = get_bed_files("figures_ALL/output_csaw.chipseq/")

    #### 1 - single k27ac rep ####
    sel_marks = c("H3K27ac", "input")
    bam_cfg_dt.sel = chip_bam_cfg_dt[mark %in% sel_marks & treatment == "GFP" & rep == "rep1"]
    bam_cfg_dt.sel$All = "all"
    np_cfg_dt.sel = chip_np_cfg_dt[mark %in% sel_marks & treatment == "GFP" & rep == "rep1"]
    col_map.sel = col_map[sel_marks]

    np_grs = easyLoad_narrowPeak(np_cfg_dt.sel$file, file_names = np_cfg_dt.sel$name)
    # olap_grs = ssvOverlapIntervalSets(np_grs)
    # ssvFeatureUpset(olap_grs)
    set.seed(0)
    query_gr = sampleCap(np_grs[[1]], 1e4)
    query_gr = resize(query_gr, width = 3e3, fix = "center")

    qc_config_feat = ConfigFeatures.GRanges(query_gr, n_peaks = Inf)
    qc_config_signal = ConfigSignal(bam_cfg_dt.sel,
                                    heatmap_limit_values = c(0, 3),
                                    flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left,
                                    center_signal_at_max = TRUE,
                                    plot_value = my_plot_value,
                                    cluster_value = my_cluster_value,
                                    sort_value = SQC_SIGNAL_VALUES$RPM,
                                    run_by = "All",
                                    color_by = "mark",
                                    color_mapping = col_map.sel)

    sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
    set.seed(0)
    sts = ChIPtSNE.runTSNE(sts)
    sts$signal_config$plot_value = my_plot_value
    sts$signal_config$cluster_value = my_cluster_value
    sts$signal_config$sort_value = SQC_SIGNAL_VALUES$RPM


    p = sts$plots$signal$heatmaps$query_features$all_signal
    my_ggsave("1_peak_eval_heatmap", p)

    t_dt = sts$signal_data$query_features$all_signal$xy_data
    p = ggplot(t_dt, aes(x = tx, y = ty)) +
        geom_point(size = base_point_size) +
        coord_fixed()
    p
    my_ggsave("1_peak_eval_tsne_pointsBW", p)

    # ctPlotBinAggregates(sts, xmin = -100, xmax = 100) +
    #     facet_grid(~mark) +
    #     coord_fixed()
    # undebug(ctPlotPoints)
    obj_dir = "presentation_ChIPtSNE_objects"
    obj_file = function(f){
        out_f = file.path(obj_dir, f)
        dir.create(dirname(out_f), showWarnings = FALSE, recursive = TRUE)
        out_f
    }
    saveRDS(sts, obj_file("1_peak_eval.Rds"))
    p1 = ctPlotPoints(sts, xmin = -100, xmax = 100, point_size = base_point_size/2, profile_value = "y_RPM") +
        facet_wrap(~mark) +
        coord_fixed() +
        scale_color_viridis_c(limits = c(0, 15), na.value = "yellow") +
        labs(color = "max\nRPM\npileup")
    p1
    my_ggsave("1_peak_eval_tsne_points", p1)

    p2 = ctPlotSummaryProfiles(sts,
                               N_floor = 0,
                               N_ceiling = 150,
                               min_fraction = .25) +
        facet_null() +
        coord_fixed()
    p2
    my_ggsave("1_peak_eval_tsne_summary_profiles", p2)

    cowplot::draw_grob(p1)
    ggdraw(xlim = c(-.5, .5), ylim = c(-.5, .5)) +
        cowplot::draw_plot(p1, x = -.5, y = -.5, width = 1, height = 1) +
        cowplot::draw_plot(p2, x = -.5, y = -.5, width = 1, height = 1)


    #### 2 - k27ac rep comparison ####
    sel_marks = c("H3K27ac")
    bam_cfg_dt.sel = chip_bam_cfg_dt[mark %in% sel_marks & treatment == "GFP" ]
    bam_cfg_dt.sel$All = "all"
    np_cfg_dt.sel = chip_np_cfg_dt[mark %in% sel_marks & treatment == "GFP" ]
    # col_map.sel = col_map[sel_marks]
    col_map.sel = c(rep2 = "#377EB8", rep1 = "darkorange")
    plot(1:2, col = col_map.sel, pch = 16, cex = 4)


    np_grs = easyLoad_narrowPeak(np_cfg_dt.sel$file, file_names = np_cfg_dt.sel$name)
    olap_grs = ssvOverlapIntervalSets(np_grs)
    p = ssvFeatureUpset(olap_grs)
    my_ggsave("2_rep_eval_upset", p)
    set.seed(0)
    query_gr = sampleCap(olap_grs, 1e4)
    query_gr = resize(query_gr, width = 3e3, fix = "center")

    qc_config_feat = ConfigFeatures.GRanges(query_gr, n_peaks = Inf)
    qc_config_signal = ConfigSignal(bam_cfg_dt.sel,
                                    heatmap_limit_values = c(0, 3),
                                    flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left,
                                    center_signal_at_max = TRUE,
                                    plot_value = my_plot_value,
                                    cluster_value = my_cluster_value,
                                    sort_value = SQC_SIGNAL_VALUES$RPM,
                                    run_by = "All",
                                    color_by = "rep",
                                    color_mapping = col_map.sel)

    sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
    set.seed(0)
    sts = ChIPtSNE.runTSNE(sts)
    saveRDS(sts, obj_file("2_rep_eval.Rds"))

    t_dt = sts$signal_data$query_features$all_signal$xy_data
    p = ggplot(t_dt, aes(x = tx, y = ty)) +
        geom_point(size = base_point_size) +
        coord_fixed()
    my_ggsave("2_rep_eval_tsne_pointsBW", p)

    # ctPlotBinAggregates(sts, xmin = -100, xmax = 100) +
    #     facet_grid(~mark) +
    #     coord_fixed()

    p = ctPlotPoints(sts, xmin = -100, xmax = 100, point_size = base_point_size) +
        # facet_wrap(~mark) +
        coord_fixed()
    my_ggsave("2_rep_eval_tsne_points", p)

    peak_memb = ssvFactorizeMembTable(sts$features_config$assessment_features$query_features)
    head(peak_memb)

    peak_memb$group = factor(peak_memb$group, levels = levels(peak_memb$group)[c(1, 3, 2)])

    rep_colors = safeBrew(peak_memb$group)
    p = ctPlotPointsAnnotation(sts, as.data.table(peak_memb), anno_var = "group", point_size = base_point_size / 2) +
        facet_wrap(~group) +
        coord_fixed() +
        scale_color_manual(values = rep_colors) +
        theme(legend.position = "bottom")
    my_ggsave("2_rep_eval_tsne_points_peak_annotation", p)

    rep_colors.simple = rep_colors[-1]
    names(rep_colors.simple) = sub(".+rep", "rep", names(rep_colors.simple))
    p = ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 150, min_fraction = .25) +
        facet_null() +
        coord_fixed() +
        labs(title = "H3K27ac replicate analysis") +
        scale_color_manual(values = rep_colors.simple)
    my_ggsave("2_rep_eval_tsne_summary_profiles", p)

    #### 3 - k27ac differential ####
    sel_marks = c("H3K27ac")
    bam_cfg_dt.sel = chip_bam_cfg_dt[mark %in% sel_marks]
    bam_cfg_dt.sel$All = "all"
    np_cfg_dt.sel = chip_np_cfg_dt[mark %in% sel_marks]
    col_map.sel = IKdata:::treatment_colors
    plot(1:2, col = col_map.sel, pch = 16, cex = 4)

    np_grs = easyLoad_narrowPeak(np_cfg_dt.sel$file, file_names = np_cfg_dt.sel$name)
    olap_grs = ssvOverlapIntervalSets(np_grs)
    p = ssvFeatureUpset(olap_grs)
    my_ggsave("3_diff_eval_peak_upset", p)

    set.seed(0)
    query_gr = sampleCap(olap_grs, 1e4)
    query_gr = resize(query_gr, width = 3e3, fix = "center")

    qc_config_feat = ConfigFeatures.GRanges(query_gr, n_peaks = Inf)
    qc_config_signal = ConfigSignal(bam_cfg_dt.sel,
                                    heatmap_limit_values = c(0, 3),
                                    flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left,
                                    center_signal_at_max = TRUE,
                                    plot_value = my_plot_value,
                                    cluster_value = my_cluster_value,
                                    sort_value = SQC_SIGNAL_VALUES$RPM,
                                    run_by = "All",
                                    color_by = "treatment",
                                    color_mapping = treat_cols)

    sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
    set.seed(0)
    sts = ChIPtSNE.runTSNE(sts)

    saveRDS(sts, obj_file("3_diff_eval.Rds"))

    t_dt = sts$signal_data$query_features$all_signal$xy_data
    p = ggplot(t_dt, aes(x = tx, y = ty)) +
        geom_point(size = base_point_size) +
        coord_fixed()
    my_ggsave("3_diff_eval_tsne_pointsBW", p)

    p = ctPlotPoints(sts, xmin = -100, xmax = 100, point_size = base_point_size) +
        # facet_wrap(~mark) +
        coord_fixed()
    my_ggsave("3_diff_eval_tsne_points", p)

    csaw_grs = pbmcapply::pbmclapply(bed_files.chipseq, function(f){
        dt = fread(f)
        GRanges(dt[FDR < .1 & width > 100])
    })
    diff_gr = csaw_grs$H3K27ac.IK_vs_GFP_in_MXP5
    split(diff_gr, diff_gr$direction)

    p_gr = sts$signal_data$query_features$all_signal$query_gr
    olaps = findOverlaps(query = resize(p_gr, 50, fix = "center"), subject = diff_gr)

    p_gr$direction = "ns"
    p_gr$direction[queryHits(olaps)] = diff_gr$direction[subjectHits(olaps)]
    p_df = as.data.table(p_gr)
    p_df$id = names(p_gr)
    table(p_df$direction)

    p = ctPlotPointsAnnotation(sts, meta_data = p_df, anno_var = "direction", point_size = base_point_size) +
        scale_color_manual(values = c(ns = "gray80", down = "blue", up = "red")) +
        coord_fixed()
    my_ggsave("3_diff_eval_tsne_points_direction", p)

    p_df = subset(p_df, direction != "ns")
    xy_dt = merge(p_df, sts$signal_data$query_features$all_signal$xy_data, by = "id")

    p_summary = ctPlotSummaryProfiles(sts, N_floor = 0, N_ceiling = 150, min_fraction = .25) +
        facet_null() +
        coord_fixed() +
        labs(title = "H3K27ac IK induction")
    my_ggsave("3_diff_eval_tsne_summary_profiles", p_summary)

    p_summary.diff = p_summary +
        annotate("point",
                 x = xy_dt[direction == "up"]$tx,
                 y = xy_dt[direction == "up"]$ty,
                 alpha = .05,
                 color = 'red') +
        annotate("point",
                 x = xy_dt[direction == "down"]$tx,
                 y = xy_dt[direction == "down"]$ty,
                 alpha = .05,
                 color = "blue")
    my_ggsave("3_diff_eval_tsne_summary_profiles_with_diff_points", p_summary.diff)


    #### 3 - dev diff plot ####
    .prepare_plot_inputs = chiptsne:::.prepare_plot_inputs
    aggregate_signals = chiptsne:::aggregate_signals
    feature_name = NULL
    signal_name = NULL
    .prepare_plot_inputs(sts, feature_name = feature_name, signal_name = signal_name)
    xmin = -200
    xmax = 200
    agg_FUN = mean
    profile_value = ssvQC:::val2var[sts$signal_config$plot_value]
    profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value]

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

    # ?set
    diff_dt = agg_dt[, .(y = mean(get(profile_value))), .(group, id, tx, ty)]
    diff_dt = dcast(diff_dt, id+tx+ty~group, value.var = "y")
    diff_dt$diff = diff_dt$`MXP5 H3K27ac IK` - diff_dt$`MXP5 H3K27ac GFP`

    lim = max(abs(diff_dt$diff))
    p = ggplot(diff_dt, aes(x = tx, y = ty, color = diff)) +
        geom_point(size = base_point_size) +
        scale_color_gradientn(colours = c("blue", "gray90", "red"), limits = c(-lim, lim)) +
        coord_fixed() +
        labs(color = "IK - GFP")
    my_ggsave("3_diff_eval_tsne_points_FC", p_summary.diff)

    #### 4 - k27ac with k4me1, Ikaros, ATAC ####
    sel_marks = c("H3K27ac", "H3K4me1",  "ATAC")
    atac_bam_cfg_dt = IKdata::atac.setup_bam_files()[cell == "MXP5"]
    atac_bam_cfg_dt$mark = "ATAC"
    bam_cfg_dt.sel = chip_bam_cfg_dt[mark %in% sel_marks]

    cn = intersect(colnames(bam_cfg_dt.sel), colnames(atac_bam_cfg_dt))
    bam_cfg_dt.sel = rbind(bam_cfg_dt.sel[, cn, with = FALSE], atac_bam_cfg_dt[, cn, with = FALSE])

    bam_cfg_dt.sel$All = "all"
    np_cfg_dt.sel = chip_np_cfg_dt[mark %in% sel_marks]
    col_map.sel = IKdata:::treatment_colors
    plot(1:2, col = col_map.sel, pch = 16, cex = 4)

    atac_dar = IKdata::atac.DAR_results_load()
    atac_dar = subset(atac_dar, DA_group != "ns")


    csaw_grs = pbmcapply::pbmclapply(bed_files.chipseq, function(f){
        dt = fread(f)
        GRanges(dt[FDR < .1 & width > 100])
    })
    diff_grs = c(csaw_grs, list(ATAC = GRanges(atac_dar)))
    diff_olap_grs = ssvOverlapIntervalSets(diff_grs)

    np_grs = easyLoad_narrowPeak(np_cfg_dt.sel$file, file_names = np_cfg_dt.sel$name)
    olap_grs = ssvConsensusIntervalSets(np_grs, min_number = 2, min_fraction = 0)
    p = ssvFeatureUpset(olap_grs)
    my_ggsave("4_combinatorial_upset", p)

    olap_grs.no_diff = subsetByOverlaps(olap_grs, diff_olap_grs, invert = TRUE)

    set.seed(0)
    olap_grs.enriched = c(
        sampleCap(olap_grs.no_diff, length(diff_olap_grs)),
        diff_olap_grs
    )
    names(olap_grs.enriched) = NULL
    mcols(olap_grs.enriched) = NULL
    olap_grs.enriched = prepare_fetch_GRanges_names(olap_grs.enriched)

    set.seed(0)
    query_gr = sampleCap(olap_grs.enriched, 1e4)
    query_gr = resize(query_gr, width = 3e3, fix = "center")

    qc_config_feat = ConfigFeatures.GRanges(query_gr, n_peaks = Inf)

    bam_cfg_dt.sel[, group := paste(mark, treatment)]
    col_map.sel = safeBrew(bam_cfg_dt.sel$group)
    qc_config_signal = ConfigSignal(bam_cfg_dt.sel,
                                    heatmap_limit_values = c(0, 3),
                                    flip_signal_mode = SQC_FLIP_SIGNAL_MODES$high_on_left,
                                    center_signal_at_max = TRUE,
                                    plot_value = my_plot_value,
                                    cluster_value = my_cluster_value,
                                    sort_value = SQC_SIGNAL_VALUES$RPM,
                                    run_by = "All",
                                    color_by = "group",
                                    color_mapping = col_map.sel)

    sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
    set.seed(0)
    sts = ChIPtSNE.runTSNE(sts)
    sts$signal_config$plot_value = my_plot_value
    sts$signal_config$cluster_value = my_cluster_value
    sts$signal_config$sort_value = SQC_SIGNAL_VALUES$linearQuantile

    saveRDS(sts, obj_file("4_combinatorial.Rds"))

    p = ctPlotPoints(sts, point_size = base_point_size / 5) +
        facet_grid(mark~treatment) +
        coord_fixed() +
        scale_color_viridis_c(limits = c(0, 3), na.value = "yellow")
    my_ggsave("4_combinatorial_tsne_points", p)

    #### 4 - diff ####
    .prepare_plot_inputs = chiptsne:::.prepare_plot_inputs
    aggregate_signals = chiptsne:::aggregate_signals
    feature_name = NULL
    signal_name = NULL
    .prepare_plot_inputs(sts, feature_name = feature_name, signal_name = signal_name)
    xmin = -200
    xmax = 200
    agg_FUN = mean
    profile_value = ssvQC:::val2var[sts$signal_config$plot_value]
    profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value]

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

    # ?set
    diff_dt = agg_dt[, .(y = mean(get(profile_value))), .(mark, treatment, id, tx, ty)]
    diff_dt = dcast(diff_dt, id+tx+ty+mark~treatment, value.var = "y")
    diff_dt = diff_dt[mark != "input"]
    diff_dt$diff = diff_dt$`IK` - diff_dt$`GFP`

    lim = quantile(abs(diff_dt$diff), .99)
    diff_dt[diff > lim, diff := lim]
    diff_dt[diff < -lim, diff := -lim]
    p_combo_diff = ggplot(diff_dt, aes(x = tx, y = ty, color = diff)) +
        geom_point(size = base_point_size/10) +
        scale_color_gradientn(colours = c("blue", "gray90", "red"), limits = c(-lim, lim)) +
        coord_fixed() +
        facet_grid(mark~.) +
        labs(color = "IK - GFP") +
        theme(panel.background = element_blank(), panel.grid = element_blank())
    my_ggsave("4_combinatorial_tsne_points_FC", p_combo_diff)

    p_combo_diff.wide = ggplot(diff_dt, aes(x = tx, y = ty, color = diff)) +
        geom_point(size = base_point_size/10) +
        scale_color_gradientn(colours = c("blue", "gray90", "red"), limits = c(-lim, lim)) +
        coord_fixed() +
        facet_grid(.~mark) +
        labs(color = "IK - GFP") +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.text = element_text(size = 6))
    my_ggsave("4_combinatorial_tsne_points_FC.wide", p_combo_diff.wide, width = 9, height = 3.2)

    sts$signal_config$run_by = "treatment"
    sts$signal_config$color_by = "mark"

    mark_cols = IKdata:::mark_colors
    mark_cols["ATAC"] = "black"

    sts$signal_config$color_mapping = mark_cols[sel_marks]
    p = ctPlotSummaryProfiles(sts, extra_vars = "group")+
        coord_fixed()
    my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_mark", p)

    sts@signal_config@run_by = "mark"
    sts@signal_config@color_by = "treatment"
    sts@signal_config@color_mapping = IKdata:::treatment_colors


    sts@signal_config@color_mapping = treat_cols


    roi = list(
        xmin = .125, xmax = .225, ymin = .075, ymax = .175
    )
    roi = list(
        xmin = -.15, xmax = 0, ymin = .2, ymax = .4
    )
    if(my_cluster_value == "linearQuantile"){
        roi = list(xmin = .02, xmax = .13, ymin = -.04, ymax = .05)
    }else if(my_cluster_value == "RPM"){
        roi = list(xmin = 0, xmax = .2, ymin = -.08, ymax = .05)
    }else{
        stop()
    }


    p_summary = ctPlotSummaryProfiles(sts, extra_vars = "group", xbins = 15) +
        coord_fixed() +
        theme(panel.background = element_blank(), panel.grid = element_blank())
    p_summary
    my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment", p_summary, width = 9, height = 6)

    # .exp = .01
    #
    # roi = lapply(roi, function(x)x-.exp)

    p_summary.rect = p_summary + annotate("rect",
                                          xmin = roi$xmin, xmax = roi$xmax,
                                          ymin = roi$ymin, ymax = roi$ymax,
                                          fill = NA, color = "black")
    p_summary.rect
    my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment.rect", p_summary.rect, width = 9, height = 6)

    p_combo_diff.rect = p_combo_diff + annotate("rect",
                                                xmin = roi$xmin, xmax = roi$xmax,
                                                ymin = roi$ymin, ymax = roi$ymax,
                                                fill = NA, color = "black")
    p_combo_diff.rect
    my_ggsave("4_combinatorial_tsne_points_FC.rect", p_combo_diff.rect, width = 4, height = 9)

    p_combo_diff.wide.rect = p_combo_diff.wide + annotate("rect",
                                                          xmin = roi$xmin, xmax = roi$xmax,
                                                          ymin = roi$ymin, ymax = roi$ymax,
                                                          fill = NA, color = "black")
    my_ggsave("4_combinatorial_tsne_points_FC.wide.rect", p_combo_diff.wide.rect, width = 9, height = 3.2)

    exp = .02
    sts$signal_config$to_run = unique(as.character(sts$signal_config$meta_data$mark))
    sts$signal_config$plot_value = ssvQC:::sqc_signal_values$linearQuantile#"RPM"
    sts@signal_data$query_features$all_signal@signal_data =
        sts$signal_data$query_features$all_signal$signal_data[mark != "input"]
    p = ctPlotSummaryProfiles(
        sts,
        extra_vars = "group",
        xrng = c(roi$xmin-exp, roi$xmax+exp),
        yrng = c(roi$ymin-exp, roi$ymax+exp),
        N_floor = 1,
        N_ceiling = 10,
        xbins = 5,
        min_fraction = .2)+
        coord_fixed() #+
    # annotate("rect",
    #          xmin = roi$xmin, xmax = roi$xmax,
    #          ymin = roi$ymin, ymax = roi$ymax,
    #          fill = NA, color = "black")
    p = p +
        theme(panel.background = element_blank(), panel.grid = element_blank())
    # p + scale_y_continuous(breaks = seq(-.5, .5, .01))
    p = p + facet_grid(.~mark)
    p
    my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment.detail", p, width = 9, height = 3.2)

    p = ctPlotPoints(sts) +
        coord_fixed(xlim = c(roi$xmin-exp, roi$xmax+exp),
                    ylim = c(roi$ymin-exp, roi$ymax+exp)) +
        annotate("rect",
                 xmin = roi$xmin, xmax = roi$xmax,
                 ymin = roi$ymin, ymax = roi$ymax,
                 fill = NA, color = "black")
    my_ggsave("4_combinatorial_tsne_points.detail", p)


    sel_id = sts$signal_data$query_features$all_signal$xy_data[tx >= roi$xmin &
                                                                   tx <= roi$xmax &
                                                                   ty > roi$ymin &
                                                                   ty < roi$ymax]$id

    shrink = .02
    sel_id = sts$signal_data$query_features$all_signal$xy_data[tx >= roi$xmin + shrink &
                                                                   tx <= roi$xmax - shrink &
                                                                   ty > roi$ymin + shrink &
                                                                   ty < roi$ymax - shrink]$id

    ctPlotSummaryProfiles(sts, xbins = 5, extra_vars = "group",
                          xrng = c(roi$xmin, roi$xmax),
                          yrng = c(roi$ymin, roi$ymax)) +
        coord_fixed() +
        theme(panel.background = element_rect(fill = "gray60"))

    sel_dt = sts$signal_data$query_features$all_signal$signal_data[id %in% sel_id]
    sel_dt$cluster_id = NULL
    sel_dt.clust = ssvSignalClustering(sel_dt, facet_ = "group", fill_ = "y_linQ")
    sel_dt.clust[, c("mark", "treatment") := tstrsplit(group, " ")]
    sel_dt.agg = sel_dt.clust[, .(y_linQ = mean(y_linQ)), .(x, group, cluster_id, mark, treatment)]
    ggplot(sel_dt.agg, aes(x = x, y = y_linQ, color = treatment)) +
        geom_path() +
        scale_color_manual(values = IKdata:::treatment_colors) +
        facet_grid(cluster_id~mark)

    ssvSignalHeatmap(sel_dt.clust, facet_ = "group", fill_ = "y_linQ") +
        scale_fill_viridis_c()

    sel_dt.assign = unique(sel_dt.clust[, .(id, cluster_id)])

    tmp = as.data.table(query_gr[sel_id])
    fwrite(tmp[, .(seqnames, start, end)], res_file("selected_ids.bed"), sep = "\t", col.names = FALSE)
    # agg_dt = sts$signal_data$query_features$all_signal$signal_data[, .(y_linQ = mean(y_linQ)), .(x, group, mark, treatment)]
    # ggplot(agg_dt, aes(x = x, y = y_linQ, color = treatment)) +
    #     geom_path() +
    #     scale_color_manual(values = IKdata:::treatment_colors) +
    #     facet_grid(.~mark)
    # dev.off()

    nice_result_gr = rtracklayer::import.bed("output_presentation_figures.102622/selected_ids.bed")
    nice_hits = subsetByOverlaps(sts$signal_data$query_features$all_signal$query_gr,
                                 nice_result_gr) %>% names

    p_gr = sts$signal_data$query_features$all_signal$query_gr
    olaps = findOverlaps(query = resize(p_gr, 50, fix = "center"), subject = nice_result_gr)

    p_df = as.data.table(p_gr)
    p_df$id = names(p_gr)
    p_df[, nice_hit := id %in% nice_hits]



    p_nice_hit = ctPlotPointsAnnotation(sts, meta_data = p_df, anno_var = "nice_hit", point_size = base_point_size / 5) +
        annotate("rect",
                 xmin = roi$xmin, xmax = roi$xmax,
                 ymin = roi$ymin, ymax = roi$ymax,
                 fill = NA, color = "black")

    my_ggsave("4_combinatorial_nice_hit", p_nice_hit)

    #### 5 tracks ####


    qc_config_signal$meta_data

    ex_gr =rtracklayer::import.gff("~/gencode.v36.annotation.gtf", feature.type = "exon")

    goi = "FFAR2"
    goi = "ERG"
    goi = "MYC"
    goi = "TGFB1"
    goi = "CD5"
    goi = "S1PR4"
    goi = "ADGRG3"
    goi = "KMT2E"
    goi_query_gr = range(subset(ex_gr, gene_name == goi))
    goi_query_gr = resize(goi_query_gr, width = width(goi_query_gr)+1e5, fix = "center")
    goi_query_gr = shift(goi_query_gr, -3e4)

    goi = sel_id[2]
    goi_query_gr = sts$signal_data$query_features$all_signal$query_gr[goi]

    goi_query_gr = resize(goi_query_gr, width = 2.3e4, fix = "center")
    # goi_query_gr = GRanges("chr19", IRanges(35435e3, 35470e3))
    # goi_query_gr = GRanges("chr21", IRanges(38350e3, 38450e3), strand = "-")
    # goi_query_gr = GRanges("chr11", IRanges(61080e3, 61150e3), strand = "+")
    query_meta_df = qc_config_signal$meta_data
    # undebug(track_chip)

    # save(qc_config_signal, goi_query_gr, treat_cols, file = "tmp.save")
    # load("tmp.save")
    t_chip = track_chip(subset(qc_config_signal$meta_data, mark != "input"), nwin = 80, nspline = 5,
                        query_gr = goi_query_gr,
                        color_VAR = "treatment",
                        color_VAR_VAR = "treatment", fill_VAR = NULL,
                        facet_VAR = "mark", color_mapping = treat_cols) #+
    #theme(panel.background = element_rect(fill = "gray60"))

    t_chip = t_chip + facet_grid(mark~., scales = "free_y")
    t_chip

    all_gr = sts$signal_data$query_features$all_signal$query_gr
    strand(all_gr) = "*"
    sel_gr = all_gr[sel_id]
    distanceToNearest(sel_gr, goi_query_gr)
    # sel_gr[257]
    # undebug(track_features)
    .check_query_gr = ssvTracks:::.check_query_gr
    .apply_x_scale = ssvTracks:::.apply_x_scale
    .apply_x_lim = ssvTracks:::.apply_x_lim
    t_features = track_features(list(All = all_gr, Selected = sel_gr), query_gr = goi_query_gr)


    relevant_goi = union(goi, c("TAB1", "SYNGR1", "CACNA1I"))
    # track_gene_reference(ex_gr, goi_query_gr)
    t_ref = track_gene_reference(subset(ex_gr, gene_name %in% relevant_goi), goi_query_gr, show_tss = TRUE)
    t_ref = track_gene_reference(ex_gr, goi_query_gr, show_tss = TRUE)


    sel_gir.hit = subsetByOverlaps(sel_gr, goi_query_gr)

    pg = assemble_tracks(
        list(t_chip+labs(title = goi), t_ref, t_features),
        query_gr = goi_query_gr, rel_heights = c(3, 1, 1)
    )


    pg

    ggsave(res_file(paste0("track_", goi, ".pdf")), pg, width = 9.1, height = 8.5)

