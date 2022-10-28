#### setup ####
library(chiptsne)
library(seqsetvis)
library(ggplot2)
library(magrittr)

base_point_size = .3

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

pdf("presentation_figures.pdf")

options(mc.cores = 20)
chip_bam_cfg_dt = IKdata::chipseq.setup_bam_files()
chip_np_cfg_dt = IKdata::chipseq.setup_peak_files()

col_map = IKdata:::mark_colors
col_map["H3K27ac"] = safeBrew(8, "set1")[2]
col_map["H3K4me1"] = safeBrew(8, "set1")[4]
col_map["input"] = "gray80"
theme_update(panel.background = element_blank(), panel.grid = element_blank())

#### single k27ac rep ####
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
                                  plot_value = SQC_SIGNAL_VALUES$RPM,
                                  cluster_value = SQC_SIGNAL_VALUES$RPM,
                                  sort_value = SQC_SIGNAL_VALUES$RPM,
                                  run_by = "All",
                                  color_by = "mark",
                                  color_mapping = col_map.sel)

sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
set.seed(0)
sts = ChIPtSNE.runTSNE(sts)
sts$signal_config$plot_value = SQC_SIGNAL_VALUES$RPM
sts$signal_config$cluster_value = SQC_SIGNAL_VALUES$RPM
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

p = ctPlotPoints(sts, xmin = -100, xmax = 100, point_size = base_point_size/2) +
    facet_wrap(~mark) +
    coord_fixed() +
    scale_color_viridis_c(limits = c(0, 5), na.value = "yellow") +
    labs(color = "max\nRPM\npileup")
p
my_ggsave("1_peak_eval_tsne_points", p)

p = ctPlotSummaryProfiles(sts,
                          N_floor = 0,
                          N_ceiling = 150,
                          min_fraction = .25) +
    facet_null() +
    coord_fixed()
p
my_ggsave("1_peak_eval_tsne_summary_profiles", p)


#### k27ac rep comparison ####
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
                                  plot_value = SQC_SIGNAL_VALUES$RPM,
                                  cluster_value = SQC_SIGNAL_VALUES$RPM,
                                  sort_value = SQC_SIGNAL_VALUES$RPM,
                                  run_by = "All",
                                  color_by = "rep",
                                  color_mapping = col_map.sel)

sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
set.seed(0)
sts = ChIPtSNE.runTSNE(sts)

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

#### k27ac differential ####
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
                                  plot_value = SQC_SIGNAL_VALUES$RPM,
                                  cluster_value = SQC_SIGNAL_VALUES$RPM,
                                  sort_value = SQC_SIGNAL_VALUES$RPM,
                                  run_by = "All",
                                  color_by = "treatment",
                                  color_mapping = col_map.sel)

sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
set.seed(0)
sts = ChIPtSNE.runTSNE(sts)

t_dt = sts$signal_data$query_features$all_signal$xy_data
p = ggplot(t_dt, aes(x = tx, y = ty)) +
    geom_point(size = base_point_size) +
    coord_fixed()
my_ggsave("3_diff_eval_tsne_pointsBW", p)

p = ctPlotPoints(sts, xmin = -100, xmax = 100, point_size = base_point_size) +
    # facet_wrap(~mark) +
    coord_fixed()
my_ggsave("3_diff_eval_tsne_points", p)

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
bed_files.chipseq[1]

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


#### dev diff plot ####
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

#### k27ac with k4me1, Ikaros, ATAC ####
sel_marks = c("H3K27ac", "H3K4me1", "IKAROS", "input", "ATAC")
atac_bam_cfg_dt = IKdata::atac.setup_bam_files()[cell == "MXP5"]
atac_bam_cfg_dt$mark = "ATAC"
bam_cfg_dt.sel = chip_bam_cfg_dt[mark %in% sel_marks]

bam_cfg_dt.sel = rbind(bam_cfg_dt.sel, atac_bam_cfg_dt)

bam_cfg_dt.sel$All = "all"
np_cfg_dt.sel = chip_np_cfg_dt[mark %in% sel_marks]
col_map.sel = IKdata:::treatment_colors
plot(1:2, col = col_map.sel, pch = 16, cex = 4)

atac_dar = IKdata::atac.DAR_results_load()
atac_dar = subset(atac_dar, DA_group != "ns")

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
                                  plot_value = SQC_SIGNAL_VALUES$RPM,
                                  cluster_value = SQC_SIGNAL_VALUES$RPM,
                                  sort_value = SQC_SIGNAL_VALUES$RPM,
                                  run_by = "All",
                                  color_by = "group",
                                  color_mapping = col_map.sel)

sts = ChIPtSNE(features_config = qc_config_feat, signal_config = qc_config_signal)
set.seed(0)
sts = ChIPtSNE.runTSNE(sts)
sts$signal_config$plot_value = SQC_SIGNAL_VALUES$linearQuantile
sts$signal_config$cluster_value = SQC_SIGNAL_VALUES$linearQuantile
sts$signal_config$sort_value = SQC_SIGNAL_VALUES$linearQuantile

p = ctPlotPoints(sts, point_size = base_point_size / 2) +
    facet_grid(mark~treatment) +
    coord_fixed() +
    scale_color_viridis_c(limits = c(0, 3), na.value = "yellow")
my_ggsave("4_combinatorial_tsne_points", p)

#### diff ####
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
p = ggplot(diff_dt, aes(x = tx, y = ty, color = diff)) +
    geom_point(size = base_point_size) +
    scale_color_gradientn(colours = c("blue", "gray90", "red"), limits = c(-lim, lim)) +
    coord_fixed() +
    facet_grid(mark~.) +
    labs(color = "IK - GFP")
my_ggsave("4_combinatorial_tsne_points_FC", p)


sts$signal_config$run_by = "treatment"
sts$signal_config$color_by = "mark"

mark_cols = IKdata:::mark_colors
mark_cols["ATAC"] = "black"

sts$signal_config$color_mapping = mark_cols[sel_marks]
p = ctPlotSummaryProfiles(sts, extra_vars = "group")+
    coord_fixed()
my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_mark", p)

sts$signal_config@run_by = "mark"
sts@signal_config@color_by = "treatment"

p_summary = ctPlotSummaryProfiles(sts, extra_vars = "group") +
    coord_fixed()

my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment", p_summary)

p_summary.rect = p_summary + annotate("rect", xmin = -.15, xmax = .2, ymin = -.1, ymax = .25, fill = NA, color = "black")
my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment.rect", p_summary.rect)

p = ctPlotSummaryProfiles(sts, extra_vars = "group", xrng = c(-.15, .2), yrng = c(-.1, .25))+
    coord_fixed()
my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment.detail", p)

sts$signal_config$run_by = "treatment"
sts$signal_config$color_by = "mark"
sts@signal_data$query_features$all_signal@signal_data = sts$signal_data$query_features$all_signal$signal_data[mark != "input"]

p = ctPlotSummaryProfiles(sts, extra_vars = "group", xrng = c(-.15, .2), yrng = c(-.1, .25))+
    coord_fixed()

my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_mark.detail", p)


sts@signal_config@run_by = "mark"
sts@signal_config@color_by = "treatment"
sts@signal_config@color_mapping = IKdata:::treatment_colors

p = ctPlotSummaryProfiles(sts, extra_vars = "group", xrng = c(-.15, .2), yrng = c(-.1, .25))+
    coord_fixed()

my_ggsave("4_combinatorial_tsne_summary_profiles_colorBy_treatment.detail2", p)

sel_id = sts$signal_data$query_features$all_signal$xy_data[tx > -.05 & tx < .08 & ty > 0 & ty < .15]$id

ctPlotSummaryProfiles(sts, xbins = 5, extra_vars = "group", xrng = c(-.05, .08), yrng = c(0, .15)) +
    coord_fixed()

# ctPlotPoints(sts, point_size = 1) +
#     facet_grid(mark~treatment) +
#     coord_fixed(xlim = c(-.05, .08), ylim = c(0, .15)) +
#     scale_color_viridis_c(limits = c(0, 1), na.value = "yellow")


sel_dt = sts$signal_data$query_features$all_signal$signal_data[id %in% sel_id]
sel_dt$cluster_id = NULL
sel_dt.clust = ssvSignalClustering(sel_dt, facet_ = "group", fill_ = "y_linQ")
sel_dt.clust[, c("mark", "treatment") := tstrsplit(group, " ")]
sel_dt.agg = sel_dt.clust[, .(y_linQ = mean(y_linQ)), .(x, group, cluster_id, mark, treatment)]
ggplot(sel_dt.agg, aes(x = x, y = y_linQ, color = treatment)) +
    geom_path() +
    scale_color_manual(values = IKdata:::treatment_colors) +
    facet_grid(cluster_id~mark)

ssvSignalHeatmap(sel_dt.clust, facet_ = "group", fill_ = "y_linQ")

sel_dt.assign = unique(sel_dt.clust[, .(id, cluster_id)])

tmp = as.data.table(query_gr[sel_id])
fwrite(tmp[, .(seqnames, start, end)], res_file("selected_ids.bed"), sep = "\t", col.names = FALSE)
# agg_dt = sts$signal_data$query_features$all_signal$signal_data[, .(y_linQ = mean(y_linQ)), .(x, group, mark, treatment)]
# ggplot(agg_dt, aes(x = x, y = y_linQ, color = treatment)) +
#     geom_path() +
#     scale_color_manual(values = IKdata:::treatment_colors) +
#     facet_grid(.~mark)
dev.off()
