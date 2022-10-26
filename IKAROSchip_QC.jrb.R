#IKAROS reps- QC

library(ssvQC)
library(data.table)
library(tidyverse)
library(seqsetvis)


peak_f <-c(   "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/CTRL_IK_R2_macs2_peaks.narrowPeak",
  "/slipstream/home/arrichman/prebALL/chromatin/IKinduction/ChIPseq/AR314_EnhIkaros/MXP5_GFP_IKAROS_rep2.Aligned.sortedByCoord.out_macs2_peaks.narrowPeak",
            "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/IK1_IK_R2_macs2_peaks.narrowPeak",
            "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/IK1_IK_R1_macs2_peaks.narrowPeak"
          )
           # "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/CTRL_IK_R1_macs2_peaks.narrowPeak")



bam_f <- c("/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/CTRL_IK_R2.Aligned.sortedByCoord.out.bam",
  "/slipstream/home/arrichman/prebALL/chromatin/IKinduction/ChIPseq/AR314_EnhIkaros/MXP5_GFP_IKAROS_rep2.Aligned.sortedByCoord.out.bam",
           "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/IK1_IK_R1.Aligned.sortedByCoord.out.bam",
           "/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/IK1_IK_R2.Aligned.sortedByCoord.out.bam"
           )
           #"/slipstream/home/arrichman/prebALL/IKbinding/HP_060222_ChIP/macs2output/CTRL_IK_R1.Aligned.sortedByCoord.out.bam")



stopifnot(file.exists(peak_f))
stopifnot(file.exists(bam_f))
stopifnot(file.exists(paste0(bam_f, ".bai")))

#features config
meta_dt.feat = data.table(file = peak_f)
meta_dt.feat[, "full_name" := tstrsplit(basename(file), "/", keep = 1)]
#meta_dt.feat[, c("sample", "treatment", "mark","rep") := tstrsplit(basename(file), "_", keep = seq(4))]
#meta_dt.feat[, "rep" := tstrsplit(basename(name), ".", keep=1)]
meta_dt.feat$sample <- "MXP5"
meta_dt.feat$treatment <- c("GFP", "GFP","IK", "IK")
meta_dt.feat$mark <- "IKAROS"
meta_dt.feat$rep <- c("rep1","rep2", "rep1", "rep2")

#meta_dt.feat <- meta_dt.feat[,-5]
#name and name_split are required
meta_dt.feat[, name := paste(sample,treatment, mark,rep, sep = "_")]
meta_dt.feat[, name_split := paste(sample,treatment,mark,rep, sep = "\n")]
#order matters for keep IgG at the same spot in each plot
#meta_dt.feat = meta_dt.feat[order(factor == "IgG")]

meta_dt.feat

###seqsetvis
peak_grs = easyLoad_narrowPeak(meta_dt.feat$file, file_names = meta_dt.feat$name)
lengths(peak_grs)

#make condition groups
peak_groups = sub("_rep.+", "", names(peak_grs))
peak_grs.grouped = split(peak_grs, peak_groups)

#generate high consensus peak sets for the replicates
con_grs = ssvConsensusIntervalSets(peak_grs.grouped[["MXP5_GFP_IKAROS"]])
ik_grs = ssvConsensusIntervalSets(peak_grs.grouped[["MXP5_IK_IKAROS"]])

#generate overlap interval set for High Consensus Peaks
All_HC_peaks = list(GFP = con_grs, IK = ik_grs)
olaps = ssvOverlapIntervalSets(All_HC_peaks)

ssvFeatureVenn(olaps)
ssvFeatureEuler(olaps)

###Configure reads/signal
meta_dt.signal = data.table(file = bam_f)
meta_dt.signal[, "full_name" := tstrsplit(basename(file), "/", keep = 1)]
#meta_dt.signal[, c("sample", "treatment", "mark", "rep") := tstrsplit(basename(file), "_", keep = seq(4))]
meta_dt.signal$sample <- "MXP5"
meta_dt.signal$treatment <- c("GFP", "GFP" , "IK", "IK")
meta_dt.signal$mark <- "IKAROS"
meta_dt.signal$rep <- c("rep1","rep2", "rep1", "rep2")

meta_dt.signal
meta_dt.signal[, name := paste(sample,treatment,mark,rep, sep = "_")]
meta_dt.signal[, name_split := paste(sample,treatment,mark, rep, sep = "\n")]
meta_dt.signal[, mapped_reads := get_mapped_reads(file)]

query_gr = olaps
query_gr = prepare_fetch_GRanges_names(query_gr)
query_gr = resize(query_gr, 3e3, fix = "center")

options(mc.cores = 20)
prof_dt = ssvFetchBam(meta_dt.signal, unique_names = meta_dt.signal$name, query_gr, return_data.table = TRUE, n_region_splits = 50, fragLens = 180)
#add peak call group info
memb_df = ssvFactorizeMembTable(query_gr)
prof_dt = merge(prof_dt, memb_df, by = "id")
#sort within peak group


##### raw reads #### 
prof_dt = within_clust_sort(prof_dt, cluster_ = "group")
ssvSignalHeatmap.ClusterBars(prof_dt, 
                             fill_limits = c(0, 500), 
                             cluster_ = "group", 
                             FUN_format_heatmap = function(p)p + labs(fill = "raw pileup"))
##### RPM ####
prof_dt[, RPM := y / mapped_reads * 1e6]
ssvSignalHeatmap.ClusterBars(prof_dt, 
                             fill_ = "RPM",
                             fill_limits = c(0, 10), 
                             cluster_ = "group", 
                             FUN_format_heatmap = function(p)p + labs(fill = "RPM pileup"))


##### linear q95 capped #### 
cap_dt = calc_norm_factors(prof_dt)
prof_dt = append_ynorm(prof_dt, cap_dt = cap_dt)
ssvSignalHeatmap.ClusterBars(prof_dt, 
                             fill_ = "y_norm",
                             fill_limits = c(0, 1), 
                             cluster_ = "group", 
                             FUN_format_heatmap = function(p)p + labs(fill = "normalized pileup"))

agg_dt = prof_dt[, .(y = mean(y), RPM = mean(RPM), y_norm = mean(y_norm)), .(sample, mark, rep, treatment, group, x)]
ggplot(agg_dt, aes(x = x, y = y_norm, color = treatment, linetype = rep, group = sample)) +
  geom_path() + 
  facet_grid(group~mark) +
  scale_color_manual(values = c(GFP = "gray", IK = "purple")) +
  theme_classic()
