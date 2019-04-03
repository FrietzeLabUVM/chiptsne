## This is how test files were generated
# bw_files = dir("/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court_hg19/",
#                pattern = "pooled", full.names = TRUE) %>% dir(pattern = "FE.bw$", full.names = TRUE)
# cfg_dt = data.table(bw_files)
# cfg_dt[, c("cell", "mark") := tstrsplit(basename(bw_files), "_", keep = 1:2)]
# cfg_dt$norm_factor = 1
# cfg_dt[mark == "H3K4me3", norm_factor := .3]
#
# biv19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_bivalent_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
# biv19$group = "bivalent"
# k4me3_19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_k4m3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
# k4me3_19$group = "k4me3_alone"
# k27me3_19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_k27me3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
# k27me3_19$group = "k27me3_alone"
# query_gr = c(biv19, k4me3_19, k27me3_19)
# query_gr = resize(query_gr, 5000, fix = "center")
#
# query_gr =sample(query_gr, 30)
#
# query_gr$id = paste0("region_", seq_along(query_gr))
#
# cfg_dt = cfg_dt[cell %in% c("ESH1", "HUES48", "HUES64")]
#
# chrSizes = fread("~/hg38_chrsizes.canon.txt")
# cs = chrSizes$V2
# names(cs) = chrSizes$V1
#generate test data bigwigs
# for(bw_f in cfg_dt$bw_files){
#     outf = file.path("inst/extdata", basename(bw_f))
#     bw_gr = rtracklayer::import.bw(bw_f, which = query_gr)
#     bw_gr = GRanges(seqsetvis::viewGRangesWinSummary_dt(bw_gr, query_gr, n_tiles = 50)[, .(seqnames, start, end, score = y)])
#
#     seqlengths(bw_gr) = cs[names(seqlengths(bw_gr))]
#     rtracklayer::export.bw(bw_gr, outf)
# }
