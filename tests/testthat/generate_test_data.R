# # This is how test files were generated
# if(FALSE){
    library(data.table)
    library(magrittr)
    library(GenomicRanges)
    library(Rsamtools)
    library(S4Vectors)
    set.seed(0)
    bw_files = dir("/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court_hg19/",
                   pattern = "pooled", full.names = TRUE) %>% dir(pattern = "FE.bw$", full.names = TRUE)
    cfg_dt = data.table(bw_files)
    cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(bw_files), "_", keep = 1:2)]
    cfg_dt$norm_factor = 1
    cfg_dt[wide_var == "H3K4me3", norm_factor := .3]
    cfg_dt[, norm_factor := norm_factor / 20]
    cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]

    bam_files = dir("/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court_hg19/",
                    pattern = "pooled", full.names = TRUE) %>% dir(pattern = "bam$", full.names = TRUE)
    cfg_dt.bam = data.table(bam_files)
    cfg_dt.bam[, c("tall_var", "wide_var") := tstrsplit(basename(bam_files), "_", keep = 1:2)]
    cfg_dt.bam$norm_factor = 1
    cfg_dt.bam[wide_var == "H3K4me3", norm_factor := .3]
    cfg_dt.bam[, norm_factor := norm_factor / 20]
    cfg_dt.bam = cfg_dt.bam[tall_var %in% c("ESH1", "HUES48", "HUES64")]

    # ex_cfg_dt = cfg_dt
    # ex_cfg_dt.bam = cfg_dt.bam
    # usethis::use_data(ex_cfg_dt, overwrite = TRUE)
    # usethis::use_data(ex_cfg_dt.bam, overwrite = TRUE)
    #
    biv19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_bivalent_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
    biv19$group = "bivalent"
    k4me3_19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_k4m3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
    k4me3_19$group = "k4me3_alone"
    k27me3_19 = fread("~/R/waldron/shiny_chiptsne/PMC5354816_k27me3_hg19.tsv", col.names = c("seqnames", "start", "end")) %>% GRanges
    k27me3_19$group = "k27me3_alone"

    data("query_gr")
    # query_gr = c(biv19, k4me3_19, k27me3_19)
    # query_gr = resize(query_gr, 5000, fix = "center")
    #
    # query_gr =sample(query_gr, 30)
    # query_gr$id = paste0("region_", seq_along(query_gr))
    #
    #
    chrSizes = fread("~/hg38_chrsizes.canon.txt")
    cs = chrSizes$V2
    names(cs) = chrSizes$V1
    # generate test data bigwigs
    for(bw_f in cfg_dt$bw_files){
        outf = file.path("inst/extdata", basename(bw_f))
        bw_gr = rtracklayer::import.bw(bw_f, which = query_gr)
        sl = seqlengths(bw_gr)
        bw_gr = GRanges(seqsetvis::viewGRangesWinSummary_dt(bw_gr, query_gr, n_tiles = 50)[, .(seqnames, start, end, score = y)])
        bw_gr = GRanges(as.character(seqnames(bw_gr)), IRanges(start(bw_gr), end(bw_gr)), seqlengths = sl, score = bw_gr$score)
        seqlengths(bw_gr)
        # seqlengths(bw_gr) = cs[names(seqlengths(bw_gr))]
        rtracklayer::export.bw(bw_gr, outf)
    }


    sbp = ScanBamParam(which = query_gr, what = "qname")


    fraction = .05

    # filts <- c("peaks", "promoters")
    # filts <- list(subsample = function(rd){sample(nrow(rd), prob = fraction)})#runif(nrow(rd)) < fraction})
    filts <- list(subsample = function(rd){runif(nrow(rd)) < fraction})
    filters <- FilterRules(filts)
    # active(filters) # all TRUE

    bam_out = parallel::mclapply(cfg_dt.bam$bam_files, function(bam_f){
        outf = file.path("inst/extdata", basename(bam_f))
        Rsamtools::filterBam(bam_f, param = sbp, destination = outf, filter = filters)

        sortf = sub(".bam$", ".sorted", outf)
        sortf = Rsamtools::sortBam(outf, sortf)
        file.remove(outf)
        file.remove(paste0(outf, ".bai"))
        Rsamtools::indexBam(sortf)
        sortf
    }) %>% unlist

    lapply(bam_out, function(b){
        sum(idxstatsBam(b)$mapped)
    })

    #
    #
    # dt_raw = seqsetvis::ssvFetchBam(cfg_dt.bam, query_gr, target_strand = "both", return_data.table = TRUE)
    # cfg_dt.filtered = cfg_dt.bam
    # cfg_dt.filtered[, bam_files := file.path("inst/extdata", basename(bam_files))]
    # dt_filtered = seqsetvis::ssvFetchBam(cfg_dt.filtered, query_gr, target_strand = "both", return_data.table = TRUE)
    #
    # dt_raw$group = "raw"
    # dt_filtered$group = "filtered"
    # dt = rbind(dt_raw, dt_filtered)
    #
    # dt[id == "region_1"]
    # library(ggplot2)
    # ggplot(dt[tall_var == "ESH1" & id %in% unique(id)[1:3]], aes(x = x, y = y, color = strand)) +
    #     geom_path() + facet_grid(id~group+wide_var)
    #
    # dcast(dt[tall_var == "ESH1" & id %in% unique(id)[1] & wide_var == "H3K4me3"], group~x+strand, value.var = "y")
# }
