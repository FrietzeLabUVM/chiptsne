#functions for fetchng and formatting genomic data for t-sne then running t-sne.

fetch_bam_dt = function(qdt,
                        qgr,
                        qwin = 50,
                        qmet = "summary",
                        agg_FUN = sum,
                        cap_value = 20,
                        high_on_right = TRUE,
                        bfc = BiocFileCache::BiocFileCache(),
                        n_cores = getOption("mc.cores", 1),
                        rname = digest::digest(list(qgr, qdt[, 1:3], qwin, qmet)),
                        force_overwrite = FALSE){
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(is.null(qdt$norm_factor)){
        qdt$norm_factor = 1
    }
    if(!all(file.exists(qdt[[1]]))){
        stop(paste(sep = "\n", "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt[[1]][!file.exists(qdt[[1]])])))

    }

    bam_fetch = function(){
        seqsetvis::ssvFetchBam(qdt[, 1:3], qgr,
                               return_data.table = TRUE,
                               win_method = qmet,
                               win_size = qwin, n_cores = n_cores)
    }
    bam_dt = bfcif(bfc, rname, bam_fetch, force_overwrite = force_overwrite)
    bam_dt$sample = NULL
    bam_dt = bam_dt[, list(y = agg_FUN(y)), list(cell, id, mark, x)]
    bam_dt
}

fetch_bw_dt = function(qdt,
                       qgr,
                       qwin = 50,
                       qmet = "summary",
                       agg_FUN = mean,
                       cap_value = 20,
                       high_on_right = TRUE,
                       bfc = BiocFileCache::BiocFileCache(),
                       n_cores = getOption("mc.cores", 1),
                       rname = digest::digest(list(qgr, qdt[, 1:3], qwin, qmet)),
                       force_overwrite = FALSE){
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(!all(file.exists(qdt[[1]]))){
        stop(paste(sep = "\n", "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt[[1]][!file.exists(qdt[[1]])])))

    }

    bw_fetch = function(){
        seqsetvis::ssvFetchBigwig(qdt[, 1:3], qgr,
                                  return_data.table = TRUE,
                                  win_method = qmet,
                                  win_size = qwin, n_cores = n_cores)
    }

    bw_dt = bfcif(bfc, rname, bw_fetch, force_overwrite = force_overwrite)
    # bw_dt$sample = NULL
    bw_dt = bw_dt[, list(y = agg_FUN(y)), list(cell, id, mark, x)]
    bw_dt
}

prep_profile_dt = function(prof_dt, norm_dt, qgr,
                           cap_value = 20,
                           high_on_right = TRUE){
    if(!all(norm_dt$norm_factor == 1)){
        prof_dt = merge(prof_dt, norm_dt)
        prof_dt[, y := y * norm_factor]
        prof_dt$norm_factor = NULL
    }
    prof_dt[y > cap_value, y := cap_value]

    if(high_on_right){
        balance_dt = prof_dt[, list(right_sum = sum(y[x > 0]), left_sum = sum(y[x < 0])), by = list(cell, id)]
        balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum, cell, id)]
        most_flipped = balance_dt[, list(fraction_flipped = sum(needs_flip) / .N), by = list(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        GenomicRanges::strand(qgr) = "+"
        GenomicRanges::strand(qgr)[most_flipped$flip_strand] = "-"
        prof_dt = merge(prof_dt, balance_dt, by = c("id", "cell"))
        remove(balance_dt)
        prof_dt[needs_flip == TRUE, x := -x]
        prof_dt$needs_flip = NULL
    }
    prof_dt
}

dt2mat = function(prof_dt, marks){
    for(m in marks){
        if(!exists("...tsne_mat")){
            dt = dcast( prof_dt[mark == m], id+cell~x, value.var = "y")
            ...tsne_mat = as.matrix(dt[, -1:-2])
            rn = paste(dt$id, dt$cell)
        }else{
            dt = dcast( prof_dt[mark == m], id+cell~x, value.var = "y")
            stopifnot(all(paste(dt$id, dt$cell) == rn))
            ...tsne_mat = cbind(...tsne_mat, as.matrix(dt[, -1:-2]))


        }
    }
    rownames(...tsne_mat) = rn
    ...tsne_mat
}

#' Title
#'
#' @param qdt
#' @param qgr
#' @param qwin
#' @param qmet
#' @param cap_value
#' @param high_on_right
#' @param bfc
#' @param n_cores
#' @param rname
#' @param force_overwrite
#'
#' @return
#' @export
#' @importFrom seqsetvis ssvFetchBigwig
#'
#' @examples
stsFetchTsneInput = function(qdt,
                             qgr,
                             file_format = NULL,
                             qwin = 50,
                             qmet = "summary",
                             cap_value = 20,
                             agg_FUN = NULL,
                             high_on_right = TRUE,
                             bfc = BiocFileCache::BiocFileCache(),
                             n_cores = getOption("mc.cores", 1),
                             rname = digest::digest(list(qgr, qdt[, 1:3], qwin, qmet)),
                             force_overwrite = FALSE){
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(is.null(qdt$norm_factor)){
        qdt$norm_factor = 1
    }
    if(!all(file.exists(qdt[[1]]))){
        stop(paste(sep = "\n", "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt[[1]][!file.exists(qdt[[1]])])))

    }
    if(is.null(file_format)){#try to determine file format
        test_file = tolower(qdt[[1]][1])
        if(grepl(".bw$", test_file) | grepl(".bigwig$", test_file)){
            file_format = "bw"
        }else if(grepl(".bam$", test_file)){
            file_format = "bam"
        }else{
            stop("couldn't determine file_format of ", qdt[[1]][1],
                 "\nplease manually set file_format to bam or bw (bigWig).")
        }

    }
    stopifnot(file_format %in% c("bam", "bw"))
    if(is.null(agg_FUN)){
        agg_FUN = switch(file_format,
                         bam = sum,
                         bw = mean)

    }
    fetch_FUN = switch(file_format,
                       bam = fetch_bam_dt,
                       bw = fetch_bw_dt)
    bw_dt = fetch_FUN(qdt, qgr, qwin, qmet, agg_FUN, cap_value, high_on_right, bfc)
    # apply transforms and thresholds
    bw_dt = prep_profile_dt(bw_dt,
                            unique(qdt[, list(cell, mark, norm_factor)]),
                            qgr,
                            cap_value = cap_value,
                            high_on_right = high_on_right)
    # cast from tidy to wide matrix
    bw_mat = dt2mat(bw_dt, unique(qdt$mark))

    return(list(bw_dt = bw_dt, bw_mat = bw_mat, query_gr = qgr))
}

#' Title
#'
#' @param tsne_mat
#' @param perplexity
#' @param n_cores
#' @param high_topright
#' @param norm1
#' @param bfc
#' @param rname
#' @param force_overwrite
#'
#' @return
#' @export
#' @import Rtsne
#' @examples
stsRunTsne = function(tsne_mat, perplexity = 100,
                      n_cores = getOption("mc.cores", 1),
                      high_topright = TRUE,
                      norm1 = TRUE,
                      bfc = BiocFileCache::BiocFileCache(),
                      rname = NULL,
                      force_overwrite = FALSE){
    set.seed(0)
    if(is.null(rname)){
        rname = digest::digest(list(
            tsne_mat[sample(1:nrow(tsne_mat), 20),
                     sample(1:ncol(tsne_mat), 20)],
            perplexity
        ))
    }
    res_tsne = bfcif(bfc, rname, force_overwrite = force_overwrite,
                     FUN = function(){
                         Rtsne(tsne_mat, num_threads = n_cores, perplexity = perplexity, check_duplicates = FALSE)
                     })

    tdt = as.data.table(res_tsne$Y)
    colnames(tdt) = c("tx", "ty")
    tdt$rn = rownames(tsne_mat)
    tdt[, c("id", "cell") := tstrsplit(rn, " ", keep = 1:2)]

    if(norm1){
        tdt$tx = rescale_capped(tdt$tx)-.5
        tdt$ty = rescale_capped(tdt$ty)-.5
    }


    if(high_topright){
        rs = rowSums(tsne_mat)
        tdt$rs = rs[tdt$rn]
        x_cutoff = mean(range(tdt$tx))
        x_flip = sum(tdt[tx > x_cutoff]$rs) < sum(tdt[tx < x_cutoff]$rs)
        if(x_flip){
            tdt[, tx := max(tx) - tx + min(tx)]
        }
        y_cutoff = mean(range(tdt$ty))
        y_flip = sum(tdt[ty > y_cutoff]$rs) < sum(tdt[ty < y_cutoff]$rs)
        if(y_flip){
            tdt[, ty := max(ty) - ty + min(ty)]
        }
        tdt$rs = NULL
    }

    tdt
}








