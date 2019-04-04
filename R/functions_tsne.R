#functions for fetchng and formatting genomic data for t-sne then running t-sne.

fetch_bam_dt = function(qdt,
                        qgr,
                        qwin = 50,
                        qmet = "summary",
                        agg_FUN = sum,
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


#' prepares query_gr GRanges for using in fetch functions
#'
#' fixed widths are recommended but not required
#'
#' metadata column id must be set and all ids must be unique
#'
#' @param query_gr GRanges of regions to be fetched
#' @param region_size numeric fixed width to apply to query_gr. if NULL, width
#'   is not enforced.
#' @param id_prefix character prefix to use for ids. 'region' is default and
#'   will yield ids such as region_1, region_2, region_3 ...
#'
#' @return GRanges ready for fetching
#' @export
#' @examples
#' library(GRanges)
#' qgr = GRanges("chr1", IRanges(1:10+100, 1:10+150))
#' prep_query_gr(qgr)
#' #region_size and id_prefix can be customized
#' prep_query_gr(qgr, region_size = 70, id_prefix = "peak")
#' names(qgr) = LETTERS[seq_along(qgr)]
#' prep_query_gr(qgr)
#' qgr$name = letters[seq_along(qgr)]
#' prep_query_gr(qgr)
#' #invalid - i.e. duplicated ids will be overwritten
#' qgr$id = rep("a", length(qgr))
#' prep_query_gr(qgr)
prep_query_gr = function(query_gr,
                         region_size = NULL,
                         id_prefix = "region"){
    if(!is.null(region_size)){
        query_gr = resize(query_gr, region_size, fix = "center")
    }
    #successively check $id, $name, and names()
    if(!is.null(query_gr$id)){
        if(!any(duplicated(query_gr$id))){
            message("using existing id")
            names(query_gr) = NULL
            query_gr$name = NULL
            return(query_gr)
        }
        message("ignoring id attribute with invalid duplicates")
    }
    if(!is.null(query_gr$name)){
        if(!any(duplicated(query_gr$name))){
            message("using existing name as id")
            query_gr$id = query_gr$name
            names(query_gr) = NULL
            query_gr$name = NULL
            return(query_gr)
        }
        message("ignoring name attribute with invalid duplicates")
    }
    if(!is.null(names(query_gr))){
        if(!any(duplicated(names(query_gr)))){
            message("using existing name as id")
            query_gr$id = names(query_gr)
            names(query_gr) = NULL
            query_gr$name = NULL
            return(query_gr)
        }
        message("ignoring names() with invalid duplicates")
    }
    message("generating unique ids")
    names(query_gr) = NULL
    query_gr$name = NULL
    query_gr$id = paste(id_prefix, seq(length(query_gr)), sep = "_")
    query_gr
}

#' applies normalizations and transformations to prof_dt, a tidy data.table of
#' profiles.
#'
#' @param prof_dt data.table of ChIP-seq signal profiles.
#' @param norm_dt data.table containing norm_factor values for cell/mark combinations.
#' @param qgr GRanges
#' @param cap_value numeric, ChIP-seq data is prone to outliers, which will
#' wash out weaker signal if not properly capped.  Default is Inf, i.e. no capping.
#' @param high_on_right boolean, should profiles be flipped when most signal is
#'   on the left? Default is TRUE.  This is appropriate when there is no
#'   biological significance to the orientation of a signal, i.e. mirror image
#'   profiles are equivalent.
#'
#' @return prepared version of prof_dt
#' @export
#'
#' @examples
#'
prep_profile_dt = function(prof_dt,
                           norm_dt,
                           qgr,
                           cap_value = Inf,
                           high_on_right = TRUE){
    if(!all(norm_dt$norm_factor == 1)){
        prof_dt = merge(prof_dt, norm_dt)
        prof_dt[, y := y * norm_factor]
        prof_dt$norm_factor = NULL
    }
    if(is.finite(cap_vale)){
        prof_dt[y > cap_value, y := cap_value]
    }


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
#' @importFrom seqsetvis ssvFetchBigwig ssvFetchBam
#'
#' @examples
stsFetchTsneInput = function(qdt,
                             qgr,
                             region_size = NULL,
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
    qgr = prep_query_gr(qgr, region_size)

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
    # fetch tidy data
    bw_dt = fetch_FUN(qdt, qgr, qwin, qmet, agg_FUN, bfc)
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








