# functions to prepare data for running tsne and running tsne

#' stsRunTsne
#'
#' @param profile_dt a tidy data.table for profile data as retrieved by
#' stsFetchTsneInput
#' @param perplexity perplexity value to use when running Rtsne.
#' @param n_cores number of cores to use when running Rtsne.  Default is
#' value of mc.cores option if set or 1.
#' @param high_topright if TRUE, returned t-sne space is flipped horizontallly
#' and vertically as required to ensure maximum signal is in the top-right of
#' the space.
#' @param norm1 if TRUE, returned t-sne space is transormed to a domain of
#' [-.5, .5] in x and y.
#'
#' @param bfc BiocFileCache obect to use when caching.
#' @param rname rname to use when querying the cache. Default is a digest
#' of 20 sampled rows and columns from the tsne input matrix and perplexity.
#' @param force_overwrite if TRUE, any contents of cache are overwritten.
#'
#' @return a tidy data.table containing t-sne embedding.  variable names
#' are tx, ty, id, and cell.
#' @export
#' @importFrom Rtsne Rtsne
#' @examples
#' data("query_gr")
#' bw_files = dir(system.file('extdata', package = "seqtsne"), pattern = ".bw$", full.names = TRUE)
#' cfg_dt = data.table(file = bw_files)
#' cfg_dt[, c("cell", "mark") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[cell %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(mark == "H3K4me3", .3, 1)]
#' profile_dt = stsFetchTsneInput(cfg_dt, query_gr)
#' tsne_res = stsRunTsne(profile_dt$bw_dt, perplexity = 15)
#' tsne_res
#'
stsRunTsne = function(profile_dt,
                      perplexity = 100,
                      n_cores = getOption("mc.cores", 1),
                      high_topright = TRUE,
                      norm1 = TRUE,
                      bfc = BiocFileCache::BiocFileCache(),
                      rname = NULL,
                      force_overwrite = FALSE){

    stopifnot(c("id", "cell", "mark", "x", "y") %in% colnames(profile_dt))
    # cast from tidy to wide matrix
    tsne_mat = dt2mat(profile_dt, unique(profile_dt$mark))
    bad_col = apply(tsne_mat, 2, var) == 0
    if(any(bad_col)){
        stop("zero variance columns detected in tsne matrix input.",
             "\n", round(sum(bad_col) / length(bad_col)*100, 2), "% of columns affected.")

    }
    set.seed(0)
    if(is.null(rname)){
        nr = min(nrow(tsne_mat), 20)
        nc = min(ncol(tsne_mat), 20)
        rname = digest::digest(list(
            tsne_mat[sample(1:nrow(tsne_mat), nr),
                     sample(1:ncol(tsne_mat), nc)],
            perplexity
        ))
    }
    res_tsne = bfcif(bfc, rname,
                     FUN = function(){
                         Rtsne::Rtsne(tsne_mat,
                                      num_threads = n_cores,
                                      perplexity = perplexity,
                                      check_duplicates = FALSE)
                     },
                     force_overwrite = force_overwrite)
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
    tdt$rn = NULL
    tdt[]
}


#' applies normalizations and transformations to prof_dt, a tidy data.table of
#' profiles.
#'
#' @param prof_dt data.table of ChIP-seq signal profiles.
#' @param norm_dt data.table containing norm_factor values for cell/mark
#'   combinations.
#' @param qgr GRanges
#' @param cap_value numeric, ChIP-seq data is prone to outliers, which will wash
#'   out weaker signal if not properly capped.  Default is Inf, i.e. no capping.
#' @param high_on_right boolean, should profiles be flipped when most signal is
#'   on the left? Default is TRUE.  This is appropriate when there is no
#'   biological significance to the orientation of a signal, i.e. mirror image
#'   profiles are equivalent.
#'
#' @return list of two items.  prepared version of prof_dt and query_gr modified
#'   to reflect any flipping required by high_on_right.
#' @export
#' @examples
#' data(profile_datatable)
#' data(query_gr)
#' #typically, norm_dt is the same configuration table used to fetch prof_dt
#' #here we derive a new norm_dt that will reduce H3K4me3 to 30% of H3K27me3.
#' norm_dt = unique(profile_datatable[, list(cell, mark)])
#' norm_dt[, norm_factor := ifelse(mark == "H3K4me3", .3, 1)]
#' prep_profile_dt(profile_datatable, norm_dt, query_gr)
prep_profile_dt = function(prof_dt,
                           norm_dt,
                           qgr,
                           cap_value = Inf,
                           high_on_right = TRUE){
    if(!all(norm_dt$norm_factor == 1)){
        prof_dt = merge(prof_dt, norm_dt, by = intersect(colnames(prof_dt),
                                                         colnames(norm_dt)))
        prof_dt[, y := y * norm_factor]
        prof_dt$norm_factor = NULL
    }
    if(high_on_right){
        balance_dt = prof_dt[, list(right_sum = sum(y[x > 0]),
                                    left_sum = sum(y[x < 0])),
                             by = list(cell, id)]
        balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum,
                                       cell,
                                       id)]
        most_flipped = balance_dt[,
                                  list(fraction_flipped = sum(needs_flip) / .N),
                                  by = list(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        GenomicRanges::strand(qgr) = "+"
        GenomicRanges::strand(qgr)[most_flipped$flip_strand] = "-"
        prof_dt = merge(prof_dt, balance_dt, by = c("id", "cell"))
        remove(balance_dt)
        prof_dt[needs_flip == TRUE, x := -x]
        prof_dt$needs_flip = NULL
    }
    if(is.finite(cap_value)){
        prof_dt[y > cap_value, y := cap_value]
    }
    list(prof_dt = prof_dt[], query_gr = qgr)
}

#' dt2mat
#'
#' casts a tidy data.table, as returned from fetch_bam_dt or fetch_bw_dt to a
#' wide matrix compatible with Rtsne
#'
#' @param prof_dt tidy data.table to cast to a matrix
#' @param marks character vector of marks to cast
#'
#' @return a wide matrix with dimensions of nrows = (length(unique(id)) * number
#'   of cells) and ncols = (numbers of viewing bins * number of marks)
#' @export
#'
#' @examples
#' n_bins = 5
#' n_cells = 3
#' n_marks = 2
#' n_ids = 5
#' bins = (seq(n_bins)-1) / (n_bins-1)
#' cells = paste("cell", LETTERS[seq(n_cells)], sep = "_")
#' marks = rev(paste("mark", rev(letters)[seq(n_marks)], sep = "_"))
#' ids = paste("region", seq(n_ids), sep = "_")
#' dt = data.table(
#'     x = rep(bins, length(marks)*length(cells)),
#'     mark = rep(marks, length(ids), each = length(bins)),
#'     cell = rep(cells, length(ids), each= length(bins)*length(marks)),
#'     id = rep(ids, each = length(bins)*length(marks)*length(cells))
#'     )
#' dt$y = runif(nrow(dt))
#' mat = dt2mat(dt, marks)
#' mat
dt2mat = function(prof_dt, marks){
    for(m in marks){
        if(m == marks[1]){
            dt = dcast( prof_dt[mark == m], id+cell~x, value.var = "y")
            wide_mat = as.matrix(dt[, -1:-2])
            rn = paste(dt$id, dt$cell)
        }else{
            dt = dcast( prof_dt[mark == m], id+cell~x, value.var = "y")
            stopifnot(all(paste(dt$id, dt$cell) == rn))
            wide_mat = cbind(wide_mat, as.matrix(dt[, -1:-2]))


        }
    }
    rownames(wide_mat) = rn
    wide_mat
}










