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
#' @param verbose if TRUE logs and status are output as messages
#' @param force_overwrite if TRUE, any contents of cache are overwritten.
#'
#' @return a tidy data.table containing t-sne embedding.  variable names
#' are tx, ty, id, and tall_var.
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom stats var
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#'
#' @examples
#' data("query_gr")
#' bw_files = dir(system.file('extdata', package = "chiptsne"), pattern = ".bw$", full.names = TRUE)
#' cfg_dt = data.table(file = bw_files)
#' cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' profile_dt = stsFetchTsneInput(cfg_dt, query_gr)
#' tsne_dt = stsRunTsne(profile_dt$bw_dt, perplexity = 15)
#' tsne_dt
stsRunTsne = function(profile_dt,
                      perplexity = 100,
                      n_cores = getOption("mc.cores", 1),
                      high_topright = TRUE,
                      norm1 = TRUE,
                      bfc = BiocFileCache::BiocFileCache(),
                      rname = NULL,
                      Y_init = NULL,
                      verbose = TRUE,
                      force_overwrite = FALSE){

    stopifnot(c("id", "tall_var", "wide_var", "x", "y") %in% colnames(profile_dt))
    # cast from tidy to wide matrix
    tsne_mat = dt2mat(profile_dt, unique(profile_dt$wide_var))
    bad_col = apply(tsne_mat, 2, stats::var) == 0
    if(any(bad_col)){
        stop("zero variance columns detected in tsne matrix input.",
             "\n", round(sum(bad_col) / length(bad_col)*100, 2), "% of columns affected.")

    }
    if(is.null(rname)){
        nr = min(nrow(tsne_mat), 50)
        nc = min(ncol(tsne_mat), 50)
        sel_r = floor(seq(0, nr)/nr*(nrow(tsne_mat)-1)+1)
        sel_c = floor(seq(0, nc)/nc*(ncol(tsne_mat)-1)+1)
        seq(nrow(tsne_mat))/nr
        rname = digest::digest(list(
            tsne_mat[sel_r, sel_c],
            perplexity
        ))
    }
    if(is.data.table(Y_init)){
        Y_init = as.matrix(Y_init[, .(tx, ty)], rownames.value = paste(Y_init$id, Y_init$tall_var))
        Y_init = Y_init[rownames(tsne_mat),]
        stopifnot(ncol(Y_init) == 2)
    }else{
        stopifnot(is.null(Y_init))
    }


    res_tsne = bfcif(bfc, rname,
                     FUN = function(){
                         Rtsne::Rtsne(tsne_mat,
                                      Y_init = Y_init,
                                      num_threads = n_cores,
                                      perplexity = perplexity,
                                      check_duplicates = FALSE)
                     },
                     force_overwrite = force_overwrite, verbose = verbose)
    tdt = as.data.table(res_tsne$Y)
    colnames(tdt) = c("tx", "ty")
    tdt$rn = rownames(tsne_mat)
    tdt[, c("id", "tall_var") := tstrsplit(rn, " ", keep = seq(2))]

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




#' dt2mat
#'
#' casts a tidy data.table, as returned from fetch_bam_dt or fetch_bw_dt to a
#' wide matrix compatible with Rtsne
#'
#' @param prof_dt tidy data.table to cast to a matrix
#' @param wide_vars character vector of wide_vars to cast
#'
#' @return a wide matrix with dimensions of nrows = (length(unique(id)) * number
#'   of tall_vars) and ncols = (numbers of viewing bins * number of wide_vars)
#' @export
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' n_bins = 5
#' n_tall_vars = 3
#' n_wide_vars = 2
#' n_ids = 5
#' bins = (seq(n_bins)-1) / (n_bins-1)
#' tall_vars = paste("tall_var", LETTERS[seq(n_tall_vars)], sep = "_")
#' wide_vars = rev(paste("wide_var", rev(letters)[seq(n_wide_vars)], sep = "_"))
#' ids = paste("region", seq(n_ids), sep = "_")
#' dt = data.table(
#'     x = rep(bins, length(wide_vars)*length(tall_vars)),
#'     wide_var = rep(wide_vars, length(ids), each = length(bins)),
#'     tall_var = rep(tall_vars, length(ids), each= length(bins)*length(wide_vars)),
#'     id = rep(ids, each = length(bins)*length(wide_vars)*length(tall_vars))
#'     )
#' dt$y = runif(nrow(dt))
#' mat = dt2mat(dt, wide_vars)
#' mat
dt2mat = function(prof_dt, wide_vars){
    for(m in wide_vars){
        if(m == wide_vars[1]){
            dt = dcast( prof_dt[wide_var == m], id+tall_var~x, value.var = "y")
            wide_mat = as.matrix(dt[, -seq_len(2)])
            rn = paste(dt$id, dt$tall_var)
        }else{
            dt = dcast( prof_dt[wide_var == m], id+tall_var~x, value.var = "y")
            stopifnot(all(paste(dt$id, dt$tall_var) == rn))
            wide_mat = cbind(wide_mat, as.matrix(dt[, -seq_len(2)]))


        }
    }
    rownames(wide_mat) = rn
    wide_mat
}










