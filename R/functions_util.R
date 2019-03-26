#' Title
#'
#' @param bfc
#' @param rname
#' @param FUN
#' @param force_overwrite
#' @param check_only
#'
#' @return
#' @export
#'
#' @examples
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE, check_only = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        message("results do exist.")
        if(check_only){
            return(TRUE)
        }
        message("loading results.")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        if(!file.exists(cache_path)){
            message("results do not exists.")
        }else{
            message("results being overwritten.")
        }
        if(check_only){
            return(FALSE)
        }
        message("executing function...")
        res = FUN()
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

#' Title
#'
#' @param gr
#' @param ref_gr
#' @param gr_size
#'
#' @return
#' @export
#'
#' @examples
my_annotate = function(gr,
                       ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz",
                                                        feature.type = "transcript", format = "gtf"),
                       gr_size = 2000){
    max_dist = 5e3
    dists = distanceToNearest(ref_gr, resize(gr, gr_size, fix = "center"))
    dists = subset(dists, distance <= max_dist)
    anno_dt = cbind(as.data.table(ref_gr[queryHits(dists)])[, .(gene_name, gene_id, transcript_id)],
                    as.data.table(gr[subjectHits(dists)]), distance = mcols(dists)[[1]])
    anno_dt
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
symbol2uniprot = function(x){
    bres = bitr(x, fromType = "SYMBOL",
                toType = c("UNIPROT"),
                OrgDb = org.Hs.eg.db)
    bres[!duplicated(bres$SYMBOL),]$UNIPROT
}

#' Title
#'
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @examples
sampleCap = function(x, n = 500){
    n = min(n, length(unique(x)))
    out = sample(unique(x), n)
    if(is.factor(out)) out = as.character(out)
    out
}

#' Title
#'
#' @param x
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
norm1 = function(x, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    (x - min(xrng)) / (max(xrng) - min(xrng))
}

#' Title
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
mybin = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(norm1(x, xrng) * (n_points-.00001))+1
}

#' Title
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#' @export
#'
#' @examples
mybin_centers = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    # xrng = range(x)
    xspc = diff(xrng)/n_points/2
    xs = seq(min(xrng)+xspc, max(xrng)-xspc, diff(xrng)/(n_points))
    xs
}
