#' bfcif
#'
#' conditionally run a function or load cached results based on rname.
#'
#' @param bfc a BiocFileCache object, as from running \code{BiocFileCache::BiocFileCache()}
#' @param rname character rname to save/load from cache.  \code{digest::digest(list(YOUR_ARGS))} works.
#' @param FUN function with no args to run.
#' @param force_overwrite if TRUE, FUN will be run and cached results replaced if present.
#' @param return_path_only if TRUE, FUN will never be run but path to file will be returned.  Particularly useful if you need to predetermined caching and paths for parallel execution.
#' @param verbose if TRUE, a series of status messages are output.  Default is FALSE.
#' @return the result of running FUN or loading cached results.
#' @export
#' @import BiocFileCache
#'
#' @examples
#' library(BiocFileCache)
#' val = 1
#' my_fun1 = function(){message(val); val}
#' bfc = BiocFileCache()
#' #basic usage
#' bfcif(bfc, "fun_results_1", my_fun1, verbose = TRUE)
#' #note how function scoping works
#' #be sure to update rname!
#' val = 2
#' bfcif(bfc, "fun_results_2", my_fun1, verbose = TRUE)
#' #but cached results haven't changed
#' bfcif(bfc, "fun_results_1", my_fun1, verbose = TRUE)
#'
#' #alternatively, we could just get the filepath
#' bfcif(bfc, "fun_results_1", my_fun1, verbose = TRUE, return_path_only = TRUE)
#' #this is true if results aren't yet cached
#' bfcif(bfc, "fun_results_3", my_fun1, verbose = TRUE, return_path_only = TRUE)
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE, return_path_only = FALSE, verbose = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        if(verbose) message("results not in cache. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        if(verbose) message("previous cache results found. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    if(return_path_only){
        if(verbose) message("returning cache path.")
        return(cache_path)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        if(verbose) message("loading previous cache results...")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        if(verbose) message("running function...", appendLF = FALSE)
        res = FUN()
        if(verbose) message("caching results...")
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

#' sampleCap
#'
#' Handy wrapper to sample() that avoid having to check against exceeding size
#' of x.
#'
#' @param x vector to sample from
#' @param n number of items to sample.  automatically reduced to length of x if
#'   too high.
#'
#' @return x sampled down to n items.
#' @export
#'
#' @examples
#' x = seq(10)
#' sampleCap(x, 5)
#' #avoid having to check number of items being sampled
#' sampleCap(x, 15)
sampleCap = function(x, n = 500){
    n = min(n, length(unique(x)))
    out = sample(unique(x), n)
    if(is.factor(out)) out = as.character(out)
    out
}

#' wraps scales::rescale to enforce range of \code{to} on output
#'
#' @param x continuous vector of values to manipulate.
#' @param to output range (numeric vector of length two)
#' @param from input range (vector of length two). If not given, is calculated from the range of x
#'
#' @return x rescaled from \code{from} domain to \code{to} domain, within bounds of \code{to}
#' @importFrom scales rescale
#' @examples
#' #behaves identically to scales::rescale when x is within 'from' domain
#' rescale_capped(0:10, to = c(0, 1), c(0, 10))
#' scales::rescale(0:10, to = c(0, 1), c(0, 10))
#' #when x exceeds 'from' domain, results are still within 'to' domain
#' rescale_capped(0:10, to = c(0,1), c(0,5))
#' #not true for scales::rescale
#' scales::rescale(0:10, to = c(0,1), c(0,5))
rescale_capped = function(x, to = c(0,1), from = range(x, na.rm = TRUE, finite = TRUE)){
    y = scales::rescale(x, to, from)
    y[y > max(to)] = max(to)
    y[y < min(to)] = min(to)
    y
}

#'
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#'
#' @examples
mybin = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    floor(rescale_capped(x, 0:1, xrng) * (n_points-.00001))+1
}

#' Title
#'
#' @param x
#' @param n_points
#' @param xrng
#'
#' @return
#'
#' @examples
mybin_centers = function(x, n_points, xrng = range(x)){
    stopifnot(length(xrng) == 2)
    # xrng = range(x)
    xspc = diff(xrng)/n_points/2
    xs = seq(min(xrng)+xspc, max(xrng)-xspc, diff(xrng)/(n_points))
    xs
}

# x: the vector n: the number of samples centered: if FALSE, then average
# current sample and previous (n-1) samples if TRUE, then average
# symmetrically in past and
# future. (If n is even, use one more sample from future.)
# from http://www.cookbook-r.com/Manipulating_data/Calculating_a_moving_average/
movingAverage = function(x, n = 1, centered = TRUE) {

    if (centered) {
        before <- floor((n - 1)/2)
        after <- ceiling((n - 1)/2)
    } else {
        before <- n - 1
        after <- 0
    }

    # Track the sum and count of number of non-NA items
    s <- rep(0, length(x))
    count <- rep(0, length(x))

    # Add the centered data
    new <- x
    # Add to count list wherever there isn't a
    count <- count + (!is.na(new))
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new

    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new <- c(rep(NA, i), x[seq_len(length(x) - i)])

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new <- c(x[(i + 1):length(x)], rep(NA, i))

        count <- count + (!is.na(new))
        new[is.na(new)] <- 0
        s <- s + new

        i <- i + 1
    }

    # return sum divided by count
    s/count
}

# adds a seriees of rectangles to ggplot p
# rects is matrix, 1 rectangle per row, columns of xmin, xmax, ymin, ymax
annotate_rects = function(p,
                          rects,
                          rect_fill = NA,
                          rect_color = "black",
                          text_color = "black",
                          text_position = c("topleft", "center", "none")
){
    anns = matrix(unlist(rects), ncol = 4, byrow = TRUE)
    p_alpha = .1

    for(i in seq_len(nrow(anns))){
        p = p +
            annotate("rect",
                     xmin = anns[i, 1],
                     xmax = anns[i, 2],
                     ymin = anns[i, 3],
                     ymax = anns[i, 4],
                     color = rect_color,
                     fill = rect_fill)
        if(!is.null(text_color)){
            switch(text_position,
                   topleft = {
                       p = p +
                           annotate("label",
                                    # x = (anns[i, 1] + anns[i, 2])/2,
                                    # y = (anns[i, 3] + anns[i, 4])/2,
                                    x = anns[i, 1],
                                    y = anns[i, 4],
                                    hjust = 1,
                                    vjust = 1,
                                    label = i,
                                    color = text_color)
                   },
                   center = {
                       p = p +
                           annotate("label",
                                    # x = (anns[i, 1] + anns[i, 2])/2,
                                    # y = (anns[i, 3] + anns[i, 4])/2,
                                    x = (anns[i, 1]+anns[i, 2])/2,
                                    y = (anns[i, 3]+anns[i, 4])/2,
                                    hjust = .5,
                                    vjust = .5,
                                    label = i,
                                    color = text_color)
                   }
            )
        }
    }
    p
}


#' Title
#'
#' @param tsne_res
#' @param cell_a
#' @param cell_b
#' @param n_points
#'
#' @return
#' @export
#'
#' @examples
calc_delta = function(tsne_res, cell_a, cell_b, n_points){
    v_dt = dcast(tsne_res[cell %in% c(cell_a, cell_b)], "id~cell", value.var = c("tx", "ty"))
    colnames(v_dt) = sub(cell_a, "cell_a", colnames(v_dt))
    colnames(v_dt) = sub(cell_b, "cell_b", colnames(v_dt))
    v_dt$bx_cell_a = mybin(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    xs = mybin_centers(v_dt$tx_cell_a, n_points = n_points, xrng = range(tsne_res$tx))
    v_dt$btx_cell_a = xs[v_dt$bx_cell_a]

    v_dt$by_cell_a = mybin(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    ys = mybin_centers(v_dt$ty_cell_a, n_points = n_points, xrng = range(tsne_res$ty))
    v_dt$bty_cell_a = ys[v_dt$by_cell_a]

    av_dt = v_dt[, list(tx_cell_b = mean(tx_cell_b), ty_cell_b = mean(ty_cell_b), N = .N), list(bx_cell_a, by_cell_a)]
    av_dt$tx_cell_a = xs[av_dt$bx_cell_a]
    av_dt$ty_cell_a = ys[av_dt$by_cell_a]
    return(list(velocity_dt = v_dt, agg_velocity_dt = av_dt))
}

#' Title
#'
#' @param data_dt
#' @param qgr
#' @param qcells
#' @param tss_ids
#'
#' @return
#' @export
#'
#' @examples
plot_profiles_selected = function(data_dt, qgr, qcells, tss_ids,
                                  color_mapping = c("H3K27me3" = "firebrick",
                                                    "H3K4me3" = "forestgreen")){
    # qcells = c("H7", "CD34", "Kasumi1", "mm1s", "Nalm6")
    # tss_ids = subset(qgr, gene_name == "RUNX1")$id[1]
    stopifnot(qcells %in% unique(data_dt$cell))
    if(!tss_ids %in% data_dt$id){
        tmp = unlist(strsplit(tss_ids, " "))
        if(length(tmp) > 1){
            tss_ids = tmp[1]
            tmp = as.numeric(tmp[-1])
            tss_ids = subset(qgr, gene_name %in% tss_ids)$id[tmp]
        }else{
            tss_ids = subset(qgr, gene_name %in% tss_ids)$id
        }

    }
    names(qgr) = qgr$id
    message(paste(as.character(qgr[tss_ids]), collapse = "\n"))
    plot_dt = data_dt[id %in% tss_ids & cell %in% qcells]
    plot_dt$cell = factor(plot_dt$cell, levels = qcells)
    plot_dt$id = factor(plot_dt$id, levels = tss_ids)

    p = ggplot(plot_dt, aes(x = x, ymin = 0, ymax = y, y = y, color = mark, fill = mark)) +
        facet_grid("cell~id", switch = "y") +
        geom_ribbon(alpha = .5) +
        geom_path(show.legend = FALSE) +
        theme_classic() +
        theme(strip.background = element_blank(), strip.placement = "outside",
              strip.text.y = element_text(angle = 180)) +
        labs(x = "", y = "") +
        scale_color_manual(values = color_mapping) +
        scale_fill_manual(values = color_mapping)
    p
}
