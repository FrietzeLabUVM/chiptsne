#' stsPlotClusterProfiles
#'
#' @param profile_dt Tidy data.table of profile information. As returned by seqsetvis::ssvFetchBam.
#' @param cluster_dt Tidy data.table of cluster information. As returned by chiptsne:::nn_clust or chiptsne:::combine_most_similar
#' @param cluster_ Variable name of cluster assignment in cluster_dt.
#' @param w Relative width of profile glpyhs as fraction of total plot range.
#' @param h Relative height of profile glpyhs as fraction of total plot range.
#'
#' @return a ggplot of TSNE point clusters overlayed with summary profiles.
#' @importFrom concaveman concaveman
#' @importFrom GGally glyphs
#' @importFrom seqsetvis safeBrew
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' setalloccol(tsne_dt)
#' setalloccol(profile_dt)
#' cluster_dt = chiptsne:::nn_clust(tsne_dt, nn = 5)
#'
#' debug(chiptsne:::stsPlotClusterProfiles)
#' chiptsne:::stsPlotClusterProfiles(profile_dt, cluster_dt)
#'
#' chiptsne:::stsPlotClusterProfiles(profile_dt, cluster_dt, w = .1, h = .1) +
#'   scale_color_manual(values = c("H3K4me3" = "forestgreen", "H3K27me3" = "firebrick"))
#'
#' cluster_dt.reduced = chiptsne:::combine_most_similar(cluster_dt, profile_dt, n_times = 3)
#' chiptsne:::stsPlotClusterProfiles(profile_dt, cluster_dt.reduced, w = .1, h = .1) +
#'   scale_color_manual(values = c("H3K4me3" = "forestgreen", "H3K27me3" = "firebrick"))
stsPlotClusterProfiles = function (profile_dt,
                                   cluster_dt = profile_dt,
                                   cluster_ = "cluster_id",
                                   w = 0.05,
                                   h = 0.05,
                                   x_var = "x",
                                   y_var = "y",
                                   id_var = "id",
                                   wide_var = "wide_var",
                                   tall_var = "tall_var") {
    setalloccol(profile_dt)
    stopifnot(cluster_ %in% colnames(cluster_dt))
    cent_dt = cluster_dt[, .(tx = median(tx), ty = median(ty)),
                         by = cluster_]
    if(is.null(profile_dt[[tall_var]])){
        set(profile_dt, j = tall_var, value = "none")
    }
    set(profile_dt, j = "tid", value = paste(profile_dt[[tall_var]], profile_dt[[id_var]]))

    if(is.null(cluster_dt[["tid"]])){
        if(is.null(cluster_dt[[tall_var]])){
            set(cluster_dt, j = "tid", value = cluster_dt[[id_var]])
        }else{
            set(cluster_dt, j = "tid", value = paste(cluster_dt[[tall_var]], cluster_dt[[id_var]]))
        }
    }

    cprof_dt = merge(profile_dt[, setdiff(colnames(profile_dt), cluster_),  with = FALSE],
                     unique(cluster_dt[, c("tid", cluster_), with = FALSE]),
                     by = "tid")
    cprof_dt = cprof_dt[, .(y = mean(y)), c(wide_var, cluster_, x_var)]
    cprof_dt = merge(cprof_dt, cent_dt, by = cluster_)
    cent_dt[, `:=`(xmin, tx - w/2)]
    cent_dt[, `:=`(xmax, tx + w/2)]
    cent_dt[, `:=`(ymin, ty - h/2)]
    cent_dt[, `:=`(ymax, ty + h/2)]
    glyph_dt = as.data.table(GGally::glyphs(cprof_dt,
                                            x_major = "tx",
                                            x_minor = x_var,
                                            y_major = "ty",
                                            y_minor = y_var,
                                            width = w,
                                            height = h))
    my_chull = function(x, y) {
        ch = concaveman::concaveman(cbind(x, y), concavity = 1.4,
                                    length_threshold = 0.005)
        list(tx = ch[, 1], ty = ch[, 2])
    }
    ch_dt = cluster_dt[, my_chull(tx, ty), by = c(cluster_)]
    cols = seqsetvis::safeBrew(length(unique(ch_dt[[cluster_]])))
    if (is.numeric(ch_dt[[cluster_]])) {
        ch_dt[[cluster_]] = factor(as.character(ch_dt[[cluster_]]),
                                   levels = as.character(sort(unique(ch_dt[[cluster_]]))))
    }

    set(glyph_dt, j = "group__", value = paste(glyph_dt$gid, glyph_dt[[wide_var]]))

    p_clust_big = ggplot() +
        geom_polygon(data = ch_dt,
                     mapping = aes_string(x = "tx",
                                          y = "ty",
                                          fill = cluster_)) +
        scale_fill_manual(values = cols) +
        geom_rect(data = cent_dt,
                  aes(xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax),
                  fill = "white",
                  color = "black") +
        geom_path(data = glyph_dt,
                  aes_string(x = "gx",
                             y = "gy",
                             group = "group__",
                             color = wide_var)) +
        labs(x = "tx", y = "ty") +
        coord_fixed()
    p_clust_big
}

#' combine_most similar
#'
#' combines 2 most similar clusters, will repeat up to n_times as long as all clusters are over min_dist apart.
#'
#' @param p_dt plot dt containing cluster information.
#' @param profile_dt profile data
#' @param n_times number of times to combine
#' @param min_dist don't combined once under this distance
#' @param cluster_ variable name of cluster assignment
#' @param new_cluster_ variable name for new cluster assignment. by default original clusters will be overwritten.
#'
#' @return p_dt with clusters combined
#' @import reshape2
#'
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' setalloccol(tsne_dt)
#' clust_res = chiptsne:::nn_clust(tsne_dt, nn = 5, return_plot = TRUE)
#' clust_res[[2]]
#' chiptsne:::combine_most_similar(clust_res[[1]], profile_dt, n_times = 3)
combine_most_similar = function(p_dt,
                                profile_dt,
                                n_times = 1,
                                min_dist = Inf,
                                cluster_ = "cluster_id",
                                new_cluster_ = cluster_){
    if(n_times < 1){
        return(p_dt)
    }
    p_dt = copy(p_dt)
    profile_dt[, tid := paste(tall_var, id)]
    if(new_cluster_ != cluster_){
        p_dt[[new_cluster_]] = p_dt[[cluster_]]
        # p_dt[, meta := as.character(get(cluster_))]
    }
    if(!is.character(p_dt[[new_cluster_]])){
        p_dt[[new_cluster_]] = as.character(p_dt[[new_cluster_]])
    }

    uniq = unique(p_dt[[new_cluster_]])
    uniq = uniq[grepl("meta", uniq)]
    if(length(uniq) == 0){
        max_meta = 0
    }else{
        max_meta = max(as.numeric(sub(".+ ", "", uniq)))
    }
    new_meta = paste("meta", max_meta + 1)
    stopifnot(!new_meta %in% uniq)
    prof_cn = setdiff(colnames(profile_dt), new_cluster_)
    clust_dt = merge(profile_dt[, prof_cn, with = FALSE], p_dt[, c("tid", "tx", "ty", new_cluster_), with = FALSE], by = "tid")
    agg_dt = clust_dt[, .(
        ymean = mean(y)
    ), c(new_cluster_, "wide_var", "x")]

    agg_w = dcast(agg_dt, paste0(new_cluster_, "~wide_var+x"), value.var = "ymean")
    agg_mat = as.matrix(agg_w[,-1])
    rownames(agg_mat) = agg_w[[new_cluster_]]

    agg_dist = dist(agg_mat)

    ddt = as.data.table(reshape2::melt(as.matrix(agg_dist)))

    ddt$Var1 = factor(ddt$Var1)
    ddt$Var2 = factor(ddt$Var2)
    ddt = ddt[as.numeric(Var2) > as.numeric(Var1)]

    if(!any(min(ddt$value) < min_dist)){
        return(p_dt)
    }

    ddt = ddt[order(value)]

    ddt$Var1 = as.character(ddt$Var1)
    ddt$Var2 = as.character(ddt$Var2)


    p_dt[get(new_cluster_) %in% c(ddt$Var1[1], ddt$Var2[1]), ][[new_cluster_]] = new_meta

    if(n_times > 1){
        p_dt = combine_most_similar(p_dt,
                                    profile_dt,
                                    n_times = n_times - 1,
                                    min_dist = min_dist,
                                    cluster_ = new_cluster_,
                                    new_cluster_ = new_cluster_)
    }
    # if(new_cluster_ != "meta"){
    #     p_dt[[new_cluster_]] = p_dt$meta
    #     p_dt$meta = NULL
    # }
    p_dt[[new_cluster_]] = factor(p_dt[[new_cluster_]])
    levels(p_dt[[new_cluster_]]) = seq_along(levels(p_dt[[new_cluster_]]))

    p_dt[]
}


#' nn_clust
#'
#' performs nearest neighbor clustering on tsne coordinates.
#'
#' @param tsne_res data.table with tsne results. (needs tx, ty, tall_var, id)
#' @param nsamp Downsample profile values down to this value. Default of Inf skips downsampling.
#' @param nn Number of neighbors to use for clustering.  Default of 100.
#' @param show_plot If TRUE, print plot to graphics device. Default is FALSE.
#' @param return_plot If TRUE, return plot in list with cluster results. Default is FALSE.
#'
#' @return data.table mapping id_var entries to clusters ("cluster_id")
#'
#' @importFrom RANN nn2
#' @importFrom Matrix Matrix
#' @importFrom igraph graph.adjacency simplify cluster_walktrap
#' @examples
#' data("profile_dt")
#' data("tsne_dt")
#' setalloccol(tsne_dt)
#' setalloccol(profile_dt)
#' clust_res = chiptsne:::nn_clust(tsne_dt, nn = 5, return_plot = TRUE)
#' clust_res$plot
#' clust_res$data
#'
#' #just data by default
#' chiptsne:::nn_clust(tsne_dt, nn = 5)
nn_clust = function (tsne_res,
                     nsamp = Inf,
                     nn = 100,
                     tall_var = "tall_var",
                     id_var = "id",
                     show_plot = FALSE,
                     return_plot = FALSE){
    setalloccol(tsne_res)
    valid_vars = c(ifelse(tall_var == "tall_none", character(), tall_var), id_var)
    valid_vars = valid_vars[!is.na(valid_vars)]
    stopifnot(valid_vars %in% colnames(tsne_res))

    if(is.null(tsne_res[[tall_var]])){
        set(tsne_res, j = tall_var, value = "none")
    }

    if(nn > .2 * nrow(tsne_res)){
        nn = floor(.2 * nrow(tsne_res))
        message("Decreasing nearest-neighbors to ", nn, ".  Original value was too high for dataset.")
    }
    set(tsne_res, j = "tid", value = paste(tsne_res[[tall_var]], tsne_res[[id_var]]))

    mat = t(as.matrix(tsne_res[, .(tx, ty)]))
    colnames(mat) = tsne_res$tid
    mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
    knn.info <- RANN::nn2(t(mat), k = nn)
    knn <- knn.info$nn.idx
    colnames(knn) = c("tid", paste0("V", seq(nn - 1)))
    knn = as.data.table(knn)
    mknn = melt(knn, id.vars = "tid")
    ADJ = Matrix::Matrix(0, ncol(mat), ncol(mat))
    ADJ[cbind(mknn$tid, mknn$value)] = 1
    rownames(ADJ) = colnames(mat)
    colnames(ADJ) = colnames(mat)
    g <- igraph::graph.adjacency(ADJ, mode = "undirected")
    g <- igraph::simplify(g)
    km <- igraph::cluster_walktrap(g)
    com <- km$membership
    names(com) <- km$names
    com_dt = data.table(tid = names(com), cluster_id = com)
    p_dt = merge(tsne_res, com_dt, by = "tid")
    p_dt$cluster_id = factor(p_dt$cluster_id)

    p = ggplot(p_dt, aes(x = tx, y = ty, color = as.character(cluster_id))) +
        labs(color = "cluster_id") +
        geom_point(size = 1.5)
    if(show_plot) plot(p)
    if(return_plot){
        return(list(data = p_dt, plot = p))
    }else{
        return(p_dt)
    }
}
