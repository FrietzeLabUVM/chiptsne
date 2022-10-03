#' Title
#'
#' @param profile_dt
#' @param cluster_dt
#' @param cluster_
#' @param w
#' @param h
#'
#' @return
#' @export
#' @importFrom concaveman concaveman
#' @importFrom GGally glyphs
#' @importFrom seqsetvis safeBrew
#' @examples
stsPlotClusterProfiles = function(profile_dt, cluster_dt,
                                  cluster_ = "cluster_id",
                                  w = .05,
                                  h = .05){
    stopifnot(cluster_ %in% colnames(cluster_dt))
    cent_dt = cluster_dt[, .(tx = median(tx), ty = median(ty)), by = cluster_]

    profile_dt[, tid := paste(tall_var, id)]
    cprof_dt = merge(profile_dt, cluster_dt[, c("tid", cluster_), with = FALSE], by = "tid")
    cprof_dt = cprof_dt[, .(y = mean(y)), c("wide_var", cluster_, "x")]
    cprof_dt = merge(cprof_dt, cent_dt, by = cluster_)


    cent_dt[, xmin := tx - w / 2]
    cent_dt[, xmax := tx + w / 2]
    cent_dt[, ymin := ty - h / 2]
    cent_dt[, ymax := ty + h / 2]

    glyph_dt = as.data.table(GGally::glyphs(cprof_dt,
                                            x_major = "tx",
                                            x_minor = "x",
                                            y_major = "ty",
                                            y_minor = "y",
                                            width = w, height = h))

    my_chull = function(x, y){
        # ch = chull(x, y)
        # list(x[ch], y[ch])
        ch = concaveman::concaveman(cbind(x, y),
                                    concavity = 1.4, length_threshold = .005)
        list(tx = ch[,1], ty = ch[,2])
    }
    ch_dt = cluster_dt[, my_chull(tx, ty), by = c(cluster_)]
    cols = seqsetvis::safeBrew(length(unique(ch_dt[[cluster_]])))
    if(is.numeric(ch_dt[[cluster_]])){
        ch_dt[[cluster_]] = factor(
            as.character(ch_dt[[cluster_]]),
            levels = as.character(sort(unique(ch_dt[[cluster_]])))
        )
    }

    p_clust_big = ggplot() +
        geom_polygon(data = ch_dt,
                     mapping = aes_string(
                         x = "tx",
                         y = "ty",
                         fill = cluster_
                     )) +
        scale_fill_manual(values = cols)

    p_clust_big +
        geom_rect(data = cent_dt,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "white", color = "black") +
        geom_path(data = glyph_dt,
                  aes(gx, gy, group = paste(gid, wide_var), color = wide_var)) +
        labs(x = "tx", y = "ty") +
        coord_fixed()
}

#' Title
#'
#' @param tsne_res
#' @param nsamp
#' @param nn
#'
#' @return
#' @export
#' @importFrom RANN nn2
#' @importFrom Matrix Matrix
#' @importFrom igraph graph.adjacency simplify cluster_walktrap
#' @examples
nn_clust = function(tsne_res, nsamp = Inf, nn = 100){
    tsne_res[, tid := paste(tall_var, id)]
    mat = t(as.matrix(tsne_res[, .(tx, ty)]))
    colnames(mat) = tsne_res$tid
    mat = mat[, sampleCap(seq(ncol(mat)), nsamp)]
    knn.info <- RANN::nn2(t(mat), k=nn)
    knn <- knn.info$nn.idx
    colnames(knn) = c("tid", paste0("V", seq(nn-1)))
    knn = as.data.table(knn)
    mknn = reshape2::melt(knn, id.vars = "tid")

    # adj <- matrix(0, ncol(mat), ncol(mat))
    # rownames(adj) <- colnames(adj) <- colnames(mat)
    # for(i in seq_len(ncol(mat))) {
    #     adj[i,colnames(mat)[knn[i,]]] <- 1
    # }
    ADJ = Matrix::Matrix(0, ncol(mat), ncol(mat))
    ADJ[cbind(mknn$tid, mknn$value)] = 1
    rownames(ADJ) = colnames(mat)
    colnames(ADJ) = colnames(mat)
    g <- igraph::graph.adjacency(ADJ, mode="undirected")
    g <- igraph::simplify(g) ## remove self loops
    # V(g)$color <- rainbow(G)[group[names(V(g))]] ## color nodes by group
    # plot(g, vertex.label=NA)
    km <- igraph::cluster_walktrap(g)
    ## community membership
    com <- km$membership
    names(com) <- km$names
    com_dt = data.table(tid = names(com), cluster_id = com)

    p_dt = merge(tsne_res, com_dt, by = "tid")
    # p_dt[, coms := paste("cluster", com)]
    # ggplot(p_dt, aes(x = tx, y = ty, color = coms)) +
    #     annotate("point", x  = p_dt$tx, y = p_dt$ty, size = .2) +
    #     geom_point(size = .5) +
    #     facet_wrap("coms")
    p = ggplot(p_dt, aes(x = tx, y = ty, color = as.character(cluster_id))) +
        labs(color = "cluster_id") +
        # annotate("point", x  = p_dt$tx, y = p_dt$ty, size = .2) +
        geom_point(size = .5) #+
    # facet_wrap("coms")
    return(list(data = p_dt, plot = p))
}

#' combine_mostsimilar
#'
#' @param p_dt
#' @param profile_dt
#' @param n_times
#' @param min_dist
#'
#' @return
#' @export
#' @import reshape2
#'
#' @examples
combine_mostsimilar = function(p_dt, profile_dt,
                               n_times = 1, min_dist = Inf,
                               cluster_ = "cluster_id",
                               new_cluster_ = cluster_){
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
    clust_dt = merge(profile_dt, p_dt[, c("tid", "tx", "ty", new_cluster_), with = FALSE], by = "tid")
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
        p_dt = combine_mostsimilar(p_dt,
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
    p_dt[]
}

