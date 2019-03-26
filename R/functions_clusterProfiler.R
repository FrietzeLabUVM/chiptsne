#' Title
#'
#' @param gene_lists
#' @param bg_genes
#' @param force_overwrite
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_kegg_fromGenes = function(gene_lists, bg_genes = NULL, force_overwrite = FALSE){
    gene_lists = lapply(gene_lists, symbol2uniprot)
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(gene_lists))
    }

    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "kegg"))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "kegg", "plot"))

    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite,
                function(){
                    message("calc KEGG res")
                    tryCatch(
                        expr = {
                            ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite,
                                       function(){
                                           message("calc compareCluster")
                                           compareCluster(geneCluster = gene_lists,
                                                          universe      = bg_genes,
                                                          fun = "enrichKEGG",
                                                          # OrgDb = org.Hs.eg.db,
                                                          keyType       = 'uniprot',
                                                          # ont           = "BP",
                                                          pAdjustMethod = "BH",
                                                          pvalueCutoff  = 0.05,
                                                          qvalueCutoff  = 0.1)
                                       })
                            # ck = my_clusterProfiler(qgr, clust_assign)
                            p = ck %>% dotplot
                            list(ck, p)
                        }, error = {
                            function(e){
                                ck = NULL
                                p = ggplot() + annotate("text", x = 0, y = 0, label ="no KEGG results")
                                list(ck, p)
                            }
                        })
                })
    res
}

#' Title
#'
#' @param qgr
#' @param clust_assign
#' @param bg_genes
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_kegg = function(qgr, clust_assign, bg_genes = NULL){
    peak_dt = as.data.table(qgr)
    peak_dt$name = names(qgr)
    mdt = merge(peak_dt, clust_assign, by.x = "name", by.y = "id")
    anno_dt = my_annotate(GRanges(mdt), ref_gr)

    anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
    anno_lists = anno_dt$V1
    names(anno_lists) = paste0("c", anno_dt$cluster_id)
    anno_lists = anno_lists[order(anno_dt$cluster_id)]
    names(anno_lists)
    lengths(anno_lists)
    # each cluster vs total genes in heatmap
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(anno_lists))
    }

    # # each cluster vs all marked by k4me3
    # bg_dt = my_annotate(peaks$H7_H3K4ME3)
    # # bg_dt = my_annotate(k4_gr$`CD34-01517_H3K4me3`)
    # bg_genes = unique(bg_dt$gene_name)

    gene_lists = lapply(anno_lists, symbol2uniprot)

    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "kegg"))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "kegg", "plot"))

    res = bfcif(bfc, rname_go_plot,
                force_overwrite = TRUE,
                function(){
                    message("calc KEGG res")
                    tryCatch(
                        expr = {
                            ck = bfcif(bfc, rname_go_dat,
                                       force_overwrite = TRUE, function(){
                                           message("calc compareCluster")
                                           compareCluster(geneCluster = gene_lists,
                                                          universe      = symbol2uniprot(bg_genes),
                                                          fun = "enrichKEGG",
                                                          # OrgDb = org.Hs.eg.db,
                                                          keyType       = 'uniprot',
                                                          # ont           = "BP",
                                                          pAdjustMethod = "BH",
                                                          pvalueCutoff  = 0.05,
                                                          qvalueCutoff  = 0.1)
                                       })
                            # ck = my_clusterProfiler(qgr, clust_assign)
                            p = ck %>% dotplot
                            list(ck, p)
                        }, error = {
                            function(e){
                                ck = NULL
                                p = ggplot() + annotate("text", x = 0, y = 0, label ="no KEGG results")
                                list(ck, p)
                            }
                        })
                })
    res
}

library(clusterProfiler)
library(org.Hs.eg.db)

#' Title
#'
#' @param gene_lists
#' @param bg_genes
#' @param force_overwrite
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_fromGenes = function(gene_lists, bg_genes = NULL, force_overwrite = FALSE){
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(gene_lists))
    }
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "BP"))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "BP", "plot"))
    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite, function(){
        message("calc go res")
        tryCatch(
            expr = {
                ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite, function(){
                    message("calc compareCluster")
                    compareCluster(geneCluster = gene_lists,
                                   universe      = bg_genes,
                                   fun = "enrichGO",
                                   OrgDb = org.Hs.eg.db,
                                   keyType       = 'SYMBOL',
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05)
                })
                # ck = my_clusterProfiler(qgr, clust_assign)
                p = ck %>% simplify %>% dotplot
                list(ck, p)
            }, error = {
                function(e){
                    ck = NULL
                    p = ggplot() + annotate("text", x = 0, y = 0, label ="no GO results")
                    list(ck, p)
                }
            })
    })
    res
}

#' Title
#'
#' @param qgr
#' @param clust_assign
#' @param bg_genes
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler = function(qgr, clust_assign, bg_genes = NULL){
    peak_dt = as.data.table(qgr)
    peak_dt$name = names(qgr)
    mdt = merge(peak_dt, clust_assign, by.x = "name", by.y = "id")
    anno_dt = my_annotate(GRanges(mdt), ref_gr)

    anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
    anno_lists = anno_dt$V1
    names(anno_lists) = paste0("c", anno_dt$cluster_id)
    anno_lists = anno_lists[order(anno_dt$cluster_id)]
    names(anno_lists)
    lengths(anno_lists)
    # each cluster vs total genes in heatmap
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(anno_lists))
    }

    # # each cluster vs all marked by k4me3
    # bg_dt = my_annotate(peaks$H7_H3K4ME3)
    # # bg_dt = my_annotate(k4_gr$`CD34-01517_H3K4me3`)
    # bg_genes = unique(bg_dt$gene_name)

    gene_lists = anno_lists

    my_clusterProfiler_fromGenes(gene_lists = gene_lists, bg_genes = bg_genes)
}

#' Title
#'
#' @param cp_res
#' @param clust_id
#'
#' @return
#' @export
#'
#' @examples
cluster_table = function(cp_res, clust_id = NULL){
    cdat = cp_res@compareClusterResult
    cdat = as.data.table(cdat)

    if(is.null(clust_id)){
        clust_id = levels(cdat$Cluster)[1]
    }
    if(!clust_id %in% levels(cdat$Cluster)){
        stop("clust_id (", clust_id, ") must be one of ", paste(levels(cdat$Cluster), collapse = ", "))
    }
    if(!clust_id %in% unique(cdat$Cluster)){
        warning("clust_id (", clust_id, ") is valid but has no significant enrichment")
        return(datatable(matrix("empty")))
    }

    gs_key = unique(cdat[, .(ID, Description)])
    setkey(gs_key, ID)

    col_order = cdat[Cluster == clust_id][order(qvalue)]$Description

    if(cp_res@fun == "enrichGO"){
        cmat = cdat[, .(gene_name = tstrsplit(geneID, "/")), by = .(ID,Cluster)]
        cmat$gene_name = unlist(cmat$gene_name)
        ctab =  dcast(cmat[Cluster == clust_id],
                      "ID~gene_name",
                      value.var = "gene_name")#, value.var = length)
    }else if(cp_res@fun == "enrichKEGG"){
        cmat = cdat[, .(uniprot = tstrsplit(geneID, "/")), by = .(ID,Cluster)]
        cmat$uniprot = unlist(cmat$uniprot)
        ctab =  dcast(cmat[Cluster == clust_id],
                      "ID~uniprot",
                      value.var = "gene_name")#, value.var = length)
    }


    dmat = as.matrix(ctab[,-1])
    dmat = ifelse(is.na(dmat), 0, -1)
    rownames(dmat) = gs_key[.(ctab$ID)]$Description

    if(cp_res@fun == "enrichGO"){

    }else if(cp_res@fun == "enrichKEGG"){
        colnames(dmat) = bitr(colnames(dmat), fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL
    }

    dmat = t(dmat)
    dmat = dmat[, col_order]
    dmat = dmat[order(rowSums(dmat)),]
    library(DT)
    datatable(dmat, options = list(pageLength = nrow(dmat))) %>% formatStyle(
        T,
        backgroundColor = styleEqual(c(0, -1), c('white', 'blue')),
        color = "00000000"

    )
}

