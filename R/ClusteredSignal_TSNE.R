#' ClusteredSignal_TSNE
#'
#' adds support for xy coordinate (tsne result) data to ClusteredSignal class.
#'
#' @slot xy_data data.table.
#'
#' @export
setClass("ClusteredSignal_TSNE",
         representation = list(
             xy_data = "data.table",
             xy_method = "character"
         ), contains = "ClusteredSignal")

ClusteredSignal_TSNE.from_ClusteredSignal = function(object, sts_parent){
    stopifnot("ClusteredSignal" %in% class(object))
    object@signal_data = setalloccol(object@signal_data) #hacky fix
    fun.aggregate = mean
    tsne_dt = switch(
        sts_parent@dimreduce_method,
        tsne = {
            run_tsne(object@signal_data,
                     sts_parent@perplexity,
                     y_var = ssvQC:::val2var[sts_parent@signal_config@cluster_value],
                     fun.aggregate = fun.aggregate,
                     norm1 = TRUE)

        },
        umap = {
            run_umap(object@signal_data,
                     y_var = ssvQC:::val2var[sts_parent@signal_config@cluster_value],
                     fun.aggregate = fun.aggregate,
                     norm1 = TRUE)
        }
    )
    method_used = switch(
        sts_parent@dimreduce_method,
        tsne = {
            "tsne"
        },
        umap = {
            "umap"
        }
    )

    new("ClusteredSignal_TSNE",
        signal_data =  object@signal_data,
        query_gr = object@query_gr,
        signal_var = object@signal_var,
        facet_var = object@facet_var,
        extra_var = object@extra_var,
        xy_data = tsne_dt,
        xy_method = method_used)
}

#' ClusteredSignal_TSNE.null
#'
#' @return A "null" ClusteredSignal_TSNE object containing no information.
#' @export
#'
#' @examples
#' ClusteredSignal_TSNE.null()
ClusteredSignal_TSNE.null = function(){
    new("ClusteredSignal_TSNE",
        signal_data =  data.table(),
        query_gr = GRanges(),
        signal_var = character(),
        facet_var = character(),
        extra_var = character())
}

setMethod("names", "ClusteredSignal_TSNE",
          function(x)
          {
              c("signal_data",
                "query_gr",
                "assignment_data",
                "query_gr.cluster_list",
                "signal_data.mean_per_cluster",
                "xy_data")

          })


setMethod("$", "ClusteredSignal_TSNE",
          function(x, name)
          {
              switch (name,
                      signal_data = x@signal_data,
                      query_gr = x@query_gr,
                      assignment_data = unique(x@signal_data[, .(id, cluster_id)]),
                      query_gr.cluster_list = {
                          assign_dt = x$assignment_data
                          qgr = x@query_gr[assign_dt$id]
                          split(qgr, assign_dt$cluster_id)
                      },
                      signal_data.mean_per_cluster = {
                          sig_var = x@signal_var
                          agg_dt = x@signal_data[, .(y = mean(get(sig_var))), c("x", "cluster_id", union(x@facet_var, x@extra_var))]
                          setnames(agg_dt, "y", sig_var)
                          agg_dt
                      },
                      xy_data = {
                          x@xy_data
                      }
              )
          })

setReplaceMethod("$", "ClusteredSignal_TSNE",
                 function(x, name, value)
                 {
                     warn_msg = "This assignment is not supported.  No effect."
                     warning(warn_msg)
                     NULL
                 })
