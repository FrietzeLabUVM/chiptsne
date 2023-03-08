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
             xy_method = "character",
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

#' ClusteredSignal_TSNE
#'
#' @param signal_profile_dt Tidy data.table containing profile information.  See output of seqsetvis::ssvFetchBam.
#' @param query_gr A GRanges containing regions to retrieve signal data at.
#' @param manual_assigned NYI but should allow manual cluster assignment.
#' @param nclust Number of k-means clusters to calculate. Default is 6.
#' @param signal_var Variable name for signal information to cluster upon in signal_profile_dt. Default is "y".
#' @param signal_var.within Variable name for ranking items within clusters.  The Default is the same as signal_var.
#' @param facet_var Variable that will eventually be used in heatmap facets.  Ensures it is preserved and not aggregated away.  Default is "name_split".
#' @param extra_var Any extra variables to preserve and avaoid aggregating away.
#' @param bfc BiocFileCache to use, uses default location otherwise.
#'
#' @return A ClusteredSignal_TSNE object containing clustering and TSNE information.
#' @export
#'
#' @examples
#' signal_profile_dt = seqsetvis::CTCF_in_10a_profiles_dt
#' setnames(signal_profile_dt, "sample", "name_split")
#' query_gr = seqsetvis::CTCF_in_10a_overlaps_gr
#' clust_sig = ClusteredSignal_TSNE(signal_profile_dt, query_gr)
ClusteredSignal_TSNE = function(signal_profile_dt,
                                query_gr,
                                manual_assigned = list(),
                                nclust = 6,
                                signal_var = "y",
                                signal_var.within = "y",
                                facet_var = "name_split",
                                extra_var = character()){
    query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)

    if(is.null(signal_profile_dt[[signal_var]])) stop("signal_var \"", signal_var, "\" not found in signal_data." )
    clust_dt = seqsetvis::ssvSignalClustering(signal_profile_dt, nclust = nclust, facet_ = "name_split", max_cols = Inf, max_rows = Inf, fill_ = signal_var)
    if(signal_var != signal_var.within){
        clust_dt = within_clust_sort(clust_dt = clust_dt, facet_ = "name_split", fill_ = signal_var.within)
    }

    if(length(manual_assigned) > 0){
        stop("manual_assigned NYI")
    }
    new("ClusteredSignal_TSNE",
        signal_data =  clust_dt,
        query_gr = query_gr,
        signal_var = signal_var,
        facet_var = facet_var,
        extra_var = extra_var)
}

#' ClusteredSignal_TSNE.fromConfig
#'
#' @return A ClusteredSignal_TSNE object containing clustering and TSNE information.
#' @export
#'
#' @examples
#' options(mc.cores = 10)
#' set.seed(0)
#' feature_config_file = system.file(package = "ssvQC", "extdata/ssvQC_peak_config.csv")
#' feature_config = QcConfigFeatures.parse(feature_config_file, process_features= TRUE)
#'
#' bam_config_file = system.file(package = "ssvQC", "extdata/ssvQC_bam_config.csv")
#' bam_config = QcConfigSignal.parse(bam_config_file)
#'
#' sqc = ssvQC(feature_config, bam_config)
#'
#' query_gr = feature_config$assessment_features$CTCF_features
#' sclust = ClusteredSignal_TSNE.fromConfig(sqc$signal_config, sqc$features_config$assessment_features)
#' sclust$
ClusteredSignal_TSNE.fromConfig = function(signal_config,
                                           query_gr,
                                           manual_assigned = list(),
                                           nclust = 6,
                                           facet_var = "name_split",
                                           extra_var = character(),
                                           bfc = ssvQC:::new_cache()){
    is_bam = grepl("bam", signal_config@read_mode)

    ssvQC::bfcif(bfc, digest_args(), function(){
        if(is_bam){
            if(signal_config@cluster_value == "RPM" | signal_config@sort_value == "RPM"){
                if(is.null(signal_config@meta_data$mapped_reads)){
                    stop("Call ssvQC.prepMappedReads() on signal_config first.")
                }
            }

            if(signal_config@cluster_value == "linearQuantile" | signal_config@sort_value == "linearQuantile"){
                if(is.null(signal_config@meta_data$cap_value)){
                    stop("Call ssvQC.prepCapValue() on signal_config first.")
                }
            }

            if(signal_config@cluster_value == "RPM_linearQuantile" | signal_config@sort_value == "RPM_linearQuantile"){
                if(is.null(signal_config@meta_data$RPM_cap_value)){
                    stop("Call ssvQC.prepCapValue() on signal_config first.")
                }

            }

            query_gr = seqsetvis::prepare_fetch_GRanges_names(query_gr)

            fetch_res = fetch_signal_at_features(signal_config, query_gr, bfc)
            prof_dt = fetch_res$prof_dt
            query_gr = fetch_res$query_gr

            if(!is.null(prof_dt$mapped_reads)){
                prof_dt[, y_RPM := y / mapped_reads * 1e6]
            }
            if(!is.null(prof_dt$RPM_cap_value)){
                prof_dt[, y_RPM_linQ := y_RPM / RPM_cap_value ]
                prof_dt[y_RPM_linQ > 1, y_RPM_linQ := 1 ]
            }

            if(!is.null(prof_dt$cap_value)){
                prof_dt[, y_linQ := y / cap_value]
                prof_dt[y_linQ > 1, y_linQ := 1]
            }

            clust_dt = ClusteredSignal_TSNE(prof_dt, query_gr,
                                            manual_assigned = manual_assigned,
                                            nclust = nclust,
                                            signal_var = ssvQC:::val2var[signal_config@cluster_value],
                                            signal_var.within = ssvQC:::val2var[signal_config@sort_value],
                                            facet_var = facet_var,
                                            extra_var = extra_var)
        }else{
            stop("NYI")
        }
        list(clust_dt = clust_dt, query_gr = query_gr)
    })
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
