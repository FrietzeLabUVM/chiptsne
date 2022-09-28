#functions for fetchng and formatting genomic data for t-sne then running t-sne.

#' stsFetchTsneInput
#'
#' retrieve data from series of either bigwig or bam files to a tidy data.table
#' for later input to stsRunTsne
#'
#' validates inputs prepares query_gr fetches either bam or bigwig based on
#' extension of first file transform and normalizes profiles after fetch
#' @param qdt data.table containing, file, tall_var, wide_var, and optionally
#'   norm_factor
#' @param qgr GRanges of regions to fetch
#' @param region_size numeric, if provided will be used to set all widths on
#'   qgr.
#' @param file_format character, describing file format. must be one of bam or
#'   bw.  Will be automatically determined from file extension if NULL.
#' @param qwin number of datapoints to use to view each region
#' @param qmet strategy to use in each region, summary or sample
#' @param cap_value maximum allowed value.  useful for outlier handling.
#' @param fetch_FUN function used to load and process data from files.
#' @param agg_FUN function used to aggregate any duplicate tall_var+wide_var mappings in
#'   qdt.  Should accept a numeric vector and return a single value.
#' @param high_on_right if TRUE, profiles where highest signal is on the left
#'   are flipped.
#' @param bfc BiocFileCache object to use to cache data
#' @param n_cores number of cores to use. Defaults to value of mc.cores or 1 if
#'   mc.cores is not set.
#' @param rname rname to use with cache.  Default is a digest of arguments used
#'   to fetch data.
#' @param force_overwrite if TRUE, any contents of cache are overwritten. Default is FALSE.
#' @param verbose if TRUE outputs progress and log messages
#' @param skip_checks if TRUE, tall_var and wide_var completeness checks are skipped.
#'   tall_var/wide_var combinations must be complete to run StsRunTsne(). Default is FALSE.
#' @return a tidy data.table of profile data.
#' @export
#' @importFrom seqsetvis ssvFetchBigwig ssvFetchBam
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' data("query_gr")
#' bw_files = dir(system.file('extdata', package = "chiptsne"), pattern = ".bw$", full.names = TRUE)
#' cfg_dt = data.table(file = bw_files)
#' cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' profile_dt = stsFetchTsneInput(cfg_dt, query_gr)
#' profile_dt
stsFetchTsneInput = function(qdt,
                             qgr,
                             region_size = NULL,
                             file_format = NULL,
                             qwin = 50,
                             qmet = "summary",
                             cap_value = 20,
                             fetch_FUN = NULL,
                             agg_FUN = NULL,
                             high_on_right = TRUE,
                             bfc = BiocFileCache::BiocFileCache(),
                             n_cores = getOption("mc.cores", 1),
                             rname = digest::digest(list(qgr,
                                                         qdt,
                                                         qwin,
                                                         qmet,
                                                         fetch_FUN)),
                             force_overwrite = FALSE,
                             verbose = TRUE,
                             skip_checks = FALSE){
    stopifnot(is.data.frame(qdt))
    qgr = prep_query_gr(qgr, region_size)

    qdt = as.data.table(qdt)
    if(is.null(qdt$norm_factor)){
        qdt$norm_factor = 1
    }
    if(!"file" %in% colnames(qdt)){
        warning("attribute 'file' not found, assuming first column is 'file'.")
        colnames(qdt)[1] = "file"
    }
    if(!"tall_var" %in% colnames(qdt)){
        warning("attribute 'tall_var' not found, assuming second column is 'tall_var'.")
        colnames(qdt)[2] = "tall_var"
    }
    if(!"wide_var" %in% colnames(qdt)){
        warning("attribute 'wide_var' not found, assuming third column is 'wide_var'.")
        colnames(qdt)[3] = "wide_var"
    }
    if(!all(file.exists(qdt$file))){
        stop(paste(sep = "\n", "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt$file[!file.exists(qdt$file)])))

    }
    if(!skip_checks){
        if(!length(unique(qdt[, paste(sort(unique(wide_var)), collapse = ","),
                              by = .(tall_var)]$V1)) == 1){
            stop("Every wide_var should be mapped to every tall_var at least once.")
        }
        if(!length(unique(qdt[, paste(sort(unique(tall_var)), collapse = ","),
                              by = .(wide_var)]$V1)) == 1){
            stop("Every tall_var should be mapped to every wide_var at least once.")

        }
    }
    if(is.null(file_format)){#try to determine file format
        test_file = tolower(qdt$file[1])
        if(grepl(".bw$", test_file) | grepl(".bigwig$", test_file)){
            file_format = "bw"
        }else if(grepl(".bam$", test_file)){
            file_format = "bam"
        }else{
            stop("couldn't determine file_format of ", qdt$file[1],
                 "\nplease manually set file_format to bam or bw (bigWig).")
        }

    }
    stopifnot(file_format %in% c("bam", "bw"))
    if(is.null(fetch_FUN)){
        fetch_FUN = switch(file_format,
                           bam = {
                               if(verbose)
                                   message(
                                       "Using default fetch_FUN for bam data.",
                                       "\nfragLen will be automatically calculated",
                                       " and pileup will be sum of both strands."
                                   )
                               fetch_bam_dt
                           },
                           bw = {
                               if(verbose)
                                   message("Using default fetch_FUN for bigwig data.")
                               fetch_bw_dt
                           })
    }
    if(is.null(agg_FUN)){
        agg_FUN = switch(file_format,
                         bam = {
                             if(verbose)
                                 message("for bam input, choosing sum() ",
                                         "to aggregate any tall_var/wide_var duplicates.")
                             sum
                         },
                         bw = {
                             if(verbose)
                                 message("for bigwig input, choosing mean() ",
                                         "to aggregate any tall_var/wide_var duplicates.")
                             mean
                         })

    }



    # browser()
    # fetch tidy data
    bw_dt = fetch_FUN(qdt[, .(file)],
                      qgr,
                      qwin,
                      qmet,
                      agg_FUN,
                      bfc,
                      n_cores,
                      rname,
                      force_overwrite,
                      verbose)

    # apply transforms and thresholds
    prep_res = prep_profile_dt(bw_dt,
                               qdt,#unique(qdt[, list(file, tall_var, wide_var, norm_factor)]),
                               qgr,
                               cap_value = cap_value,
                               high_on_right = high_on_right)
    bw_dt = prep_res[[1]]
    qgr = prep_res[[2]]


    return(list(bw_dt = bw_dt, query_gr = qgr, rname = rname))
}

#' fetch_bam_dt
#'
#' fetches a tidy data.table of scores from bigwigs, automatically uses cache.
#'
#' @param qdt data.table containing, file, tall_var, wide_var, and optionally
#'   norm_factor
#' @param qgr GRanges of regions to fetch
#' @param qwin number of datapoints to use to view each region
#' @param qmet strategy to use in each region, summary or sample
#' @param agg_FUN function used to aggregate any duplicate tall_var+wide_var mappings in
#'   qdt. Should accept a numeric vector and return a single value.
#' @param bfc BiocFileCache object to use to cache data
#' @param n_cores number of cores to use. Defaults to value of mc.cores or 1 if
#'   mc.cores is not set.
#' @param rname rname to use with cache.  Default is a digest of arguments used
#'   to fetch data.
#' @param force_overwrite if TRUE, any contents of cache are overwritten.
#'
#' @return a tidy data.table of profile data.
#' @export
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom digest digest
#' @importFrom seqsetvis ssvFetchBam
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' data("query_gr")
#' bam_files = dir(system.file('extdata', package = "chiptsne"),
#'     pattern = ".bam$", full.names = TRUE)
#' cfg_dt = data.table(file = bam_files)
#' cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' profile_dt = fetch_bam_dt(cfg_dt, query_gr)
#' profile_dt
fetch_bam_dt = function(qdt,
                        qgr,
                        qwin = 50,
                        qmet = "summary",
                        agg_FUN = sum,
                        bfc = NULL,
                        n_cores = getOption("mc.cores", 1),
                        rname = NULL,
                        force_overwrite = FALSE,
                        verbose = TRUE,
                        fragLens = NULL){
    #how to handle stranded?
    #how to handle frag lens?
    if(is.null(bfc)){
        bfc = BiocFileCache::BiocFileCache()
    }
    if(is.null(rname)){
        rname = digest::digest(list(qgr,
                                    qdt,
                                    qwin,
                                    qmet))
    }
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(!all(file.exists(qdt$file))){
        missing_files = qdt$file[!file.exists(qdt$file)]
        stop(paste(sep = "\n",
                   "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", missing_files)))

    }

    bam_fetch = function(){
        seqsetvis::ssvFetchBam(qdt, qgr, names_variable = "file",
                               fragLens = fragLens,
                               return_data.table = TRUE,
                               win_method = qmet,
                               win_size = qwin,
                               n_cores = n_cores,
                               n_region_splits = 50)
    }
    bam_dt = bam_fetch()
    # bam_dt$sample = NULL
    # bam_dt = bam_dt[, list(y = agg_FUN(y)), list(tall_var, id, wide_var, x)]
    bam_dt
}

#' fetch_bam_stranded_dt
#'
#' fetches a tidy data.table of scores from bigwigs, automatically uses cache.
#' In contrast to fetch_bam_dt, fragment size is not calculated and fragments
#' are not extended and pileup per strand is stored.
#'
#' @param qdt data.table containing, file, tall_var, wide_var, and optionally
#'   norm_factor
#' @param qgr GRanges of regions to fetch
#' @param qwin number of datapoints to use to view each region
#' @param qmet strategy to use in each region, summary or sample
#' @param agg_FUN function used to aggregate any duplicate tall_var+wide_var mappings in
#'   qdt. Should accept a numeric vector and return a single value.
#' @param bfc BiocFileCache object to use to cache data
#' @param n_cores number of cores to use. Defaults to value of mc.cores or 1 if
#'   mc.cores is not set.
#' @param rname rname to use with cache.  Default is a digest of arguments used
#'   to fetch data.
#' @param force_overwrite if TRUE, any contents of cache are overwritten.
#'
#' @return a tidy data.table of profile data.
#' @export
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom digest digest
#' @importFrom seqsetvis ssvFetchBam
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' data("query_gr")
#' bam_files = dir(system.file('extdata', package = "chiptsne"),
#'     pattern = ".bam$", full.names = TRUE)
#' cfg_dt = data.table(file = bam_files)
#' cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' profile_dt = fetch_bam_stranded_dt(cfg_dt, query_gr)
#' profile_dt
fetch_bam_stranded_dt = function(qdt,
                                 qgr,
                                 qwin = 50,
                                 qmet = "summary",
                                 agg_FUN = sum,
                                 bfc = NULL,
                                 n_cores = getOption("mc.cores", 1),
                                 rname = NULL,
                                 force_overwrite = FALSE, verbose = TRUE){
    #how to handle stranded?
    #how to handle frag lens?
    if(is.null(bfc)){
        bfc = BiocFileCache::BiocFileCache()
    }
    if(is.null(rname)){
        rname = digest::digest(list(qgr,
                                    qdt,
                                    qwin,
                                    qmet))
    }
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(!all(file.exists(qdt$file))){
        missing_files = qdt$file[!file.exists(qdt$file)]
        stop(paste(sep = "\n",
                   "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", missing_files)))

    }

    bam_fetch = function(){
        seqsetvis::ssvFetchBam(qdt, qgr, names_variable = "file",
                               target_strand = "both",
                               fragLens = NA,
                               return_data.table = TRUE,
                               win_method = qmet,
                               win_size = qwin,
                               n_cores = n_cores,
                               n_region_splits = 50)
    }
    bam_dt = bfcif(bfc, rname, bam_fetch, force_overwrite = force_overwrite, verbose = verbose)
    # bam_dt$sample = NULL
    # bam_dt = bam_dt[, list(y = agg_FUN(y)), list(tall_var, id, wide_var, x)]
    bam_dt
}

#' fetch_bw_dt
#'
#' fetches a tidy data.table of read pileups from bams, automatically uses
#' cache.
#'
#' @param qdt data.table containing, file, tall_var, wide_var, and optionally
#'   norm_factor
#' @param qgr GRanges of regions to fetch
#' @param qwin number of datapoints to use to view each region
#' @param qmet strategy to use in each region, summary or sample
#' @param agg_FUN used to aggregate any duplicate tall_var+wide_var mappings in qdt.
#' @param bfc BiocFileCache object to use to cache data
#' @param n_cores number of cores to use. Defaults to value of mc.cores or 1 if
#'   mc.cores is not set.
#' @param rname rname to use with cache.  Default is a digest of arguments used
#'   to fetch data.
#' @param force_overwrite if TRUE, any contents of cache are overwritten.
#' @return a tidy data.table of bigwig profiles at qgr views
#' @export
#' @importFrom BiocFileCache BiocFileCache
#' @importFrom digest digest
#' @importFrom seqsetvis ssvFetchBigwig
#' @importFrom GenomicRanges width
#' @examples
#' data("query_gr")
#' bw_files = dir(system.file('extdata', package = "chiptsne"),
#'     pattern = ".bw$", full.names = TRUE)
#' cfg_dt = data.table(file = bw_files)
#' cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
#' cfg_dt = cfg_dt[tall_var %in% c("ESH1", "HUES48", "HUES64")]
#' cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' profile_dt = fetch_bw_dt(cfg_dt, query_gr)
#' profile_dt
fetch_bw_dt = function(qdt,
                       qgr,
                       qwin = 50,
                       qmet = "summary",
                       agg_FUN = mean,
                       bfc = NULL,
                       n_cores = getOption("mc.cores", 1),
                       rname = NULL,
                       force_overwrite = FALSE,
                       verbose = TRUE){
    if(is.null(bfc)){
        bfc = BiocFileCache::BiocFileCache()
    }
    if(is.null(rname)){
        rname = digest::digest(list(qgr,
                                    qdt,
                                    qwin,
                                    qmet))
    }
    if(!length(unique(GenomicRanges::width(qgr))) == 1){
        warning("Multiple widths for query GRanges, it is recommended",
                " that widths be uniform.  Try GenomicRanges::resize() with",
                ' fix = "center".')
    }
    stopifnot(is.data.frame(qdt))
    qdt = as.data.table(qdt)
    if(!all(file.exists(qdt$file))){
        stop(paste(sep = "\n",
                   "\nFirst variable of qdt must be valid file paths.",
                   "Files not found:",
                   paste(collapse = " ", qdt$file[!file.exists(qdt$file)])))

    }

    bw_fetch = function(){
        seqsetvis::ssvFetchBigwig(qdt, qgr, names_variable = "file",
                                  return_data.table = TRUE,
                                  win_method = qmet,
                                  win_size = qwin, n_cores = n_cores,
                                  n_region_splits = 50)
    }
    # browser()
    bw_dt = bfcif(bfc, rname, bw_fetch, force_overwrite = force_overwrite, verbose = verbose)
    # bw_dt = bw_dt[, list(y = agg_FUN(y)), list(tall_var, id, wide_var, x)]
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
#' library(GenomicRanges)
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
#' @param norm_dt data.table containing norm_factor values for tall_var/wide_var
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
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' data(profile_dt)
#' data(query_gr)
#' #typically, norm_dt is the same configuration table used to fetch prof_dt
#' #here we derive a new norm_dt that will reduce H3K4me3 to 30% of H3K27me3.
#' norm_dt = unique(profile_dt[, list(tall_var, wide_var)])
#' norm_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
#' prep_profile_dt(profile_dt, norm_dt, query_gr)
prep_profile_dt = function(prof_dt,
                           norm_dt,
                           qgr,
                           cap_value = Inf,
                           high_on_right = TRUE){
    x = y = norm_factor = tall_var = id = left_sum = right_sum =
        needs_flip = .N = flipe_strand = fraction_flipped = NULL
    prof_dt = merge(prof_dt, norm_dt, by = intersect(colnames(prof_dt),
                                                     colnames(norm_dt)))
    if(!all(norm_dt$norm_factor == 1)){
        prof_dt[, y := y * norm_factor]
        # prof_dt$norm_factor = NULL
    }
    if(high_on_right){
        balance_dt = prof_dt[, list(right_sum = sum(y[x > 0]),
                                    left_sum = sum(y[x < 0])),
                             by = list(tall_var, id)]
        balance_dt = balance_dt[, list(needs_flip = left_sum > right_sum,
                                       tall_var,
                                       id)]
        most_flipped = balance_dt[,
                                  list(fraction_flipped = sum(needs_flip) / .N),
                                  by = list(id)]
        most_flipped[, flip_strand := fraction_flipped > .5]
        GenomicRanges::strand(qgr) = "+"
        GenomicRanges::strand(qgr)[most_flipped$flip_strand] = "-"
        prof_dt = merge(prof_dt, balance_dt, by = c("id", "tall_var"))
        remove(balance_dt)
        prof_dt[needs_flip == TRUE, x := -x]
        prof_dt$needs_flip = NULL
    }
    if(is.finite(cap_value)){
        prof_dt[y > cap_value, y := cap_value]
    }
    list(prof_dt = prof_dt[], query_gr = qgr)
}
