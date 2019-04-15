#' Data from 3 hESC cell lines for two histone modifications, H3K4me3 and
#' H3K27me3.
#'
#' @description
#' This is ENCODE data that has been aligned to hg38 by STAR.  Bams
#' and fold-enrichment bigwigs are included for 30 randomly selected regions
#' from a published list of bivalent (H3K4me3 and H3K27me3, an active and
#' repressive mark co-occuring) sites in hESC cell lines.
#'
#' @details Contains: \itemize{
#' \item \code{\link{query_gr}}
#' }
#' Data from GEO series
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312}
#' \tabular{llll}{
#' ESH1 \tab input \tab GSM433179 \tab R1 \cr
#' ESH1 \tab H3K27me3 \tab GSM537683 \tab R1 \cr
#' ESH1 \tab H3K4me3 \tab GSM537681 \tab R1 \cr
#' ESH1 \tab input \tab GSM605335 \tab R2 \cr
#' ESH1 \tab H3K27me3 \tab GSM605308 \tab R2 \cr
#' ESH1 \tab H3K4me3 \tab GSM605315 \tab R2 \cr
#' ESH1 \tab input \tab GSM605339 \tab R3 \cr
#' ESH1 \tab H3K27me3 \tab GSM466734 \tab R3 \cr
#' ESH1 \tab H3K4me3 \tab GSM469971 \tab R3 \cr
#' HUES64 \tab input \tab GSM772754 \tab R1 \cr
#' HUES64 \tab H3K27me3 \tab GSM772750 \tab R1 \cr
#' HUES64 \tab H3K4me3 \tab GSM772752 \tab R1 \cr
#' HUES64 \tab input \tab GSM772807 \tab R2 \cr
#' HUES64 \tab H3K27me3 \tab GSM669974 \tab R2 \cr
#' HUES64 \tab H3K4me3 \tab GSM669967 \tab R2 \cr
#' HUES48 \tab input \tab GSM772755 \tab R1 \cr
#' HUES48 \tab H3K27me3 \tab GSM669942 \tab R1 \cr
#' HUES48 \tab H3K4me3 \tab GSM669936 \tab R1 \cr
#' HUES48 \tab input \tab GSM772794 \tab R2 \cr
#' HUES48 \tab H3K27me3 \tab GSM772766 \tab R2 \cr
#' HUES48 \tab H3K4me3 \tab GSM772797 \tab R2 \cr
#' }
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312 ENCODE Data}
#' @source \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5354816/}{An annotated list of bivalent chromatin regions in human ES cells: a new tool for cancer epigenetic research}
#' @docType data
#' @keywords datasets
#' @name hESC_data
NULL



#' GRanges region annotation of hESC bivalent sites.
#'
#' \code{\link{hESC_data}}
#' @format A GRanges with 30 rows and mcols:
#' \describe{
#'   \item{group}{a bivalency category from \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5354816/}{PMC5354816}}
#'   \item{id}{an arbitrary unique identifier}
#' }
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312}
#' @docType data
#' @keywords datasets
#' @name query_gr
NULL

#' data.table of bam pileups
#'
#' \code{\link{hESC_data}}
#' @format data.table from fetch
#' \describe{
#'   \item{id}{unique region identifier}
#'   \item{cell}{cell facet}
#'   \item{mark}{mark facet}
#'   \item{x}{relative position to center of region}
#'   \item{y}{mean read pileup}
#' }
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312}
#' @docType data
#' @keywords datasets
#' @name bam_pileup_dt
NULL

#‘bam_pileup_dt’ ‘profile_datatable’ ‘profile_dt’ ‘tsne_dt’

#' data.table of fold-enrichment
#'
#' \code{\link{hESC_data}}
#' @format data.table from fetch
#' \describe{
#'   \item{id}{unique region identifier}
#'   \item{cell}{cell facet}
#'   \item{mark}{mark facet}
#'   \item{x}{relative position to center of region}
#'   \item{y}{mean fold-enrichment}
#' }
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312}
#' @docType data
#' @keywords datasets
#' @name profile_dt
NULL

#' data.table of tsne embedding
#'
#' \code{\link{hESC_data}}
#' @format data.table from stsRunTsne
#' \describe{
#'   \item{tx}{x coordinate}
#'   \item{ty}{y coordinate}
#'   \item{id}{unique region identifier}
#'   \item{cell}{cell facet}
#' }
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17312}{GSE17312}
#' @docType data
#' @keywords datasets
#' @name tsne_dt
NULL

