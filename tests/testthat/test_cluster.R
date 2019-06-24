testthat::context("tsne cluster")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")

test_that("nn clust", {
    clust_res = nn_clust(tsne_dt, nn = 10)
    clust_dt = clust_res$data
    clust_p = clust_res$plot
    expect_s3_class(clust_p, "ggplot")
    expect_s3_class(clust_dt, "data.table")
    expect_equal(colnames(clust_dt), c("tid", "tx", "ty", "id", "tall_var", "cluster_id"))
})

test_that("nn combine", {
    clust_res = nn_clust(tsne_dt, nn = 5)
    clust_dt = clust_res$data
    clust_p = clust_res$plot

    clust_dt6 = combine_mostsimilar(p_dt = clust_dt, profile_dt = profile_dt, n_times = 2,
                                    new_cluster_ = "meta")
    expect_s3_class(clust_p, "ggplot")
    expect_s3_class(clust_dt, "data.table")
    nclust = length(unique(clust_dt$cluster_id))
    expect_equal(length(unique(clust_dt6$cluster_id)), nclust)
    expect_equal(length(unique(clust_dt6$meta)), nclust - 2)
})

test_that("clust summary combine", {
    clust_res = nn_clust(tsne_dt, nn = 5)
    clust_dt = clust_res$data
    p_summary = stsPlotClusterProfiles(profile_dt = profile_dt, cluster_dt = clust_dt, cluster_ = "cluster_id")

    expect_s3_class(p_summary, "ggplot")
})

