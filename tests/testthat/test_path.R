testthat::context("tsne path")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")

test_that("plot_path", {
    p1 = plot_path(tsne_dt, unique(tsne_dt$tall_var), query_gr$id[1:3])
    expect_s3_class(p1, "ggplot")
    p2 = plot_path(tsne_dt, rev(unique(tsne_dt$tall_var)), query_gr$id[4:7], arrow_FUN = arrow(length = unit(.04, "npc")))
    expect_s3_class(p2, "ggplot")
})

test_that("plot_outline", {
    sel_id = query_gr$id[3:5]
    p1 = plot_outline(tsne_dt,
                      unique(tsne_dt$tall_var),
                      id_to_plot = setdiff(query_gr$id, sel_id),
                      bg_points = -1)
    expect_s3_class(p1, "ggplot")
    p2 = plot_outline(tsne_dt, unique(tsne_dt$tall_var),
                      p = p1,
                      id_to_plot = sel_id,
                      bg_points = 0,
                      bg_color = "red",
                      line_color_mapping = "red", label_type = "label",
                      fill_color_mapping = "#FF000055")
    expect_s3_class(p2, "ggplot")
})
