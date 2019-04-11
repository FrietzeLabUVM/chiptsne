testthat::context("tsne velocity")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")
options("mc.cores" = 4)

vel_dt = prep_velocity(tsne_dt, "HUES48", "HUES64")

test_that("prep_velocity variables", {
    expect_equal(
        colnames(vel_dt),
        c(
            "id",
            "tx_cell_a",
            "tx_cell_b",
            "ty_cell_a",
            "ty_cell_b",
            "bx_cell_a",
            "btx_cell_a",
            "by_cell_a",
            "bty_cell_a",
            "angle",
            "foreground",
            "distance"
        )
    )
    expect_equal(nrow(vel_dt), length(query_gr))
})

test_that("plot_velocity_bins ggplot out", {
    p1 = plot_velocity_bins(vel_dt)
    expect_s3_class(p1, "ggplot")
    p2 = plot_velocity_bins(vel_dt, bin_FUN = mean)
    expect_s3_class(p2, "ggplot")
    p3 = plot_velocity_bins(vel_dt, bin_FUN = sum)
    expect_s3_class(p3, "ggplot")
})

test_that("plot_velocity_centered ggplot out", {
    p1 = plot_velocity_centered(vel_dt)
    expect_s3_class(p1, "ggplot")
})

test_that("plot_velocity_arrows ggplot out", {
    p1 = plot_velocity_arrows(vel_dt)
    expect_s3_class(p1, "ggplot")
})

test_that("plot_regional_velocity ggplot out", {
    p1 = plot_regional_velocity(tsne_dt, "HUES48", "HUES64", 4)
    expect_s3_class(p1, "ggplot")
})

test_that("plot_recentered_velocity ggplot out", {
    p1 = plot_recentered_velocity(tsne_dt, "HUES48", "HUES64", 4, angle_as_color = TRUE)
    p1
    expect_s3_class(p1, "ggplot")
})


