testthat::context("ggimage summaries")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")
options("mc.cores" = 2)

npx = 4
npy = 3
summary_dt = prep_summary(profile_dt, tsne_dt,
                          x_points = 4, y_points = 3,
                          xrng = c(-.3, .4), yrng = c(-.45, .35))
img_res = prep_images(summary_dt,
                      x_points = 4, y_points = 3,
                      xrng = c(-.3, .4), yrng = c(-.45, .35))
color_mapping = seqsetvis::safeBrew(2)
names(color_mapping) = unique(summary_dt$wide_var)
plot_summary_raster(img_res$image_dt, x_points = 4, y_points = 3,
                    xrng = c(-.3, .4), yrng = c(-.45, .35), line_color_mapping = color_mapping)

stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4)
stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, q_wide_vars = c("H3K4me3"))
stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, plot_type = "raster")
stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, q_wide_vars = c("H3K4me3"), plot_type = "raster")

test_that("prep_images names of outputs", {
    expect_equal(names(img_res), c("image_dt", "summary_profile_dt", "x_points", "y_points", "xrng", "yrng", "line_color_mapping"))
})

test_that("prep_images variables of image_dt", {
    expect_equal(colnames(img_res$image_dt), c("bx", "by", "plot_id", "png_file", "tx", "ty", "N"))
})

test_that("prep_images variables of summary_profile_dt", {
    expect_equal(colnames(img_res$summary_profile_dt), c("bx", "by", "x", "wide_var", "y", "N", "plot_id", "ynorm", "group"))
})

test_that("prep_images parameter passthrough", {
    expect_equal(img_res$x_points, 4)
    expect_equal(img_res$y_points, 3)
    expect_equal(img_res$xrng, c(-.3, .4))
    expect_equal(img_res$yrng, c(-.45, .35))
})

img_rect = set_image_rects(img_res$image_dt,
                           x_points = img_res$x_points, y_points = img_res$y_points,
                           xrng = img_res$xrng, yrng = img_res$yrng)

test_that("set_image_rects variables", {
    expect_s3_class(img_rect, "data.frame")
    expect_s3_class(img_rect, "data.table")
    expect_equal(colnames(img_rect),
                 c("bx", "by", "plot_id", "png_file", "tx", "ty", "N", "img_size", "xmin", "xmax", "ymin", "ymax"))
})



test_that("set_image_rects variables", {
    p = plot_summary_raster(img_res$image_dt,
                      x_points = img_res$x_points, y_points = img_res$y_points,
                      xrng = img_res$xrng, yrng = img_res$yrng)
    expect_s3_class(p, "gg")
})


xrng = c(-.5, 0)
yrng = c(-.5, 0)

xrng = c(-.5, .5)
yrng = c(-.5, .5)

x_points = y_points = 6

mdt = prep_summary(profile_dt, tsne_dt, xrng = xrng, yrng = yrng,
                   x_points = x_points, y_points = y_points)
plot_summary_glyph(mdt, xrng = xrng, yrng = yrng,
                   x_points = x_points, y_points = y_points,
                   min_size = 0, N_ceiling = 4)

stsPlotSummaryProfiles(profile_dt, tsne_dt, 6, plot_type = "glyph", facet_byCell = TRUE)
stsPlotSummaryProfiles(profile_dt, tsne_dt, 6, plot_type = "glyph", facet_byCell = FALSE)
