testthat::context("ggimage summaries")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")
options("mc.cores" = 4)

npx = 4
npy = 3
img_res = prep_images(profile_dt, tsne_dt,
                        x_points = 4, y_points = 3,
                        xrng = c(-.3, .4), yrng = c(-.45, .35))

test_that("stsPrepImages names of outputs", {
    expect_equal(names(img_res), c("image_dt", "summary_profile_dt", "tsne_dt", "x_points", "y_points", "xrng", "yrng"))
})

test_that("stsPrepImages variables of image_dt", {
    expect_equal(colnames(img_res$image_dt), c("bx", "by", "plot_id", "png_file", "tx", "ty", "N"))
})

test_that("stsPrepImages variables of summary_profile_dt", {
    expect_equal(colnames(img_res$summary_profile_dt), c("bx", "by", "x", "mark", "y", "plot_id", "ynorm", "group"))
})

test_that("stsPrepImages variables of tsne_dt", {
    expect_equal(colnames(img_res$tsne_dt), c("tx", "ty", "id", "cell", "bx", "by"))
})

test_that("stsPrepImages parameter passthrough", {
    expect_equal(img_res$x_points, 4)
    expect_equal(img_res$y_points, 3)
    expect_equal(img_res$xrng, c(-.3, .4))
    expect_equal(img_res$yrng, c(-.45, .35))
})

img_rect = set_image_rects(img_res$image_dt,
                x_points = img_res$x_points, y_points = img_res$y_points,
                xrng = img_res$xrng, yrng = img_res$yrng)

p = plot_tsne_img(img_res$image_dt,
              x_points = img_res$x_points, y_points = img_res$y_points,
              xrng = img_res$xrng, yrng = img_res$yrng)
