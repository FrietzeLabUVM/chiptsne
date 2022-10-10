testthat::context("ggimage summaries")

library(data.table)
library(chiptsne)
library(testthat)

print(data("sts.test"))
options("mc.cores" = 2)
sts
npx = 4
npy = 3

profile_dt = sts$signal_data$CTCF_features$CTCF_signal$signal_data
class(sts@signal_data$CTCF_features$CTCF_signal)
tsne_dt = sts@signal_data$CTCF_features$CTCF_signal@xy_data

summary_dt = prep_summary(
    profile_dt,
    tsne_dt,
    x_points = 4,
    y_points = 3,
    xrng = c(-.3, .4),
    yrng = c(-.45, .35),
    wide_var = "wide_var",
    tall_var = "tall_var"
)
debug(prep_images)
img_res = prep_images(
    summary_dt,
    x_points = 4,
    y_points = 3,
    xrng = c(-.3, .4),
    yrng = c(-.45, .35)
)
seqsetvis::safeBrew(n = summary_dt$name)

plot_summary_raster(
    img_res$image_dt,
    x_points = 4,
    y_points = 3,
    xrng = c(-.3, .4),
    yrng = c(-.45, .35),
    line_color_mapping = color_mapping
)

# stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4)
# stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, q_wide_vars = c("H3K4me3"))
# stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, plot_type = "raster")
# stsPlotSummaryProfiles(profile_dt, position_dt = tsne_dt, x_points = 4, q_wide_vars = c("H3K4me3"), plot_type = "raster")

test_that("prep_images names of outputs", {
    expect_equal(
        names(img_res),
        c(
            "image_dt",
            "summary_profile_dt",
            "x_points",
            "y_points",
            "xrng",
            "yrng",
            "line_color_mapping"
        )
    )
})

test_that("prep_images variables of image_dt", {
    expect_equal(colnames(img_res$image_dt),
                 c("bx", "by", "plot_id", "png_file", "tx", "ty", "N"))
})

test_that("prep_images variables of summary_profile_dt", {
    expect_equal(
        colnames(img_res$summary_profile_dt),
        c(
            "bx",
            "by",
            "x",
            "wide_var",
            "y",
            "N",
            "plot_id",
            "ynorm",
            "group"
        )
    )
})

test_that("prep_images parameter passthrough", {
    expect_equal(img_res$x_points, 4)
    expect_equal(img_res$y_points, 3)
    expect_equal(img_res$xrng, c(-.3, .4))
    expect_equal(img_res$yrng, c(-.45, .35))
})

img_rect = set_image_rects(
    img_res$image_dt,
    x_points = img_res$x_points,
    y_points = img_res$y_points,
    xrng = img_res$xrng,
    yrng = img_res$yrng
)

test_that("set_image_rects variables", {
    expect_s3_class(img_rect, "data.frame")
    expect_s3_class(img_rect, "data.table")
    expect_equal(
        colnames(img_rect),
        c(
            "bx",
            "by",
            "plot_id",
            "png_file",
            "tx",
            "ty",
            "N",
            "img_size",
            "xmin",
            "xmax",
            "ymin",
            "ymax"
        )
    )
})



test_that("set_image_rects variables", {
    p = plot_summary_raster(
        img_res$image_dt,
        x_points = img_res$x_points,
        y_points = img_res$y_points,
        xrng = img_res$xrng,
        yrng = img_res$yrng
    )
    expect_s3_class(p, "gg")
})

test_that("return_data returns data.table, glyph no facet", {
    r1 = stsPlotSummaryProfiles(
        profile_dt,
        position_dt = tsne_dt,
        x_points = 4,
        return_data = TRUE,
        plot_type = "glyph"
    )
    expect_s3_class(r1, "data.table")
})

test_that("return_data returns data.table, raster facet", {
    r2 = stsPlotSummaryProfiles(
        profile_dt,
        position_dt = tsne_dt,
        x_points = 4,
        return_data = TRUE,
        plot_type = "raster"
    )
    expect_s3_class(r2, "data.table")
})
test_that("return_data returns data.table, glyph with facet", {
    r1 = stsPlotSummaryProfiles(
        profile_dt,
        position_dt = tsne_dt,
        x_points = 4,
        return_data = TRUE,
        plot_type = "glyph",
        facet_byCell = TRUE
    )
    expect_s3_class(r1, "data.table")
})
test_that("return_data returns data.table, raster with facet", {
    r2 = stsPlotSummaryProfiles(
        profile_dt,
        position_dt = tsne_dt,
        x_points = 4,
        return_data = TRUE,
        plot_type = "raster",
        facet_byCell = TRUE
    )
    expect_s3_class(r2, "data.table")
})
# plot_summary_glyph(mdt, xrng = xrng, yrng = yrng,
#                    x_points = x_points, y_points = y_points,
#                    min_size = 0, N_ceiling = 4)
#
# stsPlotSummaryProfiles(profile_dt, tsne_dt, 6, plot_type = "glyph", facet_byCell = TRUE)
# stsPlotSummaryProfiles(profile_dt, tsne_dt, 6, plot_type = "glyph", facet_byCell = FALSE)
