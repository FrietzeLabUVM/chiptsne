testthat::context("Bigwig files")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")

options("mc.cores" = 2)

bw_files = dir(system.file('extdata', package = "seqtsne"), pattern = ".bw$", full.names = TRUE)
cfg_dt = data.table(file = bw_files)
cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
# save(query_gr, cfg_dt, file = "tmp2_bw.save")
# debug(prep_profile_dt)
tsne_input = stsFetchTsneInput(cfg_dt, query_gr, force_overwrite = TRUE)

test_that("stsFetchTsneInput dimensions of outputs", {
    expect_equal(names(tsne_input), c("bw_dt", "query_gr", "rname"))
    expect_equal(nrow(tsne_input$bw_dt), length(query_gr)*nrow(cfg_dt)*50)
    expect_equal(length(tsne_input$query_gr), length(query_gr))

})
# save(tsne_input, query_gr, cfg_dt, file = "tmp_bw.save")
tsne_res = stsRunTsne(tsne_input$bw_dt, perplexity = 15, force_overwrite = TRUE)

test_that("stsRunTsne dimensions of outputs", {
    expect_equal(colnames(tsne_res), c("tx", "ty", "id", "tall_var"))
    expect_equal(nrow(tsne_res), length(query_gr)*length(unique(cfg_dt$tall_var)))
})

test_that("stsRunTsne config attributes pass through", {
    cfg_dt_extra = cfg_dt
    cfg_dt_extra[, cell := tall_var]
    cfg_dt_extra[, mark := wide_var]
    tsne_input = stsFetchTsneInput(cfg_dt_extra, query_gr, force_overwrite = TRUE)
    obs = sort(colnames(tsne_input$bw_dt))
    obs = setdiff(obs, c("id", "x", "y", "start", "end", "width", "strand", "seqnames"))
    exp = sort(colnames(cfg_dt_extra))
    # exp = setdiff(exp, "file")
    expect_equal(obs, exp)
})
#
#
# #
# # bw_dt = tsne_input$bw_dt
# # setkey(bw_dt, id, tall_var, wide_var)
# #
# # np = 15
# # img_res = make_tsne_img(tsne_input$bw_dt, tsne_res, np, odir = "~/tsne_images", force_rewrite = TRUE)
# # p = plot_tsne_img(img_res$images_dt, np, min_size = .5, N_ceiling = 2)
# # p$plot
# #
# # bw_dt[.(c("region_1", "region_2"), "ESH1")]
# #
# #
# # plot_velocity_arrows(tsne_res, "ESH1", "HUES48", p = p$plot)[[1]]
# #
# # ggplot(tsne_res, aes(x = tx, y = ty)) + geom_point()
# #
# # tsne_res[, id_short := sub(".+_", "", id)]
# # ggplot(tsne_res, aes(x = tx, y = ty, color = tall_var, label = id_short)) + geom_text()
# # lab_dt = tsne_res[, .(tx = mean(tx), ty = mean(ty)), by = .(id_short)]
# # ggplot(tsne_res, aes(x = tx, y = ty, group = id_short, label = id_short)) +
# #     geom_polygon(color = "gray70", fill = "gray90") +
# #     geom_text(data = lab_dt, color = "black", vjust = .5, hjust = .5) +
# #     theme_classic()
# #
# #
#
#
#
