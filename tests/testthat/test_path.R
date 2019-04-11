testthat::context("tsne path")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")

plot_path(tsne_dt, query_gr, unique(tsne_dt$cell), query_gr$id[1:3])
plot_path(tsne_dt, query_gr, rev(unique(tsne_dt$cell)), query_gr$id[4:7], arrow_FUN = arrow(length = unit(.04, "npc")))


# tsne_dt,
qgr = query_gr
qcells = unique(tsne_dt$cell)
tss_ids = query_gr$id[1:3]
grp_var = c("id", "grp")[1]
line_type = c("curve", "spline", "straight")[2]
label_type = c("text", "label", "none")[2]
bg_points = 5000
arrow_FUN = NULL

sel_id = query_gr$id[3:5]
p = plot_outline(tsne_dt, query_gr, unique(tsne_dt$cell),
                 id_to_plot = setdiff(query_gr$id, sel_id), bg_points = -1)
plot_outline(tsne_dt, query_gr, unique(tsne_dt$cell), p = p,
             id_to_plot = sel_id, bg_points = 0,
             bg_color = "red",
             line_color_mapping = "red", label_type = "label",
             fill_color_mapping = "#FF000055")
