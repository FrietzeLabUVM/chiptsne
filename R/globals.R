.datatable.aware=TRUE
setOldClass(c('data.frame'))
setOldClass(c('data.table', 'data.frame'))
#data.table package functions
globalVariables(c(":=", ".", ".N", "as.data.table", "copy", "tstrsplit", "dcast"))
#common data.table variables
globalVariables(c(
    "distance",
    "x",
    "y",
    "id",
    "cell",
    "mark",
    "tx",
    "ty",
    "rn",
    "N",
    "angle",
    "angle_bin",
    "bx",
    "bx_cell_a",
    "by_cell_a",
    "cell_o",
    "grp",
    "grp_o",
    "ch_i",
    "img_size",
    "index",
    "norm_factor",
    "o",        "pid",
    "plot_id",
    "png_file",
    "foreground",
    "tx_cell_a",
    "tx_cell_b",
    "tx_end",
    "ty_cell_a",
    "ty_cell_b",
    "ty_end",
    "ysm",
    "flip_strand"
))
#ggplot2 package functions
globalVariables(
    c(
        "alpha",
        "annotate",
        "arrow",
        "coord_cartesian",
        "coord_polar",
        "element_blank",
        "element_text",
        "facet_grid",
        "facet_wrap",
        "geom_curve",
        "geom_path",
        "geom_point",
        "geom_rect",
        "geom_ribbon",
        "geom_segment",
        "ggplot",
        "ggsave",
        "guides",
        "is.grob",
        "labs",
        "parameters",
        "pointsGrob",
        "rectGrob",
        "resize",
        "scale_color_brewer",
        "scale_color_gradientn",
        "scale_color_manual",
        "scale_fill_gradient2",
        "scale_fill_gradientn",
        "scale_fill_manual",
        "scale_size_continuous",
        "scale_x_continuous",
        "theme",
        "theme_classic",
        "theme_void",
        "unit",
        "xmax",
        "xmin",
        "ymax",
        "ymin",
        "ynorm"
    )
)
# .N N alpha angle angle_bin annotate arrow  bin_value bx
# bx_cell_a by_cell_a cell_o ch_i coord_cartesian coord_polar copy
# dcast distance element_blank element_text facet_grid facet_wrap
# flip_strand foreground geom_curve geom_path geom_point geom_rect
# geom_ribbon geom_segment ggplot ggsave gpar grp grp_o guides img_size
# index is.grob labs matrix_file new norm_factor o parameters pid
# plot_id png_file pointsGrob quantile rectGrob regions_file resize
# scale_color_brewer scale_color_gradientn scale_color_manual
# scale_fill_gradient2 scale_fill_gradientn scale_fill_manual
# scale_size_continuous scale_x_continuous spline theme theme_classic
# theme_void tstrsplit tx_cell_a tx_cell_b tx_end ty_cell_a ty_cell_b
# ty_end unit validObject var xmax xmin ymax ymin ynorm ysm
