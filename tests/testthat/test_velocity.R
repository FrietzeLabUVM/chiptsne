testthat::context("tsne velocity")

library(data.table)
library(seqtsne)
library(testthat)

data("query_gr")
data("profile_dt")
data("tsne_dt")
options("mc.cores" = 4)

vel_dt = prep_velocity(tsne_dt, "HUES48", "HUES64", angle.max = 160)#, max_plotted = 5)
plot_velocity_bins(vel_dt[foreground == TRUE])
plot_velocity_bins(vel_dt, bin_FUN = mean)
plot_velocity_bins(vel_dt, bin_FUN = sum)
plot_velocity_centered(vel_dt)
vel_dt$foreground = TRUE
plot_velocity_arrows(vel_dt, angle_as_color = TRUE)
p_arr = plot_velocity_arrows(vel_dt, angle_as_color = FALSE)

p_arr = ggplot() + annotate(
    "segment",
    x = vel_dt$tx_cell_a,
    xend = vel_dt$tx_cell_b,
    y = vel_dt$ty_cell_a,
    yend = vel_dt$ty_cell_b,
    color = "gray",
    arrow = arrow(length = unit(.02, "npc"))
)
p_arr
plot_regional_velocity(
    tsne_dt,
    "HUES48",
    "HUES64",
    n_points = 6,
    p = p_arr,
    strategy = "by_direction"
)
plot_recentered_velocity(
    tsne_dt,
    "HUES48",
    "HUES64",
    n_points = 12,
    p = p_arr
)


plot_regional_velocity(
    tsne_dt,
    "HUES48",
    "HUES64",
    n_points = 3,
    p = p_arr,
    strategy = "by_direction",
    angle_as_color = TRUE
)
plot_regional_velocity(
    tsne_dt,
    "HUES48",
    "HUES64",
    n_points = 6,
    p = p_arr,
    strategy = "by_direction",
    angle_as_color = TRUE
)
plot_recentered_velocity(
    tsne_dt,
    "HUES48",
    "HUES64",
    n_points = 8,
    p = p_arr,
    # arrow_FUN = NULL,
    arrow_FUN = arrow(length = unit(.02, "npc")),
    angle_as_color = TRUE
)

plot_velocity_bins(velocity_dt = vel_dt, bins = 2)
plot_velocity_bins(velocity_dt = vel_dt, bins = 36, bin_FUN = sum)

