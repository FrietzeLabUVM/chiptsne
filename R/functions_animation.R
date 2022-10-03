
animate_summary_profiles = function(profile_dt,
                                    position_dt,
                                    x_points,
                                    y_points = x_points,
                                    order_var = NULL,
                                    N_ceiling = 20,
                                    min_size = 1,
                                    ## view domain
                                    xrng = range(position_dt$tx),
                                    yrng = range(position_dt$ty),
                                    ## animation
                                    transition_length = 4, state_length = 1,
                                    nframes = 150, duration = 10,
                                    ##
                                    instead_facet = FALSE

                                    ){
    p_dt = copy(profile_dt)
    if(!is.null(order_var)){
        if(!order_var %in% colnames(profile_dt)){
            stop(order_var, " not in colnames of profile_dt")
        }else{
            p_dt = p_dt[order(get(order_var))]
        }
    }
    p_dt$wide_var = factor(p_dt$wide_var, levels = unique(p_dt$wide_var))
    summ_dt = stsPlotSummaryProfiles(
        p_dt,
        position_dt,
        x_points = x_points,
        y_points = y_points,
        xrng = xrng,
        yrng = yrng,
        return_data = TRUE,
        N_ceiling = N_ceiling,
        min_size = min_size
    )
    dc_dt = dcast(summ_dt, plot_id+x~wide_var, value.var = "gy")
    dt = reshape2::melt(dc_dt, id.vars = c("plot_id", "x"))
    dt = dt[, var(value), by = .(plot_id, x)][, .(dynamic = max(V1)), .(plot_id)]

    summ_dt = merge(summ_dt, dt, by = "plot_id")

    lab_dt = unique(summ_dt[, .(wide_var)])
    lab_dt[, gx := scales::rescale(seq_len(nrow(lab_dt))-1, c(-.5, -.2))]
    lab_dt$gy = .45
    lab_dt[, hour := sub("_.+", "", wide_var)]

    summ_p = ggplot() +
        geom_path(data = summ_dt,
                  aes(
                      x = gx,
                      y = gy,
                      group = plot_id,
                      color = dynamic
                  )) +
        coord_fixed() +
        geom_label(data = lab_dt, aes(x = gx,
                                      y = gy,
                                      label = hour)) +
        theme(
            panel.background = element_rect(fill = "darkgray"),
            panel.grid.minor = element_blank()
        ) +
        # facet_wrap("wide_var") +
        scale_color_viridis_c()
    # summ_p
    if(instead_facet){
        p = summ_p + facet_wrap("wide_var")
        ggsave(paste0("across_timecourse_summary_prof_v2_", m, ".png"),
               plot = p, width = 8, height = 8)
    }else{
        anim = summ_p +
            # geom_label(data = lab_dt, aes(label = hour, group = hour), fill = "white") +
            transition_states(wide_var, transition_length = transition_length, state_length = state_length, wrap = FALSE) +
            enter_fade() +
            exit_fade()
        render = animate(anim, nframes = nframes, duration = duration)
        save_animation(render, paste0("across_timecourse_summary_prof_v2_", m, ".gif"))
    }

}




#calculate variance per subplot

#
#
#
# agg_dt = tsne_input$bw_dt[, .(y = mean(y)), .(id, cell, mark, hour)]
# agg_dt = merge(agg_dt, position_dt)
# agg_dt[, ynorm := y / quantile(y, .995), by = .(mark)]
# agg_dt[ynorm > 1, ynorm := 1]
# ggplot(agg_dt, aes(x = tx, y = ty, color = y)) +
#     facet_wrap(hour~.) + geom_point(size = .2) + scale_color_viridis_c() +
#     labs(title = m)
# ggsave(paste0("across_timecourse_points_v2_", m, ".png"), width = 8, height = 8)
#
# bxval = "tx"
# byval = "ty"
# val = "y"
# facet_ = c("cell", "mark", "hour")
# nbins = 80
# bin_met = mean
# agg_dt[, bx := bin_values(tx, n_bins = nbins, xrng = c(-.5, .5))]
# agg_dt[, by := bin_values(ty, n_bins = nbins, xrng = c(-.5, .5))]
#
# bin_dt = agg_dt[, .(y = bin_met(y)), c(facet_, "bx", "by")]
# bvc = bin_values_centers(n_bins = nbins, c(-.5, .5))
# bin_dt[, tx := bvc[bx]]
# bin_dt[, ty := bvc[by]]
#
#
#
#
# if(!is.null(bin_dt$facet)){
#     bin_dt$facet = NULL
# }
# bin_dt = bin_dt[order(cell)][order(mark)][order(hour)]
# for(fac in facet_){
#     if(is.null(bin_dt$facet)){
#         bin_dt[, facet := get(fac)]
#     }else{
#         bin_dt[, facet := paste(facet, get(fac), sep = "_")]
#     }
# }
#
#
#
# lab_dt = unique(bin_dt[, .(hour, facet)])
# lab_dt$ty = .5
# lab_dt$tx = scales::rescale(1:8, to = c(-.5, -.1))
# lab_dt$facet = factor(lab_dt$facet, levels = lab_dt$facet)
#
# bin_dt$facet = factor(bin_dt$facet, levels = levels(lab_dt$facet))
#
# p_raster = ggplot(bin_dt, aes(x = tx, y = ty, fill = y)) +
#     geom_raster() + #facet_wrap(facet_) +
#     scale_fill_viridis_c()
# ggsave(paste0("across_timecourse_raster_v2_", m, ".png"),
#        plot = p_raster + facet_wrap(facet_),
#        width = 8, height = 8)
#
#
# p_raster = ggplot(bin_dt, aes(x = tx, y = ty, fill = y)) +
#     geom_raster() + #facet_wrap(facet_) +
#     scale_fill_viridis_c()
#
# dbl_lev = levels(bin_dt$facet)
# dbl_lev = dbl_lev[seq(1, length(dbl_lev)-1)]
#
# tmp_dt = rbind(bin_dt,
#                bin_dt[facet %in% dbl_lev, .(cell, mark, hour, bx, by, y , tx, ty, facet = paste0(facet, "rev"))][order(hour, decreasing = TRUE)])
#
# lab_dt = unique(tmp_dt[, .(hour, facet)])
# lab_dt$ty = .5
# xs = scales::rescale(1:8, to = c(-.5, -.1))
# names(xs) = sort(unique(lab_dt$hour))
# lab_dt$tx = xs[as.character(lab_dt$hour)]
# lab_dt$facet = factor(lab_dt$facet, levels = lab_dt$facet)
#
#
#
#
# p_raster2 = ggplot(tmp_dt, aes(x = tx, y = ty, fill = y)) +
#     geom_raster() + #facet_wrap(facet_) +
#     scale_fill_viridis_c()
#
#
# anim = p_raster2 +
#     geom_label(data = lab_dt, aes(label = hour, group = hour), fill = "white") +
#     transition_states(facet, transition_length = 4, state_length = 1, wrap = FALSE) +
#     enter_fade() +
#     exit_fade()
# render = animate(anim, nframes = 150, duration = 10)
# anim_save(paste0("across_timecourse_raster_v2_", m, ".gif") , render)
