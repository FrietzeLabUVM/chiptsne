#' Title
#'
#' @param profiles_dt
#' @param position_dt
#' @param n_points
#' @param xrng
#' @param yrng
#' @param rname
#' @param odir
#' @param force_rewrite
#' @param apply_norm
#' @param ylim
#' @param facet_by
#' @param ma_size
#' @param n_cores
#' @param line_colors
#'
#' @return
#' @export
#' @importFrom seqsetvis applySpline
#'
#' @examples
make_tsne_img = function(profiles_dt, position_dt, n_points,
                         xrng = range(position_dt$tx),
                         yrng = range(position_dt$ty),
                         rname = digest::digest(list(
                             profiles_dt, position_dt,
                             n_points, apply_norm,
                             ylim, line_colors, facet_by)),
                         odir = file.path("tsne_images/", rname),
                         force_rewrite = FALSE,
                         apply_norm = TRUE,
                         ylim = c(0, 1),
                         facet_by = NULL,
                         # view_rect = list(),
                         ma_size = 2,
                         n_cores = getOption("mc.cores", 1),
                         line_colors = c("H3K4me3" = "forestgreen",
                                         "H3K27me3" = "firebrick1"),
                         group_mapping = NULL){
    if(!all(unique(profiles_dt$mark) %in% names(line_colors))){
        missing_colors = setdiff(unique(profiles_dt$mark), names(line_colors))
        stop("line_colors is missing assignments for: ", paste(missing_colors, collapse = ", "))
    }
    stopifnot()
    # stopifnot(is.list(view_rect))
    # if(is.null(view_rect$xmin))
    position_dt = copy(position_dt)
    position_dt = position_dt[tx >= min(xrng) & tx <= max(xrng) & ty >= min(yrng) & ty <= max(yrng)]
    #use positional info from position_dt to bin points
    position_dt[, bx := mybin(tx, n_points, xrng = xrng)]
    position_dt[, by := mybin(ty, n_points, xrng = yrng)]
    #merge binning info to profiles
    mdt = merge(profiles_dt, position_dt[, list(bx, by, cell, id)], allow.cartesian=TRUE, by = intersect(colnames(profiles_dt), c("cell", "id")))#, by = c("cell", "id"))
    if(is.null(mdt$mark)) mdt$mark = "signal"

    if(is.null(facet_by)){
        mdt = mdt[, list(y = mean(y)), list(bx, by, x, mark)]
    }else{
        mdt = mdt[, list(y = mean(y)), list(bx, by, x, mark, get(facet_by))]
        colnames(mdt)[colnames(mdt) == "get"] = facet_by
    }
    #each combination of bx and by is a unique plot_id
    mdt[, plot_id := paste(bx, by, sep = "_")]
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    if(is.null(facet_by)){
        img_dt = unique(mdt[, list(bx, by, plot_id)])
        img_dt[, png_file := file.path(odir, paste0(plot_id, ".png"))]
    }else{
        img_dt = unique(mdt[, list(bx, by, plot_id, get(facet_by))])
        colnames(img_dt)[4] = facet_by
        img_dt[, png_file := file.path(odir, paste0(get(facet_by), "_", plot_id, ".png"))]
    }

    xs = mybin_centers(tsne_res$tx, n_points, xrng = xrng)
    ys = mybin_centers(tsne_res$ty, n_points, xrng = yrng)

    # plot(expand.grid(xs, ys), xlim = xrng, ylim = yrng)
    # rect(min(xrng), min(yrng), max(xrng), max(yrng), col = rgb(0,0,1,.1))

    img_dt[, tx := xs[bx]]
    img_dt[, ty := ys[by]]


    if(apply_norm){
        mdt[, ynorm := y / quantile(y, .95), by = list(mark)]
        mdt[ynorm > 1, ynorm := 1]
    }else{
        mdt[, ynorm := y]
    }


    if(force_rewrite){
        file.remove(img_dt$png_file[file.exists(img_dt$png_file)])
    }

    if(is.null(facet_by)){
        img_dt = merge(img_dt, position_dt[, .N, list(bx, by)])
    }else{
        tmp_dt = position_dt[, .N, list(bx, by, get(facet_by))]
        colnames(tmp_dt)[colnames(tmp_dt) == "get"] = facet_by
        img_dt = merge(img_dt, tmp_dt)
    }

    if(is.null(group_mapping)){
        mdt$group = 1
    }else{
        if(!all(unique(mdt$mark) %in% names(group_mapping))){
            missing_colors = setdiff(unique(mdt$mark), names(group_mapping))
            stop("group_mapping is missing assignments for: ", paste(missing_colors, collapse = ", "))
        }
        mdt$group = group_mapping[mdt$mark]
    }


    if(any(!file.exists(img_dt$png_file))){
        plot_info = lapply(which(!file.exists(img_dt$png_file)), function(i){
            # hidden = parallel::mclapply(which(!file.exists(img_dt$png_file)), function(i){}, mc.cores = n_cores)
            fpath = img_dt$png_file[i]
            p_id = img_dt$plot_id[i]
            if(is.null(facet_by)){
                pdt = mdt[plot_id == p_id]
            }else{
                pdt = mdt[plot_id == p_id & get(facet_by) == img_dt[[facet_by]][i]]
            }

            # pdt[, ysm := seqsetvis:::movingAverage(ynorm, n = ma_size), by = list(mark)]
            pdt[, ysm := movingAverage(ynorm, n = ma_size), by = list(mark)]

            pdt = seqsetvis::applySpline(pdt, n = 10, by_ = "mark", y_ = "ysm")
            list(pdt, fpath)
        })
        # hidden = lapply(plot_info, function(x){
        # figure out how not to copy global env
        hidden = parallel::mclapply(plot_info, function(x){
            fpath = x[[2]]
            pdt = x[[1]]
            # pdt[, ysm := seqsetvis:::movingAverage(y, n = 8), by = list(mark)]

            p = ggplot(pdt, aes(x = x, y = ysm, ymin = 0, ymax = ysm, color = mark, fill = mark)) +
                geom_ribbon(alpha = .3) +
                geom_path(size = .6, alpha = 1) +
                scale_color_manual(values = line_colors) +
                scale_fill_manual(values = line_colors) +
                theme_void() + guides(color = "none", fill =
                                          'none') +
                facet_wrap("group", ncol = 1) +
                theme(strip.background = element_blank(), strip.text = element_blank()) +
                coord_cartesian(ylim = ylim, xlim = c(-.5, .5), expand = FALSE)
            ggsave(fpath, p, width = 2, height = 2, units = "cm")
            # p
            NULL
            # })
        }, mc.cores = n_cores)
    }

    img_dt[, png_file := normalizePath(png_file)]

    if(is.factor(position_dt$cell)){
        if(!is.null(img_dt$cell)){
            img_dt$cell = factor(img_dt$cell, levels = levels(position_dt$cell))
        }
        if(!is.null(mdt$cell)){
            mdt$cell = factor(mdt$cell, levels = levels(position_dt$cell))
        }
    }

    return(list(images_dt = img_dt, summary_profiles_dt = mdt, tsne_dt = position_dt, n_points = n_points, xrng = xrng, yrng = yrng))
}

#' Title
#'
#' @param simg_dt
#' @param n_points
#' @param xrng
#' @param yrng
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#'
#' @return
#' @export
#'
#' @examples
prep_tsne_img = function(simg_dt,
                         n_points,
                         xrng,
                         yrng,
                         N_floor = 0,
                         N_ceiling = NULL,
                         min_size = .3
){
    if(is.null(N_ceiling)){
        N_ceiling = max(simg_dt$N)
    }
    simg_dt[, img_size := N]
    simg_dt[img_size > N_ceiling, img_size := N_ceiling]
    simg_dt[img_size < N_floor, img_size := N_floor]

    simg_dt[, img_size := img_size - N_floor]
    simg_dt[, img_size := img_size / N_ceiling]
    simg_dt[, img_size := img_size]
    simg_dt = simg_dt[img_size >= min_size]

    xspc = diff(xrng)/n_points/2
    yspc = diff(yrng)/n_points/2

    simg_dt[, xmin := tx - xspc * img_size]
    simg_dt[, xmax := tx + xspc * img_size]
    simg_dt[, ymin := ty - yspc * img_size]
    simg_dt[, ymax := ty + yspc * img_size]
    simg_dt
}

#' Title
#'
#' @param images_dt
#' @param n_points
#' @param p
#' @param xrng
#' @param yrng
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param show_plot
#'
#' @return
#' @export
#'
#' @examples
plot_tsne_img = function(images_dt,
                         n_points,
                         p = NULL,
                         xrng = c(-.5, .5),
                         yrng = c(-.5, .5),
                         N_floor = 0,
                         N_ceiling = NULL,
                         min_size = .3,
                         show_plot = FALSE
){
    simg_dt = copy(images_dt)
    simg_dt = prep_tsne_img(simg_dt,
                            n_points = n_points,
                            xrng = xrng,
                            yrng = yrng,
                            N_floor = N_floor,
                            N_ceiling = N_ceiling,
                            min_size = min_size
    )
    if(is.null(p)) p = ggplot()
    p = p +
        geom_image.rect(data = simg_dt,
                        aes(xmin = xmin, xmax = xmax,
                            ymin = ymin, ymax = ymax,
                            image = png_file)) +#, color = rgb(0,0,1,.2)) +
        geom_rect(data = simg_dt,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, color = "black") +
        coord_cartesian(xlim = xrng, ylim = yrng)
    if(show_plot) print(p)
    invisible(list(plot = p, plot_data = simg_dt))
}

#' Title
#'
#' @param images_dt
#' @param tsne_dt
#' @param n_points
#' @param p
#' @param xrng
#' @param yrng
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param show_plot
#'
#' @return
#' @export
#'
#' @examples
plot_tsne_img_byCell = function(images_dt,
                                tsne_dt,
                                n_points,
                                p = NULL,
                                xrng = c(-.5, .5),
                                yrng = c(-.5, .5),
                                N_floor = 0,
                                N_ceiling = NULL,
                                min_size = .3,
                                show_plot = FALSE
){
    cell = bx = by = NULL;
    # tsne_dt$bx = mybin(tsne_dt$tx, n_points, xrng)
    # tsne_dt$by = mybin(tsne_dt$ty, n_points, xrng)
    simg_dt = merge(images_dt[, list(bx, by, plot_id, png_file, tx, ty)],
                    tsne_dt[, list(N = .N), list(cell, bx, by)], by = c("bx", "by"), allow.cartesian = TRUE)
    simg_dt = prep_tsne_img(simg_dt,
                            n_points = n_points,
                            xrng = xrng,
                            yrng = yrng,
                            N_floor = N_floor,
                            N_ceiling = N_ceiling,
                            min_size = min_size
    )
    if(is.null(p)) p = ggplot()
    p = p +
        geom_image.rect(data = simg_dt,
                        aes(xmin = xmin, xmax = xmax,
                            ymin = ymin, ymax = ymax,
                            image = png_file)) +#, color = rgb(0,0,1,.2)) +
        geom_rect(data = simg_dt,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, color = "black") +
        coord_cartesian(xlim = xrng, ylim = yrng) +
        facet_wrap("cell", drop = FALSE)
    if(show_plot) print(p)
    invisible(list(plot = p, plot_data = simg_dt))
}

#' Title
#'
#' @param img_results
#' @param qcell
#' @param xrng
#' @param yrng
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param facet_by
#'
#' @return
#' @export
#'
#' @examples
make_img_plots_facet = function(img_results, qcell = NULL,
                                xrng = c(-.5, .5),
                                yrng = c(-.5, .5),
                                N_floor = 0,
                                N_ceiling = NULL,
                                min_size = .3,
                                facet_by = "cell"){
    return_list = TRUE
    if(all(c("images_dt", "summary_profiles_dt", "tsne_dt") %in% names(img_results))){
        img_results = list(img_results)
        return_list = FALSE
    }
    stopifnot(is.list(img_results))
    # stopifnot(is.list(img_results$img_res))

    if(is.null(qcell))
        qcell =
        as.character(unique(
            img_results[[1]]$tsne_dt$cell
        ))

    plots = lapply(img_results, function(x){
        img_dt = copy(x$images_dt)
        # img_dt$N = NULL
        # tdt = x$tsne_dt[cell %in% qcell, list(.N), list(bx, by)]
        # img_dt = merge(img_dt, tdt, by = c("bx", "by"))
        plot_tsne_img(img_dt,
                      n_points = x$n_points,
                      N_floor = N_floor,
                      N_ceiling = N_ceiling,
                      min_size = min_size,
                      show_plot = FALSE,
                      xrng = xrng,
                      yrng = yrng)$plot +
            coord_cartesian(xlim = xrng, ylim = yrng) +
            facet_wrap(facet_by, drop = FALSE)
    })


    # pg = cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    # pg
    if(return_list){
        plots
    }else{
        plots[[1]]
    }

}

#' Title
#'
#' @param img_results
#' @param qcell
#' @param xrng
#' @param yrng
#' @param N_floor
#' @param N_ceiling
#' @param min_size
#' @param as_facet
#'
#' @return
#' @export
#'
#' @examples
make_img_plots = function(img_results, qcell = NULL,
                          xrng = c(-.5, .5),
                          yrng = c(-.5, .5),
                          N_floor = 0,
                          N_ceiling = NULL,
                          min_size = .3,
                          as_facet = TRUE){
    return_list = TRUE
    if(all(c("images_dt", "summary_profiles_dt", "tsne_dt") %in% names(img_results))){
        img_results = list(img_results)
        return_list = FALSE
    }
    stopifnot(is.list(img_results))

    # img_results = lapply(img_results, function(x){
    #     if(is.null(x$cell)){
    #         x$cell = factor("cell")
    #     }
    #     x
    # })
    if(is.null(qcell))
        qcell =
        levels(img_results[[1]]$tsne_dt$cell)


    if(as_facet){
        plots = lapply(img_results, function(x){
            pdt = x$tsne_dt[cell %in% qcell]
            pdt$cell = factor(pdt$cell, levels = qcell)
            plot_tsne_img_byCell(x$images_dt,
                                 pdt,
                                 n_points = x$n_points,
                                 N_floor = N_floor,
                                 N_ceiling = N_ceiling,
                                 min_size = min_size,
                                 show_plot = FALSE,
                                 xrng = xrng, yrng =
                                     yrng)$plot +
                coord_cartesian(xlim = xrng, ylim = yrng)
        })
    }else{
        plots = lapply(img_results, function(x){
            img_dt = copy(x$images_dt)
            img_dt$N = NULL
            tdt = x$tsne_dt[cell %in% qcell, list(.N), list(bx, by)]
            img_dt = merge(img_dt, tdt)
            plot_tsne_img(img_dt,
                          n_points = x$n_points,
                          N_floor = N_floor,
                          N_ceiling = N_ceiling,
                          min_size = min_size,
                          show_plot = FALSE,
                          xrng = xrng,
                          yrng = yrng)$plot +
                coord_cartesian(xlim = xrng, ylim = yrng)
        })
    }


    # pg = cowplot::plot_grid(plotlist = plots, nrow = length(plots))
    # pg
    if(return_list){
        plots
    }else{
        plots[[1]]
    }
}
