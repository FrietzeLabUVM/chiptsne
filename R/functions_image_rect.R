# these methods are slight alteration of ggimage::geom_image to behave like
# geom_rect

# # grid::grobName
# grobName = function (grob = NULL, prefix = "GRID")
# {
#     if (is.null(grob))
#         grobAutoName(prefix)
#     else {
#         if (!is.grob(grob))
#             stop("invalid 'grob' argument")
#         else grobAutoName(prefix, class(grob)[1L])
#     }
# }
# # grid::grobAutoName
# grobAutoName = function (prefix = "GRID", suffix = "GROB")
# {
#     index <<- index + 1
#     paste(prefix, suffix, index, sep = ".")
# }


# library(magick)
# library(grid)
# ggname = ggimage:::ggname
#' @importFrom grid grobName
ggname = function (prefix, grob)
{
    grob$name <- grid::grobName(grob, prefix)
    grob
}
# color_image = ggimage:::color_image
color_image = function (img, color, alpha = NULL)
{
    if (is.null(color))
        return(img)
    if (length(color) > 1) {
        stop("color should be a vector of length 1")
    }
    bitmap <- img[[1]]
    col <- col2rgb(color)
    bitmap[1, , ] <- as.raw(col[1])
    bitmap[2, , ] <- as.raw(col[2])
    bitmap[3, , ] <- as.raw(col[3])
    if (!is.null(alpha) && alpha != 1)
        bitmap[4, , ] <- as.raw(as.integer(bitmap[4, , ]) * alpha)
    image_read(bitmap)
}

#' @import ggimage
#' @importFrom grid gpar
NULL

draw_key_image = function (data, params, size)
{
    kt <- getOption("ggimage.keytype")
    if (is.null(kt)) {
        kt <- "point"
    }
    if (kt == "point") {
        keyGrob <-
            pointsGrob(
                0.5,
                0.5,
                pch = 19,
                gp = grid::gpar(
                    col = alpha(data$colour,
                                data$alpha),
                    fill = alpha(data$colour, data$alpha),
                    fontsize = 3 * ggplot2::.pt,
                    lwd = 0.94
                )
            )
    }
    else if (kt == "rect") {
        keyGrob <- rectGrob(gp = grid::gpar(col = NA, fill = alpha(data$colour,
                                                             data$alpha)))
    }
    else if (kt == "image") {
        img <-
            image_read(system.file("extdata/Rlogo.png", package = "ggimage"))
        grobs <- lapply(seq_along(data$colour), function(i) {
            img <- color_image(img, data$colour[i], data$alpha[i])
            rasterGrob(0.5,
                       0.5,
                       image = img,
                       width = 1,
                       height = 1)
        })
        class(grobs) <- "gList"
        keyGrob <- ggname("image_key", gTree(children = grobs))
    }
    return(keyGrob)
}

#' geom layer for visualizing image files
#'
#'
#' @title geom_image.rect
#' @param mapping aes mapping
#' @param data data
#' @param stat stat
#' @param position position
#' @param inherit.aes logical, whether inherit aes from ggplot()
#' @param na.rm logical, whether remove NA values
#' @param ... additional parameters
#' @return geom layer
#' @importFrom ggplot2 layer
#' @examples
#' library("ggplot2")
#' library("ggimage")
#' set.seed(2017-02-21)
#' d <- data.frame(
#'     x = rnorm(10),
#'     y = rnorm(10),
#'     image = sample(c("https://www.r-project.org/logo/Rlogo.png",
#'                      "https://jeroenooms.github.io/images/frink.png"),
#'     size=10, replace = TRUE)
#'                )
#' ggplot(d, aes(xmin = x, xmax = 2*x, ymin = y, ymax = 2*y, image=image)) +
#'     geom_rect(fill = "blue") +
#'     geom_image.rect()
#'
#' @author guangchuang yu
geom_image.rect <- function(mapping = NULL,
                            data = NULL,
                            stat = "identity",
                            position = "identity",
                            inherit.aes = TRUE,
                            na.rm = FALSE,
                            ...) {

    # by <- match.arg(by, c("width", "height"))

    layer(
        data=data,
        mapping=mapping,
        geom=GeomImage.rect,
        stat=stat,
        position=position,
        show.legend=NA,
        inherit.aes=inherit.aes,
        params = list(
            na.rm = na.rm,
            ...),
        check.aes = FALSE
    )
}

##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 Geom
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 draw_key_blank
##' @importFrom grid gTree
##' @importFrom grid gList
GeomImage.rect <- ggplot2::ggproto(
    "GeomImage.rect",
    Geom,
    setup_data = function(data, params) {
        if (is.null(data$subset))
            return(data)
        data[which(data$subset), ]
    },

    default_aes = ggplot2::aes(
        image = system.file("extdata/Rlogo.png", package = "ggimage"),
        #size=0.05,
        colour = NULL,
        #angle = 0,
        alpha = 1
    ),

    draw_panel = function(data,
                          panel_params,
                          coord,
                          by,
                          na.rm = FALSE,
                          .fun = NULL,
                          height,
                          image_fun = NULL,
                          # hjust=0.5,
                          # nudge_x = 0, nudge_y = 0,
                          asp = 1) {
        # data$x <- data$x + nudge_x
        # data$y <- data$y + nudge_y
        data <-
            coord$transform(data, panel_params)

        if (!is.null(.fun) &&
            is.function(.fun))
            data$image <- .fun(data$image)

        groups <-
            split(data, factor(data$image))
        imgs <- names(groups)
        grobs <-
            lapply(seq_along(groups), function(i) {
                d <- groups[[i]]
                imageGrob.rect(d$xmin,
                               d$xmax,
                               d$ymin,
                               d$ymax,
                               imgs[i],
                               #by,
                               # hjust,
                               d$colour,
                               d$alpha,
                               image_fun,
                               #d$angle,
                               asp)
            })
        grobs <- do.call("c", grobs)
        class(grobs) <- "gList"

        ggname("geom_image.rect",
               gTree(children = grobs))
    },
    non_missing_aes = c(#"size",
        "image"),
    required_aes = c("xmin", "xmax", "ymin", "ymax"),
    draw_key = draw_key_image
)



##' @importFrom magick image_read
##' @importFrom magick image_read_svg
##' @importFrom magick image_read_pdf
##' @importFrom magick image_transparent
##' @importFrom magick image_rotate
##' @importFrom grid rasterGrob
##' @importFrom grid viewport
##' @importFrom grDevices rgb
##' @importFrom grDevices col2rgb
##' @importFrom methods is
##' @importFrom tools file_ext
imageGrob.rect <- function(xmin, xmax, ymin, ymax, img, #by, hjust,
                           colour, alpha, image_fun, #angle,
                           asp=1) {
    if (!is(img, "magick-image")) {
        if (tools::file_ext(img) == "svg") {
            img <- image_read_svg(img)
        } else if (tools::file_ext(img) == "pdf") {
            img <- image_read_pdf(img)
        } else {
            img <- image_read(img)
        }
        asp <- getAR2(img)/asp
    }

    unit <- "native"
    width = xmax - xmin
    height = ymax - ymin
    # if (any(size == Inf)) {
    #     x <- 0.5
    #     y <- 0.5
    #     width <- 1
    #     height <- 1
    #     unit <- "npc"
    # } else if (by == "width") {
    #     width <- size
    #     height <- size/asp
    # } else {
    #     width <- size * asp
    #     height <- size
    # }
    #
    # if (hjust == 0 || hjust == "left") {
    #     x <- x + width/2
    # } else if (hjust == 1 || hjust == "right") {
    #     x <- x - width/2
    # }

    if (!is.null(image_fun)) {
        img <- image_fun(img)
    }

    if (is.null(colour)) {
        grobs <- list()
        grobs[[1]] <- rasterGrob(x = xmin,
                                 y = ymin,
                                 just = c(0,0),
                                 image = img,
                                 default.units = unit,
                                 height = height,
                                 width = width,
                                 interpolate = FALSE)
    } else {
        cimg <- lapply(seq_along(colour), function(i) {
            color_image(img, colour[i], alpha[i])
        })

        grobs <- lapply(seq_along(xmin), function(i) {
            img <- cimg[[i]]
            # if (angle[i] != 0) {
            #     img <- image_rotate(img, angle[i])
            #     img <- image_transparent(img, "white")
            # }
            rasterGrob(x = xmin[i],
                       y = ymin[i],
                       just = c(0,0),
                       image = img,
                       default.units = unit,
                       height = height[i],
                       width = width[i],
                       interpolate = FALSE
                       ## gp = gpar(rot = angle[i])
                       ## vp = viewport(angle=angle[i])
            )
        })
    }
    return(grobs)
}


##' @importFrom magick image_info
getAR2 <- function(magick_image) {
    info <- image_info(magick_image)
    info$width/info$height
}


compute_just <- getFromNamespace("compute_just", "ggplot2")
