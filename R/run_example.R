#' runExample
#'
#' runs an example shiny app demonstrating how the parts fits together.
#'
#' @export
#' @return runs a shiny app
#' @importFrom shiny runApp
runExample <- function() {
    appDir <- system.file("shiny-examples", "example_app", package = "chiptsne")
    if (appDir == "") {
        stop("Could not find example directory. Try re-installing `chiptsne`.", call. = FALSE)
    }

    shiny::runApp(appDir, display.mode = "normal")
}
