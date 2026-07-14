#' Launch The CASPEN Shiny Demo
#'
#' Opens a lightweight Shiny application for demonstrating CASPEN on the
#' bundled `caspen_example` dataset. The app is intended for teaching and
#' collaborator walkthroughs, not for long production analyses.
#'
#' @param ... Additional arguments passed to `shiny::runApp()`, for example
#' `port`, `host`, or `launch.browser`.
#'
#' @return Invisibly returns the result of `shiny::runApp()`.
#'
#' @examples
#' \dontrun{
#' launch_caspen_demo()
#' launch_caspen_demo(port = 3838, launch.browser = TRUE)
#' }
#'
#' @export
launch_caspen_demo <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The shiny package is required. Install it with install.packages('shiny').")
  }
  app.dir <- system.file("shiny", "caspen-demo", package = "CASPEN")
  if (!nzchar(app.dir)) {
    stop("Could not find the packaged CASPEN Shiny demo app.")
  }
  shiny::runApp(app.dir, ...)
}
