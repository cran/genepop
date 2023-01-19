# https://stackoverflow.com/questions/37830819/developing-shiny-app-as-a-package-and-deploying-it-to-shiny-server

#' @name GUI
#' @title Call an experimental GUI for Genepop 
#' @return The return value of a `shiny::runApp()` call.
GUI <- function() {
  appDir <- system.file("genepop-shiny", package = "genepop")
  if (appDir == "") {
    stop("Could not find shiny app. Try re-installing `genepop`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}