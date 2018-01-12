#' @keywords internal
"_PACKAGE"

#' @import dplyr
#' @import ggplot2
#' @importFrom stringr str_c str_interp
#' @importFrom plotly ggplotly
#' @importFrom openxlsx write.xlsx
NULL

# quiets concerns of R CMD check about . that appears in pipelines
# and some very commonly used variable names used in NSE commands
utils::globalVariables(c("."))

# release questions
release_questions <- function() {
  c(
    "Is it passing travis, appveyor and win-builder?"
  )
}
