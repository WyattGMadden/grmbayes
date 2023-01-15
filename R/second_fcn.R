#' Print a test string
#'
#' @param test A character vector to print
#'
#' @return A string
#' @export
#'
#' @examples
#' test <- "this is a test"
#' first_fcn(test)
second_fcn <- function(test) {
    print(paste("Second function says", test))
}
