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
first_fcn <- function(test) {
    print(paste("First function says", test))
}
