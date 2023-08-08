
list_rename <- function(lst, name_append = "") {
    names(lst) <- paste0(names(lst), name_append)
    return(lst)
}


