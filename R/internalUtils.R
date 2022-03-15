# Perform a multiple, mixed ordering, using intermixed alphabetical and numerical keys
# 
# intended as a local/private function only, no export
# 
mmorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=stringr::str_sort(unique(l), numeric = T))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}
