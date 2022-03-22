## Copyright (C) 2019  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the bRi package

.loggingErrHandler <- function(e){
  write(paste("Could not log event:", e), stderr())
}

.logEvent = function(event){return(TRUE)}
