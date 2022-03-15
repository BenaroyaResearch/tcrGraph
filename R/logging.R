## Copyright (C) 2019  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the bRi package

.loggingErrHandler <- function(e){
  write(paste("Could not log event:", e), stderr())
}

.logEvent <- function(event){
  tryCatch(
    {
      httr::POST(
        url = "https://researchdb.benaroyaresearch.org/shineomatic/api/public/logging",
        body = list(logData = list(logUser = Sys.info()[["user"]], event = event, app = "tcrGraph")),
        encode = "json",
        httr::timeout(0.7)
      )
    },
    error = .loggingErrHandler,
    warning = .loggingErrHandler
  )
}
