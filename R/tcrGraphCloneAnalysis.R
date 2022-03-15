## Copyright (C) 2019  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the tcrGraph package

#' Get clone information
#' 
#' Given TCR data, identify clones and count the number of occurrences for each clone. A clone is defined as a maximal subgraph
#' containing at most the number of alpha, beta, delta, and gamma chains set by `maxA`, `maxB`, `maxD`, and `maxG`. The arguments
#' `max[A,B,D,G]Lib` and `max[A,B,D,G]Graph` allow for more refined control over the steps at which chains are filtered. 
#' `Max[A,B,D,G]Lib` arguments cause libraries with more than the allowed number of chains to be removed before clones are counted.
#' `Max[A,B,D,G]Graph` arguments cause the nodes with the fewest number of counts to be removed from each maximal subgraph until 
#' the maximum number of allowed chains remain. This subgraph is then counted as a clone.
#' 
#' @param tcrData Either a tcrGraph object or a data frame containing columns named libid, v_gene, j_gene, and the column 
#' indicated by the \code{link} argument.
#' @param maxA the maximum number of alpha chains in a clone. A value of -1 will cause alphas to be removed from consideration. Default: 2.
#' @param maxB the maximum number of beta chains in a clone. A value of -1 will cause betas to be removed from consideration. Default: 2.
#' @param maxD the maximum number of delta chains in a clone. A value of -1 will cause deltas to be removed from consideration. Default: -1.
#' @param maxG the maximum number of gamma chains in a clone. A value of -1 will cause gammas to be removed from consideration. Default: -1.
#' 
#' @param maxALib the maximum number of alpha chains in a library A value of -1 will cause alphas to be removed from consideration. Default: maxA.
#' @param maxBLib the maximum number of beta chains in a library A value of -1 will cause betas to be removed from consideration. Default: maxB.
#' @param maxDLib the maximum number of delta chains in a library A value of -1 will cause deltas to be removed from consideration. Default: maxD.
#' @param maxGLib the maximum number of gamma chains in a library A value of -1 will cause gammas to be removed from consideration. Default: maxG.
#' @param maxAGraph the maximum number of alpha chains in a subgraph A value of -1 will cause alphas to be removed from consideration. Default: maxA.
#' @param maxBGraph the maximum number of beta chains in a subgraph A value of -1 will cause betas to be removed from consideration. Default: maxB.
#' @param maxDGraph the maximum number of delta chains in a subgraph A value of -1 will cause deltas to be removed from consideration. Default: maxD.
#' @param maxGGraph the maximum number of gamma chains in a subgraph A value of -1 will cause gammas to be removed from consideration. Default: maxG.
#' 
#' @param format a string, one of "compressed", "full", or "tcrGraph". In "compressed" format, the returned data frame have one row for each 
#' clone, with all of the libs where that clone was found listed together in one column. In "full" format, each lib will be on a separate row, 
#' with a column for clone ID. In "tcrGraph" format, a tcrGraph object is returned. Default is "compressed".
#' @param link a string indicating which column to use to match chains. Default is "full_nt_sequence", but only used if \code{tcrData} is
#' a data frame object.
#' @return A data frame of information about the clones detected in the TCR graph
#' 
#' @examples 
#' \dontrun{
#' tcrData <- data.frame(
#'    libid = c("lib000", "lib000"), 
#'    v_gene = c("TRAV1-1", "TRBV1"), 
#'    j_gene = c("TRAJ11", "TRBJ1-1"), 
#'    fake_nt_sequence = c("CAT", "GAT")
#' )
#' myGraph <- makeTcrGraph(tcrData, link = "fake_nt_sequence")
#' getClonesFromTcrGraph(myGraph)
#' # A tibble: 1 x 6
#' #  cloneId cloneCounts libs   vGenes         jGenes          fake_nt_sequence
#' #  <chr>         <int> <chr>  <chr>          <chr>           <chr>           
#' #  clone1            1 lib000 TRAV1-1, TRBV1 TRAJ11, TRBJ1-1 CAT, GAT  
#' }
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @export
getClonesFromTcrGraph <- function(
  tcrData, 
  maxA = 2, 
  maxB = 2, 
  maxD = -1, 
  maxG = -1, 
  maxALib = maxA,
  maxBLib = maxB,
  maxDLib = maxD,
  maxGLib = maxG,
  maxAGraph = maxA,
  maxBGraph = maxB,
  maxDGraph = maxD,
  maxGGraph = maxG,
  format = "compressed", 
  link = "full_nt_sequence"
){
  .logEvent("getClonesFromTcrGraph")
  # check for input data type. Create a tcrGraph object so that we know we have chainType and id fields.
  if(is.tcrGraph(tcrData)){
    tcrGraph <- tcrData
  } else if(is.data.frame(tcrData)){
    tcrGraph <- makeTcrGraph(tcrData, link)
  }
  # extract link here to support multipoint linkage
  link <- tcrGraph$link
  
  # make sure we have the columns we need
  if(!all(c("libid", "id", "chainType") %in% names(tcrGraph$data)))
    stop("The 'data' object in the tcrGraph class object must have columns named 'libid', 'id', and 'chainType'.")
  if(!all(tcrGraph$data$chainType %in% c("TRA", "TRB", "TRD", "TRG")))
    stop("The chain type for each data point must be one of 'TRA', 'TRB', 'TRD', or 'TRG'")
  
  cloneDf <- tcrGraph$data
  
  # filter out excluded chains or count them
  if(maxALib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRA",]
  if(maxBLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRB",]
  if(maxDLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRD",]
  if(maxGLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRG",]
  
  cloneDf <- 
    cloneDf %>%
    # count number of alpha, beta, gamma, and delta chains
    dplyr::group_by(.data$libid, .data$chainType) %>%
    dplyr::mutate(libChainTypeCounts = dplyr::n()) %>%
    dplyr::ungroup()
  
  # find libs with too many chains
  libsToRemove <- cloneDf[
    (cloneDf$chainType == "TRA" & cloneDf$libChainTypeCounts > maxALib) |
      (cloneDf$chainType == "TRB" & cloneDf$libChainTypeCounts > maxBLib) |
      (cloneDf$chainType == "TRD" & cloneDf$libChainTypeCounts > maxDLib) |
      (cloneDf$chainType == "TRG" & cloneDf$libChainTypeCounts > maxGLib),
    ]$libid
  
  # generate a graph without multiplets to find clones from the remaining libs
  cloneDf <- cloneDf[!cloneDf$libid %in% libsToRemove, ]
  #colnames(cloneDf)[colnames(cloneDf)== "id"] <- "linkField"
  cloneGraph <- makeTcrGraph(cloneDf, link = "id")
  
  # Create an igraph object and use subgraph finding tools
  cloneIgraph <- tcrGraphToIgraph(cloneGraph)
  cloneSubgraphs <- igraph::decompose(cloneIgraph)
  
  cloneSummaryData <- data.frame()
  cloneIndex <- 1
  for (currSubgraph in cloneSubgraphs){
    currSubDf <- igraph::as_data_frame(currSubgraph, what = "vertices")
    currNodeList <- unique(c(currSubDf$name))
    
    # Generate clone's summary data for the clone df
    currTcrData <- dplyr::filter(cloneGraph$data, cloneGraph$data$id %in% currNodeList)
    
    # filter nodes based on the max allowed chains of each type
    filterSubgraphNodes <- function(data, maxType, type){
      if(maxType == -1) {
        data <- data[data$chainType != type,]
      } else {
        # retrieve all nodes of given type, then select the most counted up to the max allowed
        tmpData <- 
          data[data$chainType == type, c("id", "chainCounts")] %>%
          unique() %>%
          dplyr::arrange(.data$chainCounts)
        delta <- max((nrow(tmpData) - maxType), 0)
        idsToRemove <- tmpData[0:delta, "id", TRUE]
        data <- data[!data$id %in% idsToRemove,]
      }
      return(data)
    }
    currTcrData <- filterSubgraphNodes(currTcrData, maxAGraph, "TRA")
    currTcrData <- filterSubgraphNodes(currTcrData, maxBGraph, "TRB")
    currTcrData <- filterSubgraphNodes(currTcrData, maxDGraph, "TRD")
    currTcrData <- filterSubgraphNodes(currTcrData, maxGGraph, "TRG")
    
    if(format == "compressed"){
      currCloneSummary <- dplyr::summarise(
        currTcrData,
        cloneId = paste0("clone", cloneIndex),
        cloneCounts = length(unique(.data$libid)),
        libs = paste(unique(.data$libid), collapse = ", "),
        vGenes = paste(unique(.data$v_gene), collapse = ", "), 
        jGenes = paste(unique(.data$j_gene), collapse = ", "),
        !!link := paste(unique(.data[[link]]), collapse = ", ")
      )
    } else if (format %in% c("full", "tcrGraph")){
      currCloneSummary <- dplyr::mutate(
        currTcrData,
        cloneCounts = length(unique(currTcrData$libid)),
        cloneId = paste0("clone", cloneIndex)
      )
    } else {
      stop("Invalid 'format' argument. Format must be 'compressed', 'full', or 'tcrGraph'.")
    }
    if(nrow(currTcrData) > 0){
      cloneSummaryData <- rbind(cloneSummaryData, currCloneSummary)
      cloneIndex <- cloneIndex + 1
    }
  }
  if (format == "tcrGraph"){
    cloneSummaryData <- makeTcrGraph(cloneSummaryData, link = link)
  }
  
  return(cloneSummaryData)
}