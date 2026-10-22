## Copyright (C) 2025  Mario Rosasco, Matt Dufort, and Benaroya Research Institute
##
## This file is part of the tcrGraph package

#' Get clone information
#'
#' Given TCR data, identify clones and count the number of occurrences for each clone. A clone is defined as a maximal subgraph
#' containing at most the number of alpha, beta, delta, and gamma chains set by `maxA`, `maxB`, `maxD`, and `maxG`. The arguments
#' `max[A,B,D,G]Lib` and `max[A,B,D,G]Graph` allow for more refined control over the steps at which chains are filtered.
#' `max[A,B,D,G]Lib` arguments cause libraries with more than the allowed number of chains to be removed before clones are counted.
#' `max[A,B,D,G]Graph` arguments cause libraries that support weak edges in violating subgraphs to be split into separate subgraphs
#' until all remaining subgraphs meet the chain-count criteria.
#'
#' @param tcrData Either a tcrGraph object or a data frame containing columns named libid, v_gene, j_gene, and the column
#' indicated by the \code{link} argument.
#' @param maxA the maximum number of alpha chains in a clone. A value of -1 will cause alphas to be removed from consideration. Default: 2.
#' @param maxB the maximum number of beta chains in a clone. A value of -1 will cause betas to be removed from consideration. Default: 1.
#' @param maxD the maximum number of delta chains in a clone. A value of -1 will cause deltas to be removed from consideration. Default: -1.
#' @param maxG the maximum number of gamma chains in a clone. A value of -1 will cause gammas to be removed from consideration. Default: -1.
#'
#' @param maxALib the maximum number of alpha chains in a library. A value of -1 will cause alphas to be removed from consideration. Default: maxA.
#' @param maxBLib the maximum number of beta chains in a library. A value of -1 will cause betas to be removed from consideration. Default: maxB.
#' @param maxDLib the maximum number of delta chains in a library. A value of -1 will cause deltas to be removed from consideration. Default: maxD.
#' @param maxGLib the maximum number of gamma chains in a library. A value of -1 will cause gammas to be removed from consideration. Default: maxG.
#' @param maxAGraph the maximum number of alpha chains in a subgraph. A value of -1 will cause alphas to be removed from consideration. Default: maxA.
#' @param maxBGraph the maximum number of beta chains in a subgraph. A value of -1 will cause betas to be removed from consideration. Default: maxB.
#' @param maxDGraph the maximum number of delta chains in a subgraph. A value of -1 will cause deltas to be removed from consideration. Default: maxD.
#' @param maxGGraph the maximum number of gamma chains in a subgraph. A value of -1 will cause gammas to be removed from consideration. Default: maxG.
#'
#' @param format a string, one of "compressed", "full", or "tcrGraph". In "compressed" format, the returned data frame have one row for each
#' clone, with all of the libs where that clone was found listed together in one column. In "full" format, each lib will be on a separate row,
#' with a column for clone ID. In "tcrGraph" format, a tcrGraph object is returned. Default is "compressed".
#' @param link a string indicating which column(s) to use to match chains. Default is "full_nt_sequence", but only used if \code{tcrData} is
#' a data frame object. If \code{tcrData} is a tcrGraph object, the link used to generate the tcrGraph is used here.
#'
#' @param verbose logical, whether to provide messages on status of clone identification. Default is FALSE.
#' @return A data frame of information about the clones detected in the TCR graph. Column \code{violatesLimits} is TRUE for a
#' clone that exceeds the \code{max[A,B,D,G]Graph} limits but cannot be split, because every lib in it carries every chain
#' (e.g. a single lib with two beta chains when \code{maxBLib = 2} and \code{maxBGraph = 1}). Such clones are kept and flagged,
#' not removed.
#'
#' @importFrom stringr str_replace
#' @export
getClonesFromTcrGraph <- function(
    tcrData,
    maxA = 2,
    maxB = 1,
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
    link = "full_nt_sequence",
    verbose = FALSE
){
  # .logEvent("getClonesFromTcrGraph")

  if(verbose) message("Checking input data.\n")

  # check for input data type. Create a tcrGraph object so that we know we have chainType and id fields.
  if(is.tcrGraph(tcrData)){
    tcrGraph <- tcrData
  } else if(is.data.frame(tcrData)){
    tcrGraph <- makeTcrGraph(tcrData, link)
  } else {
    stop("'tcrData' must be a tcrGraph or data.frame.")
  }
  # extract link here to support multipoint linkage
  link <- tcrGraph$link

  # make sure we have the columns we need
  if(!all(c("libid", "id", "chainType") %in% names(tcrGraph$data)))
    stop("The 'data' object in the tcrGraph class object must have columns named 'libid', 'id', and 'chainType'.")
  if(!all(tcrGraph$data$chainType %in% c("TRA", "TRB", "TRD", "TRG")))
    stop("The chain type for each data point must be one of 'TRA', 'TRB', 'TRD', or 'TRG'")

  cloneDf <- tcrGraph$data

  if(verbose) message("Filtering to cells and chains matching input parameters.\n")

  # filter out excluded chains at library level
  if(maxALib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRA",]
  if(maxBLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRB",]
  if(maxDLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRD",]
  if(maxGLib == -1) cloneDf <- cloneDf[cloneDf$chainType != "TRG",]

  cloneDf <-
    cloneDf %>%
    # count number of alpha, beta, gamma, and delta chains per lib
    dplyr::group_by(.data$libid, .data$chainType) %>%
    dplyr::mutate(libChainTypeCounts = dplyr::n()) %>%
    dplyr::ungroup()

  # find libs with too many chains at the lib level
  libsToRemove <- cloneDf[
    (cloneDf$chainType == "TRA" & cloneDf$libChainTypeCounts > maxALib) |
      (cloneDf$chainType == "TRB" & cloneDf$libChainTypeCounts > maxBLib) |
      (cloneDf$chainType == "TRD" & cloneDf$libChainTypeCounts > maxDLib) |
      (cloneDf$chainType == "TRG" & cloneDf$libChainTypeCounts > maxGLib),
  ]$libid

  # generate data without multiplets to find clones from the remaining libs
  cloneDf <- cloneDf[!cloneDf$libid %in% libsToRemove, ]

  if(verbose)
    message("Filtering complete; removed ", length(libsToRemove), " cells",
            " and ", nrow(tcrGraph$data) - nrow(cloneDf), " total chains.\n")

  # modify node naming from original makeTcrGraph call to avoid confusion / collisions
  cloneDf$idInput <-
    str_replace(cloneDf$id, "node", "nodeOrig")

  if (nrow(cloneDf) == 0L) {
    cloneSummaryData <- data.frame()
    if (format == "tcrGraph") {
      cloneSummaryData <- makeTcrGraph(cloneSummaryData, link = link)
    }
    return(cloneSummaryData)
  }

  if(verbose) message("Splitting tcrGraph object into clones.\n")

  # Each lib starts in a subgraph pool defined by the existing connected components
  # of the graph produced by makeTcrGraph() on the current cloneDf.
  # This step is necessary because the filtering above may have changed the
  # components of the graph from the original input tcrGraph object.
  # Use 'link = "idInput"' because that field has the node identifier from the upstream
  # makeTcrGraph run, so incorporates the chain matching strategy used to build nodes.
  # Leaving this as 'link = "id"' causes problems because ID gets replaced in the
  # output data frame with the new node names.
  initGraph  <- makeTcrGraph(cloneDf, link = "idInput")
  initIgraph <- tcrGraphToIgraph(initGraph)

  # Decompose into connected components
  initComponents <- igraph::decompose(initIgraph)

  # Map each libid to a component-based subgraphId.
  # A lib can touch multiple nodes, but all of its nodes should lie within one component
  # (because nodes are connected via shared libid edges in makeTcrGraph).
  cloneDf$subgraphId <- NA_integer_

  # iterate over components
  for (currComponent in seq_along(initComponents)) {
    # identify nodes in this component
    componentNodes <- igraph::V(initComponents[[currComponent]])$name
    # identify libraries in those nodes
    componentLibs  <- unique(initGraph$data$libid[initGraph$data$id %in% componentNodes])
    # set subgraphId for those libraries
    cloneDf$subgraphId[cloneDf$libid %in% componentLibs] <- currComponent
  }

  # Safety: any lib that didn't map (should be rare) gets its own subgraphId
  if (anyNA(cloneDf$subgraphId)) {
    missingLibs <- unique(cloneDf$libid[is.na(cloneDf$subgraphId)])
    startId <- ifelse(all(is.na(cloneDf$subgraphId)), 0L, max(cloneDf$subgraphId, na.rm = TRUE))
    for (i in seq_along(missingLibs)) {
      cloneDf$subgraphId[cloneDf$libid == missingLibs[i]] <- startId + i
    }
  }

  # helper: does a subGraph violate the max*Graph constraints?
  subGraphViolatesLimits <- function(currSubgraph, graphData) {
    nodeIds <- igraph::V(currSubgraph)$name
    currTcrData <- dplyr::filter(graphData, .data$id %in% nodeIds)

    # count distinct chains of each type in this subGraph
    nUniqueA <- length(unique(currTcrData$id[currTcrData$chainType == "TRA"]))
    nUniqueB <- length(unique(currTcrData$id[currTcrData$chainType == "TRB"]))
    nUniqueD <- length(unique(currTcrData$id[currTcrData$chainType == "TRD"]))
    nUniqueG <- length(unique(currTcrData$id[currTcrData$chainType == "TRG"]))

    tooManyUniqueA <- (maxAGraph != -1) && (nUniqueA > maxAGraph)
    tooManyUniqueB <- (maxBGraph != -1) && (nUniqueB > maxBGraph)
    tooManyUniqueD <- (maxDGraph != -1) && (nUniqueD > maxDGraph)
    tooManyUniqueG <- (maxGGraph != -1) && (nUniqueG > maxGGraph)

    tooManyUniqueA || tooManyUniqueB || tooManyUniqueD || tooManyUniqueG
  }

  # ---------------------------------------------------------------------------
  # Iteratively split subgraphs by moving libs that support the weakest edge
  # of violating components into new subgraph pools.
  # ---------------------------------------------------------------------------
  # this loop runs once, and reruns if any subgraph gets split
  # set a tracker so we don't redo the whole thing from start every time we split a subgraph
  maxSubgraphIdChecked <- 0L
  # subgraph pools holding one component that violates the limits but cannot be split
  # by moving libs (every lib in it carries every chain); checked once, never again
  unsplittableIds <- integer(0)
  repeat {
    if(verbose && exists("changed") && changed)
      message("Split subgraph ", currSubgraphId, " because it violated clone parameters.")

    changed <- FALSE
    subgraphIds <- sort(unique(cloneDf$subgraphId))

    if(verbose) message(max(subgraphIds) - maxSubgraphIdChecked, " remaining subgraphs to check; this may take a little time.\n")

    # this loop iterates over each subgraph, and checks and splits it
    # it should not restart from 1, but should repeat on any subgraph that gets split
    for (currSubgraphId in subgraphIds[subgraphIds > maxSubgraphIdChecked]) {
      maxSubgraphIdChecked <- currSubgraphId - 1L # this allows the process to restart from same point

      # message("Checking subgraph ", currSubgraphId, ".")
      if (currSubgraphId %in% unsplittableIds) next

      subDf <- cloneDf[cloneDf$subgraphId == currSubgraphId, , drop = FALSE]
      if (nrow(subDf) == 0L) next

      subGraph <- makeTcrGraph(subDf, link = "idInput")
      subIgraph <- tcrGraphToIgraph(subGraph)

      # if no edges, nothing to violate
      if (igraph::ecount(subIgraph) == 0L) next

      subComponents <- igraph::decompose(subIgraph)

      for (currSubgraph in subComponents) {
        if (!subGraphViolatesLimits(currSubgraph, subGraph$data))
          next

        # this component violates the chain-count limits; find weakest edge
        if (igraph::ecount(currSubgraph) == 0L)
          next

        weightsCurrSubgraph <- igraph::E(currSubgraph)$weight
        if (is.null(weightsCurrSubgraph)) {
          stop("Edge weights not found on igraph; modify tcrGraphToIgraph() so ",
               "edges have a numeric 'weight' attribute giving support per edge.")
        }
        # if multiple edges with same weight, select the first one listed
        # this approach could be changed to select randomly if there is more than one
        weakestEdgeIndex <- which(weightsCurrSubgraph == min(weightsCurrSubgraph, na.rm = TRUE))[1L]
        edgeEnds <- igraph::ends(currSubgraph, weakestEdgeIndex)
        node1 <- edgeEnds[1]
        node2 <- edgeEnds[2]

        # find libraries (cells) that contain BOTH chains involved in this edge
        edgeData <- dplyr::filter(subGraph$data, .data$id %in% c(node1, node2))
        libs1 <- unique(edgeData$libid[edgeData$id == node1])
        libs2 <- unique(edgeData$libid[edgeData$id == node2])
        libsToMove <- intersect(libs1, libs2)

        if (length(libsToMove) == 0L) {
          # paranoia guard: no libs found for this edge; skip
          next
        }

        newId <- max(cloneDf$subgraphId) + 1L

        # If every lib in this component supports the weakest edge, every lib carries
        # every chain (each edge's weight is at most the number of libs), so moving
        # libs cannot split it: the moved pool would be the same component, violate
        # again, and be moved again forever. Isolate it in its own pool, mark it
        # unsplittable, and move on.
        componentLibs <- unique(subGraph$data$libid[
          subGraph$data$id %in% igraph::V(currSubgraph)$name])
        if (all(componentLibs %in% libsToMove)) {
          cloneDf$subgraphId[
            cloneDf$libid %in% componentLibs & cloneDf$subgraphId == currSubgraphId
          ] <- newId
          unsplittableIds <- c(unsplittableIds, newId)
          if (verbose)
            message("Subgraph ", currSubgraphId, " has a component that violates clone ",
                    "parameters but cannot be split (libs: ",
                    paste(componentLibs, collapse = ", "), "); flagged with violatesLimits.")
          changed <- TRUE
          break
        }

        cloneDf$subgraphId[
          cloneDf$libid %in% libsToMove & cloneDf$subgraphId == currSubgraphId
        ] <- newId

        changed <- TRUE
        break  # break component loop; restart with updated cloneDf
      }

      if (changed) break  # restart from updated subgraph assignments
    }

    if (!changed) break  # no violating components anywhere; we're done splitting
  }

  # ---------------------------------------------------------------------------
  # Now, for each final subgraph pool, decompose into connected components and
  # call each component a clone. Each lib appears in exactly one subgraphId.
  # Chain types with max*Graph == -1 are removed from *output* (not used as
  # constraints or reported chains).
  # This separate step is a secondary check, and also generates summary objects
  # ---------------------------------------------------------------------------
  if(verbose) message("Assigning final clone IDs.")
  cloneSummaryData <- data.frame()
  cloneIndex <- 1L

  finalSubgraphIds <- sort(unique(cloneDf$subgraphId))

  for (currSubgraphId in finalSubgraphIds) {
    subDf <- cloneDf[cloneDf$subgraphId == currSubgraphId, , drop = FALSE]
    if (nrow(subDf) == 0L) next

    subGraph <- makeTcrGraph(subDf, link = "idInput")
    subIgraph <- tcrGraphToIgraph(subGraph)
    subComponents <- igraph::decompose(subIgraph)

    for (currSubgraph in subComponents) {
      currSubDf <- igraph::as_data_frame(currSubgraph, what = "vertices")
      currNodeList <- unique(c(currSubDf$name))

      currTcrData <- dplyr::filter(subGraph$data, .data$id %in% currNodeList)

      # Remove chain types that are not considered at graph level
      if (maxAGraph == -1) currTcrData <- currTcrData[currTcrData$chainType != "TRA", ]
      if (maxBGraph == -1) currTcrData <- currTcrData[currTcrData$chainType != "TRB", ]
      if (maxDGraph == -1) currTcrData <- currTcrData[currTcrData$chainType != "TRD", ]
      if (maxGGraph == -1) currTcrData <- currTcrData[currTcrData$chainType != "TRG", ]

      if (nrow(currTcrData) == 0L)
        next

      if(format == "compressed"){
        currCloneSummary <- dplyr::summarise(
          currTcrData,
          cloneId     = paste0("clone", cloneIndex),
          cloneCounts = length(unique(.data$libid)),
          violatesLimits = currSubgraphId %in% unsplittableIds,
          libs        = paste(unique(.data$libid), collapse = ", "),
          vGenes      = paste(unique(.data$v_gene), collapse = ", "),
          jGenes      = paste(unique(.data$j_gene), collapse = ", "),
          !!link      := paste(unique(.data[[link]]), collapse = ", ")
        )
      } else if (format %in% c("full", "tcrGraph")){
        currCloneSummary <- dplyr::mutate(
          currTcrData,
          cloneCounts = length(unique(currTcrData$libid)),
          cloneId     = paste0("clone", cloneIndex),
          violatesLimits = currSubgraphId %in% unsplittableIds
        )
      } else {
        stop("Invalid 'format' argument. Format must be 'compressed', 'full', or 'tcrGraph'.")
      }

      if(nrow(currTcrData) > 0){
        cloneSummaryData <- rbind(cloneSummaryData, currCloneSummary)
        cloneIndex <- cloneIndex + 1L
      }
    }
  }

  if (format == "tcrGraph"){
    cloneSummaryData <- makeTcrGraph(cloneSummaryData, link = link)
  }

  return(cloneSummaryData)
}