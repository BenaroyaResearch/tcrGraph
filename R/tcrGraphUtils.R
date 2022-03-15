## Copyright (C) 2019  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the tcrGraph package

#' Make a tcrGraph object
#' 
#' Accepts a dataframe containing TCR information, and generates a tcrGraph object
#' 
#' @param tcrDf a data frame that contains columns libid, v_gene, j_gene, and the column indicated by the \code{link} argument
#' @param link a string or vector of strings indicating which column(s) to use to match chains. Default is "full_nt_sequence"
#' @return an S3 object of class \code{tcrGraph}
#' 
#' @examples 
#' \dontrun{
#' tcrData <- data.frame(
#'    libid = c("lib000", "lib000"), 
#'    v_gene = c("TRAV1-1", "TRBV1"), 
#'    j_gene = c("TRAJ11", "TRBJ1-1"), 
#'    fake_nt_sequence = c("CAT", "GAT")
#' )
#' makeTcrGraph(tcrData, link = "fake_nt_sequence")
#' # An object of class 'tcrGraph'
#' # Nodes defined by fake_nt_sequence
#' # Number of nodes: 2
#' # Number of edges: 2
#' 
#' makeTcrGraph(tcrData, link = c("v_gene", "fake_nt_sequence"))
#' # An object of class 'tcrGraph'
#' # Nodes defined by v_gene_fake_nt_sequence
#' # Number of nodes: 2
#' # Number of edges: 2
#' }
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @export
makeTcrGraph <- function(tcrDf, link = "full_nt_sequence"){
  .logEvent("makeTcrGraph")
  # set up evaluation environment
  currEnv <- rlang::pkg_env("tcrGraph")

  # if a multi-point link is requested, set up a column for this
  if(length(link) > 1){
    collapsedLink <- paste0(link, collapse = "_")
    tcrDf[[collapsedLink]] <- apply(tcrDf[,link], 1, function(s){paste(s, collapse = "_")})
    link <- collapsedLink
  }
  
  # collect chain info, node counts, and assign node IDs
  tcrDf <-
    tcrDf %>%
    dplyr::mutate(
      chainType = stringr::str_extract(.data$v_gene, "TR[A-Z]")
    ) %>%
    dplyr::group_by(!!rlang::parse_quo(link, currEnv)) %>%
    dplyr::mutate(
      id = paste0("node_", dplyr::cur_group_id()),
      chainCounts = dplyr::n()
    ) %>%
    dplyr::ungroup()

  # set up a node data frame
  nodes <-
    tcrDf %>%
    dplyr::group_by(.data$id) %>%
    # get unique chains, and summarize info we want to show in the tooltip
    dplyr::summarise(
      vGenes = paste(unique(.data$v_gene), collapse = ", "),
      jGenes = paste(unique(.data$j_gene), collapse = ", "),
      libs = paste(unique(.data$libid), collapse = ", "),
      chainTypes = paste(unique(.data$chainType), collapse = ", "),
      groupField = !!link,
      groupValue = paste(unique(.data[[link]])),
      value = dplyr::n()
    ) %>%
    # create html "title" to show in tool tip
    dplyr::mutate(
      title = paste0("<b>", .data$id, " - ", .data$chainTypes, "</b><br>",
        "Counts: ", .data$value, "<br>",
        "V Gene: ", .data$vGenes, "<br>",
        "J Gene: ", .data$jGenes, "<br>",
        "Libraries: ", stringr::str_trunc(.data$libs, width = 20), "<br>",
        "Defined by ", .data$groupField, " with value ", .data$groupValue, "<br>"
      ),
      !!link := .data$groupValue
    ) %>%
    dplyr::select(
      .data$id, 
      .data$value, 
      group = .data$chainTypes,
      .data$vGenes,
      .data$jGenes,
      !!link,
      .data$libs,
      .data$title
    )

  edges <-
    dplyr::full_join(tcrDf, tcrDf, by = "libid") %>%
    dplyr::select(from = .data$id.x, to = .data$id.y) %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::mutate(width = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(title = paste0("<b>",.data$from, "<br>", .data$to, "<br></b>Number of co-occurences: ", .data$width)) %>%
    dplyr::filter(.data$from != .data$to) %>%
    unique()

  tcrGraphResult <- tcrGraph(data = tcrDf, nodes = nodes, edges = edges, link = link)
  return(tcrGraphResult)
}

#' Convert a tcrGraph obect to a graphNEL object 
#' 
#' From the graphNEL class description: "This is a class of graphs that are 
#' represented in terms of nodes and an edge list. This is a suitable 
#' representation for a graph with a large number of nodes and relatively few edges."
#' 
#' Converting a tcrGraph to a graphNEL object enables the use of utility and 
#' statistics functions from the \code{graph} package.
#' 
#' Note that tcrGraph does not rely on the \code{graph} package, because it is unavailable in CRAN.
#' However, the graph package may be useful for performing operations on TCR network graphs, and so this
#' function is included for portability.
#' 
#' @param tcrGraph the tcrGraph object to convert
#' @return an object of class 'graphNEL'
#' 
#' @seealso \code{\link[graph]{graphNEL}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @export
tcrGraphToGraphNEL <- function(tcrGraph){
  if(!is.tcrGraph(tcrGraph)){stop("Object must be a tcrGraph object.")}

  .logEvent("tcrGraphToGraphNEL")
  # use igraph as a bridge to mitigate dependence on graph package
  tcrIgraph <- tcrGraphToIgraph(tcrGraph)
  graphNel <- igraph::as_graphnel(tcrIgraph)
  return(graphNel)
}

#' Convert a tcrGraph obect to an igraph object 
#' 
#' Converting a tcrGraph to an igraph object enables the use of utility and 
#' statistics functions from the \code{link[igraph]{igraph}} package.
#' 
#' @param tcrGraph the tcrGraph object to convert
#' @return an igraph graph object
#' 
#' @seealso \code{link[igraph]{igraph}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#' @export
tcrGraphToIgraph <- function(tcrGraph){
  if(!is.tcrGraph(tcrGraph)){stop("Object must be a tcrGraph object.")}
  
  .logEvent("tcrGraphToIgraph")
  # create an igraph object
  tcrIgraph <- igraph::graph_from_data_frame(
    d = tcrGraph$edges[,!names(tcrGraph$edges) %in% "title"], 
    directed = FALSE, 
    vertices = tcrGraph$nodes[,!names(tcrGraph$nodes) %in% "title"]
  )
  # simplify (remove duplicate edges and loops)
  tcrIgraph <- igraph::simplify(
    graph = tcrIgraph,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = function(x){unique(x)}
  )
  
  return(tcrIgraph)
}

#' Calculate the degree for each node of a tcrGraph obect
#' 
#' The degree of a node is the number of connections it has to other nodes.
#' 
#' @param tcrGraph the tcrGraph object to calculate the degree of
#' @return a data frame containing the degree for each node id
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
#' tcrGraphDegree(myGraph) 
#'        degree     id
#' node_1      1 node_1      
#' node_2      1 node_2        
#' }
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
tcrGraphDegree <- function(tcrGraph){
  .logEvent("tcrGraphDegree")
  tcrIgraph <- tcrGraphToIgraph(tcrGraph)
  tcrIgraphDeg <- igraph::degree(tcrIgraph)
  return(data.frame(degree = tcrIgraphDeg, id = names(tcrIgraphDeg)))
}

