## Copyright (C) 2019  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the tcrGraph package

#' tcrGraph S3 class
#' 
#' Functions to create and check for tcrGraph objects.
#' 
#' Note that it is usually not advised to manually construct a tcrGraph object. To generate a tcrGraph from a data frame of TCR data, this
#' package provides the \code{\link{makeTcrGraph}} function, which should be used in most cases.
#' 
#' Note also that when a field named "title" is present in tcrGraph nodes and edges data, this field is assumed  
#' to contain summary data as HTML that can be rendered during a mouseover event in a browser. This is useful 
#' for some interactive plotting packages such as \code{\link[visNetwork]{visNetwork}}, but it means you should
#' avoid using 'title' as a field/column name in your TCR data frame.
#' 
#' @usage
#' tcrGraph(data, nodes, edges, link)
#' 
#' ## S3 methods for class 'tcrGraph'
#' is.tcrGraph(x)
#' 
#' @param data data frame containing TCR data used to construct the nodes and edges. This must include columns named \code{libid}, 
#' \code{id}, \code{chainType}, and the column indicated by the \code{link} argument.
#' @param nodes data frame containing columns \code{id} and \code{group}, where group is one of "TRA", "TRB", "TRD", or "TRG"
#' @param edges data frame containing columns \code{from} and \code{to}, where both from and to values appear as node ids
#' @param link string indicating the column name in \code{data} used to define nodes (ie: link between libraries)
#' @param x object to be tested for class tcrGraph or coerced into a data frame
#' 
#' @aliases is.tcrGraph as.data.frame
#' 
#' @seealso \code{\link{plot.tcrGraph}}, \code{\link{print.tcrGraph}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
tcrGraph <- function(data, nodes, edges, link){
  
  # check for required data
  if(!all(c("libid", "id", "chainType", link) %in% names(data)))
    stop(paste("The 'data' object must contain columns named 'libid', 'id', 'chainType', and ", link, "."))
  
  if(!all(c("id", "group") %in% names(nodes)))
    stop("The 'nodes' object must contain columns named 'id' and 'group'.")
  
  if(!all(nodes$group %in% c("TRA", "TRB", "TRD", "TRG")))
    stop("The group assignments for each node must be one of 'TRA', 'TRB', 'TRD', or 'TRG'")
  
  if(!all(c("from", "to") %in% names(edges)))
    stop("The 'edges' object must contain columns named 'from' and 'to'.")
  
  if(!all(c(edges$from, edges$to) %in% nodes$id))
    stop("The 'from' and 'to' values in the 'edges' object must appear in the 'id' field of 'data.' ")
  
  .logEvent("tcrGraph")
  return(structure(list(data = data, nodes = nodes, edges = edges, link = link), class = "tcrGraph"))
}

#' @export
is.tcrGraph <- function(x){ inherits(x, "tcrGraph") }

#' @method as.data.frame tcrGraph
#' @export
as.data.frame.tcrGraph <- function(x, row.names = NULL, ...){
  if(!is.null(row.names)){ warning("Argument row.names is being ignored.")}
  if(!is.null(c(...))){ warning("Additional arguments are being ignored.")}
  if(!is.tcrGraph(x)){ 
    stop("Object must be a tcrGraph object.") 
  }
  if(!"id" %in% names(x$nodes)){
    stop("Node data is missing an id field and the object is no longer a valid tcrGraph class object.")
  }
  if(!"from" %in% names(x$edges)){
    stop("Edge data is missing a from field and the object is no longer a valid tcrGraph class object.")
  }
  .logEvent("as.data.frame.tcrGraph")
  nodesToMerge <- x$nodes[, !(names(x$nodes) %in% "title")]
  edgesToMerge <- x$edges[, !(names(x$edges) %in% "title")]
  names(nodesToMerge) <- paste0(names(nodesToMerge), "OfFromNode")
  names(edgesToMerge) <- paste0(
    names(edgesToMerge), 
    ifelse(grepl("^from$|^to$", names(edgesToMerge)), "", "OfEdge")
  )
  return(
    merge( 
      edgesToMerge,
      nodesToMerge,
      by.x = "from",
      by.y = "idOfFromNode"
    )
  )
}

#' Print method for tcrGraph objects
#' 
#' Prints summary information about the TCR graph structure.
#' 
#' @param x tcrGraph object
#' @param ... other arguments to be passed to generic print method
#' 
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @method print tcrGraph
#' @export
print.tcrGraph <- function(x, ...){
  if(!is.tcrGraph(x))
    stop("The argument to print.tcrGraph must be an object of type tcrGraph")
  cat(
    sprintf(
      "An object of class 'tcrGraph'\nNodes defined by %s\nNumber of nodes: %d\nNumber of edges: %d",
      x$link,
      nrow(x$nodes),
      nrow(x$edges)
    )
  )
  if(length(list(...)) > 0){
    warning("Additional arguments were passed to print.tcrGraph, and were ignored.")
  }
}

#' Plot method for tcrGraph objects
#' 
#' Creates simple plot of TCR graph structure
#' 
#' @param x object of class 'tcrGraph' to be printed
#' @param nodeScalingFactor a number to multiply the number of node occurrences by to generate the node size in the plot. Default 1.
#' @param chainPalette a list of 4 colors to be assigned to chains of type TRA, TRB, TRD, and TRG. Any names will be overridden. 
#' Default rainbow(4)
#' @param ... additional arguments to be passed to igraph.plot(...)
#' 
#' @importFrom grDevices rainbow
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
#' plot(myGraph, nodeScalingFactor = 30)
#' }
#' 
#' @seealso \code{\link[igraph]{igraph.plotting}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @method plot tcrGraph
#' @export
plot.tcrGraph <- function(x, nodeScalingFactor = 1, chainPalette = rainbow(4), ...){
  if(!is.tcrGraph(x))
    stop("The argument to plot.tcrGraph must be an object of type tcrGraph")
  
  .logEvent("plot.tcrGraph")
  names(chainPalette) <- c("TRA", "TRB", "TRD", "TRG")
  
  # leverage the plot.igraph function for simplicity
  igraphToPlot <- tcrGraphToIgraph(x)

  igraph::V(igraphToPlot)$size <- x$nodes$value[match(x$nodes$id, names(igraph::V(igraphToPlot)))] * nodeScalingFactor
  igraph::V(igraphToPlot)$label <- x$nodes$group[match(x$nodes$id, names(igraph::V(igraphToPlot)))]
  igraph::V(igraphToPlot)$color <- chainPalette[match(igraph::V(igraphToPlot)$label, names(chainPalette))]
  igraph::plot.igraph(
    igraphToPlot,
    ...
  )
  
}
