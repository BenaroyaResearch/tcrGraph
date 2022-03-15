## Copyright (C) 2020 Mario Rosasco and Benaroya Research Institute
##
## This file is part of the tcrGraph package

#' Generates an interactive visNetwork plot
#' 
#' Generates a visNetwork plot with reasonable defaults for
#' data exploration. Investigating a TCR network in this way will
#' probably be more useful than using the base plotting function
#' for networks with mrore than a few nodes. Hovering over each
#' node will provide additional information about the chain.
#' 
#' @param tcrGraph the tcrGraph object to generate a plot from
#' @return an interactive visNetwork plot
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
#' makeVisNetwork(myGraph)     
#' }
#' @seealso \code{\link[visNetwork]{visNetwork}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @export
makeVisNetwork <- function(tcrGraph){
  .logEvent("makeVisNetwork")
  
  visNetwork::visNetwork(nodes = tcrGraph$nodes, edges = tcrGraph$edges) %>% 
    visNetwork::visPhysics(enabled = FALSE) %>% 
    visNetwork::visIgraphLayout(layout = "layout_components")
}

#' Generates a circos-style plot
#' 
#' Generates a circos-style plot with reasonable default aesthetics.
#' 
#' @param tcrGraph the tcrGraph object to generate a plot from
#' @param title a string that will be displayed as the title of the plot
#' @param annoFields a string or vector of strings indicating the names of fields in the tcrGraph data 
#' that will be used to generate annotation tracks in the circos plot
#' @param annoPalettes a named list of color palettes, where names correspond to the values in \code{annoFields} (see examples).
#' If a given annotation field doesn't have a palette specified, a random palette will be assigned.
#' @param sortField a string indicating the field in \code{annoFields} to sort the sections in the circos plot by
#' @param sortMixed a logical indicating whether the sortField should be sorted by numerical as well as alphabetical characters
#' @param arcColorField a string indicating which of the fields in \code{tcrGraph} to use to color the arcs connecting sections
#' of the plot. Note that the "from" node sets the color of the arc.
#' @param sampleField a string indicating the field used to represent individual samples/cells in the tcrGraph object. 
#' For BRI data from the Research Database this should almost always be "libid"
#' @param chainField a string indicating the field used to represent a node in the tcrGraph object. For data prepared by
#' @param showLegend a logical indicating whether to render the legend given the \code{annoFields}
#' \code{makeTcrGraph()} this should be left as the default "id".
#' @param trackHeight a number representing the fraction of a unit circle that should be taken up by each track 
#' (eg: '0.1' means each track takes up 1/10 of the full circle radius)
#' #' @return a plot object
#' 
#' @examples 
#' \dontrun{
#' # self-contained example with user-defined palettes
#' tcrData <- data.frame(
#'    libid = c("lib000", "lib000", "lib001", "lib001"), 
#'    v_gene = c("TRAV1-1", "TRBV1", "TRAV1-1", "TRBV2"), 
#'    j_gene = c("TRAJ11", "TRBJ1-1", "TRAJ11", "TRBJ2-2"),
#'    project = c("P0-1", "P0-1", "P0-2", "P0-2"),
#'    fake_nt_sequence = c("CAT", "GAT", "CAT", "GAC")
#' )
#' myGraph <- makeTcrGraph(tcrData, link = "fake_nt_sequence")
#' makeCircos(
#'   myGraph,
#'   annoPalettes = list(
#'     # unnamed vectors of colors with the right length are allowed
#'     project = c("blue", "green"), 
#'     # named vectors of colors are also allowed
#'     chainType = c(TRB = "red", TRA = "purple")
#'   )
#' )
#' 
#' # We can also directly visualize the circos plot for clones data
#' makeCircos(getClonesFromTcrGraph(myGraph, format = "tcrGraph"))
#' }
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange
#' @importFrom grDevices colors
#' @import graphics
#' @import circlize
#' @export
makeCircos <- function(
  tcrGraph, 
  title = "",
  annoFields = c("project", "chainType"),
  annoPalettes = NULL,
  sortField = "project",
  sortMixed = FALSE,
  sampleField = "libid",
  chainField = "id",
  arcColorField = sortField,
  showLegend = TRUE,
  trackHeight = 0.1
){
  # housekeeping
  if(!is.tcrGraph(tcrGraph)){
    stop("The 'tcrGraph' argument must be of class 'tcrGraph'.")
  }
  if(!all(annoFields %in% names(tcrGraph$data))){
    stop("The data in 'tcrGraph' must contain fields for each of the names listed in 'annoFields'.")
  }
  if(!all(c(chainField, sampleField) %in% names(tcrGraph$data))){
    stop("The data in 'tcrGraph' must contain fields for each of the names listed in 'chainField' and 'sampleField'.")
  }
  if(!all(arcColorField %in% names(tcrGraph$data))){
    stop("The data in 'tcrGraph' must contain the field named by 'arcColorField'.")
  }
  if(!sortField %in%  names(tcrGraph$data)){
    stop("The data in 'tcrGraph' must contain the field named by 'sortField'.")
  }
  
  # setup
  df <- as.data.frame(tcrGraph$data)
  annoLegendList <- list()
  
  # set up color palettes
  for (currAnnoField in unique(c(annoFields, arcColorField))) {
    currAnnoValues <- unique(df[[currAnnoField]])
    
    currPaletteSet <- FALSE
    # check if a palette has already been set up for the current annotation
    if(currAnnoField %in% names(annoPalettes)){
      currPalette <- annoPalettes[[currAnnoField]]
      if(all(currAnnoValues %in% names(currPalette))){ # palette set
        currPaletteSet <- TRUE
      } else if(length(currAnnoValues) == length(currPalette)){ # palette set but not named
        names(currPalette) <- currAnnoValues
        currPaletteSet <- TRUE
      }
    } 
    if(!currPaletteSet){ # no palette set - assign random
      message(paste("No palette set for annotation", currAnnoField, "- assigning random colors."))
      currPalette <- sample(colors())[1:length(currAnnoValues)]
      names(currPalette) <- currAnnoValues
    }
    
    # merge colors onto df to make plotting easier
    currPaletteDf <- as.data.frame(currPalette, stringsAsFactors = FALSE)
    currPaletteDf[[currAnnoField]] <- rownames(currPaletteDf)
    colnames(currPaletteDf) <- c(paste0(currAnnoField, "Palette"), currAnnoField)
    df <- merge(df, currPaletteDf, all.x = TRUE, by = currAnnoField)
    
    annoLegendList[[currAnnoField]] <- list(legend=currAnnoValues, fill=currPalette)
  }
  
  # assign circos node IDs and remove duplicate nodes (eg: from second contig with same link values)
  df$.trackCellId <- paste0(df[[sampleField]], "-", df[[chainField]])
  df <- df[!duplicated(df$.trackCellId),]
  # sort data for plotting
  if (sortMixed) {
    df <- df[mmorder(df[[sortField]], df[[sampleField]], df[[chainField]]),]
  } else {
    df <- df[order(df[[sortField]], df[[sampleField]], df[[chainField]]),]
  }
  df$.trackCellId <- factor(df$.trackCellId, levels = df$.trackCellId)
  
  # set up the plotting environment
  message("Generating circos plot. This may take a moment, please wait...")
  
  if (showLegend) {
    nAnno <- length(annoFields)
    layoutMatrix <- t(matrix(c(1,1,2:(nAnno+1), 1,1,2:(nAnno+1)), ncol = 2))
  } else {
    layoutMatrix <- t(matrix(c(1,1,1, 1,1,1)))
  }
  layout(layoutMatrix)
  par(mar=c(0,0,0,0))
  
  # Set up the circos plot parameters
  circos.clear()
  circos.par(
    "gap.degree" = 0,
    "cell.padding" = c(0,0,0,0),
    "track.margin" = c(0, 0),
    "track.height" = trackHeight
  )
  circos.initialize(factors = df$.trackCellId, xlim = c(0,1))
  
  # plot tracks
  for (currAnnoField in annoFields){
    circos.trackPlotRegion(
      factors = df$.trackCellId, 
      bg.col = df[,paste0(currAnnoField, "Palette")],
      bg.border = df[,paste0(currAnnoField, "Palette")],
      ylim = c(0,1)
    )
  }
  
  # draw links
  while(nrow(df) > 0){
    currRow <- df[1,,FALSE]
    df <- df[-1,,FALSE]
    # get nodes that should be linked to this one
    connectedNodes <- df[df[[chainField]] == currRow[[chainField]], ".trackCellId"]
    for(currNode in connectedNodes){
      circos.link(currRow$.trackCellId, 0.5, currNode, 0.5, col = currRow[[paste0(arcColorField, "Palette")]])
    }
  }
  
  # draw the legends
  if (showLegend) {
    # draw the legends
    for (currLegendInfo in annoLegendList){
      plot.new()
      par(mar=c(0,0,0,0), new = TRUE)
      legend(
        "left", 
        legend=currLegendInfo$legend,
        # make sure name order matches
        fill=currLegendInfo$fill[match(currLegendInfo$legend, names(currLegendInfo$fill))], 
        bty = "n"
      )
    }
  }

  # reset plotting window to add title
  if (showLegend) {
    par(
      mfrow = c(1, 1),
      mar = c(5, 4, 4, 2)
    )
  } else {
    par(
      mfrow = c(1, 1),
      mar = c(0, 0, 1, 0)
    )
    

  }
  title(title)
}
