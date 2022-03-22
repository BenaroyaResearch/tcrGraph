## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----tcrGraphSetup------------------------------------------------------------
library(tcrGraph)
myTcrData <- data.frame(
  libid = c("lib000", "lib000"), 
  v_gene = c("TRAV1-1", "TRBV1"),
  j_gene = c("TRAJ11", "TRBJ1-1"), 
  fake_nt_sequence = c("CAT", "GAT")
)
myTcrGraph <- makeTcrGraph(myTcrData, link = "fake_nt_sequence")

## ----showTcrGraph-------------------------------------------------------------
print(myTcrGraph)

plot(myTcrGraph)

## ----displayNicerGraph--------------------------------------------------------
plot(
  myTcrGraph, 
  nodeScalingFactor = 30,
  chainPalette = c("orange", "cyan", "purple", "magenta"),
  main = "TCR Graph Plot"
)

## ----tcrClonesGraphSetup------------------------------------------------------
myTcrCloneData <- data.frame(
  libid = c("lib000", "lib000", "lib001", "lib001", "lib002", "lib002"), 
  v_gene = c("TRAV1-1", "TRBV1", "TRAV1-1", "TRBV1", "TRAV1-2", "TRBV2"),
  j_gene = c("TRAJ11", "TRBJ1-1", "TRAJ11", "TRBJ1-1", "TRAJ12", "TRBJ1-2"), 
  fake_nt_sequence = c("CAT", "GAT", "CAT", "GAT", "CAC", "GAG")
)
myTcrCloneGraph <- makeTcrGraph(myTcrCloneData, link = "fake_nt_sequence")
plot(
  myTcrCloneGraph,
  nodeScalingFactor = 30,
  chainPalette = c("orange", "cyan", "purple", "magenta"),
  main = "TCR Clones Plot"
)

## ----showCloneTable-----------------------------------------------------------
getClonesFromTcrGraph(myTcrCloneGraph)

## ----addSingleChainLib--------------------------------------------------------
myTcrCloneData <- rbind(
  myTcrCloneData, 
  data.frame(libid = "lib003", v_gene = "TRBV2", j_gene = "TRBJ1-2", fake_nt_sequence = "GAG")
)
myTcrCloneGraph <- makeTcrGraph(myTcrCloneData, link = "fake_nt_sequence")
plot(
  myTcrCloneGraph,
  nodeScalingFactor = 30,
  chainPalette = c("orange", "cyan", "purple", "magenta"),
  main = "TCR Clones Plot"
)
getClonesFromTcrGraph(myTcrCloneGraph)

## ----addMultipleCellLib-------------------------------------------------------
myTcrCloneData <- rbind(
  myTcrCloneData, 
  data.frame(
    libid = c("lib004", "lib004", "lib004", "lib004", "lib004"),
    v_gene = c("TRBV2", "TRAV1-2", "TRBV1", "TRAV1-1", "TRAV1-3"),
    j_gene = c("TRBJ1-2", "TRAJ12", "TRBJ1-1", "TRAJ11", "TRAJ13"),
    fake_nt_sequence = c("GAG", "CAC", "GAT", "CAT", "TAC")
  )
)
myTcrCloneGraph <- makeTcrGraph(myTcrCloneData, link = "fake_nt_sequence")
plot(
  myTcrCloneGraph,
  nodeScalingFactor = 30,
  chainPalette = c("orange", "cyan", "purple", "magenta"),
  main = "TCR Clones Plot"
)

## ----clonalAnalysisWithMultipleCellLib----------------------------------------
getClonesFromTcrGraph(myTcrCloneGraph)

## ----filterAtLibLevel---------------------------------------------------------
myTcrCloneData <- rbind(
  myTcrCloneData, 
  data.frame(
    libid = c("lib005", "lib005", "lib005"),
    v_gene = c("TRBV3", "TRAV1-2", "TRBV4"),
    j_gene = c("TRBJ1-3", "TRAJ12", "TRBJ1-4"),
    fake_nt_sequence = c("GTG", "CAC", "GTT")
  )
)
myTcrCloneGraph <- makeTcrGraph(myTcrCloneData, link = "fake_nt_sequence")
getClonesFromTcrGraph(
  myTcrCloneGraph,
  maxALib = 2, 
  maxBLib = 2,
  maxAGraph = nrow(myTcrCloneData), # allow all observed chains to occupy a graph
  maxBGraph = nrow(myTcrCloneData)
)

## ----clonalAnalysisWithFullOutput---------------------------------------------
getClonesFromTcrGraph(myTcrCloneGraph, format = "full")

## ----clonalAnalysisWithGraphOutput--------------------------------------------

clonesToPlot <- getClonesFromTcrGraph(
  myTcrCloneGraph,
  maxALib = 2, 
  maxBLib = 2,
  maxAGraph = nrow(myTcrCloneData), # allow all observed chains to occupy a graph
  maxBGraph = nrow(myTcrCloneData),
  format = "tcrGraph"
)

plot(
  clonesToPlot,
  nodeScalingFactor = 30,
  chainPalette = c("orange", "cyan", "purple", "magenta"),
  main = "TCR Clones Plot"
)

