library(tcrGraph)
library(testthat)

context("Clone Counting")

with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  
  test_that("getClonesFromTcrGraph() counts the correct number of clones when filtering isn't necessary", {
    # GIVEN some TCR data
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001"),
      v_gene = c("TRAV1-1", "TRBV1", "TRBV2", "TRAV1"),
      j_gene = c("TRAJ11", "TRBJ1-1", "TRBJ1", "TRAJ11"),
      full_nt_sequence = c("CAT", "GAT", "GAT", "CAT")
    )
    
    # GIVEN a tcrGraph object describing one clone observed twice 
    # with one alpha and one beta chain in the graph
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 1
    testClonesCts = 2
    # WHEN clones are counted with default filtering parameters
    testClones = getClonesFromTcrGraph(testGraph)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
  
    # GIVEN a valid tcrGraph object with two clones observed once
    # with one alpha and one beta chain in each graph
    testGraph = makeTcrGraph(testData, link = "v_gene")
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 2
    testClonesCts = c(1,1)
    # WHEN clones are counted with default filtering parameters
    testClones = getClonesFromTcrGraph(testGraph)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
  }),

  test_that("getClonesFromTcrGraph() counts the correct number of clones when filtering at the graph level", {
    # GIVEN a tcrGraph object describing one clone observed twice
    # with one alpha and two beta chains
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001"),
      v_gene = c("TRAV1-1", "TRBV1", "TRBV1", "TRAV1"),
      j_gene = c("TRAJ11", "TRBJ1-1", "TRBJ1", "TRAJ11"),
      full_nt_sequence = c("CAT", "GAT", "GAA", "CAT")
    )
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 1
    testClonesCts = 2
    # WHEN clones are counted allowing only one beta chain per clone
    testClones = getClonesFromTcrGraph(testGraph, maxBGraph = 1)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
    #############################################################
    
    # GIVEN a tcrGraph object describing two beta chains observed once 
    # in two independent libraries and both together in one library 
    testData = data.frame(
      libid = c("lib000", "lib001", "lib002", "lib002"),
      v_gene = c("TRBV1", "TRBV2", "TRBV1", "TRBV2"),
      j_gene = c("TRBJ1-1", "TRBJ1", "TRBJ1-1", "TRBJ1"),
      full_nt_sequence = c("GAT", "GAA", "GAT", "GAA")
    )
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 2
    testClonesCts = c(1,1)
    # WHEN clones are counted allowing only one beta chain per library
    testClones = getClonesFromTcrGraph(testGraph, maxBLib = 1)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
    #############################################################
    
    # GIVEN a tcrGraph object describing two alpha/beta clones, 
    # observed once each and connected by a common delta chain 
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001", "lib000", "lib001"),
      v_gene = c("TRAV1-1", "TRBV1", "TRBV1", "TRAV1", "TRDV1", "TRDV1"),
      j_gene = c("TRAJ11", "TRBJ1-1", "TRBJ1", "TRAJ12", "TRDJ1", "TRDJ1"),
      full_nt_sequence = c("CAT", "GAT", "GAA", "TAC", "TTT", "TTT")
    )
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 2
    testClonesCts = c(1,1)
    # WHEN clones are counted with default values, which exclude delta chains
    testClones = getClonesFromTcrGraph(testGraph)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
    #############################################################
    
  }),
  
  # set env for mocking private function
  .env = "tcrGraph"
)
  