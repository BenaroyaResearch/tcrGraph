library(tcrGraph)
library(testthat)

# GIVEN some TCR data in a data frame with fields named
# libid, v_gene, j_gene, full_nt_sequence, and chainType
testData = data.frame(
  libid = c("lib000", "lib000", "lib001"), 
  id = c("node_1", "node_2", "node_2"),
  full_nt_sequence = c("CAT", "GAT", "GAT"),
  chainType = c("TRA", "TRB", "TRB"),
  stringsAsFactors = FALSE
)
# AND data frames describing the node and edge structure for the TCR graph
testNodes = data.frame(
  id = c("node_1", "node_2"),
  group = c("TRA", "TRB"),
  stringsAsFactors = FALSE
)
testEdges = data.frame(
  from = c("node_1", "node_2"),
  to = c("node_2", "node_1"),
  stringsAsFactors = FALSE
)

context("S3 Class Constructor")
with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  
  test_that("The S3 constructor generates a valid tcrGraph object when given valid data", {
    # GIVEN some valid TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN the tcrGraph S3 constructor is called
    testGraph = tcrGraph(
      data = testData, 
      nodes = testNodes, 
      edges = testEdges,
      link = testLinkField
    )
    # THEN the object returned will be of class tcrGraph
    expect_s3_class(testGraph, "tcrGraph")
  }),
  
  test_that("The S3 constructor throws an error when given data with missing fields", {
    # GIVEN some TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN required fields are missing from the data
    # AND the tcrGraph S3 constructor is called 
    # THEN an error is thrown
    for (i in 1:ncol(testData)){
      expect_error(
        tcrGraph(
          data = testData[,-i], 
          nodes = testNodes, 
          edges = testEdges,
          link = testLinkField
        )
      )
    }
  }),
  
  test_that("The S3 constructor throws an error when given incomplete node information", {
    # GIVEN some TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN required fields are missing from the node information
    # AND the tcrGraph S3 constructor is called 
    # THEN an error is thrown
    for (i in 1:ncol(testNodes)){
      expect_error(
        tcrGraph(
          data = testData, 
          nodes = testNodes[,-i], 
          edges = testEdges,
          link = testLinkField
        )
      )
    }
  }),
  
  test_that("The S3 constructor throws an error when given incomplete edge information", {
    # GIVEN some TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN required fields are missing from the edge information
    # AND the tcrGraph S3 constructor is called 
    # THEN an error is thrown
    for (i in 1:ncol(testEdges)){
      expect_error(
        tcrGraph(
          data = testData, 
          nodes = testNodes, 
          edges = testEdges[,-i],
          link = testLinkField
        )
      )
    }
  }),
  
  test_that("The S3 constructor throws an error when given node groups aren't TCR genes", {
    # GIVEN some TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN a non-standard TCR gene group name in the node information
    testNodesWithIssue = testNodes
    testNodesWithIssue[1,"group"] = "OtherGene"
    # AND the tcrGraph S3 constructor is called 
    # THEN an error is thrown
    expect_error(
      tcrGraph(
        data = testData, 
        nodes = testNodesWithIssue, 
        edges = testEdges,
        link = testLinkField
      )
    )
  }),
  
  test_that("The S3 constructor throws an error when nodes and edges are discordant", {
    # GIVEN some TCR, node, and edge data in data frames
    expect_s3_class(testData, "data.frame")
    expect_s3_class(testNodes, "data.frame")
    expect_s3_class(testEdges, "data.frame")
    # AND a chain link field name that's present in the TCR data
    testLinkField = "full_nt_sequence"
    expect_true(testLinkField %in% names(testData))
    # WHEN the edges data references nodes absent from the node data
    testEdgesWithFromIssue = testEdges
    testEdgesWithFromIssue[1,"from"] = "OtherNode"
    # AND the tcrGraph S3 constructor is called 
    # THEN an error is thrown
    expect_error(
      tcrGraph(
        data = testData, 
        nodes = testNodes, 
        edges = testEdgesWithFromIssue,
        link = testLinkField
      )
    )
  }),
  
  # set env for mocking private function
  .env = "tcrGraph"
)

# GIVEN some TCR data in a data frame with fields named 
# libid, v_gene, j_gene, and full_nt_sequence
testData = data.frame(
  libid = c("lib000", "lib000", "lib001"), 
  v_gene = c("TRAV1-1", "TRBV1", "TRBV2"), 
  j_gene = c("TRAJ11", "TRBJ1-1", "TRBJ2"), 
  full_nt_sequence = c("CAT", "GAT", "GAT"),
  stringsAsFactors = FALSE
)  

context("S3 Class Generic Methods")
with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  
  test_that("The S3 as.data.frame method converts a tcrGraph to a data frame", {
    # GIVEN a valid tcrGraph object
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    # WHEN a generic method is used to coerce the tcrGraph object to a data frame
    testDf = as.data.frame(testGraph)
    # THEN the returned object will be of class data.frame
    expect_s3_class(testDf, "data.frame")
  }),
  
  test_that("The S3 print method prints tcrGraph summary data", {
    # GIVEN a valid tcrGraph object
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    # WHEN a generic method is called to print the tcrGraph object
    # THEN summary data will be printed for the object
    expect_output(print(testGraph), "An object of class 'tcrGraph'")
  }),
  
  # set env for mocking private function
  .env = "tcrGraph"
)