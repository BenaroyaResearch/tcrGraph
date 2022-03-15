library(tcrGraph)
library(testthat)

# GIVEN some TCR data in a data frame with fields named 
# libid, v_gene, j_gene, and full_nt_sequence
testData = data.frame(
  libid = c("lib000", "lib000", "lib001"), 
  v_gene = c("TRAV1-1", "TRBV1", "TRBV2"), 
  j_gene = c("TRAJ11", "TRBJ1-1", "TRBJ2"), 
  full_nt_sequence = c("CAT", "GAT", "GAT"),
  stringsAsFactors = FALSE
)

context("Graph Construction Utility")
with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},

  test_that("makeTcrGraph() generates a valid tcrGraph object with default chain link", {
    # GIVEN some valid TCR data in a data frame
    expect_is(testData, "data.frame")
    # WHEN the data are passed to makeTcrGraph()
    testGraph = makeTcrGraph(testData)
    # THEN the return object is of type 'tcrGraph'
    expect_is(testGraph, "tcrGraph")
  }),
  
  test_that("makeTcrGraph() throws an error when required TCR data is missing", {
    # GIVEN some TCR data in a data frame
    expect_is(testData, "data.frame")
    for(i in 1:ncol(testData)){
      # AND the data are missing one of the required fields
      testDataMissingField = testData[,-i]
      # WHEN the data are passed to makeTcrGraph()
      # THEN the function will throw an error 
      expect_error(makeTcrGraph(testDataMissingField))
    }
  }),
  
  test_that("makeTcrGraph() generates a valid tcrGraph object with user-defined chain link", {
    # GIVEN some valid TCR data in a data frame
    expect_is(testData, "data.frame")
    # AND a user-defined chain link field in the data frame
    testLinkField = "v_gene"
    expect_true(testLinkField %in% names(testData))
    # WHEN the data and the link field are passed to makeTcrGraph()
    testGraph = makeTcrGraph(testData, link = testLinkField)
    # THEN the return object is of type 'tcrGraph'
    expect_is(testGraph, "tcrGraph")
  }),
  
  test_that("makeTcrGraph() generates a valid tcrGraph object with a multipoint chain link", {
    # GIVEN some valid TCR data in a data frame
    expect_is(testData, "data.frame")
    # AND a set of user-defined chain link fields in the data frame
    testLinkFields = c("v_gene", "j_gene", "full_nt_sequence")
    expect_true(all(testLinkFields %in% names(testData)))
    # WHEN the data and the link fields are passed to makeTcrGraph()
    testGraph = makeTcrGraph(testData, link = testLinkFields)
    # THEN the return object is of type 'tcrGraph'
    expect_is(testGraph, "tcrGraph")
  }),
  # set env for mocking private function
  .env = "tcrGraph"
)
  
context("Graph Conversion")
with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  
  test_that("tcrGraphToIgraph() generates a valid igraph object", {
    # GIVEN a valid tcrGraph object
    testGraph = makeTcrGraph(testData)
    expect_is(testGraph, "tcrGraph")
    # WHEN the tcrGraph object is passed to tcrGraphToIgraph()
    testIGraph = tcrGraphToIgraph(testGraph)
    # THEN the return object is of type igraph
    expect_is(testIGraph, "igraph")
  }),
  
  test_that("tcrGraphToGraphNEL() generates a valid graphNEL object", {
    # GIVEN a valid tcrGraph object
    testGraph = makeTcrGraph(testData)
    expect_is(testGraph, "tcrGraph")
    # WHEN the tcrGraph object is passed to tcrGraphToGraphNEL()
    testGraphNel = tcrGraphToGraphNEL(testGraph)
    # THEN the return object is of type igraph
    expect_is(testGraphNel, "graphNEL")
  }),
  # set env for mocking private function
  .env = "tcrGraph"
)
  
context("Graph Statistics")
with_mock(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  
  test_that("tcrGraphDegree() produces the correct node degree with a default chain link", {
    # GIVEN a valid tcrGraph object containing two chains defined with the default chain link
    testGraph = makeTcrGraph(testData)
    expect_is(testGraph, "tcrGraph")
    # WHEN the graph is passed to tcrGraphDegree()
    testDegrees = tcrGraphDegree(testGraph)
    # THEN the return object will be a data frame with a "degree" column
    expect_is(testDegrees, "data.frame")
    expect_true("degree" %in% names(testDegrees))
    # AND the degree of both chain nodes will be correct
    expect_equal(testDegrees$degree, c(1,1))
  }),
  
  test_that("tcrGraphDegree() produces the correct node degree with a user-defined chain link", {
    # GIVEN a valid tcrGraph object containing three chains defined with a user-set chain link
    testLinkField = "v_gene"
    testGraph = makeTcrGraph(testData, link = testLinkField)
    expect_is(testGraph, "tcrGraph")
    # WHEN the graph is passed to tcrGraphDegree()
    testDegrees = tcrGraphDegree(testGraph)
    # THEN the return object will be a data frame with a "degree" column
    expect_is(testDegrees, "data.frame")
    expect_true("degree" %in% names(testDegrees))
    # AND the degree of each chain nodes will be correct
    expect_equal(testDegrees$degree, c(1,1,0))
  }),
  # set env for mocking private function
  .env = "tcrGraph"
)