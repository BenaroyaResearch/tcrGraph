library(tcrGraph)
library(testthat)

context("Clone Counting")

with_mocked_bindings(
  # mock the internal function that calls the API to prevent network requirement
  .logEvent = function(event){return(TRUE)},
  .package = "tcrGraph",
  code = {

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
  })

  test_that("getClonesFromTcrGraph() counts the correct number of clones when filtering at the graph level", {
    # GIVEN a tcrGraph object describing one clone observed twice
    # with two alpha and one beta chains
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001"),
      v_gene = c("TRAV1-1", "TRBV1", "TRBV1", "TRAV1"),
      j_gene = c("TRAJ11", "TRBJ1", "TRBJ1", "TRAJ11"),
      full_nt_sequence = c("CAT", "GAT", "GAT", "CAA")
    )
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 2
    testClonesCts = c(1,1)
    # WHEN clones are counted allowing only one beta chain per clone
    testClones = getClonesFromTcrGraph(testGraph, maxAGraph = 1)
    # THEN the return object is a data frame
    expect_s3_class(testClones, "data.frame")
    # AND the correct number of clones and clone counts is reported
    expect_equal(nrow(testClones), testNClones)
    expect_equal(testClones$cloneCounts, testClonesCts)
    #############################################################

    # GIVEN a tcrGraph object describing two alpha chains observed once
    # in two independent libraries and both together in one library
    testData = data.frame(
      libid = c("lib000", "lib001", "lib002", "lib002"),
      v_gene = c("TRAV1", "TRAV2", "TRAV1", "TRAV2"),
      j_gene = c("TRAJ1-1", "TRAJ1", "TRAJ1-1", "TRAJ1"),
      full_nt_sequence = c("GAT", "GAA", "GAT", "GAA")
    )
    testGraph = makeTcrGraph(testData)
    expect_s3_class(testGraph, "tcrGraph")
    testNClones = 1
    testClonesCts = 3
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

  })

  test_that("getClonesFromTcrGraph() terminates and flags a component that cannot be split", {
    # GIVEN one clone seen in two libs (alpha + beta), and a third lib with that clone's
    # chains plus a second beta; libs may hold two betas but clones only one
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001", "lib002", "lib002", "lib002"),
      v_gene = c("TRAV3", "TRBV2", "TRAV3", "TRBV2", "TRAV3", "TRBV2", "TRBV3"),
      j_gene = c("TRAJ3", "TRBJ2", "TRAJ3", "TRBJ2", "TRAJ3", "TRBJ2", "TRBJ3"),
      full_nt_sequence = c("CCC", "DDD", "CCC", "DDD", "CCC", "DDD", "EEE")
    )
    # WHEN clones are counted (this looped forever before the fix)
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData), maxBLib = 2)
    # THEN the two-chain clone keeps lib000 and lib001, unflagged
    # AND lib002 is split off on its own and flagged
    expect_equal(nrow(testClones), 2)
    expect_equal(sort(testClones$cloneCounts), c(1, 2))
    expect_equal(testClones$violatesLimits[testClones$cloneCounts == 2], FALSE)
    expect_equal(testClones$violatesLimits[testClones$cloneCounts == 1], TRUE)
    expect_equal(testClones$libs[testClones$violatesLimits], "lib002")

    # GIVEN two identical libs that both carry two betas, and nothing else
    testData = data.frame(
      libid = rep(c("lib000", "lib001"), each = 3),
      v_gene = rep(c("TRAV3", "TRBV2", "TRBV3"), 2),
      j_gene = rep(c("TRAJ3", "TRBJ2", "TRBJ3"), 2),
      full_nt_sequence = rep(c("CCC", "DDD", "EEE"), 2)
    )
    # WHEN clones are counted in "full" format
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData), maxBLib = 2, format = "full")
    # THEN both libs form one flagged clone
    expect_equal(length(unique(testClones$cloneId)), 1)
    expect_true(all(testClones$violatesLimits))
  })

  test_that("getClonesFromTcrGraph() flags clones whose extra chains never occur together in one lib", {
    # GIVEN two alphas that share a beta, but no lib carries both alphas
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001", "lib002"),
      v_gene = c("TRAV1", "TRBV1", "TRAV2", "TRBV1", "TRBV1"),
      j_gene = c("TRAJ1", "TRBJ1", "TRAJ2", "TRBJ1", "TRBJ1"),
      full_nt_sequence = c("AAA", "BBB", "CCC", "BBB", "BBB")
    )
    # WHEN clones are counted with default parameters (up to two alphas per clone)
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData))
    # THEN they form one clone, flagged as unpaired on the alpha chain
    expect_equal(nrow(testClones), 1)
    expect_true(testClones$unpairedMultiChain)
    expect_equal(testClones$unpairedChainTypes, "TRA")

    # GIVEN the same data plus one lib that carries both alphas
    testData = rbind(testData, data.frame(
      libid = c("lib003", "lib003", "lib003"),
      v_gene = c("TRAV1", "TRAV2", "TRBV1"),
      j_gene = c("TRAJ1", "TRAJ2", "TRBJ1"),
      full_nt_sequence = c("AAA", "CCC", "BBB")
    ))
    # WHEN clones are counted
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData))
    # THEN the dual-alpha clone is not flagged
    expect_equal(nrow(testClones), 1)
    expect_false(testClones$unpairedMultiChain)
    expect_equal(testClones$unpairedChainTypes, "")

    # GIVEN three alphas on one beta, where only two of them occur together in a lib
    testData = data.frame(
      libid = c("lib000", "lib000", "lib000", "lib001", "lib001"),
      v_gene = c("TRAV1", "TRAV2", "TRBV1", "TRAV3", "TRBV1"),
      j_gene = c("TRAJ1", "TRAJ2", "TRBJ1", "TRAJ3", "TRBJ1"),
      full_nt_sequence = c("AAA", "CCC", "BBB", "DDD", "BBB")
    )
    # WHEN clones are counted allowing three alphas, in full format
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData), maxA = 3, format = "full")
    # THEN every row of the clone carries the flag
    expect_equal(length(unique(testClones$cloneId)), 1)
    expect_true(all(testClones$unpairedMultiChain))
    expect_true(all(testClones$unpairedChainTypes == "TRA"))

    # GIVEN a clone with one alpha and one beta
    testData = data.frame(
      libid = c("lib000", "lib000", "lib001", "lib001"),
      v_gene = c("TRAV1", "TRBV1", "TRAV1", "TRBV1"),
      j_gene = c("TRAJ1", "TRBJ1", "TRAJ1", "TRBJ1"),
      full_nt_sequence = c("AAA", "BBB", "AAA", "BBB")
    )
    # THEN it is not flagged
    testClones = getClonesFromTcrGraph(makeTcrGraph(testData))
    expect_false(testClones$unpairedMultiChain)
  })
  }
)
