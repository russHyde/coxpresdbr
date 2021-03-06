###############################################################################
# Tests for statistics functions in `coxpresdbr` package
#
###############################################################################

context("Tests for geneset-combination statistics in `coxpresdbr` package")

###############################################################################

test_that("evaluate_coex_partners-data.frame: invalid input", {
  pval_empty <- tibble::tibble(
    gene_id = character(0),
    p_value = numeric(0),
    direction = integer(0)
  )

  coex_empty <- tibble::tibble(
    source_id = character(0),
    target_id = character(0)
  )

  # `x` should inherit from data.frame
  # ... (or on generalising be a DGELRT/MArrayLM object)
  expect_error(
    object = evaluate_coex_partners(
      x = "NOT A DATA FRAME",
      coex_partners = coex_empty
    ),
    info = "`x` should be a data-frame in evaluate_coex_partners"
  )

  expect_error(
    object = evaluate_coex_partners(
      x = tibble::tibble(
        GENE_ID = character(0),
        P_VALUE = numeric(0),
        DIRECTION = integer(0)
      )
    ),
    info = "`x` should have columns `gene_id`, `p_value` and `direction`"
  )

  # `coex_partners` should inherit from data.frame
  expect_error(
    object = evaluate_coex_partners(
      x = pval_empty,
      coex_partners = "Not a data-frame"
    ),
    info = "`coex_partners` should be a data-frame"
  )

  expect_error(
    object = evaluate_coex_partners(
      x = pval_empty,
      coex_partners = tibble::tibble(
        SOURCE_ID = character(0),
        TARGET_ID = character(0)
      )
    ),
    info = "`coex_partners` should have columns `source_id` and `target_id`"
  )
})

###############################################################################

test_that("evaluate_coex_partners-data.frame: valid input, few partners", {
  # source gene has no partners: no-row data-frame returns
  pval <- data.frame(
    gene_id = letters[1:3],
    p_value = c(1e-12, 0.001, 0.9),
    direction = c(1L, -1L, 1L),
    stringsAsFactors = FALSE
  )
  coex_no_partners <- tibble::tibble(
    source_id = "d",
    target_id = "e"
  )
  expect <- tibble::tibble(
    gene_id = character(0),
    n_partners = integer(0),
    z_score = numeric(0),
    p_value = numeric(0)
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_no_partners),
    expected = expect,
    info = "source gene has no partners => source gene is absent from output"
  )

  # source gene has a single partner: the returned p-value should match the
  # input p-value
  coex_one_partner <- tibble::tibble(
    source_id = "d",
    target_id = "a"
  )
  expect <- tibble::tibble(
    gene_id = "d",
    n_partners = 1L,
    z_score = 0.5 * (
      qnorm(5e-13, lower.tail = FALSE) -
        qnorm(1 - 5e-13, lower.tail = FALSE)
    ),
    p_value = 1e-12
  )
  expect_equal(
    object = as.data.frame(
      evaluate_coex_partners(pval, coex_one_partner)
    )[, 1:2],
    expected = as.data.frame(expect)[, 1:2],
    info = paste(
      "source gene has a single partner => gene-id/#partners for the result",
      "should match that of the target gene (non-numeric cols)"
    )
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_one_partner)$z_score,
    expected = expect$z_score,
    info = paste(
      "source gene has a single partner => z-score should match that of the",
      "target gene (z-score col)"
    ),
    tolerance = 1e-4,
    scale = expect$z_score
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_one_partner)$p_value,
    expected = expect$p_value,
    info = paste(
      "source gene has a single partner => p-value should match that of the",
      "target gene (p_value col)"
    ),
    tolerance = 1e-4,
    scale = expect$p_value
  )

  # p-value returned should be independent of `direction` and
  # z-scores returned should be reversed if `direction` is switched
  expect_equal(
    object = pval %>%
      dplyr::mutate(direction = -1 * direction) %>%
      evaluate_coex_partners(coex_one_partner) %>%
      as.data.frame(),

    expected = pval %>%
      evaluate_coex_partners(coex_one_partner) %>%
      dplyr::mutate(z_score = -1 * z_score) %>%
      as.data.frame(),

    info = paste(
      "if all directions change, the p-value should remain the same and the",
      "z-score should reverse"
    )
  )
})

test_that("evaluate_coex_partners: multiple partners", {
  # Highly significant input genes that all have the same direction
  # threw a bug (2019-07; all with negative sign) where the output p-value
  # was zero (this resulted from passing rows with p_value = 0 into
  # metap::sumz; something which happened because of numerical rounding while
  # converting two-tailed p-values to z-scores via one-tailed p-values within
  # metap).

  # this is a two-tailed p-value
  common_p_value <- 1e-17

  # the z-score corresponding to the equivalent one-tailed p-value is
  common_z_score <- qnorm(common_p_value / 2)

  # If all neighbours of a given node have this common z-score
  # then the neighbourhood z-score will be sqrt(# neighbours) * z
  # ... because we use the sum-of-z-scores method to summarise p-values, and
  # this has the form sum(z_i : i = 1..n) / sqrt(n)

  p_vals <- data.frame(
    gene_id = letters,
    p_value = common_p_value,
    direction = -1L,
    stringsAsFactors = FALSE
  )

  coex_partners <- data.frame(
    source_id = "a",
    target_id = letters[-1],
    stringsAsFactors = FALSE
  )

  expect_equal(
    object = evaluate_coex_partners(p_vals, coex_partners)$z_score,
    expected = common_z_score * sqrt(25),
    info = paste(
      "when all targets have the same z_score and direction (all -ve), the",
      "neighbourhood z-score should be common_z_score * sqrt(# targets)"
    )
  )

  expect_equal(
    object = evaluate_coex_partners(
      dplyr::mutate(p_vals, direction = 1L),
      coex_partners
    )$z_score,
    expected = -1 * common_z_score * sqrt(25),
    info = paste(
      "when all targets have the same z_score and direction (all +ve), the",
      "neighbourhood z_score should be common_z_score * sqrt(# targets)"
    )
  )

  expect_equal(
    object = evaluate_coex_partners(
      p_vals, coex_partners
    )$p_value,
    expected = evaluate_coex_partners(
      dplyr::mutate(p_vals, direction = 1L), coex_partners
    )$p_value,
    info = paste(
      "p-values should be two-tailed; hence should be the same upon swapping",
      "directions"
    )
  )
})

###############################################################################

test_that(".format_coex_edges_for_tidygraph: invalid input", {
  # - Input should be a CoxpresDbPartners object
  # - Input should have a non-empty `partners` data-frame
  expect_error(
    object = .format_coex_edges_for_tidygraph(
      coex_partners = "NOT A CoxpresDbPartners object",
      cluster_source_nodes_only = TRUE
    ),
    info = "`coex_partners` should be a `CoxpresDbPartners` object"
  )

  expect_error(
    object = .format_coex_edges_for_tidygraph(
      coex_partners = new("CoxpresDbPartners"),
      cluster_source_nodes_only = TRUE
    ),
    info = "`coex_partners` should have a valid/non-empty `partners` entry"
  )
})

###############################################################################

test_that(".format_coex_edges_for_tidygraph: valid input", {
  # output should have columns `from`, `to`
  test_partners <- tibble::tibble(
    source_id = c("a", "c"),
    target_id = c("c", "b")
  )

  expect_equal(
    object = .format_coex_edges_for_tidygraph(
      coex_partners = new("CoxpresDbPartners", partners = test_partners),
      cluster_source_nodes_only = TRUE
    ),
    expected = tibble::tibble(
      from = "a",
      to = "c"
    ),
    info = paste(
      "if cluster_source_nodes_only only keep rows with a source_id in the",
      "target_id"
    )
  )

  expect_equal(
    object = .format_coex_edges_for_tidygraph(
      coex_partners = new("CoxpresDbPartners", partners = test_partners),
      cluster_source_nodes_only = FALSE
    ),
    expected = tibble::tibble(
      from = c("a", "c"),
      to = c("c", "b")
    ),
    info = paste(
      "if cluster_source_nodes_only is FALSE, keep all rows"
    )
  )
})

###############################################################################

test_that(".format_unsorted_nodes_for_tidygraph: invalid input", {
  # input should be a CoxpresDbPartners with a non-empty `gene_statistics`
  # output should have columns `name`, `z_score`, `p_value`, `direction`

  expect_error(
    object = .format_unsorted_nodes_for_tidygraph(
      coex_partners = "NOT A CoxpresDbPartners object"
    ),
    info = "`coex_partners` should be a `CoxpresDbPartners` object"
  )

  expect_error(
    object = .format_unsorted_nodes_for_tidygraph(
      coex_partners = new("CoxpresDbPartners")
    ),
    info = "`coex_partners` have a non-empty `gene_statistics` entry"
  )
})

###############################################################################

test_that(".format_unsorted_nodes_for_tidygraph: valid input", {
  z <- qnorm(0.25, lower.tail = FALSE)

  test_statistics <- tibble::tibble(
    gene_id = "a", p_value = 0.5, direction = 1, z_score = z
  )
  test_statistics_no_z <- test_statistics[, 1:3]

  expect_equal(
    object = .format_unsorted_nodes_for_tidygraph(
      coex_partners = new(
        "CoxpresDbPartners",
        gene_statistics = test_statistics
      )
    ),
    expected = tibble::tibble(
      name = "a", z_score = z, p_value = 0.5, direction = 1
    ),
    info = paste(
      ".format_unsorted_nodes... replaces gene_id with name and reorders",
      "columns"
    )
  )

  expect_equal(
    object = .format_unsorted_nodes_for_tidygraph(
      coex_partners = new(
        "CoxpresDbPartners",
        gene_statistics = test_statistics_no_z
      )
    ),
    expected = tibble::tibble(
      name = "a", z_score = z, p_value = 0.5, direction = 1
    ),
    info = paste(
      ".format_unsorted_nodes... adds z_score if it's missing"
    )
  )
})

###############################################################################

test_that(".add_direction_parities_to_coex_edges: invalid input", {
  # coex_partners should have a non-empty partners and gene_statistics
  #
  expect_error(
    object = .add_direction_parities_to_coex_edges(
      coex_partners = "NOT A CoxpresDbPartners object"
    ),
    info = "`coex_partners` should be a `CoxpresDbPartners` object"
  )

  expect_error(
    object = .add_direction_parities_to_coex_edges(
      coex_partners = new(
        "CoxpresDbPartners",
        partners = tibble::tibble(
          source_id = "a",
          target_id = "b"
        )
      )
    ),
    info = "`coex_partners` have a non-empty `gene_statistics` entry"
  )
  expect_error(
    object = .add_direction_parities_to_coex_edges(
      coex_partners = new(
        "CoxpresDbPartners",
        gene_statistics = tibble::tibble(
          gene_id = "abc",
          p_value = 0.2,
          direction = 1
        )
      )
    ),
    info = "`coex_partners` have a non-empty `partners` entry"
  )
})

###############################################################################

test_that(".add_direction_parities_to_coex_edges: valid input", {
  test_coex_partners <- new(
    "CoxpresDbPartners",
    gene_statistics = tibble::tibble(
      gene_id = letters[1:3],
      p_value = c(0.2, 0.5, 0.9),
      direction = c(1, -1, 1)
    ),
    partners = tibble::tibble(
      source_id = c("a", "b"),
      target_id = c("c", "a"),
      some_other_col = c(1, 2)
    )
  )

  expect_equal(
    object = .add_direction_parities_to_coex_edges(
      test_coex_partners
    )@partners,
    expected = tibble::tibble(
      source_id = c("a", "b"),
      target_id = c("c", "a"),
      some_other_col = c(1, 2),
      direction_parity = c(TRUE, FALSE)
    ),
    info = "add_direction_parities_to_coex_edges"
  )
})

###############################################################################

test_that("cluster_by_coex_partnership: invalid input", {
  # - Input should be a CoxpresDbPartners object
  # - The CoxpresDbPartners object should have a valid `gene_statistics`
  # - The CoxpresDbPartners object should have a valid `partners`
  # - `drop_disparities` should be Boolean

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = "NOT A CoxpresDbPartners object",
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should be a `CoxpresDbPartners` object"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = new(
        "CoxpresDbPartners",
        partners = tibble::tibble(source_id = "a", target_id = "b")
      ),
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should have a valid `gene_statistics` entry"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = new(
        "CoxpresDbPartners",
        gene_statistics = tibble::tibble(
          gene_id = "abc", p_value = "0.2", direction = "1"
        )
      ),
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should have a valid/non-empty `partners` entry"
  )

  # This is a valid CoxpresDbPartners object. But there should be no edges in
  # cluster graph for the dataset (when drop_disparities = TRUE), since the
  # direction of change for the two nodes differs. cluster_by_coex_partnership
  # should be able to deal with a graph with no edges
  test_coex_partners <- new(
    "CoxpresDbPartners",
    gene_statistics = tibble::tibble(
      gene_id = letters[1:2], p_value = c(0.1, 0.2), direction = c(-1, 1),
      z = c(-0.5, 1)
    ),
    partners = tibble::tibble(
      source_id = "a",
      target_id = "b"
    )
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = test_coex_partners,
      drop_disparities = "NOT A BOOLEAN"
    ),
    info = "`drop_disparities` should be Boolean"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = test_coex_partners,
      drop_disparities = logical(0)
    ),
    info = "`drop_disparities` should be a Boolean of length 1"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = test_coex_partners,
      drop_disparities = c(TRUE, FALSE)
    ),
    info = "`drop_disparities` should be a Boolean of length 1"
  )
})

###############################################################################

# NOTE: problems with comparing two seemingly identical graphs
# expect_equal(result_graph, expected_graph) fails with uninformative message
# : re differences in component 9: component 1 ...
# : the stated components can't be accessed however, so there is no way to
#   debug the test

expect_equal_graph <- function(object, expected, info) {
  expect_equal(
    object = igraph::vertex.attributes(object),
    expected = igraph::vertex.attributes(expected),
    info = paste(info, ": vertex-equality")
  )

  expect_equal(
    object = igraph::as_data_frame(object),
    expected = igraph::as_data_frame(expected),
    info = paste(info, ": edge-equality")
  )
}

###############################################################################

test_that("cluster_by_coex_partnership: valid input", {
  # TODO: test to use drop_disparities = FALSE
  z_05 <- qnorm(0.5 / 2, lower.tail = FALSE)
  z_02 <- qnorm(0.2 / 2, lower.tail = FALSE)

  # compare cluster graph
  test_statistics_no_disparity <- tibble::tibble(
    gene_id = c("a", "b"),
    p_value = c(0.5, 0.2),
    direction = c(1, 1),
    z_score = c(z_05, z_02)
  )

  test_statistics_with_disparity <- tibble::tibble(
    gene_id = c("a", "b"),
    p_value = c(0.5, 0.2),
    direction = c(-1, 1),
    z_score = c(-1 * z_05, z_02)
  )

  test_partners <- tibble::tibble(
    source_id = c("a", "b"),
    target_id = c("b", "a")
  )

  empty_graph <- tidygraph::as_tbl_graph(
    list(
      nodes = tibble::tibble(name = character(0)),
      edges = tibble::tibble(from = numeric(0), to = numeric(0))
    )
  )

  # --- #
  nodes_no_disparity <- test_statistics_no_disparity %>%
    dplyr::rename(name = .data[["gene_id"]]) %>%
    magrittr::extract(c("name", "z_score", "p_value", "direction"))

  edges_no_disparity <- tibble::tibble(
    from = c(1, 2),
    to = c(2, 1),
    direction_parity = c(TRUE, TRUE)
  )
  expected_graph_no_disparity <- tidygraph::as_tbl_graph(
    list(nodes = nodes_no_disparity, edges = edges_no_disparity)
  )
  result_graph_no_disparity <- cluster_by_coex_partnership(
    new(
      "CoxpresDbPartners",
      gene_statistics = test_statistics_no_disparity,
      partners = test_partners
    )
  )@cluster_graph

  expect_equal_graph(
    object = result_graph_no_disparity,
    expected = expected_graph_no_disparity,
    info = "Extraction of the partnership graph for a pair of nodes"
  )

  # --- #
  # With a single disparity in the expression data
  result_graph_with_disparity <- cluster_by_coex_partnership(
    new(
      "CoxpresDbPartners",
      gene_statistics = test_statistics_with_disparity,
      partners = test_partners
    )
  )@cluster_graph

  expect_equal_graph(
    object = result_graph_with_disparity,
    expected = empty_graph,
    info = "The cluster graph should be empty for a genepair with expression
disparity"
  )

  # --- #
  nodes_ignoring_disparity <- test_statistics_with_disparity %>%
    dplyr::rename(name = .data[["gene_id"]]) %>%
    magrittr::extract(c("name", "z_score", "p_value", "direction"))

  edges_ignoring_disparity <- tibble::tibble(
    from = c(1, 2),
    to = c(2, 1),
    direction_parity = c(FALSE, FALSE)
  )
  expected_graph_ignoring_disparity <- tidygraph::as_tbl_graph(
    list(nodes = nodes_ignoring_disparity, edges = edges_ignoring_disparity)
  )
  result_graph_ignoring_disparity <- cluster_by_coex_partnership(
    new(
      "CoxpresDbPartners",
      gene_statistics = test_statistics_with_disparity,
      partners = test_partners
    ),
    drop_disparities = FALSE
  )@cluster_graph

  expect_equal_graph(
    object = result_graph_ignoring_disparity,
    expected = expected_graph_ignoring_disparity,
    info = "Make a cluster graph but ignore disparity of gene expression"
  )

  # Results should be the same whether z-scores are passed in with the
  # coex_partners object, or whether they are computed wihin
  # cluster_by_coex_partnership
  expect_equal_graph(
    object = cluster_by_coex_partnership(
      new(
        "CoxpresDbPartners",
        gene_statistics = test_statistics_no_disparity,
        partners = test_partners
      )
    )@cluster_graph,
    expected = cluster_by_coex_partnership(
      new(
        "CoxpresDbPartners",
        gene_statistics = test_statistics_no_disparity[
          ,
          c("gene_id", "p_value", "direction")
        ],
        partners = test_partners
      )
    )@cluster_graph,
    info = "Results of cluster_by_.. should be the same regardless of whether
    z-scores are computed within the cluster_by... workflow"
  )
})

###############################################################################
