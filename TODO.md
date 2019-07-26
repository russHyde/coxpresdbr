----

# TODO list for package `coxpresdbr`

<!-- "Use '-' as the bullet-point icon" -->
<!-- "Use 4 spaces to indent sublists" -->
<!-- "Use 4 dashes for horizontal break" -->
<!-- "Use 4 tildes to wrap code blocks" -->

----

# LONG-TERM

- Access coxpresdb via SPARQL

----

# README.md

- explain how to pull out the partners of a single gene

- writes up an example workflow for use of coxpresdbr

- adds a disclaimer re time-to-lookup and suggestion re reducing archive size
  before running an analysis

----

# R/

## `coxpresdbr_classes.R`

## `coxpresdbr_workflows.R`

- `run_coex_partner_workflow(
        gene_ids, gene_statistics, importer,
        gene_universe = NULL, n_partners = 100, ...)`

    - runs `cluster_by_coex_partnership` over a CoxpresDbPartners object and
      appends the resulting graph etc

## `coxpresdbr_io.R`

- `get_all_coex_partners(gene_id, importer)`

    - uses more than one gene in `gene_id` argument

## `coxpresdbr_parse.R`

- Warn the user if any of the source genes is absent from the coexpression
  database

- `get_coex_partners`

    - rewrite to return a CoxpresDbPartners object

## `coxpresdbr_stats.R`

- `.format_coex_edges_for_tidygraph`

    - split filtering source/target-ness into a separate function

- `cluster_by_coex_partnership`

    - add `cluster_source_nodes_only = BOOLEAN` to formals

    - append to a CoxpresDbPartners object:

        - @partners: data-frame of `source_id` to `target_id` where both
          `source_id` and `target_id` are in `gene_ids` (subset of data
          returned by get_coex_partners)

- `evaluate_coex_partners(x, coex_partners)`

    - rewrite to take a CoxpresDbPartners object and return an appended
      CoxpresDbPartners object.

    - rewrite to take `coex_partners` as first argument, so that I can pipe
      from get_coex_partners() %>% evaluate_coex_partners() %>%
      cluster_by_coex_partnership()

    - set `evaluate_coex_partners` to a generic function over DGELRT /
      MArrayLM etc

    - add an alternative-analysis-method switch, eg, `evaluate_coex_partners(x,
      coex_partners, method = c("sumz", "enrichment"))`

- `evaluate_coex_partners(x: DGELRT/DGEExact, coex_partners)`

- `evaluate_coex_partners(x: MArrayLM, coex_partners)`

----

- notes

    - What would user want to compare against?

        - Assume the user has performed a separate experiment for which they
          either have p-values for each gene, or correlation values for each
          gene against a comparator `gene_id`

        - Enrichment of gene-partners of a given gene in `significant` list

        - meta-analysis of p-values / z-scores for the partners of `gene_id`

        - average correlation of the partners of `gene_id` with `gene_id`

        - want some idea of statistical significance, by inverting the averaged
          z-scores, or by randomly sampling across the gene-universe

    - What data-types would the user have access to?

        - edgeR::DGELRT/DGEExact objects

            - would need to indicate gene_id column of `genes` component

        - limma::MArrayLM/TestResults objects

            - would need to indicate mappings between rows and gene_id

            - if multiple copies of a gene are present, take the average of
              it's z-scores before running analysis

        - data-types that can be passed to metap::sumz (p is a vector of
          one-tailed p-values)

        - data-frame: `gene_id`, `p_value` (two-tailed), `direction`

            - the above datasets should be converted into this form

----

# tests/

## `test_coxpresdbr_parse.R`

- tests both applying filters and combining several source-genes datasets
  together

## `test_coxpresdbr_stats.R`

- `cluster_by_coex_partnership`

    - tests for `drop_disparities = FALSE`
