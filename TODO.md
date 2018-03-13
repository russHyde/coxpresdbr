----

# TODO list for package `coxpresdbr`

<!-- "Use 4-spaces to indent sublists" -->
<!-- "Use '-' as the bullet-point icon" -->
<!-- "Use 4-dashes for horizontal break" -->

----

# README.md

- explain how to pull out the partners of a single gene

- writes up an example workflow for use of coxpresdbr

- adds a disclaimer re time-to-lookup and suggestion re reducing archive size
  before running an analysis

----

# R/

## `coxpresdbr_classes.R`

- `CoxpresDbPartners`

    - adds this class

- `CoxpresDbImporter`

   - moves this class from `coxpresdbr_io.R` into `coxpresdbr_classes.R`

## `coxpresdbr_workflows.R`

- `run_coex_partner_workflow`

    - runs `get_coex_partners` then `evaluate_coex_partners` for each gene in
      `gene_ids`

    - returns partners for each gene in `gene_ids` (all data returned by
      `get_coex_partners`) and returns z-statistics over the set of partners
      for each gene

    - results should be a named class: CoxpresDbPartners(@partners,
      @proximity_stats)

## `coxpresdbr_io.R`

- `import_all_coex_partners(gene_id, importer)`

    - uses more than one gene in `gene_id` argument

## `coxpresdbr_parse.R`

- Warn the user if any of the source genes is absent from the coexpression
  database

- `get_coex_partners`

    - rewrite to return a CoxpresDbPartners object

## `coxpresdbr_stats.R`

- `cluster_by_coex_partnership(coex_partners, source_nodes_only = TRUE,
  drop_disparities = TRUE)`

    - runs on results from `run_coex_partner_workflow`

    - links pairs (a, b) of genes from the `gene_ids` if the `a` is in the
      `coex_partners` of `b` (or vice versa) returns clusters by connectivity

    - appends to a CoxpresDbPartners object:

        - @cluster_graph: tidygraph::tbl_graph for connectivity between genes

        - @partners: data-frame of `source_id` to `target_id` where both
          `source_id` and `target_id` are in `gene_ids` (subset of data
          returned by get_coex_partners)

- `evaluate_coex_partners(x, coex_partners, ...)`

    - rewrites to take a CoxpresDbPartners object and return an appended
      CoxpresDbPartners object.

    - rewrites to take `coex_partners` as first argument, so that I can pipe
      from get_coex_partners() %>% evaluate_coex_partners() %>%
      cluster_by_coex_partnership()

    - [note x is a dataframe containing p-values / gene-ids / directions and
      must come first so that I can add a generic function for DGELRT and
      MArrayLM objects]

    - sets `evaluate_coex_partners` to a generic function over DGELRT /
      MArrayLM etc

    - adds an alternative-analysis-method switch, eg,
      `evaluate_coex_partners(x, coex_partners, method = c("sumz",
      "enrichment"), ...)`

- `evaluate_coex_partners(x: DGELRT/DGEExact, coex_partners)`

- `evaluate_coex_partners(x: MArrayLM, coex_partners)`

- What would user want to compare against?

    - Assume the user has performed a separate experiment for which they either
      have p-values for each gene, or correlation values for each gene against
      a comparator `gene_id`

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

        - if multiple copies of a gene are present, take the average of it's
          z-scores before running analysis

    - data-types that can be passed to metap::sumz (p is a vector of one-tailed
      p-values)

    - data-frame: `gene_id`, `p_value` (two-tailed), `direction`

        - the above datasets should be converted into this form

----

# tests/

## `test_coxpresdbr_parse.R`

- tests both applying filters and combining several source-genes datasets
  together

## `test_coxpresdbr_stats.R`

- tests that the same (two-sided) p-value returns on swapping the `direction`
  of input data-changes

----

# TIMELINE

----

# 2018-03-12

## `coxpresdbr_io.R`

- `import_coex_db(gene_id, db_archive)`

    - removes this function since it has been replaced by the method
      `import_all_coex_partners(gene_id, importer)`

## `coxpresdbr_parse.R`

- `get_coex_partners` | `.filter_coex_partners`

    - rewrites these functions to work with a CoxpresDbImporter object rather
      than a db_archive

    - replaces the call to `import_coex_db` with one to
      `import_all_coex_partners`

# 2018-03-09

## `coxpresdbr_io.R`

- rewrites the IO section to create a CoxpresdbImporter object

    - automatically extracts the gene-universe and filepaths from the coxpresdb
      archive

    - allows exploding a compressed archive into a temp directory

    - allows import from an exploded archive

    - rewrites all import functions to use the CoxpresdbImporter rather than a
      `db_archive` parameter

# 2018-03-06

## `coxpresdbr_stats.R`

- `evaluate_coex_partners(x, coex_partners, ...)`

    - with x a dataframe containing p-values/gene-ids/directions

    - adds function that takes a data-frame with columns `gene_id`, `p_value`
      and `direction` and a coexpression database (as a data-frame)

    - adds call to metap::twoToOne to convert two-sided tests to one-sided
      tests

    - joins `pval_df` and `coex_partners` on `gene_id` and `target_id`

    - for each source-gene in `coex_partners`:

        - computes average z-score (symmetrically wrt two-tailed test to fix
          numerical inaccuracies in metap implementation) using metap::sumz
          across the target-genes

        - counts the number of target-genes

        - converts z-score to p-value

- Note that there are some numerical inaccuracies in `metap` implementation of
  `sumz`

    - this is due to more accurate storage of p-values near 0 than near to 1

    - hence we took average of the z-scores:, ie, sumz(p)$z and -sumz(1-p)$z
      and return two-tailed p-values

# 2018-03-05

## `coxpresdbr_parse.R`

- adds function `get_coex_partners(gene_ids, db_archive, gene_universe = NULL,
  n_partners=100, mr_threshold = 0, cor_threshold = -1)`

    - returns the k-most coexpressed genes wrt each gene in `gene_ids` with
    mutual rank at-most `mrank_threshold`, correlation at-least `cor_threshold`
    and with each returned gene within the set `gene_universe`

    - extracts these genes from `db_archive`, using `import_coex_db`. So
      `db_archive` is a tar.bz2 as downloaded from coxpresdb.jp

- adds function `.filter_coex_partners` that takes a coex dataframe as input
  and applies all requested filters

    - Currently this function only works for the partnerset of a single gene
      (to save on manipulating large dataframes)

## `test_coxpresdbr_parse.R`

- tests subsetting a coex dataframe to just the top k genepartners for a given
  gene

- tests subsetting a coex dataframe to keep genepairs that have mutual rank or
  correlation above some threshold

- tests combining the coex dataframes for multiple genes together

## `coxpresdbr_io.R`

- adds function `import_coex_db(gene_id, db_archive)`

    - returns details for all partners of a given gene `gene_id` that are
      present in the tar.bz2 file `db_archive`

    - adds code to extract a file from a `*.tar.bz2` using
    `data.table::fread()`; returning a tibble

- adds function `import_coex_db_universe(db_archive)`

    - returns the ids for all genes that are present in the archive
    `db_archive`

    - Considered including a `gene_subset` so that the coxpresdb-annotated
    universe could be automatically restricted to a subset of the user-defined
    gene universe. Decided against this as it would give the function too many
    responsibilities

## `test_coxpresdbr_io.R`

- Adds a subset of the yeast coxpresdb dataset into tests/testthat for use in
  io testing
    - This was generated from the fission yeast dataset
      `Spo.v14-08.G4881-S224.rma.mrgeo.d.tar.bz2`

    - The original data was manipulated by selecting 10 genes and keeping the
      files for just those 10 genes; and with each file being restricted to
      just the selected 10 genes:
        - `tar -jtf Spo.v14*.tar.bz2 | grep -v "/$" | head > selected.txt`
        - `tar -xjf Spo.v14*.tar.bz2 $(cat selected.txt)`
        - `REGEX=$(perl -e 'chomp(@F=<>); $f=join "\\|", map {/.*\/(.*)$/; $1}
          @F; print($f)' < selected.txt)`
        - `grep -e "${REGEX}" <each_file>; <and replace the file>`

    - The data is stored in the `./tests/testthat` directory, rather than in
    `inst/extdata` or `data`

