----

# TODO list for package `coxpresdbr`

<!-- "Use 4-spaces to indent sublists" -->
<!-- "Use '-' as the bullet-point icon" -->
<!-- "Use 4-dashes for horizontal break" -->

----

# R/

## `coxpresdbr_io.R`

- `import_coex_db(gene_id, db_archive)`

    - uses more than one gene in `gene_id` argument

## `coxpresdbr_stats.R`

- What would user want to compare against?

    - Assume the user has performed a separate experiment for which they either
      have p-values for each gene, or correlation values for each gene against
      a comparator `gene_id`

    - Enrichment of gene-partners of a given gene in `significant` list

    - meta-analysis of p-values / z-scores for the partners of `gene_id`

    - average correlation of the partners of `gene_id` with `gene_id`

    - want some idea of statistical significance, by inverting the averaged
      z-scores, or by randomly sampling across the gene-universe

----

# tests/

## `test_coxpresdbr_parse.R`

- tests both applying filters and combining several source-genes datasets
  together

----

# TIMELINE

----

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

