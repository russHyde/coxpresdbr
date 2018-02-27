----

TODO list for package `coxpresdbr`

----

[comment]: <> "Use 4-spaces to indent sublists"
[comment]: <> "Use '-' as the bullet-point icon"
[comment]: <> "Use 4-dashes for horizontal break"

----

# R

## `coxpresdbr_io.R`

- adds function `import_coexpresdb(gene_id, db_archive)`

    - returns details for all partners of a given gene `gene_id` that are
      present in the tar.bz2 file `db_archive`

- adds function `get_coexpresdb_universe(db_archive, gene_universe)`

    - returns the ids for all genes that are present in the archive
    `db_archive`

    - If `gene_universe` is specified, returns the intersection of this
      with those genes that have `coxpresdb` information

## `coxpresdbr_parse.R`

- adds function `get_coex_partners(gene_id, k=100, mrank_threshold=0,
  db_archive, gene_universe)`

    - returns the k-most coexpressed genes wrt gene `gene_id` with mutual rank
    greater than mrank_threshold

    - extracts these genes from `db_archive`, a tar.bz2 as downloaded from
    coxpresdb

## `coxpresdbr_stats.R`

- What would user want to compare against?

    - Assume the user has performed a separate experiment for which they either
      have p-values for each gene, or correlation values for each gene against
      a comparator `gene_id`

    - meta-analysis of p-values / z-scores for the partners of `gene_id`

    - average correlation of the partners of `gene_id` with `gene_id`

    - want some idea of statistical significance, by inverting the averaged
      z-scores, or by randomly sampling across the gene-universe

----

# tests

----
