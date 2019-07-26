
<!-- README.md is generated from README.Rmd. Please edit that file -->

-----

# coxpresdbr: a package for importing data from coxpresdb.jp and for comparing your own datasets against coxpresdb

-----

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/russHyde/coxpresdbr.svg?branch=master)](https://travis-ci.org/russHyde/coxpresdbr)

[![Codecov test
coverage](https://codecov.io/gh/russHyde/coxpresdbr/branch/master/graph/badge.svg)](https://codecov.io/gh/russHyde/coxpresdbr?branch=master)
<!-- badges: end -->

-----

[`coxpresdb.jp`](https://coxpresdb.jp/) is a website / database that
summarises the correlation between pairs of genes across thousands of
different gene-expression studies. As of 2019, it contains downloadable
coexpression data from humans and a range of model organisms (rodents,
nematodes, flies, yeast etc). In humans, for example, the dataset may
contain all correlations amongst \~ 20k genes, as derived from
summarising \~ 160k samples.

This package provides some data-access and statistical tools for working
with `coxpresdb.jp` in *R*.

Why?

Several tools exist for combining network analysis with gene expression
data or for finding some other higher-level context for changes observed
in gene expression data (*GSEA* / *GSVA*). These may construct a network
*de novo* (eg, *WGCNA*) or identify interesting subcomponents of an
existing network after superimposing some expression statistics upon
that network (eg, *jActiveModules* in *Cytoscape*).

To use *WGCNA* requires tons of expression data. And no-one I know has
that many good quality samples.

*jActiveModules* and related tools are great. But they don’t provide
much statistical context for the submodules that they select. Would they
find an *interesting* subnetwork for white-noise input? How many
subnetworks are considered before selecting the small number of optimal
modules. And when applied, how relevant is network construction to the
results obtained: if there are non-uniformities in data-acquisition (eg,
literature searches) or sampling across the network (eg, funding /
publication biases, sub-genome-scale studies), if technical issues
influence network structure, and if the network is constructed from a
data source that may be only loosely-coupled to gene expression (eg,
protein-protein interaction data) could a good network analysis tool
still return meaningless results?

As an alternative, to put gene expression data in some network context,
one could use coexpression databases. This provides a network structure
upon which gene expression data can be overlaid (so can still work when
samples are limited in your experiment), that is built from genome-scale
studies and where the approach is closely connected to gene-expression
studies and where we can perform both gene-focussed and network-wide
analyses. (admittedly the tissue types and cellular treatments used in
the datasets from which coexpression databases are constructed may still
be biased).

Here’s what I might do:

  - compute fold-changes and p-values from an RNA-Seq experiment

  - for each significant gene, have a see what genes it tends to be
    coexpressed with / correlated with

  - have a look if those coexpression partners are also differentially
    expressed … in the same direction … and have shared ChIP marks.

Let’s formalise that.

## Importing data

### Data formats

<!-- TODO: obtain data for a specific gene via web-search:
  eg, see
  https://coxpresdb.jp/api/v3/datasets/Hsa-u.c1-0/genes/ncbigene:1/interactions?top=10
-->

<!-- TODO: import data from a collapsed single-file all-gene database -->

`coxpresdbr` currently only works with downloaded copies of the
`coxpresdb.jp` databases. A given database comes as a compressed archive
and within that archive, a file is present for each gene. In recent
memory there have been two related file formats and two compression
types (`.tar.bz2` and, more recently, `.zip`).

The first file format had three columns (tab-separated)

(`target_gene, mutual_rank, correlation_coefficient`)

but has been replaced by a two column (`target_gene, mutual_rank`)
format. Column names are not present and the `source_gene` (the name of
the file) and `target_gene` (contents of the first column) are numeric
Entrez IDs.

When working with the coxpresdb files, I tend to make a smaller archive
setting up the file for each source-gene to contain only the top 100
most correlated genes.

### The `CoxpresDbAccessor` Class

We first define an importer object that permits access to the
gene-specific files stored in the coxpresdb archive.

    my_importer <- CoxpresDbAccessor(path_to_coxpresdb_archive)

### Importing the coexpression data & filtering

Then to pull out the top 10 coexpression partners for a given set of
genes, you can call:

    get_coex_partners(my_gene_set, importer = my_importer, n_partners = 10)

This returns a data-frame with columns: ( `source_id, target_id,
mutual_rank`; correlation coefs, if present in the coxpresdb archive,
are disregarded since they aren’t present in all coxpresdb releases)

You can specify the universe of genes from which the top-N coexpression
partners are chosen by using the argument `universe=...` and you can
filter to ensure that the returned gene-pairs have a mutual-rank score
below some threshold using `mr_threshold`.

<!-- TODO: Annotating an imported dataset -->

## Gene-focussed statistical summaries

<!-- `metap` as the basis of the summaries -->

<!--
  Single comparison workflow:
  - overlaying differential-expression results on a predefined network
  - summarising over nodes
-->

<!--
  Correlation comparison workflow
  - overlaying correlation data on a predefined network
  - summarising over edges
-->

## Network-level analysis workflows

<!--
  Using coxpresdbr-derived networks as input to other network analysis tools
-->
