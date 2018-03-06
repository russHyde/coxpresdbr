# coxpresdbr 0.0.0.9000

- Added a `NEWS.md` file to track changes to the package.

- Added functions for importing CoxpresDB datasets from an existing archive.
  The functions work without having to decompress/recompress the archive.
  The user must have `tar`, `grep` and `bzip2` installed on a unix-alike
  system.

- Added a function for selecting subsets of the coexpression partners of a
  given gene (or genes): `get_coex_partners`. The user can select to keep the
  top `k` coexpression partners of a given gene, or any partners with a mutual
  rank at-most-equal to some threshold, or any partners with a correlation
  coefficient at-least-equal to some threshold

