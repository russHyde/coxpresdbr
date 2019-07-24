# coxpresdbr 0.0.0.9000

- Correlation coefficients, if present in the original coxpresdb archive, are
  no longer reported in the coexpression partners since these coefficients are
  not consistently reported across the different releases of coxpresdb (and so
  the user can no longer filter based on these values in `get_coex_partners`)

- User can use `.tar.bz2`, `.tar` or `.zip` archives as the source of CoxpresDB
  data and the files can be of two-column (`target_id, mutual_rank`) or
  three-column (`target_id, mutual_rank, correlation_coefficient`) format

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
