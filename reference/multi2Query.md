# Querying multiple slices of a contact matrix

These functions are the workhorse internal functions used to extract
counts from multiple genomic coordinates in a Hi-C contact matrix.

## Usage

``` r
.multi2DQuery(
  file,
  resolution,
  pairs,
  maxDistance = NULL,
  bed = NULL,
  BPPARAM = BiocParallel::bpparam()
)
```

## Arguments

- file:

  path to a Hi-C contact file (can be any format, (m)cool, .hic, or
  HiC-Pro-derived)

- resolution:

  resolution to use to import matrix over specified targets

- pairs:

  slices to read, provided as a Pairs object

- maxDistance:

  Maximum distance to use when compiling distance decay

- bed:

  associated bed file for HiC-Pro derived contact matrix.

- BPPARAM:

  BiocParallel parameters

## Value

a GInteractions object with `count`, `balanced`, `detrended` and
`expected` scores
