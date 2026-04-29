# Parsing (m)cool files

These functions are the workhorse internal functions used to import a
`.(m)cool` file as GInteractions (wrapped into a `HiCExperiment` object
by [`HiCExperiment()`](HiCExperiment.md) function).

## Usage

``` r
.getCoolAnchors(file, resolution = NULL, balanced = "cooler")

.getCountsFromPair(file, pair, anchors, resolution = NULL)

.getCounts(file, coords, anchors, resolution = NULL)

.fetchCool(file, path, resolution = NULL, idx = NULL, ...)

.dumpCool(file, resolution = NULL)

.lsCoolFiles(file, verbose = FALSE)

.lsCoolResolutions(file, verbose = FALSE)

.cool2seqinfo(file, resolution = NULL)

.cool2gi(file, coords = NULL, resolution = NULL)
```

## Arguments

- file:

  path to a Hi-C contact file (in (m)cool format)

- resolution:

  resolution of the contact matrix

- balanced:

  import balancing scores

- pair:

  Genomic coordinates to extract contacts for, stored as a Pairs of
  GRanges (e.g. S4Vectors::Pairs(GRanges("II:200000-300000"),
  GRanges("II:70000-100000"))).

- anchors:

  anchors from .getCoolAnchors()

- coords:

  Genomic coordinates to extract contacts for, stored as a GRanges
  object

- path:

  Internal path of the cool file to check

- idx:

  Index to extract from the cool (HDF5) file

- ...:

  Other arguments passed to .fetchCool

- verbose:

  Print resolutions in the console

## Value

Silently, a numerical vector of resolutions stored in the cool file
