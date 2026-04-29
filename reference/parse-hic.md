# Parsing hic files

These functions are the workhorse internal functions used to import a
`.hic` file as GInteractions (wrapped into a `HiCExperiment` object by
[`HiCExperiment()`](HiCExperiment.md) function).

## Usage

``` r
.hic2gi(file, coords = NULL, resolution = NULL)

.lsHicResolutions(file, verbose = FALSE)

.getHicAnchors(file, resolution = NULL)

.hic2seqinfo(file)

.dumpHic(file, resolution = NULL)
```

## Arguments

- file:

  path to a Hi-C contact file in .hic format

- coords:

  NULL, character, or GRanges. Can also be a Pairs object of paired
  GRanges (length of 1).

- resolution:

  resolution of the contact matrix to use

- verbose:

  Print resolutions in the console

## Value

a GInteractions object

vector
