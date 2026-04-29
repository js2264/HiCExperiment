# Utils functions

Utilities to facilitate parsing/handling of coordinates, GInteractions,
Pairs, ... These functions are not exported.

## Usage

``` r
splitCoords(coords)

coords2char(coords, big.mark = ",")

char2coords(char)

sortPairs(pairs)

asGInteractions(df)

sdiag(A, k = 0)

.df2symmmat(diag, score)

distanceDecay(dump, threshold = NULL)

detrendingModel(file, resolution)

.fixRegions(gis, bins, coords)
```

## Arguments

- coords:

  A set of genomic coordinates (either as a GRanges object or as a
  character string)

- big.mark:

  Separator for thousands when printing out genomic coordinates as
  character

- char:

  char (e.g. "II:30001-50000" or "II:30001-50000\|II:60001-80000")

- pairs:

  Pairs object

- df:

  a data.frame to turn into a GInteraction object.

- A:

  Numerical matrix

- k:

  secondary diagonal k

- diag:

  vector of distances to diagonal

- score:

  scores to parse into symmetrical matrix

- dump:

  dumped contacts as GInteractions, e.g. from .dumpCool

- threshold:

  maximum distance to compute distance decay for

- file:

  path to a HiC contact matrix file

- resolution:

  Resolution to use with the HiC contact matrix file

- gis:

  GInteractions object

- bins:

  Larger set of regions (usually bins from HiCExperiment)

## Value

Reformatted coordinates or GInteractions.
