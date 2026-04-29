# Pairs parsing functions

Pairs parsing functions

## Usage

``` r
.pairs2gi(
  file,
  chr1.field = NULL,
  start1.field = NULL,
  chr2.field = NULL,
  start2.field = NULL,
  strand1.field = NULL,
  strand2.field = NULL,
  frag1.field = NULL,
  frag2.field = NULL,
  nThread = 1,
  nrows = Inf
)
```

## Arguments

- file:

  pairs file. Default formatting is
  `<readname>\t<chr1>\t<start1>\t<chr2>\t<start2>`.

- chr1.field, start1.field, chr2.field, start2.field, strand1.field,
  strand2.field, frag1.field, frag2.field:

  Index of the column in which each field is contained in the pairs
  file.

- nThread:

  Number of CPUs to use to import the `pairs` file in R

- nrows:

  Number of pairs to import

## Value

a GInteractions object
