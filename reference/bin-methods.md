# HiCExperiment binning methods

HiCExperiment binning methods

## Usage

``` r
# S4 method for class 'GInteractions,numeric'
bin(x, resolution, seqinfo = NULL)

# S4 method for class 'PairsFile,numeric'
bin(x, resolution, seqinfo = NULL)
```

## Arguments

- x:

  A PairsFile or GInteractions object

- resolution:

  Which resolution to use to bin the interactions

- seqinfo:

  Seqinfo object

## Examples

``` r
pairsf <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
pf <- PairsFile(pairsf)
```
