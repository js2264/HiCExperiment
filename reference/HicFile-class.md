# `HicFile` S4 class

The `HicFile` class describes a `BiocFile` object, pointing to the
location of a .hic file (usually created with juicer) and containing 3
additional slots:

1.  resolution: at which resolution the associated .hic file should be
    parsed;

2.  pairsFile: the path (in plain character) to an optional pairs file
    (stored as a `PairsFile` object);

3.  metadata: a list metadata

HicFile methods.

## Arguments

- path:

  String; path to a .hic file

- resolution:

  numeric; resolution to use with mcool file

- pairsFile:

  String; path to a pairs file

- metadata:

  list.

- object:

  A `HicFile` object.

## See also

[`CoolFile()`](CoolFile-class.md), [`HicproFile()`](HicproFile-class.md)

## Examples

``` r
hicPath <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
pairsPath <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
hic <- HicFile(
  hicPath, 
  resolution = 16000, 
  pairsFile = pairsPath, 
  metadata = list(type = 'example')
)
hic
#> HicFile object
#> .hic file: /github/home/.cache/R/ExperimentHub/3e265ac673a5_7836 
#> resolution: 16000 
#> pairs file: /github/home/.cache/R/ExperimentHub/3e2648419160_7753 
#> metadata(1): type
resolution(hic)
#> [1] 16000
pairsFile(hic)
#>                                                  EH7703 
#> "/github/home/.cache/R/ExperimentHub/3e2648419160_7753" 
metadata(hic)
#> $type
#> [1] "example"
#> 
```
