# Example datasets provided in `HiCExperiment` & `HiContactsData`

Example datasets provided in `HiCExperiment` & `HiContactsData`

## Usage

``` r
data(centros_yeast)

contacts_yeast(full = FALSE)

contacts_yeast_eco1(full = FALSE)
```

## Format

An object of class `"GRanges"`.

## Source

HiContacts

## Arguments

- full:

  Whether to import all interactions

## Examples

``` r
data(centros_yeast)
centros_yeast
#> GRanges object with 16 ranges and 0 metadata columns:
#>        seqnames        ranges strand
#>           <Rle>     <IRanges>  <Rle>
#>    [1]        I 151583-151641      +
#>    [2]       II 238361-238419      +
#>    [3]      III 114322-114380      +
#>    [4]       IV 449879-449937      +
#>    [5]        V 152522-152580      +
#>    ...      ...           ...    ...
#>   [12]      XII 151366-151424      +
#>   [13]     XIII 268222-268280      +
#>   [14]      XIV 628588-628646      +
#>   [15]       XV 326897-326955      +
#>   [16]      XVI 556255-556313      +
#>   -------
#>   seqinfo: 17 sequences (1 circular) from R64-1-1 genome
contacts_yeast()
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
#> `HiCExperiment` object with 8,757,906 contacts over 763 regions 
#> -------
#> fileName: "/github/home/.cache/R/ExperimentHub/3e261f048867_7752" 
#> focus: "whole genome" 
#> resolutions(5): 1000 2000 4000 8000 16000
#> active resolution: 16000 
#> interactions: 267709 
#> scores(2): count balanced 
#> topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) centromeres(16) 
#> pairsFile: N/A 
#> metadata(0):
```
