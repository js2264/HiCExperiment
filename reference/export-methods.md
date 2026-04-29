# HiCExperiment export methods

Export methods to save a HiCExperiment object into a set of
HiC-Pro-style files (matrix & regions files)

## Usage

``` r
# S4 method for class 'HiCExperiment,missing,character'
export(object, prefix, format, ...)
```

## Arguments

- object:

  A HiCExperiment object

- prefix:

  Prefix used when generating output file(s).

- format:

  File format. Available: `cool` and `HiC-Pro`.

- ...:

  Extra arguments to use when exporting to `cool`. Can be
  `metadata <string>` or `chunksize <integer>`.

## Value

Path to saved files

## Examples

``` r
################################################################
## ----------- Importing .(m)cool contact matrices ---------- ##
################################################################

mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
hic <- import(mcoolPath, format = 'mcool', resolution = 16000)
export(hic["II"], prefix = 'subset_chrII', format = 'cool')
#> Writing /bins table to subset_chrII.cool file...
#> Writing /chroms table to subset_chrII.cool file...
#> Writing /pixels table to subset_chrII.cool file...
#> Writing /indexes table to subset_chrII.cool file...
#> CoolFile object
#> .mcool file: subset_chrII.cool 
#> resolution: 16000 
#> pairs file: 
#> metadata(0):
export(hic["II"], prefix = 'subset_chrII', format = 'HiC-Pro')
#> Writing regions to subset_chrII_regions.bed file...
#> Writing raw count matrix to subset_chrII_matrix.mtx file...
#> HicproFile object
#> HiC-Pro files:
#>   $ matrix:   subset_chrII_matrix.mtx 
#>   $ regions:  subset_chrII_regions.bed 
#> resolution: 16000 
#> pairs file: 
#> metadata(0):
```
