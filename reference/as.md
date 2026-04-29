# Coercing functions

Coercing functions available for HiCExperiment objects.

## Usage

``` r
# S4 method for class 'HiCExperiment'
as.matrix(x, use.scores = "balanced", sparse = FALSE)

# S4 method for class 'HiCExperiment'
as.data.frame(x)

gi2cm(gi, use.scores = "score")

cm2matrix(cm, replace_NA = NA, sparse = FALSE)

df2gi(
  df,
  seqnames1 = "seqnames1",
  start1 = "start1",
  end1 = "end1",
  seqnames2 = "seqnames2",
  start2 = "start2",
  end2 = "end2"
)
```

## Arguments

- x:

  HiCExperiment object

- use.scores:

  Which scores to use to inflate GInteractions

- sparse:

  Whether to return the contact matrix as a sparse matrix

- gi:

  GInteractions object

- cm:

  A `ContactMatrix` object

- replace_NA:

  Replace NA values

- df:

  A `data.frame` object

- seqnames1, start1, end1, seqnames2, start2, end2:

  Names (as strings) of columns containing corresponding information in
  a data.frame parsed into GInteractions (default: FALSE)

## Examples

``` r
mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
contacts <- import(mcoolPath, focus = 'XVI', resolution = 16000, format = 'cool')
gis <- interactions(contacts)
cm <- gi2cm(gis, 'balanced')
cm
#> class: ContactMatrix 
#> dim: 60 60 
#> type: dgCMatrix 
#> rownames: NULL
#> colnames: NULL
#> metadata(0):
#> regions: 60
cm2matrix(cm)[1:10, 1:10]
#>       [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
#>  [1,]  NaN        NaN        NaN        NaN        NaN        NaN        NaN
#>  [2,]  NaN 1.39502102 1.89623124 0.39560223 0.14820676 0.07426750 0.05980846
#>  [3,]  NaN 1.89623124 2.35751077 1.58129178 0.33905391 0.14339992 0.10552853
#>  [4,]  NaN 0.39560223 1.58129178 1.20960316 0.54168680 0.18242094 0.10413483
#>  [5,]  NaN 0.14820676 0.33905391 0.54168680 0.52832563 0.49867956 0.16116603
#>  [6,]  NaN 0.07426750 0.14339992 0.18242094 0.49867956 0.67745097 0.41256806
#>  [7,]  NaN 0.05980846 0.10552853 0.10413483 0.16116603 0.41256806 0.87399239
#>  [8,]  NaN 0.03340782 0.05614098 0.06121680 0.10343933 0.26074267 0.48711735
#>  [9,]  NaN 0.02780230 0.03659704 0.03372027 0.05460712 0.11214948 0.18101311
#> [10,]  NaN 0.02204195 0.03877495 0.02539418 0.02878660 0.04762856 0.08132232
#>             [,8]       [,9]      [,10]
#>  [1,]        NaN        NaN        NaN
#>  [2,] 0.03340782 0.02780230 0.02204195
#>  [3,] 0.05614098 0.03659704 0.03877495
#>  [4,] 0.06121680 0.03372027 0.02539418
#>  [5,] 0.10343933 0.05460712 0.02878660
#>  [6,] 0.26074267 0.11214948 0.04762856
#>  [7,] 0.48711735 0.18101311 0.08132232
#>  [8,] 0.59450224 0.52204112 0.16856414
#>  [9,] 0.52204112 0.61170845 0.48810265
#> [10,] 0.16856414 0.48810265 0.64767893
df2gi(data.frame(
    chr1 = 'I', start1 = 10, end1 = 100, 
    chr2 = 'I', start2 = 40, end2 = 1000, 
    score = 12, 
    weight = 0.234, 
    filtered = TRUE
), seqnames1 = 'chr1', seqnames2 = 'chr2')
#> GInteractions object with 1 interaction and 3 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2 |     score    weight
#>           <Rle> <IRanges>         <Rle> <IRanges> | <numeric> <numeric>
#>   [1]         I    10-100 ---         I   40-1000 |        12     0.234
#>        filtered
#>       <logical>
#>   [1]      TRUE
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
