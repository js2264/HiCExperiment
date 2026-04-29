# HiCExperiment import methods

Import methods to parse Hi-C files (`.(m)cool`, `.hic`, HiC-Pro derived
matrices, pairs files) into data structures implemented in the
HiCExperiment package.

## Usage

``` r
import(con, format, text, ...)

# S4 method for class 'ANY'
availableResolutions(x, ...)

# S4 method for class 'CoolFile'
availableResolutions(x)

# S4 method for class 'HicFile'
availableResolutions(x)

# S4 method for class 'HicproFile'
availableResolutions(x)

# S4 method for class 'ANY'
availableChromosomes(x, ...)

# S4 method for class 'CoolFile'
availableChromosomes(x)

# S4 method for class 'HicFile'
availableChromosomes(x)

# S4 method for class 'HicproFile'
availableChromosomes(x)
```

## Arguments

- ...:

  Extra parameters to pass to format-specific methods. A list of
  possible arguments is provided in the next section.

- con, x:

  Path or connection to a cool, mcool, .hic or HiC-Pro derived files.
  Can also be a path to a pairs file.

- format:

  The format of the output. If missing and 'con' is a filename, the
  format is derived from the file extension. This argument is
  unnecessary when files are directly provided as `CoolFile`, `HicFile`,
  `HicproFile` or `PairsFile`.

- text:

  If 'con' is missing, this can be a character vector directly providing
  the string data to import.

## Value

A `HiCExperiment` or `GInteractions` object

## import arguments for ContactFile class

`ContactFile` class gathers `CoolFile`, `HicFile` and `HicproFile`
classes. When importing a `ContactFile` object in R, two main arguments
can be provided besides the `ContactFile` itself:

- `resolution`: Resolutions available in the disk-stored contact matrix
  can be listed using `availableResolutions(file)`

- `focus`: A genomic locus (or pair of loci) provided as a string. It
  can be any of the following string structures:

  - `"II"` or `"II:20001-30000"`: this will extract a symmetrical square
    HiCExperiment object, of an entire chromosome or an portion of it.

  - `"II|III"` or `"II:20001-30000|III:40001-90000"`: this will extract
    a non-symmetrical HiCExperiment object, with an entire or portion of
    different chromosomes on each axis.

## Examples

``` r
################################################################
## ----------- Importing .(m)cool contact matrices ---------- ##
################################################################

mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
availableResolutions(mcoolPath)
#> resolutions(5): 1000 2000 4000 8000 16000
#> 
availableChromosomes(mcoolPath)
#> Seqinfo object with 16 sequences from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   I            230218       <NA>   <NA>
#>   II           813184       <NA>   <NA>
#>   III          316620       <NA>   <NA>
#>   IV          1531933       <NA>   <NA>
#>   V            576874       <NA>   <NA>
#>   ...             ...        ...    ...
#>   XII         1078177       <NA>   <NA>
#>   XIII         924431       <NA>   <NA>
#>   XIV          784333       <NA>   <NA>
#>   XV          1091291       <NA>   <NA>
#>   XVI          948066       <NA>   <NA>
import(mcoolPath, resolution = 16000, focus = 'XVI', format = 'cool')
#> `HiCExperiment` object with 535,350 contacts over 60 regions 
#> -------
#> fileName: "/github/home/.cache/R/ExperimentHub/3e261f048867_7752" 
#> focus: "XVI" 
#> resolutions(5): 1000 2000 4000 8000 16000
#> active resolution: 16000 
#> interactions: 1731 
#> scores(2): count balanced 
#> topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
#> pairsFile: N/A 
#> metadata(0):

################################################################
## ------------ Importing .hic contact matrices ------------- ##
################################################################

hicPath <- HiContactsData::HiContactsData('yeast_wt', 'hic')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
availableResolutions(hicPath)
#> resolutions(5): 1000 2000 4000 8000 16000
#> 
availableChromosomes(hicPath)
#> Seqinfo object with 17 sequences from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   I            230218       <NA>   <NA>
#>   II           813184       <NA>   <NA>
#>   III          316620       <NA>   <NA>
#>   IV          1531933       <NA>   <NA>
#>   IX           439888       <NA>   <NA>
#>   ...             ...        ...    ...
#>   XIII         924431       <NA>   <NA>
#>   XIV          784333       <NA>   <NA>
#>   XV          1091291       <NA>   <NA>
#>   XVI          948066       <NA>   <NA>
#>   M             85780       <NA>   <NA>
import(hicPath, resolution = 16000, focus = 'XVI', format = 'hic')
#> `HiCExperiment` object with 838,222 contacts over 60 regions 
#> -------
#> fileName: "/github/home/.cache/R/ExperimentHub/3e265ac673a5_7836" 
#> focus: "XVI" 
#> resolutions(5): 1000 2000 4000 8000 16000
#> active resolution: 16000 
#> interactions: 1732 
#> scores(2): count balanced 
#> topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
#> pairsFile: N/A 
#> metadata(0):

################################################################
## ------- Importing HiC-Pro derived contact matrices ------- ##
################################################################

hicproMatrixPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
hicproBedPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
#> see ?HiContactsData and browseVignettes('HiContactsData') for documentation
#> loading from cache
availableResolutions(hicproMatrixPath, hicproBedPath)
#> [1] 1000
availableChromosomes(hicproMatrixPath, hicproBedPath)
#> Seqinfo object with 17 sequences from an unspecified genome:
#>   seqnames seqlengths isCircular genome
#>   I            230218       <NA>   <NA>
#>   II           813184       <NA>   <NA>
#>   III          316620       <NA>   <NA>
#>   IV          1531933       <NA>   <NA>
#>   IX           439888       <NA>   <NA>
#>   ...             ...        ...    ...
#>   XII         1078177       <NA>   <NA>
#>   XIII         924431       <NA>   <NA>
#>   XIV          784333       <NA>   <NA>
#>   XV          1091291       <NA>   <NA>
#>   XVI          948066       <NA>   <NA>
import(hicproMatrixPath, bed = hicproBedPath, format = 'hicpro')
#> Registered S3 methods overwritten by 'readr':
#>   method                    from 
#>   as.data.frame.spec_tbl_df vroom
#>   as_tibble.spec_tbl_df     vroom
#>   format.col_spec           vroom
#>   print.col_spec            vroom
#>   print.collector           vroom
#>   print.date_names          vroom
#>   print.locale              vroom
#>   str.col_spec              vroom
#> 
#> `HiCExperiment` object with 9,503,604 contacts over 12,165 regions 
#> -------
#> fileName: "/github/home/.cache/R/ExperimentHub/3e2639eaefe7_7837" 
#> focus: "whole genome" 
#> resolutions(1): 1000
#> active resolution: 1000 
#> interactions: 2686250 
#> scores(2): count balanced 
#> topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) 
#> pairsFile: N/A 
#> metadata(1): regions
```
