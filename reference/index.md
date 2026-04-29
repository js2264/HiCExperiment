# Package index

## Importing Hi-C files

- [`AllGenerics`](AllGenerics.md)
  [`availableResolutions`](AllGenerics.md)
  [`availableChromosomes`](AllGenerics.md)
  [`resolutions`](AllGenerics.md) [`resolution`](AllGenerics.md)
  [`focus`](AllGenerics.md) [`focus<-`](AllGenerics.md)
  [`scores`](AllGenerics.md) [`scores<-`](AllGenerics.md)
  [`topologicalFeatures`](AllGenerics.md)
  [`topologicalFeatures<-`](AllGenerics.md)
  [`pairsFile`](AllGenerics.md) [`pairsFile<-`](AllGenerics.md)
  [`metadata<-`](AllGenerics.md) [`bins`](AllGenerics.md)
  [`slices`](AllGenerics.md) [`zoom`](AllGenerics.md)
  [`refocus`](AllGenerics.md) [`cis`](AllGenerics.md)
  [`trans`](AllGenerics.md) [`bin`](AllGenerics.md) : Generic functions
- [`import()`](import-methods.md)
  [`availableResolutions(`*`<ANY>`*`)`](import-methods.md)
  [`availableResolutions(`*`<CoolFile>`*`)`](import-methods.md)
  [`availableResolutions(`*`<HicFile>`*`)`](import-methods.md)
  [`availableResolutions(`*`<HicproFile>`*`)`](import-methods.md)
  [`availableChromosomes(`*`<ANY>`*`)`](import-methods.md)
  [`availableChromosomes(`*`<CoolFile>`*`)`](import-methods.md)
  [`availableChromosomes(`*`<HicFile>`*`)`](import-methods.md)
  [`availableChromosomes(`*`<HicproFile>`*`)`](import-methods.md) :
  HiCExperiment import methods

## Exporting Hi-C files

- [`export(`*`<HiCExperiment>`*`,`*`<missing>`*`,`*`<character>`*`)`](export-methods.md)
  : HiCExperiment export methods

## Methods for `HiCExperiment` S4 class

- [`HiCExperiment()`](HiCExperiment.md)
  [`makeHiCExperimentFromGInteractions()`](HiCExperiment.md)
  [`resolutions(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`resolution(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`focus(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`` `focus<-`( ``*`<HiCExperiment>`*`,`*`<character>`*`)`](HiCExperiment.md)
  [`zoom(`*`<HiCExperiment>`*`,`*`<numeric>`*`)`](HiCExperiment.md)
  [`refocus(`*`<HiCExperiment>`*`,`*`<character>`*`)`](HiCExperiment.md)
  [`scores(`*`<HiCExperiment>`*`,`*`<missing>`*`)`](HiCExperiment.md)
  [`scores(`*`<HiCExperiment>`*`,`*`<character>`*`)`](HiCExperiment.md)
  [`scores(`*`<HiCExperiment>`*`,`*`<numeric>`*`)`](HiCExperiment.md)
  [`` `scores<-`( ``*`<HiCExperiment>`*`,`*`<character>`*`,`*`<numeric>`*`)`](HiCExperiment.md)
  [`topologicalFeatures(`*`<HiCExperiment>`*`,`*`<missing>`*`)`](HiCExperiment.md)
  [`topologicalFeatures(`*`<HiCExperiment>`*`,`*`<character>`*`)`](HiCExperiment.md)
  [`topologicalFeatures(`*`<HiCExperiment>`*`,`*`<numeric>`*`)`](HiCExperiment.md)
  [`` `topologicalFeatures<-`( ``*`<HiCExperiment>`*`,`*`<character>`*`,`*`<GRangesOrGInteractions>`*`)`](HiCExperiment.md)
  [`pairsFile(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`` `pairsFile<-`( ``*`<HiCExperiment>`*`,`*`<character>`*`)`](HiCExperiment.md)
  [`` `metadata<-`( ``*`<HiCExperiment>`*`,`*`<list>`*`)`](HiCExperiment.md)
  [`subsetByOverlaps(`*`<HiCExperiment>`*`,`*`<numeric>`*`)`](HiCExperiment.md)
  [`subsetByOverlaps(`*`<HiCExperiment>`*`,`*`<logical>`*`)`](HiCExperiment.md)
  [`subsetByOverlaps(`*`<HiCExperiment>`*`,`*`<GRanges>`*`)`](HiCExperiment.md)
  [`subsetByOverlaps(`*`<HiCExperiment>`*`,`*`<GInteractions>`*`)`](HiCExperiment.md)
  [`subsetByOverlaps(`*`<HiCExperiment>`*`,`*`<Pairs>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<numeric>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<GRanges>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<logical>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<GInteractions>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<Pairs>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`` `[`( ``*`<HiCExperiment>`*`,`*`<character>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](HiCExperiment.md)
  [`fileName(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`interactions(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`` `interactions<-`( ``*`<HiCExperiment>`*`,`*`<GInteractions>`*`)`](HiCExperiment.md)
  [`length(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`` `$<-`( ``*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`` `$`( ``*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`seqinfo(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`bins(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`anchors(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`regions(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`cis(`*`<HiCExperiment>`*`)`](HiCExperiment.md)
  [`trans(`*`<HiCExperiment>`*`)`](HiCExperiment.md) :

  `HiCExperiment` S4 class

- [`AggrHiCExperiment()`](AggrHiCExperiment.md)
  [`slices(`*`<AggrHiCExperiment>`*`,`*`<missing>`*`)`](AggrHiCExperiment.md)
  [`slices(`*`<AggrHiCExperiment>`*`,`*`<character>`*`)`](AggrHiCExperiment.md)
  [`slices(`*`<AggrHiCExperiment>`*`,`*`<numeric>`*`)`](AggrHiCExperiment.md)
  [`show(`*`<AggrHiCExperiment>`*`)`](AggrHiCExperiment.md) :

  `AggrHiCExperiment` S4 class

## S4 Classes for Hi-C contact matrix file formats

- [`ContactsFile-class`](ContactsFile-class.md)
  [`ContactsFile`](ContactsFile-class.md)
  [`ContactsFile-methods`](ContactsFile-class.md)
  [`pairsFile,ContactsFile-method`](ContactsFile-class.md)
  [`resolution,ContactsFile-method`](ContactsFile-class.md)
  [`metadata<-,ContactsFile-method`](ContactsFile-class.md)
  [`metadata<-,ContactsFile,list-method`](ContactsFile-class.md) :

  `ContactsFile` S4 class

- [`CoolFile-class`](CoolFile-class.md)
  [`McoolFile-class`](CoolFile-class.md) [`CoolFile`](CoolFile-class.md)
  [`CoolFile-methods`](CoolFile-class.md)
  [`show,CoolFile-method`](CoolFile-class.md) :

  `CoolFile` S4 class

- [`HicFile-class`](HicFile-class.md) [`HicFile`](HicFile-class.md)
  [`HicFile-methods`](HicFile-class.md)
  [`show,HicFile-method`](HicFile-class.md) :

  `HicFile` S4 class

- [`HicproFile-class`](HicproFile-class.md)
  [`HicproFile`](HicproFile-class.md)
  [`HicproFile-methods`](HicproFile-class.md)
  [`show,HicproFile-method`](HicproFile-class.md) :

  `HicproFile` S4 class

- [`PairsFile-class`](PairsFile-class.md)
  [`PairsFile`](PairsFile-class.md)
  [`pairsFile,PairsFile-method`](PairsFile-class.md) :

  `PairsFile` S4 class

## Coercing methods

- [`bin(`*`<GInteractions>`*`,`*`<numeric>`*`)`](bin-methods.md)
  [`bin(`*`<PairsFile>`*`,`*`<numeric>`*`)`](bin-methods.md) :
  HiCExperiment binning methods
- [`as.matrix(`*`<HiCExperiment>`*`)`](as.md)
  [`as.data.frame(`*`<HiCExperiment>`*`)`](as.md) [`gi2cm()`](as.md)
  [`cm2matrix()`](as.md) [`df2gi()`](as.md) : Coercing functions

## Example datasets

- [`centros_yeast`](data.md) [`contacts_yeast()`](data.md)
  [`contacts_yeast_eco1()`](data.md) :

  Example datasets provided in `HiCExperiment` & `HiContactsData`
