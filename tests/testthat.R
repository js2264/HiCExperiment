# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(HiCExperiment)

coolPath <- HiContactsData::HiContactsData('yeast_wt', 'cool')
mcoolPath <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
hicPath <- HiContactsData::HiContactsData('yeast_wt', 'hic')
hicproMatrixPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
hicproBedPath <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
pairsPath <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
pairsPath <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')

test_check("HiCExperiment")
