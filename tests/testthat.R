# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(HiCExperiment)

cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
hicpro_matrix_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
hicpro_bed_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')

test_check("HiCExperiment")
