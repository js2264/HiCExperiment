test_that("PairsFile works", {
    pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')

    expect_s4_class(PairsFile(pairs_path), "PairsFile")
    expect_s4_class(import(PairsFile(pairs_path)), "GenomicInteractions")

})