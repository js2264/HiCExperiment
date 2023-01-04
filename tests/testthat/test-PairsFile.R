test_that("PairsFile works", {
    pf <- PairsFile(pairs_path)

    expect_s4_class(pf, "PairsFile")
    expect_no_warning(pairsFile(pf))
    expect_no_error(pairsFile(pf))
    expect_s4_class(import(PairsFile(pairs_path)), "GInteractions")

})
