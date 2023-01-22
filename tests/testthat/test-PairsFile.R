test_that("PairsFile works", {
    pf <- PairsFile(pairsPath)

    expect_s4_class(pf, "PairsFile")
    expect_no_warning(pairsFile(pf))
    expect_no_error(pairsFile(pf))
    expect_s4_class(import(PairsFile(pairsPath)), "GInteractions")

})
