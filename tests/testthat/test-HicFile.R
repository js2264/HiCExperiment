test_that("HicFile works", {
    
    expect_no_error(availableChromosomes(hicPath))
    expect_no_error(availableResolutions(hicPath))
    # Prepare paths and HicFiles
    hf0 <- HicFile(hicPath)
    hf2 <- HicFile(hicPath, pairsFile = pairsPath, resolution = 16000)
    hf3 <- HicFile(hicPath, resolution = 16000)
    expect_no_error(availableChromosomes(hf0))
    expect_no_error(availableResolutions(hf0))

    # No errors/warnings when printing HicFile
    expect_no_warning(show(hf0))
    expect_no_warning(hf0)
    expect_no_error(hf0)
    expect_no_warning(hf2)
    expect_no_error(hf2)
    expect_no_warning(hf3)
    expect_no_error(hf3)

    # Resolution(s) are properly handled
    expect_equal(resolution(hf0), 1000L)
    expect_equal(resolution(hf2), 16000L)
    expect_equal(resolution(hf3), 16000L)
    expect_error(HicFile(coolPath, 1000))
    expect_error(HicFile(hicPath, 2))

    # pairFile is properly handled
    expect_null(pairsFile(hf0))
    expect_equal(pairsFile(hf2), pairsPath)
    expect_null(pairsFile(hf3))

    # import is properly handled for HicFile
    expect_no_error(ihf0 <- import(hf0, focus = 'I'))
    expect_no_error(ihf2 <- import(hf2))
    expect_no_error(ihf4 <- import(hf3, focus = 'I'))
    expect_no_error(ihf5 <- import(hf3, resolution = 4000, focus = 'I', metadata = list(test = 'test')))
    expect_equal(resolution(ihf0), 1000L)
    expect_null(pairsFile(ihf0))
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(pairsFile(ihf2), pairsPath)
    expect_equal(metadata(ihf4), list())
    expect_equal(resolution(ihf4), 16000L)
    expect_null(pairsFile(ihf4))
    expect_equal(focus(ihf4), "I")
    expect_equal(resolution(ihf5), 4000L)
    expect_null(pairsFile(ihf5))
    expect_equal(metadata(ihf5), list(test = 'test'))
    expect_equal(focus(ihf5), "I")

    # import is properly handled for files provided as string
    ihf0 <- import(hicPath, focus = 'I', format = 'hic')
    ihf2 <- import(hicPath, pairsFile = pairsPath, resolution = 16000, format = 'hic')
    ihf3 <- import(hicPath, focus = 'I', resolution = 16000, format = 'hic')
    ihf4 <- import(hicPath, focus = 'I', format = 'hic')
    ihf5 <- import(hicPath, resolution = 4000, focus = 'I', metadata = list(test = 'test'), format = 'hic')
    expect_equal(resolution(ihf0), 1000L)
    expect_null(pairsFile(ihf0))
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(pairsFile(ihf2), pairsPath)
    expect_equal(resolution(ihf3), 16000L)
    expect_null(pairsFile(ihf3))
    expect_equal(resolution(ihf4), 1000L)
    expect_null(pairsFile(ihf4))
    expect_equal(focus(ihf4), "I")
    expect_equal(resolution(ihf5), 4000L)
    expect_null(pairsFile(ihf5))
    expect_equal(metadata(ihf5), list(test = 'test'))
    expect_equal(focus(ihf5), "I")

    # import is properly handled with HiCExperiment function
    expect_error(HiCExperiment(hicPath))
    expect_error(HiCExperiment(hicPath, focus = 'I'))
    ihf2 <- HiCExperiment(hicPath, focus = 'I', resolution = 16000)
    ihf3 <- HiCExperiment(hicPath, focus = 'I', resolution = 16000, pairsFile = pairsPath)
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(focus(ihf2), 'I')
    expect_null(pairsFile(ihf2))
    expect_equal(resolution(ihf3), 16000L)
    expect_equal(pairsFile(ihf3), pairsPath)
    expect_equal(resolution(ihf3), 16000L)

    # Import works for all possible combinations of ranges 
    expect_no_error(import(hicPath, format = 'hic', focus = 'I'))
    expect_no_error(import(hicPath, format = 'hic', focus = 'I:1-40000'))
    expect_no_error(import(hicPath, format = 'hic', focus = 'I|II'))
    expect_no_error(import(hicPath, format = 'hic', focus = 'III|II'))
    expect_no_error(import(hicPath, format = 'hic', focus = 'III:30000-40000|II:10000-40000'))

})
