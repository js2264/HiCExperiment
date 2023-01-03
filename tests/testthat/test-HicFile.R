test_that("HicFile works", {
    
    # Prepare paths and HicFiles
    hic_path <- HiContactsData::HiContactsData('yeast_wt', 'hic')
    pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
    hf0 <- HicFile(hic_path)
    hf2 <- HicFile(hic_path, pairsFile = pairs_path, resolution = 16000)
    hf3 <- HicFile(hic_path, resolution = 16000)

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
    expect_error(HicFile(cool_path, 1000))
    expect_error(HicFile(hic_path, 2))

    # pairFile is properly handled
    expect_null(pairsFile(hf0))
    expect_equal(pairsFile(hf2), pairs_path)
    expect_null(pairsFile(hf3))

    # import is properly handled for HicFile
    ihf0 <- import(hf0, focus = 'I')
    ihf2 <- import(hf2)
    ihf3 <- import(hf3, focus = 'I')
    ihf4 <- import(hf3, focus = 'I')
    ihf5 <- import(hf3, resolution = 4000, focus = 'I', metadata = list(test = 'test'))
    expect_equal(resolution(ihf0), 1000L)
    expect_null(pairsFile(ihf0))
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(pairsFile(ihf2), pairs_path)
    expect_equal(resolution(ihf3), 16000L)
    expect_null(pairsFile(ihf3))
    expect_equal(metadata(ihf3), list())
    expect_equal(resolution(ihf4), 16000L)
    expect_null(pairsFile(ihf4))
    expect_equal(focus(ihf4), "I")
    expect_equal(resolution(ihf5), 4000L)
    expect_null(pairsFile(ihf5))
    expect_equal(metadata(ihf5), list(test = 'test'))
    expect_equal(focus(ihf5), "I")

    # import is properly handled for files provided as string
    ihf0 <- import(hic_path, focus = 'I', format = 'hic')
    ihf2 <- import(hic_path, pairsFile = pairs_path, resolution = 16000, format = 'hic')
    ihf3 <- import(hic_path, focus = 'I', resolution = 16000, format = 'hic')
    ihf4 <- import(hic_path, focus = 'I', format = 'hic')
    ihf5 <- import(hic_path, resolution = 4000, focus = 'I', metadata = list(test = 'test'), format = 'hic')
    expect_equal(resolution(ihf0), 1000L)
    expect_null(pairsFile(ihf0))
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(pairsFile(ihf2), pairs_path)
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
    expect_error(HiCExperiment(hic_path))
    expect_error(HiCExperiment(hic_path, focus = 'I'))
    ihf2 <- HiCExperiment(hic_path, focus = 'I', resolution = 16000)
    ihf3 <- HiCExperiment(hic_path, focus = 'I', resolution = 16000, pairsFile = pairs_path)
    expect_equal(resolution(ihf2), 16000L)
    expect_equal(focus(ihf2), 'I')
    expect_null(pairsFile(ihf2))
    expect_equal(resolution(ihf3), 16000L)
    expect_equal(pairsFile(ihf3), pairs_path)
    expect_equal(resolution(ihf3), 16000L)

})
