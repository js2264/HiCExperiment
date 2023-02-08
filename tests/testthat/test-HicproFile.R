test_that("HicproFile works", {
    
    expect_no_error(availableChromosomes(hicproMatrixPath, hicproBedPath))
    expect_no_error(availableResolutions(hicproMatrixPath, hicproBedPath))
    hicpro0 <- expect_no_error(HicproFile(hicproMatrixPath))
    hicpro1 <- HicproFile(hicproMatrixPath, bed = hicproBedPath)
    hicpro2 <- HicproFile(hicproMatrixPath, bed = hicproBedPath, pairs = pairsPath)
    expect_no_error(availableChromosomes(hicpro2))
    expect_no_error(availableResolutions(hicpro2))
    
    # # import is properly handled for HicproFile
    ihicpro1 <- import(hicpro1)
    expect_no_error(show(hicpro1))

    # # import is properly handled for files provided as string
    ihicpro1 <- import(hicproMatrixPath, bed = hicproBedPath, format = 'hicpro')
    expect_error(import(hicproMatrixPath))

    # # import is properly handled with HiCExperiment function
    ihicpro2 <- HiCExperiment(hicproMatrixPath, bed = hicproBedPath)
    expect_error(HiCExperiment(hicproMatrixPath))
    expect_error(HiCExperiment(hicproMatrixPath, bed = hicproBedPath, focus = 'II'))
    expect_error(HiCExperiment(hicproMatrixPath, bed = hicproBedPath, resolution = 160000))

})
