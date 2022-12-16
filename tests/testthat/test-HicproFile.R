test_that("HicFile works", {
    
    # Prepare paths and HicproFiles
    hicpro_matrix_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_matrix')
    hicpro_bed_path <- HiContactsData::HiContactsData('yeast_wt', 'hicpro_bed')
    pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')

    hicpro0 <- HicproFile(hicpro_matrix_path)
    hicpro1 <- HicproFile(hicpro_matrix_path, bed = hicpro_bed_path)
    hicpro2 <- HicproFile(hicpro_matrix_path, bed = hicpro_bed_path, pairs = pairs_path)
    
    # import is properly handled for HicproFile
    expect_error(import(hicpro0))
    ihicpro1 <- import(hicpro1)

    # import is properly handled for files provided as string
    expect_error(import(hicpro_matrix_path))
    ihicpro1 <- import(hicpro_matrix_path, bed = hicpro_bed_path, format = 'hicpro')

    # import is properly handled with HiCExperiment function
    expect_error(HiCExperiment(hicpro_matrix_path))
    ihicpro2 <- HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path)
    expect_error(HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path, focus = 'II'))
    expect_error(HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path, resolution = 160000))

})
