test_that("HicproFile works", {
    
    hicpro0 <- expect_no_error(HicproFile(hicpro_matrix_path))
    hicpro1 <- HicproFile(hicpro_matrix_path, bed = hicpro_bed_path)
    hicpro2 <- HicproFile(hicpro_matrix_path, bed = hicpro_bed_path, pairs = pairs_path)
    
    # import is properly handled for HicproFile
    ihicpro1 <- import(hicpro1)
    expect_no_error(show(hicpro1))

    # import is properly handled for files provided as string
    ihicpro1 <- import(hicpro_matrix_path, bed = hicpro_bed_path, format = 'hicpro')
    expect_error(import(hicpro_matrix_path))

    # import is properly handled with HiCExperiment function
    ihicpro2 <- HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path)
    expect_error(HiCExperiment(hicpro_matrix_path))
    expect_error(HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path, focus = 'II'))
    expect_error(HiCExperiment(hicpro_matrix_path, bed = hicpro_bed_path, resolution = 160000))

})
