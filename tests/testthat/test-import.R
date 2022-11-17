test_that("imports work", {
    
    cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
    mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    pairs_path <- HiContactsData::HiContactsData('yeast_wt', 'pairs.gz')
    
    expect_s4_class(HiCExperiment(cool_path), 'HiCExperiment')
    expect_s4_class(HiCExperiment(mcool_path, 1000), 'HiCExperiment')
    expect_equal(
        resolution(HiCExperiment(mcool_path, 2000)), 
        2000
    )
    expect_error(HiCExperiment(cool_path, 1000))

    expect_s4_class(HiCExperiment(CoolFile(cool_path)), 'HiCExperiment')
    expect_s4_class(HiCExperiment(CoolFile(mcool_path, 1000)), 'HiCExperiment')
    expect_s4_class(HiCExperiment(CoolFile(mcool_path)), 'HiCExperiment')
    expect_error(HiCExperiment(CoolFile(cool_path, 1000)))

    expect_s4_class(import(CoolFile(cool_path), format = 'cool'), 'HiCExperiment')
    expect_s4_class(import(cool_path, format = 'cool'), 'HiCExperiment')
    expect_s4_class(import(CoolFile(mcool_path), format = 'mcool'), 'HiCExperiment')
    cf <- import(mcool_path, format = 'mcool', pairsFile = pairs_path)
    expect_s4_class(cf, 'HiCExperiment')
    expect_equal(pairsFile(cf), pairs_path)
    expect_equal(
        resolution(import(CoolFile(mcool_path, 2000), format = 'mcool')), 
        2000
    )

    expect_equal(
        resolution(import(mcool_path, format = 'mcool', resolution = 8000)), 
        8000
    )
    expect_error(import(cool_path, format = 'cool', resolution = 8000))

    expect_s4_class(import(pairs_path, format = 'pairs'), 'GenomicInteractions')
    expect_s4_class(import(PairsFile(pairs_path)), 'GenomicInteractions')

})
