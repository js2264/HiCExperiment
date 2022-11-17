test_that("CoolFile works", {
    cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
    mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    cf <- CoolFile(cool_path)
    mcf <- CoolFile(mcool_path, resolution = 2000)

    expect_s4_class(CoolFile(cool_path), "CoolFile")
    expect_error(CoolFile(cool_path, 1000))

    expect_s4_class(CoolFile(mcool_path), "CoolFile")
    expect_s4_class(CoolFile(mcool_path, 1000), "CoolFile")
    expect_s4_class(CoolFile(mcool_path, 2000), "CoolFile")
    expect_error(CoolFile(mcool_path, 2))

    expect_s4_class(import(CoolFile(mcool_path, 8000)), "HiCExperiment")
    expect_equal(resolution(import(CoolFile(mcool_path, 8000))), 8000L)

    expect_equal(resolution(cf), NULL)
    expect_equal(resolution(mcf), 2000L)
    expect_no_error(cf)
    expect_no_warning(cf)

})