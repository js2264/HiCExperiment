test_that("CoolFile works", {
    cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
    mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')

    expect_s4_class(CoolFile(cool_path), "CoolFile")
    expect_error(CoolFile(cool_path, 1000))

    expect_error(CoolFile(mcool_path))
    expect_s4_class(CoolFile(mcool_path, 1000), "CoolFile")
    expect_s4_class(CoolFile(mcool_path, 2000), "CoolFile")
    expect_error(CoolFile(mcool_path, 2))

    expect_s4_class(import(CoolFile(mcool_path, 8000)), "HiCExperiment")

})