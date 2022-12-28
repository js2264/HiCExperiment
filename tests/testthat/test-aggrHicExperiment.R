test_that("AggrHiCExperiment methods work", {
    microC <- HiCExperiment_example()
    CTCF_peaks <- topologicalFeatures(microC, 'CTCF_peaks')
    x <- AggrHiCExperiment(
        file = fileName(microC), 
        resolution = resolution(microC),
        targets = CTCF_peaks
    )
    expect_s4_class(x, "AggrHiCExperiment")
    expect_s4_class(slices(x), 'SimpleList')
    expect_true(is.array({
        slices(x, 1)
    }))
    expect_true(is.array({
        slices(x, 'count')
    }))
})
