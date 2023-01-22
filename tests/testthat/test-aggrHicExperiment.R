test_that("AggrHiCExperiment methods work", {
    data(centros_yeast)
    x <- AggrHiCExperiment(
        file = mcoolPath, 
        resolution = 8000,
        targets = centros_yeast[c(4, 7)]
    )
    expect_s4_class(x, "AggrHiCExperiment")
    expect_no_error(show(x))
    expect_s4_class(slices(x), 'SimpleList')
    expect_true(is.array({
        slices(x, 1)
    }))
    expect_true(is.array({
        slices(x, 'count')
    }))
})
