test_that("AggrHiCExperiment methods work", {
    contacts <- contacts_yeast(full = TRUE)
    centros <- topologicalFeatures(contacts, 'centromeres')
    x <- AggrHiCExperiment(
        file = fileName(contacts), 
        resolution = 2000,
        targets = centros[3:8]
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
