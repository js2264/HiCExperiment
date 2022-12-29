test_that("AggrHiCExperiment methods work", {
    contacts <- full_contacts_yeast()
    centros <- topologicalFeatures(contacts, 'centromeres')
    x <- AggrHiCExperiment(
        file = fileName(contacts), 
        resolution = 1000,
        targets = centros
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
