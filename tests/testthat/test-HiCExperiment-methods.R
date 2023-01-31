contacts_yeast <- contacts_yeast()

test_that("HiCExperiment methods work", {
    expect_identical(length(contacts_yeast), 267709L)
    expect_s4_class(contacts_yeast[seq_len(10)], 'HiCExperiment')
    expect_s4_class({
        sub <- c(
            rep(TRUE, length(contacts_yeast)/2), 
            rep(FALSE, length(contacts_yeast)/2)
        )
        contacts_yeast[sub]
    }, 'HiCExperiment')
    expect_s4_class({
        contacts_yeast['II:1-20000']
    }, 'HiCExperiment')
    expect_s4_class({
        contacts_yeast['II:1-30000|II:30000-60000']
    }, 'HiCExperiment')
    expect_type(fileName(contacts_yeast), 'character')
    expect_s4_class(seqinfo(contacts_yeast), 'Seqinfo')
    expect_type(resolutions(contacts_yeast), 'integer')
    expect_equal(resolution(contacts_yeast), 16000L)
    expect_s4_class(bins(contacts_yeast), 'GRanges')
    expect_null(focus(contacts_yeast))
    expect_no_error(focus(contacts_yeast) <- 'II')
    expect_no_error(refocus(contacts_yeast, 'III'))
    expect_s4_class(zoom(contacts_yeast, 16000), 'HiCExperiment')
    expect_s4_class(interactions(contacts_yeast), 'GInteractions')
    expect_no_error(interactions(contacts_yeast) <- interactions(contacts_yeast))
    expect_s4_class(scores(contacts_yeast), 'SimpleList')
    expect_type(scores(contacts_yeast, 1), 'double')
    expect_type(scores(contacts_yeast, 'balanced'), 'double')
    expect_no_error(scores(contacts_yeast, 'test') <- runif(length(contacts_yeast)))
    expect_s4_class(topologicalFeatures(contacts_yeast), 'SimpleList')
    expect_s4_class(topologicalFeatures(contacts_yeast, 2), 'GRanges')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'loops'), 'GRangesOrGInteractions')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'borders'), 'GRanges')
    expect_no_error(topologicalFeatures(contacts_yeast, 'test') <- GenomicRanges::GRanges())
    expect_type(pairsFile(contacts_yeast), 'NULL')
    expect_no_error(pairsFile(contacts_yeast) <- pairsPath)
    expect_type(anchors(contacts_yeast), 'list')
    expect_s4_class(regions(contacts_yeast), 'GRanges')
    expect_type(metadata(contacts_yeast), 'list')
    expect_no_error(metadata(contacts_yeast) <- list(test = 'OK'))
    expect_no_error(show(contacts_yeast))
})

test_that("checks work", {
    expect_true(check_resolution(contacts_yeast, 2000))
    expect_error(check_resolution(contacts_yeast, 3000))
    expect_true(is_square(S4Vectors::Pairs(
        first = 'I:10000-20000', 
        second = 'I:10000-20000'
    )))
})

test_that("utils works", {

    expect_type(splitCoords(GRanges('II:30-40')), 'list')
    expect_type(splitCoords('II:30-40'), 'list')

    expect_type(coords2char('II:30-40'), 'character')
    expect_type(coords2char('II:30-40|II:40-50'), 'character')
    expect_type(coords2char(GRanges('II:30-40')), 'character')

    expect_s4_class(char2coords("II:30000-50000|II:60000-80000"), 'Pairs')

    expect_equal(
        sortPairs(char2coords("II:30000-50000|II:60000-80000")), 
        char2coords("II:30000-50000|II:60000-80000")
    )
    expect_equal(
        sortPairs(char2coords("II:90000-100000|II:60000-80000")), 
        char2coords("II:60000-80000|II:90000-100000")
    )

    expect_s4_class(
        asGInteractions(
            data.frame(
                seqnames1 = 'II', 
                start1 = 23, 
                end1 = 34, 
                seqnames2 = 'II', 
                start2 = 12, 
                end2 = 54
            )
        ), 
        'GInteractions'
    )
})

test_that("coerce works", {
    expect_s4_class({
        as(contacts_yeast, 'GInteractions')
    }, 'GInteractions')
    expect_s4_class({
        as(contacts_yeast, 'ContactMatrix')
    }, 'ContactMatrix')
    expect_true(is.matrix(base::as.matrix(as.matrix(contacts_yeast))))
    expect_s4_class({
        as(contacts_yeast, 'matrix')
    }, 'dgTMatrix')
    expect_s4_class({
        as.matrix(contacts_yeast)
    }, 'dgTMatrix')
    expect_s3_class({
        as(contacts_yeast, 'data.frame')
    }, 'data.frame')
})
