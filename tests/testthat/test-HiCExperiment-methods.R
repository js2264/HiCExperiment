test_that("HiCExperiment methods work", {
    cool_path <- HiContactsData::HiContactsData('yeast_wt', 'cool')
    mcool_path <- HiContactsData::HiContactsData('yeast_wt', 'mcool')
    contacts_yeast <- contacts_yeast()
    expect_identical(length(contacts_yeast), 2945402L)
    expect_s4_class(contacts_yeast[seq_len(10)], 'HiCExperiment')
    expect_s4_class({
        sub <- c(
            rep(TRUE, length(contacts_yeast)/2), 
            rep(FALSE, length(contacts_yeast)/2)
        )
        contacts_yeast[sub]
    }, 'HiCExperiment')
    expect_s4_class({
        contacts_yeast['II:1-10000']
    }, 'HiCExperiment')
    expect_s4_class({
        contacts_yeast['II:1-10000|II:20000-40000']
    }, 'HiCExperiment')
    expect_type(fileName(contacts_yeast), 'character')
    expect_s4_class(seqinfo(contacts_yeast), 'Seqinfo')
    expect_type(resolutions(contacts_yeast), 'integer')
    expect_equal(resolution(contacts_yeast), 1000L)
    expect_s4_class(bins(contacts_yeast), 'GRanges')
    expect_null(focus(contacts_yeast))
    expect_s4_class(interactions(contacts_yeast), 'GInteractions')
    expect_s4_class(scores(contacts_yeast), 'SimpleList')
    expect_type(scores(contacts_yeast, 1), 'double')
    expect_type(scores(contacts_yeast, 'balanced'), 'double')
    expect_s4_class(topologicalFeatures(contacts_yeast), 'SimpleList')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'loops'), 'Pairs')
    expect_s4_class(topologicalFeatures(contacts_yeast, 'borders'), 'GRanges')
    expect_type(pairsFile(contacts_yeast), 'NULL')
    expect_type(anchors(contacts_yeast), 'list')
})

test_that("checks work", {
    contacts_yeast <- contacts_yeast()
    full_contacts_yeast <- full_contacts_yeast()
    expect_true(check_resolution(contacts_yeast, 2000))
    expect_error(check_resolution(contacts_yeast, 3000))
    expect_error(is_comparable(contacts_yeast, full_contacts_yeast))
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
    contacts_yeast <- contacts_yeast()
    expect_s4_class({
        as(contacts_yeast, 'GInteractions')
    }, 'GInteractions')
    expect_s4_class({
        as(contacts_yeast, 'ContactMatrix')
    }, 'ContactMatrix')
    expect_true(is.matrix({
        as(contacts_yeast, 'matrix')
    }))
    expect_s3_class({
        as(contacts_yeast, 'data.frame')
    }, 'data.frame')
})
