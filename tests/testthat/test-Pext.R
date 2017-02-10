library(vortexR)
context("test Pextinct")


test_that("Pextinct", {
    # Using Pacioni et al. example data. See ?pac.run.lhs for more details.
    data(pac.run.lhs)
    Pext <- Pextinct(pac.run.lhs[[2]], project="Pacioni_et_al", scenario="ST_Classic",
                     ST=TRUE, save2disk=FALSE, dir_out="DataAnalysis/Pextinct")

    expect_is(Pext, "data.frame")
    expect_equal(dim(Pext), c(24, 8))
})

