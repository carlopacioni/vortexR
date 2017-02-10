library(vortexR)
context("test fit_regression")


test_that("fit-regressin", {
    # Using Pacioni et al. example data. See ?pac.clas for more details.
    data(pac.clas)
    recov <- rRec(pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
                  ST=TRUE, runs=3, yr0=1, yrt=120, save2disk=FALSE)
    expect_is(reg, "glmulti")
})

