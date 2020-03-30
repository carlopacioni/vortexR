library(vortexR)
library(vortexRdata)
context("test fit_regression")


test_that("fit-regressin", {
    # Using Pacioni et al. example data. See ?pac.run.lhs and ?pac.lhs for more
    # details.
    data(pac.run.lhs, pac.lhs)

    # Remove base scenario from .stdat data
    pac.lhs.no.base <- pac.lhs[!pac.lhs$scen.name == "ST_LHS(Base)", ]

    # Use function lookup_table to obtain correct parameter values at year 0
    lkup.ST_LHS <- lookup_table(data=pac.lhs.no.base, project="Pacioni_et_al",
                                scenario="ST_LHS",
                                pop="Population 1",
                                SVs=c("SV1", "SV2", "SV3", "SV4", "SV5", "SV6",
                                      "SV7"),
                                save2disk=FALSE)

    # Remove base scenario from .run output in long format
    lrun.ST_LHS.no.base <- pac.run.lhs[[2]][
        !pac.run.lhs[[2]]$Scenario == "ST_LHS(Base)",
        ]

    reg <- fit_regression(data=lrun.ST_LHS.no.base, lookup=lkup.ST_LHS,
                          census=FALSE,
                          project="Pacioni_et_al", scenario="ST_LHS", popn=1,
                          param="N", vs=c("SV1", "SV2", "SV3"), l=2,  ncand=30,
                          save2disk=FALSE)

    expect_is(reg, "glmulti")
    coefs <- coef(reg@objects[[1]])
    expect_true(round(0.0038677, 6)==round(coefs["SV1"], 6))


    reg.prop <- fit_regression(data=lrun.ST_LHS.no.base, lookup=lkup.ST_LHS,
                          census=FALSE,
                          project="Pacioni_et_al", scenario="ST_LHS", popn=1,
                          param="GeneDiv", vs=c("SV1", "SV2", "SV3"), l=2,  ncand=30,
                          save2disk=FALSE)
    expect_is(reg.prop, "glmulti")
    coefs <- coef(reg.prop@objects[[1]])
    expect_true(round(1.663474e-04, 6)==round(coefs["SV1"], 6))

})

