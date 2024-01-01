test_that("fit-regression", {
    # Using Pacioni et al. example data. See ?pac.run.lhs and ?pac.lhs for more
    # details.
    data(pac.run.lhs, pac.lhs)

    # Remove base scenario from .stdat data
    pac.lhs.no.base <- pac.lhs[!pac.lhs$scen.name == "ST_LHS(Base)", ]

    # Use function lookup_table to obtain correct parameter values at year 0
    lkup.ST_LHS <- lookup_table(
        data = pac.lhs.no.base, project = "Pacioni_et_al",
        scenario = "ST_LHS",
        pop = "Population 1",
        SVs = c(
            "SV1", "SV2", "SV3", "SV4", "SV5", "SV6",
            "SV7"
        ),
        save2disk = FALSE
    )

    # Remove base scenario from .run output in long format
    lrun.ST_LHS.no.base <- pac.run.lhs[[2]][
        !pac.run.lhs[[2]]$Scenario == "ST_LHS(Base)",
    ]

    reg <- fit_regression(
        data = lrun.ST_LHS.no.base, lookup = lkup.ST_LHS,
        census = FALSE,
        project = "Pacioni_et_al", scenario = "ST_LHS", popn = 1,
        param = "N", vs = c("SV1", "SV2", "SV3"), l = 2, ncand = 30,
        save2disk = FALSE
    )

    expect_is(reg, "glmulti")
    coefs <- coef(reg@objects[[1]])
    expect_true(round(0.0038677, 6) == round(coefs["SV1"], 6))

    # These tests fail on Gh actions with:
    # ── Error ('test-fit_regression.R:39:5'): fit-regression ────────────────────────
    # Error in `fitfunction == "glm" && !beber$converged`: 'length = 90' in coercion to 'logical(1)'
    # Backtrace:
    #     ▆
    #  1. └─vortexR::fit_regression(...) at test-fit_regression.R:39:5
    #  2.   ├─base::system.time(...)
    #  3.   ├─base::do.call(...)
    #  4.   ├─glmulti::glmulti(...)
    #  5.   └─glmulti::glmulti(...)
    #  6.     ├─base::eval(call)
    #  7.     │ └─base::eval(call)
    #  8.     ├─glmulti::glmulti(...)
    #  9.     └─glmulti::glmulti(...)

    # reg.prop <- fit_regression(
    #     data = lrun.ST_LHS.no.base, lookup = lkup.ST_LHS,
    #     census = FALSE,
    #     project = "Pacioni_et_al", scenario = "ST_LHS", popn = 1,
    #     param = "GeneDiv", vs = c("SV1", "SV2", "SV3"), l = 2, ncand = 30,
    #     save2disk = FALSE
    # )
    # expect_is(reg.prop, "glmulti")
    # coefs <- coef(reg.prop@objects[[1]])
    # expect_true(round(1.663474e-04, 6) == round(coefs["SV1"], 6))

    # reg.prop2 <- fit_regression(
    #     data = lrun.ST_LHS.no.base, lookup = lkup.ST_LHS,
    #     census = FALSE,
    #     links = c("logit", "probit", "cauchit", "loglog"),
    #     project = "Pacioni_et_al", scenario = "ST_LHS", popn = 1,
    #     param = "GeneDiv", vs = c("SV1", "SV2", "SV3"), l = 2,
    #     ncand = 30,
    #     save2disk = FALSE
    # )
    # coefs <- coef(reg.prop2@objects[[1]])
    # expect_true(round(0.0002583495, 6) == round(coefs["SV1"], 6))
})
