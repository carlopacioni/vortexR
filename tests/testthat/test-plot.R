library(vortexR)
context("test plots")


test_that("test dot_plot", {
    data(pac.clas)
    dot <- dot_plot(data=pac.clas, project="Pacioni_et_al", scenario="ST_Classic",
                    yrs=c(80, 120),
                    params=c("PExtinct", "Nextant", "Het", "Nalleles"),
                    save2disk=FALSE)

  expect_length(dot , 4)
  expect_is(dot,"list")
  expect_is(dot[[1]], c("gg", "ggplot"))
})

test_that("test line_plot", {
    data(pac.clas)
    lineplot.st.classic <- line_plot_year(data=pac.clas, project="Pacioni_et_al",
                                          scenario="ST_Classic",
                                          params=c("Nextant"),
                                          save2disk=FALSE)
    expect_length(lineplot.st.classic , 1)
    expect_is(lineplot.st.classic,"list")
    expect_is(lineplot.st.classic[[1]], c("gg", "ggplot"))
})

test_that("test line_plot_mid", {
    data(pac.clas)
    lineMidPlot.st.classic <- line_plot_year_mid(data=pac.clas,
                                                 project="Pacioni_et_al",
                                                 scenario="ST_Classic",
                                                 yrmid=50,
                                                 params=c("Nextant"),
                                                 save2disk=FALSE)
    expect_length(lineMidPlot.st.classic , 1)
    expect_is(lineMidPlot.st.classic,"list")
    expect_is(lineMidPlot.st.classic[[1]], c("gg", "ggplot"))
})
