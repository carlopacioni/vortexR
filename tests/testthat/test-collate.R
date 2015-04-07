library(vortexR)
context("test collate")


test_that("test collate_dat", {
  # dir
  camp.dir <- system.file("extdata", "campbell", package="vortexR")
	pac.dir <- system.file("extdata", "pacioni", package="vortexR")

  # load data
  data(sta.main)
	starling <- collate_dat("Starlingv3PopBased", 10000, dir_in=camp.dir,
                        save2disk=FALSE, verbose=FALSE)

  data(pac.clas)
	woylie.st.classic <- collate_dat("Pacioni_et_al", 3,
	                                 scenario = "ST_Classic", dir_in = pac.dir,
	                                 save2disk=FALSE, verbose=FALSE)

  expect_equal(starling , sta.main)
  expect_equal(woylie.st.classic, pac.clas)
})

test_that("test collate_run", {
  # dir
  pac.dir <- system.file("extdata", "pacioni", package="vortexR")

  # load data
  data(pac.run.lhs)
  woylie.st.run <- collate_run("Pacioni_et_al", scenario = "ST_LHS",
                               dir_in = pac.dir, save2disk=FALSE, verbose=FALSE)

  expect_equal(woylie.st.run, pac.run.lhs)
})

test_that("test collate_yr", {
  # dir
  pac.dir <- system.file("extdata", "pacioni", package="vortexR")

  # load data
  data(pac.yr)
  woylie.yr <- collate_yr("Pacioni_et_al", scenario = "ST_Classic",
                               dir_in = pac.dir, save2disk=FALSE, verbose=FALSE)

  expect_equal(woylie.yr, pac.yr)
})
