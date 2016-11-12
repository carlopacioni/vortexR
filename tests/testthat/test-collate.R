library(vortexR)
context("test collate")


test_that("test collate_dat", {
  # dir
  camp.dir <- system.file("extdata", "campbell", package="vortexR")
  example_file <- system.file("extdata", "pacioni.zip", package="vortexR")
  pac.dir <- tempdir()
  unzip(example_file, exdir = pac.dir)


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
  unlink(pac.dir, recursive = TRUE)

})

test_that("test collate_run", {
  # dir
  example_file <- system.file("extdata", "pacioni.zip", package="vortexR")
  pac.dir <- tempdir()
  unzip(example_file, exdir = pac.dir)


  # load data
  data(pac.run.lhs)
  woylie.st.run <- collate_run("Pacioni_et_al", scenario = "ST_LHS",
                               dir_in = pac.dir, save2disk=FALSE, verbose=FALSE)

  expect_equal(woylie.st.run, pac.run.lhs)
  unlink(pac.dir, recursive = TRUE)

})

test_that("test collate_yr", {
  # dir
    example_file <- system.file("extdata", "pacioni.zip", package="vortexR")
    pac.dir <- tempdir()
    unzip(example_file, exdir = pac.dir)

  # load data
  data(pac.yr)
  woylie.yr <- collate_yr("Pacioni_et_al", scenario = "ST_Classic",
                               dir_in = pac.dir, save2disk=FALSE, verbose=FALSE)

  expect_equal(woylie.yr, pac.yr)
  unlink(pac.dir, recursive = TRUE)

})
