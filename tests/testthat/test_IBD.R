cat("test_IBD.R:\n")

testthat::test_that("ibd() output", { 
  locinfile <- genepopExample('w2.txt')
  outfile <- ibd(locinfile,"w2.txt.ISO", geographicScale = "Log", statistic="e")
 
  if (sessionInfo()[["R.version"]][["arch"]]=="i386") { # i386 vs x64
    testthat::expect_equal(readLines(outfile)[229], "0.0151487 [ -0.0143973 , 0.0374648 ]") 
  } else testthat::expect_equal(readLines(outfile)[229], "0.0151487 [ -0.0143977 , 0.0374648 ]")
  
  locinfile <- genepopExample('PEL1600withCoord.txt')
  outfile <- ibd(locinfile,"PEL1600withCoord.ISO", statistic = "SingleGeneDiv",
                geographicScale = "1D")
  nums <- as.numeric(unlist(strsplit(readLines(outfile)[59], "[^0-9e.-]+")))
  testthat::expect_equal(nums, c(2.92606e-06, 6.28251e-07, 6.25044e-06))
  
})

testthat::test_that("Whether BCa differs from BC as expected", { 
  outfile <- ibd(genepopExample('w2.txt'), 'w2.txt.ISO', 
                 geographicScale = 'Log', 
                 statistic='e', 
                 bootstrapMethod = "BCa", 
                 bootstrapNsim = 9)
  numsBCa <- as.numeric(unlist(strsplit(readLines(outfile)[c(229,243)], "[^0-9e.-]+")))
  
  outfile <- ibd(genepopExample('w2.txt'), 'w2.txt.ISO', 
                 geographicScale = 'Log', 
                 statistic='e', 
                 bootstrapMethod = "BC", 
                 bootstrapNsim = 9)
  numsBC <- as.numeric(unlist(strsplit(readLines(outfile)[c(229,243)], "[^0-9e.-]+")))
  # Elements 1 and 4 are the point estimate, they must be identical by each method
  # Other elements are CI bounds, small differences appear for increasing bootstrapNsim. Two are already visible for bootstrapNsim=9
  # For for bootstrapNsim=99 they are all different (but that would also be the case 
  #  if ABC were run instead of the BC/BCa, so the presnet check is stronger)
  testthat::expect_true(  all((numsBCa==numsBC)==c(TRUE,TRUE,FALSE,TRUE,FALSE,TRUE)))
})




