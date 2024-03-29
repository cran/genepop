cat("test_diff.R:\n")

testthat::test_that("test_diff() output", { 
  locinfile <- genepopExample('sample.txt')
  outfile <- test_diff(locinfile,outputFile="sample.txt.GE")
  testthat::expect_equal(readLines(outfile)[116],"mtDNA          0.2106   ")
  outfile <- test_diff(locinfile,pairs=TRUE,outputFile="sample.txt.GE2")
  testthat::expect_equal(readLines(outfile)[157],"Tertre Rotebo & last pop      14.52354   12    0.268532")
  
# Was: 
#     testthat::expect_equal(readLines(outfile)[158],"Bonneau 05    & last pop      48.07912   12    3.03e-006") # check use of scientific format.
# Now adjust to correct the variable printing of scientific notation
  nums <- as.numeric(tail(unlist(strsplit(readLines(outfile)[158], " +")), 3))
  refs <- c(48.07912, 12, 3.03e-006)
  testthat::expect_equal(nums, refs)
  
  #
  outfile <- test_diff(locinfile,genic=FALSE,outputFile="sample.txt.G")
  nums <- as.numeric(unlist(strsplit(readLines(outfile)[112], "[^0-9e.-]+"))[c(3,4,7)])
  testthat::expect_equal(nums,c(81.0679,10,3.09949e-013))
  # : avoiding potential e-11/e-011 mismatch
  #
  outfile <- test_diff(locinfile,genic=FALSE,pairs=TRUE,outputFile="sample.txt.2G2")
  testthat::expect_equal(readLines(outfile)[146],"Bonneau 05    & last pop      32.50364   10    0.000330")
})
