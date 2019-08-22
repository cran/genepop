cat("tests_HW.R:\n")

locinfile <- genepopExample('sample.txt')
outfile <- test_HW(locinfile,"deficit","sample.txt.D")
testthat::expect_equal(readLines(outfile)[133],
             "ADH-5       1.0000  0.0000  -0.3333 -0.3750  19903 switches")
outfile <- test_HW(locinfile,"global deficit")
testthat::expect_equal(readLines(outfile)[48], " 0.3329  0.0048 24797.00")

outfile <- nulls(locinfile, "sample.txt.NUL") ## FR->FR no testthat check !

# Example in Guo & Thompson 1992 Table 5
locinfile <- genepopExample('Rhesus.txt')
outfile <- HWtable_analysis(locinfile,which="Proba",batches = 1000,iterations = 1000)
testthat::expect_equal(readLines(outfile)[21],"P-value=0.684292; S.E=0.00830887 (241565 switches)")

