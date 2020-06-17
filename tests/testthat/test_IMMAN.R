library(IMMAN)
library(testthat)

test_that("IMMAN reconstruct the interlog protein network", {

 result <- IMMAN(ProteinLists=list(as.character(Celegance[1:25,]), as.character(FruitFly[1:25,])),
                 fileNames=NULL, Species_IDs  = c(6239, 7227),
        identityU = 30,
        substitutionMatrix = "BLOSUM62",
        gapOpening = -8,
        gapExtension = -8,
        NetworkShrinkage = FALSE,
        coverage = 1,
        BestHit = TRUE,
        score_threshold = 400,
        STRINGversion="11",
        InputDirectory = getwd())

 expect_that( result$IPNEdges, is_a("data.frame") )
 expect_that( result$IPNNodes, is_a("data.frame") )

})
