library(testthat)
library(IMMAN)

test_check("IMMAN")

data("H.sapiens")
data("R.norvegicus")
test_that("IMMAN reconstruct the interlog protein network", {
  IMMAN(ProteinLists=list(as.character(H.sapiens$V1), as.character(R.norvegicus$V1)), fileNames=NULL, Species_IDs,
        identityU = 30,
        substitutionMatrix = "BLOSUM62",
        gapOpening = -8,
        gapExtension = -8,
        NetworkShrinkage = FALSE,
        coverage = 1,
        BestHit = TRUE,
        score_threshold = 400,
        STRINGversion="10",
        InputDirectory = getwd())
})
