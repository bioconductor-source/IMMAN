## ------------------------------------------------------------------------

library(IMMAN)


## ------------------------------------------------------------------------

data(Celegance)
data(FruitFly)


## ------------------------------------------------------------------------

ProteinLists = list(as.character(Celegance$V1), as.character(FruitFly$V1))

List1_Species_ID = 6239  # taxonomy ID Celegance
List2_Species_ID = 7227  # taxonomy ID FruitFly

Species_IDs  = c(List1_Species_ID, List2_Species_ID)


## ------------------------------------------------------------------------

identityU = 30
substitutionMatrix = "BLOSUM62"
gapOpening = -8
gapExtension = -8
NetworkShrinkage = FALSE
coverage = 1
BestHit = TRUE
score_threshold = 400
STRINGversion="10"


## ------------------------------------------------------------------------

output = IMMAN(ProteinLists, fileNames=NULL, Species_IDs,
              identityU, substitutionMatrix,
              gapOpening, gapExtension, BestHit,
              coverage, NetworkShrinkage,
              score_threshold, STRINGversion,
              InputDirectory = getwd())


## ------------------------------------------------------------------------
output$IPNEdges
output$IPNNodes
output$Networks
output$Networks[[1]]
output$maps
output$maps[[2]]

