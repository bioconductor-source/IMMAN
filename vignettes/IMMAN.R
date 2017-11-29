## ------------------------------------------------------------------------

library(IMMAN)


## ------------------------------------------------------------------------

data(H.sapiens)
data(R.norvegicus)
data(Celegance)
data(FruitFly)


## ------------------------------------------------------------------------

ProteinLists = list(as.character(H.sapiens$V1), as.character(R.norvegicus$V1),
 as.character(Celegance$V1), as.character(FruitFly$V1))

List1_Species_ID = 9606  # taxonomy ID List1 Homo sapiens
List2_Species_ID = 10116 # taxonomy ID List2 Rat
List3_Species_ID = 6239  # taxonomy ID List3 Celegans
List4_Species_ID = 7227  # taxonomy ID List4 FruitFly

Species_IDs  = c(List1_Species_ID, List2_Species_ID, List3_Species_ID, List4_Species_ID)


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


