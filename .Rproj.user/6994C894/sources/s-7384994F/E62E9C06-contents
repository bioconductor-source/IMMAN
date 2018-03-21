#setwd("C:\\Users\\user\\Desktop\\PEIMAN R")
##data = read.table("ptm.txt" , header = TRUE, sep ="\t")

attach("1.rda")
names(data) <- c("AC" , "OS", "PTM")

x <- read.table("list.txt" , header = FALSE)
UserList <- x[,1]
UserList <- as.character(UserList)

OS <- "Homo sapiens (Human)."
OSindx <- which( data$"OS" == OS)
data.OS <- data[OSindx,]

indx <- which(UserList %in% data.OS$"AC" )

ValidList   <- UserList[indx]

P = unique(data.OS$"AC")
N <- length(P)
n <- length(ValidList)


indx <- which(data.OS$"AC" %in% ValidList)

u <- as.character(data.OS[indx,]$"PTM")
v <- as.character(data.OS[indx,]$"AC")
Rep <- data.frame(v,u)
names(Rep) <- c("AC" , "PTM")

List <- unique(Rep$"PTM")
List <- as.character(List)
freq = c()
for( i in 1:length(List)){
freq = c ( freq , length( which( Rep$"PTM" %in% List[i]) ) )
}
user <- data.frame(List , freq)

List <- unique(data.OS$"PTM")
List <- as.character(List)
freq = c()
for( i in 1:length(List)){
  freq = c ( freq , length(which( data.OS$"PTM" %in% List[i]  )) )
}
uniprot <- data.frame(List , freq)


indx <- which(uniprot$"List" %in% user$"List" , arr.ind = TRUE )
uniprot <- uniprot[indx,]

user = user[ order(user$"List") , ]
uniprot = uniprot[ order(uniprot$"List") , ] 

out <- data.frame( PTM = uniprot$"List" , FreqinList = user$"freq" ,
                  FreqinUniProt = uniprot$"freq", Sample = n, Population = N)

x = out$FreqinList
m = out$FreqinUniProt
nn = N - m
k = n
out$pvalue <- 1 - phyper(x , m , nn , k)
out$Corrected.Pvalue <- p.adjust(out$pvalue, method = "fdr")
out



