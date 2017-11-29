#' @title Interlog protein network reconstruction by Mapping and Mining ANalysis
#'
#' @description
#' A function for reconstructing Interlog Protein Network (IPN) integrated
#' from  Protein-protein Interaction Networks (PPIN) from different species. Users can overlay different
#' PPINs to mine conserved common network between diverse species. It helps to retrieve
#' IPN with different degrees of conservation to have better protein function
#' prediction and PPIN analysis.
#' @param ProteinLists a list in which each element contains protein names of a species as a character vector.
#' If it was NULL then the protein lists file name should be addressed in fileNames parameter.
#' @param fileNames a character vector, containing names of text files containing protein list for each species.
#' The protein list of each species must be in a column without header and rownames in seperate ".txt" files.
#' The ProteinLists argument should be include at least two text file names addressing the protein list of each species
#' which are in UniProt accession IDs format.
#' @param Species_IDs a numeric vector; taxonomy ID for each organism which are provided in fileNames
#' @param identityU numeric; value for selecting proteins whose alignment score is greater or equal than identityU
#' @param substitutionMatrix a scoring substitution matrix to be used for alignment setting.
#' @param gapOpening numeric; indicating the cost for opening a gap in the alignment
#' @param gapExtension The incremental cost incurred along the length of the gap in the alignment
#' @param BestHit logical; if TRUE describes a pair protein sequence among two different species
#' which is the reciprocal best hit in sequence similarity analysis, whilst,
#' if it is FALSE, indicates a nonreciprocal best hit
#' @param coverage Number of connected proteins pairs in each Ortholog Protein Set (OPS) pair
#' (termed as “coverage”) to
#' reconstruct an edge of OPS pair in the IPN (Interlog Protein Network)
#' @param NetworkShrinkage logical; if TRUE OPSs that are similar to each other would be merged.
#' @param score_threshold numeric; STRINGdb score for protein protein interaction (PPI) selection in STRING database
#' @param STRINGversion character; indicating which version of STRING database should program
#' search in for the score of PPIs.
#' @param InputDirectory By default is getwd(). You can set this parameter to indicate where
#' the downloaded file from STRING should be saved.
#' @seealso \code{\link[Biostrings]{PairwiseAlignments}}
#' @return
#' a list contaning four elements:
#'
#' IPNEdges : data.frame; Edges of resulted interlog protein network.
#'
#' IPNNodes : data.frame; Nodes of resulted interlog protein network.
#' Each node represents an OPS which is a set of ortholog proteins.
#'
#' Network : list; Retrived PPINs of each input species.
#'
#' maps : list; It includes data.frames indicating STRING_id data base matched to their corresponding
#' UNIPROT_AC. The number of data.frames is according to the the number of species.
#'
#' IPN : an igraph object representing the interlog protein network.
#'
#' @author Payman Nickchi, Abdollah Safari, Minoo Ashtiani, Mohieddin Jafari
#' @examples
#' data(H.sapiens)
#' data(R.norvegicus)
#'
#' ProteinLists = list(as.character(H.sapiens$V1), as.character(R.norvegicus$V1))
#'
#' List1_Species_ID = 9606  # taxonomy ID List1 Homo sapiens
#' List2_Species_ID = 10116 # taxonomy ID List2 Rat
#'
#' Species_IDs  = c(List1_Species_ID, List2_Species_ID)
#'
#'identityU = 30
#'substitutionMatrix = "BLOSUM62"
#'gapOpening = -8
#'gapExtension = -8
#'NetworkShrinkage = FALSE
#'coverage = 1
#'BestHit = TRUE
#'score_threshold = 400
#'STRINGversion="10"
#'
#'# Run the IMMAN function for the parameters
#' output = IMMAN(ProteinLists, fileNames=NULL, Species_IDs,
#'               identityU, substitutionMatrix,
#'               gapOpening, gapExtension, BestHit,
#'               coverage, NetworkShrinkage,
#'               score_threshold, STRINGversion,
#'               InputDirectory = getwd())
#'
#'output$IPNEdges
#'output$IPNNodes
#'output$Networks
#'output$Networks[[1]]
#'output$maps
#'output$maps[[2]]
#'
#' @export
#' @importFrom Biostrings pairwiseAlignment
#' @importFrom Biostrings pid
#' @importFrom igraph graph_from_data_frame
#' @importFrom graphics plot
#' @importFrom STRINGdb STRINGdb
#' @importFrom utils read.csv
#' @importFrom igraph layout_in_circle
#' @importFrom seqinr read.fasta

IMMAN <- function(ProteinLists, fileNames = NULL, Species_IDs,
                 identityU, substitutionMatrix,
                 gapOpening, gapExtension, BestHit,
                 coverage, NetworkShrinkage,
                 score_threshold, STRINGversion,
                 InputDirectory = getwd()){

  align <- function(pattern, subject, type){

    p = pairwiseAlignment(pattern, subject, scoreOnly = FALSE,
                          substitutionMatrix = substitutionMatrix,
                          gapOpening = gapOpening,
                          gapExtension = gapExtension)

    return( pid(p  , type = "PID1") )
  }

  combination3 <- function(L1, L2, L3) {

    # res1_tem <- unlist(Lists[1])
    # res2_tem <- unlist(Lists[2])
    # res3_tem <- unlist(Lists[3])

    res1_tem <- L1
    res2_tem <- L2
    res3_tem <- L3

    # combination
    x.names <- character(0)
    y.names <- character(0)
    z.names <- character(0)
    for (i in 1 : nrow(res1_tem)) {

      if (sum(res1_tem[i, ]) == 1) {
        vec_tem <- res2_tem[i, ] * res3_tem[(res1_tem[i,] == 1), ]
        if (sum(vec_tem) >= 1) {
          x.names <- c(x.names, rep(rownames(res1_tem)[i], sum(vec_tem)))
          y.names <- c(y.names, rep(rownames(res3_tem)[(res1_tem[i,] == 1)], sum(vec_tem)))
          z.names <- c(z.names, colnames(res3_tem)[vec_tem == 1])
        }
      }

      if (sum(res1_tem[i, ]) > 1) {
        list_tem <- apply(cbind(res3_tem[(res1_tem[i,] == 1), ], which(res1_tem[i,] == 1)), 1,
                           function(x) {
                             vec_tem <- res2_tem[i, ] * x[-length(x)]
                             if (sum(vec_tem) >= 1) {
                               x.names_tem <- rep(rownames(res1_tem)[i], sum(vec_tem))
                               y.names_tem <- rep(rownames(res3_tem)[x[length(x)]], sum(vec_tem))
                               z.names_tem <- colnames(res3_tem)[vec_tem == 1]
                               list(xx = x.names_tem,
                                    yy = y.names_tem,
                                    zz = z.names_tem)
                             }
                           })

        x.names <- c(x.names, unlist(lapply(list , function(x) unlist(x[1]))))
        y.names <- c(y.names, unlist(lapply(list_tem, function(x) unlist(x[2]))))
        z.names <- c(z.names, unlist(lapply(list_tem, function(x) unlist(x[3]))))
      }
    }

    return(list(x.names, y.names, z.names))
  }

  if(is.null(ProteinLists)){
  if(is.element(FALSE,file.exists(fileNames)) %in% FALSE) {
    files<-list()
    ProteinLists<-list()
    for (i in 1:length(fileNames)) {
      files[[i]] <- read.csv(fileNames[i],
                             header = F)
      ProteinLists[[i]]<-as.character(as.data.frame(files[[i]])$V1)
    }
  }
  else stop("Error: The file name does not exist in the diretory path")
}

  Species_IDs<-sapply(Species_IDs,list)

  list_num <- length(ProteinLists)

  for (i in 1 : list_num) {
    if( is.character(ProteinLists[[i]]) == "FALSE" ) {
      stop(paste("Error: ProteinLists", i," should be a vector type of character", sep = ""))
    }
  }

  if( !(coverage %in% 1:4) ){
    stop("Error: Coverage should be between 1 up to 4")
  }


  ## Print status
  message("Step 1/4:Downloading amino acid sequences...")
  PS_list <- list()
  for (i in 1 : list_num) {
    message(paste("Downloading amino acid sequences of List", i, sep = ""))
    PS_list <- c(PS_list, list(apply(data.frame(Protein = ProteinLists[[i]]), 1,
                                     function (x) {
                                       as.character(read.fasta(file =
                                                                 paste("http://www.uniprot.org/uniprot/",
                                                                       x,".fasta",sep=""),
                                                               seqtype ="AA", as.string = TRUE,
                                                               set.attributes = FALSE))
                                     })))
  }

  names(PS_list) <- c(paste("PS", seq(1 : list_num), sep = ""))


  message("Step 2/4: Alignment...")

  tem_list <- list()
  res_list <- list()
  for (i in 1 : (list_num - 1)) {
    for (j in (i + 1) : list_num) {
      message(paste("Align List", i," with List", j, sep = ""))
      unbinres = t(apply(as.matrix(PS_list[[i]], ncol = 1), 1, function (x) {
        apply(as.matrix(PS_list[[j]], ncol = 1), MARGIN = 1, FUN=align, x)
      }))
      res = matrix(as.numeric(unbinres > identityU), nrow(unbinres), ncol(unbinres))
      rownames(res) <- ProteinLists[[i]]
      colnames(res) <- ProteinLists[[j]]
      indx = res == 1
      tem_list = c(tem_list, list(indx * unbinres))
      res[res == 1] <- 0
      res_list <- c(res_list, list(res))
    }
  }

  pair_num <- list_num * (list_num - 1) / 2
  names(tem_list) <- paste("tem", seq(1 : pair_num), sep = "")
  names(res_list) <- paste("res", seq(1 : pair_num), sep = "")

  # tem_list_backup <- tem_list
  # res_list_backup <- res_list
  #
  # tem_list <- tem_list_backup
  # res_list  <- res_list_backup

  if (BestHit == TRUE) {
    for (i in 1 : pair_num) {
      tem_list[[i]] <- apply(tem_list[[i]], 1, function(x) {
        ind_tem <- which(x == max(x))
        if (length(ind_tem) > 1) colnames(tem_list[[i]])[ind_tem][apply(tem_list[[i]][, ind_tem], 2, max) == x[ind_tem]]
        else colnames(tem_list[[i]])[ind_tem][max(tem_list[[i]][, ind_tem]) == x[ind_tem]]
      })

      for (j in 1 : length(tem_list[[i]])) {
        if (length(unlist(tem_list[[i]][j])) > 0) res_list[[i]][names(tem_list[[i]][j]), unlist(tem_list[[i]][j])] <- 1
      }
    }
  }
  if( BestHit == FALSE){
    for (i in 1 : pair_num) {
      tem_list[[i]] <- apply(tem[[i]], 1, function(x) {
        colnames(tem[[i]])[x != 0]
      })

      for (j in 1 : length(tem_list[[i]])) {
        if (length(unlist(tem_list[[i]][j])) > 0) res_list[[i]][names(tem_list[[i]][j]), unlist(tem_list[[i]][j])] <- 1
      }
    }
  }

  message("Step 3/4: Detection in STRING...")

  string_db_list <- list()
  map_list <- list()
  for (i in 1 : list_num) {
    string_db_list <- c(string_db_list, list(STRINGdb$new(version = STRINGversion,
                                                          species= Species_IDs[[i]],
                                                          score_threshold = score_threshold,
                                                          input_directory = getwd())))
  }
  if (list_num == 4) {
    # 1-2
    x_tem1 <- rep(names(tem_list[[1]]), unlist(lapply(tem_list[[1]], length)))
    y_tem1 <- unlist(tem_list[[1]])
    # 1-3
    x_tem2 <- rep(names(tem_list[[2]]), unlist(lapply(tem_list[[2]], length)))
    z_tem1 <- unlist(tem_list[[2]])
    # 1-4
    x_tem3 <- rep(names(tem_list[[3]]), unlist(lapply(tem_list[[3]], length)))
    w_tem1 <- unlist(tem_list[[3]])
    # 2-3
    y_tem2 <- rep(names(tem_list[[4]]), unlist(lapply(tem_list[[4]], length)))
    z_tem2 <- unlist(tem_list[[4]])
    # 2-4
    y_tem3 <- rep(names(tem_list[[5]]), unlist(lapply(tem_list[[5]], length)))
    w_tem2 <- unlist(tem_list[[5]])
    # 3-4
    z_tem3 <- rep(names(tem_list[[6]]), unlist(lapply(tem_list[[6]], length)))
    w_tem3 <- unlist(tem_list[[6]])

    # x-y-z
    list.names1 <- combination3(res_list[[1]], res_list[[2]], res_list[[4]])
    x.names1 <- unlist(list.names1[1])
    y.names1 <- unlist(list.names1[2])
    z.names1 <- unlist(list.names1[3])
    mat.xyz1 <- cbind(x.names1, y.names1, z.names1)

    # x-y-w
    list.names2 <- combination3(res_list[[1]], res_list[[3]], res_list[[5]])
    x.names2 <- unlist(list.names2[1])
    y.names2 <- unlist(list.names2[2])
    w.names2 <- unlist(list.names2[3])
    mat.xyw1 <- cbind(x.names2, y.names2, w.names2)

    x.inters1 <- intersect(unique(x.names1), unique(x.names2))
    y.inters1 <- lapply(unique(x.inters1), function(x) {
      y.inters1 <- intersect(unique(y.names1[x.names1 == x]), unique(y.names2[x.names2 == x])) })

    mat.xyz2 <- matrix(NA, ncol = 3, nrow = 1)
    mat.xyw2 <- matrix(NA, ncol = 3, nrow = 1)
    for (i in 1 : length(x.inters1)) {
      mat.xyz2 <- rbind(mat.xyz2, mat.xyz1[mat.xyz1[ , 2] %in% unlist(y.inters1[i]), ])
      mat.xyw2 <- rbind(mat.xyw2, mat.xyw1[(mat.xyw1[, 1] == x.inters1[i]) & (mat.xyw1[ , 2] %in% unlist(y.inters1[i])), ])
    }
    mat.xyz2 <- mat.xyz2[-1, ]
    mat.xyw2 <- mat.xyw2[-1, ]

    x.inters2 <- intersect(unique(mat.xyz2[,1]), unique(unique(mat.xyw2[,1])))
    y.inters2 <- lapply(unique(x.inters2), function(x) {
      y.inters2 <- intersect(unique(mat.xyz2[mat.xyz2[,1] == x, 2]), unique(mat.xyw2[mat.xyw2[, 1] == x, 2])) })

    # Final lists
    mat.xyzw <- matrix(NA, nrow = 1, ncol = 4)
    for (i in 1 : length(y.inters2)) {
      for (j in 1 : length(y.inters2[[i]])) {
        z_tem <- mat.xyz2[(mat.xyz2[,1] == x.inters2[i]) & (mat.xyz2[,2] == y.inters2[[i]][j]), 3]
        w_tem <- mat.xyw2[(mat.xyw2[,1] == x.inters2[i]) & (mat.xyw2[,2] == y.inters2[[i]][j]), 3]
        res6_tem.sub <- matrix(res_list[[6]][z_tem,w_tem], ncol = length(w_tem), nrow = length(z_tem), T)
        sum_tem <- sum(res6_tem.sub)
        if (sum_tem > 0) {
          tem.mat <- matrix(NA, nrow = 1, ncol = 2)
          for (k in 1 : length(z_tem)) {
            tem.mat <- rbind(tem.mat, t(rbind(rep(z_tem[k], sum(res6_tem.sub[k,])),
                                                w_tem[res6_tem.sub[k,] == 1])))
          }
          tem.mat <- matrix(tem.mat[-1, ], ncol = 2, nrow = nrow(tem.mat) - 1)
          mat.xyzw <- rbind(mat.xyzw,
                            cbind(matrix(c(rep(x.inters2[i], sum_tem),
                                           rep(y.inters1[[i]][j], sum_tem)), ncol = 2, nrow = sum_tem),
                                  tem.mat))
        }
      }
    }
    mat.xyzw <- mat.xyzw[-1, ]


    message("Detecting List1 in STRING")
    map1 = string_db_list[[1]]$map( data.frame(UNIPROT_AC = unique(mat.xyzw[, 1])) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)
    if( nrow(map1) == 0 ) {
      print(ProteinLists[[1]])
      stop("Error: none of the proteins in list1 mapped to STRING ID")
    }

    message("Detecting List2 in STRING")
    map2 = string_db_list[[2]]$map(data.frame(UNIPROT_AC = unique(mat.xyzw[, 2])) ,
                                   "UNIPROT_AC" , removeUnmappedRows = TRUE)
    if( nrow(map2) == 0 ) {
      print(ProteinLists[[2]])
      stop("Error: none of the proteins in list2 mapped to STRING ID")
    }

    message("Detecting List3 in STRING")
    map3 = string_db_list[[3]]$map( data.frame(UNIPROT_AC = unique(mat.xyzw[, 3])) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)
    if( nrow(map3) == 0 ) {
      print(ProteinLists[[3]])
      stop("Error: none of the proteins in list3 mapped to STRING ID")
    }

    message("Detecting List4 in STRING")
    map4 = string_db_list[[4]]$map( data.frame(UNIPROT_AC = unique(mat.xyzw[, 4])) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)
    if( nrow(map4) == 0 ) {
      print(ProteinLists[[4]])
      stop("Error: none of the proteins in list4 mapped to STRING ID")
    }

    OPS = data.frame(node1 = mat.xyzw[, 1], node2 = mat.xyzw[, 2], node3 = mat.xyzw[, 3], node4 = mat.xyzw[, 4])
    OPS <- merge(OPS, map1, by.y = "UNIPROT_AC", by.x = "node1")
    colnames(OPS)[5] <- "STRING_id_1"
    OPS <- merge(OPS, map2, by.y = "UNIPROT_AC", by.x = "node2")
    colnames(OPS)[6] <- "STRING_id_2"
    OPS <- merge(OPS, map3, by.y = "UNIPROT_AC", by.x = "node3")
    colnames(OPS)[7] <- "STRING_id_3"
    OPS <- merge(OPS, map4, by.y = "UNIPROT_AC", by.x = "node4")
    colnames(OPS)[8] <- "STRING_id_4"

    OPS <- data.frame(node1 = OPS[, 5], node2 = OPS[, 6], node3 = OPS[, 7], node4 = OPS[, 8])

    OPSLabel = c()
    flag_tem <- TRUE
    if (nrow(OPS) > 10) {OPSLabel = c(OPSLabel, paste("OPS000", c(1 : 9), sep=""))
    } else {
      OPSLabel = c(OPSLabel, paste("OPS000", c(1 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 100) {OPSLabel = c(OPSLabel, paste("OPS00", c(10 : 99), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel, paste("OPS00", c(10 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 1000) {OPSLabel = c(OPSLabel, paste("OPS00", c(100 : 999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel, paste("OPS00", c(100 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 10000) {OPSLabel = c(OPSLabel, paste("OPS00", c(1000 : 9999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel, paste("OPS00", c(1000 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }

    OPS <- cbind(OPS, OPSLabel = OPSLabel)


    message("Step 4/4: Retrieving String Network...")

    message("Retrieving List1")
    Network1 = data.frame(from = string_db_list[[1]]$get_interactions(OPS$node1)$from,
                          to = string_db_list[[1]]$get_interactions(OPS$node1)$to)
    if(nrow(Network1) == 0){
      print(ProteinLists[[1]])
      stop("Error: No interaction was detected for ProteinLists1")
    }

    message("Retrieving List2")
    Network2 = data.frame(from = string_db_list[[2]]$get_interactions(OPS$node2)$from,
                          to = string_db_list[[2]]$get_interactions(OPS$node2)$to)
    if(nrow(Network2) == 0){
      print(ProteinLists[[2]])
      stop("Error: No interaction was detected for ProteinLists2")
    }

    message("Retrieving List3")
    Network3 = data.frame(from = string_db_list[[3]]$get_interactions(OPS$node3)$from,
                          to = string_db_list[[3]]$get_interactions(OPS$node3)$to)
    if(nrow(Network3) == 0){
      print(ProteinLists[[3]])
      stop("Error: No interaction was detected for ProteinLists3")
    }

    message("Retrieving List4")
    Network4 = data.frame(from = string_db_list[[4]]$get_interactions(OPS$node4)$from,
                          to = string_db_list[[4]]$get_interactions(OPS$node4)$to)
    if(nrow(Network4) == 0){
      print(ProteinLists[[4]])
      stop("Error: No interaction was detected for ProteinLists4")
    }

    message("Producing IPN...")

    node1 = c()
    node2 = c()
    l = nrow(OPS)
    # i = 1
    for (i in 1 : (l - 1)) {
      node_tem <- apply(OPS[c((i + 1) : l), ], 1, function(x){
        a = c(as.character(OPS[i,1]) , as.character(x[1]))
        b = c(as.character(OPS[i,2]) , as.character(x[2]))
        c = c(as.character(OPS[i,3]) , as.character(x[3]))
        d = c(as.character(OPS[i,4]) , as.character(x[4]))

        cond1 = ifelse(nrow(Network1[((Network1$from == a[1]) & (Network1$to) == a[2]), ]) != 0 |
                         nrow(Network1[((Network1$from == a[2]) & (Network1$to) == a[1]), ]) != 0, TRUE, FALSE)
        cond2 = ifelse(nrow(Network2[((Network2$from == b[1]) & (Network2$to) == b[2]), ]) != 0 |
                         nrow(Network2[((Network2$from == b[2]) & (Network2$to) == b[1]), ]) != 0, TRUE, FALSE)
        cond3 = ifelse(nrow(Network3[((Network3$from == c[1]) & (Network3$to) == c[2]), ]) != 0 |
                         nrow(Network3[((Network3$from == c[2]) & (Network3$to) == c[1]), ]) != 0, TRUE, FALSE)
        cond4 = ifelse(nrow(Network4[((Network4$from == d[1]) & (Network4$to) == d[2]), ]) != 0 |
                         nrow(Network4[((Network4$from == d[2]) & (Network4$to) == d[1]), ]) != 0, TRUE, FALSE)

        if (((cond1 + cond2 + cond3 + cond4 == 1) & (coverage == 1)) |
            ((cond1 + cond2 + cond3 + cond4 == 2) & (coverage == 2)) |
            ((cond1 + cond2 + cond3 + cond4 == 3) & (coverage == 3)) |
            ((cond1 + cond2 + cond3 + cond4 == 4) & (coverage == 4))) {
          return(c(as.character(OPS[i, 5]), as.character(x[5])))
        }

        if ((NetworkShrinkage == FALSE)) {
          t1 = as.character(OPS[i, 1]) == as.character(x[1])
          t2 = as.character(OPS[i, 2]) == as.character(x[2])
          t3 = as.character(OPS[i, 3]) == as.character(x[3])
          t4 = as.character(OPS[i, 4]) == as.character(x[4])

          mycond = cond1 + cond2 + cond3 + cond4
          TT = t1 + t2 + t3 + t4
          if (TT + mycond >= coverage){
            return(c(as.character(OPS[i, 5]), as.character(x[5])))
          }

        }
      })

      if (! is.null(node_tem)) {
        node1 <- c(node1, unlist(node_tem)[seq(1, length(unlist(node_tem)), 2)])
        node2 <- c(node2, unlist(node_tem)[seq(2, length(unlist(node_tem)), 2)])
      }
    }

    EdgeList = data.frame(node1 , node2)  #, node3, node4)

    map_list <- list(map1, map2, map3, map4)
    network_list <- list(Network1, Network2, Network3, Network4)
  }
  if (list_num == 3) {
    # 1-2
    x_tem1 <- rep(names(tem_list[[1]]), unlist(lapply(tem_list[[1]], length)))
    y_tem1 <- unlist(tem_list[[1]])
    # 1-3
    x_tem2 <- rep(names(tem_list[[2]]), unlist(lapply(tem_list[[2]], length)))
    z_tem1 <- unlist(tem_list[[2]])
    # 2-3
    y_tem2 <- rep(names(tem_list[[3]]), unlist(lapply(tem_list[[3]], length)))
    z_tem2 <- unlist(tem_list[[3]])

    # combination
    x.names <- character(0)
    y.names <- character(0)
    z.names <- character(0)
    for (i in 1 : nrow(res_list[[1]])) {

      if (sum(res_list[[1]][i, ]) == 1) {
        vec_tem <- res_list[[2]][i, ] * res_list[[3]][(res_list[[1]][i,] == 1), ]
        if (sum(vec_tem) >= 1) {
          x.names <- c(x.names, rep(rownames(res_list[[1]])[i], sum(vec_tem)))
          y.names <- c(y.names, rep(rownames(res_list[[3]])[(res_list[[1]][i,] == 1)], sum(vec_tem)))
          z.names <- c(z.names, colnames(res_list[[3]])[vec_tem == 1])
        }
      }

      if (sum(res_list[[1]][i, ]) > 1) {
        list_tem <- apply(cbind(res_list[[3]][(res_list[[1]][i,] == 1), ], which(res_list[[1]][i,] == 1)), 1,
                           function(x) {
                             vec_tem <- res_list[[2]][i, ] * x[-length(x)]
                             if (sum(vec_tem) >= 1) {
                               x.names_tem <- rep(rownames(res_list[[1]])[i], sum(vec_tem))
                               y.names_tem <- rep(rownames(res_list[[3]])[x[length(x)]], sum(vec_tem))
                               z.names_tem <- colnames(res_list[[3]])[vec_tem == 1]
                               list(xx = x.names_tem,
                                    yy = y.names_tem,
                                    zz = z.names_tem)
                             }
                           })

        x.names <- c(x.names, unlist(lapply(list_tem, function(x) unlist(x[1]))))
        y.names <- c(y.names, unlist(lapply(list_tem, function(x) unlist(x[2]))))
        z.names <- c(z.names, unlist(lapply(list_tem, function(x) unlist(x[3]))))
      }
    }


    message("Detecting List1 in STRING")

    map1 = string_db_list[[1]]$map( data.frame(UNIPROT_AC = unique(x.names)) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)


    if( nrow(map1) == 0 ) {
      stop("Error: none of the proteins in list1 mapped to STRING ID")
    }


    message("Detecting List2 in STRING")

    map2 = string_db_list[[2]]$map( data.frame(UNIPROT_AC = unique(y.names)) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)

    if( nrow(map2) == 0 ) {
      stop("Error: none of the proteins in list2 mapped to STRING ID")
    }


    message("Detecting List3 in STRING")

    map3 = string_db_list[[3]]$map(data.frame(UNIPROT_AC = unique(z.names)) ,
                                   "UNIPROT_AC" , removeUnmappedRows = TRUE)

    if( nrow(map3) == 0 ) {
      stop("Error: none of the proteins in list3 mapped to STRING ID")
    }

    OPS = data.frame(node1 = x.names, node2 = y.names, node3 = z.names)
    OPS <- merge(OPS, map1, by.y = "UNIPROT_AC", by.x = "node1")
    colnames(OPS)[4] <- "STRING_id_1"
    OPS <- merge(OPS, map2, by.y = "UNIPROT_AC", by.x = "node2")
    colnames(OPS)[5] <- "STRING_id_2"
    OPS <- merge(OPS, map3, by.y = "UNIPROT_AC", by.x = "node3")
    colnames(OPS)[6] <- "STRING_id_3"

    OPS <- data.frame(node1 = OPS[, 4], node2 = OPS[, 5], node3 = OPS[, 6])

    OPSLabel = c()
    flag_tem <- TRUE
    if (nrow(OPS) > 10) {OPSLabel = c(OPSLabel,paste("OPS000", c(1 : 9), sep=""))
    } else {
      OPSLabel = c(OPSLabel,paste("OPS000", c(1 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 100) {OPSLabel = c(OPSLabel,paste("OPS00", c(10 : 99), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel,paste("OPS00", c(10 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 1000) {OPSLabel = c(OPSLabel,paste("OPS00", c(100 : 999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel,paste("OPS00", c(100 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 10000) {OPSLabel = c(OPSLabel,paste("OPS00", c(1000 : 9999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel,paste("OPS00", c(1000 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }

    OPS <- cbind(OPS, OPSLabel = OPSLabel)


    message("Step 4/4: Retrieving String Network...")

    message("Retrieving List1")
    Network1 = data.frame(from = string_db_list[[1]]$get_interactions(OPS$node1)$from,
                          to = string_db_list[[1]]$get_interactions(OPS$node1)$to)
    if(nrow(Network1) == 0){
      print(ProteinLists[[1]])
      stop("Error: No interaction was detected for ProteinLists1")
    }

    message("Retrieving List2")
    Network2 = data.frame(from = string_db_list[[2]]$get_interactions(OPS$node2)$from,
                          to = string_db_list[[2]]$get_interactions(OPS$node2)$to)
    if(nrow(Network2) == 0){
      print(ProteinLists[[2]])
      stop("Error: No interaction was detected for ProteinLists2")
    }

    message("Retrieving List3")
    Network3 = data.frame(from = string_db_list[[3]]$get_interactions(OPS$node3)$from,
                          to = string_db_list[[3]]$get_interactions(OPS$node3)$to)
    if(nrow(Network3) == 0){
      print(ProteinLists[[3]])
      stop("Error: No interaction was detected for ProteinLists3")
    }


    message("Producing IPN...")

    node1 = c()
    node2 = c()
    l = nrow(OPS)

    # i = 1
    for (i in 1 : (l - 1)) {
      node_tem <- apply(OPS[c((i + 1) : l), ], 1, function(x){
        a = c(as.character(OPS[i,1]) , as.character(x[1]))
        b = c(as.character(OPS[i,2]) , as.character(x[2]))
        c = c(as.character(OPS[i,3]) , as.character(x[3]))

        cond1 = ifelse(nrow(Network1[((Network1$from == a[1]) & (Network1$to) == a[2]), ]) != 0 |
                         nrow(Network1[((Network1$from == a[2]) & (Network1$to) == a[1]), ]) != 0, TRUE, FALSE)
        cond2 = ifelse(nrow(Network2[((Network2$from == b[1]) & (Network2$to) == b[2]), ]) != 0 |
                         nrow(Network2[((Network2$from == b[2]) & (Network2$to) == b[1]), ]) != 0, TRUE, FALSE)
        cond3 = ifelse(nrow(Network3[((Network3$from == c[1]) & (Network3$to) == c[2]), ]) != 0 |
                         nrow(Network3[((Network3$from == c[2]) & (Network3$to) == c[1]), ]) != 0, TRUE, FALSE)

        if (((cond1 + cond2 + cond3 == 1) & (coverage == 1)) |
            ((cond1 + cond2 + cond3 == 2) & (coverage == 2)) |
            ((cond1 + cond2 + cond3 == 3) & (coverage ==3))) {
          return(c(as.character(OPS[i,4]), as.character(x[4])))
        }

        if ((NetworkShrinkage == FALSE)) {
          t1 = as.character(OPS[i,1]) == as.character(x[1])
          t2 = as.character(OPS[i,2]) == as.character(x[2])
          t3 = as.character(OPS[i,3]) == as.character(x[3])

          mycond = cond1 + cond2 + cond3
          TT = t1 + t2 + t3
          if (TT+mycond >= coverage){
            return(c(as.character(OPS[i,4]), as.character(x[4])))
          }

        }
      })

      if (! is.null(node_tem)) {
        node1 <- c(node1, unlist(node_tem)[seq(1, length(unlist(node_tem)), 2)])
        node2 <- c(node2, unlist(node_tem)[seq(2, length(unlist(node_tem)), 2)])
      }
    }

    EdgeList = data.frame(node1 , node2) #, node3)

    map_list <- list(map1, map2, map3)
    network_list <- list(Network1, Network2, Network3)
  }
  if (list_num == 2) {
    x <- rep(names(tem_list[[1]]), unlist(lapply(tem_list[[1]], length)))
    xperim = unlist(tem_list[[1]])

    message("Detecting List1 in STRING")
    map1 = string_db_list[[1]]$map(data.frame(UNIPROT_AC = unique(x)) ,
                                   "UNIPROT_AC" , removeUnmappedRows = TRUE)

    if ( nrow(map1) == 0 ) {
      stop("Error: none of the proteins in list1 mapped to STRING ID")
    }

    message("Detecting List2 in STRING")
    map2 = string_db_list[[2]]$map( data.frame(UNIPROT_AC = unique(xperim)) ,
                                    "UNIPROT_AC" , removeUnmappedRows = TRUE)

    if (nrow(map2) == 0) {
      stop("Error: none of the proteins in list2 mapped to STRING ID")
    }

    OPS = data.frame(node1 = x, node2 = xperim)
    OPS <- merge(OPS, map1, by.y = "UNIPROT_AC", by.x = "node1")
    colnames(OPS)[3] <- "STRING_id_1"
    OPS <- merge(OPS, map2, by.y = "UNIPROT_AC", by.x = "node2")
    colnames(OPS)[4] <- "STRING_id_2"

    OPS <- data.frame(node1 = OPS[, 3], node2 = OPS[, 4])

    OPSLabel = c()
    flag_tem <- TRUE
    if (nrow(OPS) > 10) {OPSLabel = c(OPSLabel, paste("OPS000", c(1 : 9), sep=""))
    } else {
      OPSLabel = c(OPSLabel,paste("OPS000", c(1 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 100) {OPSLabel = c(OPSLabel, paste("OPS00", c(10 : 99), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel, paste("OPS00", c(10 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 1000) {OPSLabel = c(OPSLabel, paste("OPS00", c(100 : 999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel,paste("OPS00", c(100 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }
    if (nrow(OPS) > 10000) {OPSLabel = c(OPSLabel,paste("OPS00", c(1000 : 9999), sep=""))
    } else {
      if (flag_tem) OPSLabel = c(OPSLabel,paste("OPS00", c(1000 : nrow(OPS)), sep=""))
      flag_tem <- FALSE
    }

    OPS <- cbind(OPS, OPSLabel = OPSLabel)

    message("Step 4/4: Retrieving String Network...")

    message("Retrieving List1")
    Network1 = data.frame(from = string_db_list[[1]]$get_interactions(OPS$node1)$from,
                          to = string_db_list[[1]]$get_interactions(OPS$node1)$to)

    if (nrow(Network1) == 0) {
      print(ProteinLists[[1]])
      stop("Error: No STRING network was detected for ProteinLists1")
    }

    message("Retrieving List2")
    Network2 = data.frame(from = string_db_list[[2]]$get_interactions(OPS$node2)$from,
                          to = string_db_list[[2]]$get_interactions(OPS$node2)$to)

    if (nrow(Network2) == 0) {
      print(ProteinLists[[2]])
      stop("Error: No STRING network was detected for ProteinLists2")
    }

    message("Producing IPN...")

    node1 = c()
    node2 = c()
    l = nrow(OPS)

    # i = 1
    for (i in 1 : (l - 1)) {
      node_tem <- apply(OPS[c((i + 1) : l), ], 1, function(x){
        a = c(as.character(OPS[i,1]) , as.character(x[1]))
        b = c(as.character(OPS[i,2]) , as.character(x[2]))

        cond1 = ifelse(nrow(Network1[((Network1$from == a[1]) & (Network1$to) == a[2]), ]) != 0 |
                         nrow(Network1[((Network1$from == a[2]) & (Network1$to) == a[1]), ]) != 0, TRUE, FALSE)
        cond2 = ifelse(nrow(Network2[((Network2$from == b[1]) & (Network2$to) == b[2]), ]) != 0 |
                         nrow(Network2[((Network2$from == b[2]) & (Network2$to) == b[1]), ]) != 0, TRUE, FALSE)

        if (((cond1+cond2 == 1) & (coverage == 1)) | ((cond1+cond2 == 2) & (coverage == 2))) {
          return(c(as.character(OPS[i,3]), as.character(x[3])))
        }

        if ((NetworkShrinkage == FALSE)) {
          t1 = as.character(OPS[i,1]) == as.character(x[1])
          t2 = as.character(OPS[i,2]) == as.character(x[2])

          mycond = cond1 + cond2
          TT = t1 + t2
          if (TT+mycond >= coverage){
            return(c(as.character(OPS[i,3]), as.character(x[3])))
          }

        }
      })

      if (! is.null(node_tem)) {
        node1 <- c(node1, unlist(node_tem)[seq(1, length(unlist(node_tem)), 2)])
        node2 <- c(node2, unlist(node_tem)[seq(2, length(unlist(node_tem)), 2)])
      }
    }

    EdgeList = data.frame(node1 , node2)

    map_list <- list(map1, map2)
    network_list <- list(Network1, Network2)
  }


  if (nrow(EdgeList) != 0) {
    IPN = graph_from_data_frame(d = EdgeList, directed = FALSE)
    reslist = list(IPNEdges = EdgeList , IPNNodes = OPS,
                   Networks = network_list,
                   maps = map_list)

    plot(IPN, layout = layout_in_circle(IPN))
    message("DONE!")
    return(reslist)
  }

  if (nrow(EdgeList) == 0) {
    reslist = list(IPNEdges = EdgeList , IPNNodes = OPS,
                   Networks = network_list,
                   maps = map_list)
    message("DONE! But EdgeList is empty!")
    return(reslist)
  }
}
