#' Extract vital rates from time series of Comadre/Compadre matrices
#'
#' @param cadre A Compadre or Comadre object, subset to contain two or more "Individual" matrices for a single population
#'
#' @return A list with two elements
#' \describe{
#'  \item{vrmat}{A matrix with a row for each vital rate and a column for each matrix from cadre}
  #' \item{matstruct}{A list with lots of structural information about the vital rate decomposition of the matrix, which can be used by \code{\link{vrtomat}} to reconstruct the matrix.}
#' }
#' 
#' @export
#'
extractVRs <- function(cadre) {
  require(popbio)
  # Get lists of A, F, U matrices
  matlists <- make_mat_lists(cadre)
  nmat <- nrow(cadre$metadata)
  matrank <- nrow(cadre$mat[[1]]$matF)

  # Find structural nonzeros from mean matrices
  Fmean <- mean(matlists$F)
  Fkeep <- Fmean > 0
  Fnames <- make_mat_element_names("F", matrank)[Fkeep]
  Findices <- make_index_list(Fnames)

  Umean <- mean(matlists$U)
  Skeep <- colSums(Umean) > 0
  Snames <- paste("S", seq(1, matrank)[Skeep], sep = ".")
  Sindices <- make_index_list(Snames)

  gkeep <- Umean > 0
  Gmultnames <- make_mat_element_names("g", matrank)[gkeep]
  Gmultindices <- make_index_list(Gmultnames)
  gcount <- plyr::count(Gmultindices, vars = "j")
  growStages <- subset(gcount, freq > 1, select = j)
  keepRows <- Gmultindices$j %in% growStages$j
  Gmultnames <- Gmultnames[keepRows]
  Gmultindices <- Gmultindices[keepRows,]
  Gindices <- ddply(Gmultindices, "j", subset, i > min(i))
  Gindices$element <- sub("g", "G", Gindices$element)
  Gnames <- Gindices$element
  Gbase <- apply(Umean > 0, 2, function(x) match(TRUE, x))

  # Classify the growth types.
  # Note that this classification is incomplete and will likely result in errors
  #   if a stage is skipped
  gtype <- vector("character", matrank)
  for (jj in 1:matrank) {
    if (gcount$freq[jj] == 1) {
      if (Umean[jj, jj] > 0) {
        gtype[jj] <- "stasis"
      } else if (Umean[jj+1, jj] > 0) {
        gtype[jj] <- "age"
      }
    } else if (gcount$freq[jj] == 2) {
      gtype[jj] <- "stage"
      if (length(setdiff(subset(Gmultindices, j == jj)$i, c(jj, jj+1))) > 0) {
        gtype[jj] <- "size"
      }
    } else if (gcount$freq[jj] > 2) {
      gtype[jj] <- "size"
    }
  }

  matstruct <- list(matrank = matrank, Findices = Findices, Sindices = Sindices,
                    gtype = gtype, Gindices = Gindices, Gbase = Gbase,
                    Amean = mean(matlists$A))

  # Set up the matrix to store the vrs in
  vrnames <- c(Fnames, Snames, Gnames)
  vecF <- 1:length(Fnames)
  vecS <- (length(Fnames) + 1):(length(Fnames) + length(Snames))
  vecG <- (length(Fnames) + length(Snames) + 1):length(vrnames)
  vrmat <- matrix(0, nrow = length(vrnames), ncol = nmat, dimnames = list(vrnames))

  for (id in 1:nmat) {
    mats <- cadre$mat[[id]]
    vrmat[vecF, id] <- mats$matF[Fkeep]
    svec <- apply(mats$matU, 2, sum)
    vrmat[vecS, id] <- svec[Skeep]
    gmat <- mats$matU / matrix(svec, matrank, matrank, byrow = TRUE)
    Gvr <- NULL
    for (jj in growStages$j) {
      if (gtype[jj] == "stage") {
        Gvr <- c(Gvr, gmat[jj+1, jj])
      } else if (gtype[jj] == "size") {
        Gloc <- numeric(gcount$freq[jj] - 1)
        gloc <- gmat[Umean[, jj] > 0, jj]
        Gloc[1] <- 1 - gmat[Gbase[jj], jj]
        if (length(Gloc) > 1) {
          for (ii in 2:length(Gloc)) {
            Gloc[ii] <- Gloc[ii-1] - gloc[ii]
          }
        }
        Gvr <- c(Gvr, Gloc)
      }
    }
    vrmat[vecG, id] <- Gvr
  }
  vrmat[vrmat < .Machine$double.eps] <- 0

  return(list(vrmat = vrmat, matstruct = matstruct))
}


make_mat_lists <- function(cadre) {
  require(plyr)
  mats <- cadre$mat
  matlist <- list(A = llply(mats, function(x) x$matA),
                  F = llply(mats, function(x) x$matF),
                  U = llply(mats, function(x) x$matU),
                  C = llply(mats, function(x) x$matC))
  return(matlist)
}

make_start_date <- function(cadre, id) {

}

make_mat_element_names <- function(prefix, matrank) {
  matrix(paste(prefix, ".",  rep(1:matrank, times = rep(matrank, matrank)),
               ".", 1:matrank, sep=''), nrow = matrank, byrow = TRUE)
}

make_index_list <- function(Names) {
  Names <- as.vector(Names)
  nn <- length(strsplit(Names[1], ".", fixed = TRUE)[[1]])
  indices <- matrix(unlist(strsplit(Names, ".", fixed = TRUE)),
                    ncol = nn, byrow = TRUE)[,-1]
  indices <- data.frame(indices, stringsAsFactors = FALSE)
  indices <- data.frame(apply(indices,2,as.numeric))
  if (ncol(indices) == 1) {
    names(indices) <- "j"
  } else {
    names(indices) <- c("i","j")
  }
  indices <- cbind(data.frame(element = Names), indices)
  indices
}