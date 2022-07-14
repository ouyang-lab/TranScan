#' Hi-C inter-chromosomal translocation detection via Scan clustering
#'
#' This function allows you to detect the inter-chromosomal translocation
#' from Hi-C data. It requires two intra-regions and the corresponding
#' inter-region as inputs. We suggest rounded log counts as inputs.
#' @param chr1 First chromosome
#' @param chr2 Second chromosome
#' @param chr1_2 Inter-chromosomal region of chr1 and chr2
#' @param outputpath The path to store all results
#' @param sizes Chrom sizes of chr1_2; vector of length 2
#' @param resolution Resolution of Hi-C data
#' @param quantitest Grid of testing, must be 2^k, default 256, meaning 256 * 256
#' @param alpha FDX, default 0.05
#' @param gammath Gamma threshold, default 0.1
#' @param number_of_bandwidth Number of bandwidths used, default 20
#' @param diag_remove Removed short-range interactions, must be k * resolution, default 2e+06
#' @param add_gap Added gap between two chromosomes, must be k * resolution, defult 2500000
#' @return A heatmap and a rejected region plot will be generated under the output folder. The matrix form of rejected region is saved as RR.txt. All intermediate objects will also be saved.
#' @examples
#' data(chr10)
#' data(chr20)
#' data(chr10_20)
#' p_sizes = "../hg38.sizes"
#' p_outdir = "output"
#' sizes = read.table(p_sizes, col.name=c("chrom", "size"))
#' chr10_l = subset(sizes, chrom=="chr10")$size
#' chr20_l = subset(sizes, chrom=="chr20")$size
#' size_dim = c(chr10_l, chr20_l)
#' TranScan(chr10, chr20, chr10_20, p_outdir, size_dim)
#' @export
TranScan <- function(chr1, chr2, chr1_2, outputpath, sizes, resolution = 500000,
                     quantitest = 256, alpha = 0.05, gammath = 0.1,
                     number_of_bandwidth = 20,
                     diag_remove = 2000000, add_gap = 2500000){
  library(ggplot2); library(dplyr); library(ggpubr)
  dir.create(outputpath, F)
  setwd(outputpath)
  
  bins1 <- (sizes[1] %/% resolution)+1
  bins2 <- (sizes[2] %/% resolution)+1
  chr1[, 1:2] <- chr1[, 1:2] %/% resolution + 1
  chr2[, 1:2] <- chr2[, 1:2] %/% resolution + 1
  chr1_2[, 1:2] <- chr1_2[, 1:2] %/% resolution + 1
  #bin_seg <- max(max(chr1[, 1:2]), max(chr1_2[, 1]))
  bin_seg <- bins1
  add_line <- add_gap / resolution
  chr2[, 1:2] <- chr2[, 1:2] + add_line + bin_seg
  chr1_2[, 2] <- chr1_2[, 2] + add_line + bin_seg
  sep_line <- bin_seg + add_line %/% 2 + 1
  dat <- rbind(chr1, chr2, chr1_2)
  colnames(dat) <- c("V1", "V2", "V3")
  dims <- c(bins1+add_line+bins2, bins1+add_line+bins2)
  dat <- dat %>% filter(V1 <= V2 - diag_remove / resolution)
  dat <- Matrix::sparseMatrix(dat[, 1], dat[, 2], x = dat[, 3], dims=dims, symmetric = T)
  dat <- reshape2::melt(as.matrix(dat))
  rm(list = c("chr1", "chr2", "chr1_2"))
  colnames(dat) <- c("V1", "V2", "V3")
  
  p1 <- ggplot(dat, aes(V1, -V2)) +
    geom_tile(aes(fill = V3)) +
    scale_fill_gradient2(low = "blue", mid = "white",
                         high = "red", midpoint = mean(dat$V3)) +
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"))+
    geom_hline(yintercept = -sep_line) + geom_vline(xintercept = sep_line)
  png("Heatmap.png")
  print(p1)
  dev.off()

  dati <- dat
  xmin <- min(dati[, 1])
  xmax <- max(dati[, 1])
  ymin <- min(dati[, 2])
  ymax <- max(dati[, 2])
  dati[,1] <- (dati[, 1] - min(dati[, 1])) / (max(dati[, 1]) - min(dati[, 1]))
  dati[,2] <- (dati[, 2] - min(dati[, 2])) / (max(dati[, 2]) - min(dati[, 2]))
  colnames(dati) <- c("V1", "V2", "V3")
  sep_ratio <- (sep_line - xmin) / (xmax - xmin)

  ################################
  # Prepare grid point positions #
  ################################
  Repeat_it <- function(row){
    rep(1, row[3]) %*% t(row[1:2])
  }
  newdat <- do.call("rbind", apply(dati, 1, Repeat_it))
  dati <- as.matrix(newdat)

  #############################
  # Kernel density estimation #
  #############################
  grigliay <- rep(1, quantitest) %*% t(seq(0, 1, length= quantitest))
  grigliax <- t(grigliay)
  # stabilisci i valori di n, alpha e gamma
  n <- dim(dati)[1]   # total number of events

  # KERNEL SMOOTHING
  # input: vector of 2 bandwidths, 2 column matrix of data
  # output: matrix of density estimate vettorizzata

  library("MASS")
  kernel.fun <- function(bandwidths, dati){
    zeta <- kde2d(dati[,1], dati[,2], h= bandwidths, lims = c(0,1,0,1), n=quantitest)$z
    return(c(zeta))
  }

  kernel.spezza.fun <- function(bandwidths, dati){
    ampiezza <- dim(dati)[1] /10
    ris <- 0
    for (blocco in 1:10) {
      datib <- dati[(blocco-1)*ampiezza + 1:ampiezza , ]
      zeta <- kde2d(datib[,1], datib[,2], h= bandwidths, lims = c(0,1,0,1), n=quantitest)$z
      ris <- ris + zeta
    }
    ris <- ris/10
    return(c(ris))
  }

  # TEST STATISTICS
  # input: vector of 2 bandwidths, vector of density estimate
  # output: vector of standardized values

  teststat.fun <- function(bandwidths, fhat) {
    mu0 <- 1
    sigma0 <- sqrt(1/(2*pi*4*bandwidths[1]*bandwidths[2]) - 1)
    teststat <- sqrt(n)*(fhat - mu0)/sigma0
    return(teststat)
  }

  # PETERBARG (it is "almost" it: we do not consider here the area of the set)
  # input: vector of 2 bandwidths, vector of standardized test statistics
  # output: vector of "almost" tail probability

  peterbarg.fun <- function(bandwidths, teststat){
    detC <- (4*bandwidths[1]*bandwidths[2])*(1- 4*bandwidths[1]*bandwidths[2])
    peterbarg <- ((pi*detC)^(-1)) * (teststat^2) * (1-pnorm(teststat))
    return(peterbarg)
  }

  # SUPERSET
  # input: vector of vector of density estimate
  # output: binary vector of superset

  piterpolin <- function(z){(z^2)*(1-pnorm(z))}
  minroot <- optimize (f=piterpolin, interval=c(0, 10), maximum=TRUE)$maximum  # it's the mimimum value for piterbargh

  superset.fun <- function(bandwidths, fhat){
    teststat <- teststat.fun(bandwidths, fhat)
    teststat <- pmax(teststat, minroot)
    sorttest <- sort(teststat, decreasing=TRUE)
    detC <- (4*bandwidths[1]*bandwidths[2])*(1- 4*bandwidths[1]*bandwidths[2])
    piter <- ((pi*detC)^(-1))* (sorttest^2)*(1-pnorm(sorttest))   # piterbarg formula (without area of set)
    piter <- rev(cummin(rev(piter)))      # force pvalues to be monotone increasing
    piter <- piter*(seq(1, 1/length(piter), length=length(piter)))      # include area of the set
    if(max(piter) < alpha){thresh <- min(sorttest - 1)} else{
      thresh <- sorttest[min((1:length(piter))[piter>=alpha])]
    }
    superset <- 1*(teststat <= thresh)
    return(superset)
  }

  onearg.superset.fun <- function(vect){
    superset.fun(vect[1:2], vect[-(1:2)])
  }

  # AUGMENTATION
  # input: binary vector of familywise rejected region, vector of density estimate
  # output: binary vector of augmented region


  augment.fun <- function(fwise, fhat){
    areafwise <- sum(fwise)
    quantinuovi <- floor(areafwise*gammath/(1-gammath)) + length(fwise)*(areafwise<((1-gammath)/gammath))
    scegli <- (1-fwise)*fhat
    thresh <- sort(scegli, decreasing=TRUE)[quantinuovi] + (1+max(fhat))*(areafwise<((1-gammath)/gammath))
    nuovi <- 1*(scegli>=thresh)
    return(nuovi)
  }

  onearg.augment.fun <- function(vect){
    lungh <- length(vect)/2
    augment.fun(vect[1:lungh], vect[lungh+(1:lungh)])
  }

  # SHAVING
  # input: binary vector of rejected region, vector of 2 bandwidths
  # output: binary vector of shaved region

  shave.fun <- function(bands, reject){
    ellypse <- 1*(((grigliax-grigliax[quantitest/2, quantitest/2])/(1*bands[1]))^2 +
                    ((grigliay-grigliay[quantitest/2, quantitest/2])/(1*bands[2]))^2 <= 1)
    ellypse <- rbind(row(ellypse)[ellypse==1], col(ellypse)[ellypse==1])  - quantitest/2
    xmin <- abs(min(ellypse[1,]))
    xmax <- abs(max(ellypse[1,]))
    ymin <- abs(min(ellypse[2,]))
    ymax <- abs(max(ellypse[2,]))
    bigmatrix <- matrix(1, nrow=quantitest + xmin + xmax, ncol=quantitest + ymin + ymax)
    bigmatrix[xmin + (1:quantitest), ymin + (1:quantitest)] <- matrix(reject, nrow=quantitest, ncol=quantitest)

    shift.fun <- function(pair){
      matr <- bigmatrix[xmin + pair[1] + (1:quantitest), ymin + pair[2] + (1:quantitest)]
      return(c(matr))
    }

    shaved <- reject
    for (index in (1:dim(ellypse)[2])){shaved <- shaved*shift.fun(ellypse[,index])}
    #shaved <- apply(ellypse, 2, shift.fun)
    #shaved <- apply(shaved, 1, prod)
    return(shaved)
  }


  # input: a vector where the first 2 entries are the bandwidths and the others are reject
  onearg.shave.fun <- function(vect){
    shave.fun(vect[1:2], vect[-(1:2)])
  }

  oversmooth <-  1.1* sd(c(dati)) /(n^(1/6))
  hvect <- c(seq(1.5/quantitest, oversmooth, length=number_of_bandwidth))
  # bandwidth for x and y
  h1vect <- hvect*sd(dati[,1])/min(sd(dati))
  h2vect <- hvect*sd(dati[,2])/min(sd(dati))
  bandmatr <- rbind(h1vect, h2vect) # 2 * 20 matrix
  fhat <- apply(bandmatr, 2, kernel.spezza.fun, dati=dati) # 65536 * 20 matrix
  bandfhat <- rbind(bandmatr, fhat) # Each column: 2 bandwidth + 256^2 test point

  Ucompl <- 1- apply(bandfhat, 2, onearg.superset.fun)
  Ucomplfhat <- rbind(Ucompl, fhat)
  FDXaugment<- apply(Ucomplfhat, 2, onearg.augment.fun)
  FDXreject <- Ucompl + FDXaugment
  print("2")
  bandUcompl <- rbind(bandmatr, Ucompl)
  Ushaved <- apply(bandUcompl, 2, onearg.shave.fun)

  Ushavedfhat <- rbind(Ushaved, fhat)
  augshaved <- apply(Ushavedfhat, 2, onearg.augment.fun)
  augshaved <- Ushaved + augshaved

  #
  # Drop the 1 pixel situation
  #
  sep_rej <- quantitest * sep_ratio
  inter_sum <- function(vec){
    mat <- matrix(vec, quantitest, quantitest)
    mat[1:floor(sep_rej), 1:floor(sep_rej)] <- 0
    mat[ceiling(sep_rej):quantitest, ceiling(sep_rej):quantitest] <- 0
    return(sum(mat))
  }
  selected_index <- which.max(apply(augshaved, 2, inter_sum))
  finalset <- augshaved[, selected_index]
  toplot <- matrix(finalset, quantitest, quantitest)
  dat <- reshape2::melt(toplot)
  p1 <- ggplot(dat, aes(Var1, -Var2)) +
    geom_tile(aes(fill = value), show.legend = F) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_classic() +
    theme(aspect.ratio = 1,
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,0,0,0), "lines"))+
    geom_hline(yintercept = -sep_rej) + geom_vline(xintercept = sep_rej)
  png("RR.png")
  print(p1)
  dev.off()

  # extract kde matrix with the best selected bandwidth for breakpoint selection
  kde_mat = matrix(fhat[,selected_index], quantitest, quantitest)
  kde_long = reshape2::melt(kde_mat)

  # write out rejection regions and kde smoothed matrix
  if (inter_sum(finalset) >= 3){
    write.table(toplot, file = "RR.txt", row.names = F, col.names = F)
    write.table(kde_mat, file = "KDE.txt", row.names = F, col.names = F)
  }
  save(list=ls(all.names = TRUE), file="all_objects.Rdata")
}

