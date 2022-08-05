
# library (aricode)
# https://cran.r-project.org/web/packages/aricode/aricode.pdf
#rm(list = ls())

library(SpatialPack)
library(magick)
library(aricode)
# library for purity
library(funtimes)

library(flexclust)
library(scales)
library(cluster)
# Library for ploting similarity matrix
library(klic)
library(plotrix)
library(rapportools)
library(survival)
library(SNFtool)

# ****************************** IMG core functions begin *********************************
library(geometry)



#' bins_fd
#' @description returns the number of bins according to the Freedman-Diaconis rule
#' @param vec numeric vector
#' @return bin width
bins_fd <- function(vec) {
  #diff(range(vec)) / (2 * IQR(vec) / length(vec)^(1/3))
  fd<- (2 * IQR(vec)) / (length(vec)^(1/3))
  fd
}

bins_scotts <- function(vec) {
  bins <- nclass.scott(vec)
  bin_width <- (max(vec) - min(vec))/bins
  bin_width
  
}

bins_sturges <- function(vec) {
  bins <- nclass.Sturges(vec)
  bin_width <- (max(vec) - min(vec))/bins
  bin_width
}

bins_rice <- function(vec) {
  bins <- 2 * (length(vec)^(1/3))
  bin_width <- (max(vec) - min(vec))/bins
  bin_width
}

get_bin_width <- function(vec, bin_rule){
  bin_width <- 0
  if(bin_rule == "fd"){
    #print("using fd bin width estimator")
    bin_width <- bins_fd(vec)
  }else if(bin_rule == "scotts"){
    #print("using scotts bin width estimator")
    bin_width <- bins_scotts(vec)
  }else if (bin_rule == "sturges"){
    #print("using sturges bin width estimator")
    bin_width <- bins_sturges(vec)
  }else if (bin_rule == "rice"){
    #print("using rice bin width estimator")
    bin_width <- bins_rice(vec)
  }else {
    print("Error: no bin rule is defined!! using fd as default")
    bin_width <- bins_fd(vec)
  }
  
  bin_width
}

# *** IM RBF Kernel Starts Here... ****
IM_RBF <- function(X, Y, sgma=0.5, bin_rule){
  
  wts <- c(rep(1, dim(X)[2]))
  mx <- c(rep(1, dim(X)[2]))
  
  for (k in 1:dim(X)[2]) {
    k_width <-  get_bin_width(X[,k], bin_rule)
    #print(xmat[k,])
    if(is.empty(k_width)){
      k_width <- 1
    }
    #print(paste("bin_width = ", k_width, sep = ""))
    wts[k] <- k_width
    # get max bucket number or the number of buckets
    mx[k] <- (max(X[,k]) - min(X[,k]))/k_width
    if(mx[k] == 0){
      print("bin_width is zero...")
      mx[k] = 1
    }
    
    #print(paste("mx[k] = ", mx[k], sep = ""))
  }

  
  xwts = wts * mx
  X = X/xwts
  Y = Y/xwts
  
  xtrns1 <- matrix(1, nrow = dim(X)[1], ncol = 1)
  xtrns2 <- matrix(1, nrow = dim(Y)[1], ncol = 1)
  for( i in 1:dim(X)[1]){
    xtrns1[i,1] = X[i,] %*% X[i,]
  }
  for( i in 1:dim(Y)[1]){
    xtrns2[i,1] = Y[i,] %*% Y[i,]
  }
  
  k1 =  xtrns1 %*% t(matrix(1, nrow = dim(X)[1], ncol = 1))
  k2 =  matrix(1, nrow = dim(Y)[1], ncol = 1) %*% t(xtrns2)
  k = k1 + k2
  
  k = k - (2 * (X %*% t(Y)))

  k = k * (-1. / (2 * (sgma**2)))
  
  exp(k)
}
# *** IM Kernel ends *****



########### ROMDEX FUNCTION #################################
RomdexML <- function(xdf, sigval = 0.5, bin_rule = "fd"){
  
  xdfnorm = standardNormalization(xdf)
  
  t1=Sys.time()
  
  W1 = IM_RBF(xdfnorm, xdfnorm, sigval, bin_rule)
  #W1 = IM_RBF(xdf, xdf, sigval, bin_rule)
  
  
  an.error.occured <- FALSE
  estClus <- NULL
  estGroups <- NULL
  tryCatch( {
    
    #Groups with just geneExpression
    #estClus = estimateNumberOfClustersGivenGraph(W1, 2:10)[[1]]
    estGroups = spectralClustering(W1,7);
    
  } , error = function(e) {
    an.error.occured <<- TRUE
    estGroups = c(1:1474)
  })
  if(an.error.occured){
    print('error occured.')
    estGroups = c(1:1474)
  }
  
  t3 = Sys.time() - t1
  #print(paste('time taken: ', t3))
  #Cox log rank test GE: Gene Expression
  factor(estGroups);
}


########### Main #####################

#xdf <- xdf[, which(colMeans(!is.empty(xdf)) > 0.65)]

dataDir = "./caltech101-7/"
clabels = c("Faces", "Motorbikes", "dollar_bill", "garfield", "snoopy", "stop_sign", "windsor_chair")
#clabels = c("LFaces", "LMotorbikes")

obsr = 0
for (cats in clabels) {
  obsr = obsr + length(list.files(paste(dataDir, cats, sep = "")))
}
print(paste("observations: ", obsr))

true_labels = c(1:obsr)
lbl = 0

dtMat = matrix(0, obsr, 26400)
rowIdx = 0
for (lb in clabels) {
  cDir = paste(dataDir, lb, sep = "")
  fildir = list.files(cDir)
  
  lbl = lbl+1
  for (f in fildir){
    rowIdx = rowIdx+1;
    
    imagePath = paste(cDir, "/", f, sep ="")
    #print(imagePath)
    img_small_gray = image_convert(image_scale(image_read(imagePath), "220x120!"), type = 'grayscale')
    x <- as.array(image_data(img_small_gray)[1,,])
    #dtMat[rowIdx,] = as.double(as.vector(x))/255
    dtMat[rowIdx,] = as.integer(as.vector(x))
    true_labels[rowIdx] = lbl
  }
}

print(paste("labels: ", length(true_labels), ", instances: ", rowIdx))


#ari = 0

dims = 1
vari = numeric(dims)
vnmi = numeric(dims)
vpurity = numeric(dims)
vsigma = numeric(dims)

strt = 0
for (i in 1:dims){
  
  vsigma[i] = strt + 15.75
  
  groups = RomdexML(dtMat, sigval = vsigma[i], bin_rule = "fd")
  
  vari[i] = ARI(groups, true_labels)
  vnmi[i] = NMI(groups, true_labels)
  vpurity[i] = purity(groups, true_labels)[1]$pur
  
  print(paste('sigma: ', vsigma[i], 'ARI: ', vari[i], 'NMI: ', vnmi[i], 'Purity: ', vpurity[i], '\n'))
  
  
}

ml_results <- data.frame('sigma' = vsigma, 'ARI' = vari, 'NMI' = vnmi, 'PURITY' = vpurity)
write.csv(ml_results, "./caltech-7-norm-2-fd-full.csv", row.names = FALSE)

