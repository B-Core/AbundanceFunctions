################################################################################
# Functions involved in reducing technical bias
#
# normMatrix       --- apply bias reduction methods to data matrices
# normalizeLoess   --- fork of normalize.loess with NA tolerance
#
# Authors: Theresa Lusardi, Julja Burchard, and Mark Fisher
# Started - July 2016
################################################################################

normMatrix <-
function(tag, raw.mat, expt.design, 
                normvec=c("loess","qspln","quant"), 
                normFUN=c("normalizeLoess","normalize.qspline","preprocessCore:::normalize.quantiles"), 
                normarg=list( loess=list(method="loess"),
                             qspln=list(samples=c( max(round(nrow(raw.mat)/1000), 100),12*nrow(raw.mat)^(-.7)),na.rm=TRUE),
                             quant=list(copy=c(TRUE,FALSE))), 
                normPkg=c("affy","affy","affy"), 
                depth.est = list(upper.quartile=0.75, max=1), 
                bkgdFUN=function(x,probs=.75){quantile(x,probs=probs,na.rm=TRUE)/10^(max(c(trunc(log10(quantile(x,probs=probs,na.rm=TRUE)))-1,1)) )} 
                ){
  # Add lograw to normvec ##Feedback ##TS ##MF
  # apply bias reduction functions to abundance data
  # tag: string identifier for output files. suggest including data type.
  # raw.mat: a matrix of raw data, rownames are unique feature identifiers, column names are unique sample identifiers
  # expt.design: a named list of character vectors, each vector having one element per sample in raw.mat column order, each string representing a design factor
  # it has optional attributes: 
  #   normBy: element name(s) by which to split data before norm
  #   lmBy: for later, formula for regression using expt.design elements
  # the following 4 vectors _MUST have the same # elem in same order_
  # normvec: vector of tags for bias reduction methods
  # normFUN: vector of string names for these bias reduction functions 
  # normarg: list of arguments to pass to each bias reduction function
  #          if 0 elements, no extra args are passed in do.call
  #          if 1 value for an argument, will be passed in do.call
  #          if 2 values, 2nd will be used as alt in tryCatch
  #          >2 values per argument are currently not supported
  #  notes: R lazy eval allows function args to ref other args not yet assigned
  #         qspline args are from function default and help page, respectively
  # normPkg: vector of string names of packages containing these functions
  # depth.est: named list of quantiles per sample to report to console
  # bkgdFUN: string specifying function to use to calculate bkgd to add

  # require package(s) with norm functions in them
  uPkg = unique(normPkg)
  for(i in 1:length(uPkg) ){
    eval(parse(text=paste('require(',uPkg[i],')')))
  }
  
  # depth estimate(s)
  for(mystat in names(depth.est)){
    cat(mystat,'of total units per sample\n')
    cat( sprintf(signif(sapply(1:ncol(raw.mat),function(x){quantile(raw.mat[,x],probs=depth.est[[mystat]],na.rm=TRUE)}),2),fmt='%1.1e'), "\n")
  }
  # set background adjustment
  bkgd = bkgdFUN(raw.mat)
  message(sprintf("Set values: min = %1.2f, background = %1.2f", min(raw.mat), bkgd))

  # normalizations used for microarray RNA expression data
  # add background to raw data for logging and low-end noise amelioration
  LoM.norm = vector(mode='list',length=length(normvec)+1 )
  names(LoM.norm) = c("alograw", normvec)
  # Create a log2 transform of the raw matrix. Do not add background!
  mynorm = "alograw"
  message(sprintf("mynorm = %s", mynorm))
  LoM.norm[[mynorm]] = log2(raw.mat+1)

  for( i in 1:length(normvec) ){
    mynorm = normvec[i]
    myargs = normarg[[i]]
    if(grepl(':::',normFUN[i]) ){
      ans = strsplit(normFUN[i],':::')
      myf = getFromNamespace(x=ans[[1]][2],ns=ans[[1]][1])
    } else {
      myf = get(normFUN[i])
    }
    if( length(myargs)==0 ){ 
      LoM.norm[[mynorm]] = do.call(myf,list(raw.mat))
      cat(mynorm,": no args given\n")
    } else if( any(lengths(myargs)>1) ){
      myflag = FALSE
      LoM.norm[[mynorm]] = tryCatch(
      log2(do.call( myf,c(list(raw.mat+bkgd),lapply(myargs,FUN=function(x){x[1]})) )),
      error=function(e){myflag=TRUE; log2(do.call( myf,c(list(raw.mat+bkgd),lapply(myargs,FUN=function(x){x[min(c(2,length(x)))]} )) )) } 
      )
      cat(mynorm,": ran with 2 args :",if(myflag){unlist(lapply(myargs,FUN=function(x){x[min(c(2,length(x)))]} ))
      }else{ unlist(lapply(myargs,FUN=function(x){x[1]})) },"\n")
    } else {
      LoM.norm[[mynorm]] = log2(do.call( myf,c(list(raw.mat+bkgd),myargs) ))
      cat(mynorm,": ran with 1 arg :",unlist(myargs),"\n")
    }
    # some norms remove string annotation!! 
    # verify that matrix is not reordered before using new norms
    colnames(LoM.norm[[mynorm]]) = colnames(raw.mat)
    rownames(LoM.norm[[mynorm]]) = rownames(raw.mat)
  }

  # write normalized data to files
  for( mynorm in normvec ){
    my.dt = data.table(LoM.norm[[mynorm]],keep.rownames=TRUE) ##feedback note this is not a data table..
    # data.table names the rowname column "rn"
    names(my.dt)[names(my.dt)=="rn"] = "Identifier"
    write.csv(my.dt,file=paste(mynorm,tag,"csv",sep='.'),quote=F) 
  }

  # return list of normalized matrices
  return(LoM.norm)
}
normalizeLoess <-
function (rawmat, rowlim = 
              if(method=="loess"){ min(20000, nrow(rawmat), na.rm=T )
                }else{ min(2000, nrow(rawmat), na.rm=T) },
              subset = order(rowMeans(log2(rawmat),na.rm=T))
                [round(seq(from=1,to=nrow(rawmat),by=nrow(rawmat)/rowlim))],
              method = "loess", epsilon = 10^-2, maxit = 1, log.it = TRUE, 
              span = 2/3, family.loess = c("symmetric","gaussian"),
              verbose = TRUE
                            ) {
  # rawmat:  a matrix with one column per abundance profile to normalize
  # rowlim:  loess: maximum number of rows of data to include in fit
  #          lowess: denominator, fraction of x range by which to space fit pts
  # subset:  a subset of the data to which to fit a loess 
  # epsilon: a tolerance value (supposed to be a small value 
  #          - used as a stopping criterion).
  # maxit:   maximum number of iterations
  # log.it:  logical. If TRUE, the log2 of rawmat is fitted
  # verbose: logical. If TRUE, reports current pair of chips being fitted
  # method:  loess (standard formulaic) or lowess (fast interpolated) regression
  # span:    loess/lowess parameter alpha which controls the degree of smoothing
  # family.loess: if "gaussian" loess fitting is by least-squares
  #               if "symmetric" loess uses a re-descending M estimator with
  #               Tukey's biweight function

  cat(method,"regression with",rowlim,"points","\n")
  # arguments
  # rawmat dimensions
  J = dim(rawmat)[2]
  II = dim(rawmat)[1]
  # fit loess to log MA plot? (default)
  if (log.it) {
    rawmat = log2(rawmat)
  }
  # loess controls
  family.loess = family.loess[1] # use first option if several are given
  change = epsilon + 1
  iter = 0
  w = c(0, rep(1, length(subset)), 0)

  # iteratively fit loess/lowess to virtual pairwise MA plots & save adjustments
  #  based on assumption that mean ratio is zero at all abundances
  while (iter < maxit) {
    # loop parameter
    iter = iter + 1

    # initialize matrix of adjusments for this iteration
    means = matrix(0, II, J)
    # loop over all unique pairs of columns
    for (j in 1:(J - 1)) {
      for (k in (j + 1):J) {

        # calculate ratio (y) and mean abundance (x)
        #  use non-NA value as x if one of two columns has missing value
        y = rawmat[, j] - rawmat[, k]
        x = rowSums( cbind(rawmat[, j], rawmat[, k]), na.rm=T )/2

        # select subset to fit
        #   default subset is selected for range and reproducibility
        # original pads subset with extrema with zero weights
        index = c(order(x)[1], subset, order(-x)[1])

        # extract subset to fit
        xx = x[index]
        yy = y[index]

        # fit loess/lowess to subset
        if( method=="loess") {
          aux = loess(yy ~ xx, span = span, degree = 1, 
            weights = w, family = family.loess)
          # predict adjustments to abundances. note that predict uses na.pass
          aux = predict(aux, data.frame(xx = x))/J
        } else if( method=="lowess") {
          # no prediction -- fit all points
          aux = lowess(x=x, y=y, f=span, 
                       delta=diff(range(x,na.rm=T))/rowlim )$y/J
        } else { stop( paste("No method", method) ) }
        # add adjustments to adjustment matrix
        means[, j] = rowSums( cbind(means[, j], aux), na.rm=T )
        means[, k] = rowSums( cbind(means[, k], -aux), na.rm=T )
        if (verbose) {
          cat("Done with", j, "vs", k, "in iteration", iter, "\n")
        }
      }
    }
    rawmat = rawmat - means
    change = max(colMeans((means[subset, ])^2))
    if (verbose) { cat(iter, change, "\n") }
  }
  if ((change > epsilon) & (maxit > 1)) { 
    warning(paste("No convergence after", maxit, "iterations.\n"))
  }
  if (log.it) {
    return(2^rawmat)
  } else { return(rawmat) }
}
