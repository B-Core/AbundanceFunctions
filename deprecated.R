###############################################################################
# Parking area for AbundanceFunctions not currently in use.
#
# norm_matrix  --- Apply bias reduction functions to abundance data
#                  <== Removed from RNAseq_abundance_functions.R
#                  ==> Replaced by normMatrix in processData.R
#                      15-Oct-2016
#
###############################################################################


norm_matrix = function(tag, raw.mat, expt.design, 
            normvec=c("loess","qspln","quant"), 
            normFUN=c("normalize.loess","normalize.qspline","preprocessCore:::normalize.quantiles"), 
            normarg=list(
                  loess=list(family=c("symmetric","gaussian")),
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
    cat(mystat,'of total read counts per sample\n')
    cat( sprintf(signif(sapply(1:ncol(raw.mat),function(x){quantile(raw.mat[,x],probs=depth.est[[mystat]],na.rm=TRUE)}),2),fmt='%1.1e'), "\n")
  }
  # set background adjustment
  bkgd = bkgdFUN(raw.mat)
  message(sprintf("Set values: min = %1.2f, background = %1.2f", min(raw.mat), bkgd))

  # not yet implemented:
  # if indicated, normalize by selected annCol

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
    my.dt = data.frame(LoM.norm[[mynorm]],keep.rownames=TRUE) ##feedback note this is not a data table..
    write.csv(my.dt,file=paste(mynorm,tag,"csv",sep='.'),quote=F) 
  }

  # return list of normalized matrices
  return(LoM.norm)
}
