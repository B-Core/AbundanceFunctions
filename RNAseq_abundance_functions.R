# draft generalized abundance functions

read_STAR = function( useme.cols, label.from.colname, annCol.label, annCol.names, annCol.normBy, annCol.lmBy, readir = '.', readpattern = '.', outypes = list(gene.counts='Gene.out.tab',SJ.counts='SJ.out.tab'), outcol = list(gene.counts = c("Gene","Reads","FwdReads","RevReads"),SJ.counts = c("Chr","IntronStart","IntronStop","Strand","IntronMotif","Annotated","UniqueReads","MultiReads","MaxOverhang")), outkey = list(gene.counts = 1, SJ.counts = 1:3), outsum = list(gene.counts = c("N_ambiguous","N_multimapping","N_noFeature","N_unmapped"),SJ.counts = c()), stranded.col = list(gene.counts = c(1,3,4),SJ.counts = c(1:3,7)), unstranded.col = list(gene.counts = c(1,2),SJ.counts = c(1:3,7)), stranded=F, filesep="/", write2file=TRUE){
  # read in STAR output from processing of RNAseq data
  # format 08/26/2016: one directory per FASTQ file, named after FASTQ file
  # contains gene count and splice-junction count files
  # arguments:
  # useme.cols: '(C|D)' -- custom regex for processed data to retain and normalize
  # label.from.colname: '^.*AB_([0-9]+_[^_]+).*$' -- regex to extract sample labels from FASTQ file names
  # annCol.label: c('^.*AB_[0-9]+_[0-9]+-([^_]+).*$') -- regex(es) to extract exptl design factors from FASTQ filenames
  # annCol.names: c('condition') -- vector of titles for exptl design factors
  # annCol.normBy: for later, exptl design factor indicating batches to normalize separately, or NULL if no batch normalization is intended
  # annCol.lmBy: for later, formula for regression testing of normalized data
  # readdir: directory in which to find STAR output directories
  # readpattern: pattern by which to recognize the STAR output directories
  # !!the following 4 arguments MUST have the _same lenghts & element names_!!
  # outypes: suffixes of STAR output file types to crunch
  # outcol: columns of each STAR output format
  # outkey: number(s) of columns to use as keys
  # outsum: rownames of summary stats saved by STAR as data rows
  # stranded.col: col to save for stranded library prep
  # currently saving both Fwd and Rev, but Fwd may suffice -- test this
  # unstranded.col: col to save for unstranded library prep
  # stranded: logical flag indicating whether library prep was strand-specific
  # filesep: delimiter for pathnames

  # imports
  require(data.table)
  
  # collate raw gene counts
  samps = dir(path=readir, pattern=readpattern)
  # trimmed sample labels
  samp.labels = gsub(label.from.colname,'\\1', samps)

  # intermediate structure for results
  tab = vector(mode="list",length=length(samps))
  for(i in 1:length(samps)){
    for(ctype in names(outcol) ){
      mysuffix = outypes[[ctype]]
      mycol = outcol[[ctype]]
      mykey = outkey[[ctype]]
      myfile = paste(samps[i],dir(samps[i],mysuffix),sep=filesep)
      tab[[i]][[ctype]] = fread(myfile, verbose=F)
      # type-specific processing
      names(tab[[i]][[ctype]])[mykey] = mycol[mykey]
      names(tab[[i]][[ctype]])[setdiff(1:length(mycol),mykey)] = paste(mycol[setdiff(1:length(mycol),mykey)],samps[i],sep='.')
      setkeyv(tab[[i]][[ctype]],mycol[mykey])
    }
  }
  # combined data.tables
  comb.dat = vector(mode="list",length=length(outypes))
  names(comb.dat) = names(outypes)
  for(ctype in names(outypes) ){
    for(i in 1:length(samps)){
      # set up column indices as needed
      if(stranded){ mycol = names(tab[[i]][[ctype]])[ stranded.col[[ctype]] ]
      } else { mycol = names(tab[[i]][[ctype]])[ unstranded.col[[ctype]] ] }
      if( length(comb.dat[[ctype]])>0 ){
        comb.dat[[ctype]] = merge(comb.dat[[ctype]], tab[[i]][[ctype]][, mget(mycol)],all=T )
      } else {
        comb.dat[[ctype]] = tab[[i]][[ctype]][,mget(mycol)] }
    }
  } 
  
  # data matrices
  LoM.raw = vector( mode='list',length=length(names(outypes)) )
  names(LoM.raw) = names(outypes)
  for(ctype in names(outypes) ){
    if(stranded){ mycol = outcol[[ctype]][ setdiff(stranded.col[[ctype]],outkey[[ctype]]) ]
    } else { mycol = outcol[[ctype]][ setdiff(unstranded.col[[ctype]],outkey[[ctype]]) ] }
    # prepare row labels
    vec = apply(MARGIN=1, X=comb.dat[[ctype]][,mget(outcol[[ctype]][outkey[[ctype]] ])], FUN=paste0, collapse='.'); vec = gsub('\\s+','',vec)
    # split count types into separate matrices
    LoM.raw[[ctype]] = vector( mode='list',length=length(mycol) )
    names(LoM.raw[[ctype]]) = mycol
    for(colx in mycol){
      LoM.raw[[ctype]][[colx]] = as.matrix(comb.dat[[ctype]][,mget(names(comb.dat[[ctype]])[ grepl(useme.cols,names(comb.dat[[ctype]])) & grepl(paste0('^',colx),names(comb.dat[[ctype]])) ])])
      rownames(LoM.raw[[ctype]][[colx]]) = vec
      colnames(LoM.raw[[ctype]][[colx]]) = samp.labels
      #colnames(LoM.raw[[ctype]]) = paste(mycol,rep(samp.labels, each=round(ncol(LoM.raw[[ctype]])/length(samp.labels)) ), sep='.')
      # remove any STAR summary stats rows
      if( length(outsum[[ctype]]) > 0 ){
        rowmask = !rownames(LoM.raw[[ctype]][[colx]]) %chin% outsum[[ctype]]
        LoM.raw[[ctype]][[colx]] = LoM.raw[[ctype]][[colx]][rowmask,]
      }
    }
  }

  # annotations for columns, as specified in arguments
  annCol = NULL
  for(i in 1:length(annCol.label) ){annCol = c(annCol,list(gsub(annCol.label[i],'\\1', samps )) )}
  names(annCol) = annCol.names
  attr(annCol,'normBy') = annCol.normBy
  attr(annCol,'lmBy') = annCol.lmBy

  # write to file
  if(write2file){
    for(ctype in names(LoM.raw) ){
      my.dt = data.frame(LoM.raw[[ctype]],keep.rownames=TRUE)
      write.csv(my.dt,file=paste('raw',ctype,"csv",sep='.'),quote=F)
    }
    cat(sapply(annCol, toString), file='expt.design.txt', sep="\n")
  }

  # return list with data matrix and exptl design
  return(list(LoM.raw = LoM.raw, expt.design = annCol))
}


norm_matrix = function(tag, raw.mat, expt.design, normvec=c("loess","qspln","quant"), normFUN=c("normalize.loess","normalize.qspline","preprocessCore:::normalize.quantiles"), normarg=list(loess=list(family=c("symmetric","gaussian")),qspln=list(samples=c( max(round(nrow(raw.mat)/1000), 100),12*nrow(raw.mat)^(-.7)),na.rm=TRUE),quant=list(copy=c(TRUE,FALSE))), normPkg=c("affy","affy","affy"), depth.est = list(upper.quartile=0.75, max=1), bkgdFUN=function(x,probs=.75){quantile(x,probs=probs,na.rm=TRUE)/10^(max(c(trunc(log10(quantile(x,probs=probs,na.rm=TRUE)))-1,1)) )} ){
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


  # not yet implemented:
  # if indicated, normalize by selected annCol


  # normalizations used for microarray RNA expression data
  # add background to raw data for logging and low-end noise amelioration
  LoM.norm = vector(mode='list',length=length(normvec) )
  names(LoM.norm) = normvec
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


