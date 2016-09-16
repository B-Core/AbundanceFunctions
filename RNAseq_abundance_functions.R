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

regressMatrix = function(normmat, expt.design, lm_expression, response_var="y", contr_list=NULL, source_me=c("/home/users/burchard/R-utils/lm.pval.R")){
  # function to run linear regression on multiresponse data matrix with given model and optional contrasts
  # Arguments
  # source_me: character vector of R function file names to source
  # normmat: a matrix of normalized (sample-bias-reduced) abundance data
  #  rownames are unique identifiers for features measured per row
  #  colnames are unique identifiers for samples
  # expt.design: a named list of character vectors, each vector having one element per sample in raw.mat column order, each string representing a design factor OR covariate in the conduct of the experiment
  #   examples: Treatment, Gender, Date, Contaminant_marker_gene_abundance, Sample_purity
  # lm_expression: a string comprising an expression of the form "y  ~ model"
  # response_var gives the string to be used as the name of the response matrix
  #   by default, y is used to represent the multiresponse data matrix
  #   the predictors in the model MUST be named elements of expt.design
  #   example: "y ~ Treatment + Gender"
  # contr_list: a named list, with one element per expt.design element name used in lm_expression, whose contents are either a valid contrasts matrix or a list with elements "baseline" and "contr.FUN"
  #   if contrasts is NULL, default coding will be used
  #      R's default: dummy or treatment contrast, comparing each level to base
  #        by default, alphanumerically first factor level is base
  #        regression coefficients are adjusted group mean ratios to base
  #      useful alternative: summed contrast, comparing levels to grand mean
  #        by default, alphanumerically last factor level is omitted
  #        regression coefficients are adjusted group mean ratios to grand mean
  #   example matrix: dummy or treatment contrast with non-default baseline
  #                   this matrix causes Tx abundances to be compared to Sham
  #                     14d TxA  14d TxB
  #      14d Sham          0        0
  #      14d TxA           1        0
  #      14d TxB           0        1
  #   example list: this list will cause creation of the matrix above
  #   list(baseline="Sham", contr.FUN="contr.treatment")

  # imports
  require(qvalue)
  source(source_me) # very trusting!

  # constants
  contr_list_names = c("baseline","contr.FUN"); fundx=2; basedx=1

  # calculate contrasts matrices if contrasts input is provided
  if( !is.null(contr_list) ){
    contr_lsmat = vector(mode='list',length=length(contr_list))
    for( fac in names(contr_list) ){
      # test input assumptions
      if( !any(grepl(fac, names(expt.design))) ){
        stop(paste("Factor",fac,"in contrasts list is not in exptl design list"))
      }
      if( !typeof(expt.design[[fac]]) %in% c("character","double","integer") | 
        !is.null(attr(expt.design[[fac]],"dim")) ){
        stop(paste("Factor",fac,"is not a string or numeric vector"))
      }
      # use tested input
      if(is.list(contr_list[[fac]]) ){
        if(sum(names(contr_list[[fac]]) %in% contr_list_names)<length(contr_list_names) ){
          stop("Missing",paste(contr_list_names,collapse=" or "),"for",fac)
        }
        # create matrix for list input
        myfun = contr_list[[fac]][[contr_list_names[fundx]]]
        mybase = contr_list[[fac]][[contr_list_names[basedx]]]
        if(myfun=="contr.treatment"){
          contr_lsmat[[fac]] = contr.treatment( levels(as.factor(expt.design[[fac]])), base = which(levels(as.factor(expt.design[[fac]])) == mybase) )
        } else if(myfun=="contr.sum"){
          fac_v = c( setdiff(levels(as.factor(expt.design[[fac]])), mybase), mybase) # put baseline last for contr.sum
          contr_lsmat[[fac]] = contr.sum(fac_v)
        } else {
          stop(paste("Contrast function",myfun,"is not yet supported"))
        }
        # recover matrix input
      } else if(is.matrix(contr_list[[fac]]) ){
        contr_lsmat[fac] = contr_list[fac]
      } else if(class(contr_list[[fac]])=="AsIs" & length(dim(contr_list[[fac]]))==2 ){ # matrix wrapped in I()
        class(contr_lsmat[[fac]]) = "matrix"
        contr_lsmat[fac] = contr_list[fac]
      } else {
        stop(paste("Contrasts specification for",fac,",is not a list or matrix"))
      }
    }
  }
  
  # set up regression formula with environment
  # test input
  if( typeof(lm_expression) != "character" | !is.null(attr(lm_expression,"dim")) ){
    stop("Regression formula was not supplied as a string") 
  }
  # check supplied factors for presence in formula
  lm_fac = regmatches(lm_expression,gregexpr('(\\w+)',lm_expression))[[1]]
  lm_fac = setdiff(lm_fac, response_var)
  n = length(lm_fac) - sum(lm_fac %in% names(expt.design)) 
  if( n != 0 ){
    stop(n,"regression factors not found in exptl design list")
  }
  # pull out factors supplied in exptl design for use in lm, & assign contrasts
  lm_list = lapply(expt.design[names(expt.design) %in% lm_fac],as.factor)
  # add contrasts if contrasts input was given
  if( !is.null(contr_list) ){
    for(fac in names(lm_list) ){
      if( any(grepl( fac, names(contr_lsmat) )) ){
        contrasts(lm_list[[fac]]) = contr_lsmat[[fac]]
      }
    }
  }

  # set up environment for regression, with desired factors present
  # use parent of globalenv as parent of regression env to avoid confounding
  lm_env = list2env( lm_list, parent=parent.env(globalenv()) )
  assign(x=response_var, value=t(normmat), envir=lm_env)
  lm_formula = as.formula( lm_expression, env=lm_env)

  # run regression 
  obj_lm = lm( lm_formula)
  b_mat = t(obj_lm$coefficients)
  p_mat = t(lm.pval(obj_lm)$pval); class(p_mat) = "matrix"
  # replace NA p-values (bad for qvalue), at risk of spike at 1 (bad for qvalue)
  p_mat[is.na(p_mat)] = 1

  # convert p-values to q-values for evaluation of factor impact on abundances
  q_list = NULL
  for( fac in colnames(p_mat) ){
    myflag = TRUE
    obj_qvalue = tryCatch( qvalue(p_mat[,fac]), 
          error=function(e){myflag=FALSE; cat("No qvalues for",fac,"\n")} )
    if( myflag) { q_list[[fac]] = I(obj_qvalue) }
  }

  # return values
  return( list( b_mat=b_mat, p_mat=p_mat, q_list=q_list ) )

}

