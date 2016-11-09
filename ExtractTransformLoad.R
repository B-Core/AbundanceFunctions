################################################################################
# Functions specific to extraction, transformation, loading
#
# read_STAR       --- Read in STAR gene and splice junction count data
# processCELfiles --- Read in Affymetrix expression microarray data
#
################################################################################


read_STAR <-
function( useme.cols, label.from.colname, annCol.label, 
          annCol.names, annCol.normBy=NULL, annCol.lmBy=NULL, 
          readir = '.', readpattern = '.', 
          outypes = list(gene.counts='Gene.out.tab',SJ.counts='SJ.out.tab'), 
          outcol = list(gene.counts = c("Gene","Reads","FwdReads","RevReads"),
          SJ.counts = c("Chr","IntronStart","IntronStop","Strand","IntronMotif","Annotated","UniqueReads","MultiReads","MaxOverhang")), 
          outkey = list(gene.counts = 1, SJ.counts = 1:3), 
          outsum = list(gene.counts = c("N_ambiguous","N_multimapping","N_noFeature","N_unmapped"), SJ.counts = c()), 
          stranded.col = list(gene.counts = c(1:4),SJ.counts = c(1:3,7:8)), 
          unstranded.col = list(gene.counts = c(1:4),SJ.counts = c(1:3,7)), 
          stranded=F, filesep="/", write2file=TRUE
          ){
  # read in STAR output from processing of RNAseq data
  # return list : LoM.raw contains gene and SJ count matrics
  #               expt.design is a list with expt info parsed from FASTQ names
  #               myreads indicates the properly stranded version of gene counts
  #       (list(LoM.raw = LoM.raw, expt.design = annCol, myreads=myreads))
  # input format 08/26/16: one directory per FASTQ file, named after FASTQ file
  #                        contains gene count and splice-junction count files
  # arguments:
  # useme.cols: '(C|D)' -- custom regex for processed data to retain and normalize
  # label.from.colname: '^.*AB_([0-9]+_[^_]+).*$' -- regex to extract sample labels from FASTQ file names
  # annCol.label: c('^.*AB_[0-9]+_[0-9]+-([^_]+).*$') -- regex(es) to extract exptl design factors from FASTQ filenames
  # annCol.names: c('condition') -- vector of titles for exptl design factors
  # annCol.normBy: for later, exptl design factor indicating batches to normalize separately, or NULL if no batch normalization is intended
  # annCol.lmBy: for later, formula for regression testing of normalized data
  # readdir: directory in which to find STAR output directories
  # readpattern: pattern by which to recognize the STAR output directories
  # !!the following 4 arguments MUST have the _same lengths & element names_!!
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
      myfile = file.path(readir,samps[i],dir(file.path(readir,samps[i]),mysuffix))
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

  # test to confirm strandedness of gene-based data type
  i = grep('gene', names(LoM.raw),ignore.case=T)[1]
  fwdmk = grepl('fwd|forward|sense', names(LoM.raw[[i]]), ignore.case=T)
  revmk = grepl('rev|anti', names(LoM.raw[[i]]), ignore.case=T)
  f = which(fwdmk)[1]; r = which(revmk)[1]; u = which(!fwdmk & !revmk)[1]
  if( any(is.na( c(f,r,u) )) ){
    message("Forward: ",f,"; Reverse: ",r,"; Unstranded: ",u)
    stop('Gene-based data must provide forward, reverse and unstranded types')
  }
  testv = log2(colSums(LoM.raw[[i]][[u]],na.rm=T)+1) - log2(colSums(LoM.raw[[i]][[f]],na.rm=T)+colSums(LoM.raw[[i]][[r]],na.rm=T)+1)
  if( sum(round(testv,0))>1 ){
      stop("Forward and reverse reads don't add up to total reads!")
  } else {
    testv =  log2(colSums(LoM.raw[[i]][[f]],na.rm=T)) - log2(colSums(LoM.raw[[i]][[r]],na.rm=T)+1)
    if( sum(round(testv,0))>1){            # fwd strand excess over rev
      myreads = names(LoM.raw[[i]])[f]
      cat(" Using", myreads, "gene counts","\n")
    } else if( sum(round(testv,0))<=-1){   # rev strand excess over fwd
      myreads = names(LoM.raw[[i]])[r]
      cat(" Using", myreads, "gene counts","\n")
    } else {
      myreads = names(LoM.raw[[i]])[u] # no strand inequality
      cat(" Using", myreads, "gene counts","\n")
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
      for(colx in names(LoM.raw[[ctype]]) ){
        my.dt = data.table(LoM.raw[[ctype]][[colx]],keep.rownames=TRUE)
        # data.table names the rowname column "rn"
        names(my.dt)[names(my.dt)=="rn"] = "Identifier"
        write.csv(my.dt,file=paste('raw',ctype,colx,"csv",sep='.'),
                  quote=F, row.names=F)
      }
    }
    cat(sapply(annCol, toString), file='expt.design.txt', sep="\n")
  }

  # return list with data matrix and exptl design
  return(list(LoM.raw = LoM.raw, expt.design = annCol, myreads=myreads))
}
Rm_quantile_of_very_different_intensities <-
function(lin_scale_abundance_mat, quant_val){
  #Removes the quant_val...th quantile of the most differing abundance metric (e.g., probe intensity) across an expression matrix; returns the remainder of the matrix
  #lin_scale_abundance_mat is a matrix of linear scale abundance data
  #quant_val is the cutoff you'd like to use (so, only retain rows where the diff is less than the quant_val the quantile of the data)
  #returns a linear-scale matrix subset of lin_scale_abundance_mat
  log_v_mat = log2(lin_scale_abundance_mat)
  z_tmp = sapply(1:nrow(log_v_mat), FUN=function(x){diff(range(log_v_mat[x,], na.rm=T))})
  log_vec_for_include = z_tmp<quantile(z_tmp,probs=quant_val) #make argument
  final_mat = lin_scale_abundance_mat[log_vec_for_include,]
  rownames(final_mat) = rownames(lin_scale_abundance_mat)[log_vec_for_include]
  return(lin_scale_abundance_mat[log_vec_for_include,])
}
mapSJ2feature <-
function (SJmat, gtf.file, 
  gtf.key = "gene", gtf.feature = "feature", 
  gtf.orig.col = c("gene_id","gene_name","seqname","start","end"), 
  gtf.col = c("Gene","Symbol","Chr","start","stop"), 
  ikey = 1, isym = 2, ic = 3, i0 = 4, i1 = 5, #chr, start, stop
  pos.col = "Pos", pos.delim ="\\.", delim.sub="_",
  source.me="~burchard/git.R/AbundanceFunctions/readENSgtf.R"
  ) {
  # genomic annotation of splice junction files
  # returns data.table with SJmat values and annotation for genomic features
  #         overlapping start and stop of the splice junction
  # TO DO: generalize for any input genomic interval type
  # gtf.file : file of genome annotation in gtf or gff format
  #           usually best to parse GTF used for data processing
  # gtf.feature : column with biotype AND rowtype information
  # gtf.key : key rowtype to keep, given as value found in feature column
  # gtf.orig.col : column names to keep in extracted gtf.file rows 
  # ikey, isym, ic, i0, i1 : short indices for key feature, gene symbol and 
  #                          SJ position columns chr, start, stop
  # pos.col : column for SJ identifiers comprising chr, start, stop
  # pos.delim : delimiter for splitting SJ identifiers into chr, start, stop
  # delim.sub : replacement for delimiter if found chr field
  # source.me : R file with readENSgtf function


  # imports
  require(data.table)
  require(IRanges)
  source(source.me) # very trusting!
  
  
  # read in genome annotation
  # whole file -- use for biotype annotation
  gtf = readENSgtf(filename=gtf.file)
  message("Genome annotation file ",gtf.file," read in with ",nrow(gtf)," rows")
  feature.gtf = gtf[get(gtf.feature)==gtf.key, mget(gtf.orig.col)]
  names(feature.gtf) = gtf.col
  setorderv(feature.gtf,cols=gtf.col[c(ic,i0,i1)],na.last=TRUE)
  # note setkeyv MUST be run 2nd as setorderv wipes key
  setkeyv(feature.gtf,gtf.col[ikey]) 
  
  
  # collect all SJ positions and assign to features as possible
  
  # collect SJ positions
  SJs_v =  unique(rownames(SJmat))
  # defend against delimiters in chromosome names by keeping last 2 delimiters
  mypattern = paste0('^(.*?',pos.delim,'.*?)(\\.[0-9]+\\.[0-9]+)$')
  idx = grep(mypattern, SJs_v)
  if( length(idx)>0 ){ # some chr have delimiter
    tmp1 = gsub(mypattern, '\\1', SJs_v[idx])
    tmp1 = gsub(pos.delim, delim.sub, tmp1)
    tmp2 = gsub(mypattern, '\\2', SJs_v[idx])
    SJs_v[idx] = paste0(tmp1,tmp2)
  }
  # split out fields
  SJs_dt = data.table(do.call(rbind,
                      lapply(SJs_v,function(x){strsplit(x,pos.delim)[[1]]})))
  names(SJs_dt) = gtf.col[c(ic,i0,i1)]
  # put chr names back as they were
  if( length(idx)>0 ){
    set(SJs_dt, i=as.integer(idx), j=gtf.col[ic], 
        value=gsub(delim.sub,pos.delim,SJs_dt[idx,get(gtf.col[ic])]) )
  }
  # convert start, stop to numeric
  set(SJs_dt, j=gtf.col[i0], value=as.numeric(SJs_dt[,get(gtf.col[i0])]) )
  set(SJs_dt, j=gtf.col[i1], value=as.numeric(SJs_dt[,get(gtf.col[i1])]) )
  setorderv(SJs_dt,cols=gtf.col[c(ic,i0,i1)],na.last=TRUE)
  # SJ start feature; initiate to NA of type character
  set(SJs_dt,j=gtf.col[ikey],value="NA"); SJs_dt[,(gtf.col[ikey])]=NA
  # SJ end feature if different
  set(SJs_dt,j=paste0(gtf.col[ikey],1),value="NA");SJs_dt[,(paste0(gtf.col[ikey],1))]=NA
  # feature assignment: first, split by chromosomes
  for(chr in unique(feature.gtf[,get(gtf.col[ic])]) ){
    gdx = which(feature.gtf[,get(gtf.col[ic])]==chr)
    sjdx = which(SJs_dt[,get(gtf.col[ic])]==chr)
    # find feature containing start of SJ
    #  return vector positions are positions in sjdx; 
    #  values are positions in gdx
    gdx0 = findOverlaps(subject=IRanges(start=feature.gtf[gdx,get(gtf.col[i0])],
                                  end=feature.gtf[gdx,get(gtf.col[i1])]), 
                        query=IRanges(start=SJs_dt[sjdx,get(gtf.col[i0])], 
                                  width=1), select="first")
    sjdx0 = which(!is.na(gdx0)); gdx0 = gdx0[sjdx0]
    # find feature containing end of SJ
    gdx1 = findOverlaps(subject=IRanges(start=feature.gtf[gdx,get(gtf.col[i0])],
                                  end=feature.gtf[gdx,get(gtf.col[i1])]), 
                        query=IRanges(start=SJs_dt[sjdx,get(gtf.col[i1])], 
                                  width=1), select="first")
    sjdx1 = which(!is.na(gdx1)); gdx1 = gdx1[sjdx1]
    # annotate SJs with feature identifiers where mapped
    set(SJs_dt,i=as.integer(sjdx[sjdx0]), j=gtf.col[ikey],
        value=feature.gtf[gdx[gdx0],get(gtf.col[ikey])])
    set(SJs_dt,i=as.integer(sjdx[sjdx0]), j=gtf.col[isym],
        value=feature.gtf[gdx[gdx0],get(gtf.col[isym])])
    set(SJs_dt,i=as.integer(sjdx[sjdx1]), j=paste0(gtf.col[ikey],1),
        value=feature.gtf[gdx[gdx1],get(gtf.col[ikey])])
    set(SJs_dt,i=as.integer(sjdx[sjdx1]), j=paste0(gtf.col[isym],1),
        value=feature.gtf[gdx[gdx1],get(gtf.col[isym])])
  }
  # add column build like SJmat rownames to translation data.table
  # using non-escaped version of delimiter 
  my.delim = gsub("\\\\","",pos.delim)
  set(SJs_dt, j=pos.col, value=SJs_dt[,
      do.call(paste,c(mget(gtf.col[c(ic,i0,i1)]),list(sep=my.delim)))] )
  setkeyv(SJs_dt,pos.col)

  # combine ID translations with data
  SJmat_dt = data.table(SJ_mat, keep.rownames=T)
  names(SJmat_dt)[names(SJmat_dt)=="rn"] = pos.col
  setkeyv(SJmat_dt, pos.col)
  SJs_dt = merge(SJs_dt, SJmat_dt)

  # return data.table
  return(SJs_dt)
  
}
processCELfiles <-
function(pathToCELfiles){
  #' Read in CEL files from a directory argument, using the oligo package first and then, failing that, trying the affy package instead. Requires the relevent array package to be installed from BioConductor before run. Not sure of a good strategy for installing these if they haven't already been. Also want more sophistocated error handling here, but SO recommends against nested tryCatches(), because they don't work as expected:
  #' http://stackoverflow.com/questions/35716394/nested-try-catch-in-r
  #' 
  #' @param pathToCELfiles a string containing the path to the directory containing 
  #' @return A list with two elements. The first element is the affy or oligo object (GeneFeatureSet). The second is a matrix with rownames and colnames that contains non-normalized intensities for each probe in a given array
  #' @examples
  #' path2CELfiles1 = "/Users/fishema/Desktop/CEL_files_for_testing/HTA/"
  #' test1_ls = convertCELtoNonNormMat(path2CELfiles1)
  #' path2CELfiles2 = "/Users/fishema/Desktop/CEL_files_for_testing/Rhesus/"
  #' test2_ls = convertCELtoNonNormMat(path2CELfiles2)
  #' path2CELfiles3 = "/Users/fishema/Desktop/CEL_files_for_testing/SNP6/"
  #' test3_ls = convertCELtoNonNormMat(path2CELfiles3)
  #' test3_mat = test3_ls[[2]]
  #' @export
  
  require(oligo)
  require(affy)
  tryCatch({
    eCELs = list.celfiles(pathToCELfiles,full.names=T)
    Data_obj = read.celfiles(eCELs)
    return_val = list(Data_obj, exprs(Data_obj))
    names(return_val) = c("MicroarrayObject", "NonNormalizedMatrix")
  }, warning = function(war) {
    print(paste("Warning in convertCELtoNonNormMat function:  ",war))
    return_val = NULL
    #return(return_val)
  }, error = function(err) {
    # Might be affy then
    print(paste("Array type not suitable for oligoClasses package. Trying affy package:  ",err))
    Data_obj=ReadAffy(celfile.path=pathToCELfiles)
    non_norm_mat = exprs(Data_obj)
    return_val = list(Data_obj, non_norm_mat)
    names(return_val) = c("MicroarrayObject", "NonNormalizedMatrix")
    #return(return_val)
  }, finally = {
    print(paste("Done with convertCELtoNonNormMat"))
    return(return_val)
  }) # END tryCatch
}
