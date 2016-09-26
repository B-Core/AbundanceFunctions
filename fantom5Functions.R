fetchTssRangesWithXFoldOverUnderExpressionInFFOntologyID = function(ffID_v, ff_mat, x, onlyUpReg=T){
  #check the numeric section of the ff_mat for 0s,NAs???? and replace with col-wise non-zero min.
  numericPartOfFf_mat = data.matrix(ff_mat[,grep('tpm', colnames(ff_mat))])
  numericPartOfFf_mat = apply(numericPartOfFf_mat, 2, function(x) as.numeric(x))
  numericPartOfFf_mat_test = apply(numericPartOfFf_mat, 2, function(x) replace(x,x %in% c(0), min(x[x>0]))) #include NA?
  #min(numericPartOfFf_mat[numericPartOfFf_mat > 0],na.rm=T)
  
  #log2 transform the tpm count containing columns
  lg2FF_mat = log2(numericPartOfFf_mat_test)#, na.omit=T)
  #take rowmeans of the same columns
  rMeans = rowMeans(lg2FF_mat,na.rm=T)
  #find ratio to average for each of the columns (a difference in log2 scale)
  lg2Ratio_mat= lg2FF_mat - rMeans #can I do this?
  #search for those rows with >abs(x) fc for one of any columns containing ids in ffID_v and return those log ratios, along with the chr, start, and stop
  hit_ls = list()
  for (ffID in ffID_v){
    if(onlyUpReg==TRUE){
      HighExprsInColOfInterestLogic_v = lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))] > x
    }else{
      HighExprsInColOfInterestLogic_v = abs(lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))]) > abs(x)
    }
    if(sum(HighExprsInColOfInterestLogic_v, na.rm=T)>0){
      logRatio_mat = lg2Ratio_mat[HighExprsInColOfInterestLogic_v,grep(ffID, colnames(lg2Ratio_mat)),drop=F]
      splits_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v,1,drop=F], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3])
      splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
      entrezIDsplit_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v, 5, drop=F], function(x) unlist(strsplit(x,':'))[2])
      entrezID_mat = matrix(unlist(entrezIDsplit_ls), nrow=length(entrezIDsplit_ls), byrow = T)
      description_mat = ff_mat[HighExprsInColOfInterestLogic_v, "description", drop=F]
      logRatioJoin_mat = cbind(splits_mat,entrezID_mat,description_mat, logRatio_mat)
      colnames(logRatioJoin_mat) = c("Chr","Start","Stop","EntrezID","Description", colnames(logRatio_mat))
      hit_ls[[length(hit_ls)+1]] = logRatioJoin_mat=logRatioJoin_mat
      names(hit_ls)[length(hit_ls)] = colnames(lg2Ratio_mat)[grep(ffID, colnames(lg2Ratio_mat))]
    }
  }
  return(hit_ls)
  # logic_mat = lg2Ratio_mat > abs(x)
  # logic_mat[,grep()]
  #colsToKeep = grep(ffID_v[i], colnames(ff_mat))
}

getTopXCellTypes = function(ff_mat=NULL, chromoNum, range_v, x, colStrWithChromLocations_str){
  if(is.null(ff_mat)){
    require(data.table)
    ff_mat = as.matrix(fread("grep -v '^#' /home/exacloud/lustre1/CompBio/genomic_resources/fantom5/human/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt"))
  }
  #check that range_v is valid
  if (length(range_v)!=2){
    print("Too few or two many items in range_v. Should be a vector of length 2")
    return(NULL)
  }else{
    #make sure it goes in correct order
    range_v = range(range_v)
    #parse the first column of ff_mat, which has format of, chr10:100150986..100150988,+
    annotationParsed_ls = strsplit(ff_mat[,colStrWithChromLocations_str], split='chr|:|\\.\\.|,')
    #find first list element with 5
    FIRST=1
    while (length(annotationParsed_ls[[FIRST]])<5) {
      FIRST = FIRST+1
    }
    matches_ls = list()
    for (locCounter in FIRST:length(annotationParsed_ls)){
      if(range_v[1] <= as.numeric(annotationParsed_ls[[locCounter]][4]) & range_v[2] >= as.numeric(annotationParsed_ls[[locCounter]][3]) & chromoNum==annotationParsed_ls[[locCounter]][2]){
        #print(locCounter)
        #browse()
        #inRangeEntry_ls = sort(as.numeric(ff_mat[locCounter,grep('tpm', colnames(ff_mat))]), decreasing=T, index.return=T) #this is a super time consuming step
        inRangeEntry_ls = sort(ff_mat[locCounter,grep('tpm', colnames(ff_mat))], decreasing=T, index.return=T) #this is a super time consuming step
        names(inRangeEntry_ls[[1]]) = colnames(ff_mat)[grep('tpm', colnames(ff_mat))[inRangeEntry_ls[[2]]]]
        matches_ls[[length(matches_ls)+1]] = inRangeEntry_ls[[1]][1:x]
        names(matches_ls)[length(matches_ls)]=paste0(annotationParsed_ls[[locCounter]][2],":", annotationParsed_ls[[locCounter]][3],"-", annotationParsed_ls[[locCounter]][4])
        #return(inRangeEntry_ls[[1]][1:x])
      }#else search nearby for closest tss
    }
    return(matches_ls)
  }
}