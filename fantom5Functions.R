getTopXCellTypes = function(ff_dt=NULL, CHROMONUM, range_v, x){
  if(is.null(ff_dt)){
    require(data.table)
    ff_dt = fread("grep -v '^#' /home/exacloud/lustre1/CompBio/genomic_resources/fantom5/human/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")
  }
  #check that range_v is valid
  if (length(range_v)!=2){
    print("Too few or two many items in range_v. Should be a vector of length 2")
    return(NULL)
  }else{
    #make sure it goes in correct order
    range_v = range(range_v)
    #parse the first column of ff_dt, which has format of, chr10:100150986..100150988,+
    annotationParsed_ls = strsplit(ff_dt[["00Annotation"]], split='chr|:|\\.\\.|,')
    #find first list element with 5
    FIRST=1
    while (length(annotationParsed_ls[[x]])<5) {
      FIRST = FIRST+1
    }
    matches_ls = list()
    for (locCounter in FIRST:length(annotationParsed_ls)){
      if(range_v[1] <= annotationParsed_ls[[locCounter]][4] & range_v[2] >= annotationParsed_ls[[locCounter]][3] & CHROMONUM==annotationParsed_ls[[locCounter]][2]){
        inRangeEntry_ls = sort(as.matrix(ff_dt)[locCounter,grep('tpm', colnames(ff_dt))], decreasing=T, index.return=T) #this is a super time consuming step
        matches_ls = list(matches_ls, inRangeEntry_ls[[1]][1:x])
        #return(inRangeEntry_ls[[1]][1:x])
      }
    }
    return(matches_ls)
  }
}