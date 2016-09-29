getTopNCellTypesDtVersionWithOutOfRangeReturns_dtWrap = function(listOfListOfDataTables_input,N){
  master_dt = NULL
  listOfListOfDataTables = copy(listOfListOfDataTables_input)
  for(m in 1:length(listOfListOfDataTables)){
    for(s in 1:length(listOfListOfDataTables[[m]])){
      contributing_dt[,query_info := names(listOfListOfDataTables)[m]]
      contributing_dt[,match_info := names(listOfListOfDataTables[[m]][s])]
      for (tpmcol in 1:N){
        contributing_dt [ ,paste0("CellType_",tpmcol) := names(listOfListOfDataTables[[m]][[s]])[tpmcol]]
      }
      #setcolorder(contributing_dt,c("query_info","match_info",names(contributing_dt)[! names(contributing_dt)  %in% c("match_info", "query_info")]))
      master_dt = rbind(master_dt,contributing_dt) #rbindlist(contributing_dt,contributing_dt,use.names = F,fill=F)
      rm(contributing_dt)
    }
  }
  return(master_dt)
}

getSortedDTwithSplitOfFirstCol = function(ff_dt,nameOfCol1="00Annotation"){
  #Function to return a data table of the fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #ff_dt is a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.
  #nameOfCol1 is is a string containing the name() of the column in ff_dt containing the position information. In the ff_dt, note that the position information needs to be parsed out, which this function('s dependency) does automatically. Default is "00Annotation".
  splits_ls = lapply(ff_dt[,get(nameOfCol1)], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3])
  splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
  cNamesTmp = copy(names(ff_dt))
  ff_dt[,chr := gsub('chr','',splits_mat[,1])]
  ff_dt[,chr := gsub('X','25',ff_dt[,chr])]
  ff_dt[,chr := gsub('Y','26',ff_dt[,chr])]
  ff_dt[,Start := splits_mat[,2]]
  ff_dt[,Stop := splits_mat[,3]]
  #ff_dt[,get(nameOfCol1):=NULL]
  #ff_dt[,1,with=F] =NULL
  setcolorder(ff_dt,c("chr", "Start", "Stop",cNamesTmp))
  #coerce chr, Start, and Stop into numeric (the idea being that things that turn into NAs *should* be NAs)
  ff_dt$chr = as.numeric(ff_dt$chr)
  ff_dt$Start = as.numeric(ff_dt$Start)
  ff_dt$Stop = as.numeric(ff_dt$Stop)
  ff_dt_return = ff_dt[order(chr,Start)]
  ff_dt_return[,chr := gsub('26','Y',ff_dt_return[,chr])]
  ff_dt_return[,chr := gsub('25','X',ff_dt_return[,chr])]
  return(ff_dt_return)
}

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

getTopNCellTypesDtVersionWithOutOfRangeReturns = function(ff_dt,chromoNum, range_v, N, colStrWithChromLocations_str="00Annotation") {
  #Takes as input a genome position (chr,Start,Stop), looks for entries in ff_dt with an intersect with the query, and returns the top N tpm columns FOR EACH hit as a separate data.table. 
  #In other words, the output is a LIST of data.tables.
  #If no matching entries are found, the function will search for the nearest entry and return the top N tpm columns of that. If there are two nearest neighbors, both will be returned.
  #ff_dt is a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt. If no data.table is provided, the function will look for the hard-coded path in exacloud (see below).
  #chromoNum is a numeric() chromosome number (nb NOT "chr1") of the query position.
  #range_v is a numeric() vector containing the Start and Stop bp position of the querty (e.g., range_v = c(100,123))
  #N is a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #colStrWithChromLocations_str is a string containing the name() of the column in ff_dt containing the position information. In the ff_dt, note that the position information needs to be parsed out, which this function('s dependency) does automatically. Default is "00Annotation".
  if(is.null(ff_dt)){ #unchecked
    require(data.table)
    ff_dt = fread("grep -v '^#' /home/exacloud/lustre1/CompBio/genomic_resources/fantom5/human/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")
  }
  ffSorted_dt = getSortedDTwithSplitOfFirstCol(ff_dt=copy(ff_dt), nameOfCol1=colStrWithChromLocations_str)
  #check that range_v is valid
  if (length(range_v)!=2){
    print("Too few or two many items in range_v. Should be a vector of length 2")
    return(NULL)
  }else{
    #make sure it goes in correct order
    range_v = range(range_v)
    #Find the first row where Start of our query is >= Start in dt
    chromHits = which(ffSorted_dt[,chr]==as.character(chromoNum))
    startRowChrom = chromHits[1]
    endRowChrom = chromHits[length(chromHits)]
    directHit = which(ffSorted_dt[startRowChrom:endRowChrom,Stop]>=range_v[1] & range_v[2]>= ffSorted_dt[startRowChrom:endRowChrom,Start])
    #I think this has to be a for loop, because I have to sort each row independently, and there's no way to retain colnames under those circumstances
    #topCells = t(apply(ffSorted_dt[startRowChrom:endRowChrom,][directHit,grep('tpm', names(ffSorted_dt)),with=F],1,sort))
    #topCells=sapply(ffSorted_dt[startRowChrom:endRowChrom,][directHit,grep('tpm', names(ffSorted_dt)),with=F],function(x)order(x,decreasing=TRUE))
    #topCells[1:N]
    if(length(directHit)>0){
      matches_ls = list()
      for (dh in 1:length(directHit)){
        topCells_dt = sort(ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],grep('tpm', names(ffSorted_dt)),with=F], decreasing = TRUE)[,1:N,with=F]
        #topCellsName_v = names(topCells_dt)
        matches_ls[[length(matches_ls)+1]] = topCells_dt
        names(matches_ls)[length(matches_ls)] = paste(chromoNum,max(range_v[1],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Start]),min(range_v[2],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Stop]),"matched_region",sep="_")
      }
      return(matches_ls)
    }
    if(length(directHit)==0){
      ##Look for whether the query start is beyond the very last Stop in ffSorted_dt for this chromosome
      if (range_v[1] > ffSorted_dt[startRowChrom:endRowChrom,][nrow(ffSorted_dt[startRowChrom:endRowChrom,]),Stop]){
        nearestUpstreamNeighborIdx = nrow(ffSorted_dt[startRowChrom:endRowChrom,])
        nearestDownstreamNeighborIdx = nearestUpstreamNeighborIdx
      }else{
        ##Look for first occassion where Start in the ffSort_dt is greater than Stop in the query
        nearestDownstreamNeighborIdx = which(ffSorted_dt[startRowChrom:endRowChrom,Start]>range_v[2])[1]
        if (nearestDownstreamNeighborIdx>1){
          nearestUpstreamNeighborIdx = nearestDownstreamNeighborIdx-1
        }else{
          nearestUpstreamNeighborIdx=nearestDownstreamNeighborIdx
        }
      }
      #is the Stop of the nearestUpstreamNeighbor closer to start of the query, or is the Start of nearestDownstreamNeighbor closer to the stop of the query?
      QstopDiffDownStart = ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Start] - range_v[2]
      QstartDiffUpStop= range_v[1] - ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Stop]
      if (QstopDiffDownStart < QstartDiffUpStop){
        newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Start], ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Stop])
        return(getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str))
      }
      if(QstartDiffUpStop < QstopDiffDownStart){
        newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Start], ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Stop])
        return(getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str))
      }
      if(QstopDiffDownStart == QstartDiffUpStop){
        newPositions = range(ffSorted_dt[startRowChrom:endRowChrom,][nearestDownstreamNeighborIdx,Stop],ffSorted_dt[startRowChrom:endRowChrom,][nearestUpstreamNeighborIdx,Start])
        return(getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ff_dt,chromoNum=chromoNum, range_v=newPositions, N=N, colStrWithChromLocations_str=colStrWithChromLocations_str))
      }
    }
  }
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

# ##Create dummy data set and test getTopNCellTypesDtVersionWithOutOfRangeReturns function
# ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0))
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(1,10000), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(31,33), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(1,3), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(10000,10030), N=2, colStrWithChromLocations_str="Annotation")
