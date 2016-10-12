#####################################################################################################
# Set of functions useful for cell specificity analysis using the Fantom5 phase1 database.
#
# FFTopNCellTypesWrapper    --- Wrap the output of getTopNCellTypesDtVersionWithOutOfRangeReturns over a data.table of queries instead of one at a time and convert the output of that function into a data.table of cell type/sample names instead of the list of list of data.tables using the ConvertFFlsOflsOfDtToDt function.
# ConvertFFlsOflsOfDtToDt   --- Generate a data.table with N+2 columns. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt).
# getSortedDTwithSplitOfFirstCol       --- Generate and return a data table of the fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
# generateLg2RatioToAveAllMat --- Generate and return a matrix containing log2-scale ratios to average all for the numerical portion of the fantom5 data
# fetchTssRangesWithXFoldOverUnderExpressionInFFOntologyID    --- C
# getTopNCellTypesDtVersionWithOutOfRangeReturns      --- P
# getTopXCellTypes         --- P
#
# Author: Mark Fisher
# Started - October, 2016
#####################################################################################################


FFTopNCellTypesWrapper = function(ff_dt, query_dt, ChrColName, StartColName, StopColName, N, savePath, colStrWithChromLocations_str="00Annotation"){ #, queryDtENSRidColName="ENSRid", queryDtRegDescriptionColName="gwasEnsemblRegDescription", queryDtVarLenAcrossENSRidColName="VarLenAcrossENSRid_x"){
  #' Wrap the output of getTopNCellTypesDtVersionWithOutOfRangeReturns over a data.table of queries instead of one at a time and convert the output of that function into a data.table of cell type/sample names instead of the list of list of data.tables using the ConvertFFlsOflsOfDtToDt function.
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user. If no data.table is provided, the function will look for the hard-coded path in exacloud.
  #' @param query_dt a data.table with at least three columns: one corresponding to a chromosome number, one corresponding to the start position, and one correpsonding to the stop position. Rows populated with positions (e.g., regulatory element positions) of interest to the user.
  #' @param ChrColName a string containing the name() of the column in query_dt containing the chromosome position information.
  #' @param StartColName a string containing the name() of the column in query_dt containing the start position information.
  #' @param StopColName a string containing the name() of the column in query_dt containing the stop position information.
  #' @param N a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #' @param savePath a string that specifies the destination directory for the two return items (a data table and list of list of data tables)
  #' @return A list containing two elements: 1) a data table with ncol() N+2. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt). 2) A list of list of data tables. The length() of the list will be nrow(query_dt). Each element of the list will contain a list of data tables, where each element of this list will represent one match (as previously mentioned, there will always be at least one "match"--direct or nearest--and sometimes more than one match). Each match will correspond to a data.table with nrow()=1 and ncol()=N. The names() of the data.table will correspond to the cell types/samples in ff_dt, and the entries in the table will be tpm counts from ff_dt, sorted from highest to lowest. If a chromosome number in the query data table doesn't match a chromosome number in ff_dt, it simply gets skipped.
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0))
  #' q_dt = data.table(chrP = c(10,10,10,10,11), startP = c(1,31,1,10000,1), stopP=c(10000,33,3,10030,100))
  #' output_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation")
  
  
  #If ff_dt was not provided, try to load it from its home in exacloud #Attn hardcoded
  if(is.null(ff_dt)){ #unchecked
    require(data.table)
    ff_dt = fread("grep -v '^#' /home/exacloud/lustre1/CompBio/genomic_resources/fantom5/human/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")
  }
  
  #Check that ff_dt and query_dt are data.tables. If not, try to coerce them to data.tables
  if(! (is.data.table(ff_dt) & is.data.table(query_dt))){
    ff_dt = as.data.table(ff_dt)
    query_dt = as.data.table(query_dt)
  }
  
  #Initiate and populate the list of list of data.tables
  FFtopNcells_lsOfdt_ls = list()
  for(r in 1:nrow(query_dt)){
    if(r%%50==0) {print(r)} #print r every 50th iteration of the loop to reassure an eager analyst
    FFtopNcells_lsOfdt_ls[[r]] = getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ff_dt,chromoNum=as.numeric(query_dt[r,get(ChrColName)]), range_v=c(as.numeric(query_dt[r,get(StartColName)]),as.numeric(query_dt[r,get(StopColName)])), N=N, colStrWithChromLocations_str=colStrWithChromLocations_str)
    
    #if the r-th element of the list was successfully populated (it wouldn't be if, say, there was no chromosome number in ff_dt that matched the query), let's go ahead and give that element a name
    if(length(FFtopNcells_lsOfdt_ls) ==r){
      names(FFtopNcells_lsOfdt_ls)[r] = paste(query_dt[r,get(ChrColName)],query_dt[r,get(StartColName)],query_dt[r,get(StopColName)],"query",sep="_")
    }
  }
  
  #Save the list of list of data.tables to a text file. Not for client viewing, but useful for analyst to check.
  sink(file=paste0(savePath,"FFtop",N,"cells_lsOfdt_ls.txt"))
  print(FFtopNcells_lsOfdt_ls)
  sink()
  
  #Convert the list of list of data.tables to one data.table more amenable to downstream analyses
  Converted_dt=ConvertFFlsOflsOfDtToDt(FFtopNcells_lsOfdt_ls,N)
  write.table(Converted_dt, file=paste0(savePath, "Fantom5Top",N,"CellTypes.txt"), sep="\t", row.names = F)
  
  #Function will return both the data.table Converted_dt and the list of list of data.tables.
  output_ls = list(Converted_dt, FFtopNcells_lsOfdt_ls)
  names(output_ls) = c(paste0("Data table of query and match output for top ",N," cell types" ), "List of list of data tables output containing actual tpms")
  return(output_ls)
} #end FFTopNCellTypesWrapper

ConvertFFlsOflsOfDtToDt = function(listOfListOfDataTables_input,N){
  #' Generate a data.table with N+2 columns. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt).
  #' 
  #' @param listOfListOfDataTables_input a list of list of data.tables generated by the repeated calls to the getTopNCellTypesDtVersionWithOutOfRangeReturns function in the FFTopNCellTypesWrapper function (i.e., the FFtopNcells_lsOfdt_ls object from the wrapper function).
  #' @param N a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #' @return A list containing two elements: 1) a data table with ncol() N+2. The first column ("query_info") contains the chromosome, start, and stop positions of the query, separated by underscore delimiters. The second column ("match_info") contains the chromosome, start, and stop position of the INTERSECT region of the query position and any match within ff_dt. It also includes the original column name of its match in ff_dt following an "OriginalColName" string and an underscore delimiter. Failing a match, the match_info column will report the entire position of the nearest ff_dt TSS. If there are two nearest ff_dt TSSs, two separate matches will be reported, each reporting the entire position of its ff_dt TSS. Note that the query position may intersect more than one ff_dt TSS. In this case, each intersect will be reported as a separate row in the output. The output data.table will therefore have nrow() >= nrow(query_dt). 2) A list of list of data tables. The length() of the list will be nrow(query_dt). Each element of the list will contain a list of data tables, where each element of this list will represent one match (as previously mentioned, there will always be at least one "match"--direct or nearest--and sometimes more than one match). Each match will correspond to a data.table with nrow()=1 and ncol()=N. The names() of the data.table will correspond to the cell types/samples in ff_dt, and the entries in the table will be tpm counts from ff_dt, sorted from highest to lowest. If a chromosome number in the query data table doesn't match a chromosome number in ff_dt, it simply gets skipped.
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' q_dt = data.table(chrP = c(10,10,10,10,11), startP = c(1,31,1,10000,1), stopP=c(10000,33,3,10030,100))
  #' N = 2
  #' output_ls = FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = q_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation")
  #' listOfListOfDataTables_ls = output_ls[[2]]
  #' Converted_dt = ConvertFFlsOflsOfDtToDt(listOfListOfDataTables_ls,N)
   
  #Start with empty data.table
  master_dt = NULL
  
  #Fresh copy of listOfListOfDataTables_input to avoid manipulating it by reference
  listOfListOfDataTables = copy(listOfListOfDataTables_input)
  
  #Loop through each data table in the list of list of data tables
  for(m in 1:length(listOfListOfDataTables)){
    for(s in 1:length(listOfListOfDataTables[[m]])){
      
      #And create a "row" which has to be a data.table for rbind to work containing a query_info column
      contributing_dt = data.table(query_info = names(listOfListOfDataTables)[m] )
      
      #Add a match_info column to that same "row"
      contributing_dt[,match_info := names(listOfListOfDataTables[[m]][s])]
      
      for (tpmcol in 1:N){
        #And then add a new column for each of the N most-higlhy-ranked tpm counts, and populate with the CAGE peak IDs corresponding to these highly-ranked counts
        contributing_dt [ ,paste0("CellType_",tpmcol) := names(listOfListOfDataTables[[m]][[s]])[tpmcol]]
      }
      
      #Tack the "row" onto the master data.table
      master_dt = rbind(master_dt,contributing_dt)
      
      #And reset the "row"
      rm(contributing_dt)
    }
  }
  
  #When you're done looping through the list of lists of data.tables, return the resulting master data.table that's been generated
  return(master_dt)
} #END ConvertFFlsOflsOfDtToDt

getSortedDTwithSplitOfFirstCol = function(ff_dt,nameOfCol1="00Annotation"){
  #' Generate and return a data table of the fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @param nameOfCol1 a string matching the column in ff_dt that contains the TSS ID information (e.g., "chr10:10..20,-"). This function will parse out the chromosome, start, and stop.
  #' @return A data.table of fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' sorted_dt = getSortedDTwithSplitOfFirstCol(ff_dt = ffDummy_dt,nameOfCol1="Annotation")
  
  
  #Parse the fantom5 annotation (TSS identifier) column
  splits_ls = lapply(ff_dt[,get(nameOfCol1)], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3])
  splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
  
  #Get a placeholder vector of column names of ff_dt to later arrange the column order. You'll have to get a deep copy of this because data.table does everything by reference
  cNamesTmp = copy(names(ff_dt))
  
  #Remove chr string and substitute X and Y chromosome for 23 and 24, respectively, temporarily while the rows get sorted.
  ff_dt[,chr := gsub('chr','',splits_mat[,1])]
  ff_dt[,chr := gsub('X','23',ff_dt[,chr])]
  ff_dt[,chr := gsub('Y','24',ff_dt[,chr])]
  
  #Populate start and stop columns with the parsed start and stop locations. Re-arrange the column order.
  ff_dt[,Start := splits_mat[,2]]
  ff_dt[,Stop := splits_mat[,3]]
  setcolorder(ff_dt,c("chr", "Start", "Stop",cNamesTmp))
  
  #coerce chr, Start, and Stop into numeric with the idea being that things that turn into NAs *should* be NAs (e.g., see the first few rows of ff_dt)
  ff_dt$chr = as.numeric(ff_dt$chr)
  ff_dt$Start = as.numeric(ff_dt$Start)
  ff_dt$Stop = as.numeric(ff_dt$Stop)
  
  #Sort ff_dt by chromosome and then by start position
  ff_dt_return = ff_dt[order(chr,Start)]
  
  #Restore X and Y chromosomes to their former glory (i.e., their previous chromosome identifiers)
  ff_dt_return[,chr := gsub('23','X',ff_dt_return[,chr])]
  ff_dt_return[,chr := gsub('24','Y',ff_dt_return[,chr])]
  return(ff_dt_return)
} #END getSortedDTwithSplitOfFirstCol

generateLg2RatioToAveAllMat = function (ff_DtOrMat){
  #' Generates and returns a matrix containing log2-scale ratios to average all for the numerical portion of the fantom5 data
  #' 
  #' @param ff_DtOrMat a data.table, matrix, or data.frame of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @return A matrix of log2-scale ratios to average all for the numerical portion of the fantom5 data. Zeroes in the original tpm count matrix are replaced with the smallest non-zero value within the same CAGE peak ID to avoid log2(0)-related errors/issues.
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  #' lgRatio_mat = generateLg2RatioToAveAllMat(ffDummy_dt)
  
  #Coerce ff_mat to matrix if not already and possible
  ff_mat = as.matrix(ff_DtOrMat)
  
  #Subset ff_mat by the portion that is numeric (the tpm columns), coerce to numeric, check for 0s,NAs???? and replace with col-wise non-zero min.
  numericPartOfFf_mat = data.matrix(ff_mat[,grep('tpm', colnames(ff_mat))])
  numericPartOfFf_mat = apply(numericPartOfFf_mat, 2, function(x) as.numeric(x))
  numericPartOfFf_mat_test = apply(numericPartOfFf_mat, 2, function(x) replace(x,x %in% c(0), min(x[x>0]))) #include NA?
  #min(numericPartOfFf_mat[numericPartOfFf_mat > 0],na.rm=T)
  
  #log2 transform the tpm count containing columns
  lg2FF_mat = log2(numericPartOfFf_mat_test)
  
  #take rowmeans of the same columns
  rMeans = rowMeans(lg2FF_mat)
  
  #find ratio to average for each of the columns (a difference in log2 scale)
  lg2Ratio_mat= lg2FF_mat - rMeans
  
  return(lg2Ratio_mat)
}

fetchTssRangesWithXFoldOverUnderExpressionInFFOntologyID = function(ffID_v, ff_mat, x, onlyOverExprs=T){
  #' For each fantom5 ontology ID or query string of interest in ffID_v, generates a matrix of just entries with log2-fold ratios greater than x (or abs(log2(ratio)>x)), decorated with chromosome, start position, stop position, entrezID, and description. If no entries match at the >x (or abs()>x) threshold, populate that matrix with NULL. Assembles these matrices into a list and returns the list.
  #' 
  #' @param ffID_v a vector containing keywords or -even better- fantom5 ontology IDs by which to narrow the search for 
  #' @param ff_mat a matrix
  #' @param onlyOverExprs a boolean specifying whether the user wants to inspect only over-expressed
  #' @param x a number representing the desired log2-fold over- and possibly under-expression desired
  #' @return a list of matrices, one matrix per list element, of just entries with log2-fold ratios greater than x (or abs(log2(ratio)>x)), decorated with chromosome, start position, stop position, entrezID, and description.
  #' 
  #' @examples
  #' #Attn add entrezID to 5th column of ffDummy_mat
  #' ffDummy_mat = as.matrix(data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3), entrezID=c("ID1","ID2","ID3","ID4","ID5","ID6","ID7","ID8","ID9","ID10","ID11"),description=c("Description1","Description2","Description3","Description4","Description5","Description6","Description7","Description8","Description9","Description10","Description11")))
  #' ffID_v1 = c("Two", "Three")
  #' x_var = 3
  #' hits_ls = fetchTssRangesWithXFoldOverUnderExpressionInFFOntologyID(ffID_v=ffID_v1, ff_mat=ffDummy_mat, x=x_var, onlyOverExprs=T)
  
  lg2Ratio_mat = generateLg2RatioToAveAllMat(ff_mat)
  
  #search for those rows with >abs(x) fc for one of any columns containing ids in ffID_v and return those log ratios, along with the chr, start, and stop
  hit_ls = list()
  
  #loop through the fantom5 ontology IDs or query strings and if you find entries with log2-scale ratios >x, assemble a matrix of just these entries, decorated with chromosome, start position, stop position, entrezID, and description
  for (ffID in ffID_v){
    if(onlyOverExprs==TRUE){
      HighExprsInColOfInterestLogic_v = lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))] > x
    }else{
      HighExprsInColOfInterestLogic_v = abs(lg2Ratio_mat[,grep(ffID, colnames(lg2Ratio_mat))]) > abs(x)
    }
    
    #If there were TSSs in our cell type/sample of interest >x
    if(sum(HighExprsInColOfInterestLogic_v, na.rm=T)>0){
      #create a matrix of that cell type sample column
      logRatio_mat = lg2Ratio_mat[HighExprsInColOfInterestLogic_v,grep(ffID, colnames(lg2Ratio_mat)),drop=F]
      
      #parse the identifier column of ff_mat
      splits_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v,1,drop=F], function(x) unlist(strsplit(x,':|\\.+|,'))[1:3]) #Attn don't use number
      splits_mat = matrix(unlist(splits_ls), nrow=length(splits_ls), byrow = T)
      
      entrezIDsplit_ls = lapply(ff_mat[HighExprsInColOfInterestLogic_v, 5, drop=F], function(x) unlist(strsplit(x,':'))[2]) #Attn don't use number
      entrezID_mat = matrix(unlist(entrezIDsplit_ls), nrow=length(entrezIDsplit_ls), byrow = T)
      description_mat = ff_mat[HighExprsInColOfInterestLogic_v, "description", drop=F]
      logRatioJoin_mat = cbind(splits_mat,entrezID_mat,description_mat, logRatio_mat)
      colnames(logRatioJoin_mat) = c("Chr","Start","Stop","EntrezID","Description", colnames(logRatio_mat))
      hit_ls[[length(hit_ls)+1]] = logRatioJoin_mat=logRatioJoin_mat
      names(hit_ls)[length(hit_ls)] = colnames(lg2Ratio_mat)[grep(ffID, colnames(lg2Ratio_mat))]
    }
  }
  return(hit_ls)
} #END fetchTssRangesWithXFoldOverUnderExpressionInFFOntologyID

getTopNCellTypesDtVersionWithOutOfRangeReturns = function(ff_dt,chromoNum, range_v, N, colStrWithChromLocations_str="00Annotation") {
  #' LEFT OFF HERE
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @param nameOfCol1 a string matching the column in ff_dt that contains the TSS ID information (e.g., "chr10:10..20,-"). This function will parse out the chromosome, start, and stop.
  #' @return A data.table of fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
  
  #Takes as input a genome position (chr,Start,Stop), looks for entries in ff_dt with an intersect with the query, and returns the top N tpm columns FOR EACH hit as a separate data.table. 
  #In other words, the output is a LIST of data.tables.
  #If no matching entries are found, the function will search for the nearest entry and return the top N tpm columns of that. If there are two nearest neighbors, both will be returned.
  #ff_dt is a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt. If no data.table is provided, the function will look for the hard-coded path in exacloud (see below).
  #chromoNum is a numeric() chromosome number (nb NOT "chr1") of the query position.
  #range_v is a numeric() vector containing the Start and Stop bp position of the query (e.g., range_v = c(100,123))
  #N is a numeric() argument the retreives the top N largest tpm-containing columns in a given row of ff_dt
  #colStrWithChromLocations_str is a string containing the name() of the column in ff_dt containing the position information. In the ff_dt, note that the position information needs to be parsed out, which this function('s dependency) does automatically. Default is "00Annotation".
  
  #Check that N isn't greater than the number of cell types/samples available
  if(N > length(grep('tpm', names(ff_dt)))){
    return(NULL)
  }else{
    
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
      #make sure the range goes in correct order
      range_v = range(range_v)
      
      #Find the first row where Start of our query is >= Start in dt
      chromHits = which(ffSorted_dt[,chr]==as.character(chromoNum))
      if(length(chromHits)<1){
        print(paste0("Invalid or non-matching chromosome number: ",chromoNum))
        #browse()
        return(NULL)
      }else{
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
            names(matches_ls)[length(matches_ls)] = paste(chromoNum,max(range_v[1],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Start]),min(range_v[2],ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],Stop]),"matched_region_OriginalColName",ffSorted_dt[startRowChrom:endRowChrom,][directHit[dh],][[colStrWithChromLocations_str]],sep="_")
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
  }
}

getTopXCellTypes = function(ff_mat=NULL, chromoNum, range_v, x, colStrWithChromLocations_str){
  #' Generate and return a data table of the fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' 
  #' @param ff_dt a data table of hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt or any version or subset thereof as defined by user.
  #' @param nameOfCol1 a string matching the column in ff_dt that contains the TSS ID information (e.g., "chr10:10..20,-"). This function will parse out the chromosome, start, and stop.
  #' @return A data.table of fantom5 hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt data SORTED numerically by numeric() chromosome number and THEN by Start position
  #' @examples
  #' ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
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
# ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-", "chr12:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86,0), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0,3),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0,3))
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(1,10000), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(31,33), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(1,3), N=2, colStrWithChromLocations_str="Annotation")
# getTopNCellTypesDtVersionWithOutOfRangeReturns(ff_dt=ffDummy_dt,chromoNum=10, range_v=c(10000,10030), N=2, colStrWithChromLocations_str="Annotation")
# ffDummy_dt = data.table(Annotation=c("chr10:10..20,-", "chr10:25..30,-","chr10:35..100,-","chr10:106..205,-","chr10:223..250,-","chr10:269..478,-","chr10:699..1001,-","chr10:2000..2210,-","chr10:2300..2500,-","chr10:2678..5678,-"),tpmOne=c(0,0,0.213,1,1.2,0.5,0.7,0.9,0.8,0.86), tpmTwo=c(100,1000,1001,1500,900,877,1212,1232,1312,0),tpmThree=c(0.2138595,0,0,0,0,0,0.6415786,0,0,0))
# query_dt = data.table(chrP = c(10,10,10,10,11), startP = c(1,31,1,10000,1), stopP=c(10000,33,3,10030,100))
# FFTopNCellTypesWrapper(ff_dt = ffDummy_dt, query_dt = query_dt, ChrColName = "chrP", StartColName ="startP", StopColName = "stopP", N = 2, savePath = "~/Desktop/", colStrWithChromLocations_str="Annotation")

# numericPartOfFf_mat = data.matrix(ffDummy_dt[,grep('tpm', colnames(ffDummy_dt)),with=F])
# ffRanks_mat = matrix(nrow=(nrow(numericPartOfFf_mat)), ncol=length(grep('tpm', colnames(ffDummy_dt))))
# colnames(ffRanks_mat) = colnames(ffDummy_dt)[grep('tpm', colnames(ffDummy_dt))]
# for (r in 1:nrow(numericPartOfFf_mat)){
#   ffRanks_mat[r,] = rank(-numericPartOfFf_mat[r,], na.last = T)
# }
# rankCount_mat = matrix(nrow=1, ncol=length(grep('tpm', colnames(ffDummy_dt))))
# colnames(rankCount_mat) = colnames(ffRanks_mat)
# N=2
# for (cNum in 1:ncol(ffRanks_mat)){
#   rankCount_mat[1,cNum] = sum(ffRanks_mat[,cNum]<=N)
#   #colnames(rankCount_mat)[cNum] = colnames(ffRanks_mat)[cNum]
# }
 
 #ffDummy_dt[,sum(test_v%in%.I), by= .I] #1:nrow(ffDummy_dt)]
 # test_v = c(0,0,0.86)
 # which(ffDummy_dt[, rowSums(.SD), .SDcols=2:ncol(ffDummy_dt)] == sum(test_v))
 #ffDummy_dt[, if (all(.SD %in% test_v)) .SD, Annotation]
