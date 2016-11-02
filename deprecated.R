###############################################################################
# Parking area for AbundanceFunctions not currently in use.
#
# norm_matrix  --- Apply bias reduction functions to abundance data
#                  <== Removed from RNAseq_abundance_functions.R
#                  ==> Replaced by normMatrix in processData.R
#                      15-Oct-2016
# make_heatmap_versatile --- Makes a heat map using the NMF aheatmap function, using only some columns
#                  <== Removed from abundance_functions.R
#                  ==> Replaced by makeHeatmap currently in Proteomics.Tools.x1.R
#                      28-Oct-2016
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

make_heatmap_versatile = function(matrix_with_samples_for_sig_IDs, sample_string_vec, gsub_remove_str='', annotation_mat=NULL, annotate_with_gene_names = F, ID_colname=NULL, Symbol_colname=NULL, save_image_file=F, string_to_lead_file_name_with = "", c.lim=NULL, annCol=NULL, main=NULL, reso=600){
  #Makes a heat map using the NMF aheatmap function, using only the columns of a matrix in sample_string_vec
  #matrix_with_samples_for_sig_IDs is a matrix containing all rows with with you want to make a heatmap. It can contain some columns you don't want, because you specify which columns you want in sample_string_vec. Assumes that matrix_with_samples_for_sig_IDs already has rownames of ID (e.g., transcript cluster ID).
  #sample_string_vec is a vector of column names that you want to use to contstruct the heatmap. The most permissive thing you can do is sample_string_vec=colnames(matrix_with_samples_for_sig_IDs).
  #gsub_remove_str: let's say that your column names are very long, and you don't want to label the columns of your heat map with them. You would provide a regular expression pattern in gsub_remove_str, and it will REMOVE things that match that pattern when naming the columns in the heat map.
  #annotation_mat is a matrix containing at the very least a column of IDs (can be longer than rownames(matrix_with_samples_for_sig_IDs)), and a column of something else (e.g., gene symbols) that match those IDs. To be used in case, say, you want to annotate the rows of the heat map with gene symbols instead of transcript cluster IDs.
  #annotate_with_gene_names a boolean deciding whether you want to annotate the heat map using annotation_mat
  #ID_colname is the column name in annotation_mat that contains the IDs of the same type as in rownames(matrix_with_samples_for_sig_IDs), e.g. transcript cluster ID
  #Symbol_colname is the column name in annotation_mat that contains the values (e.g, gene symbols) that the user actually wants to annotate the heat map with
  #save_image_file is a boolean specifying whether 600 dpi png files should be saved or not
  #string_to_lead_file_name_with is a string the represents the first part of the file name to be saved. Can also be used to specify the destination path of the save file if not current working directory.
  #c.lim are the limits in which most of the color variation of the heat map will occur. The min. and max values will also be included in the heat map color range, but everything between them and the c.lim range will be the darkest two colors.
  #Needs Ann_col argument
  #Needs color list arguments for the different tracks
  if(is.null(c.lim)){
    c.lim = quantile(matrix_with_samples_for_sig_IDs[,sample_string_vec], probs=.95,na.rm=T)
  }
  subset_of_mat_of_interest = matrix_with_samples_for_sig_IDs[,sample_string_vec]
  class(subset_of_mat_of_interest) = "numeric"
  if(is.null(annCol)){
    annCol=gsub(gsub_remove_str, '',sample_string_vec)
  }
  #c.lim = quantile(subset_of_mat_of_interest, probs=.95,na.rm=T)
  if(class(subset_of_mat_of_interest) != "matrix"){
    subset_of_mat_of_interest = t(as.matrix(subset_of_mat_of_interest))
  }
  if(ncol(subset_of_mat_of_interest)>2 & nrow(subset_of_mat_of_interest)>1){
    if(annotate_with_gene_names ==F){
      tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=annCol,labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
      #print(tmp)
      indices_in_heat_map_order = rev(tmp$rowInd)
      ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
      ordered_actual_values = subset_of_mat_of_interest[indices_in_heat_map_order,]
      if(save_image_file==T){
        png(filename=paste0(string_to_lead_file_name_with,"_heatmap.png"),width=5,height=5.4,units="in",res=reso)
        #print(tmp)
        tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=annCol,labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = rep(" ", nrow(subset_of_mat_of_interest)))
        dev.off()
      }
      return(list(ordered_probeset_ids, ordered_actual_values,tmp))
    } else{
      #then annotate_with_gene_names =T
      ##Find the gene names that correpond to rownames of subset_of_mat_of_interest
      gn_names = vec_of_gene_symbols_given_vec_of_IDs(vec_of_IDs=rownames(subset_of_mat_of_interest), master_mat=annotation_mat, gene_symbol_colname_in_master=Symbol_colname, ID_colname_in_master=ID_colname)
      tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=annCol,labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = gn_names)
      #print(tmp)
      indices_in_heat_map_order = rev(tmp$rowInd)
      col_indices_in_heat_map_order = tmp$colInd
      ordered_probeset_ids = rownames(subset_of_mat_of_interest) [indices_in_heat_map_order]
      ordered_gene_ids = vec_of_gene_symbols_given_vec_of_IDs(vec_of_IDs=ordered_probeset_ids, master_mat=annotation_mat, gene_symbol_colname_in_master=Symbol_colname, ID_colname_in_master=ID_colname)
      ordered_actual_values = subset_of_mat_of_interest[indices_in_heat_map_order,col_indices_in_heat_map_order]
      if(save_image_file==T){
        png(filename=paste0(string_to_lead_file_name_with,"_heatmap.png"),width=5,height=5.4,units="in",res=reso)
        tmp=aheatmap(subset_of_mat_of_interest,color=paste0("-PiYG:", nrow(matrix_with_samples_for_sig_IDs)), breaks=c(min(subset_of_mat_of_interest),seq(from=-c.lim,to=c.lim, by=c.lim/((nrow(matrix_with_samples_for_sig_IDs)-2)/2)), max(subset_of_mat_of_interest)),annCol=annCol,labCol=gsub(gsub_remove_str, '',sample_string_vec),labRow = gn_names)
        dev.off()
      }
      return(list(ordered_probeset_ids, ordered_gene_ids, ordered_actual_values,tmp))
    }
  } else{
    return("")
  }
}

convertEnsemblToBed = function(Ensembl_mat, justThreeCols=F){
  chromName = gsub('^ *', '', Ensembl_mat[,1])
  chromName = gsub('^chr', '', Ensembl_mat[,1])
  start = as.numeric(Ensembl_mat[,2])-1
  end = as.numeric(Ensembl_mat[,3])
  if(justThreeCols == T){
    final_mat = cbind(chromName, start, end)
  }else{
    final_mat = cbind(chromName, start, end, Ensembl_mat[,4:ncol(Ensembl_mat)])
  }
  return(final_mat)
}
