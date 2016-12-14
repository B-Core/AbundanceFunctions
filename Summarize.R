##############################################################################################
# Functions specific to summarization of abundance data (e.g., microarray probes to probesets)
#
# Summarize_by_some_custom_ID     --- Apply oligo::summarize to arbitrary probe-->probeset mapping
# Remove_mulimappers_and_return_probe_IDs --- Remove multimapping probes
# complex_process_probes          --- Wraps the above and more; needs work
# 
#
##############################################################################################



Summarize_by_some_custom_ID <-
function(oligoFeatureSetObj=NULL, NormWithRownames_mat=NULL, featureID_v=NULL, customID_v=NULL, meth="medianpolish"){
  #' Please refer to the README for comments about this function in its current implementation
  #' Summarizes either a _FeatureSet object from the oligo package or a NormWithRownames_mat.
  #' @description
  #' If the summarization occurs using a _FeatureSet (e.g., ExonFeatureSet, ExpressionFeatureSet, GeneFeatureSet, etc.) object from oligo, summarization is accomplished by applying rma() with background and normalized set to "FALSE". If using NormWithRownames_mat, the function maps the probe IDs in featureID_v to the same rownames in the matrix. If both arguments are supplied, it will return the summarization based on NormWithRownames_mat but will throw a warning if the two summarizations produce different nrows(). Assumes that the featureID_v is ALREADY MAPPED CORRECTLY to the customID_v (i.e., they should be the same length and correspond to one another). customID_v could contain probeset IDs or something else entirely.
  #' @param oligoFeatureSetObj a FeatureSet object generated from oligo's read.celfiles function.
  #' @param NormWithRownames_mat a matrix of normalized data on a linear scale. It's rownames should be feature IDs (e.g., probe IDs for microarrays)
  #' @param featureID_v a vector of strings of IDs of the same type/ilk as rownames(NormWithRownames_mat)
  #' @param customID_v a vector of strings of IDs to which those probes/IDs from featureID_v are to be collapsed/summarized
  #' @param meth a string passed to oligo::summarize that specifies the summarization method.
  #' @return a matrix of nrow <= length(unique(customID_v))
  #' @example 
  #' Summarize_by_some_custom_ID(microarrayProcessing_ls$NonNormalizedMatrix,featureID_v=probeset_df$fid, customID_v=probeset_df$fsetid, meth="medianpolish")
  #' @export
  if (is.null(oligoFeatureSetObj) & is.null(NormWithRownames_mat)){
    stop("_FeatureSet and/or matrix not provided")
  }
  if (!(is.null(oligoFeatureSetObj))){
    tmp = oligo::rma(oligoFeatureSetObj, background=F, normalize=F)
    sData = exprs(tmp)
  }
  if(!(is.null(NormWithRownames_mat))){
    if(!(nrow(NormWithRownames_mat) >= length(unique(customID_v)))){
      stop("The matrix is alredy smaller than the number of IDs you're trying to collapse it to.")
    }
    if(!(length(customID_v) == length(featureID_v))){
      stop("The vectors are not the same length as one another.")
    }
    
    ##NormWithRownames_mat gets logged again by the summarize function, so we bring it back to linear space here##
    NormWithRownames_mat = 2^NormWithRownames_mat
    
    idx2=match(featureID_v,rownames(NormWithRownames_mat)) #position in y where x is
    idx1=which(!is.na(idx2))
    idx2=idx2[idx1]
    ready_for_summarization=matrix(NA, ncol=ncol(NormWithRownames_mat), nrow=length(customID_v))  #nrow=min(nrow(NormWithRownames_mat), length(customID_v))) #probably need to test this when ncol(matrix) is smaller than length(customID_v)
    ready_for_summarization[idx1,]=NormWithRownames_mat[idx2,]
    rownames(ready_for_summarization)=customID_v #[idx1] #should I index by idx1 here? #there will be NAs in here if NormWithRownames_mat doesn't contain all probes in the array 
    colnames(ready_for_summarization) = colnames(NormWithRownames_mat)
    sDataTmp = oligo::summarize(ready_for_summarization, method=meth)
    if(exists("sData")){
      if (nrow(sData) != nrow(sDataTmp)){
        warning("Summarization using oligo's rma function and summarizing based on user-supplied vectors don't agree on matrix size. This likely means that you're using the inappropriate featureID_v and customID_v")
      }
    }
    sData = sDataTmp
    #out_mat = as.data.table(sData)
    #out_dt$probeset_id = rownames(sData)
    #rownames(out_dt) = rownames(sData)
  }
  return(sData)
} #END Summarize_by_some_custom_ID
Remove_mulimappers_and_return_probe_IDs <-
function(vec_of_probe_IDs, vec_of_probeset_IDs){
  #Takes as arguments a vector of probesetID and a vector of probes. Asks how many probesets are using each of the probes in that vector. Returns only those probes that are not used by more than one probeset.
  #Assumes that vec_of_probe_IDs and vec_of_probeset_IDs are ordered the same way and correspond.
  num_probesets_using_ea_probe=ave(vec_of_probeset_IDs,vec_of_probe_IDs, FUN=function(x){length(unique(x[!is.na(x)]))})
  sorted_match_counts<-table(num_probesets_using_ea_probe)[order(table(num_probesets_using_ea_probe),decreasing=T)]#table basically assigns them to bins with different counts of each occurrence
  non_multiple_mappers=vec_of_probe_IDs[!(num_probesets_using_ea_probe>1)]
  return(non_multiple_mappers)
}
complex_process_probes <-
function(probe_mat,probe_ID_vec, customID_v, quantile_to_remove=0.99999){
  ##For arrays like the HTA array, takes a matrix of probe instensities with rownames and colnames and summarizes it, removed extremely variariable probes and summarizes again, removed multi-mappers except in the case where that means removing all probes mapping to a probe set and summarizes again, then takes log ratio to average all for some of these. Returns a list of many different kinds of output.
  #probe_mat is the matrix of feature intensities, with rownames and colnames
  #probe_ID_vec is a vector of probe IDs corresponding EXACTLY to the probe set IDs in customID_v (see below)
  #customID_v is a vector of custom probe set IDs (e.g., for HTA, probes can map to transcript clusters or to exons. Those probe set IDs will be different)
  #quantile_to_remove is essentially the percentile of wildly variable probes you want to remove from the data set before summarization.
  
  summarized_including_all = Summarize_by_some_custom_ID(probe_mat, probe_ID_vec, customID_v)
  summarized_including_all = summarized_including_all[order(probeset_id)]
  crazies_removed = Rm_quantile_of_very_different_intensities(probe_mat,quantile_to_remove)
  probe_mats_ls = list(probe_mat,crazies_removed)
  summarized_crazies_removed = Summarize_by_some_custom_ID(crazies_removed, probe_ID_vec, customID_v)
  summarized_crazies_removed = summarized_crazies_removed[order(probeset_id)]
  #Find multi-mappers
  non_multiple_mapper_probes = Remove_mulimappers_and_return_probe_IDs(probe_ID_vec, customID_v)
  #Remove multi-mappers
  idx_in_arg2_with_no_multi = match(non_multiple_mapper_probes,probe_ID_vec)
  idx1=which(!is.na(idx_in_arg2_with_no_multi)) 
  idx_in_arg2_with_no_multi=idx_in_arg2_with_no_multi[idx1]
  custom_probeset_name_no_multi = customID_v[idx_in_arg2_with_no_multi]
  no_multi_summarized =  Summarize_by_some_custom_ID(crazies_removed, non_multiple_mapper_probes, custom_probeset_name_no_multi)
  #patch those missing from no_multi_summarized with summarized_crazies_removed
  only_multi_idx = !(summarized_crazies_removed$probeset_id %chin% no_multi_summarized$probeset_id)
  rows_missing_from_no_multi_summarized = summarized_crazies_removed[only_multi_idx,]
  #rownames(missing_multi_only_non_norm) = missing_multi_only_non_norm$probeset_id
  merged_no_multi_summarized = rbind(no_multi_summarized, rows_missing_from_no_multi_summarized)
  merged_no_multi_summarized = merged_no_multi_summarized[order(probeset_id)]
  print(sum(merged_no_multi_summarized$probeset_id==summarized_crazies_removed$probeset_id & merged_no_multi_summarized$probeset_id==merged_no_multi_summarized$probeset_id))
  E_av_all_min_multi = rowMeans(merged_no_multi_summarized[,-ncol(merged_no_multi_summarized), with=F], na.rm=TRUE)
  E_lr_av_all_min_multi = merged_no_multi_summarized[,-ncol(merged_no_multi_summarized), with=F]- E_av_all_min_multi
  rownames(E_lr_av_all_min_multi) = merged_no_multi_summarized[[ncol(merged_no_multi_summarized)]]
  E_av_all_all_multi = rowMeans(summarized_crazies_removed[,-ncol(summarized_crazies_removed), with=F], na.rm=TRUE)
  E_lr_av_all_all_multi= summarized_crazies_removed[,-ncol(summarized_crazies_removed), with=F] - E_av_all_all_multi
  rownames(E_lr_av_all_all_multi) = summarized_crazies_removed[[ncol(summarized_crazies_removed)]]
  E_lr_av_all_ls = list(E_lr_av_all_min_multi, E_lr_av_all_all_multi)
  names(E_lr_av_all_ls) = c("E_lr_av_all_min_multi", "E_lr_av_all_all_multi")
  rtrn_ls = list(summarized_including_all,probe_mats_ls,summarized_crazies_removed,merged_no_multi_summarized, E_lr_av_all_ls)
  names(rtrn_ls) = c("summarized_untouched", "probe matrix list", "summarized after extreme probes removed", "summarized after multi-mappers are minimized", "list of log ratio to average all")
  return(rtrn_ls)
}
