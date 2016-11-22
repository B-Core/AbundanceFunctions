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
function(normalized_matrix_with_rownames, feature_ID_vec, custom_ID_vec, meth="medianpolish"){
  #' Summarizes the rows in normalized_matrix_with_rownames on the basis of those rownames that are found in feature_ID_vec.
  #' @description  hasn't been troubleshooted yet for what happens when nrow(normalized_matrix_with_rownames) is smaller than length(custom_ID_vec)
  #' Maps the probe IDs in feature_ID_vec to the same rownames in the matrix. Assumes that the feature_ID_vec is ALREADY MAPPED CORRECTLY to the custom_ID_vec (i.e., they should be the same length and correspond to one another). custom_ID_vec could contain probeset IDs or something else entirely.
  #' @param normalized_matrix_with_rownames is a matrix of normalized data on a linear scale. It's rownames should be feature IDs (e.g., probe IDs for microarrays)
  #' @param feature_ID_vec is a vector of strings of IDs of the same type/ilk as rownames(normalized_matrix_with_rownames)
  #' @param custom_ID_vec is a vector of strings of IDs to which those probes/IDs from feature_ID_vec are to be collapsed/summarized
  #' @return a matrix of nrow length(unique(custom_ID_vec)) with rownames %in% unique(custom_ID_vec)
  #' @examples 
  #' @export
  #browser()
  if(!(nrow(normalized_matrix_with_rownames) == length(custom_ID_vec) & length(custom_ID_vec) == length(feature_ID_vec))){
    stop("The matrix doesn't have the same number of rows as the vectors")
  }
  idx2=match(feature_ID_vec,rownames(normalized_matrix_with_rownames)) #position in y where x is
  idx1=which(!is.na(idx2))
  idx2=idx2[idx1]
  ready_for_summarization=matrix(NA, ncol=ncol(normalized_matrix_with_rownames), nrow=length(custom_ID_vec))  #nrow=min(nrow(normalized_matrix_with_rownames), length(custom_ID_vec))) #probably need to test this when ncol(matrix) is smaller than length(custom_ID_vec)
  ready_for_summarization[idx1,]=normalized_matrix_with_rownames[idx2,]
  rownames(ready_for_summarization)=custom_ID_vec #[idx1] #should I index by idx1 here? #there will be NAs in here if normalized_matrix_with_rownames doesn't contain all probes in the array 
  colnames(ready_for_summarization) = colnames(normalized_matrix_with_rownames)
  sData = oligo::summarize(ready_for_summarization, method=meth)
  #out_mat = as.data.table(sData)
  #out_dt$probeset_id = rownames(sData)
  #rownames(out_dt) = rownames(sData)
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
function(probe_mat,probe_ID_vec, custom_ID_vec, quantile_to_remove=0.99999){
  ##For arrays like the HTA array, takes a matrix of probe instensities with rownames and colnames and summarizes it, removed extremely variariable probes and summarizes again, removed multi-mappers except in the case where that means removing all probes mapping to a probe set and summarizes again, then takes log ratio to average all for some of these. Returns a list of many different kinds of output.
  #probe_mat is the matrix of feature intensities, with rownames and colnames
  #probe_ID_vec is a vector of probe IDs corresponding EXACTLY to the probe set IDs in custom_ID_vec (see below)
  #custom_ID_vec is a vector of custom probe set IDs (e.g., for HTA, probes can map to transcript clusters or to exons. Those probe set IDs will be different)
  #quantile_to_remove is essentially the percentile of wildly variable probes you want to remove from the data set before summarization.
  
  summarized_including_all = Summarize_by_some_custom_ID(probe_mat, probe_ID_vec, custom_ID_vec)
  summarized_including_all = summarized_including_all[order(probeset_id)]
  crazies_removed = Rm_quantile_of_very_different_intensities(probe_mat,quantile_to_remove)
  probe_mats_ls = list(probe_mat,crazies_removed)
  summarized_crazies_removed = Summarize_by_some_custom_ID(crazies_removed, probe_ID_vec, custom_ID_vec)
  summarized_crazies_removed = summarized_crazies_removed[order(probeset_id)]
  #Find multi-mappers
  non_multiple_mapper_probes = Remove_mulimappers_and_return_probe_IDs(probe_ID_vec, custom_ID_vec)
  #Remove multi-mappers
  idx_in_arg2_with_no_multi = match(non_multiple_mapper_probes,probe_ID_vec)
  idx1=which(!is.na(idx_in_arg2_with_no_multi)) 
  idx_in_arg2_with_no_multi=idx_in_arg2_with_no_multi[idx1]
  custom_probeset_name_no_multi = custom_ID_vec[idx_in_arg2_with_no_multi]
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
