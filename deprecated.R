# norm_matrix     --- Please use normMatrix() in instead.

get_pairwise_log_ratio <-
function(cntrl, other_treatment_str, E_treatment_wise_intensities){
  cntrl_col = E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""), colnames(E_treatment_wise_intensities))]
  other_treatment_col = E_treatment_wise_intensities[, grep(paste("^",other_treatment_str,"$", sep=""), colnames(E_treatment_wise_intensities))]
  vec_of_lfc = other_treatment_col - cntrl_col
  return(vec_of_lfc)
}
get_lm_obj_stats_two_factors <-
function (lm.obj, dpath){
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  p_vals_we_care_about = p_vals[,grep(paste("^","Tx.*","$",sep=""),colnames(p_vals))]
  class(p_vals_we_care_about) = "matrix"
  q_vals = list()
  for (i in 1:ncol(p_vals_we_care_about)){
    if(sum(is.na(p_vals_we_care_about[,i]))>0){
      rows_with_NAs = which(is.na(p_vals_we_care_about[,i]))
      print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
      #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
      #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
      #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
    }
    plot_p_val_hist_two_factor(p_vals_we_care_about[,i], dpath, colnames(p_vals_we_care_about)[i])
    q_vals[[i]] = qvalue(p_vals_we_care_about[,i])
    names(q_vals)[i] = colnames(p_vals_we_care_about)[i]
    #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
    print(plot(q_vals[[i]], rng=c(0,0.8)))
    #dev.off()
  }
  return(list(betas, p_vals, q_vals))
}
get_lm_obj_stats <-
function (lm.obj, cntrl_str, other_treatment_str, dpath, filename_str){
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  p_vals_we_care_about = as.vector(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))])
  if(sum(is.na(p_vals_we_care_about))>0){
    rows_with_NAs = which(is.na(p_vals_we_care_about))
    print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
    #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
    #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
    #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
  }
  plot_p_val_hist(p_vals_we_care_about, cntrl_str, other_treatment_str, dpath, filename_str)
  q_vals = qvalue(p_vals_we_care_about)
  #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
  print(plot(q_vals, rng=c(0,0.8)))
  #dev.off()
  return(list(betas, p_vals, q_vals))
}
subset_and_lm_by_pair <-
function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("aa",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx))
}
subset_and_lm_by_pair_with_MDS_covariates <-
function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix, MDS_x, MDS_y){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("aa",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx+MDS_x+MDS_y))
}
subset_and_lm_by_pair_with_second_treatment_covariate <-
function (cntrl_treatment_str_tr1, cntrl_treatment_str_tr2, tx, exprs_matrix){
  Tx_1 = tx$Treatment_text
  Tx_1[Tx_1==cntrl_treatment_str_tr1] = paste("AA",cntrl_treatment_str_tr1, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx_1 = as.factor(Tx_1)
  Tx_2 = tx$Treatment_text_2
  Tx_2[Tx_2==cntrl_treatment_str_tr2] = paste("AA",cntrl_treatment_str_tr2, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx_2 = as.factor(Tx_2)
  #print(Tx)
  return(lm(t(exprs_matrix)~Tx_1+Tx_2))
}
subset_and_lm_by_pair_with_percentP_covariate <-
function (cntrl_treatment_str, other_treatment_str, tx, exprs_matrix, percentP_vec){
  Tx = tx$Treatment_text; Tx[Tx==cntrl_treatment_str] = paste("AA",cntrl_treatment_str, sep="")#lm() default contrast setting is dummy contrasts. The default dummy is the first factor in alphabetical order. This is a hack to make sure that's the one we want
  Tx = as.factor(Tx)
  print(Tx)
  return(lm(t(exprs_matrix)~Tx+percentP_vec))
}
generate_lr_p_beta_q_lfc_aveIntensities <-
function(cntrl, other_treatment_str, E_treatment_wise_intensities, E, tx, lm.obj, dpath, filename_str){
  all_stats_list = get_lm_obj_stats (lm.obj, cntrl, other_treatment_str, dpath, filename_str)
  ##begin output matrix with pairwise log ratios/log fold change
  output_matrix = get_pairwise_log_ratio(cntrl, other_treatment_str, E_treatment_wise_intensities) 
  output_matrix = as.matrix(output_matrix)
  ##add betas
  output_matrix = cbind(output_matrix, all_stats_list[[1]][,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(all_stats_list[[1]]))])
  ##add p-vals and q-vals
  output_matrix= cbind(output_matrix, all_stats_list[[2]][,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(all_stats_list[[2]]))], all_stats_list[[3]]$qvalues)
  ##add treatmentwise average intensities
  output_matrix= cbind(output_matrix, E_treatment_wise_intensities[,grep(paste("^",other_treatment_str,"$",sep=""),colnames(E_treatment_wise_intensities))], E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""),colnames(E_treatment_wise_intensities))])
  ##add column names
  colnames(output_matrix)=c(paste(other_treatment_str,"lr",cntrl,sep="_"), paste("beta", other_treatment_str,"vs",cntrl, sep="_"), paste("p_val", other_treatment_str,"vs",cntrl, sep="_"), paste("q_val", other_treatment_str,"vs",cntrl, sep="_"), paste("avg_intensity",other_treatment_str,sep="_"), paste("avg_intensity", cntrl, sep="_"))
  rownames(output_matrix) = rownames(E)
  return(output_matrix)
}
Noramlization_wrap <-
function(wd_path_str=params$proj_path, array_type_str=params$array_type, custom_probeset_mapping_needed = F, probe_ID_vec=NULL, custom_ID_vec=NULL){
  setwd(wd_path_str)
  eCELs = list.celfiles('.',full.names=T)
  switch(array_type_str,
         "PrimeView_Human" = {
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Loess##
           normalized_loess_Data_obj_precursor=normalize.loess(exprs(Data_obj))
           colnames(normalized_loess_Data_obj_precursor) = colnames(Expression_rma_normalized_Data_obj)
           loess_Data_obj = Data_obj
           loess_Data_obj = `exprs<-`(loess_Data_obj, normalized_loess_Data_obj_precursor)
           sampleNames(phenoData(loess_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(loess_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           normalized_loess_Data_obj = computeExprSet(loess_Data_obj, pmcorrect.method="pmonly", summary.method="medianpolish")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj - E_av_all_loess
           
           ##Combine##
           Expression_lr_to_av_all_list = list(E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("RMA_quantile", "LOESS")
           Expression_all_list = list(Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("RMA_quantile", "LOESS")
           
           pckg_name="primeView_human"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj, Expression_non_normalized_affy_obj))
         },
         "HG-U133_Plus_2" ={
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##MAS5##
           mas5_normalized_Data_obj=mas5(Data_obj, normalize=T,sc=500, analysis="absolute")
           Expression_mas5_normalized_Data_obj = log2(exprs(mas5_normalized_Data_obj)) 
           E_av_all_mas5 = rowMeans(Expression_mas5_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_mas5 = Expression_mas5_normalized_Data_obj - E_av_all_mas5
           qc_Data_obj = qc(Data_obj)
           percent_p = percent.present(qc_Data_obj)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Loess##
           normalized_loess_Data_obj=expresso(Data_obj,bgcorrect.method="none", normalize.method="loess", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj- E_av_all_loess
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Combine##
           Expression_lr_to_av_all_list = list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           Expression_all_list = list(Expression_mas5_normalized_Data_obj, Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           
           library(hgu133plus2.db)
           pckg_name="hgu133plus2"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj,Expression_non_normalized_affy_obj))
         },
         "Mouse430_2" ={
           Data_obj=ReadAffy(celfile.path=wd_path_str)
           
           ##MAS5##
           mas5_normalized_Data_obj=mas5(Data_obj, normalize=T,sc=500, analysis="absolute")
           Expression_mas5_normalized_Data_obj = log2(exprs(mas5_normalized_Data_obj))
           E_av_all_mas5 = rowMeans(Expression_mas5_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_mas5 = Expression_mas5_normalized_Data_obj - E_av_all_mas5
           qc_Data_obj = qc(Data_obj)
           percent_p = percent.present(qc_Data_obj)
           
           ##RMA##
           normalized_rma_Data_obj=expresso(Data_obj, bgcorrect.method="rma", normalize.method="quantiles", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_rma_normalized_Data_obj = exprs(normalized_rma_Data_obj)
           E_av_all_rma = rowMeans(Expression_rma_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_rma = Expression_rma_normalized_Data_obj - E_av_all_rma
           
           ##Loess##
           normalized_loess_Data_obj=expresso(Data_obj,bgcorrect.method="none", normalize.method="loess", summary.method="medianpolish", pmcorrect.method="pmonly")
           Expression_loess_normalized_Data_obj = exprs(normalized_loess_Data_obj)
           E_av_all_loess = rowMeans(Expression_loess_normalized_Data_obj, na.rm=TRUE)
           E_lr_av_all_loess = Expression_loess_normalized_Data_obj- E_av_all_loess
           
           ##Unnormalized##
           non_norm_Data_obj = Data_obj
           sampleNames(phenoData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           sampleNames(protocolData(non_norm_Data_obj)) = colnames(Expression_rma_normalized_Data_obj)
           non_normalized_summarized = expresso(Data_obj, bg.correct = F, normalize = F, pmcorrect.method = "pmonly", summary.method = "medianpolish")
           Expression_non_normalized_affy_obj = exprs(non_normalized_summarized) #maybe needs to run through bg subtract and summarization?
           
           ##Combine##
           Expression_all_list = list(Expression_mas5_normalized_Data_obj, Expression_rma_normalized_Data_obj, Expression_loess_normalized_Data_obj)
           names(Expression_all_list) = c("MAS5", "RMA_quantile", "LOESS")
           Expression_lr_to_av_all_list = list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           names(Expression_lr_to_av_all_list) = c("MAS 5", "RMA quantile", "LOESS")
           
           library(mouse4302.db)
           pckg_name="mouse4302"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj,Expression_non_normalized_affy_obj))
         },
         "HTA-2_0" ={ 
           Data_obj = read.celfiles(eCELs)
           probe_mats_ls_all = list()
           
           ##MAS5##
           #No mismatch probes for ST arrays, so MAS5 not a thing
           
           ##Unnormalized##
           #summarize with the multimappers
           non_norm_Data_obj_precursor = exprs(Data_obj)
           Complex_list = complex_process_probes(probe_mat=non_norm_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_non_norm = Complex_list[[1]]
           probe_mats_ls_non_norm = Complex_list[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all, probe_mats_ls_non_norm)
           summarized_crazies_removed_non_norm = Complex_list[[3]]
           merged_no_multi_summarized_non_norm = Complex_list[[4]]
           E_lr_av_all_ls_non_norm = Complex_list[[5]]
           
           non_norm_ls = list(summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm)
           names(non_norm_ls) = c("All_multi_non_norm", "Min_multi_non_norm")
           #probe_level_abundance_list = probe_mats_ls_non_norm
           #Note that some of these objects will have a non-numeric column
           
           ##RMA##
           normalized_rma_Data_obj_precursor=preprocessCore:::normalize.quantiles(exprs(Data_obj))
           rownames(normalized_rma_Data_obj_precursor) = rownames(exprs(Data_obj))
           colnames(normalized_rma_Data_obj_precursor) = colnames(exprs(Data_obj))
           Complex_list_rma = complex_process_probes(probe_mat=normalized_rma_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_rma = Complex_list_rma[[1]]
           probe_mats_ls_rma = Complex_list_rma[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all,probe_mats_ls_rma)
           summarized_crazies_removed_rma = Complex_list_rma[[3]]
           merged_no_multi_summarized_rma = Complex_list_rma[[4]]
           E_lr_av_all_ls_rma = Complex_list_rma[[5]]
           
           ##Loess##
           #summarize with the multimappers but without offenders
           normalized_loess_Data_obj_precursor=normalize.loess(exprs(Data_obj))
           normalized_loess_Data_obj_precursor_log_first=normalize.loess(log2(exprs(Data_obj)))
           Complex_list_loess = complex_process_probes(probe_mat=normalized_loess_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_loess = Complex_list_loess[[1]]
           probe_mats_ls_loess = Complex_list_loess[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all, probe_mats_ls_loess)
           probe_mats_ls_loess[[length(probe_mats_ls_loess)+1]] = normalized_loess_Data_obj_precursor_log_first
           names(probe_mats_ls_loess)[length(probe_mats_ls_loess)] = c("Log_first_Loess")
           summarized_crazies_removed_loess = Complex_list_loess[[3]]
           merged_no_multi_summarized_loess = Complex_list_loess[[4]]
           E_lr_av_all_ls_loess = Complex_list_loess[[5]]
           
           ##Combine##
           Expression_all_list = list(summarized_crazies_removed_rma, merged_no_multi_summarized_rma, summarized_crazies_removed_loess, merged_no_multi_summarized_loess)
           names(Expression_all_list) = c("RMA_all_multi", "RMA_min._multi","Loess_all_multi", "Loess_min._multi")
           Parallel_non_norm_ls = list(summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm,summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm)
           names(Parallel_non_norm_ls) = c("Non_norm_all_multi", "Non_norm_min_multi","Non_norm_all_multi", "Non_norm_min_multi")
           #Expression_all_list = match_up_row_order(Expression_minimally_multi_mapping_no_norm_Data_obj, Expression_all_list, is_second_arg_a_list=T)
           #Expression_lr_to_av_all_list = list(E_lr_av_all_rma, E_lr_av_all_loess) #E_lr_av_all_no_norm,  #list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           Expression_lr_to_av_all_list = append(E_lr_av_all_ls_rma, E_lr_av_all_ls_loess)
           names(Expression_lr_to_av_all_list) = c("RMA_min._multi","RMA_all_multi","Loess_min._multi","Loess_all_multi")
           
           #library(hta20sttranscriptcluster.db)
           pckg_name="hta20sttranscriptcluster"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj, non_norm_ls, probe_mats_ls_all,Parallel_non_norm_ls))
         },
         "MTA-1_0" ={ 
           Data_obj = read.celfiles(eCELs)
           probe_mats_ls_all = list()
           
           ##MAS5##
           #No mismatch probes for ST arrays, so MAS5 not a thing
           
           ##Unnormalized##
           #summarize with the multimappers
           non_norm_Data_obj_precursor = exprs(Data_obj)
           Complex_list = complex_process_probes(probe_mat=non_norm_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_non_norm = Complex_list[[1]]
           probe_mats_ls_non_norm = Complex_list[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all, probe_mats_ls_non_norm)
           summarized_crazies_removed_non_norm = Complex_list[[3]]
           merged_no_multi_summarized_non_norm = Complex_list[[4]]
           E_lr_av_all_ls_non_norm = Complex_list[[5]]
           
           non_norm_ls = list(summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm)
           names(non_norm_ls) = c("All_multi_non_norm", "Min_multi_non_norm")
           #probe_level_abundance_list = probe_mats_ls_non_norm
           #Note that some of these objects will have a non-numeric column
           
           ##RMA##
           normalized_rma_Data_obj_precursor=preprocessCore:::normalize.quantiles(exprs(Data_obj))
           rownames(normalized_rma_Data_obj_precursor) = rownames(exprs(Data_obj))
           colnames(normalized_rma_Data_obj_precursor) = colnames(exprs(Data_obj))
           Complex_list_rma = complex_process_probes(probe_mat=normalized_rma_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_rma = Complex_list_rma[[1]]
           probe_mats_ls_rma = Complex_list_rma[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all,probe_mats_ls_rma)
           summarized_crazies_removed_rma = Complex_list_rma[[3]]
           merged_no_multi_summarized_rma = Complex_list_rma[[4]]
           E_lr_av_all_ls_rma = Complex_list_rma[[5]]
           
           ##Loess##
           #summarize with the multimappers but without offenders
           normalized_loess_Data_obj_precursor=normalize.loess(exprs(Data_obj))
           normalized_loess_Data_obj_precursor_log_first=normalize.loess(log2(exprs(Data_obj)))
           Complex_list_loess = complex_process_probes(probe_mat=normalized_loess_Data_obj_precursor,probe_ID_vec=probe_ID_vec, custom_ID_vec=custom_ID_vec, quantile_to_remove=0.99999)
           summarized_including_all_loess = Complex_list_loess[[1]]
           probe_mats_ls_loess = Complex_list_loess[[2]]
           probe_mats_ls_all = append(probe_mats_ls_all, probe_mats_ls_loess)
           probe_mats_ls_loess[[length(probe_mats_ls_loess)+1]] = normalized_loess_Data_obj_precursor_log_first
           names(probe_mats_ls_loess)[length(probe_mats_ls_loess)] = c("Log_first_Loess")
           summarized_crazies_removed_loess = Complex_list_loess[[3]]
           merged_no_multi_summarized_loess = Complex_list_loess[[4]]
           E_lr_av_all_ls_loess = Complex_list_loess[[5]]
           
           ##Combine##
           Expression_all_list = list(summarized_crazies_removed_rma, merged_no_multi_summarized_rma, summarized_crazies_removed_loess, merged_no_multi_summarized_loess)
           names(Expression_all_list) = c("RMA_all_multi", "RMA_min._multi","Loess_all_multi", "Loess_min._multi")
           Parallel_non_norm_ls = list(summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm,summarized_crazies_removed_non_norm, merged_no_multi_summarized_non_norm)
           names(Parallel_non_norm_ls) = c("Non_norm_all_multi", "Non_norm_min_multi","Non_norm_all_multi", "Non_norm_min_multi")
           #Expression_all_list = match_up_row_order(Expression_minimally_multi_mapping_no_norm_Data_obj, Expression_all_list, is_second_arg_a_list=T)
           #Expression_lr_to_av_all_list = list(E_lr_av_all_rma, E_lr_av_all_loess) #E_lr_av_all_no_norm,  #list(E_lr_av_all_mas5, E_lr_av_all_rma, E_lr_av_all_loess)
           Expression_lr_to_av_all_list = append(E_lr_av_all_ls_rma, E_lr_av_all_ls_loess)
           names(Expression_lr_to_av_all_list) = c("RMA_min._multi","RMA_all_multi","Loess_min._multi","Loess_all_multi")
           
           #library(hta20sttranscriptcluster.db)
           #pckg_name="hta20sttranscriptcluster"
           return(list(Expression_all_list, Expression_lr_to_av_all_list, Data_obj, non_norm_ls, probe_mats_ls_all,Parallel_non_norm_ls))
         },
         stop("Enter something that switches me!")
  )
}
Generate_mean_by_treatment_matrix <-
function(Normalized_abundance_list, vec_of_treatments){
  #Function to average abundance values across treatment and return a list of matrices with same rownames and only one col. per unique treatment type.
  #Normalized_abundance_list is a list of abundance matrices, data frames, or data tables (with names() of the list populated). Each matrix/data frame/data table has rownames() and colnames()
  #vec_of_treatments is a vector with a treatment assigned to each sample/column. E.g., if .CEL files were called 001.CEL, 002.CEL, 003.CEL, and 004.CEL, and 001 and 002 are treatment x and 003 and 004 are treatment y, vector would be c("x", "x", "y", "y").
  #Can deal with treatments that only have one column of data/one sample
  #Returns a list of matrices in the same order as the original list. Each colname with be a treatment name from unique(vec_of_treatments) concatenated to, "_Mn"
  final_mat_list =list()
  final_mat = NULL
  for (i in 1:length(Normalized_abundance_list)){
    for (a in 1:length(unique(vec_of_treatments))){
      cols_of_interest = which(vec_of_treatments==unique(vec_of_treatments)[a])
      if (length(cols_of_interest)==1){
        final_mat = cbind(final_mat, as.data.frame(Normalized_abundance_list[[i]])[,cols_of_interest])
      } else{
        final_mat = cbind(final_mat, rowMeans(as.data.frame(Normalized_abundance_list[[i]])[,cols_of_interest]))
      }
      colnames(final_mat)[a] =paste0(unique(vec_of_treatments)[a],"_Mn")
    }
    rownames(final_mat) = rownames(Normalized_abundance_list[[i]])
    final_mat_list[[i]] = final_mat
    names(final_mat_list)[i] = names(Normalized_abundance_list)[i]
    rm(final_mat)
    final_mat=NULL
  }
  return(final_mat_list)
}
generate_lr_p_q_aveIntensities_wilcox <-
function(cntrl, other_treatment_str, E_treatment_wise_intensities, E, tx, all_stats_list){
  ##begin output matrix with pairwise log ratios/log fold change
  output_matrix = get_pairwise_log_ratio(cntrl, other_treatment_str, E_treatment_wise_intensities) 
  output_matrix = as.matrix(output_matrix)
  
  ##add p-vals and q-vals
  output_matrix= cbind(output_matrix, all_stats_list[[1]], all_stats_list[[2]]$qvalues)
  
  ##add treatmentwise average intensities
  output_matrix= cbind(output_matrix, E_treatment_wise_intensities[,grep(paste("^",other_treatment_str,"$",sep=""),colnames(E_treatment_wise_intensities))], E_treatment_wise_intensities[,grep(paste("^",cntrl,"$",sep=""),colnames(E_treatment_wise_intensities))])
  
  ##add column names
  colnames(output_matrix)=c(paste(other_treatment_str,"lr",cntrl,sep="_"), paste("p_val", other_treatment_str,"vs",cntrl, sep="_"), paste("q_val", other_treatment_str,"vs",cntrl, sep="_"), paste("avg_intensity",other_treatment_str,sep="_"), paste("avg_intensity", cntrl, sep="_"))
  rownames(output_matrix) = rownames(E)
  
  return(output_matrix)
}
plot_p_val_hist_two_factor <-
function (pval_vec, dpath, filename_str){
  print(hist(pval_vec, breaks = round(length(pval_vec)/150), main = filename_str, xlab="P-value"))
}
get_wilcox_based_qval_vec <-
function (wilcox_pval_vec, cntrl_str, other_treatment_str, dpath){
  q_vals = qvalue(wilcox_pval_vec);
  #png(filename=paste(dpath,"wilcox_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144);
  print(plot(q_vals));
  #dev.off();
  return(q_vals);
}
plot_p_val_hist <-
function (pval_vec, cntrl_str, other_treatment_str, dpath, filename_str){
  #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_pDist.png", sep=""),width=5,height=5.4,units="in",res=144)
  print(hist(pval_vec, breaks = round(length(pval_vec)/150), main = paste(filename_str, other_treatment_str, "vs.", cntrl_str, sep=" "), xlab="P-value"))
  #dev.off()
}
wilcox_loop <-
function (cntrl_str, other_treatment_str, tx, exprs_matrix){
  wilcox_obj_ls = NULL
  for (x in 1:nrow(exprs_matrix)){
    w_obj = wilcox.test(exprs_matrix[x,grep(cntrl_str,colnames(exprs_matrix))],exprs_matrix[x,grep(other_treatment_str,colnames(exprs_matrix))])$p.value
    wilcox_obj_ls = c(wilcox_obj_ls, w_obj)
    names(wilcox_obj_ls)[x] = paste(other_treatment_str, "vs", cntrl_str, rownames(exprs_matrix)[x], sep="_")
  }
  return (wilcox_obj_ls)
}
lm_factor_any_number_of_covariates_dummy_contrasts <-
function(vec_of_covariate_strings, vec_of_control_strings, phenotype_dt, exprs_matrix){
  #Return an object of type lm() using dummy contrasts
  #Could be made less hacky using factor(tx$Treatment_text, levels=c("Onex","Twox")) instead of adding AA in front of control
  #A generalization of subset_and_lm_by_pair_with_second_treatment_covariate and subset_and_lm_by_pair_with_MDS_covariates with any number of covariates
  #vec_of_covariate_strings is a vector of strings corresponding to colnames(phenotype_dt) of interest to include as covariates
  #vec_of_control_strings is a vector of strings, each element corresponding to the control treatment in the column referred to in the corresponding element in vec_of_covariate_strings
  #phenotype_dt is a data table with the each row corresponding to phenotype of of a sample. Name of phenotype is the colname() of that column. Example below
  # Inv_sample_name Treatment_text Treatment_text_2 Treatment_text_3 Treatment_text_4 Treatment_text_5
  #1                9           Twox              Low          Jan2716            X0077               OS
  #2               10           Onex             High          Oct1215            X1308               OD
  #3               11           Twox              Low           Feb116            X0108               OD
  #4               12           Twox             High          Oct1215            X1308               OS
  #exprs_matrix is any matrix of abundance data with rownames() and colnames()
  
  covariate_list =list()
  for (i in 1:length(vec_of_covariate_strings)){
    col_of_interest = which(colnames(phenotype_dt)==vec_of_covariate_strings[i])
    #assign(paste0("TX_",i), phenotype_dt$vec_of_covariate_strings[i])
    covariate_list[[i]] = phenotype_dt[[col_of_interest]]
    covariate_list[[i]] [covariate_list[[i]]==vec_of_control_strings[i]] = paste0("AA", vec_of_control_strings[i])
    covariate_list[[i]] = as.factor(unlist(covariate_list[[i]], use.names = F))
    names(covariate_list)[i] = paste0("TX_",i)
  }
  master_str = "lm(t(exprs_matrix)~"
  for (i in 1:length(covariate_list)){
    new_str = paste0("unlist(covariate_list[[",i,"]])")
    master_str = paste0(master_str, new_str)
    if (i<length(covariate_list)){ #add a plus sign to everything but the last one
      master_str = paste0(master_str, "+")
    } 
  }
  master_str = paste0(master_str,")")
  return(eval(parse(text=master_str)))
}
get_lm_obj_stats_any_num_covariates_dummy <-
function (lm.obj, vec_of_other_treatment_strings, vec_of_control_strings, dpath, filename_str){
  #WARNING: can't remember whether this still needs troubleshooting
  #Returns a list with element 1 being a matrix of regression coefficients (betas)
  #Element 2 is a matrix of p-values
  #Element 3 is a list produced by the qvalue() function
  #lm.obj is an object created by an lm() call to a matrix. The colnames() of that matrix were temporarily assigned strings where treatment corresponded to the particular sample (e.g., 1, 2, 3, 4 were reassigned High Low High Low) before being passed to one of the lm
  #vec_of_other_treatment_strings is a vector of characters ...
  betas = t(lm.obj$coefficients)
  p_vals = t(lm.pval(lm.obj)$pval)
  q_val_list = list()
  for(other_treatment_string_cnt in 1: length(vec_of_other_treatment_strings)){
    p_vals_we_care_about = as.vector(p_vals[,grep(paste("^","unlist",vec_of_other_treatment_strings[other_treatment_string_cnt],"$",sep=""),colnames(p_vals))])
    #Report when there was an NA p-value, which has potential to gum up downstream works
    #Need to standardize what to do here
    if(sum(is.na(p_vals_we_care_about))>0){
      rows_with_NAs = which(is.na(p_vals_we_care_about))
      print(paste0("NA p-values in rows ", rows_with_NAs, "converting to mean p-val and mean intercept across the array"))
      #p_vals_we_care_about[rows_with_NAs] = mean(p_vals_we_care_about, na.rm=T) #this may or may not mess with q-value's pi0 assignment; investigate the p-val distribution of q-value QC plots
      #p_vals[rows_with_NAs,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("^","Tx",other_treatment_str,"$",sep=""),colnames(p_vals))],na.rm=T)
      #p_vals[rows_with_NAs,grep(paste("(Intercept)",sep=""),colnames(p_vals))] = mean(p_vals[,grep(paste("(Intercept)",sep=""),colnames(p_vals))],na.rm=T)
    }
    plot_p_val_hist(p_vals_we_care_about, cntrl_str, other_treatment_str, dpath, filename_str)
    q_val_list[[other_treatment_string_cnt]] = qvalue(p_vals_we_care_about)
    names(q_val_list)[other_treatment_string_cnt] = paste0("qvalue_of_", vec_of_other_treatment_strings[other_treatment_string_cnt])
    #png(filename=paste(dpath, filename_str,"_", other_treatment_str, "_vs_", cntrl_str,"_q_val_data.png", sep=""),width=5,height=5.4,units="in",res=144)
    print(plot(q_val_list[[other_treatment_string_cnt]], rng=c(0,0.8)))
    #dev.off()
  }
  return(list(betas, p_vals, q_val_list))
}
generate_adj_log_ratio_p_val_q_val_mn_abundance_general <-
function(lm.obj, dpath, abundance_mat_w_rownames,name_of_abundance_mat_w_rownames=NULL, grep_str_for_cols_we_want = '^.*$', vec_of_strings_you_want_to_rename_by=NULL, feature_ID = "Probe_Set_ID",annotation_status=F,annotation_path = getwd(), feature_ID_col_name_in_annotation="transcript_cluster_id", phenotype_dt=tx, vec_of_phenotype_colnames_you_want_to_take_treatment_means_of = colnames(tx)[-1]){
  betas = t(lm.obj$coefficients)
  rownames(betas) = rownames(abundance_mat_w_rownames)
  beta_cols_we_care_about = grep(grep_str_for_cols_we_want, colnames(betas))
  if(is.null(vec_of_strings_you_want_to_rename_by)){
    vec_of_strings_you_want_to_rename_by = colnames(betas)[beta_cols_we_care_about]
  }
  if(length(beta_cols_we_care_about)==0){
    print("No columns matched your search string for regression coefficients")
    return(NULL)
  }
  if(length(beta_cols_we_care_about)==1){
    betas_we_want = as.matrix(betas[,beta_cols_we_care_about])
    colnames(betas_we_want) = gsub('^','Adj_lg_ratio_',vec_of_strings_you_want_to_rename_by)
  }else{
    betas_we_want = betas[,beta_cols_we_care_about]
    colnames(betas_we_want) = gsub('^','Adj_lg_ratio_',vec_of_strings_you_want_to_rename_by)
  }
  p_vals = t(lm.pval(lm.obj)$pval)
  class(p_vals) = "matrix"
  rownames(p_vals) = rownames(abundance_mat_w_rownames)
  p_val_cols_we_want = grep(grep_str_for_cols_we_want,colnames(p_vals))
  if(length(p_val_cols_we_want)==0){
    print("No columns matched your search string for p-values")
    return(NULL)
  }
  if(length(p_val_cols_we_want)==1){
    p_vals_we_want = as.matrix(p_vals[,p_val_cols_we_want])
    colnames(p_vals_we_want) = gsub('^','P_val_',vec_of_strings_you_want_to_rename_by)
  }else{
    p_vals_we_want = p_vals[,p_val_cols_we_want]
    colnames(p_vals_we_want) = gsub('^','P_val_',vec_of_strings_you_want_to_rename_by)
  }
  q_vals = list()
  for (i in 1:ncol(p_vals_we_want)){
    if(sum(is.na(p_vals_we_want[,i]))>0){
      rows_with_NAs = which(is.na(p_vals_we_want[,i]))
      print(paste0("NA p-values in rows ", rows_with_NAs)) #, "converting to mean p-val and mean intercept across the array"))
      browser()
    }
    plot_p_val_hist_two_factor(p_vals_we_want[,i], dpath, colnames(p_vals_we_want)[i])
    q_vals[[i]] = qvalue(p_vals_we_want[,i])
    names(q_vals)[i] = gsub('^P_val_','Q_val_', colnames(p_vals_we_want)[i])
  }
  all_stats_list = list(betas_we_want, p_vals_we_want, q_vals)
  names(all_stats_list) = c("betas", "p_vals", "q_vals")
  
  ##begin output matrix with transcript cluster IDs
  output_matrix = rownames(all_stats_list[[1]]) 
  output_matrix = as.matrix(output_matrix)
  colnames(output_matrix) = c(feature_ID)
  
  ##optionally add gene symbols
  if(annotation_status==T){
    annotation_dt = fread(annotation_path, header=T)
    #idx1 = match(annotation_dt[["transcript_cluster_id"]], output_matrix[,1])
    idx1 = match(output_matrix[,1],annotation_dt[["transcript_cluster_id"]])
    output_matrix = cbind(output_matrix, annotation_dt[["Gene_symbol"]][idx1])
    colnames(output_matrix)[ncol(output_matrix)] = "Gene_Symbol"
  }
  
  ##add betas
  output_matrix = cbind(output_matrix, all_stats_list[[1]])
  
  ##add p-vals
  output_matrix= cbind(output_matrix, all_stats_list[[2]])
  
  ##add q-values
  for (i in 1:length(all_stats_list[[3]])){
    output_matrix= cbind(output_matrix,all_stats_list[[3]][[i]]$qvalues)
    colnames(output_matrix)[ncol(output_matrix)] = names(all_stats_list[[3]])[i]
  }
  
  ##add treatment-wise mean expression
  Abundance_final_ls = list(abundance_mat_w_rownames)
  names(Abundance_final_ls)[1] = name_of_abundance_mat_w_rownames
  for(p in 1:length(vec_of_phenotype_colnames_you_want_to_take_treatment_means_of)){
    Tx_means = Generate_mean_by_treatment_matrix(Normalized_abundance_list=Abundance_final_ls, vec_of_treatments=phenotype_dt[[vec_of_phenotype_colnames_you_want_to_take_treatment_means_of[p]]])[[1]]
    output_matrix= cbind(output_matrix,Tx_means)
  }
  
  final_rtrn = list(all_stats_list,output_matrix)
  names(final_rtrn) = c("list_of_beta_p_and_q", "master_matrix_with_betas_ps_qs_treatmnt_mn")
  return(final_rtrn)
}
rqd.libraries <-
function (path) {
# path indicates the location of all other scripts used for analysis

  # Required libraries
  require(data.table)
  require(oligo)
  require(affy)
  require(limma)
  require(gdata)
  require(RColorBrewer)
  require(qvalue)
  require(limma)
  require(AnnotationDbi)
  library(NMF) # use aheatmap() function
  library(hexbin)
  library(grid)
  library(stats)
#  source("https://bioconductor.org/biocLite.R")
  biocLite("marray")

  # Custom files all located in directory indicated in 'path'
#  RFiles = "/Users/lusardi/Documents/RFiles"
  source(paste(path, '/lm.pval.R', sep = ''))
  source(paste(path, '/gs.pval.R', sep = ''))
  source(paste(path, '/gs.wrapper.R', sep = ''))
  source(paste(path, '/array_functions_mf.R', sep = ''))
  
  cat("Made it to the libraries call!\n")

  # Basic color map settings for line plots.
  # colormap for plotting
  colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7)) 
  colmap = colorRampPalette( c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7)) )

}
vec_of_gene_symbols_given_vec_of_IDs <-
function(vec_of_IDs, master_mat, gene_symbol_colname_in_master, ID_colname_in_master){
  #returns a vector of strings that map to the input vector of strings in a given matrix
  #vec_of_IDs is a vector of strings
  #master_mat is the matrix containing at least two columns, one of which is a superset of vec_of_IDs and the other of which is a superset of the matched strings
  #gene_symbol_colname_in_master is the column name in master_mat corresponding to the target strings that you want as output
  #ID_colname_in_master is the column name in master_mat corresponding to the input vector of strings
  master_mat_match = match(vec_of_IDs,master_mat[,ID_colname_in_master])
  idx1=which(!is.na(master_mat_match))
  master_mat_match=master_mat_match[idx1]
  return(master_mat[master_mat_match, gene_symbol_colname_in_master])
}
norm_matrix <-
function(...){
  # Description of function:
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

  warning("Function norm_matrix() has been deprecated and moved to deprecated.R. Please use normMatrix() in processData.R instead.")
}
merge_affy_provided_annotation_with_big_stats_mat <-
function (big_stats_mat, annotation_csv_path){
  annotation_dt = fread(annotation_csv_path)
  big_stats_dt = as.data.table(big_stats_mat) #retains colnames, loses rownames
  setkey(annotation_dt, "Probe Set ID")
  setkey(big_stats_dt, Probeset_IDs)
  final_dt = annotation_dt[big_stats_dt]
  return(final_dt)
}
merge_annotation_dt_with_big_stats_mat <-
function(big_stats_mat, big_stats_mat_key, annotation_dt, annotation_dt_key){
  ##assumes that the annotation_dt has two columns: probeset and gene symbol and that big_stat_mat has a probeset ID column. The string representing the probeset colname in each case gets passed to the function
  big_stats_dt = as.data.table(big_stats_mat)
  setkeyv(big_stats_dt, c(big_stats_mat_key))
  setkeyv(annotation_dt, c(annotation_dt_key))
  return(annotation_dt[big_stats_dt])
}
merge_affy_pck_annotation_with_big_stats_mat <-
function (big_stats_mat, pckg_name){
  colnum_of_interest = which(colnames(big_stats_mat)=="Probeset_IDs")
  pckg_probesetIDs_df = as.data.frame(big_stats_mat[,colnum_of_interest])
  colnames(pckg_probesetIDs_df)=c(paste(pckg_name,"probesetIDs", sep="_"))
  
  x <- get(paste(pckg_name,"ENTREZID", sep="")) # Get the probe identifiers that are mapped to an ENTREZ Gene ID 
  mapped_probes <- mappedkeys(x) # Convert to a list 
  pckg_probesetID_to_pckg_entrez_ls <- as.list(x[mapped_probes])
  pckg_probesetIDs = names(pckg_probesetID_to_pckg_entrez_ls)
  pckg_entrez_IDs = unlist(pckg_probesetID_to_pckg_entrez_ls, use.names=F)
  pckg_probesetID_to_pckg_entrez_df=as.data.frame(cbind(pckg_probesetIDs, pckg_entrez_IDs))
  colnames(pckg_probesetID_to_pckg_entrez_df) = c(paste(pckg_name,"probesetIDs", sep="_"), paste(pckg_name,"entrez_IDs", sep="_"))
  
  pckg_merged = merge(pckg_probesetIDs_df, pckg_probesetID_to_pckg_entrez_df, all.x=TRUE, by.x = paste(pckg_name,"probesetIDs",sep="_"), by.y = paste(pckg_name,"probesetIDs",sep="_"))
  
  x <- get(paste(pckg_name,"SYMBOL", sep="")) # Get the probe identifiers that are mapped to an ENTREZ Gene ID 
  mapped_probes <- mappedkeys(x) # Convert to a list 
  pckg_probesetID_to_pckg_symbol_ls <- as.list(x[mapped_probes])
  pckg_probesetIDs = names(pckg_probesetID_to_pckg_symbol_ls)
  pckg_symbol_IDs = unlist(pckg_probesetID_to_pckg_symbol_ls, use.names=F)
  pckg_probesetID_to_pckg_symbol_df=as.data.frame(cbind(pckg_probesetIDs, pckg_symbol_IDs))
  colnames(pckg_probesetID_to_pckg_symbol_df) = c(paste(pckg_name,"probesetIDs", sep="_"), paste(pckg_name, "symbol_IDs", sep="_"))
  
  pckg_merged = merge(pckg_merged, pckg_probesetID_to_pckg_symbol_df, all.x=TRUE, by.x =paste(pckg_name,"probesetIDs",sep="_"), by.y = paste(pckg_name,"probesetIDs",sep="_"))
  big_stats_mat = rename_nonunique_cols(big_stats_mat)
  pckg_merged = merge(pckg_merged, big_stats_mat, all=T, by.x=paste(pckg_name,"probesetIDs",sep="_"), by.y="Probeset_IDs")
  return(pckg_merged)
}
map_hsa_entrez_to_mogene20st_probesetIDs <-
function(vec_of_hsa_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst){
  rm(vec_of_mogene20st_probesetIDs)
  ##translate hsa entrez IDs to mmu entrez IDs
  ls_of_mmu_entrez_IDs = inpIDMapper(vec_of_hsa_entrez_IDs, srcSpecies="HOMSA", destSpecies="MUSMU",srcIDType="EG",destIDType="EG", keepMultDestIDMatches = FALSE)
  ##this is a list where names() corresponds to human entrez IDs and the items in the list correspond to mouse entrez IDs. One problem will be if multiple mouse entrez IDs map to the same human entrez ID, in which case the length of that element in the list will be >1. I don't think that these are ones we want, but let's check whether we even have them##
  xs=lapply(ls_of_mmu_entrez_IDs, function(e) length(e)>1)
  xs=unlist(xs)
  print(paste("there are ",length(names(ls_of_mmu_entrez_IDs)[xs]), " mouse entrez IDs that have more than one human entrez ID mapped to them", sep=""))
  ##for now, subset the mouse entrez ids by those not included in this list##
  xsub=lapply(ls_of_mmu_entrez_IDs, function(e) e[1])
  xsub=unlist(xsub)
  vec_of_mmu_entrez_IDs=xsub[!xs]
  return (map_mmu_entrez_to_mogene20st_probesetIDs(vec_of_mmu_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst))
}
map_mmu_entrez_to_mogene20st_probesetIDs <-
function(vec_of_mmu_entrez_IDs, mogene20st_trans_clust_ENTREZID_lst){
  if (length(vec_of_mmu_entrez_IDs) != length(unique(vec_of_mmu_entrez_IDs))){
    print ("your mouse entrez id vector contains non-unique items")
  }
  print(vec_of_mmu_entrez_IDs[1])
  vec_of_mogene20st_probesetIDs=c()
  for (mmus_eID in vec_of_mmu_entrez_IDs){
    vec_of_mogene20st_probesetIDs=c(vec_of_mogene20st_probesetIDs, names(which(mogene20st_trans_clust_ENTREZID_lst==mmus_eID)))	
  }
  return (vec_of_mogene20st_probesetIDs)
}
get_mat_of_greater_than_ave_sd_probesets <-
function(E_lr_av_all_data_frame){
  #takes a data frame of abundance to av all ratios with proper rownames and returns a subset of rows where the standard dev. in ratio is greater than global standard dev.
  sd_by_row = sapply(1:nrow(E_lr_av_all_data_frame),function(x){sd(E_lr_av_all_data_frame[x,],na.rm=T)})
  sd_by_row = as.matrix(sd_by_row)
  rownames(sd_by_row)=rownames(E_lr_av_all_data_frame)
  global_sd = sd(as.matrix(E_lr_av_all_data_frame))
  sd_greater_log_mat = sd_by_row>global_sd
  E_lr_av_changing_genes = E_lr_av_all_data_frame[sd_greater_log_mat[,1],]
  return(E_lr_av_changing_genes)
}
