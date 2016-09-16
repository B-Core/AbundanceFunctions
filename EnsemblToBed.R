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
