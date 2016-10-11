#####################################################################################################
# Functions useful for microarray analysis.
#
# convertCELtoNonNormMat           --- read CEL files and output matrix of intensities
#
# Authors: Mark Fisher
# Started - October, 2016
#####################################################################################################

#Depends on rqd.libraries from Proteomics.Tools.x1.R
library(oligo)
library(affy)

convertCELtoNonNormMat = function(pathToCELfiles){
  #' Read in CEL files from a directory argument, using the oligo package first and then, failing that, trying the affy package instead. Requires the relevent array package to be installed from BioConductor before run.
  #' 
  #' @param pathToCELfiles a string containing the path to the directory containing 
  #' @return A matrix with rownames and colnames that contains non-normalized intensities for each probe in a given array
  #' @examples
  #' path2CELfiles1 = "/Users/fishema/Desktop/CEL_files_for_testing/HTA/"
  #' test1_mat = convertCELtoNonNormMat(path2CELfiles1)
  #' path2CELfiles2 = "/Users/fishema/Desktop/CEL_files_for_testing/Rhesus/"
  #' test2_mat = convertCELtoNonNormMat(path2CELfiles2)
  #' path2CELfiles3 = "/Users/fishema/Desktop/CEL_files_for_testing/SNP6/"
  #' test3_mat = convertCELtoNonNormMat(path2CELfiles3)
  
  tryCatch({
    eCELs = list.celfiles(pathToCELfiles,full.names=T)
    Data_obj = read.celfiles(eCELs)
    return_val = exprs(Data_obj)
  }, warning = function(war) {
    print(paste("Warning in convertCELtoNonNormMat function:  ",war))
    return_val = NULL
    return(return_val)
  }, error = function(err) {
    # Might be affy then
    print(paste("Array type not suitable for oligoClasses package. Trying affy package:  ",err))
    Data_obj=ReadAffy(celfile.path=pathToCELfiles)
    non_norm_mat = exprs(Data_obj)
    return_val = non_norm_mat
    return(return_val)
  }, finally = {
    print(paste("Done with convertCELtoNonNormMat"))
  }) # END tryCatch
} #END convertCELtoNonNormMat