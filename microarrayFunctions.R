#####################################################################################################
# Functions useful for microarray analysis.
#
# convertCELtoNonNormMat           --- read CEL files and output matrix of intensities
#
# Authors: Mark Fisher
# Started - October, 2016
#####################################################################################################

##Proteomics.Tools.x1.R for calling proper libraries
repoPath = "~/Git_repos/AbundanceFunctions/"
source(paste0(repoPath,"Proteomics.Tools.x1.R"))

convertCELtoNonNormMat = function(pathToCELfiles){
  #' Read in CEL files from a directory argument, using the oligo package first and then, failing that, trying the affy package instead. Requires the relevent array package to be installed from BioConductor before run. Not sure of a good strategy for installing these if they haven't already been. Also want more sophistocated error handling here, but SO recommends against nested tryCatches(), because they don't work as expected:
  #' http://stackoverflow.com/questions/35716394/nested-try-catch-in-r
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
  
  require(oligo)
  require(affy)
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