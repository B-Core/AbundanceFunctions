#######################################################################################
# 
# Plot wrappers... Allow a tidy line of code designed to loop through multiple simple plots
#
#  qcQ_wrapper    --- Plots p-value and q matrix for all factors of a linear regression 
#                     Optionally plots only selected factors
#
# Authors:  Theresa Lusardi
# Started:  October 2016
#######################################################################################


# ********Q-value QC Wrapper *********************************************************
qcQ_wrapper = function(norm_x, p_mat, q_list, metadata, plotdata, plot2file = TRUE,
                       bonusMDS = FALSE, oneclass = NULL, qcut = 0.1, 
                       histbins=40, colorspec = c("#7b3294", "#008837"),
                       facSel = NULL, filesep="/") {

# Wrapper for qcQvalue function to verify q-value analysis terms
# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for even descent to pi0 in the top left tuning plot
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data
#   3) Optional MDS plot (bonusMDS=TRUE) at q-value cut (default 0.1, see below)
#
# Arguments
#   norm_x: abundance data input into regression (eg, regressMatrix()). nrow(norm_x)==length(p_mat)
#   p_mat: matrix of p-values for each factor. Returned by regressMatrix()
#   q_list: S3 object; list of p-value lists for each factor. Returned by regressMatrix() #Feedback supposed to be q-value lists?
#   metadata: list of factors; each list element is a factor (eg, sex = c("m", "f"...)) #Feedback perhaps as.factor(c("m","f"))?
#             Vector order same as norm_x column order
#             was attribs
#   plotdata is a list of info relevant to labeling and saving the plot #Feedback with the following elements:
#     plotdir:  plot destination directory
#     plotbase:  base filename for the plot. Suggest: bias reduction method
#     plottitle:  title for all plots
#     plotSubtitle:  optional subtitle; suggest lm_expression
#   plot2file if TRUE, otherwise to studio...
#   bonusMDS: if TRUE plots an MDS of data at qcut (see below)
#   oneclass: metadata factor for labeling bonusMDS data; #Feedback is this actually a factor, or is it a vector?
#             if NULL or not specified, defaults to the first element of metadata
#   qcut: either a number between 0 and 1 used as an upper qvalue limit for MDS,
#     OR a number >1, assumed to be the top n qvalues to use in MDS
#   histbins: # bins for p-value histogram; may want to change if there is suspicious behavior
#   colorspec: optional range of colors for MDS plot. Defaults to hotpink/green!
#   facSel: optional vector of experimental factors to plot (rather than all factors) #Feedback to plot for the MDS plot? What if you want to plot based on second element of metadata list?
#   #Feedback what about filesep?
# Returns... Error message or... #Feeback LOL! Thank you for the good humor!!
#
  
  # Define factors to plot
  if(is.null(facSel)) {
    facUse = colnames(p_mat)
  } else {
    # Confirm elements of facSel in the model
    facCheck = facSel %in% colnames(p_mat)
    if (sum(facCheck) == 0) {
      message(paste0("facSel value ", facSel, " not in model.\n"))
      message("Plotting defaults.")
      facUse = colnames(p_mat)
    } else
    if (sum(facCheck) != length(facSel)) {
      facBad = facSel[which(!facCheck)]
      message(paste0("facSel value ", facBad, " not in model.\n"))
      facUse = facSel[which(facCheck)]
      message(paste0("Plotting ", facUse, "\n"))
    }
  }

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  # Establish the number of elements in the normalized data
  if(is.null(dim(norm_x))) {  # TRUE - norm_x is a vector
    norm_len = length(norm_x)
  } else { norm_len = nrow(norm_x) }

  ##### Plot data for each factor in facUse
  message("facUse length = ", length(facUse))
  orig_plotbase = plotdata$plotbase
  orig_plottitle = plotdata$plottitle
  for (factor in facUse) {
    skipFac = FALSE
    plotdata$plotbase = paste(orig_plotbase, factor, sep = "_")
    plotdata$plottitle = paste(orig_plottitle, factor, sep = ", ")
    # Check length of p-value list against norm_x
    if(nrow(p_mat) != norm_len) {
      message(sprintf("WARNING: Length of %s p_vals does not match norm_x",factor))
      skipFac = TRUE
    }
    # Check length of qvalues
    if(length(q_list[[factor]]$qvalues) != norm_len) {
      message(sprintf("WARNING: Length of %s q_list$qvalues does not match norm_x", factor))
      skipFac = TRUE
    }

    if (!skipFac) {
      message(sprintf("Factor %s to qcQvalues!", factor))
      qcQvalues(norm_x = norm_x, pvalue_v = p_mat[,factor], obj_qvalue = q_list[[factor]],
                qcut = qcut, attribs = metadata, oneclass = oneclass, plotdata = plotdata,
                colorspec = colorspec, histbins=histbins, plot2file = plot2file,
                p_hist=TRUE, q_plots=TRUE, MDS_plot=bonusMDS, filesep=filesep) #Feedback where does qcQvalues actually get defined?
    }
  }
  return("qcQ_wrapper done!")
}  # End of qcQ_wrapper function

