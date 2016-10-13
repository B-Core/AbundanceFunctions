#####################################################################################################
# Set of functions useful in proteomics analysis.
#
# rqd.libraries    --- A tidy list of the libraries used for proteomics analysis.
# maplot           --- A plain vanilla MA plot with loess fit line
# hex.maplot       --- A 2D MA plot. Needs refinement.
# summary.plots    --- Creates box, scatter, density and SD/mean plots
# qc.clusters      --- Plot heatmap and MDS for abundance data
#                  --- Makes avg-all ratios for heatmap, drops small SD rows
# qcQvalue         --- Plot p-value histogram, q-value auto-qc & MDS w q-cut
# plotRatios       --- Plot heatmap and MDS for ratio data
# scatterplot      --- X-Y Scatter plot
# makeMDSplot      --- MDS plot child function called by qc functions
# makeHeatmap      --- heatmap child function called by qc & design functions
# getCidFromHeatmap--- Extract index and cluster ID from heatmap based on a number of branches by which to cut
#
# Authors: Theresa Lusardi, Julja Burchard, and Mark Fisher
# Started - July 2016
#####################################################################################################

# ******** Libraries ****************************************************************
rqd.libraries = function (path) {
# path indicates the location of all other scripts used for analysis

  # Required libraries
  require(data.table)
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

# ******** Plain Vanilla MA Plot ****************************************************************
maplot = function (v1, v2, v12names, plotdata, yrange = NULL, plot2file = FALSE, plotIDOffset = 0) {
#  v1, v2:  vectors of equal length to compare by ma plot
#  v12names:  list of strings describing vectors v1 and v2
#    v1, v2 string names; will appear on the plot!
#  plotdata is a list of info relevant to labeling and saving the plot
#    plotdir:  plot destination directory
#    plotbase:  base filename for the plot
#    plottitle:  title for all plots
#  yrange:  if specified, sets the y-axis range
#    Note:  fivenumber or quantile vector will be included in the y-axis label, so it will be clear if the plot has been limited.
#  plot2file:  if true, plot to file as specified in plotdata
#  plotIDOffset:  offset value for the plot ID... default to zero

  plotID = 4 + plotIDOffset
  plotDesc = sprintf('MA_%s.vs.%s', v12names$v1, v12names$v2) 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=5.4,units="in",res=144)
  }

  M = v1 - v2
  A = 0.5*(v1 + v2)
  ranges = fivenum(M)
  if (!is.null(yrange)) {
     ylimits = c(-max(abs(yrange), na.rm = TRUE), max(abs(yrange), na.rm = TRUE))
    plot(A, M, pch = '.', ylim = ylimits)
  } else plot(A, M, pch = '.')
  title(xlab = 'A', ylab = 'M', main = plotdata$plottitle)
  mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                  ranges[1], ranges[2], ranges[3], ranges[4], ranges[5]),
        side = 2, line = 2, cex = 0.75)
  mtext (sprintf('%s\n%s', v12names$v1, v12names$v2),
         side = 1, line = 3, adj = 0, padj = 0, cex = .75)
  lfit = loess(M ~ A)  # Calculate a loess fit of the MA plot to show variance
  MC = predict (lfit, A)
  abline(h = 0, v = 0, col = 'grey60')  # Add an x axis at 0
  lines(A[order(A)], MC[order(A)], col = 'red')

  if(plot2file) dev.off()
} # end of maplot function

# ******** Hex Binned MA Plot ****************************************************************
hex.maplot = function (comp1, comp2, c12names, yrange = FALSE, plot2file = FALSE, plotdir = NULL, baseplotname = NULL) {
#  comp1, comp2:  vectors of equal length to compare by ma plot
#  c12names:  list of three strings describing plot title, comp1, and comp2
#  yrange:  if specified, sets the y-axis range
#    Note:  fivenumber or quantile vector will be included in the y-axis label, so it will be clear if the plot has been limited.
#  plot2file:  if true, plot to file as specified in plotdir and baseplotname

# This comes from https://cran.r-project.org/web/packages/hexbin/vignettes/hexagon_binning.pdf
# It needs some customization!

  M = comp1 - comp2
  A = 0.5*(comp1 + comp2)
  ranges = fivenum(M)

  hb1 = hexbin(A, M, xbins = 40)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1,2)))
  pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1))
  nb = plot(hb1, type = 'n', xlab = 'A', ylab = 'M', main = 'MA with pts', newpage = F)
  pushHexport(nb$plot.vp)

  grid.points(A, M, pch = 16, gp = gpar(cex = 0.4))
  popViewport()
  nb$hbin = hb1
  hexVP.abline(nb$plot.vp, h = 0, col = gray(0.6))
  hexMA.loess(nb)
  popViewport()
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
  plotMAhex(matrix(c(comp1, comp2), nrow = length(comp1)), newpage = F, main = "MA plot with Hex!", legend = 0)
  hexVP.abline(hb$plot.vp, h=0, col = gray(0.6))
  hexMA.loess(hb)
  popViewport()
  plot.new()

}

# ******** Summary *************************************************************
summary.plots = function (rawmat, normmat, mynorm, samp.labels, samp.classes, colorspec, plotdata, plot2file = FALSE, histbins = 40, expand.2D = 5, filesep="/", plotIDOffset = 0) {
# rawmat is the matrix of raw data. Should have any background addition and
#  be scaled according to normalized data
# normmat is a matrix of experimental data in columns (with headers!)
# mynorm is a label for the current normalization method
# samp.labels is a vector of brief, pithy display labels for each sample
# samp.classes is a vector of tags for each sample, used to pick colors
#  a unique tag per sample colors by sample
#  a shared tag by experimental groupings colors by experimental group
# colorspec is a vector of color specifiers for colorRampPalette  
# plotdata is a list of info relevant to labeling and saving the plot
#  plotdir:  plot destination directory
#  plotbase:  base filename for the plot
#  plottitle:  title for all plots
# histbins:  default # of bins for density plot
# expand.2D = multiplier of histbins for 2D histogram
# plotIDOffset:  a number to offset the plotID values... Can be used to organize
# Consider adding color! Maybe shade by experimental group

  # imports
  require(data.table)

# test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

# set up plotting colors
  u.samp.classes = unique(samp.classes)
  colmap = colorRampPalette(colorspec)
  colvec = colmap(length(u.samp.classes))
  # map colors back to samples
  colvec = colvec[as.numeric(as.factor(samp.classes))]
  # R doesn't resort on unique: saves order of occurrence
  u.col.classes = unique(colvec) 

# Boxplot of all data
  plotID = 1 + plotIDOffset
  plotDesc = 'boxplot'
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=5.4,units="in",res=144)
  }
  boxplot(normmat, main = plotdata$plottitle, border=colvec, las = 2)
  if(plot2file) dev.off()

# Scatterplot of all data
  plotID = 2 + plotIDOffset
  plotDesc = 'scatter'
  if(plot2file) {
    png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=5,height=5.4,units="in",res=144)
  }
  yl = range(normmat,na.rm=T); xl = range(rawmat,na.rm=T)
  plot(x=0, y=0, type="n", xlim=xl, ylim=yl, pch='.', xlab='', ylab = '')
  title(xlab='log2(raw + background)',
        ylab = paste(mynorm, 'normalization'),
        main = plotdata$plottitle)
  for(j in 1:ncol(rawmat) ){
    points(rawmat[,j], normmat[,j], pch='.', col=colvec[j])
  }
  legend(x="topleft", legend=u.samp.classes, col=u.col.classes, pch=16, cex=.5)
  if(plot2file) dev.off()

# Density plot
  plotID = 3 + plotIDOffset
  plotDesc = 'density'
  if(plot2file) {
    png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=5,height=5.4,units="in",res=144)
  }

  dl = NULL; yl = c(0,0); xl = range(normmat, na.rm=T)
  for( j in 1:ncol(normmat) ) {
    dl = c(dl,list(density(normmat[,j],from=xl[1],to=xl[2],n=histbins,na.rm=T)))
    dl[[j]]$y = dl[[j]]$y/sum(dl[[j]]$y)
    yl = range(yl,dl[[j]]$y)
  }
    plot(x=0,y=0,type="n",xlim=xl,ylim=yl,
         xlab=paste(mynorm,' normalized intensity'),
         ylab='Fraction of features (smoothed)',
         main = plotdata$plottitle)

  for(j in 1:ncol(normmat) ){
    lines(dl[[j]]$x,dl[[j]]$y,lty=(j%%ncol(normmat)) + 1, col=colvec[j])
  }
  legend(x="topright", legend=u.samp.classes, col=u.col.classes, lty=1, lwd=2, cex=.6)
                                                                                  
  if(plot2file) dev.off()

# SD plot
  plotID = 7 + plotIDOffset
  plotDesc = 'spread'
  if(plot2file) {
    png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=5,height=5.4,units="in",res=144)
  }

  # calculate group-wise SDs per gene
  sdmat = NULL; int_mat = NULL
  # find rows informative on SD for all classes
  row_mk = logical(nrow(normmat)); row_mk[] = TRUE
  for( myclass in u.samp.classes ){
    row_mk = row_mk & rowSums(!is.na(normmat[,samp.classes==myclass]))>1 }
  # calculate SDs and intensities
  for( myclass in u.samp.classes ){
    int_mat = cbind(int_mat, rowMeans(normmat[row_mk,samp.classes==myclass],na.rm=T))
    if(sum(samp.classes==myclass)>1){
      sdmat = cbind(sdmat, sapply(which(row_mk),
        function(x){ sd( normmat[x,samp.classes==myclass] ,na.rm=T) }) )
    } else {
      sdmat = cbind(sdmat, rep(0,nrow(normmat)) )
    }
  }
  colnames(sdmat) = u.samp.classes; colnames(int_mat) = u.samp.classes
  rownames(sdmat) = rownames(normmat)[row_mk]

  # colors
  # bright version of each color: increase to cmult% of distance to 100%
  cmult = 0.8
  clower = rgb(t( cmult*(255 - col2rgb(colvec)) + col2rgb(colvec) )/255 )
  # build color list
  collist = NULL; ncols = 64
  for(j in 1:ncol(sdmat) ){
    # color axis for each group, lighter at low end & mostly transparent 
    col1 = unique(clower[samp.classes==colnames(sdmat)[j] ])[1]
    col2 = unique(colvec[samp.classes==colnames(sdmat)[j] ])[1]
    collist = c(collist, list(colorRampPalette( colors=c( adjustcolor(col1, alpha.f=.85), adjustcolor( col2, alpha.f=0.85) ), alpha=T)(ncols) ) )
  }
  # 2D histogram
  nbins = histbins*expand.2D
  yl = range(sdmat, na.rm=T)
  ybin = seq(from=yl[1]*(1-1/nbins),to=yl[2]*(1-1/nbins),length.out=nbins+1)
  yplt = filter(ybin,filter=c(.5,.5),sides=1)[2:(nbins+1)]
  xl = range(normmat[row_mk,], na.rm=T)
  xbin = seq(from=xl[1]*(1-1/nbins),to=xl[2]*(1-1/nbins),length.out=nbins+1)
  xplt = filter(xbin,filter=c(.5,.5),sides=1)[2:(nbins+1)]
  cmax = log2(nrow(sdmat)/nbins)
  # calculate x, y, color
  freq = NULL;
  for(i in 1:ncol( sdmat ) ){
    # frequency of abundance vs SD as binned
    freq[[i]] = as.data.table(table( findInterval(int_mat[,i],xbin), findInterval(sdmat[,i],ybin) ))
    freq[[i]] = freq[[i]][,V1:=as.numeric(V1)]
    freq[[i]] = freq[[i]][,V2:=as.numeric(V2)]
    freq[[i]] = freq[[i]][,N:=log2(N+1)]
    # transform to color numbers
    freq[[i]] = freq[[i]][,Nn := round((ncols-1)*N/cmax)+1]
    freq[[i]][Nn>ncols,Nn := ncols]
    # trim empty bins
    freq[[i]] = freq[[i]][N>0,]
  }
  # make single x, y, color vecs -- alpha seems to have limits
  myx = NULL; myy = NULL; myc = NULL
  for(i in 1:ncol( sdmat ) ){
    myx = c(myx,freq[[i]]$V1)
    myy = c(myy,freq[[i]]$V2)
    myc = c(myc,collist[[i]][freq[[i]]$Nn])
  }
  randx = sample(length(myx))
  # plot
    plot(x=0,y=0,type="n",xlim=xl,ylim=yl,
         xlab=paste(mynorm,' normalized intensity'),
         ylab='Standard deviation',
         main = plotdata$plottitle)

  points(xplt[myx[randx]],yplt[myy[randx]],col=myc[randx],pch=15,cex=.6)

  # add loess fits as in Cope et al. Bioinformatics 20:323-331 (2004), Fig.2
  for(i in 1:ncol(sdmat) ){
    lfit = lowess( y = sdmat[,i], x = int_mat[,i], f=0.2, delta=(1/histbins)*diff(range(int_mat[,i])) )
    # plot heavy lines a bit darker than regular colors
    # lowess returns list elements x and y in plotting order
    lines( lfit$x, lfit$y, col=adjustcolor(u.col.classes[i],red.f=.75,green.f=.75,blue.f=.75), lwd=3 )
  }

  # add legend
  legend(x="topright", legend=u.samp.classes, col=u.col.classes, pch=16, cex=.6)

  if(plot2file) dev.off()

}

# ******** Abundance QC MDS & Heatmap ******************************************
qc.clusters = function (normmat, rawmat, attribs, oneclass, plotdata,
                        colorspec, center = 'norm', clim.pct, mask = NULL, 
                        plot2file = FALSE, filesep='/') {
# This function considers only a single clustering variable... Will need to be gneralized
# Data will be centered around the average of all data in the normalized data set
# As a QC plot, want to filter out non-changing data, otherwise the alg gets slow
# Uses the attributes to cluster the data; could make this a longer list...
#  normmat:  data matrix 
#  rawmat:  raw data matrix (optionally used to center normmat)
#  attribs:  list of sample classifications to be tracked in clustering
  # each list element contains a string vector with one label per sample
#  oneclass: string name of attribs element to be used in MDS plot
#  plotdata is a list of info relevant to labeling and saving the plot
#    plotdir:  plot destination directory
#    plotbase:  base filename for the plot
#    plottitle:  title for all plots
#  colorspec is a vector of color specifiers for colorRampPalette  
#  center:  center data on 'norm' for normalized average
#           any other value will trigger centering on the average of rawmat
#  clim.pct:  0:1 - fraction of data to limit max color
#  mask:  NULL default masks rows with SD > SD of entire data set; 
#         Numeric vector of length 1 masks rows with SD > value
#         Vector of length nrows(normmat) is a custom mask   ***** TBImplemented
#  plot2file:  if true, plot to file as specified in plotdata

  # imports
  require(NMF)
  require(limma)

# test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }


  plotID = 5
  plotDesc = 'Heatmap' 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=7,units="in",res=300)
  }

  # Center data
  if(center == 'norm') {
    ratiomat = sweep(normmat, 1, rowMeans(normmat,na.rm=T), '-')
    message(sprintf('centered on norm, max=%1.3f, min=%1.3f', max(ratiomat,na.rm=T), min(ratiomat,na.rm=T)))
  } else ratiomat = sweep(normmat, 1, rowMeans(rawmat,na.rm=T), '-')
  
  # Determine a lower threshold for data of interest
  # If there is too much non-changing data, the cluster calculation will stall
  # As this is a QC plot, limit to genes more variable than average
  mySD = sd(ratiomat, na.rm = T)
  mymask = sapply(1:nrow(ratiomat),function(x){sd(ratiomat[x,],na.rm=T)}) > mySD
  message(sprintf('mySD = %1.3f, mymask rows = %s, !mymask rows = %s',
                   mySD, sum(mymask), sum(!mymask)))

  ah_ls = makeHeatmap(normmat=normmat[mymask,], ratiomat=ratiomat[mymask,], attribs=attribs, plottitle = plotdata$plottitle, clim.pct=clim.pct)

  if(plot2file) dev.off()


  plotID = 6
  plotDesc = 'MDS' 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5.4,height=5.4,units="in",res=300)

  obj_MDS = makeMDSplot(normmat=normmat, attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, ngenes=sum(mymask))

  if(plot2file) dev.off()

  invisible(list(ah_ls=ah_ls, obj_MDS=obj_MDS))
  }
}  # end of qc.clusters function


# ******** Q-value QC **********************************************************
qcQvalues = function (norm_x, pvalue_v, obj_qvalue, qcut, attribs, oneclass, 
                      plotdata, colorspec, histbins=40, plot2file = FALSE, 
                      filesep='/') {
# This is a check on the proper running of q-value analysis
# plots:
#  1) histogram of p-values input to q-value. 
#     these must NOT have a right peak or qvalue() will return spurious results
#  2) qvalue's default plots, with full qvalue range c(0,1) plotted
#     these should show an even descent to pi0 in the top left tuning plot,
#     and a slow or steep rise in q-values in remaining plots depending on 
#     resolving power of data
#  3) MDS plot restricted by q-value cut.
#  norm_x:  abundance data, vector or matrix. nrow(normmat)==length(pvalue_v) 
#  pvalue_v: vector of p-values previously used as input to qvalue()
#  obj_qvalue: S3 object of qvalue class (a list!) returned by qvalue()
#  qcut: either a number between 0 and 1 used as an upper qvalue limit for MDS,
#     OR a number >1, assumed to be the top n qvalues to use in MDS
#  attribs:  list of sample classifications to be tracked in clustering
#     each list element contains a string vector with one label per sample
#  oneclass: string name of attribs element to be used in MDS plot
#  plotdata is a list of info relevant to labeling and saving the plot
#     plotdir:  plot destination directory
#     plotbase:  base filename for the plot. Suggest: bias reduction method
#     plottitle:  title for all plots
#  colorspec is a vector of color specifiers for colorRampPalette  

  # constants
  fudgefac = 2 # fold excess permitted in rows recovered by top n qcut when specified as a quantile of qvalues

  # imports
  require(qvalue)
  require(limma)

  # argument tests
  if( is.null(dim(norm_x)) ){ # not matrix
    norm_len = length(norm_x)
  } else { norm_len = nrow(norm_x) 
  }
  if( !exists("qvalues", obj_qvalue) ){
    stop("qvalue object does not contain qvalues element")
  }
  if( norm_len != length(pvalue_v) | norm_len != length(obj_qvalue$qvalues) ){
    stop("Data matrix has ",norm_len," rows, p-value vector has ",length(pvalue_v)," elements and qvalue object has ",length(obj_qvalue$qvalues)," q-values, but all must be the same size")
  }

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  plotID = '8p'
  plotDesc = 'p.value_histogram' 
  if(plot2file) {
  png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=5.4,units="in",res=300)
  }
  hist(pvalue_v, nclass=histbins, main=plotdata$plottitle)

  if(plot2file) dev.off()


  plotID = '8q'
  plotDesc = 'q.value_QC.plot.array' 
  if(plot2file) {
  png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=7,height=5.4,units="in",res=300)
  }
  plot(obj_qvalue, rng=c(0,1), cex.axis=.6)

  if(plot2file) dev.off()
  

  # make first-draft MDS plot if abundance data are a matrix 
  #  (multiple samples to compare)
  obj_MDS = NULL
  if( !is.null(dim(norm_x)) ){
    plotID = '6q0'
    plotDesc = 'MDS_q.value_QC' 
    if(plot2file) {
    png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
        width=5.4,height=5.4,units="in",res=300)
    }

    # data to plot
    if( !is.numeric(qcut) | qcut<0 | qcut>length(obj_qvalue$qvalues) ){
      stop(paste(qcut," must be a number > 0 and <= nrow(data)"))
    } else if( qcut <=1 ) { # assume this is a q-value on which to cut
      mymask = obj_qvalue$qvalues < qcut & !is.na(obj_qvalue$qvalues)
      titlestr= paste("qvalue <",signif(qcut,2) )
    } else { # assume this is a number of top genes by qvalue on which to cut
      q_quan = quantile(obj_qvalue$qvalues, probs=qcut/length(obj_qvalue$qvalues), na.rm=T)
      mymask = obj_qvalue$qvalues <= q_quan & !is.na(obj_qvalue$qvalues)
      titlestr= paste("top",sum(mymask),"q-values <=",signif(q_quan,3) )
      # test for lack of resolution in q-values; if so use p-values instead
      if( sum(mymask)>(qcut*fudgefac) ) { #identical qvalues at cut increase sum
        p_quan = quantile( pvalue_v, probs=qcut/length(pvalue_v), na.rm=T)  
        mymask = pvalue_v <= p_quan & !is.na(obj_qvalue$qvalues)
        titlestr= paste("top",sum(mymask),"p-values <=",signif(p_quan,3) )
      }
    }

    obj_MDS = makeMDSplot(normmat=norm_x[mymask,], attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, subtitle=titlestr, ngenes=sum(mymask))

    if(plot2file) dev.off()
  }

  invisible( list(rowmask=mymask,obj_MDS=obj_MDS) )

} # end of qcQvalue()


# ******** MDS & Heatmap plots for ratios **************************************
plotRatios = function (ratiomat, attribs, oneclass, plotdata, colorspec, 
                       rowmask=NULL, tag="selected",
                       NAtol=ncol(ratiomat)/2, SDtol=0.01, 
                       clim.pct=0.99, clim_fix=NULL, 
                       plot2file = FALSE, filesep='/', 
                       cexRow=0.00001, png_res=300) {
# This is intended to display differential expression 
# It will create a heatmap of ratios and an MDS plot from the same ratios
# ratiomat: matrix of ratio data in columns (with headers!)
# rowmask: optional mask of rows to plot
#   default is to create an all-TRUE rowmask of length nrow(ratiomat)
# NAtol: max number of NAs tolerated; rows with more are masked out
# SDtol: min SD within row tolerated; rows with less are masked out
# attribs: list of sample classifications to be tracked in clustering
#   each list element contains a string vector with one label per sample
# tag: word or phrase to indicate what's special about _these_ ratios
# oneclass: string name of attribs element to be used in MDS plot
# plotdata: list of info relevant to labeling and saving the plot
#   plotdir:  plot destination directory
#   plotbase:  base filename for the plot
#   plottitle:  title for all plots
#  colorspec is a vector of color specifiers for colorRampPalette  
#  clim.pct:  0:1 - fraction of data to limit max color
#  clim_fix: if set, max abs ratio to show; data>clim_fix are shown as clim_fix
#  cexRow: rowlabel size for heatmap, set to ~invisible by default
#  png_res: resolution of saved png files in dpi; default 300

  # imports
  require(NMF)

  # arguments

  # test plotdir for filesep; add if absent
  if( !grepl(paste0(filesep,'$'),plotdata$plotdir) ){
    plotdata$plotdir = paste0(plotdata$plotdir,filesep)
  }

  # check for rowmask; create if NULL
  if( is.null(rowmask) ){
    rowmask = logical(nrow(ratiomat)); rowmask[]= TRUE
  }
  # adjust rowmask if needed to keep problem rows out of plots
  NAmask = rowSums(is.na(ratiomat)) <= NAtol
  SDmask = sapply(1:nrow(ratiomat), 
                  function(x){sd(ratiomat[x,],na.rm=T)}) >= SDtol
  SDmask[ is.na(SDmask) ] = FALSE # in case of all-NA rows
  rowmask = rowmask & NAmask & SDmask
  # report selection size
  message(sprintf('selected rows = %s, excess NA rows = %s, low SD rows = %s',
                   sum(rowmask), sum(!NAmask), sum(!SDmask) ))

  plotID = '5q1'
  plotDesc = paste('Heatmap', tag, sep="_")
  if(plot2file) {
  png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc), width=5,height=7,units="in",res=png_res)
  }

  # plot heatmap
  ah_ls = makeHeatmap( ratiomat=ratiomat[rowmask,], attribs=attribs, plottitle = plotdata$plottitle, clim.pct=clim.pct, clim_fix=clim_fix, cexRow=cexRow)
  
  if(plot2file) dev.off()


  plotID = '6q1'
  plotDesc = paste('MDS', tag, sep="_") 
  if(plot2file) {
  png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc), width=5.4,height=5.4,units="in",res=png_res)
  }

  # plot MDS
  obj_MDS = makeMDSplot(normmat=ratiomat[rowmask,], attribs=attribs, oneclass=oneclass, colorspec=colorspec, plottitle=plotdata$plottitle, ngenes=sum(rowmask))

  if(plot2file) dev.off()


  # return processed data
  invisible( list( ah_ls=ah_ls, obj_MDS=obj_MDS) )
  
}


# ******** Scatter ****************************************************************
scatterplot = function (normmat, attribs, plotdata, plot2file = FALSE, plotIDOffset = 0) {
# This is intended to show individual replicates vs. average of experimental group
# It will create one plot per column, and an average vs average
# normmat is a matrix of experimental data in columns (with headers!)
# attribs is a vector of experimental categories (one entry per column in normmat)
# plotdata is a list of info relevant to labeling and saving the plot
#  plotdir:  plot destination directory
#  plotbase:  base filename for the plot
#  plottitle:  title for all plots
#

  plotID = ifelse(plotIDOffset == 0, '2a', paste0(plotIDOffset + 2, "a"))
  plotlims = c(floor(min(normmat, na.rm=TRUE)), ceiling(max(normmat, na.rm=TRUE)))
  expts = unique(attribs)
  n.expts = length(expts)

  for (expt in expts)  {
    avg = rowMeans(normmat[, attribs == expt])
    for (i in 1:ncol(normmat))  {

      if (attribs[i] == expt) {  # Only plot the data in the experimental group
        if(plot2file) {
          plotDesc = sprintf('Scatter.average_%s.vs.replicate_%s', expt, colnames(normmat)[i]) 
          png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID,
              plotdata$plotbase, plotDesc), width=5,height=5.4,units="in",res=144)
        }  # end of plotfile setup
      
        plot(avg, normmat[,i], pch='.', xlab='', ylab = '', xlim = plotlims, ylim = plotlims)
        avgranges = fivenum(avg)
        repranges = fivenum(normmat[, i])
        title(xlab = sprintf('Average %s', expt),
              ylab = sprintf('Replicate %s', colnames(normmat)[i]),
              main = plotdata$plottitle)
        mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                      repranges[1], repranges[2], repranges[3], repranges[4], repranges[5]),
               side = 2, line = 2, cex = 0.75)
        mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                      avgranges[1], avgranges[2], avgranges[3], avgranges[4], avgranges[5]),
               side = 1, line = 2, cex = 0.75)

        if(plot2file) dev.off()
      }  # End of plot replicates test
    }  # end of i loop through all columns
  }  # end of expt loop through experimental groups

  # Now plot averages against each other
  for (i in 1:(n.expts-1)) {
    for (j in (i+1):n.expts) {
      avgx = rowMeans(normmat[, attribs == expts[i] ])
      avgy = rowMeans(normmat[, attribs == expts[j] ])
      xrange = fivenum(avgx)
      yrange = fivenum(avgy)
      
      if(plot2file) {
        plotDesc = sprintf('Scatter.average_%s.vs.average%s', expts[i], expts[j])
        png(filename = sprintf('%s%s_%s_%s.png', plotdata$plotdir, plotID,
                               plotdata$plotbase, plotDesc),
            width=5,height=5.4,units="in",res=144)
      }  # end of plotfile setup

      plot(avgx, avgy, pch='.', xlab='', ylab='', xlim = plotlims, ylim=plotlims)
      title(xlab = sprintf('Average %s', expts[i]),
            ylab = sprintf('Average %s', expts[j]),
            main = plotdata$plottitle)
      mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                    xrange[1], xrange[2], xrange[3], xrange[4], xrange[5]),
            side = 2, line = 2, cex = 0.75)
      mtext(sprintf('(%1.2f, %1.2f, %1.2f, %1.2f, %1.2f)',
                    yrange[1], yrange[2], yrange[3], yrange[4], yrange[5]),
            side = 1, line = 2, cex = 0.75)
      if(plot2file) dev.off()
    }  # end of j loop
  }  # end of i loop
}  # end of scatterplot function


##### **** Child functions that perform subtasks **************************#####


makeMDSplot = function (normmat, attribs, oneclass, colorspec, plottitle, 
                        subtitle=NULL, ngenes=NULL) {
# This function considers only a single clustering variable in coloring MDS plot... Will need to be gneralized
# Uses the matrix values to cluster the data
#  normmat:  data matrix, with unique & informative colnames 
#  attribs:  list of sample classifications to be tracked in clustering
#    each list element contains a string vector with one label per sample
#  oneclass: string name of attribs element to be used in MDS plot
#  colorspec is a vector of color specifiers for colorRampPalette  
#    colors are generated per sample by their clustering variable values
#  plottitle:  title for all plots
#  subtitle: optional subtitle to add below title on this plot
#  ngenes: number of "top" genes for plotMDS to select for distance calculation

  # imports
  require(limma)

  # factor to plot
  samp.classes = attribs[[oneclass]]
  # set up plotting colors
  u.samp.classes = unique(samp.classes)
  colmap = colorRampPalette(colorspec)
  colvec = colmap(length(u.samp.classes))
  # map colors back to samples
  colvec = colvec[as.numeric(as.factor(samp.classes))]
  # R doesn't re-sort on unique: saves order of occurrence
  u.col.classes = unique(colvec) 

  # number of "top" genes used for sample-sample distance calculation
  if( is.null(ngenes) ) { ngenes = nrow(normmat) }

  # plot MDS and capture numeric results
  par( mar=c(5,4,4,4)+0.1 ) #default mar=c(5,4,4,2)+0.1; margin in lines
  obj_MDS = plotMDS(normmat, col=colvec, labels=colnames(normmat),
          top=ngenes, main=plottitle, cex=.6 )
  if( !is.null(subtitle) ){ mtext(subtitle, cex=.8) }

  # add legend
  axl = par("usr") # c(x1,x2,y1,y2) == extremes of user coords in plot region
  # allow plotting anywhere on device, and widen right margin
  par(xpd=NA) 
  legend(x=axl[2]-.025*abs(diff(axl[1:2])), y=axl[4], legend=u.samp.classes, col=u.col.classes, pch=16, cex=.6)

  return(obj_MDS)

} # end makeMDSplot


makeHeatmap = function (ratiomat, attribs, plottitle, normmat=NULL,
                        clim.pct=.99, clim_fix=NULL, colorbrew="-PiYG:64", 
                        cexRow=0.00001, 
                        cexCol=min(0.2 + 1/log10(ncol(ratiomat)), 1.2),
                        labcoltype=c("colnames","colnums") ) {
# This function makes a heatmap of ratio data, displaying experimental design values as tracks
# Uses the matrix values to cluster the data
#  normmat:  abundance data matrix, optionally used to set color limits
#   set to null to use ratiomat for color limits
#  ratiomat:  ratio data matrix, optimally derived from normmat,
#    with unique & informative colnames 
#  attribs:  list of sample classifications to be tracked in clustering
#    each list element contains a string vector with one label per sample
#    set to NA to omit
#  colorbrew is a colorBrewer string specifying a heatmap color scale 
#    colors are generated per sample by their clustering variable values
#  plottitle:  title for all plots
#  clim.pct:  0:1 - fraction of data to limit max color
#  clim_fix: if set, max abs ratio to show; data>clim_fix are shown as clim_fix
#  cexRow: rowlabel size for heatmap, set to ~invisible by default
#  cexCol: collabel size for heatmap, set to aheatmap default by default
#  labcoltype: colnames to show ratiomat column names, colnums to show col #s

  # imports
  require(NMF)

  # set column labels
  if( any(grepl("colnames",labcoltype,ignore.case=T)) ){
    labCol = colnames(ratiomat)
  } else {
    labCol = 1:ncol(ratiomat)
  }
  # calculate the number of colors in each half vector
  if(length(colorbrew)>1){ # color vector given
    halfbreak = length(colorbrew)/2
  } else if(any(grepl('^[^:]+:([0-9]+)$',colorbrew)) ){
    halfbreak = as.numeric(sub('^[^:]+:([0-9]+)$','\\1',colorbrew))/2
  } else { # aheatmap also takes a single colorspec and makes 2-color map
    halfbreak = 1
  }
  if(halfbreak>1){halfbreak=halfbreak-1} # make room for outer limits!

  # calculate the color limits 
  # if normmat is provided: clim.pct percentile of data, divided in half (~dist to med) gets color scale
  # if no normmat: clim.pct percentile of ratio data gets color scale
  #    values outside this range get max color
  # inner limit
  if( !is.null(normmat) ){
    c.lim = quantile( sapply(1:nrow(normmat),
                           function(x){ diff(range( normmat[x,], na.rm=T )) }),
                    probs=clim.pct, na.rm=T)/2
  } else { 
    c.lim = quantile( ratiomat, probs=clim.pct, na.rm=T )
  }
  # outer limits
  if( !is.null(clim_fix) ){
    ratiomat[ratiomat>clim_fix] = clim_fix
    ratiomat[ratiomat<=-clim_fix] = -clim_fix
  }
  c.lim0 = max(abs(ratiomat),na.rm=T)
  c.min0 = min(ratiomat,na.rm=T) # true up lower limit -- make this optional?
  
  # make vector of colorbreaks
  if(halfbreak>1) { # use inner and outer limits
    colorbreaks=c(-c.lim0, seq(from=-c.lim, to=c.lim, by=c.lim/halfbreak), c.lim0)
  } else { # use outer limits only
    colorbreaks=c(c.min0,(c.min0+c.lim0)/2, c.lim0)
  }

  # Plot
  ah_ls = aheatmap(ratiomat, cexRow=cexRow, 
           color=colorbrew, breaks=colorbreaks,
           annCol=attribs, labCol=colnames(ratiomat),
           main=plottitle)

  return(ah_ls)
} # end makeHeatmap


getCidFromHeatmap = function(aheatmapObj, numSubTrees, cutByRowLogic=T){
  # Return a list of length 2 containing 1) a named vector of cluster IDs when cutting the dendrogram by numSubTrees and 2) a vector of elements in the order as they appear in the dendrogram (left to right for columns and top to bottom for rows).
  # 
  # aheatmapObj: an object returned by the aheatmap() function from the NMF package
  # numSubTrees: a numerical argument specifying the number of subclades/subtrees to split the dendrogram row or column into.
  # cutByRowLogic: a boolean specifying whether you want to return cluster IDs and orders for rows (TRUE) or columns (FALSE).
  # returns a list of length 2 containing 1) a named vector of cluster IDs when cutting the dendrogram by numSubTrees and 2) a vector of elements in the order as they appear in the dendrogram (left to right for columns and top to bottom for rows).
  # examples
  # normmat = NULL
  # ratiomat = structure(c(-2.96...
  # attribs = list(ARPEcellStatus=c(rep("No",1),rep("ARPE",1),rep("No",1)))
  # plottitle = "CommonCellTypesInOurImportantTSSs"
  # clim.pct = .99
  # clim_fix = NULL
  # colorbrew = "-PiYG:64"
  # cexRow=0.01
  # labcoltype=c("colnames","colnums")
  # source("abundance_functions.R")
  # aheatmapObj = makeHeatmap (normmat=NULL, ratiomat=ratiomat, attribs=attribs, plottitle=plottitle, clim.pct=clim.pct, colorbrew=colorbrew, cexRow=cexRow, labcoltype=labcoltype)
  # numSubTrees = 4
  # cutByRowLogic = T
  # test_ls = getCidFromHeatmap(aheatmapObj, numSubTrees, cutByRowLogic)
  # test_ls$clusterIDByElement
  #  #1  2  3  4  5  6  7  8  9 10 
  #  #1  2  2  3  3  3  2  3  3  4 
  # test_ls$elementOrderByDendrogram
  # #[1]  3  2  7  4  9  8  6  5  1 10

  if(cutByRowLogic == TRUE){
    clusterIDByElement = cutree(as.hclust(aheatmapObj$Rowv), k = numSubTrees)
    elementOrderByDendrogram = rev(as.hclust(aheatmapObj$Rowv)$order)
    return_ls = list(clusterIDByElement, elementOrderByDendrogram)
    names(return_ls) = c("clusterIDByElement", "elementOrderByDendrogram")
  }else{
    clusterIDByElement = cutree(as.hclust(aheatmapObj$Colv), k = numSubTrees)
    elementOrderByDendrogram = as.hclust(aheatmapObj$Colv)$order
    return_ls = list(clusterIDByElement, elementOrderByDendrogram)
    names(return_ls) = c("clusterIDByElement", "elementOrderByDendrogram")
  }
  return(return_ls)
}


# ******** Code Snippets! ****************************************************************
  # some code from Mark to colorize plots!
#  require(RColorBrewer)
#  cool_cols = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(tx$Inv_sample_name))) #or however else you get the number of samples
#  for (i in 1:length(Expression_all_list)){
#    mymat = as.matrix(Expression_all_list[[i]])
#    axl = range(c(range(Expression_non_normalized_affy_obj, na.rm = T),range(2^mymat, na.rm = T)),na.rm=T)
#    xl = axl
#    tst = as.matrix(Expression_non_normalized_affy_obj)
#    plot(2^tst, mymat,log='x',pch='.',xlim=xl,ylim=log2(axl), main=paste("Raw vs. norm plot: ", names(Expression_all_list)[i], sep=""), xlab="Non-normalized expression", ylab=paste(names(Expression_all_list)[i], " expression", sep=""))
#    for(j in 1:length(unique(tx$Inv_sample_name))){
#      current_tmnt = unique(tx$Inv_sample_name)[j]
#      cols_of_treatment = which(tx$Inv_sample_name==current_tmnt)
#      for (k in cols_of_treatment){
#        points(2^tst[,k], mymat[,k], col=cool_cols[j], pch=".")
#      }
#    }
#  }
