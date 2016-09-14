#####################################################################################################
# Set of functions useful in proteomics analysis.
#
# rqd.libraries    --- A tidy list of the libraries used for proteomics analysis.
# maplot           --- A plain vanilla MA plot with loess fit line
# hex.maplot       --- A 2D MA plot. Needs refinement.
# summary.plots    --- Creates box, scatter, and density plots
# qc.clusters      --- Plot heatmap and MDS - cuts data with small SD
#                  --- Needs some generalization for multiple variables
# scatterplot      --- X-Y Scatter plot
# mds.plot         --- MDS plot
#
# Author: Theresa Lusardi
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
maplot = function (v1, v2, v12names, plotdata, yrange = NULL, plot2file = FALSE) {
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

  plotID = 4
  plotDesc = sprintf('MA_%s.vs.%s', v12names$v1, v12names$v2) 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=5.4,units="in",res=144)
  }

  M = v1 - v2
  A = 0.5*(v1 + v2)
  ranges = fivenum(M)
  if (!is.null(yrange)) {
     ylimits = c(-max(abs(yrange)), max(abs(yrange)))
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

# ******** Summary ****************************************************************
summary.plots = function (rawmat, normmat, mynorm, samp.labels, samp.classes, colorspec, plotdata, plot2file = FALSE, histbins = 40, expand.2D = 5, filesep="/") {
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
  plotID = 1
  plotDesc = 'boxplot'
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=5,height=5.4,units="in",res=144)
  }
  boxplot(normmat, main = plotdata$plottitle, border=colvec, las = 2)
  if(plot2file) dev.off()

# Scatterplot of all data
  plotID = 2
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
  plotID = 3
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
  plotID = 7
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
         ylab='Standard deviation per gene',
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

# ******** Heatmap ****************************************************************
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
    normmat = sweep(normmat, 1, rowMeans(normmat,na.rm=T), '-')
    message(sprintf('centered on norm, max=%1.3f, min=%1.3f', max(normmat,na.rm=T), min(normmat,na.rm=T)))
  } else normmat = sweep(normmat, 1, rowMeans(rawmat,na.rm=T), '-')
  
  # calculate the color limits 
  # color limit for plotting: clim.pct percentile of data, divided in half (~dist to med)
  # try two different normalization methods...
#   c.lim1 = max(abs(quantile(normmat, c(1-clim.pct, clim.pct))))
  # inner limit
  c.lim = quantile( sapply(1:nrow(normmat),
                           function(x){ diff(range( normmat[x,], na.rm=T )) }),
                    probs=clim.pct, na.rm=T)/2
  # outer limit
  c.lim0 = max(abs(normmat),na.rm=T)
  
  # Determine a lower threshold for data of interest
  # If there is too much non-changing data, the histogram calculation will stall
  mySD = sd(normmat, na.rm = T)
  mymask = sapply(1:nrow(normmat),function(x){sd(normmat[x,],na.rm=T)}) > mySD
  message(sprintf('c.lim = %1.3f, mymask rows = %s, !mymask rows = %s',
                   c.lim, sum(mymask), sum(!mymask)))

  # Plot
  aheatmap(normmat[mymask,], color="-PiYG:64",
           breaks=c(-c.lim0, seq(from=-c.lim, to=c.lim, by=c.lim/31), c.lim0),
           annCol=attribs, labCol=colnames(normmat),
           main=plotdata$plottitle)

  if(plot2file) dev.off()

  plotID = 6
  plotDesc = 'MDS' 
  if(plot2file) {
  png(filename = sprintf('%s%i_%s_%s.png', plotdata$plotdir, plotID, plotdata$plotbase, plotDesc),
      width=4.5,height=4.9,units="in",res=300)
  }
  # factor to plot
  samp.classes = attribs[[oneclass]]
  # set up plotting colors
  u.samp.classes = unique(samp.classes)
  colmap = colorRampPalette(colorspec)
  colvec = colmap(length(u.samp.classes))
  # map colors back to samples
  colvec = colvec[as.numeric(as.factor(samp.classes))]
  # R doesn't resort on unique: saves order of occurrence
  u.col.classes = unique(colvec) 

  # plot MDS and capture numeric results
  obj_MDS = plotMDS(normmat, col=colvec, labels=colnames(normmat),
          top=sum(mymask), main=plotdata$plottitle, cex=.6 )

  # add legend
  axl = par("usr") # c(x1,x2,y1,y2) == extremes of user coords in plot region
  par(xpd=NA) # allow plotting anywhere on device
  legend(x=axl[2]-.025*abs(diff(axl[1:2])), y=axl[4], legend=u.samp.classes, col=u.col.classes, pch=16, cex=.6)

  if(plot2file) dev.off()
}  # end of qc.clusters function


# ******** Scatter ****************************************************************
scatterplot = function (normmat, attribs, plotdata, plot2file = FALSE) {
# This is intended to show individual replicates vs. average of experimental group
# It will create one plot per column, and a average vs average
# normmat is a matrix of experimental data in columns (with headers!)
# attribs is a vector of experimental categories (one entry per column in normmat)
# plotdata is a list of info relevant to labeling and saving the plot
#  plotdir:  plot destination directory
#  plotbase:  base filename for the plot
#  plottitle:  title for all plots
#

  plotID = '2a'
  plotlims = c(floor(min(normmat)), ceiling(max(normmat)))
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
