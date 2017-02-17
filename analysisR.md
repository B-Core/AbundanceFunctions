
analysis.R.md
=====================
author: Janice
date: 10-25-2016
output: html_document

analysis.R runs prebuilt functions on RNAseq alignments

R session should be run in /proj/<projectname>/results/STAR/ directory or setwd("/proj/<projectname>/results/STAR/") while in a current R session.

### Import libraries and CompBio Repos:
- AbundanceFunctions
- AssociationFunctions

This should be a separate updated but static copy of the repo so you are running and testing your script on a static instance. That way drastic updates will not effect your current analysis.

```{r}
require(data.table)
require(NMF)
require(affy)
require(limma)
require(AnnotationDbi)
source("<pathtoproj>/code/AbundanceFunctions/RNAseq_abundance_functions.R")
source("<pathtoproj>/code/AbundanceFunctions/Proteomics.Tools.x1.R")
source("<pathtoproj>/code/AbundanceFunctions/readENSgtf.R")
source("<pathtoproj>/projs/RSchweitzer_RNA160803RS/code/AbundanceFunctions/processData.R")
source("<pathtoproj>/projs/RSchweitzer_RNA160803RS/code/AssociationFunctions/gs.wrapper.R")
```
### Constants
Set some constants 
for later use. 


```{r}
hclust.limit = 2^16 # max data rows for hclust
hc.frac.cut = 0.75; # quantile of data distribution reqd in one group's worth of data, if too many rows for hclust() show minimal expression
# regression parameter
na.lim = 0 # max NAs per row tolerated by lm() at least in some cases
# plotting colors
colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7))

```

### Project annotation
```
# custom project title
proj.title = "RNA160803RS"
# custom annotation
md.file = "/home/exacloud/lustre1/CompBio/projs/RSchweitzer_RNA160803RS/data/expt.design.txt"
```


### Gene file annotation.
Make sure you have the proper GTF used for alignment and the correct release.
TaxID can be found www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
```
taxID=10090
gene2ENSfile="/<>/CompBio/genomic_resources/anno_sources/ncbi/gene2ensembl.gz"
gene2ENS.col = c("taxID","EntrezID","Gene","RefSeqTranscript","EnsemblTranscript","RefSeqProtein","EnsemblProtein")
gtfFile = "/<>/CompBio/genomic_resources/gtf/mm10/release-85/Mus_musculus.GRCm38.85.gtf"
gtf.feature = "gene"
gtf.orig.col = c("gene_id","gene_name","gene_biotype")
gtf.col = c("Gene","Symbol","biotype")
```



### Setup for reading in STAR output

project directory where star output was stored.
```
readir = '/<pathtoproj>/<projname>/results/STAR'
```


Extract the names of the samples from the output names.
Use `readpattern` to specifiy the star output to read in usually all the files, depends on groups.
Use `label.from.colname` as *unique* identifiers of each sample.

```
readpattern='0336.*Ovation'
useme.cols = 'RS'
label.from.colname='^.*RS_([^R]+).*$'![][]
```

This helps annotate your STAR.data object. These are the "conditions"/"labels"/"groups" for your sample names. Use 'annCol.label' as group identifier. These can be multiple conditions. I'm not sure how to 
add more than one yet.
```
annCol.label ='^.*RS_([^R]+).*$'
annCol.names='group'
annCol.lmBy='group'
```

### Read in STAR output *finally*
Creates STAR.data "list" and writes to current directory (results/STAR) the following csv files:
- raw.gene.counts.FwdReads.csv
- raw.gene.counts.Reads.csv
- raw.gene.counts.RevReads.csv
- raw.SJ.counts.UniqueReads.csv


```
STAR.data = read_STAR(useme.cols=useme.cols,label.from.colname=label.from.colname,annCol.label=annCol.label,annCol.names=annCol.names,annCol.normBy=NULL,annCol.lmBy=annCol.lmBy,readpattern=readpattern,unstranded.col = list(gene.counts = c(1:4),SJ.counts = c(1:3,7)))

```
* STAR.data structure=> list of list of list of list of matrices*

  <blockquote> STAR.data[
    expt.design,
    myreads,
    LoM.raw[
      gene.counts[
        Reads(matrix),
        FwdReads(matrix),
        RevReads(matrix)
      ],
      SJ.counts[
      UniqueReads(matrix)
      ]
    ]
  ]</blockquote>

> STAR.data$expt.design
> $sample
>[1] "IY27_P7_1_Ovation_" "IY27_P7_1_Ovation_" "IY27_P7_1_Tru"
> [4] "IY27_P7_1_Tru"      "IY27_P7_2_Ovation_" "IY27_P7_2_Ovation_"
> [7] "IY27_P7_2_Tru"      "IY27_P7_2_Tru"      "IY27_P7_3_Ovation_"
>[10] "IY27_P7_3_Ovation_" "IY27_P7_3_Tru"      "IY27_P7_3_Tru"
>
>attr(,"lmBy")
>[1] "sample"
>



Object that holds the labels/groups/conditions of your project. Might need parsing from <proj>/data/expt.design.txt. This was done manually for now.

```
annCol=STAR.data[["expt.design"]]
```

**Test to confirm readpattern='^RNA16'**

```
i = grep('gene', names(STAR.data$LoM.raw),ignore.case=T)[1]
fwdmk = grepl('fwd|forward|sense', names(STAR.data$LoM.raw[[i]]), ignore.case=T)
revmk = grepl('rev|anti', names(STAR.data$LoM.raw[[i]]), ignore.case=T)
f = which(fwdmk)[1]; r = which(revmk)[1]; u = which(!fwdmk & !revmk)[1]
if( any(is.na( c(f,r,u) )) ){
  stop('Gene-based data must provide forward, reverse and unstranded types')
}
testv = log2(colSums(STAR.data$LoM.raw[[i]][[u]],na.rm=T)+1) - log2(colSums(STAR.data$LoM.raw[[i]][[f]],na.rm=T)+colSums(STAR.data$LoM.raw[[i]][[r]],na.rm=T)+1)
if( sum(round(testv,0))>1 ){
    stop("Forward and reverse reads don't add up to total reads!")
} else {
  testv =  log2(colSums(STAR.data$LoM.raw[[i]][[f]],na.rm=T)) - log2(colSums(STAR.data$LoM.raw[[i]][[r]],na.rm=T)+1)
  if( sum(round(testv,0))>1){            # fwd strand excess over rev
    myreads = names(STAR.data$LoM.raw[[i]])[f]
    cat(" Using", myreads, "counts","\n")
  } else if( sum(round(testv,0))<=-1){   # rev strand excess over fwd
    myreads = names(STAR.data$LoM.raw[[i]])[r]
    cat(" Using", myreads, "counts","\n")
  } else {
    myreads = names(STAR.data$LoM.raw[[i]])[u] # no strand inequality
    cat(" Using", myreads, "counts","\n")
  }
}
```

### Parse of Ensembl GTF used for processing
```
gtf = readENSgtf(filename=gtfFile)
genes.gtf = gtf[feature==gtf.feature, mget(gtf.orig.col)]
names(genes.gtf) = gtf.col
setkeyv(genes.gtf,gtf.col[1])
##### parse of NCBI's EntrezID to Ensembl translation
Entrez2Ensembl = fread(paste("zgrep",paste0("-E '^",taxID,"'"),gene2ENSfile))
names(Entrez2Ensembl) = gene2ENS.col
setkeyv(Entrez2Ensembl,gene2ENS.col[3])
##### add EntrezIDs to genes.gtf
tmp = Entrez2Ensembl[,mget(gene2ENS.col[2:3])]; tmp=tmp[!duplicated(tmp),]
genes.gtf = merge(genes.gtf,tmp,all.x=TRUE)
rm(tmp)
```





# Bias Reduction
uses AbundanceFunctions/processData.R
output normalized csv files in results/STAR directory
- lowess.gene.counts.csv
- lowess.SJ.counts.csv
- qspln.gene.counts.csv
- qspln.SJ.counts.csv
- quant.gene.counts.csv
- quant.SJ.counts.csv
```
LoM.norms = vector(mode='list',length=length(STAR.data$LoM.raw))
names(LoM.norms) = names(STAR.data$LoM.raw)
for(tag in names(STAR.data$LoM.raw) ){
  if( length(STAR.data$LoM.raw[[tag]])>1 ){
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[myreads]],expt.design=STAR.data$expt.design)
  } else {
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[1]],expt.design=STAR.data$expt.design)
  }
}
```



# produce QC plots
```
mytypes = names(LoM.norms)
mynorms = names(LoM.norms[[1]])
samp.labels = gsub('^[^.]+\\.','',colnames(LoM.norms[[1]][[1]]))
grpBy = annCol [[annCol.lmBy]]
clim.pct=0.96
histbins=20
for(i in 1:length(mytypes) ){
  for(j in 1:length(mynorms)){
    # set up matrices and config for plotting
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 ){
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    } else {
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    plotdata = list(plotdir='./',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    if(sum(rowmask)>hclust.limit){ #hack! build optimization!
      rowmask = rowSums(rawmat>=quantile(rawmat,probs=hc.frac.cut,na.rm=T),na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    }
    ans = summary.plots(rawmat=log2(rawmat +1), normmat=normmat, mynorm=mynorms[j], samp.labels=samp.labels, samp.classes=grpBy, plotdata=plotdata,plot2file=TRUE,histbins=histbins, colorspec=colors.rgb)
    ans = qc.clusters(rawmat=log2(rawmat[rowmask,] +1), normmat=normmat[rowmask,], attribs=annCol, oneclass=annCol.lmBy, colorspec=colors.rgb, plotdata=plotdata, plot2file=TRUE, clim.pct=clim.pct)
  }
}
```


# run regression based on Group
```
regress_lsls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(regress_lsls) = names(STAR.data$LoM.raw)
contr_ls = list(group=list(baseline="hiNa",contr.FUN="contr.treatment"))
lm_expr = "y ~ group"
rowmask_ls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(rowmask_ls) = names(STAR.data$LoM.raw)
for(i in 1:length(mytypes) ){
  for(j in 1:length(mynorms)){
  # set up data matrices
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 ){
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    } else {
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    # prepare rowmask for heatmap/MDS (remove non-expr or low expr>hclustlim)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) <= na.lim
    rowmask_ls[[ mytypes[i] ]][[ mynorms[j] ]] = rowmask # save for later
    # regression
    regress_lsls[[mytypes[i]]][[mynorms[j]]] = regressMatrix(normmat[rowmask,], expt.design=annCol, lm_expression=lm_expr, contr_list = contr_ls)
  }
}
```
# preview results

```
topn=1500
for( i in 1:length(mytypes) ){
  for( j in 1:length(mynorms) ){
    ans = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    cat("1-pi0",mytypes[i],mynorms[j],signif(unlist(lapply(2:length(ans),function(x){1-ans[[x]]$pi0})),2),"\n")
    cat("top",mytypes[i],mynorms[j],signif(unlist(lapply(2:length(ans),function(x){quantile(ans[[x]]$qvalues,probs=topn/length(ans[[x]]$qvalues))})),3),topn,"out of",length(ans[[2]]$qvalues),"\n")
}}

# screen features and plot regression QC for treatment groups
ngene_v = c(500,1500,3000)
lmBy = attr(annCol,"lmBy")
histbins=20
screen_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(screen_lsmk) = names(STAR.data$LoM.raw)
for( i in 1:length(mytypes) ){
  for( j in 1:length(mynorms) ){
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]
    mymk= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    # pull out regression design and set up plot config
    tmp = unlist(lapply(contr_ls,function(x){x$baseline}))
    tmp=paste(names(tmp),tmp,sep=".")
    plotdata = list(plotdir='./',plotbase=paste(mynorms[j],mytypes[i],tmp,sep='.'),plottitle=proj.title)
    for( ngenes in ngene_v ){ # number of top genes by q- or p-value to plot
      for(fac in names(reg_ls$q_list)[grep(lmBy,names(reg_ls$q_list))] ){
        pdata = plotdata
        pdata$plotbase = paste(plotdata$plotbase,make.names(fac),ngenes,sep='.')
        ans = qcQvalues(norm_x=normmat[mymk,], pvalue_v=reg_ls$p_mat[,fac], obj_qvalue=reg_ls$q_list[[fac]], qcut=ngenes, attribs=annCol, oneclass=lmBy, plotdata=pdata, colorspec=colors.rgb, histbins=histbins, plot2file=TRUE)
        screen_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep=":")] = ans['rowmask']
      }
    }
  }
}

# heatmaps with ratios by experimental design
# parameters
ngene_v = c(500,1500,3000)
ratioby_ls = list(group=contr_ls$group$baseline)
ratio_fold = 1.3
intensity_fold = 2
cut_ls = list(q_combine="OR", rcut_fold=ratio_fold, icut_fold=intensity_fold)
# save ratios and selections
select_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(select_lsmk) = names(STAR.data$LoM.raw)
ratio_lsmat = vector(mode='list',length=length(STAR.data$LoM.raw))
names(ratio_lsmat) = names(STAR.data$LoM.raw)
for( i in 1:length(mytypes) ){
  for( j in 1:length(mynorms) ){
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels
    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
		reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
		reg_ls = reg_ls[!grepl('Intercept',names(reg_ls))]
		plotdata = list(plotdir='./',plotbase=paste(mynorms[j],sub('counts','ratios',mytypes[i]),'minratio',paste0(ratio_fold,'x'),'minexpr',round(min(normmat,na.rm=T)+log2(intensity_fold),1),sep='.'),plottitle=proj.title)
		for( ngenes in ngene_v ) {
      # calculate and plot ratios
			cut_ls$qcut = ngenes
      pdata = plotdata
      pdata$plotbase = paste(plotdata$plotbase,'ratio.vs',ratioby_ls$strain,ngenes,sep='.')
      ans = designRatios(normmat[mymask,], attribs=annCol, ratioby_ls=ratioby_ls, plotdata=pdata, q_list=reg_ls, cut_ls=cut_ls, colorspec=colors.rgb, plot2file=TRUE)
      # save selections
      select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep=":")] = ans['rowmask']
    }
    # save ratiomat (not dependent on ngene cut) after ngene loop
    ratio_lsmat[[mytypes[i]]][[mynorms[j]]] = ans$ratiomat
	}
}
[]: 
