# coding=utf-8
import os
import sys
import textwrap
import argparse
import subprocess
import stat
from operator import itemgetter


def generate_meta_file(read_dir, sample_meta_data_list, select_meta_data_list, split_hyphen=None):
    """
    Generates a meta_file
    by column orientation
    """
    sample_meta_data_list = sample_meta_data_list.strip().split(',')
    select_meta_data_list = select_meta_data_list.strip().split(',')

    check = set(sample_meta_data_list).issuperset(set(select_meta_data_list))
    if check:

        files = [f for f in os.listdir(read_dir) if f.startswith('RNA')]
        files.sort()
        project = files[0].split('_')[0]
        out_f = open(project + '_metadata.txt', 'w')

        if split_hyphen:
            out_f.write('\t'.join(select_meta_data_list) + '\t' + 'ID_Group' + '\n')
        else:
            out_f.write('\t'.join(select_meta_data_list) + '\n')

        idx = [sample_meta_data_list.index(i) for i in select_meta_data_list]
        for f in files:
            metadata = f.split('_')
            pertinent = list(itemgetter(*idx)(metadata))
            if split_hyphen:
                hyphen_idx = [x for x in pertinent if '-' in x]
                if len(hyphen_idx) > 1:
                    print 'Multiple meta fields containing hyphens!'
                else:
                    pertinent.append(hyphen_idx[0].split('-')[0])
                    out_f.write('\t'.join(pertinent) + '\n')
            else:
                out_f.write('\t'.join(pertinent) + '\n')
        out_f.close()
    else:
        for i in select_meta_data_list:
            if i not in sample_meta_data_list:
                print (i)
        sys.exit(1)


def generate_slurm_submit_script(project_title, code_dir):
    log = '/'.join(code_dir.split('/')[:-2])+'/logs'
    submit = '/'.join(code_dir.split('/')[:-2])
    out_f = open(os.path.join(submit, project_title + '_analysis.submit'), 'w')
    cmd = """
#!/bin/sh
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name={project_title}
#SBATCH --error={log}/{project_title}.AutoGenerate.stdout
#SBATCH --output={log}/{project_title}.AutoGenerate.stderr
#SBATCH -p debug

Rscript {submit}/{project_title}_analysis.R
"""
    reformatted_cmd = textwrap.dedent(cmd).strip()
    context = {"project_title": project_title, "log" : log, "submit":submit}
    out_f.write(reformatted_cmd.format(**context))
    out_f.close()


def launch_slurm_submit_script(project_title, code_dir):
    log = '/'.join(code_dir.split('/')[:-2])
    sub_script = "%s_analysis.submit" % (project_title)
    analysis_script = "%s_analysis.R" % (project_title)
    scripts = [sub_script, analysis_script]
    cmd = "sbatch %s/%s_analysis.submit" % (log, project_title)
    # Set permissions on scripts to be launched to scheduler

    current_permissions = stat.S_IMODE(os.lstat(log).st_mode)
    for s in scripts:
        os.chmod(os.path.join(log, s), current_permissions)

    subprocess.call(cmd, shell=True)


def generate_abundance_script(read_dir, meta_file, code_dir, taxID, gtfFile, project_title,
                              baseline, SampleID="SampleID", mart_dataset="mmusculus_gene_ensembl",
                              lmBy="ID_Group", gtf_feature="gene", read_pattern="^RNA1", useme_cols="RNA1",
                              label_from_colname="^.*?RNA[^_]*?_([^_]+)_.*", path_type="gene.counts",
                              path_norms="loess", covariate=False):
    gtf_read_dir = '/'.join(gtfFile.split('/')[:-1])
    log = '/'.join(code_dir.split('/')[:-2])
    out_f = open(os.path.join(log, project_title + '_analysis.R'), 'w')

    code_context = {"code_dir": code_dir, "meta_file": meta_file, "SampleID": SampleID, "taxID": taxID,
                    "gtfFile": gtfFile, "gtf_feature": gtf_feature, "project_title": project_title,
                    "gtf_read_dir": gtf_read_dir, "read_dir": read_dir, "read_pattern": read_pattern,
                    "useme_cols": useme_cols, "lmBy": lmBy, "baseline": baseline, "path_type": path_type,
                    "path_norms": path_norms, "mart_dataset": mart_dataset, "label_from_colname": label_from_colname}

    if not covariate:
        contrast_str = """contr_ls = list("{lmBy}" = list(baseline="{baseline}", contr.FUN = "contr.treatment"))""".format(**code_context)
        lm_expr = "y ~ {lmBy}".format(**code_context)
        code_context['lm_expr'] = lm_expr
        code_context['contr_ls'] = contrast_str
        code_context['annColplotme'] = '"%s"' % lmBy
        code_context['annCollmBy'] = '"%s"' % lmBy
        code_context['oneclass'] = '"%s"' % lmBy
    else:
        if len(covariate.split(',')) != len(baseline.split(',')):
            print """Provide co-variates and their associated baselines in ordered comma delimited list .ie

            covariate = 'Time,Pressure'
            baseline = '0Hr,0mmHg'

            resulting in:

            contr_ls = list(Time = list(baseline="0Hr", contr.FUN = "contr.treatment"), Pressure = list(baseline="0mmHg", contr.FUN = "contr.treatment"))

            &

            lm_expr = 'y ~ Time + Pressure

            ***Note***
            Both Time and Pressure need to be column headers in the table provided as an argument to expt.design in regressMatrix

            covariate list provided:
            """
            print covariate.split(',')
            print "baseline list provided:"
            print baseline.split(',')
        else:
            code_context['oneclass'] = '"%s"' % lmBy
            covariate_list = covariate.split(',')
            baseline_list = baseline.split(',')
            lm_expr = "y ~ " +' + '.join(covariate_list)
            code_context['lm_expr'] = lm_expr
            cov_str = ""
            for i in range(len(zip(covariate_list,baseline_list))):
                cov_str += """'%s' = list(baseline='%s', contr.FUN = 'contr.treatment'),"""%(covariate_list[i],baseline_list[i])
            covariate_str = "list(" + cov_str[:-1] + ")"
            code_context['contr_ls'] = covariate_str
            reformat_covariate = '"'+'","'.join(covariate_list)+'"'
            code_context['annColplotme'] = 'c(%s)'%reformat_covariate
            code_context['annCollmBy'] = 'c(%s)'%reformat_covariate


    code = """
require(data.table)
require(NMF)
require(affy)
require(limma)
#require(AnnotationDbi)
require(biomaRt)
source("{code_dir}AbundanceFunctions/BiasReduce.R")
source("{code_dir}AbundanceFunctions/ExtractTransformLoad.R")
source("{code_dir}AbundanceFunctions/DifferentialAnalysis.R")
source("{code_dir}AbundanceFunctions/NonVisualOutput.R")
source("{code_dir}GenomicsFunctions/ReadAndParse.R")
source("{code_dir}AssociationFunctions/gs.wrapper.R")
source("{code_dir}AssociationFunctions/PathwayAnalysis.R")
source("{code_dir}BcorePlotting/SummaryPlots.R")
source("{code_dir}BcorePlotting/MultipleTestingCorrection.R")
source("{code_dir}BcorePlotting/ClusteringPlots.R")


setwd("{read_dir}")

dir.create(file.path(getwd(),'{project_title}'),showWarnings=FALSE)
oneclass = {oneclass}
# constants
# max data rows for hclust
hclust.limit = 2^16
annCol.lmBy = {annCollmBy}
# quantile of data distribution reqd in one group's worth of data if too many rows for hclust()
hc.frac.cut = 0.75;
SJ.counts.na.frac = 0.25;
# max fraction of samples not having detected a splice junction for the splice
# junction to be retained in raw data


# regression parameters
na.lim = 0 # max NAs per row tolerated by lm() at least in some cases
do.not.regress = "alograw" # control norm not to be used for regression stats
# plotting colors
colors.rgb = c(rgb(0,0,0),rgb(0.1,0.1,1),rgb(0,.7,.7),rgb(0,.7,0),rgb(.7,1,0),rgb(.7,0,.7))
md.file = "{meta_file}"
md.orientation = "byRow" # sampleIDs are in @ row. alt:byCol (IDs in @ col)
md.IDcol = "{SampleID}" # reqd if md.orientation is byRow; byCol==headers are IDs

# gene annotation
taxID = {taxID}
gene2ENSfile = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/ncbi/gene2ensembl.gz"
gene2ENS.col = c("taxID","EntrezID","Gene","RefSeqTranscript","EnsemblTranscript","RefSeqProtein","EnsemblProtein")
gtfFile = "{gtfFile}"
gtf.feature = "{gtf_feature}"
gtf.orig.col = c("gene_id","gene_name","gene_biotype")
gtf.col = c("Gene","Symbol","biotype")
proj.title = "{project_title}"
gtf.Rdir = "{gtf_read_dir}"
genome.func = "{code_dir}GenomicsFunctions"


readdir = "{read_dir}"
readpattern = "{read_pattern}"
useme.cols = "{useme_cols}"
label.from.colname = "{label_from_colname}"
samps = dir(path=readdir, pattern=readpattern)
samp.labels = gsub(label.from.colname,'\\\\1', samps)

annCol.names = "group"
annCol.label = "{label_from_colname}"

# read in STAR alignments

print("Read in Abundance data")

STAR.data = read_STAR(useme.cols=useme.cols,label.from.colname=label.from.colname,annCol.label=annCol.label,annCol.names=annCol.names,annCol.normBy=NULL,annCol.lmBy=annCol.lmBy, readpattern=readpattern, unstranded.col=list(gene.counts=c(1:4), SJ.counts=c(1:3,7)))

print("Read in Abundance data : Complete")

print("Filter SJ.counts data")
# filter out SJs with too many NAs
{{
    if(any( names(STAR.data$LoM.raw)=="SJ.counts" & exists("SJ.counts.na.frac")))
    {{
      STAR.data$SJ.counts.orig = STAR.data$LoM.raw$SJ.counts
      for(tag in names(STAR.data$LoM.raw$SJ.counts))
      {{
        STAR.data$LoM.raw$SJ.counts[[tag]] = STAR.data$LoM.raw$SJ.counts[[tag]][rowSums(is.na(STAR.data$LoM.raw$SJ.counts[[tag]]))<= (SJ.counts.na.frac*ncol(STAR.data$LoM.raw$SJ.counts[[tag]])), ]
      }}
    }}
}}
print("Filtering SJ.counts : Complete")
save.image("./{project_title}.RData")

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

myreads = STAR.data$myreads

print("Filter gene.counts data")

{{
    if(any( names(STAR.data$LoM.raw)=="gene.counts"))
    {{
      STAR.data$LoM.orig.raw = STAR.data$LoM.raw$gene.counts
      for(tag in names(STAR.data$LoM.raw$gene.counts))
      {{
        STAR.data$LoM.raw$gene.counts[[tag]] = STAR.data$LoM.raw$gene.counts[[tag]][rowSums(STAR.data$LoM.raw$gene.counts[[tag]])>0, ]
      }}
    }}
}}

print("Filtering gene.counts : Complete")

save.image("./{project_title}.RData")

LoM.norms = vector(mode='list',length=length(STAR.data$LoM.raw))
names(LoM.norms) = names(STAR.data$LoM.raw)


print("Bias reduction with normMatrix")

for(tag in names(STAR.data$LoM.raw))
{{
  if( length(STAR.data$LoM.raw[[tag]])>1 )
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[myreads]])
  }}
  else
  {{
    LoM.norms[[tag]] = normMatrix(tag=tag,raw.mat=STAR.data$LoM.raw[[tag]][[1]])
  }}
}}

print("Bia reduction with normMatrix : Complete")

save.image("./{project_title}.RData")

## custom: reannotate from metadata
## run if necessary metadata is in file and not also in FASTQ file names
# map annotation to read matrix
md.dt = fread("{meta_file}")

print("Evaluate meta-data file for sample name consistency")

{{
    if( md.orientation == "byCol" ){{ # samples are one per column
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), names(md.dt) )
    }}
    else
    {{ # samples are one per row
      idx2 = match( colnames(STAR.data$LoM.raw[[1]][[myreads]]), md.dt[,get(md.IDcol)] )
    }}
}}
idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
{{
    if( sum(idx1==1:ncol(STAR.data$LoM.raw[[1]][[myreads]]))==ncol(STAR.data$LoM.raw[[1]][[myreads]]))
    {{
      if( md.orientation == "byCol" )
      {{ # samples are one per column
        annCol = NULL
        namecol = setdiff( 1:ncol(md.dt), idx2 )
        md.factors = as.character(md.dt[,namecol,with=F])
        for(k in 1:nrow(md.dt) )
        {{
          annCol[[ md.dt[k,get(names(md.dt)[namecol])] ]][idx1] = as.vector( md.dt[ k, mget(names(md.dt)[idx2]) ] )
        }}
      }}
      else
      {{ # samples are one per row
        annCol = NULL
        namecol = setdiff(names(md.dt), md.IDcol)
        md.factors = as.character(namecol)
        for( k in setdiff(names(md.dt),md.IDcol))
        {{
          annCol[[ k ]][idx1] = as.vector( md.dt[ idx2, get(k) ] )
        }}
      }}
    }}
    else
    {{
      stop(paste("Some samples have no annotation in",md.file))
    }}
}}

print("Evaluation of meta-data : Complete")

mytypes = names(LoM.norms)
mynorms = names(LoM.norms[[1]])

#grpBy = annCol[[annCol.lmBy]]
grpBy = annCol[[oneclass]]
annCol.plotme = {annColplotme}
clim.pct=0.96
histbins=20
dir.create(file.path(getwd(),'{project_title}/summary.plotsPlots'),showWarnings=FALSE)
dir.create(file.path(getwd(),'{project_title}/qc.clustersPlots'),showWarnings=FALSE)

print("Generate summary and qc.cluster plots")

for(i in 1:length(mytypes) )
{{
  for(j in 1:length(mynorms))
  {{
    # set up matrices and config for plotting
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}

    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels
    plotdata = list(plotdir='./{project_title}/summary.plotsPlots',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)))

    if(sum(rowmask)>hclust.limit)
    {{
      rowmask = rowSums(rawmat>=quantile(rawmat,probs=hc.frac.cut,na.rm=T),na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) < (ncol(normmat)/length(unique(grpBy)) )
    }}

    ans = summary.plots(rawmat=log2(rawmat +1), normmat=normmat, mynorm=mynorms[j], samp.labels=samp.labels, samp.classes=grpBy, plotdata=plotdata,plot2file=TRUE,histbins=histbins, colorspec=colors.rgb)
    plotdata = list(plotdir='./{project_title}/qc.clustersPlots',plotbase=paste(mynorms[j],mytypes[i],sep='.'),plottitle=proj.title)
    ans = qc.clusters(rawmat=log2(rawmat[rowmask,] +1), normmat=normmat[rowmask,], attribs=annCol[annCol.plotme], oneclass=oneclass, colorspec=colors.rgb, plotdata=plotdata, plot2file=TRUE, clim.pct=clim.pct)
  }}
}}

print("Generate summary and qc.cluster plots : Complete")

save.image("./{project_title}.RData")

regress_lsls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(regress_lsls) = names(STAR.data$LoM.raw)

#contr_ls = list("{lmBy}"=list(baseline="{baseline}",contr.FUN="contr.treatment"))
contr_ls = {contr_ls}

# set baseline for regression in parameter(s) of interest
# contr.treatment generates regression coefficients that are like (adjusted mean) ratios of other groups to the baseline group.
# contr.sum generates coefficients that are like (adjusted mean) ratios to average all for all but the mandadory ommitted treatment group (because there is always one fewer independent pairwise comparison than there are pairs).


lm_expr = "{lm_expr}"

rowmask_ls = vector(mode='list',length=length(STAR.data$LoM.raw))
names(rowmask_ls) = names(STAR.data$LoM.raw)
dir.create(file.path(getwd(),'{project_title}/regressMatrixPlots'),showWarnings=FALSE)

print("RegressMatrix with specified contrasts")

for(i in 1:length(mytypes))
{{
  for(j in which(!mynorms %in% do.not.regress))
  {{
    # set up data matrices
    if( length(STAR.data$LoM.raw[[ mytypes[i] ]])>1 )
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[myreads]]
    }}
    else
    {{
      rawmat = STAR.data$LoM.raw[[ mytypes[i] ]][[1]]
    }}
    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./{project_title}/regressMatrixPlots/',plotbase=paste(mynorms[j],mytypes[i],'vs',tmp,sep='.'),plottitle=proj.title)
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(rawmat)=samp.labels; colnames(normmat)=samp.labels

    # prepare rowmask for heatmap/MDS (remove non-expr or low expr>hclustlim)
    rowmask = rowSums(rawmat>1,na.rm=T) > (ncol(rawmat)/length(unique(grpBy)) ) & rowSums(is.na(normmat)) <= na.lim
    rowmask_ls[[ mytypes[i] ]][[ mynorms[j] ]] = rowmask # save for later

    # regression

    regress_lsls[[mytypes[i]]][[mynorms[j]]] = regressMatrix(normmat[rowmask,], expt.design=annCol[annCol.lmBy], lm_expression=lm_expr, contr_list = contr_ls, plot2file = TRUE, plotdata = plotdata)
  }}
}}

print("RegressMatrix with specified contrasts : Complete")

save.image("./{project_title}.RData")

topn=500 # number with which to play

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress) )
  {{
    ans = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    cat(mytypes[i],mynorms[j],"1-p0:",signif(unlist(lapply(2:length(ans),function(x){{1-ans[[x]]$pi0}})),2),"\n")
    cat(mytypes[i],mynorms[j],"qcut:",signif(unlist(lapply(2:length(ans),function(x){{quantile(ans[[x]]$qvalues,probs=topn/length(ans[[x]]$qvalues))}})),3),
    "returns ~",topn,"out of",length(ans[[2]]$qvalues),"\n")
  }}
}}

# Plots for each design factor (default) or for factors specified in facSel
#   1) Histogram of p values that were included in the design
#      Look for a peak on the left, and no peaks in the middle or on the right
#   2) qvalue's default plots, with full qvalue range c(0,1) plotted
#      Look for descent to pi0 in the top left tuning plot with good asymptote
#      The slow/steep rise in q-values in remaining plots depends on resolving power of data


lmBy = annCol.lmBy
histbins=20
select_lsmk = vector(mode='list',length=length(STAR.data$LoM.raw))
names(select_lsmk) = names(STAR.data$LoM.raw)
dir.create(file.path(getwd(),'{project_title}/qcQvaluePlots'),showWarnings=FALSE)

print("Perform qcQvalue & Generate masks")

for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]
    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    # pull out regression design and set up plot config

    tmp = unlist(lapply(contr_ls,function(x){{x$baseline}}))
    tmp=paste(names(tmp),tmp,sep=".")

    plotdata = list(plotdir='./{project_title}/qcQvaluePlots',plotbase=paste(mynorms[j],mytypes[i],tmp,sep='.'),plottitle=proj.title)

    ##for(fac in names(reg_ls$q_list)[grep(lmBy,names(reg_ls$q_list))])

    for(fac in names(reg_ls$q_list))
    {{
      if(!fac %in% '(Intercept)')
      {{
        pdata = plotdata
        pdata$plotbase = paste(plotdata$plotbase,make.names(fac),sep='.')
        ans = qcQvalues(norm_x=normmat[mymask,], pvalue_v=reg_ls$p_mat[,fac], obj_qvalue=reg_ls$q_list[[fac]], attribs=annCol[lmBy],
        oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, histbins=histbins, plot2file=TRUE)
        select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac)] = ans['rowmask']
      }}
    }}
  }}
}}

print("Perform qcQvalue & Generate masks : Complete")

save.image("./{project_title}.RData")

# Build of ratios to baselines for desired design factors
# Also build masks selecting genes based on q-values per factor
# and follow-up considerations such as expression level and fold change
# Loop over q-value cuts to assist with final cut selection.
# Plot heatmaps and MDS plots based on selected genes and designed ratios


ngene_v = c(200,500,1000) # q-value cuts by number; can also cut by q-value
ratioby_ls = list("{lmBy}"=contr_ls${lmBy}$baseline)
ratio_fold = 1.3
intensity_fold = 2

cut_ls = list(q_combine="OR", rcut_fold=ratio_fold, icut_fold=intensity_fold)

# settings for heatmaps
#annCol.plotme = annCol.lmBy # heatmap tracks

clustrowmin = 10 # min data rows for heatmap and MDS plots

# save ratios and selections

ratio_lsmat = vector(mode='list',length=length(STAR.data$LoM.raw))
names(ratio_lsmat) = names(STAR.data$LoM.raw)

print("Perform PlotRatios")

dir.create(file.path(getwd(),'{project_title}/plotRatiosPlots'),showWarnings=FALSE)
for( i in 1:length(mytypes))
{{
  for( j in which(!mynorms %in% do.not.regress))
  {{
    normmat= LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    colnames(normmat)=samp.labels

    if( grepl('SJ',mytypes[i]))
    {{#not the most robust way to find SJs
      # map SJ positions to overlapping genes
      # specifying defaults to fix "Gene" and Pos as colnames
      gtf.col = c("Gene","Symbol","Chr","start","stop")
      pos.col = "Pos"
      mySJ_dt = mapSJ2feature(STAR.data$LoM.raw[[ mytypes[i] ]][[1]], pos.col=pos.col, gtf.col=gtf.col, gtf.file=gtfFile, gtf.Rdir=genome.func)
      # min data rows for heatmap and MDS plots
      # map gene IDs back to data matrix using mapping built earlier
      # this enables include_ID to select rows

      idx2 = match( rownames(STAR.data$LoM.raw[[ mytypes[i] ]][[1]]), mySJ_dt[,get(pos.col)] )
      idx1 = which(!is.na(idx2)); idx2 = idx2[idx1]
      rownames(normmat)[idx1] = mySJ_dt[idx2,get(gtf.col[1])]
    }}

    mymask= rowmask_ls[[mytypes[i]]][[mynorms[j]]]
    reg_ls = regress_lsls[[mytypes[i]]][[mynorms[j]]]$q_list
    reg_ls = reg_ls[!grepl('Intercept',names(reg_ls))]
    plotdata = list(plotdir='./{project_title}/plotRatiosPlots',plotbase=paste(mynorms[j],sub('counts','ratios',mytypes[i]),'minratio',paste0(ratio_fold,'x'),'minexpr',round(min(normmat,na.rm=T)+log2(intensity_fold),1),sep='.'),plottitle=proj.title)

    #for(fac in names(reg_ls)[grep(lmBy,names(reg_ls))])
    for(fac in names(reg_ls))
    {{
        for( ngenes in ngene_v )
        {{
          # calculate and plot ratios
          cut_ls$qcut = ngenes
          pdata = plotdata
          pdata$plotbase = paste(plotdata$plotbase,'ratio.vs',ratioby_ls${lmBy},ngenes,sep='.')

          # this function makes the ratios and cuts
          ans = designRatios(normmat[mymask,], q_list=reg_ls[[fac]], attribs=annCol, ratioby_ls=ratioby_ls, cut_ls=cut_ls)

          # this function creates the heatmap and MDS plot
          if( sum(ans$rowmask)>clustrowmin )
          {{
            ans2= plotRatios( ratiomat=ans$ratiomat, attribs=annCol[annCol.plotme], oneclass=oneclass, plotdata=pdata, colorspec=colors.rgb, rowmask=ans$rowmask, plot2file=TRUE)
          }}
          # save selections
          select_lsmk[[mytypes[i]]][[mynorms[j]]][paste(fac,ngenes,sep="_")] = ans['rowmask']
        }}
        # save ratiomat (not dependent on ngene cut) after ngene loop
        ratio_lsmat[[mytypes[i]]][[mynorms[j]]] = ans$ratiomat
    }}
  }}
}}

print("Perform PlotRatios : Complete")

save.image("./{project_title}.RData")

dir.create(file.path(getwd(),'{project_title}/PairwiseScatter'),showWarnings=FALSE)

print("Perform Pairwise Scatter")

for( i in 1:length(mytypes))
{{
  for( j in 1:length(mynorms))
  {{
    normmat = LoM.norms[[ mytypes[i] ]][[ mynorms[j] ]]
    attribs = annCol[[oneclass]]
    plotdir = paste(getwd()[1], "/{project_title}/PairwiseScatter/",sep='')
    plotdata = list(plotdir=plotdir,plotbase=paste(mynorms[j],mytypes[i],'Scatter',lmBy,sep=':'),plottitle=proj.title)
    scatterplot(normmat, attribs=attribs, plotdata=plotdata, plot2file=TRUE)
  }}
}}

print("Perform Pairwise Scatter : Complete")

save.image("./{project_title}.RData")

print("Perform biomaRt gene annotation")

ratio_mat = ratio_lsmat[["{path_type}"]][["{path_norms}"]]
reg_ls = regress_lsls[["{path_type}"]][["{path_norms}"]]
ID = rownames(ratio_mat)


project_species = useMart("ENSEMBL_MART_ENSEMBL", dataset = "{mart_dataset}", host = "dec2016.archive.ensembl.org")
human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "dec2016.archive.ensembl.org")

print("biomaRt gene annotations : Loaded")

hsa_entrezID = getLDS(attributes = "ensembl_gene_id",mart = project_species,  attributesL = "entrezgene", martL = human)
Entrez.ID = character(length(ID))
Entrez.ID[] = NA
idx2 = match(ID,hsa_entrezID$Gene.ID)
idx1 = which(!is.na(idx2)); idx2= idx2[idx1]
Entrez.ID[idx1] = hsa_entrezID$EntrezGene.ID[idx2];
Entrez.ID = as.numeric(Entrez.ID)

###Human gene symbols

hsa_symbol = getLDS(attributes = "ensembl_gene_id",mart = project_species,  attributesL = "hgnc_symbol", martL = human)
Anno.Symbol = character(length(ID))
Anno.Symbol = NA
idx2 = match(ID,hsa_symbol$Gene.ID)
idx1 = which(!is.na(idx2)); idx2= idx2[idx1]
Anno.Symbol[idx1] = unlist(hsa_symbol$HGNC.symbol[idx2]);
# assemble IDs for annotation

print("Species event ID's mapped to Hsapiens_gene_ensembl")


backgroundset = as.data.table(cbind(ID, Entrez.ID, Anno.Symbol))

# assemble signatures

reg_ls = regress_lsls$gene.counts${path_norms}

masks = names(select_lsmk[["{path_type}"]][["{path_norms}"]])
dir.create(file.path(getwd(),'{project_title}/tables'),showWarnings=FALSE)

print("Perform pathway enrichment for all factor levels relative to baseline")

for (k in masks)
{{
    sig_gmt = NULL
    rowmask = select_lsmk[["{path_type}"]][["{path_norms}"]][[k]]
    sig_gmt[['all']] =  ID[rowmask]
    fileSettings = list(directory=file.path(getwd(),'{project_title}/tables'),baseFilename=paste(k,'with_Uni',sep='_'))
    path_ls = ags.wrapper(setlist=sig_gmt, backgroundset, include.identifiers=TRUE, anno.uni=TRUE,
        fileSettings = fileSettings,
        functiondir = "{code_dir}AssociationFunctions/",
        resourcedir = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/old_association_resources",
        return.OM=TRUE, ecut=0.05, ocut=5)

    fileSettings = list(directory=file.path(getwd(),'{project_title}/tables'),baseFilename=paste(k,'with_out_Uni',sep='_'))
    path_ls = ags.wrapper(setlist=sig_gmt, backgroundset, include.identifiers=TRUE, anno.uni=FALSE,
        fileSettings = fileSettings,
        functiondir = "{code_dir}AssociationFunctions/",
        resourcedir = "/home/exacloud/lustre1/BioCoders/DataResources/AnnotationSources/old_association_resources",
        return.OM=TRUE, ecut=0.05, ocut=5)
}}

print("Perform pathway enrichment for all factor levels relative to baseline : Complete")

print("Generate Abundance and Ratio table with associated q-value and p-values")

gtf.Rdir = "{gtf_read_dir}"

out_norm_mat = LoM.norms$gene.counts$loess[rowmask_ls${path_type}${path_norms},]
out_table = outputTable(normmat= out_norm_mat, gtf.file = "{gtfFile}",ratiomat = ratio_lsmat${path_type}${path_norms}, q_list=reg_ls$q_list, gtf.Rdir=genome.func, gtf.key='transcript')
write.csv(out_table, row.names = FALSE, file=file.path(getwd(), './{project_title}/tables', paste("{project_title}","{path_norms}","Normed_with_Ratio_and_Abundance.csv", sep="_")),quote=FALSE)

print("Generate Abundance and Ratio table with associated q-value and p-values : Complete")

save.image("./{project_title}.RData")
"""
    reformatted_code = textwrap.dedent(code).strip()
    out_f.write(reformatted_code.format(**code_context))
    out_f.close()


def main():
    parser = argparse.ArgumentParser(description='Generate Abundance Workflow Wrapper')
    parser.add_argument("-d", "--read_dir", type=str, help="Absolute path to abundance data directory i.e ../STAR/")
    metagroup = parser.add_argument_group('Mandatory Metadata arguments')
    metagroup.add_argument("-mb", "--meta", action="store_true", dest="meta", help="Boolean to determine whether to make a meta file", default=False)
    metagroup.add_argument("-md", "--sample_meta_data_list", type=str, help="comma delimited list of meta information provided by underscore delimited sample directory name i.e RNA160606DM_294-1_S2_L001_R1_001 = 'Project','ID','sample','lane','r1','01'")
    metagroup.add_argument("-ms", "--select_meta_data_list", type=str, help="subset of the comma delimited list of meta information provided by underscore delimited sample directory name i.e 'ID','sample': These fields will be incorporated into metadata.txt")
    metagroup.add_argument("-sh", "--split_hyphen", action="store_true", dest="split_hyphen", help="split and incorporate hyphenated meta data field into meta data i.e. 294-1 in RNA160606DM_294-1_S2_L001_R1_001 will be incorporated into metadata.txt under the column ID_Group as 294", default=False)

    abundancegroup = parser.add_argument_group('Mandatory Abundance script generation arguments')
    abundancegroup.add_argument("-lj", "--launch_job", action="store_true", dest="launch_job", help="Generate SLURM submission script and launch job", default=False)
    abundancegroup.add_argument("-mf", "--meta_file", type=str, help="Absolute path to metafile.txt generated via .generate_meta_file")
    abundancegroup.add_argument("-c", "--code_dir", type=str, help="Absolute path to code directory i.e. ProjectDirectory/code/")
    abundancegroup.add_argument("-t", "--taxID", type=str, help="TaxID can be found www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi")
    abundancegroup.add_argument("-g", "--gtfFile", type=str, help="Absolute path to gtf used for alignment")
    abundancegroup.add_argument("-p", "--project_title", type=str, help="Project title associated with abundance dataset. (Will be incorporated into base file name of plots)")
    abundancegroup.add_argument("-b", "--baseline", type=str, help="Baseline to generate contrasts against when generating lm",default=False)

    optabundancegroup = parser.add_argument_group('Optional Abundance script generation arguments')
    optabundancegroup.add_argument("-co", "--covariate", type=str, help="Covariates to generate contrasts against when generating lm")
    optabundancegroup.add_argument("-id", "--SampleID", type=str, help="Column in metadata.txt that identifies samples read in by label_from_colname", default="SampleID")
    optabundancegroup.add_argument("-bm", "--mart_dataset", type=str, help="biomaRt dataset to use. Datasets include: mmusculus_gene_ensembl | hsapiens_gene_ensembl | aplatyrhynchos_gene_ensembl | drerio_gene_ensembl | ggallus_gene_ensembl | oaries_gene_ensembl | rnorvegicus_gene_ensembl | sscrofa_gene_ensembl", default="mmusculus_gene_ensembl")
    optabundancegroup.add_argument("-lm", "--lmBy", type=str, help="Column name in metadata.txt that contains [baseline] value", default="ID_Group")
    optabundancegroup.add_argument("-gf", "--gtf_feature", type=str, help="GTF feature to annotate abundance and ratio tables with", default="gene")
    optabundancegroup.add_argument("-rp", "--read_pattern", type=str, help="Read pattern expression provided to R to read in sample associated abundance information", default="^RNA1")
    optabundancegroup.add_argument("-uc", "--useme_cols", type=str, help="Read pattern expression provided to R to select data to be incorporated in STAR.data ", default="RNA1")
    optabundancegroup.add_argument("-lc", "--label_from_colname", type=str, help="Read pattern expression provided to R to select for unique sample label identifiers", default="^.*?RNA[^_]*?_([^_]+)_.*")
    optabundancegroup.add_argument("-pt", "--path_type", type=str, help="gene.counts | SJ.counts", default="gene.counts")
    optabundancegroup.add_argument("-pn", "--path_norms", type=str, help=" alograw | loess | lowess | qspln | quant", default="loess")

    args = parser.parse_args()
    if args.meta:
        generate_meta_file(args.read_dir, args.sample_meta_data_list, args.select_meta_data_list, args.split_hyphen)
    else:
        generate_abundance_script(args.read_dir, args.meta_file, args.code_dir, args.taxID, args.gtfFile,
                                  args.project_title,
                                  args.baseline, args.SampleID, args.mart_dataset,
                                  args.lmBy, args.gtf_feature, args.read_pattern, args.useme_cols,
                                  args.label_from_colname, args.path_type,
                                  args.path_norms, args.covariate)
        if args.launch_job:
            generate_slurm_submit_script(args.project_title, args.code_dir)
            launch_slurm_submit_script(args.project_title, args.code_dir)

if __name__ == '__main__':
    main()

