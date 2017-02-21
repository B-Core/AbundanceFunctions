################################################################################
# Functions involved in generating tabular output or other non-plot output
#
# outputTable     --- a data.table for reporting, with abundance data and identifiers
#                     and optionally also with ratio data and regression results
#
# Authors: Theresa Lusardi, Julja Burchard, and Mark Fisher
# Started - July 2016
################################################################################


outputTable <-
function(normmat, gtf.file, rowmask=NULL, ratiomat=NULL, 
  q_list=NULL, p_flag=TRUE, b_x=NULL,  is_a ="gene", 
  gtf.key = "gene", gtf.feature = "feature", 
  gtf.orig.col = c("gene_id","gene_name","seqname","start","end"), 
  gtf.col = c("Gene","Symbol","Chr","start","stop"), 
  ikey = 1, isym = 2, ic = 3, i0 = 4, i1 = 5, #chr, start, stop
  pos.col = "Pos", pos.delim ="\\.", delim.sub="_",
  gtf.Rfile="ReadAndParse.R", gtf.Rdir="GenomicFunctions"
  ){
  # return a data.table for reporting, with abundance data and identifiers
  #     and optionally also with ratio data and regression results:
  #       q-values, optionally p-values, optionally coefficients
  # normmat: matrix of abundance data with identifiers as rownames
  # ratiomat: matrix of ratio data with identifiers as rownames
  # rowmask: logical vector used to shrink normmat to q_list nrow, if applicable
  # q_list: either an object returned by qvalue(), or a list of same
  # p_flag: if TRUE, include p-values as well as q-values in the output
  # b_x: vector or matrix of regression coefficients to include in the output
  #      assumed to be the same nrow as q_list, optimally from same regression
  # is_a : "SJ" for splice junction normmat, otherwise same as gtf.key
  # gtf.file : file of genome annotation in gtf or gff format
  #           usually best to parse GTF used for data processing
  # gtf.feature : column with biotype AND rowtype information
  # gtf.key : key rowtype to keep, given as value found in feature column
  # gtf.orig.col : column names to keep in extracted gtf.file rows 
  # ikey, isym, ic, i0, i1 : short indices for key feature, gene symbol and 
  #                          SJ position columns chr, start, stop
  # pos.col : column for SJ identifiers comprising chr, start, stop
  # pos.delim : delimiter for splitting SJ identifiers into chr, start, stop
  # delim.sub : replacement for delimiter if found chr field
  # gtf.Rfile : R file with readENSgtf function
  # gtf.Rdir : /Path/To/Directory containing R file with readENSgtf function
  #            trailing slash is currently NOT enabled

   
  # imports
  require(data.table)
  source(paste(gtf.Rdir,gtf.Rfile,sep="/")) # very trusting!

  # check required arguments
  if( is.null(dim(normmat)) ){
    stop("Argument normmat is not a matrix")
  }
  if( !file.exists(gtf.file) ){
    stop("GTF file does not exist: ",gtf.file)
  }
  # check optional arguments
  if( !is.null(rowmask) ){
    if(length(rowmask) != nrow(normmat) ){
      stop("rowmask length must match nrow(normmat)")
    }
    rdx = as.integer(which(rowmask)) # save positions for later
  } else {
    rdx = as.integer(1:nrow(normmat))# save positions anyway
  }

  # convert normmat to data.table
  out_dt = data.table(normmat, keep.rownames=T)

  # add ratio data if given
  if( !is.null(ratiomat) ){
    if( is.null(dim(ratiomat)) ){
      warning("Argument ratiomat is not a matrix. Skipping...")
    } else if( nrow(ratiomat)==nrow(normmat) & 
              sum(rownames(ratiomat)==rownames(normmat),na.rm=T) ==
                  nrow(normmat) ){
      colnames(ratiomat)=paste0(colnames(ratiomat),".ratio")
      out_dt = cbind(out_dt,ratiomat)
    } else if( !is.null(rowmask) & 
              nrow(ratiomat) == sum(rowmask) &
              sum(rownames(ratiomat)==rownames(normmat)[rdx],na.rm=T) ==
                  nrow(ratiomat) ){
      set( out_dt, j=paste0(colnames(ratiomat),".ratio"), value=as.double(NA) )
      set( out_dt, i=rdx, j=paste0(colnames(ratiomat),".ratio"),
          value=data.table(ratiomat) )
    }
  }

  # add q_list data if given
  if( !is.null(q_list) ){
    if( any(grepl("qvalues",names(q_list) )) ){ # one qvalue object supplied
      # rearrange for easier looping
      tmp = q_list; q_list=NULL; q_list[[1]] = tmp
      names(q_list)[1] = "Q"
    }
    for( i in 1:length(q_list) ){
      qname = make.names(names(q_list)[i])
      if( !any(grepl("qvalues",names(q_list[[i]]) )) ){
        warning("No qvalues in q_list element ",i)
      } else if((!is.null(rowmask) & sum(rowmask)==length(q_list[[i]]$qvalues))|
              (is.null(rowmask) & nrow(normmat)==length(q_list[[i]]$qvalues)) ){
        if( p_flag ){
          set(out_dt, i=rdx, j=paste0(qname,".pval"), value=q_list[[i]]$pvalues)
        }
        set(out_dt, i=rdx, j=paste0(qname,".qval"), value=q_list[[i]]$qvalues)
      } else {
        warning("Dimensions of q_list do not match normmat. Skipping...")
      }
    }
  }

  # add coefficient (beta) data if given
  if( !is.null(b_x) ){
    if( is.null(dim(b_x)) ){ # vector data
      if((!is.null(rowmask) & sum(rowmask)==length(b_x))|
         (is.null(rowmask) & nrow(normmat)==length(b_x)) ){
        set(out_dt, i=rdx, j="beta", value=b_x)
      }
    } else {                 # matrix data
      if( nrow(b_x)==nrow(normmat) &
         sum(rownames(b_x)==rownames(normmat),na.rm=T) ==
                  nrow(normmat) ){
        colnames(b_x)=paste0(colnames(b_x),".beta")
        out_dt = cbind(out_dt,b_x)
      } else if( nrow(b_x)==sum(rowmask) &
                !is.null(rowmask) & 
              sum(rownames(b_x)==rownames(normmat)[rdx],na.rm=T) ==
                  nrow(b_x) ){
        set( out_dt, j=paste0(colnames(b_x),".beta"), value=as.double(NA) )
        set( out_dt, i=rdx, j=paste0(colnames(b_x),".beta"),
          value=data.table(b_x) )
      }
    }
  }

  # map data to identifiers
  if( any(grepl('SJ|splice',is_a,ignore.case=T)) ){
    map_dt = mapSJ2feature(normmat=normmat, gtf.file=gtf.file, 
              gtf.key=gtf.key, gtf.orig.col=gtf.orig.col, gtf.col=gtf.col, 
              ikey=ikey, isym=isym, ic=ic, i0=i0, i1=i1, pos.col=pos.col, 
              pos.delim=pos.delim, delim.sub=delim.sub, source.me=source.me )
    # merge with out_dt
    setkeyv(out_dt,"rn")
    out_dt = merge(map_dt[,mget(setdiff(names(map_dt),colnames(normmat)))],
                   out_dt, by.x=pos.col, by.y="rn")
  } else {
    # read in genome annotation
    gtf = readENSgtf(filename=gtf.file)
    message("Genome annotation file ",gtf.file," read in with ",nrow(gtf)," rows")
    feature.gtf = gtf[get(gtf.feature)==gtf.key, mget(gtf.orig.col)]
    names(feature.gtf) = gtf.col
    setorderv(feature.gtf,cols=gtf.col[c(ic,i0,i1)],na.last=TRUE)
    # note setkeyv MUST be run 2nd as setorderv wipes key
    setkeyv(feature.gtf,gtf.col[ikey]) 
    # merge
    setkeyv(out_dt,"rn") # data.table names rowname column "rn"
    out_dt = merge(feature.gtf, out_dt, by.x=gtf.col[ikey], by.y="rn", all.y=T)
  }

  # return data.table
  return(out_dt)

}
