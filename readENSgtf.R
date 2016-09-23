# read ENSEMBL GTF format
#
# specifications: http://uswest.ensembl.org/info/website/upload/gff.html
#
# track lines:
# Although not part of the formal GFF specification, Ensembl uses track lines to further configure sets of features (thus maintaining compatibility with UCSC). Track lines should be placed at the beginning of the list of features they are to affect.  The track line consists of the word 'track' followed by space-separated key=value pairs - see the example below. Valid parameters used by Ensembl are
# *name - unique name to identify this track when parsing the file
# *description - Label to be displayed under the track in Region in Detail
# *priority - integer defining the order in which to display tracks, if multiple tracks are defined.
#
# fields:
# 1. seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
# 2. source - name of the program that generated this feature, or the data source (database or project name)
# 3. feature - feature type name, e.g. Gene, Variation, Similarity
# 4. start - Start position of the feature, with sequence numbering starting at 1.
# 5. end - End position of the feature, with sequence numbering starting at 1.
# 6. score - A floating point value.
# 7. strand - defined as + (forward) or - (reverse).
# 8. frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
# 9. attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature. [Additionally, some attributes are themselves comma-delimited key-value pairs.]

readENSgtf = function (filename, gtf.colnames=c('seqname','source','feature','start','end','score','strand','frame','attribute'), feature.col=9, comma.sub="|", curly.open.sub='<',curly.end.sub='>'){
  # arguments
  # filename = required pathname to file to be read
  # gtf.colnames = ENSEMBL GTF format field names. please check for accuracy

  # imports
  require(data.table)
  require(jsonlite)

  # int for indexing
  a = as.integer(feature.col)
  # read in data
  # default fread settings skip tracklines, identify typeof each column
  gtf = fread(filename)
  names(gtf) = gtf.colnames
  # parse attribute column: transform to JSON and use JSON parsers
  # first, protect commas or curly braces if any within fields
  mymk = grepl(',',gtf[,get(gtf.colnames[9])])
  if( sum(mymk)>0 ){
    set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub(',',comma.sub,gtf[,get(gtf.colnames[a])]))
  }
  mymk = grepl('\\{',gtf[,get(gtf.colnames[9])])
  if( sum(mymk)>0 ){
    set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub('\\{',curly.open.sub,gtf[,get(gtf.colnames[a])]))
  }
  mymk = grepl('\\}',gtf[,get(gtf.colnames[9])])
  if( sum(mymk)>0 ){
    set(x=gtf,i=as.integer(which(mymk)),j=a,value=gsub('\\}',curly.open.sub,gtf[,get(gtf.colnames[a])]))
  }
  # next, adapt ENSEMBL GTF collapsed fields to JSON format
  set(x=gtf,j=a,value=paste("{",gtf[,get(gtf.colnames[a])],"}",sep=''))
  set(x=gtf,j=a,value=gsub('; ?',',',gtf[,get(gtf.colnames[a])]))
  set(x=gtf,j=a,value=gsub('([,{}])([A-Za-z0-9_.-]+) ','\\1"\\2" : ',gtf[,get(gtf.colnames[a])]))
  set(x=gtf,j=a,value=gsub(',}','},',gtf[,get(gtf.colnames[a])]))
  # begin and end properly
  gtf[1,(gtf.colnames[a]) := gsub('^','[',get(gtf.colnames[a]))]
  gtf[nrow(gtf),(gtf.colnames[a]) := gsub(',$',']',get(gtf.colnames[a]))]
  # JSON
  gtf.attributes = as.data.table(fromJSON(gtf[,get(gtf.colnames[a])]))
  # merge
  gtf = cbind(gtf,gtf.attributes)

  return( gtf ) 
}


