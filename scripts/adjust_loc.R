#!/usr/bin/env Rscript

# Input Parametres
# ---------------- #
avinput.path=as.character(commandArgs(TRUE)[1])
loc.path=as.character(commandArgs(TRUE)[2])
# ---------------- #

df.loc = read.delim(loc.path,stringsAsFactors = F,header = F)
df.avi = read.delim(avinput.path,stringsAsFactors = F,header = F)

df.ann=NULL
for (i in 1:nrow(df.loc)) {
  if ( nchar(df.loc[i,3])>nchar(df.loc[i,4]) ) {
    alt="-"
    ref=substr(df.loc[i,3],start = 2,stop = nchar(df.loc[i,3]))
    start=df.loc[i,2]+1
    end=df.loc[i,2]+nchar(df.loc[i,3])-1
    chr=df.loc[i,1]
  } else if (nchar(df.loc[i,3]) == nchar(df.loc[i,4])) {
    alt=df.loc[i,4]
    ref=df.loc[i,3]
    start=df.loc[i,2]
    end=df.loc[i,2]
    chr=df.loc[i,1]
  } else {
    alt=substr(df.loc[i,4],start = 2,stop = nchar(df.loc[i,4]))
    ref="-"
    start=df.loc[i,2]
    end=df.loc[i,2]
    chr=df.loc[i,1]
  }
  df.ann=rbind(df.ann,data.frame("chr"=chr,"start"=start,"end"=end,"ref"=ref,"alt"=alt,stringsAsFactors = F))
}

rownames(df.avi)=paste(df.avi[,1],df.avi[,2],df.avi[,3],df.avi[,4],df.avi[,5],sep = ".")
rownames(df.loc)=paste(df.ann[,1],df.ann[,2],df.ann[,3],df.ann[,4],df.ann[,5],sep = ".")
df.loc = df.loc[rownames(df.avi),]
write.table(df.loc,file = loc.path,row.names = F,sep = "\t",col.names = F,quote = F)

