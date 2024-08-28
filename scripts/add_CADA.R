#!/usr/bin/env Rscript

CADA_dict.path = "~/opt/rare/Rdata/CADA_dict.rds"
df.dict = readRDS(CADA_dict.path)

rare_res.path = as.character(commandArgs(TRUE)[1])
cada.res.path = as.character(commandArgs(TRUE)[2])

df.rare = read.delim(rare_res.path,stringsAsFactors = F)
df.cada = read.delim(cada.res.path,stringsAsFactors = F)

# Assign min score to negative values and normalize between 0 and 1
df.cada$score[df.cada$score<0] = min(df.cada$score[df.cada$score>0])
df.cada$score = df.cada$score/100

# assign scores
df.dict$score = df.cada$score[match(df.dict$entrezID,df.cada$gene_id)]
df.dict = df.dict[order(df.dict$score,decreasing = T),]

# merge CADA res
df.rare$CADA.score = df.dict$score[match(df.rare$Symbol,df.dict$symbol)]
df.rare$RARE.score = apply(cbind(df.rare$Maverick.Score,df.rare$CADA.score),MARGIN = 1,function(x) mean(x,na.rm = T))
df.rare = df.rare[order(df.rare$RARE.score,decreasing = T),]

cols = c("CHROM","POS","REF","ALT","QUAL","VarType","Symbol","AD","DP","GQ","GT","VAF","Maverick.Score","CADA.score","RARE.score","AlphaMissense","ESM1b","EVE","REVEL","Orphanet.id","Orphanet.Phenotype","OMIM.id","OMIM.title","OMIM.Phenotype")
df.rare = cbind(df.rare[,cols],df.rare[,!colnames(df.rare)%in%cols])

write.table(df.rare,file = rare_res.path,quote = F,sep = "\t",row.names = F)

