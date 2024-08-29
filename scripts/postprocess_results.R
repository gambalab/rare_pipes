#!/usr/bin/env Rscript

orphanet.path = "~/opt/rare/Rdata/orphanet.rds"
omim.path = "~/opt/rare/Rdata/omim.rds"
df.orphanet = readRDS(orphanet.path)
df.omim = readRDS(omim.path)

tsv.path = as.character(commandArgs(TRUE)[1])
maverick.res.path = as.character(commandArgs(TRUE)[2])
out.res.path = as.character(commandArgs(TRUE)[3])

# Load tsv with input mutations
df.vcf.tsv = read.delim(file = tsv.path,skip = 1,stringsAsFactors = F)
rownames(df.vcf.tsv) = paste(df.vcf.tsv$CHROM,df.vcf.tsv$POS,df.vcf.tsv$REF,df.vcf.tsv$ALT,sep=".")
colnames(df.vcf.tsv) = sub(x = colnames(df.vcf.tsv),pattern = "dbNSFP_",replacement = "",fixed = T)
colnames(df.vcf.tsv) = sub(x = colnames(df.vcf.tsv),pattern = "_score",replacement = "",fixed = T)

# load Maverick results, fix chrs names and compute id
df.maverick = read.delim(file = maverick.res.path,stringsAsFactors = F)
df.maverick = subset(df.maverick,!is.na(hg38_chr))
df.maverick$hg38_chr[df.maverick$hg38_chr==23] = "X"
df.maverick$hg38_chr[df.maverick$hg38_chr==24] = "Y"
df.maverick$hg38_chr[df.maverick$hg38_chr==25] = "MT"
df.maverick$hg38_chr = paste0("chr",df.maverick$hg38_chr)
rownames(df.maverick) = paste(df.maverick$hg38_chr,df.maverick$hg38_pos.1.based.,df.maverick$ref,df.maverick$alt,sep = ".")

# Make final report, merging infos from orphanet and omim
df.vcf.tsv = df.vcf.tsv[rownames(df.vcf.tsv)%in%rownames(df.maverick),]
df.maverick = df.maverick[rownames(df.vcf.tsv),]
df.vcf.tsv$Maverick.Score=df.maverick$TotalScore

# first exstract vartype and symbol
df.vcf.tsv$VarType=sapply(strsplit(df.vcf.tsv$ANN,split = "|",fixed = T),function(x) x[[2]])
df.vcf.tsv$Symbol=sapply(strsplit(df.vcf.tsv$ANN,split = "|",fixed = T),function(x) x[[4]])

# Now merge orphanet and OMIM
df.vcf.tsv$Orphanet.id = df.orphanet$Orphanet.id[match(df.vcf.tsv$Symbol,df.orphanet$Symbol)]
df.vcf.tsv$Orphanet.Phenotype = df.orphanet$Phenotype[match(df.vcf.tsv$Symbol,df.orphanet$Symbol)]
df.vcf.tsv$Orphanet.inheritance = df.orphanet$inheritance_type[match(df.vcf.tsv$Symbol,df.orphanet$Symbol)]
df.vcf.tsv$OMIM.id = df.omim$OMIM.id[match(df.vcf.tsv$Symbol,df.omim$Symbol)]
df.vcf.tsv$OMIM.title = df.omim$OMIM.title[match(df.vcf.tsv$Symbol,df.omim$Symbol)]
df.vcf.tsv$OMIM.Phenotype = df.omim$Phenotype[match(df.vcf.tsv$Symbol,df.omim$Symbol)]

# sort by maverick ranks
df.vcf.tsv = df.vcf.tsv[order(df.vcf.tsv$Maverick.Score,decreasing = T),]

# Order coloumns
cols = c("CHROM","POS","REF","ALT","QUAL","VarType","Symbol","AD","DP","MED_DP","MIN_DP","GQ","GT","PL","VAF","Maverick.Score","AlphaMissense","ESM1b","EVE","REVEL","Orphanet.id","Orphanet.Phenotype","Orphanet.inheritance","OMIM.id","OMIM.title","OMIM.Phenotype")
df.vcf.tsv = cbind(df.vcf.tsv[,colnames(df.rare)%in%cols],df.vcf.tsv[,!colnames(df.vcf.tsv)%in%cols])

# save the results
write.table(df.vcf.tsv,file = out.res.path,quote = F,sep = "\t",row.names = F)




