mimTitles = read.delim("~/src/rare_pipes/data/mimTitles.txt",stringsAsFactors = F,skip = 2,header = T)
genemap2 = read.delim("~/src/rare_pipes/data/genemap2.txt",stringsAsFactors = F,skip = 3,header = T)
morbidmap = read.delim("~/src/rare_pipes/data/morbidmap.txt",stringsAsFactors = F,skip = 3,header = T)


mim2gene = read.delim("~/src/rare_pipes/data/mim2gene.txt",stringsAsFactors = F,skip = 4,header = T)
colnames(mim2gene) = c("OMIM.id","OMIM.entry","Entrez.id","Symbol","ENS.id")
mim2gene = subset(mim2gene,ENS.id!="")

# add the title for each entry and clean
mim2gene$OMIM.gene.title = mimTitles$Preferred.Title..symbol[match(mim2gene$OMIM.id,mimTitles$MIM.Number)]
mim2gene$OMIM.gene.title = sapply(strsplit(x = mim2gene$OMIM.gene.title,split = ";",fixed = T), function(x) x[[1]])

# add Phenotype
morbidmap = morbidmap[!duplicated(morbidmap$MIM.Number),]
mim2gene$Phenotype = morbidmap$X..Phenotype[match(mim2gene$OMIM.id,morbidmap$MIM.Number)]

# take only genes with a phenotype and rm dup
mim2gene = subset(mim2gene,!is.na(Phenotype))
mim2gene = mim2gene[!duplicated(mim2gene$ENS.id),]

saveRDS(mim2gene,"~/src/rare_pipes/Rdata/omim.rds")
