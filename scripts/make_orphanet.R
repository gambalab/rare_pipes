library(xml2)
library(dplyr)
library(xmlconvert)

ensToSymbol = function(df,col,organism,verbose=T)
{
  require(AnnotationHub)
  require(fastmatch)
  organism = tolower(organism)
  organism = base::match.arg(arg = organism,choices = c("human","mouse"),several.ok = F)
  org.map = c("Homo Sapiens","Mus Musculus")
  names(org.map) = c("human","mouse")
  
  message("... Retrieving gene annotation from AnnotationHub()")
  ah <- AnnotationHub::AnnotationHub()
  ahDb <- AnnotationHub::query(ah,pattern = c(org.map[organism],"EnsDb"), ignore.case = TRUE)
  id <- tail(rownames(AnnotationHub::mcols(ahDb)),n=1)
  edb <- ah[[id]]
  ens.map <- subset(genes(edb,return.type = "data.frame"),seq_name%in%c(as.character(1:22),"X","Y","MT") & !gene_biotype%in%"LRG_gene")
  
  if(organism %in% "human")
  {
    message(".. Start converting human ensamble id to symbols")
    df$symb = NA
    df$symb = ens.map$gene_name[fastmatch::fmatch(df[,col],ens.map$gene_id)]
    message("Done!")
  }
  
  if (organism %in% "mouse")
  {
    message(".. Start converting mouse ensamble id to symbols")
    df$symb = NA
    df$symb = ens.map$gene_name[fastmatch::fmatch(df[,col],ens.map$gene_id)]
    message("Done!")
  }
  
  return(df)
}


# Download last version from https://www.orphadata.com/genes/
data <- xml2::read_xml("~/src/rare_docker/data/en_product6.xml")

# Point locations
point <- xml_find_all(x = data,xpath = "//Disorder")
res <- NULL
for (i in 1:length(point)) {
  print(i)
  d <- xml2::xml_children(point[i])
  orphaCode <- NA
  disease <- NA
  for (j in 1:length(d)) {
    nm <- xml2::xml_name(x = d[j])
    if (nm=="OrphaCode"){
      orphaCode <- xml2::xml_integer(x = d[j])
    }
    if (nm=="Name"){
      disease <- xml2::xml_text(x = d[j])
    }
    
    if (nm=="DisorderGeneAssociationList"){
      g <- xml2::xml_children(xml2::xml_children(d[j]))
      ix <- which(xml_name(g) %in% "Gene")
      g <- xml2::xml_children(xml2::xml_children(g[ix]))
      ix <- grepl(pattern = "Ensembl",ignore.case = F,x = xml_text(g),fixed = T)
      ens <- gsub(pattern = "Ensembl",replacement = "",x = xml_text(g)[ix],fixed = T)
      if(length(ens)>0)
        res <- rbind(res,data.frame(orphaID=orphaCode,disease_name=disease,ensID=ens,stringsAsFactors = F))
    }
  }
}

res$orphaID <- paste0("ORPHA:",res$orphaID)
res <- ensToSymbol(df = res,col = "ensID",organism = "human",verbose = T)

# collapse by Ens
u = unique(res$ensID)
df = data_frame(ENS.id=u,Orphanet.id=NA,Phenotype=NA,Symbol=NA)
for (i in 1:length(u)) {
  tmp = subset(res,ensID%in%u[i])
  df$Orphanet.id[i] = paste(tmp$orphaID,collapse = "|")
  df$Phenotype[i] = paste(tmp$disease_name,collapse = "|")
  df$Symbol[i] = paste(unique(tmp$symb),collapse = "|")
}

saveRDS(df,"~/src/rare_docker/Rdata/orphanet.rds")



