entrezToSymbol = function(organism="human")
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
  f = genes(edb,return.type = "data.frame")
  f = subset(f,!is.na(entrezid) & !gene_biotype%in%"LRG_gene")
  ens = unlist(f$entrezid)
  
  ens.map <- data.frame(ensID=names(ens),entrezID=ens,stringsAsFactors = F)
  ens.map$symbol = f$symbol[fastmatch::fmatch(ens.map$ensID,f$gene_id)]
  ens.map = subset(ens.map,!is.na(symbol) & symbol!="" & ensID!="LRG_1231")
  return(ens.map)
}

ens.map = entrezToSymbol()
ens.map$entrezID = paste0("Entrez:",ens.map$entrezID)
df = read.csv("~/src/rare_pipes/data/CADA.dict.csv",stringsAsFactors = F)
ens.map = subset(ens.map,entrezID%in%df$gene_id)

saveRDS(ens.map,"~/src/rare_pipes/Rdata/CADA_dict.rds")
