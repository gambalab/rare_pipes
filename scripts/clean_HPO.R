#!/usr/bin/env Rscript

HPO.dict.path="~/opt/rare/Rdata/CADA_HPO_dict.rds"
df.HPO = readRDS(HPO.dict.path)

HPO.list = as.character(commandArgs(TRUE)[1])

HPO.list = unlist(strsplit(x = HPO.list,split = ",",fixed = T))
HPO.list = intersect(HPO.list,df.HPO$HPO.id)
cat(paste(HPO.list,collapse = ","))
