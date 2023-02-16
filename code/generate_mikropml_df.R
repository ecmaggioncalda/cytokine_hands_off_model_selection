#mikropml input file generation
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#PHENO ----
paste0(snakemake@input[["pheno"]])
pheno <- readr::read_delim(file = snakemake@input[["pheno"]],
                           delim = "\t")

pheno_merge <- as.data.frame(pheno)
rownames(pheno_merge) <- pheno_merge$genome_id
pheno_merge <- pheno_merge[, -1, drop = FALSE]

#GENO ----
paste0(snakemake@params[['path']])
geno <- read.delim(file = snakemake@params[['path']],
                   row.names = 1)

#FEATURE SELECTION ----
feature <- read_csv(snakemake@input[["feature_file"]])

feature_sub <- feature %>%
  select(full_names,
         all_of(colnames(pheno_merge))) %>%
  filter(if_any(where(is.numeric), ~.x == 1))

index <- sapply(feature_sub$full_names,
                function(x){
                  
                  which(x == geno$variant)
                  
                })

if(length(index) == length(feature_sub$full_names)){
  stop()
}
if(!any(is.na(index))){
  stop()
}

geno_sub <- geno[index,] %>%
  select(variant,
         all_of(rownames(pheno_merge))) %>%
  mutate(variant = gsub("_$", "", variant)) %>%
  column_to_rownames("variant")

geno_merge <- t(geno_sub)

if(sum(rownames(pheno_merge) %in% rownames(geno_merge)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

index <- match(rownames(pheno_merge), rownames(geno_merge))

geno_ordered <- geno_merge[index, , drop = FALSE]

if(sum(rownames(pheno_merge) == rownames(geno_ordered)) != length(rownames(pheno_merge))){
  stop("mismatch between pheno and geno contents")
}

complete_frame <- cbind(pheno_merge,
                        geno_ordered)

print("complete frame generated with pheno:geno, export to csv for mikropml preprocessing")

#GENERATE FILES ----
#patient and genome factors
write_csv(complete_frame,
          file = snakemake@output[['file_name']],
          col_names = TRUE)
