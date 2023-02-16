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

# if(snakemake@wildcards[['genome']] == "pan"){
#   
#   print(paste0("using pan genome path:", snakemake@params[['pan_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['pan_path']],
#                      row.names = 1)
#   
# }else if(snakemake@wildcards[['genome']] == "core"){
#   
#   print(paste0("using core genome path:", snakemake@params[['core_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['core_path']],
#                      row.names = 1)
#   
# }else if(snakemake@wildcards[['genome']] == "gene"){
#   
#   print(paste0("using gene genome path:", snakemake@params[['gene_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['gene_path']],
#                      row.names = 1)
#   
# }else if(snakemake@wildcards[['genome']] == "struct"){
#   
#   print(paste0("using gene genome path:", snakemake@params[['struct_path']]))
#   
#   geno <- read.delim(file = snakemake@params[['struct_path']],
#                      row.names = 1)
#   
# }

print("Geno matrix has completed read in")

geno_merge <- t(geno)

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
