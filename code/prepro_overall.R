# Builds individual phenotype file and generates necessary directories for snakemake workflow
# This file would need to be updated based on the measures that the investigator is taking to adjust the continuous measure
# In this example, cytokines of interest are adjusted based on association with weighted elix score
source("code/log_smk.R") #this assigns the log file for the run
#LIBRARIES ----
library(tidyverse)
library(doFuture)

#These two lines assign proper variables for using the cluster
doFuture::registerDoFuture()
future::plan(future::multicore, workers = snakemake@resources[["ncores"]])

#FILES ----
pheno <- readr::read_csv(file = snakemake@input[["metadata"]])
phenotype_cols <- snakemake@wildcards[["phenotype"]]

paths <- paste0(c("data/pheno/",
                  "results/",
                  "benchmarks/",
                  "figures/",
                  "log/",
                  "data/mikropml/"),
                phenotype_cols)

for(i in 1:length(paths)){
  
  if(dir.exists(paths[i]) == FALSE){
    
    dir.create(paths[i])
    
  }else{
    
    print(paste0("directory ", paths[i], " already exists"))
    
  }
}

new_dir3 <- paste0("results/", phenotype_cols, "/runs")

if(dir.exists(new_dir3) == FALSE){
  
  dir.create(new_dir3)
    
  }else{
    
    print(paste0("directory ", new_dir3, " already exists"))
    
  }

#UNADJUSTED FILE ----
unadjusted_out <- pheno %>%
  select(genome_id,
         all_of(phenotype_cols),
         elix_weighted_update) %>%
  drop_na() %>%
  mutate(elix_weighted_update = NULL)
  
  write_tsv(unadjusted_out,
            file = paste0("data/pheno/", phenotype_cols, "/raw.tsv"))

#EVALUATE SIGNIFICANT ELIX ASSOCIATIONS, GENERATE ADJUSTED ----
  model <- paste0("summary(lm(",
                  phenotype_cols,
                  " ~ elix_weighted_update, data = pheno))")
  
  out_model <- eval(str2lang(model))
  
  out_coefs <- out_model$coefficients
  
  if(out_model$r.squared > 0.10){
    
    adjusted_out <- pheno %>%
      mutate(adjusted_phenotype = eval(str2lang(phenotype_cols)) - (elix_weighted_update*out_coefs[2,1] + out_coefs[1,1])) %>%
      #view()
      select(genome_id,
             adjusted_phenotype,
             elix_weighted_update) %>%
      drop_na() %>%
      mutate(elix_weighted_update = NULL) #%>%
      #view()
    
    colnames(adjusted_out) <- c("genome_id",
                                phenotype_cols)
    
    write_tsv(adjusted_out,
              file = paste0("data/pheno/", phenotype_cols, "/adjusted.tsv"))
  }
