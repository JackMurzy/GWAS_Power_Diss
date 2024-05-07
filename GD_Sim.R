library(qqman)
library(genio)
library(parallel)
library(data.table)
library(tidyverse)
setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test")
# Add PLINK to path (no need to run multiple times)
# plink_path <- "~/Desktop/University/Biomedicine/y3/Dissertation/plink_mac_20231211"
# Sys.setenv(PATH = paste(Sys.getenv("PATH"), plink_path, sep = ":"))
# Sys.getenv("PATH")
rm(list = ls())

generate_genomic_data <- function(n_individuals = 100, 
                                  n_snps = 1000, 
                                  n_chrs = 22, 
                                  heritability = 0.8, 
                                  n_causal_snps = 5,
                                  n_cases = 46,
                                  file_name = "simulated_data",
                                  parallel = FALSE,
                                  human = TRUE,
                                  chrom_data = NULL){
  # Load the timing functions
  tic <- function(){return(Sys.time())}
  toc <- function(tic){return(as.numeric(difftime(Sys.time(), tic, units = "secs")))}
  time_format <- function(seconds){
    hour_secs <- floor(seconds/3600)*3600
    min_secs <- floor((seconds - hour_secs)/60)*60
    sec_secs <- floor(seconds - (hour_secs+min_secs))
    return(paste0(hour_secs/3600,"hrs ",min_secs/60,"mins ",sec_secs,"secs"))
  }
  setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test")
  # Add PLINK to path if not already there
  plink_path <- "~/Desktop/University/Biomedicine/y3/Dissertation/plink_mac_20231211"
  if (!grepl(plink_path, Sys.getenv("PATH"))) {
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), plink_path, sep = ":"))
  }
  
  # Begin the computations
  total_start <- tic()
  cat("\014")
  
  # Load the chrom data if required
  if(human & is.null(chrom_data))
    chrom_data <- data.frame(
      chr = c(1:22),
      length = c(214066000,222889000,186938000,169035000,170954000,165022000,149414000,125148000,107440000,127894000,129193000,
                 125198000,93711000,89344000,73467000,74037000,73267000,73078000,56044000,63317000,33824000,33786000),
      snp_count = c(129931,103664,93140,84426,117882,96317,71752,57834,62013,61298,84663,
                    59245,53093,44112,37814,38735,34621,45135,25676,29478,20916,28410)
    )
  
  # If paralell in use sequester cores
  if(parallel){
    cores_avail <- detectCores()
    numCores <- cores_avail - 2  # Reserve one core for system operations
    if(n_chrs < numCores){
      numCores <- n_chrs
    }
    print(n_chrs)
    print(paste0("Parralell SNP processing: ",numCores," in use, ",cores_avail-numCores," reserved"))
    cl <- makeCluster(numCores)
    # Define and export necessary variables and functions to all worker cores
    vars_to_export <- list(n_chrs = n_chrs,
                           n_snps = n_snps, 
                           n_individuals = n_individuals,
                           chrom_data = chrom_data, 
                           tic = tic, 
                           toc = toc,
                           time_format = time_format)
    list_names <- names(vars_to_export)
    clusterExport(cl, list_names, envir = list2env(vars_to_export))
  }
  
  # Function to process each chromosome - note this needs severe work !!!
  process_chromosome <- function(chr) {
    print(paste0("Gathering data for chrom: ",chr,"/",n_chrs))
    
    # Base positions and number on chrom data if availible
    if(!is.null(chrom_data)){
      if(n_snps > chrom_data[chr,3]){
        n_snps <- chrom_data[chr,3]
      }
      # Randomly generate position based of chrom spread
      positions <- sort(sample(1:chrom_data[chr,2],n_snps))
    } else {
      positions <- sort(sample(1:9000000, n_snps))
    }
    snp_ids <- paste0("rs", chr, positions)
    
    
    # Generate SNPs, IDs and .map data
    map_data <- data.frame(CHR = chr, SNP = snp_ids, GD = 0, BPP = positions)
    
    # Generate allele data - this is separate from genotype data since it can be random
    bases <- c("A","T","C","G")
    minor_alleles <- sample(bases,n_snps,replace = TRUE)
    major_alleles <- sapply(minor_alleles, function(x) sample(bases[bases != x], 1))
    # Calculate allele vectors
    allele_vec <- mapply(function(min,maj){c(paste0(min, min), paste0(min, maj), paste0(maj, maj))},
                         minor_alleles, major_alleles, SIMPLIFY = FALSE)
    
    # Generate HWE vectors
    maf <- runif(n_snps, min = 0.05, max = 0.5)
    p <- maf
    q <- 1 - p
    p_aa <- p^2
    p_Aa <- 2 * p * q
    p_AA <- q^2
    
    first_prob <- c(p_aa[1], p_Aa[1], p_AA[1])
    # Generate random genotypes for the first SNP
    genotypes_mx <- sample(c(0, 1, 2), n_individuals, prob = first_prob, replace = TRUE)
    # Generate the alleles for first SNP
    alleles_df <- as.data.frame(allele_vec[[1]][genotypes_mx+1])
    genotypes_df <- as.data.frame(genotypes_mx)
    
    # Calculate LD metric based on distances between SNPs
    distances <- c(0, diff(positions))
    ld_metric <- 0.999^distances
    # Pre-calculate average MAF outside the loop to save computation
    average_maf <- (maf[1:(length(maf)-1)] + maf[2:length(maf)]) / 2
    
    
    run_end <- NULL
    for (i in 2:n_snps) {
      start <- tic()
      print(paste0("Generating SNP data: ",i,"/",n_snps))
      prev_genotypes <- genotypes_df[, i - 1]
      ld <- ld_metric[i]
      
      # Adjust LD weight considering MAF - not sure if this is correct, this needs work
      # Idea is that when LD is high, maf input becomes much lower
      # However if one MAF is high then LD of the other should not matter since it cant be in LD ? 
      # INVESTIGATE ME. Pretty sure this can be changed to get rid of the ld weight section
      ld_weight <- ld * (1 - average_maf[i-1])
      maf_weight <- 1 - ld_weight
      # Base probabilities on MAF
      base_probs <- c(p_aa[i], p_Aa[i], p_AA[i])
      probs <- matrix(0, nrow = n_individuals, ncol = 3)
      
      for (j in 1:n_individuals) {
        probs[j, ] <- maf_weight * base_probs
        if (prev_genotypes[j] == 0) {
          # Increase the probability of repeating the same genotype (0)
          probs[j, 1] <- probs[j, 1] + ld_weight * ld
        } else if (prev_genotypes[j] == 1) {
          # Increase the probability of repeating the same genotype (1)
          probs[j, 2] <- probs[j, 2] + ld_weight * ld
        } else {
          # Increase the probability of repeating the same genotype (2)
          probs[j, 3] <- probs[j, 3] + ld_weight * ld
        }
      }
      
      # Sample genotypes based on calculated probabilities of MAF and ld weight
      genotypes_sampled <- apply(probs, 1, function(p) sample(c(0, 1, 2), 1, prob = p))
      genotypes_df <- cbind(genotypes_df, genotypes_sampled)
      
      # Map genotype to the allele vector list and store
      allele_sampled <- allele_vec[[i]][genotypes_sampled+1]
      alleles_df <- cbind(alleles_df,allele_sampled)
      cat("\014")
      print(paste0("Completed: ",i/n_snps*100,"%"))
      run_end <- rbind(run_end,toc(start))
      mean(run_end)
      print(paste0("Estimated time remaining: ",time_format((n_snps-i)*mean(run_end))))
    }
    colnames(genotypes_df) <- snp_ids
    colnames(alleles_df) <- snp_ids
    
    print("Collating results")
    # Checks for similarity and stores allele info. Similarity is a proxy for indentical by descnet score
    similarity_scores <- numeric(ncol(genotypes_df) - 1)
    observed_maf <- numeric(ncol(genotypes_df))
    for (i in 1:ncol(genotypes_df)) {
      current_snp <- genotypes_df[, i]
      observed_maf[i] <- sum(current_snp == 1) / (2 * n_individuals) + sum(current_snp == 0) / n_individuals
      if (i > 1) {
        prev_snp <- genotypes_df[, i-1]
        similarity <- sum(current_snp == prev_snp) / n_individuals
        similarity_scores[i-1] <- similarity
      }
    }
    # Merge data into a useful data to be exported frame
    allele_info_df <- data.frame(chr = chr,id = snp_ids,pos = positions,minor = minor_alleles, major = major_alleles,
                                 MAF = maf,LD_Prev = ld_metric,Similarity = c(1,similarity_scores),Observed_maf = observed_maf)
    
    # Return all useful values
    return(list(map_data = map_data, allele_data = alleles_df, genotype_data = genotypes_df,allele_info = allele_info_df))  # Return a list of data
  }
  
  # Use parallel or sequential processing
  results <- if (parallel) {
    parLapply(cl, 1:n_chrs, process_chromosome)
  } else {
    lapply(1:n_chrs, process_chromosome)
  }
  
  # Stop and close the cluster
  if (parallel) {
    stopCluster(cl)
  }
  
  # Assign Phenotypes
  print("Assigning phenotypes")
  allele_infos <- lapply(results, `[[`, "allele_info")
  total_allele_info <- do.call(rbind,allele_infos)
  # How many causal SNPs are there
  filt_snps <- total_allele_info %>%
    filter(Observed_maf <= 0.15,
           Observed_maf >= 0.05,
           LD_Prev <= 0.8,
           LD_Prev >= 0.2)
  snps_remain <- nrow(filt_snps)
  if(snps_remain < n_causal_snps){
    n_causal_snps_adj <- snps_remain
    print(paste0("WARNING: Filters mean only ",snps_remain," potentially causal SNPs remain, less than the ",n_causal_snps," desired"))
  } else {
    n_causal_snps_adj <- n_causal_snps
  }
  # Causal SNPs are sampled based on a probability vector which conssiders their LD to previous SNP (high LD = more likely causal)
  causal_snp_data <- filt_snps %>%
    arrange(desc(LD_Prev)) %>%
    group_by(chr) %>%
    mutate(rank_weight = 1 /row_number()) %>%
    ungroup() %>%
    slice_sample(n = n_causal_snps_adj, weight_by = rank_weight, replace = FALSE) %>%
    mutate(Causal = T)
  # Function above needs lots of work. Look into modifying causal probability distribution
  
  # Get the right information to adjust the phenotypes accordingly
  indexed_causal_data <- total_allele_info %>%
    mutate(index_col = row_number()) %>%
    right_join(causal_snp_data)
  genotype_total <- do.call(cbind,lapply(results, `[[`, "genotype_data"))
  # Assign phenotype weights based on the normalised sum of genotypes and the heritability weight
  geno_index_data <- genotype_total[,indexed_causal_data$index_col] %>%
    mutate(person_index = row_number())
  # Note the bottom filters out all cases where the genotype score is max, and we square scores to make differnces bigger
  geno_pheno_data_filt <- geno_index_data %>%
    mutate(phenotype_vals = rowSums(.)-person_index,
           weight_genetic = (1 / phenotype_vals)^2,
           weight_genetic = weight_genetic / sum(weight_genetic),
           weight_uniform = 1 / n(),
           weight = heritability * weight_genetic + (1 - heritability) * weight_uniform) %>%
    filter(phenotype_vals < n_causal_snps_adj*2) # Filter out any individuals where all rows are 2
  # Note the code above breaks if n_causal = 1 since a LOT of individuals shoudl have 2 genotype assuming HWE  
                                 
  # Send an error messsage if we need to return fewer cases
  potential_cases <- nrow(geno_pheno_data_filt)
  if(n_cases > potential_cases){
    n_cases_adj <- potential_cases
    print(paste0("WARNING: Filters mean only ",potential_cases," potential cases remain, less than the ",n_cases," desired"))
  } else {
    n_cases_adj <- n_cases
  }
  
  # Assign and arrange out phenotypes
  geno_pheno_data <- geno_pheno_data_filt %>%
    slice_sample(n = n_cases_adj, weight_by = weight, replace = FALSE) %>%
    mutate(phenotype = 2) %>%
    right_join(geno_index_data) %>%
    arrange(person_index) %>%
    mutate(phenotype = if_else(is.na(phenotype), 1, phenotype))
  print("Phenotypes assigned !")
  
  # Extract and combine data from results
  print("Merging map file")
  map_datas <- lapply(results, `[[`, "map_data")
  plink_map_data <- do.call(rbind, map_datas)
  fwrite(plink_map_data, file = paste0(file_name, ".map"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  print("Merging ped file")
  ped_easy <- data.frame(FID = 1:n_individuals, IID = 1:n_individuals, PID = 0, MID = 0, Sex = sample(1:2,n_individuals,replace = TRUE))
  alleles_total <- lapply(results, `[[`, "allele_data")
  plink_ped_data <- cbind(ped_easy,geno_pheno_data$phenotype,alleles_total)
  fwrite(plink_ped_data, file = paste0(file_name, ".ped"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Write out the final bfile
  print(paste0("Making final bfile with name: '",file_name,"'"))
  system(paste0("plink --file ",file_name," --make-bed --out ",file_name))
  genomic_data <- read_plink("simulated_data")
  print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
  print(paste0("File name: '",file_name,"'"))
  
  # Collect the allele info too and add causality
  allele_infos <- do.call(rbind,lapply(results, `[[`, "allele_info")) %>%
    left_join(causal_snp_data) %>%
    select(-rank_weight)
  genomic_data$allele_info <- allele_infos
  genomic_data$causal_SNP_weights 
  
  return(genomic_data)
}

## Quick test of results
genomic_data <- generate_genomic_data(n_individuals = 200, 
                                      n_snps = 10000, 
                                      n_chrs = 22,
                                      heritability = 1, 
                                      n_causal_snps = 3,
                                      n_cases = 46,
                                      file_name = "simulated_data",
                                      parallel = TRUE,
                                      human = TRUE,
                                      chrom_data = NULL)


system("plink --bfile simulated_data --assoc --ci 0.95 --maf 0.05 --nonfounders --out sim_assoc_results")
sim_results_assoc <- fread("sim_assoc_results.assoc", head=TRUE)

genomic_data$allele_info %>%
  filter(Causal) %>%
  select("SNP" = id) %>%
  inner_join(sim_results_assoc)

sim_results_assoc %>%
  arrange(P) %>%
  slice(1:5) 

qq(sim_results_assoc$P, main = "Q-Q plot of GWAS p-values : log")
manhattan(sim_results_assoc,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc", col = c("lightblue","blue"))
head(sim_results_assoc)
sum(!sim_results_assoc[!is.finite(sim_results_assoc$P)])


causal_SNPs <- genomic_data$allele_info %>%
  filter(Causal) %>%
  select("SNP" = id)

sim_results_assoc %>%
  arrange(P)
  inner_join(causal_SNPs)
  arrange(P) %>%
  slice(1:5)


sim_results_assoc_check  <- fread("simulated_data.ped", head=TRUE)
head(sim_results_assoc_check$rs13036)

sim_results_assoc_check[1:100,1191887]

head(genomic_data$bim)

genomic_data$bim %>%
  mutate(row = row_number()) %>%
  filter(id == "rs1820925840")

1191881





# 
# 
#   map_datas <- list()
#   alleles_total <- list()
#   for(chr in 1:n_chrs){
#     print(paste0("Gathering data for chrom: ",chr,"/",n_chrs))
#     
#     positions <- sort(sample(1:9000000, n_snps))
#     snp_ids <- paste0("rs",chr,positions)
#     genetic_distance <- 0 # No need to set this, its not used in --assoc
#     map_datas[[chr]] <- data.frame(CHR = chr, SNP = snp_ids,GD = 0, BPP = positions)
# 
# 
#     # Generate allele data
#     bases <- c("A","T","C","G")
#     minor_alleles <- sample(bases,n_snps,replace = TRUE)
#     major_alleles <- sapply(minor_alleles, function(x) sample(bases[bases != x], 1))
#     # Calculate allele vectors
#     allele_vec <- mapply(function(min,maj){c(paste0(min, min), paste0(min, maj), paste0(maj, maj))},
#                          minor_alleles, major_alleles, SIMPLIFY = FALSE)
#     
#     # Generate HWE vectors
#     maf <- runif(n_snps, min = 0.05, max = 0.5)
#     p <- maf
#     q <- 1 - p
#     p_aa <- p^2
#     p_Aa <- 2 * p * q
#     p_AA <- q^2
#     
#     first_prob <- c(p_aa[1], p_Aa[1], p_AA[1])
#     # Generate random genotypes for the first SNP
#     genotypes_mx <- sample(c(0, 1, 2), n_individuals, prob = first_prob, replace = TRUE)
#     # Generate the alleles for first SNP
#     alleles_df <- as.data.frame(allele_vec[[1]][genotypes_mx+1])
#     genotypes_df <- as.data.frame(genotypes_mx)
#     
#     # Calculate LD metric based on distances between SNPs
#     distances <- c(0, diff(positions))
#     ld_metric <- 0.999^distances
#     # Pre-calculate average MAF outside the loop to save computation
#     average_maf <- (maf[1:(length(maf)-1)] + maf[2:length(maf)]) / 2
#     
#     
#     run_end <- NULL
#     for (i in 2:n_snps) {
#       start <- tic()
#       print(paste0("Generating SNP data: ",i,"/",n_snps))
#       prev_genotypes <- genotypes_df[, i - 1]
#       ld <- ld_metric[i]
#       
#       # Adjust LD weight considering MAF
#       ld_weight <- ld * (1 - average_maf[i-1])
#       maf_weight <- 1 - ld_weight
#       # Base probabilities on MAF
#       base_probs <- c(p_aa[i], p_Aa[i], p_AA[i])
#       
#       for (j in 1:n_individuals) {
#         probs[j, ] <- maf_weight * base_probs
#         if (prev_genotypes[j] == 0) {
#           # Increase the probability of repeating the same genotype (0)
#           probs[j, 1] <- probs[j, 1] + ld_weight * ld
#         } else if (prev_genotypes[j] == 1) {
#           # Increase the probability of repeating the same genotype (1)
#           probs[j, 2] <- probs[j, 2] + ld_weight * ld
#         } else {
#           # Increase the probability of repeating the same genotype (2)
#           probs[j, 3] <- probs[j, 3] + ld_weight * ld
#         }
#       }
#       
#       # Sample genotypes based on calculated probabilities
#       genotypes_sampled <- apply(probs, 1, function(p) sample(c(0, 1, 2), 1, prob = p))
#       genotypes_df <- cbind(genotypes_df, genotypes_sampled)
#       
#       # Map genotype to the allele vector list and store
#       allele_sampled <- allele_vec[[i]][genotypes_sampled+1]
#       alleles_df <- cbind(alleles_df,allele_sampled)
#       cat("\014")
#       print(paste0("Completed: ",i/n_snps*100,"%"))
#       run_end <- rbind(run_end,toc(start))
#       mean(run_end)
#       print(paste0("Estimated time remaining: ",time_format((n_snps-i)*mean(run_end))))
#     }
#     colnames(genotypes_df) <- snp_ids
#     colnames(alleles_df) <- snp_ids
#     alleles_total[[chr]] <- alleles_df
#   }
#   
#   # Combine and make the map file
#   print("Merging map file")
#   plink_map_data <- do.call(rbind, map_datas)
#   fwrite(plink_map_data, file = paste0(file_name, ".map"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#   
#   # Combine and make the ped file 
#   print("Merging ped file")
#   ped_easy <- data.frame(FID = 1:n_individuals, IID = 1:n_individuals, PID = 0, MID = 0, Sex = sample(1:2,n_individuals,replace = TRUE))
#   phenotypes <- data.frame(P = sample(1:2,n_individuals, replace = TRUE))
#   combined_alleles <- do.call(cbind, alleles_total)
#   plink_ped_data <- cbind(ped_easy,phenotypes,combined_alleles)
#   fwrite(plink_ped_data, file = paste0(file_name, ".ped"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#   
#   # Write out the final bfile
#   print(paste0("Making final bfile with name: '",file_name,"'"))
#   system(paste0("plink --file ",file_name," --make-bed --out ",file_name))
#   genomic_data <- read_plink("simulated_data")
#   print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
#   print(paste0("File name: '",file_name,"'"))
#   return(genomic_data)
# }






