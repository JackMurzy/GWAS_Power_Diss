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
                                  n_snps = 100000, 
                                  chrs = 1:22, 
                                  heritability = 0.8, 
                                  n_causal_snps = 5,
                                  n_cases = 46,
                                  file_name = "simulated_data",
                                  parallel = FALSE,
                                  human = TRUE,
                                  chrom_data = NULL,
                                  coverage = 1,
                                  seed = NULL){
  # Set the seed
  set.seed(seed)
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
  
  # How many chroms are we dealing with
  n_chrs <- length(chrs)
  
  # Load the chrom data if required
  if(human & is.null(chrom_data)){
    chrom_data <- data.frame(
      chr = c(1:22),
      length = c(214066000,222889000,186938000,169035000,170954000,165022000,149414000,125148000,107440000,127894000,129193000,
                 125198000,93711000,89344000,73467000,74037000,73267000,73078000,56044000,63317000,33824000,33786000),
      snp_count = c(129931,103664,93140,84426,117882,96317,71752,57834,62013,61298,84663,
                    59245,53093,44112,37814,38735,34621,45135,25676,29478,20916,28410)
    )
  }
  
  if(!is.null(chrom_data)){
    # Scale data based off the coverage factor
    chrom_data$length <- coverage*chrom_data$length
    # Create snp matrix based off chromosome size and count
    snp_matrix <- round(n_snps/sum(chrom_data$length)*chrom_data$length)
    chrom_data$snp_count <- snp_matrix
  } else {
    # Set snp matrix as equal across all chromosomes
    snp_matrix <- rep(n_snps,n_chrs)
  }
  
  # If paralell in use sequester cores
  numCores <- 1
  if(parallel){
    cores_avail <- detectCores()
    numCores <- cores_avail - 2  # Reserve one core for system operations
    if(n_chrs < numCores){
      numCores <- n_chrs
    }
    print(paste0("Parralell SNP processing: ",numCores," cores in use, ",cores_avail-numCores," reserved"))
    
    # Run a pid check to see how big clusters are
    pid_check <- function(chr){data.frame(chr = chr, id = Sys.getpid())}
    cl <- makeCluster(numCores)
    pid_check_res <- do.call(rbind,parLapply(cl, chrs, pid_check))
    pid_df <- as.data.frame(table(pid_check_res$id))
    cluster_sizes <- rep(pid_df$Freq)
    
    # What is the optimal sum of values for each group
    target_sum <- sum(snp_matrix) / length(cluster_sizes)
    
    # Use a greedy algorithm to assign optimal cluster groups
    cluster_groups <- vector("list", length(cluster_sizes))
    remaining_values <- snp_matrix
    remaining_chrs <- chrs
    for (i in seq_along(cluster_sizes)) {
      cluster_sum <- 0
      cluster_values <- c()
      cluster_chroms <- c()
      # Find value from remaining_values that (when added to the group) minimises the difference 
      # between the group sum and the target sum defined as the optimal length. 
      while (length(cluster_values) < cluster_sizes[i] && length(remaining_values) > 0) {
        min_diff_index <- which.min(abs(remaining_values + cluster_sum - target_sum))
        cluster_values <- c(cluster_values, remaining_values[min_diff_index])
        cluster_sum <- sum(cluster_values)
        remaining_values <- remaining_values[-min_diff_index]
        
        # Store chroms for the actual output
        cluster_chroms <- c(cluster_chroms, remaining_chrs[min_diff_index])
        remaining_chrs <- remaining_chrs[-min_diff_index]
      }
      cluster_groups[[i]] <- cluster_chroms
    }
    # Store the ordered chroms
    ordered_chrs <- do.call(c,cluster_groups)
    
    # Define and export necessary variables and functions to all worker cores
    clusterEvalQ(cl, library("tidyverse"))
    vars_to_export <- list(n_chrs = n_chrs,
                           snp_matrix = snp_matrix, 
                           n_individuals = n_individuals,
                           chrom_data = chrom_data, 
                           seed = seed,
                           tic = tic, 
                           toc = toc,
                           time_format = time_format)
    list_names <- names(vars_to_export)
    clusterExport(cl, list_names, envir = list2env(vars_to_export))
  }
  
  # Function to process each chromosome
  process_chromosome <- function(chr) {
    set.seed(seed)
    print(paste0("Gathering data for chrom: ",chr,"/",n_chrs))
    n_snps <- snp_matrix[chr]
    # Base positions and number on chrom data if availible
    if(!is.null(chrom_data)){
      # The code below ensures the number of SNPs is based off the actual count, but 
      # it makes more sense to base the number of SNPs based off the size of the chromosome
      # if(n_snps > chrom_data[chr,3]){
      #   n_snps <- chrom_data[chr,3]
      # }
      positions <- sort(sample(1:chrom_data$length[chr],n_snps))
    } else {
      # Set positions equally across 125,000,000 SNPs (set based off a rough average of chrom sizes in human data)
      positions <- sort(sample(1:150000000, n_snps))
    }
    snp_ids <- paste0("rs",chr,"_", positions)
    
    
    # Generate SNPs, IDs and .map data
    map_data <- data.frame(CHR = chr, SNP = snp_ids, GD = 0, BPP = positions)
    
    # Generate allele data
    bases <- c("A","T","C","G")
    minor_alleles <- sample(bases,n_snps,replace = TRUE)
    major_alleles <- sapply(minor_alleles, function(x) sample(bases[bases != x], 1))
    # Calculate allele vectors
    allele_vec <- mapply(function(min,maj){c(paste0(min, min), paste0(min, maj), paste0(maj, maj))},
                         minor_alleles, major_alleles, SIMPLIFY = FALSE)
    
    # Assign MAF based off a beta distribution (ensuring we have a large enough pool to filter from)
    maf_initial <- rbeta(n_snps*100, 1, 10)
    maf_filtered <- maf_initial[maf_initial >= 0.05 & maf_initial <= 0.5]
    maf <- sample(maf_filtered, n_snps)

    # Generate HWE vectors
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
    # average_maf <- (maf[1:(length(maf)-1)] + maf[2:length(maf)]) / 2
    
    
    run_end <- NULL
    for (i in 2:n_snps) {
      start <- tic()
      print(paste0("Generating Chr",chr," SNP data: ",i,"/",n_snps))
      prev_genotypes <- genotypes_df[, i - 1]
      ld <- ld_metric[i]

      # Since maf is random and LD is "correct' based off distance, adjust the
      # maf weight based on LD so that when LD is high, it is considered more
      # maf_weight <- 1-ld
      
      # Adjust LD weight considering MAF
      # ld_weight <- ld * (1 - average_maf[i-1])
      # maf_weight <- 1 - ld_weight
      
      # Base probabilities on MAF
      this_probs <- c(p_aa[i], p_Aa[i], p_AA[i])
      prev_probs <- c(p_aa[i-1], p_Aa[i-1], p_AA[i-1])
      # Since MAF was decided randomly, when LD is high it should be closer to the prev maf instead
      base_probs <- ld * prev_probs + (1 - ld) * this_probs
      
      # # Create a probability matrix for all individuals across 0,1,2 genotypes
      # probs <- matrix(0, nrow = n_individuals, ncol = 3)
      # for (j in 1:n_individuals) {
      #   # Note this used to be: maf_weight * base_probs -- when this was in use !!
      #   probs[j, ] <- base_probs
      #   # Note this used to be: + ld * ld_weight -- when this was in use !!
      #   if (prev_genotypes[j] == 0) {
      #     # Increase the probability of repeating the same genotype (0) 
      #     probs[j, 1] <- probs[j, 1] + ld
      #   } else if (prev_genotypes[j] == 1) {
      #     # Increase the probability of repeating the same genotype (1)
      #     probs[j, 2] <- probs[j, 2] + ld
      #   } else {
      #     # Increase the probability of repeating the same genotype (2)
      #     probs[j, 3] <- probs[j, 3] + ld
      #   }
      # }
      
      # Create a probability list for all base_probs for when prev genotype is 0,1,2 (stored as 1,2,3 respectively) 
      # This does the exact same as the code above but is more efficent
      prob_list <- lapply(1:3, function(i) {
        base_probs[i] <- base_probs[i] + ld
        base_probs
      })
      # Call these probabilities based off the previous genotypes
      probs <- do.call(rbind, prob_list[prev_genotypes+1])
      
      # Sample genotypes based on calculated probabilities
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
    # Checks for similarity and stores allele info
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
    parLapply(cl, ordered_chrs, process_chromosome)
  } else {
    lapply(chrs, process_chromosome)
  }
  

  # Work out which SNPs are causal
  print("Assigning causal SNPs")
  
  allele_infos <- lapply(results, `[[`, "allele_info")
  total_allele_info <- do.call(rbind,allele_infos)
  all_snps <- total_allele_info %>%
    mutate(index_col = row_number())
  
  # Get genotypes in
  genotype_total <- do.call(cbind,lapply(results, `[[`, "genotype_data"))
  
  # Work out all the weightings for each genotype
  print(paste0("Calculating probability metrics for causal SNP assignment across all ",nrow(all_snps)," SNPs"))
  split_snp_data <- split.default(genotype_total, rep(1:numCores, each = ceiling(ncol(genotype_total) / numCores))[1:ncol(genotype_total)])
  weight_snps <- function(snp_data_input) {
    data.frame(sd = apply(snp_data_input, 2, sd)) %>%
      mutate(sd_norm = (sd - min(sd)) / (max(sd) - min(sd))) %>%
      mutate(two_count = colSums(snp_data_input == 2)) %>%
      mutate(two_norm = 1 - (two_count - min(two_count)) / (max(two_count) - min(two_count))) %>%
      mutate(sd_count_weight = sd_norm + two_norm) %>%
      rownames_to_column("id")
  }
  # If parallel then run over all cores for speed if data is large
  weighted_snp_data <- if (parallel & ncol(genotype_total) >= 100000) {
    parLapply(cl, split_snp_data, weight_snps)
  } else {
    lapply(split_snp_data, weight_snps)
  }
  # Collect the data and join to all_snps
  snp_data <- bind_rows(weighted_snp_data) %>%
    left_join(all_snps) %>%
    select(chr, id, sd, sd_norm, two_count, two_norm, sd_count_weight,
           pos, minor, major, "maf_assigned" = MAF, "maf_observed" = Observed_maf, 
           "ld_toprev" = LD_Prev, "similarity_toprev" = Similarity, index_col)
  
  # Close the connections now they are un-needed
  if(parallel){
    stopCluster(cl)
  }
  
  # Apply some categorical filters
  print(paste0("Applying categorical filter for causal SNPs based on: ",0.05," ≤ MAF ≤ ",0.15," & ",0.2," ≤ LD ≤ ",0.8))
  # Filter SNPs based on their maf and LD to ensure they aren't overly rare and have a mid-range LD structure
  filt_snps <- all_snps %>%
    filter(Observed_maf <= 0.15,
           Observed_maf >= 0.05)
           # LD_Prev <= 0.8,
           # LD_Prev >= 0.2)
  # Send a warning message if we are left with too few SNPs (likely means we are underpowered)
  snps_remain <- nrow(filt_snps)
  if(snps_remain < n_causal_snps){
    n_causal_snps_adj <- snps_remain
    warning(paste0("Filters mean only ",snps_remain," potentially causal SNPs remain, less than the ",n_causal_snps," desired"))
  } else {
    n_causal_snps_adj <- n_causal_snps
  }
  
  print("Applying categorical filter for causal SNPs to eliminate SNPs where all genotypes are 2 & where no genotype is 0")
  # Filter SNPs to ensure no SNP full of 2 is causal and no SNP without a 0 can be causal (ensures max range)
  pot_causal <- genotype_total[, filt_snps$index_col, drop = FALSE] %>%
    select(where(~ !all(. == 2))) %>% 
    select(where(~ any(. == 0)))
  # If we've filtered even further than before then say so now !
  pot_causal_num <- ncol(pot_causal)
  if(pot_causal_num < n_causal_snps_adj){
    warning(paste0("Further filters mean only ",pot_causal_num," potentially causal SNPs remain, less than the ",n_causal_snps_adj," desired"))
    n_causal_snps_adj <- pot_causal_num
  }
  
  # Pick the causal SNP from our subset of SNPs
  causal_snp_data <- snp_data %>%
    filter(id %in% colnames(pot_causal)) %>%
    group_by(chr) %>%
    arrange(desc(sd_count_weight)) %>%
    mutate(chr_weight = 1/row_number()) %>%
    ungroup() %>%
    mutate(rank_weight = sd_count_weight + chr_weight*0.25) %>%
    mutate(rank_weight_norm = (rank_weight - min(rank_weight)) / (max(rank_weight) - min(rank_weight))) %>%
    {
      causal_ids <- sample(.$id,size = n_causal_snps_adj, prob = .$rank_weight_norm, replace = FALSE)
      mutate(., causal = id %in% causal_ids)
    }
  
  causal_snps <- causal_snp_data %>% 
    filter(causal)
  
  print(paste0("Assigning probabilities for each individual being a case based on genotype data (",heritability*100,"%) and environmental factors (",(1-heritability)*100,"%)"))
  # Note the bottom filters out all cases where the genotype score is max, and we square scores to make differnces bigger
  geno_index_data <- pot_causal[,causal_snps$id,drop = FALSE] %>%
    mutate(index_row = row_number())
  geno_pheno_weights <- geno_index_data %>%
    mutate(phenotype_vals = rowSums(.)-index_row) %>%
    mutate(weight_genetic = (1- (phenotype_vals - min(phenotype_vals)) / (max(phenotype_vals) - min(phenotype_vals)))^2) %>%
    mutate(weight_uniform = 0.5) %>%
    mutate(weight = heritability * weight_genetic + (1 - heritability) * weight_uniform) %>%
    mutate(weight = if_else(weight<=0, 0.000000001, weight))
  
  # Send an error message if too many 0 cases are seen (means we needed more individuals or fewer cases)
  good_cases <- n_individuals-length(which(geno_pheno_weights$weight==0.000000001))
  if(n_cases > good_cases){
    warning(paste0("Filters mean only ",good_cases," good potential cases remain, less than the ",n_cases," desired"))
  }
  
  print(paste0("Assigning ",n_cases," cases based on the genotype data of ",n_causal_snps_adj," SNPs across ",n_individuals," individuals"))
  # Assign and arrange out phenotypes
  geno_pheno_data <- geno_pheno_weights %>%
    {
      causal_ids <- sample(.$index_row,size = n_cases, prob = .$weight, replace = FALSE)
      mutate(., phenotype = ifelse(index_row %in% causal_ids,2,1))
    }
  
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
  allele_info_fnal <- snp_data %>%
    left_join(causal_snp_data) %>%
    arrange(chr)
  genomic_data$allele_info <- allele_info_fnal
  genomic_data$causal_genotype <- geno_pheno_data
  genomic_data$causal_snps <- as.data.frame(causal_snps)
  return(genomic_data)
}

## Quick test of results
genomic_data <- generate_genomic_data(n_individuals = 400, 
                                      n_snps = 22000, 
                                      chrs = 1:22,
                                      heritability = 0.8, 
                                      n_causal_snps = 5,
                                      n_cases = 166,
                                      file_name = "simulated_data",
                                      parallel = TRUE,
                                      human = TRUE,
                                      chrom_data = NULL,
                                      coverage = 0.01,
                                      seed = 123)
system("plink --bfile simulated_data --assoc --ci 0.95 --maf 0.05 --nonfounders --out sim_assoc_results")
sim_results_assoc <- fread("sim_assoc_results.assoc", head=TRUE)

genomic_data$allele_info %>%
  filter(causal) %>%
  select("SNP" = id) %>%
  inner_join(sim_results_assoc)
sim_results_assoc %>%
  arrange(P) %>%
  slice(1:5) 

qq(sim_results_assoc$P, main = "Q-Q plot of GWAS p-values : log")
manhattan(sim_results_assoc,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc", col = c("lightblue","blue"))
head(sim_results_assoc)
sum(!sim_results_assoc[!is.finite(sim_results_assoc$P)])



