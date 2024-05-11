library(qqman)
library(genio)
library(parallel)
library(data.table)
library(tidyverse)
setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test")
rm(list = ls())

#' Simulate genomic data
#' 
#' @description
#' `generate_genomic_data` Creates bed, bim, fam, map and ped files containing simulated genomic data
#' @details
#' This function will generate varying sets of genomic data according to user requirements. The length of time the function takes 
#' to run is proportional to the number of individuals and SNPs it is asked to simulate as well as the speed of the laptop. Parallel processing
#' has been implemented for speed, specifically this will use all but 2 cores when parallel is set to true. The output files are designed to work with
#' PLINK 1.9. NOTE that this function cant simulate chrom 23 data since this would have different genotype values depending on sex.
#' @param n_individuals (numeric) The number of individuals to simulate SNPs across
#' @param n_snps (numeric) The total number of SNPs to simulate across individuals, this is distributed across chromosomes in correspondence to their size to 
#' ensure similar coverage
#' @param chrs (numeric vector) Vector of chromosomes to simulate
#' @param n_causal_snps (numeric) The number of SNPs assigned as causal
#' @param n_cases (numeric) The number of individuals who will have case phenotypes (the rest will be controls)
#' @param heritability (numeric) A weight between 0 and 1 for how much of the disease is decided by genetics or the environment where 1 = purely genetics. 
#' @param coverage (numeric) A weight between 0 and 1 to scale the length of chromosomes from chrom_data to simulate SNPs across. 
#' E.g. if chrom 1 is 214066000bp long and you wanted to simulate 8000 SNPs then these would be highly spread apart. But if coverage is 0.001 then the chromosome
#' is only treated as being 214066bp long. 
#' @param file_name (character) File name for bed, bim, fam, map and ped files
#' @param directory (character) Path to store the output files
#' @param parallel (boolean) Enable parallel processing, FALSE by default but HIGHLY recommended to be TRUE if the users computer has multiple cores
#' @param human (boolean) Default TRUE, if so then this uses pre-loaded chrom_data containing the length of chromosomes 1:22 (23 excluded since this is the sex
#' chromosome and can't be simulated by this function). These lengths can be altered by the coverage parameter and determine the scale SNP positions are sampled over
#' @param chrom_data (df) Default NULL, if a dataframe is provided this must contain 'chr' and 'length' columns. If the human parameter is TRUE and chrom_data
#' is NULL then preloaded data is used, this can be overridden if chrom_data is specified
#' @param seed (numeric) Sets the seed for random number generation to make results repeatable. NOTE outputs of the same seed when the parallel parameter differs
#' still differ due to environmental constraints, but within each TRUE/FALSE parallel group the results are identical. 
#' @outputs Creates bed, bim, fam, map and ped files. Returns a list containing three dataframes for allele_info, genotype_info and causal_snps information. 
generate_genomic_data <- function(n_individuals = 200, 
                                  n_snps = 10000, 
                                  chrs = 1:22, 
                                  n_causal_snps = 5,
                                  n_cases = 46,
                                  heritability = 0.8, 
                                  coverage = 1,
                                  file_name = "simulated_data",
                                  directory = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test",
                                  parallel = FALSE,
                                  human = TRUE,
                                  chrom_data = NULL,
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
  setwd(directory)
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
                    59245,53093,44112,37814,38735,34621,45135,25676,29478,20916,28410)) %>%
      mutate(coverage = snp_count/length*100)
  }
  
  if(!is.null(chrom_data)){
    #Create snp matrix based off chromosome size and count
    snp_matrix <- round(n_snps/sum(chrom_data$length)*chrom_data$length)
    chrom_data$snp_count <- snp_matrix
    # Scale data based off the coverage factor
    chrom_data$length <- coverage*chrom_data$length
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
    set.seed(seed+chr)
    print(paste0("Gathering data for chrom: ",chr,"/",n_chrs))
    n_snps <- snp_matrix[chr]
    # Base positions and number on chrom data if availible
    if(!is.null(chrom_data)){
      # max_coverage <- chrom_data$coverage[chr]
      # check_coverage <- n_snps/chrom_data$length[chr]
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
    
    # Calculate LD metric based on distances between SNPs
    distances <- c(0, diff(positions))
    ld_metric <- 0.999^distances
    complement_ld <- 1 - ld_metric
    
    # The probability of each snp being aa,aA,AA is ld * prev_prob + (1-ld) * this_prob
    adjusted_aa <- ld_metric * c(p_aa[1],head(p_aa,-1)) + (1-ld_metric)*p_aa
    adjusted_Aa <- ld_metric * c(p_Aa[1],head(p_Aa,-1)) + (1-ld_metric)*p_Aa
    adjusted_AA <- ld_metric * c(p_AA[1],head(p_AA,-1)) + (1-ld_metric)*p_AA
    
    # Calculate the sums of probabilities for each combination
    sum_aa_Aa <- adjusted_aa + adjusted_Aa
    sum_aa_AA <- adjusted_aa + adjusted_AA
    sum_Aa_AA <- adjusted_Aa + adjusted_AA
    
    # Calculate the proportional adjustments for each genotype combination
    prop_adj_AA_Aa <- adjusted_Aa / sum_Aa_AA
    prop_adj_Aa_AA <- adjusted_AA / sum_Aa_AA
    prop_adj_AA_aa <- adjusted_aa / sum_aa_AA
    prop_adj_aa_AA <- adjusted_AA / sum_aa_AA
    prop_adj_Aa_aa <- adjusted_aa / sum_aa_Aa
    prop_adj_aa_Aa <- adjusted_Aa / sum_aa_Aa
    
    # Create a list of lists. Containing all the baseline probabilities for each SNP. Regarding whether the previous
    # snp was aa,aA or AA (in each case we add ld)
    genotype_prob_bases <- lapply(1:length(adjusted_AA), function(i) {
      list(
        c(ld_metric[i], complement_ld[i] * prop_adj_AA_Aa[i], complement_ld[i] * prop_adj_Aa_AA[i]),
        c(complement_ld[i] * prop_adj_AA_aa[i], ld_metric[i], complement_ld[i] * prop_adj_aa_AA[i]),
        c(complement_ld[i] * prop_adj_Aa_aa[i], complement_ld[i] * prop_adj_aa_Aa[i], ld_metric[i])
      )
    })
    
    # # The probability of each snp given the previous snp matched it is
    # adjusted_aa_ld <- adjusted_aa + ld_metric # CHANGE ME - make weight = ld for one it is and renormalise rest !
    # adjusted_Aa_ld <- adjusted_Aa + ld_metric
    # adjusted_AA_ld <- adjusted_AA + ld_metric
    # 
    # # Create a list of lists. Containing all the baseline probabilities for each SNP. Regarding whether the previous
    # # snp was aa,aA or AA (in each case we add ld)
    # genotype_prob_bases <- lapply(1:length(adjusted_AA), function(i) {
    #   list(
    #     c(adjusted_aa_ld[i], adjusted_Aa[i], adjusted_AA[i]),
    #     c(adjusted_aa[i], adjusted_Aa_ld[i], adjusted_AA[i]),
    #     c(adjusted_aa[i], adjusted_Aa[i], adjusted_AA_ld[i])
    #   )
    # })
    
    
    # Initialise the genotypes_df and alleles_df tables
    first_prob <- c(adjusted_aa[1], adjusted_Aa[1], adjusted_AA[1])
    genotypes_mx <- sample(c(0, 1, 2), n_individuals, prob = first_prob, replace = TRUE)
    genotypes_df <- as.data.frame(genotypes_mx)
    alleles_df <- as.data.frame(allele_vec[[1]][genotypes_mx+1])
    
    # For each snp work out what individual genotypes should be based on the probabilities
    run_end <- NULL
    for(i in 2:n_snps){
      start <- tic()
      print(paste0("Generating Chr",chr," SNP data: ",i,"/",n_snps))
      prev_genotypes <- genotypes_df[, i - 1]
      
      # Get pre-made list of probabilities
      prob_list <- genotype_prob_bases[[i]]
      
      # Call these probabilities based off the previous genotypes
      probs <- do.call(rbind, prob_list[prev_genotypes+1])
      
      # Sample genotypes based on calculated probabilities
      genotypes_sampled <- apply(probs, 1, function(p) sample(c(0, 1, 2), 1, prob = p))
      genotypes_df <- cbind(genotypes_df,genotypes_sampled)
      
      # Map genotype to the allele vector list and store
      allele_sampled <- allele_vec[[i]][genotypes_sampled+1]
      alleles_df <- cbind(alleles_df,allele_sampled)
      
      # Output timing calculation
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
  if (parallel) {
    print(paste0("Using parallel processing (",numCores," cores) to assign genotypes"))
    results <- parLapply(cl, ordered_chrs, process_chromosome)
    stopCluster(cl)
  } else {
    results <- lapply(chrs, process_chromosome)
  }
  
  
  # Work out which SNPs are causal
  print("Assigning causal SNPs")
  allele_infos <- lapply(results, `[[`, "allele_info")
  total_allele_info <- do.call(rbind,allele_infos)
  all_snps <- total_allele_info %>%
    mutate(index_col = row_number())
  
  # Get genotypes in
  genotype_total <- do.call(cbind,lapply(results, `[[`, "genotype_data"))
  
  # Calculate the number of desired controls the user wants
  n_controls <- n_individuals-n_cases
  n_control_prob <- n_controls/n_individuals
  
  
  # # Recalculate the number of cores available in case this has changed
  # if (parallel & ncol(genotype_total) >= 100000) {
  #   cores_avail <- detectCores()
  #   numCores <- cores_avail - 2
  # } else {
  #   numCores <- 1
  # }
  # 
  # # Work out all the weightings for each genotype
  # print(paste0("Calculating probability metrics for causal SNP assignment across all ",nrow(all_snps)," SNPs"))
  # print(paste0("Splitting genotype data into ",numCores," groups to process probability metrics"))
  # split_snp_data <- split.default(genotype_total, rep(1:numCores, each = ceiling(ncol(genotype_total) / numCores))[1:ncol(genotype_total)])
  # investigate_snps <- function(snp_data_input) {
  #   data.frame(sd = apply(snp_data_input, 2, sd)) %>%
  #     mutate(zero_count = colSums(snp_data_input == 0)) %>%
  #     mutate(one_count = colSums(snp_data_input == 1)) %>%
  #     mutate(two_count = colSums(snp_data_input == 2)) %>%
  #     mutate(two_coverage = two_count/n_individuals) %>%
  #     mutate(control_similarity = 1-abs(two_coverage-n_control_prob)) %>%
  #     mutate(geno_total = colSums(snp_data_input)) %>%
  #     rownames_to_column("id")
  # }
  # 
  #
  # # If parallel then run over all cores for speed if data is large
  # if (parallel & ncol(genotype_total) >= 100000) {
  #   cl2 <- makeCluster(numCores)
  #   clusterEvalQ(cl2, library("tidyverse"))
  #   vars_to_export <- list(n_individuals = n_individuals,
  #                          n_control_prob = n_control_prob)
  #   list_names <- names(vars_to_export)
  #   clusterExport(cl2, list_names, envir = list2env(vars_to_export))
  #   
  #   print(paste0("Processing SNP data across ",numCores," cores"))
  #   investigated_snp_data <- parLapply(cl2, split_snp_data, investigate_snps)
  #   stopCluster(cl2)
  # } else {
  #   print("Processing SNP data singularly")
  #   investigated_snp_data <- lapply(split_snp_data, investigate_snps)
  # }

  # After testing non-paralell processing for this stage seems considerably quicker
  investigated_snp_data <- data.frame(sd = apply(genotype_total, 2, sd)) %>%
      mutate(zero_count = colSums(genotype_total == 0)) %>%
      mutate(one_count = colSums(genotype_total == 1)) %>%
      mutate(two_count = colSums(genotype_total == 2)) %>%
      mutate(two_coverage = two_count/n_individuals) %>%
      mutate(control_similarity = 1-abs(two_coverage-n_control_prob)) %>%
      mutate(geno_total = colSums(genotype_total)) %>%
      rownames_to_column("id")
  
  # Collect the data and join to all_snps
  snp_data <- bind_rows(investigated_snp_data) %>%
    left_join(all_snps) %>%
    select(chr, id, sd, zero_count, one_count, two_count, two_coverage, control_similarity, geno_total,
           pos, minor, major, "maf_assigned" = MAF, "maf_observed" = Observed_maf, 
           "ld_toprev" = LD_Prev, "similarity_toprev" = Similarity, index_col)
  
  # Apply some categorical filters
  print(paste0("Applying categorical filter for causal SNPs based on: ",0.05," ≤ MAF ≤ ",0.15," & ",0.2," ≤ LD ≤ ",0.8))
  print("Applying categorical filter for causal SNPs to eliminate SNPs where all genotypes are 2 & where no genotype is 0")
  filt_snps <- snp_data %>%
    # Filter SNPs based on their maf and LD to ensure some structure
    filter(maf_observed >= 0.05,
           # maf_observed <= 0.15,
           # ld_toprev >= 0.2,
           # Causal SNPs must have at least 1 0 and must not be full of 2's
           zero_count != 0,
           two_count != 2*n())
  
  # Send a warning message if we are left with too few SNPs (likely means we are underpowered)
  filt_snp_count <- nrow(filt_snps)
  if(filt_snp_count < n_causal_snps){
    n_causal_snps_adj <- filt_snp_count
    warning(paste0("Filters mean only ",filt_snp_count," potentially causal SNPs remain, less than the ",n_causal_snps," desired"))
  } else {
    n_causal_snps_adj <- n_causal_snps
  }
  
  # Assign weights to causal SNPs
  causal_snp_weights <- filt_snps %>%
    mutate(sd_scale = (sd - min(sd)) / (max(sd) - min(sd))) %>%
    mutate(zero_scale = (zero_count - min(zero_count)) / (max(zero_count) - min(zero_count))) %>%
    mutate(one_scale = (one_count - min(one_count)) / (max(one_count) - min(one_count))) %>%
    mutate(two_scale = (two_count - min(two_count)) / (max(two_count) - min(two_count))) %>%
    mutate(ctrl_sim_scale = (control_similarity - min(control_similarity)) / (max(control_similarity) - min(control_similarity))) %>%
    mutate(geno_scale = (geno_total - min(geno_total)) / (max(geno_total) - min(geno_total))) %>%
    mutate(ctr_adj_sd = ctrl_sim_scale*sd_scale) %>%
    mutate(ctr_adj_sd_scale = (ctr_adj_sd - min(ctr_adj_sd)) / (max(ctr_adj_sd) - min(ctr_adj_sd))) %>%
    group_by(chr) %>%
    arrange(desc(ctr_adj_sd_scale)) %>%
    mutate(chr_rank = 1/row_number()) %>%
    ungroup() %>%
    mutate(causal_weight = ctr_adj_sd_scale*0.8 + chr_rank*0.2) %>%
    mutate(causal_weight_scale = (causal_weight - min(causal_weight)) / (max(causal_weight) - min(causal_weight)))
  
  # Assign n causal SNPs
  causal_snps <- causal_snp_weights
  while(nrow(causal_snps) > n_causal_snps_adj){
    is_causal <- runif(nrow(causal_snps)) < causal_snps$causal_weight_scale
    while(length(which(is_causal)) < n_causal_snps_adj){
      is_causal <- runif(nrow(causal_snps)) < causal_snps$causal_weight_scale
    }
    causal_snps <- causal_snps[is_causal,]
  }
  
  # Create final allele info for later
  allele_info_fnal <- causal_snps %>% 
    right_join(snp_data) %>%
    mutate(causal = (id %in% causal_snps$id)) %>%
    arrange(chr)
  
  print(paste0("Assigning probabilities for each individual being a case based on genotype data (",heritability*100,"%) and environmental factors (",(1-heritability)*100,"%)"))
  # Note the bottom filters out all cases where the genotype score is max, and we square scores to make differnces bigger
  geno_index_data <- genotype_total[,causal_snps$id,drop = FALSE] %>%
    mutate(index_row = row_number())
  geno_pheno_weights <- geno_index_data %>%
    mutate(phenotype_vals = rowSums(.)-index_row) %>%
    mutate(phenotype_scale = (phenotype_vals - min(phenotype_vals)) / (max(phenotype_vals) - min(phenotype_vals))) %>%
    mutate(uniform_weight = 0.5) %>%
    mutate(herit_adj_pheno = heritability * (1-phenotype_scale) + (1 - heritability) * uniform_weight) %>%
    mutate(case_scale = (herit_adj_pheno - min(herit_adj_pheno)) / (max(herit_adj_pheno) - min(herit_adj_pheno)))
  
  # Send an error message if too many 0 cases are seen (means we needed more individuals or fewer cases)
  possible_cases <- n_individuals - sum(geno_pheno_weights$case_scale == 0)
  if(possible_cases < n_cases){
    warning(paste0("There are only ",possible_cases," possible cases, less than the ",n_cases," desired"))
  }
  
  print(paste0("Assigning ",n_cases," cases based on the genotype data of ",n_causal_snps_adj," SNPs across ",n_individuals," individuals"))
  # Assign n cases to individuals
  case_individuals <- geno_pheno_weights
  case_assigned <- data.frame()
  while(nrow(case_assigned) < n_cases){
    is_case <- runif(nrow(case_individuals)) < case_individuals$case_scale
    case_assigned <- rbind(case_assigned,case_individuals[is_case,])
    if(nrow(case_assigned) > n_cases){
      deficit <- nrow(case_assigned) - n_cases
      case_assigned <- case_assigned %>% 
        mutate(remove_weight = 1-(case_scale - min(case_scale)) / (max(case_scale) - min(case_scale))) %>%
        slice_sample(.,n = deficit,weight_by = .$remove_weight) %>%
        anti_join(case_assigned,.)
    }
    case_individuals <- anti_join(case_individuals,case_assigned)
  }
  
  # Apply the phenotype to main data set
  geno_pheno_data <- case_assigned %>%
    right_join(geno_pheno_weights) %>%
    mutate(phenotype = ifelse(index_row %in% case_assigned$index_row,2,1)) %>%
    arrange(index_row)
  
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
  
  # Write out the genomic data information
  genomic_data$allele_info <- allele_info_fnal
  genomic_data$genotype_info <- geno_pheno_data
  genomic_data$causal_snps <- as.data.frame(causal_snps)
  return(genomic_data)
}

## Quick test of results (note these are not default values)
genomic_data <- generate_genomic_data(n_individuals = 2000, 
                                      n_snps = 100000, 
                                      chrs = 1:22,
                                      n_causal_snps = 5,
                                      n_cases = 246,
                                      heritability = 0.8, 
                                      coverage = 0.01,
                                      file_name = "simulated_data",
                                      directory = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test",
                                      parallel = TRUE,
                                      human = TRUE,
                                      chrom_data = NULL,
                                      seed = 123)

# Run a GWAS
system("plink --bfile simulated_data --assoc --ci 0.95 --maf 0.05 --nonfounders --out sim_assoc_results")
sim_results_assoc <- fread("sim_assoc_results.assoc", head=TRUE)

# See what the probabilities of each of the DEFINED causal SNPs being causal are
genomic_data$allele_info %>%
  filter(causal) %>%
  select("SNP" = id) %>%
  inner_join(sim_results_assoc)
# See what the top 5 causal SNPs picked by PLINK are
sim_results_assoc %>%
  arrange(P) %>%
  slice(1:5) 

# Visualise results as a QQ and manhattan plot
qq(sim_results_assoc$P, main = "Q-Q plot of GWAS p-values : log")
manhattan(sim_results_assoc,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc", col = c("lightblue","blue"))




