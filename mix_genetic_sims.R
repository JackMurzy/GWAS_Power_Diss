# NCmisc::list.functions.in.file(filename = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/Fnal_Diss_Functions/fnal_mixed_genetic_sim.R",
#                                alphabetic = TRUE)
library(stats)
library(DirichletReg)
library(parallel)
library(qs)
library(data.table)
library(tidyverse)

## Function to make genetic data, simulate causal SNPs and assign phenotypes
#' Simulate genomic data
#' 
#' @description
#' `simulate_genotype_data_mixed` Creates distinct sets of population genotype data
#' @details
#' This function will generate a list containing multiple populations worth of genotype, allele and SNP information
#' @param n_individuals (numeric) The number of individuals to simulate SNPs across
#' @param n_snps (numeric) The total number of SNPs to simulate across individuals, this is distributed across chromosomes in correspondence to their size to 
#' ensure similar coverage
#' @param chrs (numeric vector) Vector of chromosomes to simulate
#' @param pair_ld (numeric) A value between 0 and 1 (not inclusive) that is raised to the distance between a SNP and the one before it. If the distance is 1bp then
#'  the LD value would be pair_ld^1. A value close to 1 (0.999..) means LD will be near identical, a value of 0 means there is no LD at all.  
#' @param coverage (numeric) A weight between 0 and 1 to scale the length of chromosomes from chrom_data to simulate SNPs across. 
#' E.g. if chrom 1 is 214066000bp long and you wanted to simulate 8000 SNPs then these would be highly spread apart. But if coverage is 0.001 then the chromosome
#' is only treated as being 214066bp long. 
#' @param parallel (boolean) Enable parallel processing, FALSE by default but HIGHLY recommended to be TRUE if the users computer has multiple cores
#' @param human (boolean) Default TRUE, if so then this uses pre-loaded chrom_data containing the length of chromosomes 1:22 (23 excluded since this is the sex
#' chromosome and can't be simulated by this function). These lengths can be altered by the coverage parameter and determine the scale SNP positions are sampled over
#' @param chrom_data (df) Default NULL, if a dataframe is provided this must contain 'chr' and 'length' columns. If the human parameter is TRUE and chrom_data
#' is NULL then preloaded data is used, this can be overridden if chrom_data is specified
#' @param seed (numeric) Sets the seed for random number generation to make results repeatable. NOTE outputs of the same seed when the parallel parameter differs
#' still differ due to environmental constraints, but within each TRUE/FALSE parallel group the results are identical. 
#' @outputs Creates bed, bim, fam, map and ped files suitable for PLINK. Returns a list containing six dataframes for allele data, 
#' genotype data, allele + phenotype data, genotype + phenotype data,  total SNP data, and causal SNP data. 
simulate_genotype_data_mixed <- function(
    n_individuals = 1000,
    n_snps = 100000,
    pair_ld = 0.999,
    coverage = 0.01,
    chrs = 1:22,
    n_populations = 2,
    human = TRUE,
    chrom_data = NULL,
    parallel = TRUE,
    seed = 123,
    temp_dir = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/parallel_tmp_dir"){
  
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
    scaled_chrom_length <- coverage*chrom_data$length
    if(any(snp_matrix > scaled_chrom_length)){
      warning(paste0("The number of SNPs per chromosome exceeds the adjusted (by coverage) chromosome length,
                     chrom length is therefore being adjusted to match the desried number of SNPs"))
      extend_length <- which(snp_matrix > scaled_chrom_length)
      scaled_chrom_length[extend_length] <- snp_matrix[extend_length]
    }
    chrom_data$length <- scaled_chrom_length
  } else {
    # Set SNP matrix as equal across all chromosomes
    snp_matrix <- rep(n_snps,n_chrs)
  }
  
  # If parallel in use sequester cores
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
                           n_populations = n_populations,
                           pair_ld = pair_ld,
                           seed = seed,
                           temp_dir = temp_dir,
                           tic = tic, 
                           toc = toc,
                           time_format = time_format)
    list_names <- names(vars_to_export)
    clusterExport(cl, list_names, envir = list2env(vars_to_export))
  }
  
  # Function to process each chromosome
  process_chromosome <- function(chr) {
    # Set the temporary directory
    Sys.setenv(TMPDIR = temp_dir)
    
    if(!is.null(seed)){
      set.seed(seed+chr)
    }
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
    
    # Calculate LD metric based on distances between SNPs
    distances <- c(0, diff(positions))
    ld_metric <- pair_ld^distances
    complement_ld <- 1 - ld_metric
    ld_copy <- round(n_individuals * ld_metric)
    
    # Generate allele data
    bases <- c("A","T","C","G")
    
    # For each population produce a set of genomic data
    population_data <- list()
    for(pop in 1:n_populations){
      set.seed(seed+chr+pop)
      minor_alleles <- sample(bases,n_snps,replace = TRUE)
      major_alleles <- sapply(minor_alleles, function(x) sample(bases[bases != x], 1))
      # Calculate allele vectors
      allele_vec <- mapply(function(min,maj){c(paste0(min, min), paste0(min, maj), paste0(maj, maj))},
                           minor_alleles, major_alleles, SIMPLIFY = FALSE)
      
      # Assign MAF based off a beta distribution (ensuring we have a large enough pool to filter from)
      maf_initial <- rbeta(n_snps*100, 1, 15) # Changed from 1,10 !
      maf_filtered <- maf_initial[maf_initial <= 0.5]
      maf <- sample(maf_filtered, n_snps)
      
      
      # Generate HWE vectors
      p <- maf
      q <- 1 - p
      p_aa <- p^2
      p_Aa <- 2 * p * q
      p_AA <- q^2
      
      # The probability of each snp being aa,aA,AA is ld * prev_prob + (1-ld) * this_prob
      adjusted_aa <- ld_metric * c(p_aa[1],head(p_aa,-1)) + complement_ld*p_aa
      adjusted_Aa <- ld_metric * c(p_Aa[1],head(p_Aa,-1)) + complement_ld*p_Aa
      adjusted_AA <- ld_metric * c(p_AA[1],head(p_AA,-1)) + complement_ld*p_AA
      
      # Create a probability list
      adjusted_prob_list <- lapply(1:length(adjusted_AA), function(i) {
        c(adjusted_aa[i],adjusted_Aa[i],adjusted_AA[i])
      })
      
      # When LD is high then the previous SNP should inform the next
      # When LD is low then the maf should inform the previous SNP
      # Therefore take LD as a percentage, 1 = 100% and 0 = 0%, fill in a random amount of this percentage and sample the rest
      # from the MAF
      
      # Initialise the genotypes_df and alleles_df tables
      first_prob <- adjusted_prob_list[[1]]
      genotypes_mx <- sample(c(2, 1, 0), n_individuals, prob = first_prob, replace = TRUE)
      genotypes_df <- as.data.frame(genotypes_mx)
      alleles_df <- as.data.frame(allele_vec[[1]][4-(genotypes_mx+1)])
      
      # For each snp work out what individual genotypes should be based on the probabilities and calculate metrics
      similarity_scores <- numeric(n_snps)
      observed_maf <- similarity_scores
      zero_count <- observed_maf
      one_count <- zero_count
      two_count <- one_count
      geno_total <- two_count
      similarity_scores[1] <- 1
      run_end <- NULL
      for(i in 1:n_snps){
        start <- tic()
        
        # Generate new data
        if(i > 1){
          print(paste0("Generating Chr",chr," SNP data: ",i,"/",n_snps))
          prev_genotypes <- genotypes_df[, i - 1]
          
          # Create a vector of indices to copy
          indices_to_copy <- sample(1:n_individuals, ld_copy[i])
          
          # Create a new vector and copy the specified genotypes
          new_genotypes <- prev_genotypes
          
          # Fill the remaining positions with a random sample of 0, 1, and 2
          remaining_indices <- setdiff(1:n_individuals, indices_to_copy)
          new_genotypes[remaining_indices] <- sample(2:0,length(remaining_indices),prob = adjusted_prob_list[[i]], replace = TRUE)
          genotypes_df <- cbind(genotypes_df,new_genotypes)
          
          # Map genotype to the allele vector list and store
          allele_sampled <- allele_vec[[i]][4-(new_genotypes+1)]
          alleles_df <- cbind(alleles_df,allele_sampled)
          
          # Calculate similarity score
          similarity <- sum(new_genotypes == prev_genotypes) / n_individuals
          similarity_scores[i] <- similarity
        } else {
          new_genotypes <- genotypes_df
        }
        
        # Calculate count scores
        zero_count[i] <- sum(new_genotypes == 0)
        one_count[i] <- sum(new_genotypes == 1)
        two_count[i] <- n_individuals - (zero_count[i] + one_count[i])
        geno_total[i] <- one_count[i] + 2*two_count[i]
        # Calculate the MAF
        observed_maf[i] <- (geno_total[i]) / (2 * n_individuals)
        
        # Output timing calculation
        cat("\014")
        print(paste0("Completed: ",i/n_snps*100,"%"))
        run_end <- rbind(run_end,toc(start))
        mean(run_end)
        print(paste0("Estimated time remaining: ",time_format((n_snps-i)*mean(run_end))))
      }
      colnames(genotypes_df) <- snp_ids
      colnames(alleles_df) <- snp_ids
      
      # Add in metrics
      print("Merging results")
      # Merge data into a useful data to be exported frame
      allele_info_df <- data.frame(chr = chr,id = snp_ids,pos = positions,minor = minor_alleles, major = major_alleles,
                                   maf_assigned = maf,ld_prev = ld_metric,similarity = similarity_scores,maf_observed = observed_maf,
                                   zero_count = zero_count, one_count = one_count, two_count = two_count, geno_total = geno_total)
      
      #Cleanup temp files after finished 
      temp_files <- list.files(temp_dir,full.names = TRUE)
      unlink(temp_files)
      
      # Return all useful values
      population_data[[paste0("pop_",pop)]] <- list(allele_info = allele_info_df, genotype_data = genotypes_df, allele_data = alleles_df)
    }
    return(population_data)  # Return a list of data
  }
  
  # Use parallel or sequential processing
  if (parallel) {
    print(paste0("Using parallel processing (",numCores," cores) to assign genotypes"))
    results <- parLapply(cl, ordered_chrs, process_chromosome)
    stopCluster(cl)
  } else {
    results <- lapply(chrs, process_chromosome)
  }
  set.seed(seed)
  
  genotype_data <- list()
  for(pop in 1:n_populations){
    print(paste0("Reading in results for population ",pop,"/",n_populations))
    population_results <- lapply(results, `[[`, paste0("pop_",pop))
    
    # Get allele info
    allele_infos <- lapply(population_results, `[[`, "allele_info")
    total_allele_info <- do.call(rbind,allele_infos)
    snp_info <- total_allele_info %>%
      mutate(index_col = row_number()) %>%
      arrange(chr,pos)
    # Get genotype info
    genotype_info <- do.call(cbind,lapply(population_results, `[[`, "genotype_data"))
    arrange_genotype <- genotype_info[,snp_info$index_col]
    # Get allele data (corresponds to genotype data)
    allele_data <- do.call(cbind,lapply(population_results, `[[`, "allele_data"))
    arrange_allele <- allele_data[,snp_info$index_col]
    
    # Finally change index_col to reflect new order
    arrange_snp <- snp_info %>%
      mutate(index_col = row_number())
    
    genotype_data[[paste0("pop_",pop)]] <- list(allele_info = arrange_snp, 
                                                genotype_info = arrange_genotype, 
                                                allele_data = arrange_allele)
    
  }
  
  print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
  
  return(genotype_data)
}



## Function to combine multiple population data using a haplotype resampling approach
#' Simulate genomic data
#' 
#' @description
#' `combine_populations` Mixes multiple populations together by imputing and resampling their haplotypes into a single new population
#' @details
#' This function generates genomic data for a single mixed population based on N populations. This uses a resampling approach to
#' piece together haplotype blocks in a biologically informed manner
#' @param pop_genotype_data (list) The output of the `simulate_genotype_data_mixed` function
#' @param n_individuals (numeric) The number of individuals to simulate SNPs across
#' @param n_chrs (numeric vector) Number of chromosomes being simulated
#' @param admixture_weights (numeric vector) Mixing weights for each population in the same order they are listed in within pop_genotype_data 
#' @param file_name (character) File name for all output files
#' @param directory (character) Path to store the output files
#' @param parallel (boolean) Enable parallel processing, FALSE by default but HIGHLY recommended to be TRUE if the users computer has multiple cores
#' @param seed (numeric) Sets the seed for random number generation to make results repeatable. NOTE outputs of the same seed when the parallel parameter differs
#' still differ due to environmental constraints, but within each TRUE/FALSE parallel group the results are identical. 
#' @outputs List of four dataframes containing the original basic allele info for each population, the new mixed populations genotype and allele info
#' as well as a dataframe of each individuals admixture ratio between each population
combine_populations <- function(pop_genotype_data = genotype_data,
                                n_individuals = 1000,
                                n_chrs = 22,
                                admixture_weights = c(1,1),
                                seed = 123,
                                parallel = TRUE,
                                directory = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_tests_50_50",
                                file_name = "mixed_data_50_50"){
  setwd(directory)
  split_group_allele_info <- list()
  for(pop in names(pop_genotype_data)){
    spec_pop_genotype <- pop_genotype_data[[pop]]
    split_group_allele_info[[pop]] <- spec_pop_genotype$allele_info %>%
      group_by(chr) %>%
      mutate(entries_per_chr = n(),
             threshold = ceiling(entries_per_chr / 30)) %>%
      mutate(split_point = case_when(similarity < mean(similarity) ~ TRUE, TRUE ~ FALSE)) %>%
      mutate(group = {
        group_id <- 1
        counter <- 0
        groups <- vector("integer", n())
        for (i in seq_along(similarity)) {
          counter <- counter + 1
          if (split_point[i] && counter >= threshold[1]) {
            group_id <- group_id + 1
            counter <- 0
          }
          groups[i] <- group_id
        }
        groups
      }) %>%
      ungroup() %>%
      mutate(population = pop) %>%
      group_by(chr,group) %>%
      mutate(global_group = cur_group_id()) %>%
      ungroup() %>%
      mutate(lead_id = lead(id))
  }
  
  set.seed(seed)
  individ_mix_mx <- rdirichlet(n_individuals,admixture_weights)
  individ_mix <- individ_mix_mx %>%
    as.data.frame() %>%
    mutate(indvidual = row_number())
  pop_names <- names(split_group_allele_info)
  colnames(individ_mix) <- c(pop_names,"individual")
  indvid_probs <- individ_mix %>%
    pivot_longer(cols = pop_names, names_to = "population", values_to = "prob_pick")
  total_data <- do.call(rbind,split_group_allele_info)
  n_snps <- nrow(split_group_allele_info[[1]])
  
  process_individual_group <- function(indv_groups) {
    ordered_popgen <- list()
    ordered_popgallele <- list()
    admixture_check <- list()
    picked_data_ordered <- list()
    for (indv in indv_groups) {
      cat("Processing individual", indv, "in cluster", Sys.getpid(), "\n")
      set.seed(seed)
      total_data_remain <- total_data %>%
        arrange(chr, pos, population) %>%
        mutate(individual = indv) %>%
        left_join(indvid_probs) %>%
        sample_frac()
      picked_data <- NULL
      
      while (nrow(total_data_remain) != 0) {
        cat("Remaining rows:", nrow(total_data_remain), "\n")
        row_pick <- ifelse(nrow(total_data_remain) < 1000, nrow(total_data_remain), 1000)
        pick_frame <- total_data_remain %>%
          slice_sample(n = row_pick, weight_by = prob_pick)
        
        selected_groups <- total_data_remain %>%
          semi_join(pick_frame, by = c("population", "global_group"))
        
        duplicate_ids <- selected_groups %>%
          group_by(id) %>%
          filter(n() > 1) %>%
          ungroup()
        
        non_duplicate_ids <- selected_groups %>%
          anti_join(duplicate_ids)
        
        while (nrow(duplicate_ids) != 0) {
          selected_dupe <- duplicate_ids %>%
            semi_join(duplicate_ids[1,], by = c("population", "global_group"))
          non_duplicate_ids <- bind_rows(non_duplicate_ids, selected_dupe)
          duplicate_ids <- duplicate_ids %>%
            filter(!id %in% selected_dupe$id)
        }
        
        picked_data <- bind_rows(picked_data, non_duplicate_ids)
        total_data_remain <- anti_join(total_data_remain, non_duplicate_ids, by = "id")
      }
      
      picked_data_ordered[[indv]] <- picked_data %>%
        arrange(chr, pos)
      rm(picked_data)
      gc() # Clean up RAM
    }
    setwd(directory)
    qsave(picked_data_ordered, file = paste0(file_name,"_result_chunk_",indv_groups[1], ".qs"))
    return(indv_groups)
  }
  
  # Convert pop_genotype to something much shorter
  extract_needed_data <- function(pop_list) {
    pop_list[c("genotype_info", "allele_data")]
  }
  pop_genotype_dif_data <- lapply(pop_genotype_data, extract_needed_data)
  gc()
  split_into_groups <- function(vec, n) {
    groups <- cut(seq_along(vec), breaks = n, labels = FALSE)
    split(vec, groups)
  }
  
  total_start <- tic()
  if (parallel) {
    cores_avail <- detectCores()
    numCores <- cores_avail - 1  # Reserve one core for system operations
    
    indv_groups <- split_into_groups(c(1:n_individuals),(numCores*2))
    cl <- makeCluster(numCores)
    # Export necessary packages to each cluster
    clusterEvalQ(cl, {
      library(dplyr)
      library(qs)
    })
    vars_to_export <- list(seed = seed,
                           total_data = total_data, 
                           indvid_probs = indvid_probs,
                           individ_mix = individ_mix,
                           process_individual_group = process_individual_group,
                           directory = directory,
                           file_name = file_name)
    list_names <- names(vars_to_export)
    clusterExport(cl, list_names, envir = list2env(vars_to_export))
    print(paste("Using", numCores, "cores for parallel execution"))
    numRounds <- ceiling(length(indv_groups) / numCores)
    results <- list()
    for (i in 1:numRounds) {
      start_index <- ((i - 1) * numCores) + 1
      end_index <- min(i * numCores, length(indv_groups))
      indv_groups_chunk <- indv_groups[start_index:end_index]
      results[[i]] <- parLapply(cl, indv_groups_chunk, process_individual_group)
      gc()
    }
    stopCluster(cl)
  } else {
    result <- lapply(indv_groups, process_individual)
  }
  print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
  setwd(directory)
  # Save incase function runs out of RAM
  qsave(pop_genotype_data, paste0(file_name,"_pop_genotype_data.qs"))
  # Reload back in 
  pop_genotype_data <- qread(paste0(file_name,"_pop_genotype_data.qs"))
  # Convert pop_genotype to something much shorter
  extract_needed_data <- function(pop_list) {
    pop_list[c("genotype_info", "allele_data")]
  }
  pop_genotype_dif_data <- lapply(pop_genotype_data, extract_needed_data)
  
  
  total_start <- tic()
  group_ordered_popgen <- list()
  group_ordered_popgallele <- list()
  group_admixture_check <- list()
  for (indv_group_use in indv_groups) {
    print("Reading in qs file")
    indv_pop_data <- qread(paste0(directory,"/",file_name,"_result_chunk_",indv_group_use[1],".qs"))
    ordered_popgen <- list()
    ordered_popgallele <- list()
    admixture_check <- list()
    for(individual_in_use in indv_group_use){
      individual_data <- indv_pop_data[[individual_in_use]]
      print(paste0("Getting results for individual ",individual_in_use,"/",n_individuals))
      popgen_data <- list()
      popallele_data <- list()
      for (pop_use in unique(individual_data$population)) {
        pop_data <- individual_data %>% filter(population == pop_use)
        popgen_data[[pop_use]] <- pop_genotype_dif_data[[pop_use]]$genotype_info[individual_in_use, pop_data$index_col]
        popallele_data[[pop_use]] <- pop_genotype_dif_data[[pop_use]]$allele_data[individual_in_use, pop_data$index_col]
      }
      names(popgen_data) <- NULL
      names(popallele_data) <- NULL
      ordered_popgen[[individual_in_use]] <- do.call(cbind, popgen_data) %>%
        select(individual_data$id)
      ordered_popgallele[[individual_in_use]] <- do.call(cbind, popallele_data) %>%
        select(individual_data$id)
      admixture_check[[individual_in_use]] <- individual_data %>%
        group_by(population) %>%
        summarise(count = n()) %>%
        mutate(percentage = (count / sum(count)) * 100,
               individual = individual_in_use,
               pop_weights = {
                 unique_pops <- unique(individual_data$population)
                 individ_mix %>% select(all_of(unique_pops)) %>% slice(individual_in_use) %>% as.matrix() %>% c()
               })
      
    }
    group_ordered_popgen[[indv_group_use[1]]] <- do.call(rbind, ordered_popgen)
    group_ordered_popgallele[[indv_group_use[1]]] <- do.call(rbind, ordered_popgallele)
    group_admixture_check[[indv_group_use[1]]] <- do.call(rbind,admixture_check)
    rm(indv_pop_data)
  }
  total_ordered_popgen <- do.call(rbind, group_ordered_popgen)
  total_ordered_popallele <- do.call(rbind,group_ordered_popgallele)
  total_admixture_check <- do.call(rbind,group_admixture_check)
  print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
  
  
  extract_allele_data <- function(pop_list) {
    pop_list[c("allele_info")]
  }
  pop_genotype_allele_data <- lapply(pop_genotype_data, extract_allele_data)
  total_mixed_data <- list(
    pop_alleles = pop_genotype_allele_data,
    total_popgen = total_ordered_popgen,
    total_popallele = total_ordered_popallele,
    total_admixture = total_admixture_check
  ) 
  setwd(directory)
  qsave(total_mixed_data, file = paste0(file_name,"total_mixed_results.qs"))
  return(total_mixed_data)
}

