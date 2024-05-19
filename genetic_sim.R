# NCmisc::list.functions.in.file(filename = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/Sim2.R",
#                                                                alphabetic = TRUE)
library(parallel)
library(data.table)
library(stats)
library(tidyverse)

simulate_genotype_data <- function(
    n_individuals = 200,
    n_snps = 10000,
    pair_ld = 0.999,
    coverage = 1,
    chrs = 1:22,
    human = TRUE,
    chrom_data = NULL,
    parallel = TRUE,
    seed = NULL,
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
    
    # Generate allele data
    bases <- c("A","T","C","G")
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
    
    
    # Calculate LD metric based on distances between SNPs
    distances <- c(0, diff(positions))
    ld_metric <- pair_ld^distances
    complement_ld <- 1 - ld_metric
    
    # Calcualte LD as a copy percentage
    ld_copy <- round(n_individuals * ld_metric)
    
    
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
    return(list(allele_info = allele_info_df, genotype_data = genotypes_df, allele_data = alleles_df))  # Return a list of data
  }
  
  # Use parallel or sequential processing
  if (parallel) {
    print(paste0("Using parallel processing (",numCores," cores) to assign genotypes"))
    results <- parLapply(cl, ordered_chrs, process_chromosome)
    stopCluster(cl)
  } else {
    results <- lapply(chrs, process_chromosome)
  }
  
  # Get SNP info
  print("Reading in results")
  allele_infos <- lapply(results, `[[`, "allele_info")
  total_allele_info <- do.call(rbind,allele_infos)
  snp_info <- total_allele_info %>%
    mutate(index_col = row_number()) %>%
    arrange(chr,pos)
  # Get genotype info
  genotype_info <- do.call(cbind,lapply(results, `[[`, "genotype_data"))
  arrange_genotype <- genotype_info[,snp_info$index_col]
  
  # Get allele data (corresponds to genotype data)
  allele_data <- do.call(cbind,lapply(results, `[[`, "allele_data"))
  arrange_allele <- allele_data[,snp_info$index_col]
  
  # Finally change index_col to reflect new order
  arrange_snp <- snp_info %>%
    mutate(index_col = row_number())
  
  print(paste0("Done ! -- Total time elapsed: ",time_format(toc(total_start))))
  genotype_data <- list(allele_info = arrange_snp, genotype_info = arrange_genotype, allele_data = arrange_allele)
  return(genotype_data)
}

causal_snp_assignment <- function(
    genotype_data = NULL,
    n_causal_snps = 5,
    causal_maf = 0.01,
    seed = NULL){
  
  set.seed(seed)
  
  # If the names are not expected then return an error and stop the function
  if(!any(names(genotype_data) %in% c("allele_info","genotype_info","allele_data"))){
    stop("genotype_data must be a list of dataframes called allele_info, genotype_info and allele_data")
  }
  
  # Load in the data
  allele_info <- genotype_data$allele_info
  genotype_info <- genotype_data$genotype_info
  
  # Filter for causal variants based on desired MAF (also filter out any where everything is a 0)
  causal_weights <- allele_info %>%
    filter(zero_count != nrow(genotype_info)) %>% 
    filter(maf_observed >= causal_maf) %>%
    mutate(maf_diff = abs(maf_observed - causal_maf)) %>%
    mutate(maf_sim_norm = 1-(maf_diff - min(maf_diff)) / (max(maf_diff) - min(maf_diff)))
  
  causal_snp <- NULL
  causal_reweight <- causal_weights
  for(i in 1:n_causal_snps){
    one_weighted <- causal_reweight[causal_reweight$maf_sim_norm == 1,]
    selected_snp <- one_weighted[sample(1:nrow(one_weighted), 1), ]
    causal_snp <- rbind(causal_snp, selected_snp)
    
    # Remove causal SNP +- 1000bp from remaining data
    causal_reweight <- causal_reweight %>%
      filter(!(chr == selected_snp$chr & 
                 pos >= selected_snp$pos - 2000 & 
                 pos <= selected_snp$pos + 2000)) %>%
      mutate(maf_sim_norm = (maf_sim_norm - min(maf_sim_norm)) / (max(maf_sim_norm) - min(maf_sim_norm)))
  }
  
  return(causal_snp)
}


phenotype_assigment <- function(
    genotype_data = NULL,
    causal_snps = NULL,
    n_cases = 50,
    genef_size = 1,
    envef_size = 0,
    heritability = 0.8,
    override_ncases = FALSE,
    seed = NULL){
  
  set.seed(seed)
  
  # If the names are not expected then return an error and stop the function
  if(!any(names(genotype_data) %in% c("allele_info","genotype_info"))){
    stop("genotype_data must be a list of dataframes called allele_info and genotype_info")
  }
  if(!any(names(causal_snps) %in% c('chr','id','pos','minor','major','maf_assigned','ld_prev','similarity','maf_observed','index_col','maf_diff','maf_sim_norm'))){
    stop("causal_snps must be a dataframe with columns: 'chr','id','pos','minor','major','maf_assigned','ld_prev','similarity','maf_observed','index_col','maf_diff' and 'maf_sim_norm'")
  }
  
  # If effect size is not the correct length vector stop
  if(length(genef_size) == 1){
    warning("genef_size is being assumed equal for each causal_snp")
    genef_size <- rep(genef_size,nrow(causal_snps))
  } else if(length(genef_size) != nrow(causal_snps)){
    stop("genef_size must be a vector of values which correspond to each causal_snp")
  }
  
  # Load in the data
  allele_info <- genotype_data$allele_info
  genotype_info <- genotype_data$genotype_info
  allele_data <- genotype_data$allele_data
  n_individuals <- nrow(genotype_info)
  
  # Get the genetic risk score
  causal_snp_genotype <- genotype_info[, causal_snps$index_col, drop = FALSE]
  causal_snp_alleles <- allele_data[, causal_snps$index_col, drop = FALSE]
  genetic_score <- rowSums(as.matrix(causal_snp_genotype) %*% diag(genef_size))
  
  # Normalise risk score
  if(length(unique(genetic_score)) != 1){
    normalised_gen <- (genetic_score - min(genetic_score)) / (max(genetic_score) - min(genetic_score))
  } else {
    normalised_gen <- genetic_score
  }
  
  # Get the environmental risk score
  # If effect size is not the correct length vector stop
  if(length(envef_size) == 1){
    warning("envef_size is being assumed equal for each causal_snp")
    envef_size <- rep(envef_size,n_individuals)
  } else if(length(envef_size) != n_individuals){
    stop("envef_size must be a vector of values which correspond to each individual")
  }
  if(length(unique(envef_size)) != 1){
    normalised_envo <- (envef_size - min(envef_size)) / (max(envef_size) - min(envef_size))
  } else{
    normalised_envo <- envef_size
  }
  
  if(length(heritability) == 1){
    warning("heritability is being assumed equal for each causal_snp")
    heritability <- rep(heritability,n_individuals)
  } else if(length(heritability) != n_individuals){
    stop("heritability must be a vector of values which correspond to each individual")
  }
  if(length(unique(heritability)) != 1){
    normalised_herit <- (heritability - min(heritability)) / (max(heritability) - min(heritability))
  } else{
    normalised_herit <- heritability
  }
  
  # Get the final effects
  final_effect <- normalised_gen*normalised_herit + (1-normalised_herit)*normalised_envo
  
  # Check there aren't too many 0's
  possible_cases <- sum(final_effect != 0)
  if(possible_cases < n_cases){
    warning(paste0("Only ",possible_cases," individuals could be a case, the rest have a 0% chance as stands. This usually occurs when MAF, n_individuals and heritibility is too low"))
    if(override_ncases){
      warning(paste0("Overriding the number of cases required, returning ",possible_cases," cases"))
      n_cases <- possible_cases
    } else {
      warning(paste0("Adding in small environmental effect to ensure ",n_cases," can be returned. NOTE ",n_cases - possible_cases," will are effectively 'random' cases."))
      final_effect[final_effect==0] <- 1e-5
    }
  }
  
  # Assign phenotypes
  pheno <- rep(1,n_individuals)
  pheno[sample(1:n_individuals,n_cases,prob = final_effect)] <- 2
  phenotypic_data <- data.frame(phenotype = pheno, totef = final_effect, genetic_score = genetic_score, genef = normalised_gen, environment_score = envef_size, envef = normalised_envo)
  geno_pheno_data <- cbind(causal_snp_genotype,phenotypic_data)
  allele_pheno_data <- cbind(causal_snp_alleles,phenotypic_data)
  
  # Work out a metric for each causal SNP
  case_control_stats <- geno_pheno_data %>%
    arrange(phenotype) %>%
    group_by(phenotype) %>%
    summarise(across(starts_with("rs"), sum)) %>%
    pivot_longer(cols = -phenotype, names_to = "id", values_to = "sum") %>%
    pivot_wider(names_from = phenotype, values_from = sum, names_prefix = "sum_") %>%
    rename(sum_control = sum_1, sum_case = sum_2) %>%
    mutate(case_control_ratio = sum_case/sum_control)
  
  
  # Return data
  return(list(geno_pheno_data = geno_pheno_data,allele_pheno_data = allele_pheno_data,case_control_stats = case_control_stats))
}

## Function to make genetic data
genetic_sim <- function(
    n_individuals = 2000,
    n_snps = 10000,
    chrs = 1:22,
    n_causal_snps = 5,
    causal_maf = 0.05,
    n_cases = 50,
    genef_size = 1,
    envef_size = 0,
    heritability = 0.8,
    pair_ld = 0.999,
    coverage = 0.01,
    write_plink = TRUE,
    file_name = "simulated_data",
    directory = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test",
    temp_dir = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/parallel_tmp_dir",
    human = TRUE,
    chrom_data = NULL,
    override_ncases = FALSE,
    parallel = TRUE,
    seed = NULL){
  
  # Load genotype data
  print("Generating genotype data ")
  genotype_data <- simulate_genotype_data(
    n_individuals = n_individuals,
    n_snps = n_snps,
    pair_ld = pair_ld,
    coverage = coverage,
    chrs = chrs,
    human = human,
    chrom_data = chrom_data,
    parallel = parallel,
    seed = seed,
    temp_dir = temp_dir
  )
  
  print("Assigning causal SNPs")
  causal_snps <- causal_snp_assignment(
    genotype_data = genotype_data,
    n_causal_snps = n_causal_snps,
    causal_maf = causal_maf,
    seed = seed
  )
  
  print("Assigning phenotypes")
  phenotype_indv <- phenotype_assigment(
    genotype_data = genotype_data,
    causal_snps = causal_snps,
    n_cases = n_cases,
    genef_size = genef_size,
    envef_size = envef_size,
    heritability = heritability,
    override_ncases = override_ncases,
    seed = seed
  )
  
  print("Compiling data nicely")
  total_allele_info <- genotype_data$allele_info %>%
    mutate(causal = ifelse(id %in% causal_snps$id,T,F))
  index_cols <- causal_snps$index_col
  causal_snps_better <- causal_snps %>% 
    select(-index_col) %>%
    left_join(phenotype_indv$case_control_stats) %>%
    mutate(index_col = index_cols)

  
  genetic_data <- list(causal_snp_data = causal_snps_better,
                       causal_geno_pheno_data = phenotype_indv$geno_pheno_data,
                       causal_allele_pheno_data = phenotype_indv$allele_pheno_data,
                       all_snp_data = total_allele_info,
                       all_genotype_data = genotype_data$genotype_info,
                       all_allele_data = genotype_data$allele_data)
  
  
  # Format files for PLINK
  if(write_plink){
    setwd(directory)
    # Add PLINK to path if not already there
    plink_path <- "~/Desktop/University/Biomedicine/y3/Dissertation/plink_mac_20231211"
    if (!grepl(plink_path, Sys.getenv("PATH"))) {
      Sys.setenv(PATH = paste(Sys.getenv("PATH"), plink_path, sep = ":"))
    }
    print("Creating PLINK files")
    print("Making map file")
    map_data <- data.frame(CHR = total_allele_info$chr, SNP = total_allele_info$id, GD = 0, BPP = total_allele_info$pos)
    fwrite(map_data, file = paste0(file_name, ".map"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    print("Making ped file")
    ped_easy <- data.frame(FID = 1:n_individuals, IID = 1:n_individuals, PID = 0, MID = 0, Sex = sample(1:2,n_individuals,replace = TRUE))
    phenotype <- phenotype_indv$geno_pheno_data$phenotype
    ped_data <- cbind(ped_easy,phenotype,genetic_data$all_allele_data)
    fwrite(ped_data, file = paste0(file_name, ".ped"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Write out the final bfile
    print(paste0("Making final bfile with name: '",file_name,"'"))
    system(paste0("plink --file ",file_name," --make-bed --out ",file_name))
  }
  print("Finished !")
  return(genetic_data)
}

