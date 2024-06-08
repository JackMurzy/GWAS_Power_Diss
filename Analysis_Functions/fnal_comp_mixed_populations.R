# Compare the mixed popualtions to each other
library(qs)
library(ggsignif)
library(ggnewscale)
library(tidyverse)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_sim.R")
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_assoc_stats_functions.R")

setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_tests_50_50")
mixed_sim_data_50_50 <- qread("total_mixed_results_50_50.qs")
setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_tests")
mixed_sim_data_100_10 <- qread("total_mixed_results.qs")

indv_pop_genotype_data  <- simulate_genotype_data_mixed(
  n_individuals = 1000,
  n_snps = 100000,
  pair_ld = 0.999,
  coverage = 0.01,
  chrs = 1:22,
  n_populations = 2,
  human = TRUE,
  chrom_data = NULL,
  parallel = TRUE,
  seed = 123)



genotype_data <- list()
total_mixed_data <- list(mixed_sim_data_50_50,mixed_sim_data_100_10)
for(i in 1:length(total_mixed_data)){
  mixed_sim_data <- total_mixed_data[[i]]
  admix <- mixed_sim_data$total_admixture
  allele_info_1 <- mixed_sim_data$pop_genotypes$pop_1$allele_info
  allele_info_2 <- mixed_sim_data$pop_genotypes$pop_2$allele_info
  geneotype_info <- mixed_sim_data$total_popgen
  allele_data <- mixed_sim_data$total_popallele
  zero_count <- NULL
  one_count <- zero_count
  two_count <- one_count
  geno_total <- two_count
  observed_maf <- geno_total
  n_individuals <- nrow(geneotype_info)
  for(j in 1: ncol(geneotype_info)){
    print(j)
    new_genotypes <- geneotype_info[,j]
    # Calculate count scores
    zero_count[j] <- sum(new_genotypes == 0)
    one_count[j] <- sum(new_genotypes == 1)
    two_count[j] <- n_individuals - (zero_count[j] + one_count[j])
    geno_total[j] <- one_count[i] + 2*two_count[j]
    # Calculate the MAF
    observed_maf[j] <- (geno_total[j]) / (2 * n_individuals)
  }
  allele_info <- allele_info_1
  allele_info$zero_count <- zero_count
  allele_info$one_count <- one_count
  allele_info$two_count <- two_count
  allele_info$geno_total <- geno_total
  allele_info$maf_observed <- observed_maf
  genotype_data[[i]] <- list(
    allele_info = allele_info,
    genotype_info = geneotype_info,
    allele_data = allele_data
  )
}

# Add in the single population data
genotype_data[[3]] <- indv_pop_genotype_data$pop_1
genotype_data[[4]] <- indv_pop_genotype_data$pop_2

# get all allele info filtered for just the potentially causal SNPS
pot_causal_snp <- list()
for(i in 1:length(genotype_data)){
  pot_causal_snp[[i]] <- genotype_data[[i]]$allele_info %>%
    filter(zero_count != nrow(genotype_data[[i]]$genotype_info)) %>%
    filter(maf_observed >= 0.05) %>%
    filter(maf_observed <= 0.1) %>%
    select(chr,id,pos,index_col)
}

selected_snps <- inner_join(pot_causal_snp[[1]],pot_causal_snp[[2]])

total_geno <- rbind(genotype_data[[1]]$genotype_info,genotype_data[[2]]$genotype_info,
                    genotype_data[[3]]$genotype_info,genotype_data[[4]]$genotype_info)
total_geno_mx <- as.matrix(total_geno)
library(bigstatsr)
library(bigreadr)

big_geno <- as_FBM(total_geno_mx)
pca_result <- big_SVD(big_geno, fun.scaling = big_scale(center = TRUE, scale = TRUE), k = 10)
pca_scores <- data.table(PC1 = pca_result$u[, 1], PC2 = pca_result$u[, 2])

pca_scores$Cohort <- c(rep("AB1",1000),rep("AB2",1000),rep("A",1000),rep("B",1000))
pca_scores$c_number <- c(rep("3",1000),rep("4",1000),rep("1",1000),rep("2",1000))
pca_scores_use <- pca_scores %>%
  arrange(c_number)
pca_scores_use$Cohort <- factor(pca_scores_use$Cohort, levels = unique(pca_scores_use$Cohort))
# Plot the PCA analysis
pca_analysis_results <- ggplot(pca_scores_use, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(alpha = 0.7) +
  labs(x = "Principal Component 1", y = "Principal Component 2",
       title = "PCA of Population Genotype Data",
       color = "Population") +
  scale_color_manual(values = c("A" = "blue", "B" = "red", "AB1" = "green", "AB2" = "purple")) +
  theme_minimal() + 
  theme(
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
ggsave(filename = "/Users/jackmurzynowski/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/pca_pop_analysis.png", 
       plot = pca_analysis_results, device = "png", bg = 'white',
       width = 13, height = 12)

## Use paralell processing
run_simulation <- function(i) {
  run_results_i <- list()
  n_individuals <- 1000
  for (j in 1:length(genotype_data)) {
    print(paste0("On run: ", i, "/", 50, " for data: ", j, "/", length(genotype_data)))
    genotype_data_use <- genotype_data[[j]]
    
    # Randomly select 20 SNPs as causal
    set.seed(i)
    causal_snps <- selected_snps %>%
      slice_sample(., n = 20)
    
    phenotype_indv <- phenotype_assigment(
      genotype_data = genotype_data_use,
      causal_snps = causal_snps,
      n_cases = 500,
      genef_size = 1,
      envef_size = 0.2,
      heritability = 0.8,
      override_ncases = FALSE,
      seed = 123
    )
    
    print("Compiling data nicely")
    total_allele_info <- genotype_data_use$allele_info %>%
      mutate(causal = ifelse(id %in% causal_snps$id, T, F))
    index_cols <- causal_snps$index_col
    causal_snps_better <- causal_snps %>%
      select(-index_col) %>%
      left_join(phenotype_indv$case_control_stats) %>%
      mutate(index_col = index_cols)
    
    genetic_data <- list(
      causal_snp_data = causal_snps_better,
      causal_geno_pheno_data = phenotype_indv$geno_pheno_data,
      causal_allele_pheno_data = phenotype_indv$allele_pheno_data,
      all_snp_data = total_allele_info,
      all_genotype_data = genotype_data_use$genotype_info,
      all_allele_data = genotype_data_use$allele_data
    )
    
    directory <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_tests_50_50"
    file_name <- paste0("mixed_pop_50_50_",i,"_",j)
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
    ped_easy <- data.frame(FID = 1:n_individuals, IID = 1:n_individuals, PID = 0, MID = 0, Sex = sample(1:2, n_individuals, replace = TRUE))
    phenotype <- phenotype_indv$geno_pheno_data$phenotype
    ped_data <- cbind(ped_easy, phenotype, genetic_data$all_allele_data)
    fwrite(ped_data, file = paste0(file_name, ".ped"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # Write out the final bfile
    print(paste0("Making final bfile with name: '", file_name, "'"))
    system(paste0("plink --file ", file_name, " --make-bed --out ", file_name),ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Perform the PLINK GWAS
    print("Performing GWAS")
    system(paste0("plink --file ", file_name, " --assoc --ci 0.95 --maf 0.01 --nonfounders --out ", file_name),ignore.stdout = TRUE, ignore.stderr = TRUE)
    sim_results_assoc <- fread(paste0(file_name, ".assoc"), head = TRUE) %>%
      arrange(P)
    
    plink_file_name <- file_name
    system(paste0("rm ", directory, "/", plink_file_name, "*"))
    
    causal_snps <- genetic_data$causal_snp_data %>%
      select("CHR" = chr, "SNP" = id, "BP" = pos, "Control_T" = sum_control, "Case_T" = sum_case, "CsCtrl_R" = case_control_ratio) %>%
      inner_join(sim_results_assoc) %>%
      arrange(P)
    
    assoc_stats <- generate_assoc_stats(association_data = sim_results_assoc,
                                        causal_snps = causal_snps)
    inflation <- inflation_factor(sim_results_assoc$P)
    
    # # Get bonferroni and pval thresholds
    # bon_fdr <- calc_bonfdr(pval_data = sim_results_assoc$P)
    # 
    # # Generate plots for html report
    # print("Generating plots")
    # gwas_title <- paste0("Population AB₁")
    # if(!is.infinite(bon_fdr$bonferroni)){
    #   major_threshold <- -log10(bon_fdr$bonferroni)
    # } else {
    #   major_threshold <- NULL
    # }
    # if(!is.infinite(bon_fdr$fdr)){
    #   minor_threshold <-  -log10(bon_fdr$fdr)
    # } else {
    #   minor_threshold <- NULL
    # }
    # qqmanh_plots <- generate_qqmanhat(association_data = sim_results_assoc,
    #                                   causal_snps = causal_snps,
    #                                   plot_type = c('manhattan','qq'),
    #                                   logp = TRUE,
    #                                   major_threshold = major_threshold,
    #                                   minor_threshold = NULL,
    #                                   qq_thresholds = FALSE,
    #                                   min_pval = -log10(0.01),
    #                                   highlight_causal = TRUE,
    #                                   annotate_causal = FALSE,
    #                                   title = gwas_title,
    #                                   annotate_inflation = TRUE,
    #                                   directory = directory,
    #                                   filename = plink_file_name)
    
    
    # Mutate causal SNPs
    association_data_mut <- sim_results_assoc %>%
      mutate(true_causal = ifelse(SNP %in% causal_snps$SNP, 1, 0))
    
    # Create ROC curve using p-values
    roc_curve <- roc(response = association_data_mut$true_causal, predictor = 1 - log10(association_data_mut$P))
    auc <- roc_curve$auc
    
    run_results_i[[j]] <- assoc_stats %>%
      mutate(seed = i,
             auc = auc,
             inflation_fac = inflation,
             population = j)
    files <- list.files("/Users/jackmurzynowski/Library/Application Support/CloudDocs/session/i", full.names = TRUE)
    unlink(files, recursive = TRUE)
  }
  return(run_results_i)
}

# Get the number of available cores
num_cores <- detectCores()-1

# Create a cluster with the desired number of cores
cl <- makeCluster(num_cores)

# Export necessary variables and packages to the cluster
clusterEvalQ(cl, {
  library(dplyr)
  library(data.table)
  library(pROC)
  library(tidyr)
})
clusterExport(cl, c("genotype_data", "selected_snps", "phenotype_assigment", "generate_assoc_stats","inflation_factor"))


# Run the simulation in parallel
run_results <- parLapply(cl, 1:50, run_simulation)

# Stop the cluster
stopCluster(cl)

# Flatten the list of results
run_results_total <- do.call(c, run_results)

total_run_stats <- do.call(rbind,run_results_total)
qsave(total_run_stats, file = "total_comparator_stats.qs")
setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_tests_50_50")
total_run_stats <- qread("total_comparator_stats.qs")

measure_data <- list()
for(measure in c("TP","FP","auc","inflation_fac")){
  if(measure == "TP"){
    title_use <- "Mean GWAS True Positives"
    y_use <- "Mean  Normalized True Positives"
  } else if(measure == "FP"){
    title_use <- "Mean GWAS False Positives"
    y_use <- "Mean Normalized False Positives"
  } else if(measure == "auc"){
    title_use <- "Mean GWAS AUC"
    y_use <- "Mean Normalized AUC"
  } else {
    title_use <- "Mean GWAS Inflation"
    y_use <- "Mean Normalized Inflation"
  }
  
  measure_stats <- total_run_stats %>%
    select(population, measure = !!sym(measure)) %>%
    mutate(across(measure, ~ (.-min(.)) / (max(.) - min(.))))
  
  pop_combinations <- combn(unique(measure_stats$population), 2, simplify = FALSE)
  sig_results <- list()
  for (pair in pop_combinations) {
    pop1 <- pair[1]
    pop2 <- pair[2]
    
    pop1_data <- measure_stats %>% filter(population == pop1) %>% pull(measure)
    pop2_data <- measure_stats %>% filter(population == pop2) %>% pull(measure)
    
    # Ensure the lengths are equal for paired differences
    min_length <- min(length(pop1_data), length(pop2_data))
    differences <- pop1_data[1:min_length] - pop2_data[1:min_length]
    
    # Check normality of differences
    shapiro_result <- shapiro.test(differences)
    
    if (shapiro_result$p.value > 0.05) {
      # Use paired t-test if differences are normally distributed
      test_result <- t.test(pop1_data[1:min_length], pop2_data[1:min_length], paired = TRUE)
      test_type <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test if differences are not normally distributed
      test_result <- wilcox.test(pop1_data[1:min_length], pop2_data[1:min_length], paired = TRUE, exact = FALSE)
      test_type <- "Wilcoxon test"
    }
    sig_results[[paste(pop1, "vs", pop2)]] <- list(
      test_result = test_result,
      test_type = test_type
    )
  }
  populations <- unique(total_run_stats$population)
  model_matrix <- matrix("", nrow = length(populations), ncol = length(populations),
                         dimnames = list(populations, populations))
  for (pops_comp in names(sig_results)) {
    measure_results <- sig_results[[pops_comp]]
    pops <- strsplit(pops_comp, " vs ")[[1]]
    pop1 <- pops[1]
    pop2 <- pops[2]
    p_value <- measure_results$test_result$p.value
    V_stat <- if (!is.na(measure_results$test_result$statistic)) measure_results$test_result$statistic else ""
    model_matrix[pop2, pop1] <- p_value
    model_matrix[pop1, pop2] <- V_stat
  }
  model_matrix[upper.tri(model_matrix)] <- ""
  
  # Get metrics
  model_summary <- measure_stats %>%
    group_by(population) %>%
    summarize(
      Mean = mean(measure),
      Lower_CI = mean(measure) - qt(1 - 0.05 / 2, n() - 1) * sd(measure) / sqrt(n()),
      Upper_CI = mean(measure) + qt(1 - 0.05 / 2, n() - 1) * sd(measure) / sqrt(n()),
      SE = sd(measure) / sqrt(n())
    ) %>%
    arrange(.,desc(Mean)) %>%
    mutate(rank = row_number())
  model_summary$population <- factor(model_summary$population, levels = model_summary$population)
  model_rank <- model_summary %>% select(population,Upper_CI,rank)
  
  
  
  
  signif_comparisons <- model_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "Population1") %>%
    pivot_longer(cols = -Population1, names_to = "Population2", values_to = "p_value") %>%
    filter(p_value != "") %>%
    mutate(p_value = as.numeric(p_value)) %>%
    mutate(significance = case_when(
      p_value < 0.0001 ~ "****",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*"
    )) %>%
    filter(!if_any(everything(), is.na)) %>%
    select(Population1, Population2, p_value, significance) %>%
    left_join(model_rank, by = c("Population1" = "population")) %>%
    left_join(model_rank, by = c("Population2" = "population")) %>%
    mutate(priCI = ifelse(Upper_CI.x > Upper_CI.y, Upper_CI.x, Upper_CI.y),
           secCI = ifelse(Upper_CI.x < Upper_CI.y, Upper_CI.x, Upper_CI.y)) %>%
    arrange(priCI) %>%
    group_by(secCI) %>%
    mutate(secrank = cur_group_id()) %>%
    ungroup() %>%
    arrange(priCI,secCI) 
  # Get coordinates for sig lines
  comparisons_list <- list()
  annotations <- c()
  y_positions <- c()
  base_y_position <- max(model_summary$Upper_CI) * 1
  max_ci_track <- c()
  prev_pos <- NULL
  
  range <- max(model_summary$Upper_CI)-min(model_summary$Lower_CI)
  for (i in 1:nrow(signif_comparisons)) {
    comparisons_list[[i]] <- c(signif_comparisons$Population1[i], signif_comparisons$Population2[i])
    annotations[i] <- signif_comparisons$significance[i]
    priCI <- signif_comparisons[i,]$priCI
    secrank <- signif_comparisons[i,]$secrank
    predict_pos <- priCI + (range*0.05) * secrank
    if(secrank != 1){
      if(prev_pos > priCI){
        predict_pos <- prev_pos + (range*0.1)
      }
    } else {
      if(i != 1){
        if(prev_pos + (range*0.1) > priCI){
          predict_pos <- prev_pos + (range*0.1)
        } else {
          predict_pos <- predict_pos - (range*0.05) * secrank # Initial elevation above the plot
        }
      } else {
        predict_pos <- predict_pos - (range*0.05) * secrank # Initial elevation above the plot
      }
    }
    prev_pos <- predict_pos
    y_positions[i] <- predict_pos
  }
  
  model_labels <- model_summary %>% 
    mutate(population = case_when(
      str_detect(population, "1") ~ "AB1",
      str_detect(population, "2") ~ "AB2",
      str_detect(population, "3") ~ "A",
      str_detect(population, "4") ~ "B",
      TRUE ~ population
    )) %>% select(population) %>% as.matrix() %>% as.character()
  
  y_breaks <- pretty(c(0, max(model_summary$Upper_CI) * 1.01), n = 5)
  y_labels <- y_breaks
  barsig_plot <- ggplot(model_summary, aes(x = population, y = Mean)) +
    geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    labs(title = title_use, y = paste0(y_use), x = "Population") +
    theme_minimal() +
    geom_signif(
      comparisons = comparisons_list,
      annotations = annotations,
      y_position = y_positions,
      tip_length = 0.03,
      textsize = 4,
      vjust = 0.4
    ) + 
    scale_y_continuous(
      breaks = y_breaks, 
      labels = y_labels,
      expand = expansion(mult = c(0, 0.2))
    ) +
    scale_x_discrete(labels = model_labels) + 
    theme(
      plot.margin = unit(c(3, 1, 1, 1), "lines"),
      panel.grid.major.y = element_line(color = "grey", linewidth = 0.5),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(hjust = 0.25, vjust = 2.5, size = 14),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12)
    )
  
  measure_data[[measure]] <- list(plot = barsig_plot,
                                  data = model_summary)
  
}

# Plot all three results in a graph
pop_structure_plots <- plot_grid(
  measure_data$TP$plot,
  measure_data$FP$plot,
  measure_data$auc$plot,
  measure_data$inflation_fac$plot,
  ncol = 4,
  labels = c("A","B","C", "D"),
  label_size = 20
)
ggsave(filename = "/Users/jackmurzynowski/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_analysis.png", 
       plot = pop_structure_plots, device = "png", bg = 'white',
       width = 16, height = 10)

# Get the manhattan plots
mix_50_50_man <- qqmanh_plots_50_50$manhattan
mix_50_50_qq <- qqmanh_plots_50_50$qq
reg_man <- qqmanh_plots_reg$manhattan
reg_qq <- qqmanh_plots_reg$qq

pop_structure_top_plots <- plot_grid(
  reg_man,
  mix_50_50_man,
  mix_50_50_qq,
  ncol = 3,
  labels = c("A","B","C"),
  label_size = 20
)
pop_structure_bottom_plots <- plot_grid(
  measure_data$TP$plot,
  measure_data$FP$plot,
  measure_data$auc$plot,
  measure_data$inflation_fac$plot,
  ncol = 4,
  labels = c("D","E","F", "G"),
  label_size = 20
)
total_grid <- plot_grid(
  pop_structure_top_plots,
  pop_structure_bottom_plots,
  ncol = 1)
ggsave(filename = "/Users/jackmurzynowski/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/mix_pop_analysis_total.png", 
       plot = total_grid, device = "png", bg = 'white',
       width = 20, height = 14)

