# Functions for calculating association statistics
# Check what packages are needed and load them in
# NCmisc::list.functions.in.file(filename = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_assoc_stats_functions.R",
#                                alphabetic = TRUE)
library(ComplexHeatmap)
library(pROC)
library(ggplot2)
library(stats)
library(tidyverse)

# Function to calculate the inflation factor
inflation_factor <- function(p_values) {
  print("Calculating inflation factor")
  chisq <- qchisq(1 - p_values, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  return(lambda)
}

# Function to caucluate power at significance threshold(s) - a vector can be supplied
calculate_power <- function(causal_snps, pval_thresholds) {
  print(paste0("Calculating power factor at ",length(pval_thresholds)," pval thresholds"))
  power <- sapply(pval_thresholds, function(threshold) {
    sum(causal_snps$P < threshold) / nrow(causal_snps)
  })
  data.frame(Pval_Threshold = pval_thresholds, Power = power)
}

# Calcualte bonferroni and FDR pval thresholds
calc_bonfdr <- function(pval_data){
  bonfdr <- list()
  bonfdr$bonferroni <- 0.05 / length(pval_data)
  bonfdr$fdr <- max(pval_data[p.adjust(pval_data, method = "fdr") < 0.05], na.rm = TRUE)
  return(bonfdr)
}

# Function to generate TP,FP,TN,FN,Sens,Spec and F statistic data for GWAS at either a bonferroni pval or FDR pval
generate_assoc_stats <- function(association_data, causal_snps) {
  print("Generating association statistics")
  if(!is.data.frame(association_data) | !all(c("P") %in% colnames(association_data))){
    stop("association_data must be a dataframe containing a 'P' column")
  }
  if(is.data.frame(causal_snps)){
    if(!"SNP" %in% colnames(causal_snps)){
      stop("causal_snps must contain a 'SNP' column")
    } else {
      causal_snps_vect <- causal_snps$SNP
    }
  } else {
    causal_snps_vect <- causal_snps
  }
  
  # Calculate thresholds
  bonferroni_threshold <- 0.05 / nrow(association_data)
  fdr_threshold <- max(association_data$P[p.adjust(association_data$P, method = "fdr") < 0.05], na.rm = TRUE)
  
  # Define a helper function to calculate stats
  calc_assoc_stats <- function(data, threshold) {
    tp <- sum(data$P < threshold & data$SNP %in% causal_snps_vect)
    fp <- sum(data$P < threshold & !data$SNP %in% causal_snps_vect)
    fn <- sum(data$P >= threshold & data$SNP %in% causal_snps_vect)
    tn <- sum(data$P >= threshold & !data$SNP %in% causal_snps_vect)
    sensitivity <- tp / (tp + fn)
    specificity <- tn / (tn + fp)
    precision <- tp / (tp + fp)
    f_stat <- 2*((precision*sensitivity)/(precision+sensitivity))
    return(c(TP = tp, FP = fp, TN = tn, FN = fn, Sensitivity = sensitivity, Specificity = specificity,Precision = precision,F_statistic = f_stat))
  }
  
  # Calculate stats for Bonferroni and FDR thresholds
  bonferroni_stats <- calc_assoc_stats(association_data, bonferroni_threshold)
  fdr_stats <- calc_assoc_stats(association_data, fdr_threshold)
  
  # Create assoc_stats table
  assoc_stats <- data.frame(
    Type = c("Bonferroni", "FDR"),
    Threshold = c(bonferroni_threshold, fdr_threshold),
    t(cbind(bonferroni_stats, fdr_stats))
  )
  row.names(assoc_stats) <- NULL
  
  return(assoc_stats)
}


# Function to generate custom QQ and Manhattan plots
generate_qqmanhat <- function(association_data,causal_snps,plot_type = c('manhattan','qq'),logp = TRUE, major_threshold = -log10(5e-08),
                              minor_threshold = -log10(1e-05),qq_thresholds = FALSE, highlight_causal = TRUE,annotate_causal = FALSE,
                              title = "Assoc data",min_pval = -log10(0.01),annotate_inflation = TRUE, directory = NULL, filename = NULL){
  print("Generating a qq/manhattan plot of results")
  
  # Run preliminary checks to ensure the function can execute
  if(!is.data.frame(association_data) | !all(c("CHR","BP","SNP","P") %in% colnames(association_data))){
    stop("association_data must be a dataframe containing 'CHR', 'BP', 'SNP' and 'P' columns")
  }
  if(!is.data.frame(causal_snps) | !all(c("SNP","P") %in% colnames(causal_snps))){
    stop("causal_snps must be a dataframe containing a 'SNP' and 'P' column")
  }
  if(!exists("ggmanhattan",mode = "function") | !exists("ggqq",mode = "function")){
    # Load in the modified Manhattan and qq functions
    source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/ggqqman.R")
  }
  
  # Highlight significant SNPs
  if(highlight_causal){
    if(logp){
      colour_major_threshold <-ifelse(is.null(major_threshold),1000000,major_threshold) 
      colour_minor_threshold <-ifelse(is.null(minor_threshold),1000000,minor_threshold) 
      colour_causal <- causal_snps %>%
        mutate(logp = -log10(P),
               colour = ifelse(logp > colour_minor_threshold, "orange","red"),
               colour = ifelse(logp > colour_major_threshold, "green",colour))
    } else {
      colour_major_threshold <-ifelse(is.null(major_threshold),0,major_threshold) 
      colour_minor_threshold <-ifelse(is.null(minor_threshold),0,minor_threshold) 
      colour_causal <- causal_snps %>%
        mutate(colour = ifelse(P < colour_minor_threshold, "orange","red"),
               colour = ifelse(P < colour_major_threshold, "green",colour))
    }
    highlight_snps <- list(name = colour_causal$SNP,colour = colour_causal$colour)
  } else {
    highlight_snps <- NULL 
  }
  
  # If both are requested then change titles and filenames accordingly
  if(all(c('qq','manhattan') %in% plot_type)){
    qq_title <- paste0("QQ: ",title)
    qq_filename <- paste0("qq_",filename)
    manhat_title <- paste0("Manhattan: ",title)
    manhat_filename <- paste0("manhattan_",filename)
  } else {
    qq_title <- title
    manhat_title <- title
    qq_filename <- filename
    manhat_filename <- filename
  }
  
  # Generate the manhattan plot
  qqmanplot <- list()
  for(i in plot_type){
    if(i == "manhattan"){
      print("Generating Manhattan plot")
      manhat_plot <- ggmanhattan(association_data,
                                 chr = "CHR", 
                                 bp = "BP", 
                                 p = "P", 
                                 snp = "SNP", 
                                 title = manhat_title,
                                 logp = logp,
                                 min_pval = min_pval,
                                 col = c("lightblue", "blue"), 
                                 suggestiveline = minor_threshold,
                                 genomewideline = major_threshold,
                                 suggestLinecol = 'red',
                                 genomeLinecol = 'blue',
                                 highlight = highlight_snps,
                                 annotateHighlight = annotate_causal,
                                 chrlabs = c(1, 2, 3, 4, 5, 7, 9, 12, 15, 19))
      
      # Save the plot
      if(!is.null(directory) & !is.null(filename)) {
        if(!dir.exists(directory)) {
          dir.create(directory, recursive = TRUE)
        }
        print("Saving Manhattan plot")
        ggsave(filename = file.path(directory, paste0(manhat_filename,".png")), plot = manhat_plot, device = "png", bg = 'white')
      }
      
      # Save results to list
      qqmanplot$manhattan <- manhat_plot
      
    } else if(i == 'qq'){
      print("Generating QQ plot")
      if(!qq_thresholds){
        minor_threshold_qq <- NULL
        major_threshold_qq <- NULL
      } else {
        minor_threshold_qq <- minor_threshold
        major_threshold_qq <- major_threshold
      }
      qq_plot <- ggqq(association_data,
                      p = "P",
                      snp = "SNP",
                      title = qq_title,
                      logp = logp, 
                      inflation = annotate_inflation,
                      min_pval = min_pval,
                      point_colour = "blue", 
                      diag_colour = "red", 
                      suggestiveline = minor_threshold_qq,
                      genomewideline = major_threshold_qq,
                      suggestLinecol = 'red',
                      genomeLinecol = 'blue',
                      highlight = highlight_snps,
                      annotateHighlight = FALSE)
      
      # Save the plot
      if(!is.null(directory) & !is.null(filename)) {
        if(!dir.exists(directory)) {
          dir.create(directory, recursive = TRUE)
        }
        print("Saving QQ plot")
        ggsave(filename = file.path(directory, paste0(qq_filename,".png")), plot = qq_plot, device = "png", bg = 'white')
      }
      
      # Save results to list
      qqmanplot$qq <- qq_plot
      
    }
  }
  
  # If the plot only contains one entry then just return this 
  if(length(qqmanplot) == 1){
    qqmanplot <- qqmanplot[[1]]
  }
  return(qqmanplot)
}

# Function to generate an ROC plot along with max MAF, best P threshold and dataframes containing this info
generate_roc <- function(association_data,causal_snps,colour = "viridis",title = "ROC Curve for GWAS Results using P-values",
                         directory = NULL, filename = NULL){
  print("Generating ROC plots and data")
  if(!is.data.frame(association_data) | !all(c("SNP","P") %in% colnames(association_data))){
    stop("association_data must be a dataframe containing a 'SNP' and 'P' column")
  }
  if(is.data.frame(causal_snps)){
    if(!"SNP" %in% colnames(causal_snps)){
      stop("causal_snps must contain a 'SNP' column")
    } else {
      causal_snps_vect <- causal_snps$SNP
    }
  } else {
    causal_snps_vect <- causal_snps
  }
  
  # Mutate causal SNPs
  association_data_mut <- association_data %>%
    mutate(true_causal = ifelse(SNP %in% causal_snps_vect, 1, 0))
  
  # Create ROC curve using p-values
  roc_curve <- roc(response = association_data_mut$true_causal, predictor = 1 - log10(association_data_mut$P))
  
  # Extract thresholds and corresponding TPR and FPR
  thresholds <- roc_curve$thresholds
  tpr <- roc_curve$sensitivities
  fpr <- 1 - roc_curve$specificities
  precision <- tpr / (tpr + fpr)
  f_score <- 2 * ((precision * tpr) / (precision + tpr))
  max_f_score <- max(f_score, na.rm = TRUE)
  max_f_index <- which(f_score == max_f_score)
  max_f_tpr <- tpr[max_f_index]
  max_f_fpr <- fpr[max_f_index]
  max_f_threshold <- thresholds[max_f_index]
  max_f_pval <- 1 - max_f_threshold
  max_f_point <- data.frame(
    Max_F_score = max_f_score,
    False_Positive_Rate = max_f_fpr,
    True_Positive_Rate = max_f_tpr,
    P_value = max_f_pval
  )
  
  # Combine the data into a data frame
  roc_data <- data.frame(
    pval = 1 - thresholds,
    True_Positive_Rate = tpr,
    False_Positive_Rate = fpr,
    neg_log_pval = -log10(1 - thresholds)
  )
  
  # Calculate AUC
  auc_value <- auc(roc_curve)
  
  # Get p-value points
  # pval_points <- roc_data %>%
  #   filter(True_Positive_Rate != lead(True_Positive_Rate, default = first(True_Positive_Rate))) %>%
  #   mutate(label = format(pval, digits = 3, scientific = TRUE))
  # target_tpr_values <- c(0.25, 0.5, 0.75, 1)
  # find_closest_tpr <- function(target, tpr_values) {
  #   tpr_values[which.min(abs(tpr_values - target))]
  # }
  # closest_tpr_values <- sapply(target_tpr_values, find_closest_tpr, roc_data$True_Positive_Rate)
  # pval_points <- roc_data %>%
  #   filter(True_Positive_Rate %in% closest_tpr_values) %>%
  #   group_by(True_Positive_Rate) %>%
  #   slice_min(order_by = pval) %>%
  #   ungroup() %>%
  #   mutate(label = format(pval, digits = 3, scientific = TRUE))
  # pval_points$label <- factor(pval_points$label, levels = unique(pval_points$label[order(pval_points$pval,decreasing = T)]))
  
  # Plot the ROC curve with ggplot2
  roc_curve_plot <- ggplot(roc_data, aes(x = False_Positive_Rate, y = True_Positive_Rate, color = neg_log_pval)) +
    geom_path(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") + # Add diagonal baseline
    geom_point(data = max_f_point, aes(x = False_Positive_Rate, y = True_Positive_Rate), size = 4, color = "red") +
    annotate("text", x = 0.65, y = 0.40, label = paste("Max F-score (", round(max_f_score, 2), "),"), 
             size = 4, hjust = 0, color = "red") +
    annotate("text", x = 0.65, y = 0.35, label = paste("p-value =", format(max_f_pval, digits = 3, scientific = TRUE)), 
             size = 4, hjust = 0, color = "red") +
    annotate("text", x = 0.65, y = 0.2, label = paste("AUC =", round(auc_value, 3)), 
             size = 4, hjust = 0) +
    annotate("text", x = 0.65, y = 0.15, label = paste("Causal variants =", length(causal_snps_vect)), 
             size = 4, hjust = 0) +
    labs(title = title,
         x = "False Positive Rate",
         y = "True Positive Rate",
         color = "-log10(p)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) + 
    scale_color_viridis_c(option = colour)
  
  # Save the plot
  if(!is.null(directory) & !is.null(filename)) {
    if(!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
    }
    ggsave(filename = file.path(directory, paste0(filename,".png")), plot = roc_curve_plot, device = "png", bg = 'white')
  }
  
  return(list(roc_plot = roc_curve_plot, roc_data = roc_data, max_f_data = max_f_point, roc_info = roc_curve, auc = auc_value))
}



# Calculating haplotype frequencies from genotypic data
generate_ld_metrics <- function(snp_data,genotype_data,snp_id,window_size,index = TRUE,directory = NULL,filename = NULL){
  print("Generating LD heatmap and matrix")
  if(!is.data.frame(snp_data) | !all(c("causal","chr","pos","index_col") %in% colnames(snp_data))){
    stop("snp_data must be a dataframe containing 'chr', 'pos', 'causal' and 'index_col' columns")
  }
  
  # Filter for the snp of interest
  snp_interest <- snp_data %>% 
    filter(id == snp_id)
  if(nrow(snp_interest)==0){
    stop("snp_data doesn't contain the SNP of interest")
  }
  chrom_filt <- snp_data %>%
    filter(chr == snp_interest$chr)
  
  
  # Get the range of SNPs we wish to calculate LD for
  interest_range <- if(index){
    chrom_filt %>%
      filter(index_col >= snp_interest$index_col - window_size & index_col <= snp_interest$index_col + window_size)
  } else {
    chrom_filt %>%
      filter(pos >= snp_interest$pos - window_size & pos <= snp_interest$pos + window_size)
  }
  
  if(!all(interest_range$id %in% colnames(genotype_data))){
    stop("genotype_data columns do not match names of requested SNPs from snp_data")
  } else if(!exists("gen_ld_matrix",mode = "function")){
    # Load in the ld function
    source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/gen_ld_matrix.R")
  }
  
  # Get the genotype of the interested-in SNPs
  index_genotype_data <- genotype_data[,interest_range$index_col]
  
  # Calculate LD based off genotype
  ld_data <- gen_ld_matrix(index_genotype_data)
  
  # Save the plot
  if(!is.null(directory) & !is.null(filename)) {
    if(!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
    }
    filepath <- file.path(directory, paste0(filename,".png"))
    png(filepath, width = 800, height = 800, bg = "white")
    # Draw the heatmap
    draw(ld_data$ld_heatmap)
    # Close the graphics device
    dev.off()
  }
  return(ld_data)
}