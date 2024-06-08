# Function to iterate over a set of parameters, creating graphs and storing data
rm(list = ls())
directory <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test"
setwd(directory)
temp_dir <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/parallel_tmp_dir"
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_sim.R")
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_assoc_stats_functions.R")
library(DT)
library(rlist)
library(tidyverse)

# Temp files appear in /Users/jackmurzynowski/Library/Application Support/CloudDocs/session/i
# Cached increases appear in /Users/jackmurzynowski/Library/Caches/CloudKit


n_individuals = c(200,400,800,1600,3200,4000)
n_snps = c(2500,5000,10000,25000,50000,100000)
pair_ld = 0.999
coverage = 0.01
temp_dir = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/parallel_tmp_dir"

# Start loop to get all individuals, repeating experiments for each pairing
time_results <- data.frame(
  n_snps = integer(),
  n_individuals = integer(),
  time = numeric()
)
for(snp_size in n_snps){
  for(indv in n_individuals){
    for(i in 1:5){
      print(paste0("Testing SNP: ",snp_size," for Individual: ",indv," on run ",i))
      # Generate the genetic data
      time_taken <- system.time({
        genotype_data <- simulate_genotype_data(
          n_individuals = indv,
          n_snps = snp_size,
          pair_ld = pair_ld,
          coverage = coverage,
          chrs = 1:22,
          human = TRUE,
          chrom_data = NULL,
          parallel = TRUE,
          seed = 123,
          temp_dir = temp_dir
        )
      })
      # Store the results in the dataframe
      time_results <- rbind(time_results, data.frame(
        n_snps = snp_size,
        n_individuals = indv,
        time = time_taken["elapsed"]
      ))
    }
  }
}
saveRDS(summary_df, file = file.path(directory, "timing_data_summary.rds"))

# Convert time_taken to log scale for better visualization
calculate_summary <- function(df) {
  df %>%
    group_by(n_snps, n_individuals) %>%
    summarise(
      mean_time = mean(time),
      lower_ci = mean(time) - 1.96 * sd(time) / sqrt(n()),
      upper_ci = mean(time) + 1.96 * sd(time) / sqrt(n())
    )
}
summary_df <- calculate_summary(time_results)
time_results$time_log <- log10(time_results$time)

max(summary_df$mean_time)

summary_df[9,] <- c(4.7)

# Plot using ggplot2
sim_runtime <- ggplot(summary_df, aes(x = n_individuals, y = mean_time, group = n_snps, shape = as.factor(n_snps))) +
  geom_line(aes(linetype = as.factor(n_snps)), linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = "95% CI"), width = 0.2) +
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(name = "", values = "red") +
  labs(title = "Genomic Data Simulation Time",
       x = "Number of individuals",
       y = "Time (sec)",
       shape = "Number of SNPs",
       linetype = "Number of SNPs") +
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 1), color = guide_legend(order = 2, override.aes = list(linetype = "solid"))) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, "in")
  )
directory <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test"
setwd(directory)
response <- "auc"
ggsave(filename = paste0(directory,"/",response,"_sim_runtime.png"), plot = sim_runtime, device = "png", bg = 'white')


# Test a bunch of models
generate_and_fit_models <- function(data, response, predictors) {
  # Create a list of all combinations of regular, polynomial, and square root terms
  terms <- c("", predictors, paste0("I(", predictors, "^2)"), paste0("sqrt(", predictors, ")"), paste0("log(",predictors,")"))
  formulas <- c()
  
  # Generate all possible non-empty combinations of terms
  for (i in 2:length(terms)) {
    combinations <- combn(terms, i, simplify = FALSE)
    for (combo in combinations) {
      # Ensure at least one predictor is included in the formula
      if (any(predictors %in% combo) || any(grepl("sqrt", combo)) || any(grepl("I", combo))) {
        formulas <- c(formulas, paste(response, "~", paste(combo[combo != ""], collapse = " + ")))
      }
    }
  }
  # Initialize a list to store the models and their summaries
  models <- list()
  model_summaries <- data.frame(
    formula = character(),
    adj_r_squared = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Fit models for each formula and store the results
  for (f in formulas) {
    model <- lm(as.formula(f), data = data)
    summary_model <- summary(model)
    
    models[[f]] <- model
    model_summaries <- rbind(model_summaries, data.frame(
      formula = f,
      adj_r_squared = summary_model$adj.r.squared,
      stringsAsFactors = FALSE
    ))
  }
  
  list(models = models, model_summaries = model_summaries)
}

response <- "time"
predictors <- c("n_individuals", "n_snps")

model_results <- generate_and_fit_models(data = time_results, 
                                         response = response, 
                                         predictors = predictors)

all_model_formulas_with_coefficients <- function(model_summaries, models) {
  # Initialize an empty list to store the results
  results <- list()
  
  # Loop through each model
  for (i in 1:nrow(model_summaries)) {
    formula <- model_summaries$formula[i]
    adj_r_squared <- model_summaries$adj_r_squared[i]
    
    # Extract the model
    model <- models[[formula]]
    
    # Extract coefficients
    coefficients <- summary(model)$coefficients
    coeff_str <- paste0(rownames(coefficients), " * ", round(coefficients[, "Estimate"], 4), collapse = " + ")
    
    # Store the results
    results[[i]] <- data.frame(
      formula = formula,
      coefficients = coeff_str,
      adj_r_squared = adj_r_squared,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine the list into a dataframe and sort by adjusted R-squared
  results_df <- do.call(rbind, results)
  results_df <- results_df[order(-results_df$adj_r_squared), ]
  
  return(results_df)
}
all_formulas_with_coefficients <- all_model_formulas_with_coefficients(
  model_results$model_summaries, 
  model_results$models) %>%
  unique() %>%
  filter(grepl('time \\~ log\\(n_individuals\\) \\+ log\\(n_snps\\)',formula)) %>%
  slice(1:5)

test_sqrt <- lm(time ~ log(n_individuals) + log(n_snps), data = time_results)
summary(test_sqrt)

best_formula_with_coefficients <- best_model_formula_with_coefficients(
  model_results$model_summaries, 
  model_results$models)
