# Function to iterate over a set of parameters, creating graphs and storing data
rm(list = ls())
setwd("~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test")
temp_dir <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/parallel_tmp_dir"
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_sim.R")
source(file = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/genetic_assoc_stats_functions.R")
library(DT)
library(rlist)
library(tidyverse)


# Generate the genetic data
dir_use <- "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/Sanity_test"
plink_file_name <- "total_association"
manhat_plt <- list()
for(herit in 0:1){
  genetic_data <- genetic_sim(
    n_individuals = 4000,
    n_snps = 100000,
    chrs = 1:22,
    n_causal_snps = 1,
    causal_maf = 0.05,
    n_cases = 2000,
    genef_size = 1,
    envef_size = 0,
    heritability = herit,
    pair_ld = 0.999,
    coverage = 0.01,
    write_plink = TRUE,
    file_name = plink_file_name,
    directory = dir_use,
    human = TRUE,
    chrom_data = NULL,
    override_ncases = FALSE,
    parallel = TRUE,
    seed = set_seed)
  gc()
  
  # Run PLINK GWAS
  print("Running PLINK")
  system(paste0("plink --bfile ",plink_file_name," --assoc --ci 0.95 --maf 0.01 --nonfounders --out ",plink_file_name))
  sim_results_assoc <- fread(paste0(plink_file_name,".assoc"), head=TRUE) %>%
    arrange(P)
  
  # Remove all PLINK files
  system(paste0("rm ",directory,"/",plink_file_name,"*"))
  
  # Get causal SNP data
  causal_snps <- genetic_data$causal_snp_data %>%
    select("CHR" = chr,"SNP" = id,"BP" = pos,"MAF" = maf_observed,"N_0" = zero_count,"N_1" = one_count,"N_2" = two_count,
           "Control_T" = sum_control, "Case_T" = sum_case, "CsCtrl_R" = case_control_ratio) %>%
    inner_join(sim_results_assoc) %>%
    arrange(P)
  
  total_run_data <- list(sim_results_assoc = sim_results_assoc,
                         causal_snps = causal_snps)
  
  
  # Run PLINK GWAS
  print("Running PLINK")
  system(paste0("plink --bfile ",plink_file_name," --assoc --ci 0.95 --maf 0.01 --nonfounders --out ",plink_file_name))
  sim_results_assoc <- fread(paste0(plink_file_name,".assoc"), head=TRUE) %>%
    arrange(P)
  
  # Remove all PLINK files
  system(paste0("rm ",directory,"/",plink_file_name,"*"))
  
  # Get causal SNP data
  causal_snps <- genetic_data$causal_snp_data %>%
    select("CHR" = chr,"SNP" = id,"BP" = pos,"MAF" = maf_observed,"N_0" = zero_count,"N_1" = one_count,"N_2" = two_count,
           "Control_T" = sum_control, "Case_T" = sum_case, "CsCtrl_R" = case_control_ratio) %>%
    inner_join(sim_results_assoc) %>%
    arrange(P)
  
  total_run_data <- list(sim_results_assoc = sim_results_assoc,
                         causal_snps = causal_snps)
  
  bon_pval <- 0.05 / length(total_run_data$sim_results_assoc$P)
  bon_fdr <- calc_bonfdr(pval_data = total_run_data$sim_results_assoc$P)
  
  if(herit == 0){
    title_use <- "No Genetic Effect"
  } else {
    title_use <- "Pure Genetic Effect"
  }
  manhat_plt[[herit]] <- generate_qqmanhat(association_data = total_run_data$sim_results_assoc,
                                       causal_snps = total_run_data$causal_snps,
                                       plot_type = c('manhattan'),
                                       logp = TRUE,
                                       major_threshold = -log10(bon_pval),
                                       minor_threshold = NULL,
                                       qq_thresholds = FALSE,
                                       min_pval = -log10(0.01),
                                       highlight_causal = TRUE,
                                       annotate_causal = FALSE,
                                       title = title_use,
                                       annotate_inflation = TRUE,
                                       directory = dir_use,
                                       filename = "one_assoc")
}


manhat_plt_none <- manhat_plt[[1]] + 
  labs(title = "Only Environmental Effect") +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12))

manhat_plt_one <- manhat_plt[[2]] + 
  labs(title = "Only Genetic Effect") +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16, face = "plain"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12))

# Combine the plots side by side with labels A and B
combined_plot <- plot_grid(
  plot_grid(ggdraw() + draw_label("A", fontface = 'bold', x = 0, hjust = 0, size = 20), manhat_plt_none, ncol = 1, rel_heights = c(0.05, 1)),
  plot_grid(ggdraw() + draw_label("B", fontface = 'bold', x = 0, hjust = 0, size = 20), manhat_plt_one, ncol = 1, rel_heights = c(0.05, 1)),
  ncol = 2,
  rel_widths = c(1, 1)
)

# Investigate the LD structure of the causal SNP
total_run_data$causal_snps
top_snps <- total_run_data$sim_results_assoc %>%
  filter(P <= bon_pval)
top_snps %>%
  mutate(distance = abs(BP - 968376)) %>%
  arrange(desc(distance))
ld_data <- generate_ld_metrics(snp_data = genetic_data$all_snp_data,
                               genotype_data = genetic_data$all_genotype_data,
                               snp_id = "rs11_968376",
                               window_size = 5000,
                               index = FALSE,
                               directory = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/Sanity_test",
                               filename = "testing_ld")

ld_long <- as.data.frame(ld_data$ld_matrix) %>%
  rownames_to_column(var = "SNP1") %>%
  pivot_longer(-SNP1, names_to = "SNP2", values_to = "value")

unique_snps <- unique(c(ld_long$SNP1, ld_long$SNP2))

# Create a data frame to map SNPs to their highlight status
highlight_status <- data.frame(
  SNP = unique_snps,
  fontface = ifelse(unique_snps == top_snps$SNP[1], "bold", "plain"),
  color = ifelse(unique_snps == top_snps$SNP[1], "green", "grey"),
  size = ifelse(unique_snps %in% top_snps$SNP, 12, 6) # Larger font size for highlighted SNPs
) %>% 
  mutate(color = case_when(unique_snps == top_snps$SNP[1] ~ "green", T ~ 
                             ifelse(unique_snps %in% top_snps$SNP, "black", "grey"))) 

highlighted_positions <- highlight_status %>%
  filter(size == 12) %>%
  mutate(position = match(SNP, unique_snps))

# Determine the positions for the outer lines
line_positions <- highlighted_positions$position[c(1, nrow(highlighted_positions))]
green_line_position <- match(top_snps$SNP[1], unique_snps)

min_pos <- min(line_positions)
max_pos <- max(line_positions)
highlight_status <- highlight_status %>%
  mutate(size = ifelse(size == 12, 12, 
                       ifelse(row_number() < min_pos, pmax(1, 10 - (min_pos - row_number())),
                              ifelse(row_number() > max_pos, pmax(1, 10 - (row_number() - max_pos)), 8))))

# Create the heatmap
heatmap_plot <- ggplot(ld_long, aes(SNP2, SNP1)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "LD R^2") +
  labs(title = "Local Region Around Causal SNP") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,
                               face = highlight_status$fontface[match(ld_long$SNP2, highlight_status$SNP)],
                               color = highlight_status$color[match(ld_long$SNP2, highlight_status$SNP)],
                               size = highlight_status$size[match(ld_long$SNP2, highlight_status$SNP)]),
    axis.text.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12, vjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16, face = "plain")
  ) +
  geom_segment(aes(x = line_positions[1]-0.5, xend = line_positions[1]-0.5, y = line_positions[1]-0.5, yend = line_positions[2]+0.5), linetype = "dashed", color = "blue") +
  geom_segment(aes(x = line_positions[2]+0.5, xend = line_positions[2]+0.5, y = line_positions[1]-0.5, yend = line_positions[2]+0.5), linetype = "dashed", color = "blue") +
  geom_segment(aes(x = line_positions[1]-0.5, xend = line_positions[2]+0.5, y = line_positions[1]-0.5, yend = line_positions[1]-0.5), linetype = "dashed", color = "blue") +
  geom_segment(aes(x = line_positions[1]-0.5, xend = line_positions[2]+0.5, y = line_positions[2]+0.5, yend = line_positions[2]+0.5), linetype = "dashed", color = "blue") +
  geom_segment(aes(x = green_line_position, xend = green_line_position, y = min(SNP1), yend = max(SNP1)), linetype = "dashed", color = "green") +
  geom_segment(aes(x = min(SNP2), xend = max(SNP2), y = green_line_position, yend = green_line_position), linetype = "dashed", color = "green") +
  annotate("text", x = green_line_position, y = 1.5, label = "Causal SNP", color = "green", size = 5, angle = 0, hjust = -0.1) +
  annotate("text", x = line_positions[2] + 0.5, y = line_positions[1] - 1.5, label = "Significant SNPs", color = "blue", size = 5, angle = 0, hjust = 1)

final_combined_plot <- plot_grid(
  combined_plot,
  plot_grid(ggdraw() + draw_label("C", fontface = 'bold', x = 0, hjust = 0, size = 20), heatmap_plot, ncol = 1, rel_heights = c(0.05, 1)),
  ncol = 1,
  rel_heights = c(0.8, 1)
)
ggsave(filename = paste0(directory,"/check_sim_perf.png"), plot = final_combined_plot, device = "png", bg = 'white')

