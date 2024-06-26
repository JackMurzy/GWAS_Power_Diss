# NCmisc::list.functions.in.file(filename = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/ggqqman.R",
#                                                                alphabetic = TRUE)
library(ggplot2)
library(stats)
library(tidyverse)


# Better function to generate a ggplot manhattan plot
ggmanhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", "gray60"), logp = TRUE, title = "Manhattan Plot",
                        min_pval = -log10(0.01), chrlabs = c(1, 2, 3, 4, 5, 7, 9, 12, 15, 19), suggestiveline = -log10(1e-05), suggestLinecol = 'red',
                        genomewideline = -log10(5e-08), genomeLinecol = 'blue', highlight = NULL,annotateHighlight = FALSE) {
  
  # Ensure the data is as expected
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Ensure all columns are the correct class
  plot_data <- x %>%
    mutate(CHR = as.numeric(!!sym(chr)),
           BP = as.numeric(!!sym(bp)),
           P = as.numeric(!!sym(p)),
           SNP = !!sym(snp),
           logp = case_when(logp ~ -log10(P), T ~ P)) %>%
    arrange(CHR,BP) %>%
    group_by(CHR) %>%
    mutate(chr_length = max(BP)) %>%
    ungroup()
  
  sum_plot_data <- plot_data %>%
    select(CHR,chr_length) %>%
    unique() %>%
    mutate(chr_sum = lag(cumsum(chr_length), default = 0),
           centre = chr_length/2 + chr_sum)
  
  # Join all the data into a final plot
  fnal_plot_data <- sum_plot_data %>%
    right_join(plot_data) %>%
    mutate(fnal_pos = chr_sum + BP)
  
  # Set a minimum Pval to bother plotting
  significant_data <- fnal_plot_data %>% filter(logp >= min_pval)
  insignificant_boxes <- sum_plot_data %>%
    mutate(
      xmin = chr_sum,
      xmax = chr_sum + chr_length,
      ymin = 0,
      ymax = min_pval
    )
  
  # Get the labels for the chromosomes
  chromosome_labels <- ifelse(sum_plot_data$CHR %in% chrlabs, sum_plot_data$CHR, "")
  
  plot <- ggplot() +
    geom_point(data = significant_data, aes(x = fnal_pos, y = logp, color = as.factor(CHR)), alpha = 0.75, size = 1.3) +
    geom_rect(data = insignificant_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = as.factor(CHR)), alpha = 1) +
    scale_x_continuous(label = chromosome_labels, breaks = sum_plot_data$centre) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length.x = unit(0.2, "cm")) +
    labs(title = title, x = "Chromosome", y = expression(-log[10](italic(p))))
  
  
  # Add colourings
  if(!is.null(col)){
    rep_sum <- nrow(sum_plot_data)/2
    plot <- plot + 
      scale_color_manual(values = rep(col, rep_sum)) +
      scale_fill_manual(values = rep(col, rep_sum))
  }
  
  # Add in significance lines
  if (!is.null(suggestiveline)) {
    plot <- plot + geom_hline(yintercept = suggestiveline, color = suggestLinecol, linetype = "dashed")
  }
  if (!is.null(genomewideline)) {
    plot <- plot + geom_hline(yintercept = genomewideline, color = genomeLinecol)
  }
  
  # Recolour highlights
  if (!is.null(highlight)) {
    highlight_snps <- if (is.list(highlight)) highlight[[1]] else highlight
    highlight_clrs <- if (is.list(highlight)) highlight[[2]] else rep("green3", length(highlight))
    highlight_data <- fnal_plot_data %>% 
      filter(SNP %in% highlight_snps) %>%
      arrange(match(SNP, highlight_snps))
    plot <- plot + geom_point(data = highlight_data, aes(x = fnal_pos, y = logp),color = highlight_clrs, size = 2,alpha = 0.75)
    if (annotateHighlight) {
      plot <- plot + geom_text(data = highlight_data,aes(x = fnal_pos, y = logp, label = SNP),vjust = -1,color = "black", size = 3)
    }
  }
  return(plot)
}

# Function to generate a qq plot
ggqq <- function(x, p = "P", snp = "SNP", point_colour = "blue", diag_colour = "red", logp = TRUE, title = "QQ Plot",
                 highlight = NULL, annotateHighlight = FALSE, inflation = TRUE, min_pval = -log10(0.01),
                 suggestiveline = -log10(1e-05), suggestLinecol = 'red',
                 genomewideline = -log10(5e-08), genomeLinecol = 'blue'){
  
  # Ensure the data is as expected
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Calculate expected and observed -log10(p-values)
  n <- nrow(x)
  qq_data <- x %>%
    mutate(Expected = -log10((rank(as.numeric(!!sym(p))) - 0.5) / n),
           Observed = if(logp) -log10(!!sym(p)) else !!sym(p))
  
  # Separate data based on min_pval
  line_data <- qq_data %>% filter(Observed <= min_pval)
  point_data <- qq_data %>% filter(Observed > min_pval)
  
  # Create the plot
  plot <- ggplot() +
    geom_line(data = line_data, aes(x = Expected, y = Observed), color = point_colour, size = 1.3) +
    geom_point(data = point_data, aes(x = Expected, y = Observed), color = point_colour, alpha = 0.75, size = 1.3) +
    geom_abline(slope = 1, intercept = 0, color = diag_colour, linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = expression(Expected ~ -log[10](italic(p))), y = expression(Observed ~ -log[10](italic(p)))) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Highlight points if needed
  if (!is.null(highlight)) {
    highlight_snps <- if (is.list(highlight)) highlight[[1]] else highlight
    highlight_clrs <- if (is.list(highlight)) highlight[[2]] else rep("green3", length(highlight))
    highlight_data <- qq_data %>%
      filter(!!sym(snp) %in% highlight_snps) %>%
      arrange(match(!!sym(snp), highlight_snps))
    plot <- plot + geom_point(data = highlight_data, aes(x = Expected, y = Observed), color = highlight_clrs, size = 2, alpha = 0.75)
    if (annotateHighlight) {
      plot <- plot + geom_text(data = highlight_data, aes(x = Expected, y = Observed, label = !!sym(snp)), vjust = -1, color = "black", size = 3)
    }
  }
  
  # Add inflation factor if needed
  if(inflation){
    p_values <- qq_data$P
    chisq <- qchisq(1 - p_values, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    plot <- plot + annotate("text", x = 0.65, y = 0.2, label = paste("Inflation factor (λ) =", round(lambda, 3)), 
                            size = 4, hjust = 0)
  }
  
  # Add in significance lines
  if (!is.null(suggestiveline)) {
    plot <- plot + geom_hline(yintercept = suggestiveline, color = suggestLinecol, linetype = "dashed")
  }
  if (!is.null(genomewideline)) {
    plot <- plot + geom_hline(yintercept = genomewideline, color = genomeLinecol)
  }
  
  return(plot)
}