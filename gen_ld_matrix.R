# NCmisc::list.functions.in.file(filename = "~/Desktop/University/Biomedicine/Y3/Dissertation/PLINK_Tutorials/Simulate_Data_Test/gen_ld_matrix.R",
#                                                                                               alphabetic = TRUE)
library(circlize)
library(ComplexHeatmap)
library(stats)


# Uses method from: "https://bmcgenomdata.biomedcentral.com/articles/10.1186/s12863-020-0818-9#Equ3"

# Function to get the counts of genotypes
get_counts <- function(genotype_data,snp1 = 1,snp2 = 2){
  # Ensure the genotypes are factors with levels 0, 1, and 2
  genotype_data[[snp1]] <- factor(genotype_data[[snp1]], levels = 0:2)
  genotype_data[[snp2]] <- factor(genotype_data[[snp2]], levels = 0:2)
  
  # Create a table to count occurrences of each genotype combination
  genotype_table <- table(genotype_data[[snp1]], genotype_data[[snp2]])
  counts <- c(
    genotype_table["0", "0"], # AA BB
    genotype_table["0", "1"], # AA Bb
    genotype_table["0", "2"], # AA bb
    genotype_table["1", "0"], # Aa BB
    genotype_table["1", "1"], # Aa Bb
    genotype_table["1", "2"], # Aa bb
    genotype_table["2", "0"], # aa BB
    genotype_table["2", "1"], # aa Bb
    genotype_table["2", "2"]  # aa bb
  )
  return(counts)
}

# Log-liklihood function for estimating haplotype frequencies (this is from the paper)
loglik <- function(p_AB, p_Ab, p_aB, counts) {
  # Calculate the frequency of the fourth haplotype (ab) by subtracting the other three from 1
  p_ab <- 1 - p_AB - p_Ab - p_aB
  
  # Calculate the expected frequencies of each genotype combination based on haplotype frequencies
  # A small constant (1e-10) is added to avoid taking the log of zero. These are calculated based off HWE (Table 1)
  f_AABB <- p_AB^2 + 1e-10
  f_AABb <- 2*p_AB*p_Ab + 1e-10
  f_AAbb <- p_Ab^2 + 1e-10
  f_AaBB <- 2*p_AB*p_aB + 1e-10
  f_AaBb <- 2*(p_AB*p_ab + p_Ab*p_aB) + 1e-10
  f_Aabb <- 2*p_Ab*p_ab + 1e-10
  f_aaBB <- p_aB^2 + 1e-10
  f_aaBb <- 2*p_aB*p_ab + 1e-10
  f_aabb <- p_ab^2 + 1e-10
  
  # Calculate the log-likelihood by summing the product of genotype counts and log of expected frequencies.
  # This is made from equation 7 of the paper
  loglik <- counts[1]*log(f_AABB) + 
    counts[2]*log(f_AABb) + 
    counts[3]*log(f_AAbb) + 
    counts[4]*log(f_AaBB) + 
    counts[5]*log(f_AaBb) + 
    counts[6]*log(f_Aabb) + 
    counts[7]*log(f_aaBB) + 
    counts[8]*log(f_aaBb) + 
    counts[9]*log(f_aabb)
  
  return(loglik)
}

# Transformation function for haplotype frequencies
transform_freqs <- function(u, v, w) {
  # Transform the optimization parameters (u, v, w) into haplotype frequencies
  # This is based on equation 8 which was solved for pAb, pAb and paB instead
  p_AB <- u*v*w
  p_Ab <- u*v*(1-w)
  p_aB <- u*(1-v)
  
  # Return a vector of the first three haplotype frequencies
  return(c(p_AB, p_Ab, p_aB))
}

# Function to estimate haplotype frequencies using Constrained ML
estimate_haplotype_freqs <- function(genotype_data, snp1 = 1, snp2 = 2) {
  
  # Get genotype counts using the provided function
  counts <- get_counts(genotype_data, snp1, snp2)
  
  #' optim is a function in R that performs optimisation. It finds values of parameters that minimise or maximise a given
  #' objective function.par = the starting point of the parameters which are initial guesses for u,v and w respectively. 
  #' fn=function(x) -loglik(...) specifies the objective function to be minimised. In this case, we want to maximise the log-likelihood, 
  #' but since optim() minimises by default, we put a negative sign in front of loglik() to maximise it instead.Inside loglik(), 
  #' we have transform_freqs(x[1],x[2],x[3])[1], transform_freqs(x[1],x[2],x[3])[2], and transform_freqs(x[1],x[2],x[3])[3]. 
  #' These expressions transform the optimization parameters (u, v, w) into the haplotype frequencies (p_AB, p_Ab, p_aB) using the 
  #' transform_freqs() function. Method="L-BFGS-B" specifies the optimization algorithm to be used. L-BFGS-B is a widely used algorithm for 
  #' solving optimization problems with box constraints (i.e., lower and upper bounds on the parameters).
  #' When code is run, optim searches for values of u,v and w that maximise the loglik function, iteritively adjusting the parameter values until
  #' it finds the best solution or reaches a stopping criterion. The final result is stored in sol
  sol <- optim(par=c(0.3, 0.3, 0.3), fn=function(x) -loglik(transform_freqs(x[1],x[2],x[3])[1], # p_AB
                                                            transform_freqs(x[1],x[2],x[3])[2], # p_Ab
                                                            transform_freqs(x[1],x[2],x[3])[3], # p_aB
                                                            counts), 
               method="L-BFGS-B", 
               lower=c(0,0,0), # Lower bounds for the optimization parameters
               upper=c(1,1,1)  # Upper bounds for the optimization parameters
  )
  
  # Extract the estimated haplotype frequencies from the optimization result
  p_AB_hat <- transform_freqs(sol$par[1], sol$par[2], sol$par[3])[1]
  p_Ab_hat <- transform_freqs(sol$par[1], sol$par[2], sol$par[3])[2] 
  p_aB_hat <- transform_freqs(sol$par[1], sol$par[2], sol$par[3])[3]
  p_ab_hat <- 1 - p_AB_hat - p_Ab_hat - p_aB_hat
  
  return(list(p_AB = p_AB_hat, p_Ab = p_Ab_hat, p_aB = p_aB_hat, p_ab = p_ab_hat))
}

# Function to compute linkage disequilibrium (LD) measures from haplotype frequencies
compute_ld <- function(haplotype_freqs) {
  
  # Extract the haplotype frequencies from the input list
  p_AB <- haplotype_freqs$p_AB
  p_Ab <- haplotype_freqs$p_Ab
  p_aB <- haplotype_freqs$p_aB
  p_ab <- haplotype_freqs$p_ab
  
  # Calculate the allele frequencies
  p_A <- p_AB + p_Ab
  p_a <- p_aB + p_ab
  p_B <- p_AB + p_aB
  p_b <- p_Ab + p_ab
  
  # Calculate the D statistic (deviation from linkage equilibrium)
  D <- p_AB - p_A * p_B
  
  # Calculate the maximum and minimum possible values for D
  Dmax <- min(p_A * p_b, p_a * p_B)
  Dmin <- max(-p_A * p_B, -p_a * p_b)
  
  # Calculate the D' statistic (normalized D)
  if (D < 0) {
    D_prime <- D / Dmin
  } else {
    D_prime <- D / Dmax
  }
  
  # Calculate the r-squared statistic (square of the correlation coefficient)
  r_squared <- D^2 / (p_A * p_a * p_B * p_b)
  
  # Return a list of the calculated LD measures
  return(list(D = D, D_prime = D_prime, r_squared = r_squared))
}

# Function to estimate LD measures from genotype data
estimate_ld <- function(genotype_data, snp1 = 1, snp2 = 2) {
  # Estimate haplotype frequencies from the genotype data
  haplotype_freqs <- estimate_haplotype_freqs(genotype_data, snp1, snp2)
  # Compute LD measures using the estimated haplotype frequencies
  ld <- compute_ld(haplotype_freqs)
  return(ld)
}

# Calcualte r^2 directly 
calculate_ld <- function(snp_data, snp1, snp2) {
  ld <- estimate_ld(snp_data, snp1, snp2)
  return(ld$r_squared)
}

# Function to calculate an LD matrix from genotype data
gen_ld_matrix <- function(genotype_data){
  snp_count <- ncol(genotype_data)
  ld_matrix <- matrix(0, nrow = snp_count, ncol = snp_count)
  for (i in 1:(snp_count-1)) {
    for (j in (i+1):snp_count) {
      cat("\014")
      print(paste0("Comparing SNP ",i, " with SNP ",j))
      ld_matrix[i, j] <- calculate_ld(genotype_data, i, j)
      # Fill in flip side to save on computation
      ld_matrix[j, i] <- ld_matrix[i, j]
    }
  }
  diag(ld_matrix) <- 1
  rownames(ld_matrix) <- colnames(genotype_data)
  colnames(ld_matrix) <- colnames(genotype_data)
  
  # Generate the heatmap
  color_scale <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  ld_heatmap <- Heatmap(ld_matrix,
                        name = "LD",
                        col = color_scale,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        show_column_dend = FALSE,
                        show_row_dend = FALSE,
                        heatmap_legend_param = list(
                          title = "LD Value", at = c(0, 0.25, 0.5, 0.75, 1),
                          labels = c("0", "0.25", "0.5", "0.75", "1"))
  )
  
  return(list(ld_matrix = ld_matrix, ld_heatmap = ld_heatmap))
}




