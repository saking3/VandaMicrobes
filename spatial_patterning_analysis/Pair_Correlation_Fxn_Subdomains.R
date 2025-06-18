library(spatstat)

# Parameters
sub_win_size <- 0.4
r_vals <- seq(0, 0.05, length.out = 100)
nsim <- 100
min_points <- 20

# Null hypothesis functions for PCF
pcf_csr <- function(r) rep(1, length(r))  # PCF for CSR is just 1

# Main analysis function
analyze_site_with_avg_pcf <- function(filename, output_subfolder) {
  cat("Processing:", filename, "\n")
  
  data <- read.csv(filename)
  pp <- ppp(data$X, data$Y, window = owin(c(min(data$X), max(data$X)), c(min(data$Y), max(data$Y))))
  cat("Planar point pattern:", summary(pp)$n, "points\n")
  
  bounds_x <- range(data$X)
  bounds_y <- range(data$Y)
  xseq <- seq(bounds_x[1], bounds_x[2] - sub_win_size, length.out = 5)
  yseq <- seq(bounds_y[1], bounds_y[2] - sub_win_size, length.out = 5)
  subwindows <- expand.grid(x = xseq, y = yseq)
  
  all_pcf_vals <- list()
  all_lower_bounds_pcf <- list()
  all_upper_bounds_pcf <- list()
  
  if (!dir.exists(output_subfolder)) dir.create(output_subfolder, recursive = TRUE)
  
  for (i in 1:nrow(subwindows)) {
    cat(paste0("Subwindow ", i, "/", nrow(subwindows), ": "))
    sw <- subwindows[i, ]
    win <- owin(c(sw$x, sw$x + sub_win_size), c(sw$y, sw$y + sub_win_size))
    
    sub_pp <- pp[win]
    n_points <- npoints(sub_pp)
    cat(n_points, "points\n")
    
    if (n_points < min_points) {
      cat("Skipping subwindow with too few points.\n")
      next
    }
    
    envelope_result <- envelope(sub_pp, pcf, nsim = nsim)
    r_vals <- envelope_result$r
    pcf_vals <- envelope_result$obs
    lower_bound_pcf <- envelope_result$lo
    upper_bound_pcf <- envelope_result$hi
    
    # Save individual subdomain plot
    sub_plot_filename <- file.path(output_subfolder, paste0("subdomain_", i, "_PCF.png"))
    png(sub_plot_filename, width = 800, height = 800)
    
    plot(r_vals, pcf_vals, type = "l", lwd = 2, col = "blue",
         main = paste("Subdomain", i, "-", n_points, "points"),
         xlab = "Distance (r)", ylab = "Pair Correlation Function (PCF)")
    
    polygon(c(r_vals, rev(r_vals)), c(rep(0, length(r_vals)), rev(upper_bound_pcf)),
            border = FALSE, col = rgb(116/255, 118/255, 120/255, 0.4))
    polygon(c(r_vals, rev(r_vals)), c(rep(0, length(r_vals)), rev(lower_bound_pcf)),
            border = FALSE, col = rgb(116/255, 118/255, 120/255, 0.4))
    
    lines(r_vals, pcf_csr(r_vals), col = "green", lty = 2, lwd = 2)
    dev.off()
    
    # Store results for averaging
    all_pcf_vals[[length(all_pcf_vals)+1]] <- pcf_vals
    all_lower_bounds_pcf[[length(all_lower_bounds_pcf)+1]] <- lower_bound_pcf
    all_upper_bounds_pcf[[length(all_upper_bounds_pcf)+1]] <- upper_bound_pcf
  }
  
  if (length(all_pcf_vals) == 0) {
    cat("No subwindows with sufficient points, skipping file.\n")
    return(NULL)
  }
  
  avg_pcf_vals <- Reduce("+", all_pcf_vals) / length(all_pcf_vals)
  avg_lower_bounds_pcf <- Reduce("+", all_lower_bounds_pcf) / length(all_lower_bounds_pcf)
  avg_upper_bounds_pcf <- Reduce("+", all_upper_bounds_pcf) / length(all_upper_bounds_pcf)
  
  # Save averaged plot
  png_filename <- file.path(output_subfolder, "averaged_pcf_with_envelope.png")
  png(png_filename, width = 800, height = 800)
  
  plot(r_vals, avg_pcf_vals, type = "l", lwd = 2, col = "blue", 
       main = paste("Averaged PCF for", basename(filename)), 
       xlab = "Distance (r)", ylab = "Pair Correlation Function (PCF)")
  
  polygon(c(r_vals, rev(r_vals)), c(rep(0, length(r_vals)), rev(avg_upper_bounds_pcf)), 
          border = FALSE, col = rgb(116/255, 118/255, 120/255, 0.4))
  polygon(c(r_vals, rev(r_vals)), c(rep(0, length(r_vals)), rev(avg_lower_bounds_pcf)), 
          border = FALSE, col = rgb(116/255, 118/255, 120/255, 0.4))
  
  lines(r_vals, pcf_csr(r_vals), col = "green", lty = 2, lwd = 2)
  dev.off()
}

# === RUN ANALYSIS ON ALL CSVs IN THE FOLDER ===

input_folder <- "G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/2_Data_Cleaning/3_Extracted_Maxima_and_Heights/Maxima_Points_CSVs"
output_base <- "G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/3_Statistics/PCF_Subdomain"

csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)

for (csv_file in csv_files) {
  file_base <- tools::file_path_sans_ext(basename(csv_file))
  output_folder <- file.path(output_base, file_base)
  analyze_site_with_avg_pcf(csv_file, output_folder)
}

