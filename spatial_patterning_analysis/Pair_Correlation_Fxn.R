# Load required library
library(spatstat)

# Set the working directory to the folder containing your CSV files
setwd('G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/2_Data_Cleaning/3_Extracted_Maxima_and_Heights/Maxima_Points_CSVs')

# Define the output folder for saving plots
output_folder <- "G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/2_Data_Cleaning/3_Extracted_Maxima_and_Heights/Maxima_Points_CSVs/Pairwise_Correlation_Function"

# Create the output folder if it does not exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# List all CSV files in the folder
file_list <- list.files(pattern = "\\.csv$")

# Loop through each file and calculate the pairwise correlation function with envelopes
for (file_name in file_list) {
  # Read the CSV file
  data <- read.csv(file_name, header = TRUE)
  
  # Dynamically calculate the center and radius for the circular site
  x_center <- mean(data[, 1])
  y_center <- mean(data[, 2])
  max_distance <- max(sqrt((data[, 1] - x_center)^2 + (data[, 2] - y_center)^2))
  
  # Define a circular window using the center and radius
  circular_window <- disc(radius = max_distance, centre = c(x_center, y_center))
  
  # Convert to point pattern object using the circular window
  data_ppp <- as.ppp(data[, c(1, 2)], W = circular_window)
  
  # Calculate the pairwise correlation function
  pcf_result <- pcf(data_ppp, correction = "isotropic", nlarge = Inf)
  
  # Generate an envelope for the PCF
  envelope_result <- envelope(data_ppp, fun = pcf, nsim = 100, correction = "isotropic", nlarge = Inf)
  
  # Extract data for plotting
  r <- envelope_result$r
  g <- envelope_result$obs
  g_low <- envelope_result$lo
  g_high <- envelope_result$hi
  
  # Generate the PNG filename based on the CSV file name
  png_filename <- file.path(output_folder, paste0(tools::file_path_sans_ext(file_name), "_pcf_with_envelope.png"))
  
  # Save the plot as a PNG
  png(png_filename, width = 800, height = 800)  # Open PNG device
  
  # Plot the pairwise correlation function with envelopes
  plot(r, g, type = "l", col = "blue", lwd = 2, 
       xlab = "Distance (m)", ylab = "Pairwise Correlation Function (g(r))",
       main = paste("Pairwise Correlation Function:", file_name))
  
  # Add the envelopes (light gray)
  polygon(c(r, rev(r)), c(g_high, rev(g_low)), border = NA, col = rgb(0.7, 0.7, 0.7, 0.4))
  
  # Re-plot the observed PCF on top
  lines(r, g, col = "blue", lwd = 2)
  
  dev.off()  # Close the PNG device
}


