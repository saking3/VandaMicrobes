# Load required library
library(spatstat)

# Set the working directory to the folder containing your CSV files
setwd('G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/2_Data_Cleaning/3_Extracted_Maxima_and_Heights/Maxima_Points_CSVs')

output_folder <- "G:/Official_Vanda_Organizing/Spatial_Morphologic_Data/2_Data_Cleaning/3_Extracted_Maxima_and_Heights/Maxima_Points_CSVs/Ripleys_L"

# List all CSV files in the folder
file_list <- list.files(pattern = "\\.csv$")

# Loop through each file and calculate Ripley's L function for circular sites
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
  
  # Calculate Ripley's K function
  kest <- Kest(data_ppp, correction = "isotropic", nlarge = Inf)
  
  # Create envelope for Ripley's K
  envelope_result <- envelope(data_ppp, Kest, nlarge = Inf, nsim = 100)
  
  # Prepare data for Ripley's L plot
  x <- envelope_result$r
  L <- sqrt(envelope_result$lo / pi) - envelope_result$r
  U <- sqrt(envelope_result$hi / pi) - envelope_result$r
  
  # Generate the PNG filename based on the CSV file name
  png_filename <- file.path(output_folder, paste0(tools::file_path_sans_ext(file_name), "_ripleysL.png"))
  
  # Save the plot as a PNG
  png(png_filename, width = 800, height = 800)  # Open PNG device
  
  # Plot Ripley's L function
  plot(kest, sqrt(./pi) - r ~ r, xaxs = "i", yaxs = "i", yaxt = "n", xaxt = "n",
       ylab = NA, xlab = NA, legend = FALSE, lwd = 2, main = file_name)
  axis(1, tcl = 0.4, lwd.ticks = 2, mgp = c(1, 0.4, 0))
  mtext(side = 1, line = 1.9, text = expression(bold("Distance (m)")), cex = 1.3)
  mtext(side = 2, line = 1.9, text = expression(bold("Ripley's L")), cex = 1.3)
  axis(2, tcl = 0.4, lwd.ticks = 2, mgp = c(1, 0.4, 0))
  polygon(c(x, rev(x)), c(L, rev(U)), border = FALSE, col = rgb(116 / 255, 118 / 255, 120 / 255, 0.4))
  
  dev.off()  # Close the PNG device
}

