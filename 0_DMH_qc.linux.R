

# Load defined functions
source("defined_functions.R")

# Read input parameters from a TXT file
read_input_params <- function(txt_file) {
  params <- list()
  con <- file(txt_file, "r")
  lines <- readLines(con)
  close(con)
  params$dir_path <- lines
  return(params)
}

# Main function to perform QC
perform_qc <- function(dir_paths, output_dir) {
  projects <- basename(dir_paths)
  for (j in seq_along(dir_paths)) {
    output_path <-
      file.path(output_dir, paste0(projects[j], ".output.log"))
    sink(file = output_path, type = "output")
    print(j)
    cellranger_qc(dir_paths[j], projects[j], output_dir)
    sink(type = "output")
  }
}

# Read input parameters from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

# Set working directory to control output location
setwd(output_dir)

# Read input parameters
input_params <- read_input_params(input_file)

# Extract data from input parameters
dir_paths <- input_params$dir_path

# Perform QC
start_time_qc <- Sys.time()
perform_qc(dir_paths, output_dir)
end_time_qc <- Sys.time()
print(end_time_qc - start_time_qc)
gc()
