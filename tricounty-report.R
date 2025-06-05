# Libraries --------------------------------------------------------------------
library(lubridate)
library(writexl)

# Define Helper Functions ------------------------------------------------------

#' Set up filesystem
#' 
#' `create_filesystem` creates the given folders if they don't already exist.
#' 
#' @param input Filepath. Folder where input file will be placed.
#' @param processed Filepath. Folder where processed input files are moved.
#' @param internal Filepath. Folder to hold internal reports.
#' @param public Filepath. Folder to hold public reports.
#' 
#' @returns NULL.
create_filesystem <- function(input, processed, internal, public) {
  # - Create folders if needed
  for (f in c(input, processed, internal, public)) {
    if (!dir.exists(f)) {
      dir.create(f)
    }
  }
}

#' Clear out old reports before generating new ones.
#' 
#' `clear_old_reports` deletes reports from previous runs and returns a list of 
#' the reports that were deleted.
#' 
#' @param i_folder Filepath. Folder containing internal reports.
#' @param p_folder Filepath. Folder containing public reports.
#' 
#' @returns The list of old reports that were cleared.
clear_old_reports <- function(i_folder, p_folder) {
  # - Remove old internal reports
  i_reports <- list(list.files(i_folder, full.names = TRUE))
  p_reports <- grepl(".xlsx|public_report_", 
                     list.files(p_folder, full.names = TRUE))
  p_reports <- list(list.files(p_folder, full.names = TRUE)[p_reports])
  
  old_reports <- c(i_reports, p_reports)
  
  do.call(file.remove, old_reports)
  
  old_reports
}

#' Write report CSV files
#' 
#' `write_report_csv` writes the given data to the specified folder with the 
#' given filename.
#' 
#' @param data Dataframe. Report data.
#' @param filename String. Report filename.
#' @param folder Filepath. Report destination folder.
#' 
#' @returns NULL.
write_report_csv <- function(data, filename, folder) {
  write.csv(data, file.path(folder, filename), row.names = FALSE)
}

#' Validate input EpiTrax data
#' 
#' 'validate_data' checks the data for expected columns and data types, removes
#' unneeded columns, and returns the resulting data.
#' 
#' @param data Dataframe. EpiTrax data to validate.
#' 
#' @returns The validated data with all unneeded columns removed.
validate_data <- function(data) {
  # Check column names
  expected_cols <- c("patient_mmwr_year", "patient_mmwr_week", "patient_disease")
  actual_cols <- colnames(data)
  
  if (!all(expected_cols %in% actual_cols)) {
    stop("The EpiTrax data is missing one of the following fields:\n\n\t",
         paste(expected_cols, collapse=", "), 
         "\n\nPlease add the missing fields to the file and try again.")
  }
  # Check column data types
  if (class(data$patient_mmwr_week) != "integer" ||
      class(data$patient_mmwr_year) != "integer" ||
      class(data$patient_disease) != "character") {
    stop("One or more columns in the EpiTrax dataset has an incorrect 
         data type:\n\n",
         "\t'patient_mmwr_week' should be of type 'integer'\n",
         "\t'patient_mmwr_year' should be of type 'integer'\n",
         "\t'patient_disease' should be of type 'character'\n",
         "\nPlease try again with a valid dataset."
    )
  }
  # Remove all columns we're not using
  # - Note this also rearranges the columns into the order of expected_cols
  data <- data[expected_cols]
  
  data
}

#' Format input EpiTrax data
#' 
#' 'format_week_num' formats the input EpiTrax dataset with month numbers 
#' using the field 'patient_mmwr_week' and filters rows older than five years.
#' 
#' @param data Dataframe. Data to format.
#' 
#' @returns The formatted data.
format_week_num <- function(data) {
  # Format data
  data$month <- with(data, month(
    ymd(patient_mmwr_year * 10000 + 0101) + 
      patient_mmwr_week * 7
  ))
  data$patient_mmwr_week <- NULL
  data$counts <- 1 # Makes easier to use aggregate()
  colnames(data) <- c("year", "disease", "month", "counts")
  # - Rearrange columns for easier debugging
  data <- data[c("disease", "month", "year", "counts")]
  
  # - Extract last years of data
  data <- with(data, data[year >= (max(year) - 5), ])
  
  data
}

#' Read in input EpiTrax data
#' 
#' 'read_epitrax_data' reads EpiTrax data from a CSV, validates and formats it,
#' then returns the data.
#' 
#' @param input_folder Filepath. Folder containing input data.
#' @param processed_folder Filepath. Optional folder to move the processed data
#' 
#' @returns The validated and formatted EpiTrax data from the input file.
read_epitrax_data <- function(input_folder, processed_folder = NULL) {
  # Get file name from input data folder
  fname <- list.files(input_folder)
  
  if (length(fname) == 0) {
    stop("No input file provided.")
  }
  
  # Check for only 1 file
  if (length(fname) > 1) {
    stop("Please include only 1 file in the '", input_folder, "' folder.")
  }
  
  # Check file has correct extension
  fpath <- file.path(input_folder, fname)
  if (!file.exists(fpath) || !grepl("\\.csv$", fpath)) {
    stop("Please add an EpiTrax data file (.csv) to the '", input_folder, 
         "' folder")
  }
  
  # Read data from file
  input_data <- read.csv(fpath, header = TRUE)
  
  # Validate and format data
  input_data <- validate_data(input_data)
  input_data <- format_week_num(input_data)
  
  # Move processed input file to processed_folder (if provided)
  if (!is.null(processed_folder)) {
    file.rename(fpath, file.path(processed_folder, fname))
  }
  
  # Return data from file
  input_data
}

#' Reshape data frame with each month as a separate column
#' 
#' 'reshape_monthly_wide' reshapes a given data frame with diseases for rows and
#' months for columns.
#' 
#' @param df Dataframe. Data to reshape with months as columns.
#'
#' @returns The reshaped data frame.
reshape_monthly_wide <- function(df) {
  m_df <- with(df, reshape(
    merge(
      df,
      expand.grid(
        disease = unique(disease),
        month = unique(month)
      ),
      all = TRUE
    ),
    direction = "wide",
    idvar = "disease",
    timevar = "month"
  ))
  # - Set NA values to 0
  m_df[is.na(m_df)] <- 0
  # - Update column names to more human-readable format
  colnames(m_df) <- c("disease", month.abb[1:(ncol(m_df) - 1)])
  
  m_df
}

#' Get the public disease list
#' 
#' 'get_public_disease_list' reads the public list from a given CSV file or uses
#' the default diseases if the file doesn't exist.
#' 
#' The provided public disease list file must contain two columns that map the 
#' EpiTrax disease name to a public-facing name for the public report.
#' @param filepath Filepath. Public disease list CSV file.
#' @param default_diseases String vector. List of default diseases to use if the
#' above file doesn't exist.
#' 
#' @returns A dataframe containing the diseases to include in the public report 
#' and the name to use for each disease in the public report.
get_public_disease_list <- function(filepath, default_diseases) {
  
  if (file.exists(filepath)) {
    
    d_list <- read.csv(filepath, header = TRUE)
    
    # Validate file
    if (is.null(d_list$EpiTrax_name) || is.null(d_list$Public_name)) {
      stop("File '", filepath, "' is incorrectly formatted. Please use the 
           column names: 'EpiTrax_name' and 'Public_name'.")
    }
    
    d_list
    
  } else {
    # If the file doesn't exist, use the list of diseases in the input data
    warning("File '", filepath, "' not found. Using list from EpiTrax input
            dataset instead.")
    
    default_diseases <- sort(default_diseases)
    
    d_list <- data.frame(
      EpiTrax_name = default_diseases,
      Public_name = default_diseases
    )
    
    d_list
  }
}

#' Create a public report
#' 
#' 'create_public_report' creates a public report for the given month.
#' 
#' @param cases Dataframe. Disease case counts for each month and year.
#' @param avgs Dataframe. Disease case count averages for each month.
#' @param d_list Dataframe. List of diseases to use for the report.
#' @param m Integer. The report month.
#' @param y Integer. The report year.
#' @param r_folder Filepath. Destination folder for the public report.
#' 
#' @returns List containing the report name and data.
create_public_report <- function(cases, avgs, d_list, m, y, r_folder) {
  
  m_counts <- with(cases, cases[year == y & month == m, c("disease", "counts")])
  
  # - Only take the rows with data in the final report
  m_counts <- subset(m_counts, disease %in% avgs$disease)
  
  # - Create the report data frame initializing the Current_Rate column to 0
  m_report <- data.frame(
    Disease = avgs$disease,
    Current_Rate = 0, 
    Historical_Rate = round(avgs[m + 1], digits = 1)
  )
  
  # - Update the Current_Rate column with values from m_counts
  for (i in 1:length(m_counts$disease)) {
    d <- m_counts$disease[i]
    m_report[m_report$Disease == d, ]$Current_Rate <- m_counts$counts[i]
  }
  
  # - Add Trends column
  m_report$Trend <- mapply(function(x, y) {
    ifelse(x > y, "↑", ifelse(x < y, "↓", "→"))
  }, m_report$Current_Rate, m_report[[3]])
  
  # - Wait until final step to convert disease names to public-facing versions
  m_report$Disease <- d_list$Public_name
  
  # - Write to CSV file
  r_name <- paste0("public_report_", colnames(m_report)[3], report_year)
  write_report_csv(m_report, paste0(r_name, ".csv"), r_folder)
  
  list("name" = r_name, "report" = m_report)
}


# Set up file system -----------------------------------------------------------
input_data_folder <- "input_epitrax_data"
processed_data_folder <- "processed_epitrax_data"
internal_reports_folder <- "internal_reports"
public_reports_folder <- "public_reports"

xl_files <- list() # Internal reports to combine into single .xlsx file

create_filesystem(
  input = input_data_folder,
  processed = processed_data_folder,
  internal = internal_reports_folder,
  public = public_reports_folder
)

clear_old_reports(internal_reports_folder, public_reports_folder)


# Read in EpiTrax data ---------------------------------------------------------
epitrax_data <- read_epitrax_data(input_data_folder)
epitrax_data_yrs <- sort(unique(epitrax_data$year))
epitrax_data_diseases <- unique(epitrax_data$disease)
report_year <- max(epitrax_data_yrs)


# Compute annual counts for each disease ---------------------------------------
annual_counts <- aggregate(counts ~ disease + year, 
                           data = epitrax_data, 
                           FUN = sum)
# - Reshape data to use years as columns and diseases as rows
annual_counts <- with(annual_counts, reshape(
  merge(
    annual_counts,
    expand.grid(
      disease = unique(disease),
      year = unique(year)
    ),
    all = TRUE
  ),
  direction = "wide",
  idvar = "disease",
  timevar = "year"
))
# - Set NA values to 0
annual_counts[is.na(annual_counts)] <- 0
# - Update column names to more human-readable format
colnames(annual_counts) <- c("disease", epitrax_data_yrs)
# - Write to CSV
write_report_csv(annual_counts, "annual_counts.csv", internal_reports_folder)
# - Add to Excel List
xl_files[["annual_counts"]] <- annual_counts


# Compute monthly counts for each year -----------------------------------------
month_counts <- aggregate(counts ~ disease + year + month, 
                          data = epitrax_data, 
                          FUN = sum)

for (y in epitrax_data_yrs) {
  # - Extract counts for given year
  m_df <- month_counts[month_counts$year == y, ]
  
  # - Remove year column (don't want to save it to CSV)
  m_df$year <- NULL
  
  # - Reshape data to use months as columns and disease as rows
  m_df <- reshape_monthly_wide(m_df)
  
  # - Write to CSV
  fname <- paste0("monthly_counts_", y)
  write_report_csv(m_df, paste0(fname, ".csv"), internal_reports_folder)
  
  # - Add to Excel List
  xl_files[[fname]] = m_df
}


# Compute monthly averages for all years except current year -------------------
# - Extract all previous years
epitrax_data_prev_yrs <- epitrax_data[epitrax_data$year != report_year,]
num_yrs <- length(unique(epitrax_data_prev_yrs$year))

# - Compute average counts for each month
monthly_avgs <- aggregate(counts ~ disease + month, 
                          data = epitrax_data_prev_yrs, 
                          FUN = sum)

monthly_avgs$counts <- monthly_avgs$counts / num_yrs

# - Reshape data to use months as columns and disease as rows
monthly_avgs <- reshape_monthly_wide(monthly_avgs)

# - Write to CSV
avgs_fname <- with(epitrax_data_prev_yrs,
                   paste0("monthly_avgs_", min(year), "-", max(year), ".csv"))
write_report_csv(monthly_avgs, avgs_fname, internal_reports_folder)

# - Add to Excel List
xl_files[["monthly_avgs"]] <- monthly_avgs

# - Combine internal reports into single .xlsx file
write_xlsx(xl_files, 
           file.path(internal_reports_folder, "internal_reports.xlsx"))


# Prepare Public Report --------------------------------------------------------
diseases <- get_public_disease_list(
  file.path(public_reports_folder, "disease_list_for_public_report.csv"),
  default_diseases = epitrax_data_diseases
)

# - Remove rows from monthly_avgs that aren't going into the public report
monthly_avgs <- subset(monthly_avgs, disease %in% diseases$EpiTrax_name)

# - Get diseases from public report list that weren't in the EpiTrax data
missing_diseases <- subset(diseases, !(EpiTrax_name %in% monthly_avgs$disease))

# If there are any missing diseases, add them
if (length(missing_diseases$EpiTrax_name) > 0) {
  # - Fill the missing diseases in with avg = 0
  missing_avgs <- data.frame(
    disease = missing_diseases$EpiTrax_name
  )
  
  missing_cols <- colnames(monthly_avgs)[2:length(colnames(monthly_avgs))]
  missing_avgs[, missing_cols] <- 0.0
  
  # - Combine monthly_avgs with missing_avgs
  monthly_avgs <- rbind(monthly_avgs, missing_avgs)
  
  # - Sort alphabetically so missing diseases are correctly placed
  monthly_avgs <- monthly_avgs[order(monthly_avgs$disease),]
}

# - Find the previous month of the report year
report_month <- max(month_counts[month_counts$year == report_year, ]$month) - 1

# - Create report table for most recent month and for 1 and 2 months prior
xl_files <- list()

r <- create_public_report(
  cases = month_counts, 
  avgs = monthly_avgs, 
  d_list = diseases, 
  m = report_month, 
  y = report_year,
  r_folder = public_reports_folder
)
xl_files[[r[["name"]]]] <- r[["report"]]

r <- create_public_report(
  cases = month_counts, 
  avgs = monthly_avgs, 
  d_list = diseases, 
  m = report_month - 1, 
  y = report_year,
  r_folder = public_reports_folder
)
xl_files[[r[["name"]]]] <- r[["report"]]

r <- create_public_report(
  cases = month_counts, 
  avgs = monthly_avgs, 
  d_list = diseases, 
  m = report_month - 2, 
  y = report_year,
  r_folder = public_reports_folder
)
xl_files[[r[["name"]]]] <- r[["report"]]

# - Combine public reports into single .xlsx file
write_xlsx(xl_files, 
           file.path(public_reports_folder, "public_reports.xlsx"))
