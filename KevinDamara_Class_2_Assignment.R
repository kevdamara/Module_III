# Assignment
height <- c(1.75, 1.76, 1.82, 1.67)
c(68, 78, 85, 75) -> weight
smoking_status = c("Yes", "No", "No", "Yes")

# Arithmetic
BMI <- weight/(height^2)
BMI

# Comparisons
BMI > 25
BMI < 18.5
height >= 1.75
weight <= 65
smoking_status == "No"
smoking_status != "No"

# Logical combos
(BMI > 25) & (smoking_status == "Yes")
(BMI > 25) | (smoking_status == "Yes")
!(smoking_status == "Yes")
# Vectors
num_vec <- c(1, 2, 3, 4)
chrc_vector <- c("gene1", "gene2", "gene3")
logical_vector <- c(TRUE, FALSE, TRUE)
mix_vector <- c("gene1", 1, "gene2", 2) # coerces to character
num_vec[2]; num_vec[2:4]

# Lists
all_vectors <- list(num_vec=num_vec, chrc=chrc_vector, logi=logical_vector)
all_vectors[[2]]       # second element
all_vectors$chrc       # by name

# Matrices
my_matrix <- matrix(1:9, nrow=3, ncol=3)                 # column-wise
my_matrix_byrow <- matrix(1:9, nrow=3, ncol=3, byrow=TRUE) # row-wise
my_matrix_byrow[2,3]; my_matrix_byrow[2,]

# Data frame
data <- data.frame(
  patient_id = c("P1", "P2", "P3"),
  age        = c(65, 78, NA),
  diagnosis  = c("cancer", "diabetes", "cancer"),
  stringsAsFactors = FALSE
)

str(data); head(data); tail(data, 2); dim(data); names(data)
data$patient_id
data[1:2, c(1,3)]
data$new_column <- c(1,2,3)

# Missing values
is.na(data); sum(is.na(data)); colSums(is.na(data)); rowSums(is.na(data))

clean_data1  <- na.omit(data)                            # drop any row with NA
clean_data_2 <- data[, colSums(is.na(data))==0]          # keep columns with no NA

# Safer replacements (avoid coercing the whole data.frame!)
clean_data_3 <- data;  clean_data_3$age[is.na(clean_data_3$age)] <- 0
clean_data_4 <- data;  clean_data_4$age[is.na(clean_data_4$age)] <- mean(data$age, na.rm=TRUE)

# Basic
calculate_BMI <- function(weight, height) {
  weight / (height^2)
}

calculate_BMI(weight = 60, height = 1.75)
calculate_BMI(weight = c(68,78,85,75), height = c(1.75,1.76,1.82,1.67))

# Default argument
calculate_BMI <- function(weight, height = 1.75) {
  weight / (height^2)
}
calculate_BMI(60)  # uses default height

# Lazy evaluation (unused args allowed)
calculate_BMI <- function(weight, height, age) {
  weight / (height^2)
}
calculate_BMI(60, 1.65)  # 'age' ignored

# Make sure folders exist
input_dir  <- "Raw_Data"
output_dir <- "Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# files to process
files_to_process <- c("BMI_data_1.csv", "BMI_data_2.csv")

# helper (BMI)
calculate_BMI <- function(weight, height) weight / (height^2)

result_list <- list()

for (file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n")
  input_file_path <- file.path(input_dir, file_name)
  
  if (!file.exists(input_file_path)) {
    cat("  ⚠️ File not found:", input_file_path, "\n"); next
  }
  
  data <- read.csv(input_file_path, header = TRUE)
  cat("  File imported. Checking missing values...\n")
  
  if ("height" %in% names(data)) {
    miss_h <- sum(is.na(data$height))
    cat("  Missing 'height':", miss_h, "\n")
    data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
  }
  
  if ("weight" %in% names(data)) {
    miss_w <- sum(is.na(data$weight))
    cat("  Missing 'weight':", miss_w, "\n")
    data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
  }
  
  data$bmi <- calculate_BMI(data$weight, data$height)
  cat("  ✅ BMI calculated.\n")
  
  result_list[[file_name]] <- data
  
  output_file_path <- file.path(output_dir, paste0("BMI_results_", file_name))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("  💾 Saved:", output_file_path, "\n")
}

# convenient objects
results_1 <- result_list[["BMI_data_1.csv"]]
results_2 <- result_list[["BMI_data_2.csv"]]

# Setup
input_dir  <- "Raw_Data"
output_dir <- "Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

deg_files <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")  # match exact names

# 1) classification function
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1           # treat missing padj as non-significant
  if (!is.na(logFC) && logFC >  1 && padj < 0.05) return("Upregulated")
  if (!is.na(logFC) && logFC < -1 && padj < 0.05) return("Downregulated")
  return("Not_Significant")
}

# 2) loop over files
for (fn in deg_files) {
  cat("\nProcessing:", fn, "\n")
  path_in <- file.path(input_dir, fn)
  
  if (!file.exists(path_in)) { cat("  ⚠️ Missing:", path_in, "\n"); next }
  
  df <- read.csv(path_in, header = TRUE)
  
  # basic checks
  needed <- c("Gene_Id", "padj", "logFC")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols)) {
    stop("File ", fn, " is missing columns: ", paste(missing_cols, collapse=", "))
  }
  
  # 3) replace missing padj with 1
  na_p <- sum(is.na(df$padj))
  if (na_p > 0) cat("  Replacing", na_p, "NA padj with 1\n")
  df$padj[is.na(df$padj)] <- 1
  
  # 4) add 'status' column using vectorized mapply
  df$status <- mapply(classify_gene, df$logFC, df$padj)
  
  # 5) summaries
  sig_mask  <- df$padj < 0.05
  up_mask   <- df$status == "Upregulated"
  down_mask <- df$status == "Downregulated"
  
  cat("  Summary (table of status):\n")
  print(table(df$status))
  cat("  Significant (padj < 0.05):", sum(sig_mask), "\n")
  cat("  Upregulated:              ", sum(up_mask), "\n")
  cat("  Downregulated:            ", sum(down_mask), "\n")
  
  # 6) save output
  out_name <- paste0("DEGs_processed_", fn)
  path_out <- file.path(output_dir, out_name)
  write.csv(df, path_out, row.names = FALSE)
  cat("  💾 Saved:", path_out, "\n")
}

save.image("KevinDamara_Class_2_Assignment.RData")

Results/DEGs_processed_DEGs_Data_1.csv
Results/DEGs_processed_DEGs_Data_2.csv

list.files("Results")
list.files("Results", pattern = "DEGs_processed")
df1 <- read.csv("Results/DEGs_processed_DEGs_Data_1.csv")
head(df1)
file.exists("Results/DEGs_processed_DEGs_Data_1.csv")
df1 <- read.csv("Results/DEGs_processed_DEGs_Data_1.csv")
head(df1)

df1 <- read.csv("Results/DEGs_processed_DEGs_Data_1.csv")
head(df1)

df2 <- read.csv("Results/DEGs_processed_DEGs_Data_2.csv")
head(df2)

save.image("KevinDamara_Class_2_Assignment.RData")


