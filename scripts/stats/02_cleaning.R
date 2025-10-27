source("scripts/preprocessing/00_setup.R")

# 1. Load required packages and data
# Load the feature table from data/processed

# 2. Check the structure of the table
# View number of rows and columns, column names, and missing values

# 3. Remove technical samples
# Delete columns with names like MM8, QC, ACN

# 4. Check missing values
# Count NAs per feature (row) and per sample (column)

# 5. Filter out low-quality features
# Remove features with too many missing values (e.g., >30%)

# 6. Save the cleaned table
# Export the cleaned data as feature_table_neg_clean.csv

# 7. Compare before and after cleaning
# Check how many features and NAs were before vs after