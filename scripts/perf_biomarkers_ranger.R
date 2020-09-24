# Find important nodes for MCC performance, using the 1
# stable state model dataset and random forests via `ranger`
library(dplyr)
library(tibble)
library(ranger)

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

# Read the link operator data (see `get_1ss_lo_data.R`)
lo_data = readRDS(file = "data/lo_data.rds")

# check the model names have the same order
stopifnot(all(rownames(ss_data) == rownames(lo_data)))

## MCC forests :)

# get models MCC scores
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds") # see `emba_analysis.R`
models_mcc = mcc_res$models.mcc

# check order
stopifnot(all(rownames(ss_data) == names(models_mcc)))
stopifnot(all(rownames(lo_data) == names(models_mcc)))

# drop the model names and convert to matrix for faster simulation
models_mcc = unname(models_mcc)
model_names = rownames(lo_data)
rownames(ss_data) = NULL
rownames(lo_data) = NULL
ss_data = as.matrix(ss_data)
lo_data = as.matrix(lo_data)

# run ranger for stable state data
set.seed(42)
ss_mcc_ranger_res = ranger::ranger(x = ss_data, y = models_mcc, write.forest = FALSE,
 num.trees = 500, mtry = 16, num.threads = 4, importance = 'impurity',
 verbose = TRUE)
# save result
saveRDS(object = ss_mcc_ranger_res, file = "data/ss_mcc_ranger_res.rds")

# run ranger for link operator state data
set.seed(42)
lo_mcc_ranger_res = ranger::ranger(x = lo_data, y = models_mcc, write.forest = FALSE,
  num.trees = 500, mtry = 16, num.threads = 4, importance = 'impurity',
  verbose = TRUE)
# save result
saveRDS(object = lo_mcc_ranger_res, file = "data/lo_mcc_ranger_res.rds")
