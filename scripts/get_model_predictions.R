library(dplyr)
library(stringr)
library(emba)

# point to directory with the Drabme results from the Zenodo dataset
data_dir = "/home/john/tmp/balance_paper/cascade_1.0_hsa_fixpoints"
res_dirs = list.dirs(path = data_dir, recursive = FALSE)

model_pred_data = list()
index = 1
for (res_dir in res_dirs) {
  print(paste0("Reading directory: ", res_dir, " (", index, ")"))
  model_pred_file = stringr::str_subset(
    string = list.files(res_dir, full.names = TRUE),
    pattern = "model_predictions.tab")
  model_pred_data[[index]] = emba::get_model_predictions(model.predictions.file = model_pred_file)
  index = index + 1
}

res = dplyr::bind_rows(model_pred_data)

# read the stable state data (only models who have 1 stable state are included)
ss_data = readRDS(file = "data/ss_data.rds") # see `get_ss_data.R`

# keep the model predictions only for the 1 stable state models (and in the same order)
res = res[rownames(ss_data),]

# `model_predictions` is a data.frame where each row refers to a link operator
# mutated boolean model and the columns are the 21 drug combinations that were
# tested for synergy using the HSA rule for each model.
# Each value in the data.frame is either 0 (antagonistic), 1 (synergistic) or NA
# (couldn't perform assessment of synergy because of lack of attractors)
saveRDS(res, file = "data/model_predictions.rds")
