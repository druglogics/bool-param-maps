#################
# Random Forest #
#################

library(dplyr)
library(tibble)
library(usefun)
library(randomForest)

# read data
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
lo_mat = readRDS(file = "data/lo_mat.rds")

# classification: `ss_num` as factor
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()

# Tune `mtry`
set.seed(42)
tune_mtry_list = list()
for(i in 1:100) {
  print(i)
  model_samples_indexes = sample(x = 1:nrow(lo_mat), size = 10000)
  tune_mtry_list[[i]] = randomForest::tuneRF(
    x = lo_mat[model_samples_indexes,], y = ss_num[model_samples_indexes],
    ntreeTry = 500, stepFactor = 1.3, trace = TRUE, plot = FALSE)
}

# tidy up
tune_mtry_tbl_list = lapply(tune_mtry_list, tibble::as_tibble)
mtry_data = dplyr::bind_rows(tune_mtry_tbl_list)

# save result
saveRDS(object = mtry_data, file = "data/mtry_data.rds")

# Get Importance measures from Random Forest (random samples)
set.seed(42)
rf_imp_data = list()
for(i in 1:20) { # sample 20 x 100000 models
  print(i)
  # larger model samples/more data to train is always better (100000)
  model_samples_indexes = sample(x = 1:nrow(lo_mat), size = 100000)
  # the best `mtry` from above is between 14 and 18 => we choose 16
  rf_res = randomForest::randomForest(x = lo_mat[model_samples_indexes,],
    y = ss_num[model_samples_indexes], ntree = 500, mtry = 16, importance = TRUE,
    do.trace = TRUE)
  rf_imp_data[[i]] = rf_res$importance
}

# keep only the two importance scores
rf_imp_data = lapply(rf_imp_data, function(mat) {mat[,c("MeanDecreaseAccuracy", "MeanDecreaseGini")]})

# save result
saveRDS(object = rf_imp_data, file = "data/rf_imp_data.rds")
