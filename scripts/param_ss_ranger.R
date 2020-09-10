#####################################
# Random Forest with ranger package #
#####################################
library(dplyr)
library(tibble)
library(ranger)

# read data
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
lo_mat = readRDS(file = "data/lo_mat.rds")

# classification: `ss_num` as factor
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()

# Run ranger RF
set.seed(42)
model_samples_indexes = sample(x = 1:nrow(lo_mat), size = 4000000)
ranger_res = ranger::ranger(x = lo_mat[model_samples_indexes,],
  y = ss_num[model_samples_indexes],
  num.trees = 300, mtry = 16, #save.memory = TRUE,
  num.threads = 4, importance = 'impurity')

saveRDS(object = ranger_res, file = "data/ranger_res.rds")
