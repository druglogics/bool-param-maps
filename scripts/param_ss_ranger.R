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
ranger_res = ranger::ranger(x = lo_mat, y = ss_num, write.forest = FALSE,
  num.trees = 500, mtry = 16, num.threads = 8, importance = 'impurity', verbose = TRUE)

# Save result
saveRDS(object = ranger_res, file = "data/ranger_res.rds")
