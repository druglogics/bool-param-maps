########
# UMAP #
########

library(dplyr)
library(tibble)
# Use `uwot` library (fast, multiple thread support)
library(uwot)

# NOTE: we suggest having a lot of memory to run this script! (64GB should be enough)
# Memory problems was the reason we do not run UMAP on the whole dataset

# read data
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
lo_mat = readRDS(file = "data/lo_mat.rds")

# classification: `ss_num` as factor
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()

# choose seed and model sample for unsupervised case
set.seed(42)
n_samples = 6000000
model_samples_indexes = sample(x = 1:nrow(lo_mat), size = n_samples)
saveRDS(object = model_samples_indexes, file = "data/indexes.rds")

lo_umap = uwot::umap(
  X = lo_mat[model_samples_indexes,],
  metric = "manhattan", n_threads = 4) # put more `n_threads` if you have them!
saveRDS(object = lo_umap, file = "data/lo_umap.rds")

# choose seed and model sample for the supervised case
set.seed(42)
n_samples_sup = nrow(lo_mat)/3 # 1/3 of the dataset, ~ 3 million rows
model_samples_indexes_sup = sample(x = 1:nrow(lo_mat), size = n_samples_sup)
saveRDS(object = model_samples_indexes_sup, file = "data/indexes_sup.rds")

lo_sumap = uwot::umap(
  X = lo_mat[model_samples_indexes_sup,],
  y = ss_num[model_samples_indexes_sup],
  metric = "manhattan", n_threads = 4) # put more `n_threads` if you have them!
saveRDS(object = lo_sumap, file = "data/lo_sumap.rds")
