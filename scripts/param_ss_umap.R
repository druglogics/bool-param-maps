################################################################
# Use UMAP on the whole CASCADE 1.0 dataset (parameterization) #
################################################################
library(dplyr)
library(tibble)
library(uwot)

# NOTE: we suggest having a lot of memory to run this script! (64GB should be enough)

# read data
models_ss_stats = readRDS(file = "data/models_ss_stats.rds") # see `count_models_ss.R`
lo_mat = readRDS(file = "data/lo_mat.rds") # see `get_lo_mat.R`
# order of models is alright, since `lo_mat` was produced from `models_ss_stats`

################
# Unsupervised #
################

# set neighbors
n_neighbors = c(8,14,20) # local to global

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  set.seed(42)
  lo_umap = uwot::umap(X = lo_mat, n_threads = 8, n_neighbors = i, verbose = TRUE)
  saveRDS(object = lo_umap, file = paste0("data/lo_umap/lo_umap_", i, "nn.rds"))
}

# run once with hamming distance
set.seed(42)
lo_umap = uwot::umap(X = lo_mat, n_threads = 8, n_neighbors = 14,
  metric = "hamming", min_dist = 0.5, verbose = TRUE)
saveRDS(object = lo_umap, file = paste0("data/lo_umap/lo_umap_14nn_ham_0_5_min_dist.rds"))

##############
# Supervised #
##############

# classification: `ss_num` as factor
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  set.seed(42)
  lo_sumap = uwot::umap(X = lo_mat, y = ss_num, n_threads = 8, n_neighbors = i, verbose = TRUE)
  saveRDS(object = lo_sumap, file = paste0("data/lo_umap/lo_sumap_", i, "nn.rds"))
}

# run once with `euclidean` distance, 14 `n_neighbors` and `min_dist` = 0.3
set.seed(42)
lo_sumap = uwot::umap(X = lo_mat, y = ss_num, n_threads = 8, n_neighbors = 14,
  min_dist = 0.3, verbose = TRUE)
saveRDS(object = lo_sumap, file = "data/lo_umap/lo_sumap_14nn_0_3_min_dist.rds")
