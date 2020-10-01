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
n_neighbors = c(2,8,14,20) # local to global

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  #indx = sample(x = 1:nrow(lo_mat), size = 50000)
  set.seed(42)
  lo_umap = uwot::umap(X = lo_mat, n_threads = 8, n_neighbors = i, verbose = TRUE)
  saveRDS(object = lo_umap, file = paste0("data/lo_umap/lo_umap_", i, "nn.rds"))
}

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
