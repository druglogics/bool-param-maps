library(dplyr)
library(stringr)
library(usefun)

models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
node_stats = readRDS(file = "data/node_stats.rds")
lo_nodes = node_stats %>% pull(node)
num_bits = length(lo_nodes)

# make data wider by transforming the model name from decimal
# to it's binary link operator representation
lo_mat = sapply(models_ss_stats %>% pull(model_number), function(model_name) {
  stringr::str_split(
    string = usefun::dec_to_bin(decimal_num = model_name, bits = num_bits),
    simplify = TRUE, # return matrix
    pattern = "")} %>%
    as.integer()) %>%
  t() # transpose result => columns are now the nodes

# name the columns as the link-operator nodes
colnames(lo_mat) = lo_nodes

# save data
# `lo_mat` has 23 columns for each node, describing all possible
# link-operator parameterizations (from 0...0 to 1...1). As such,
# it has 2^23 rows.
saveRDS(object = lo_mat, file = "data/lo_mat.rds")