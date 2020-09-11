# Use this script to create the link operator data for the
# models that have 1 stable state
library(dplyr)
library(stringr)
library(usefun)
library(emba)

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

# For the order and names of the link operator nodes (columns)
# (see `get_node_stats.R` script)
node_stats = readRDS(file = "data/node_stats.rds")
lo_nodes = node_stats %>% pull(node)
num_bits = length(lo_nodes)

# Make the model link operator data (binary encoding: AND-NOT => 0, OR-NOT => 1)
# from the corresponding numerical model representation of the stable state data!
lo_data = sapply(rownames(ss_data), function(model_name) {
  model_name = stringr::str_remove(model_name, "network_") %>% as.integer()
  stringr::str_split(
    string = usefun::dec_to_bin(decimal_num = model_name, bits = num_bits),
    simplify = TRUE, # return matrix
    pattern = "")} %>%
    as.integer()) %>%
  t() # transpose result => columns are now the nodes

# name the columns as the link-operator nodes
colnames(lo_data) = lo_nodes

# save result
saveRDS(object = lo_data %>% as.data.frame(), file = "data/lo_data.rds")
