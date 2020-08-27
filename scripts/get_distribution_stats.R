library(dplyr)
library(tibble)
library(emba)

# read the models stable state data
res_ss_df = readRDS(file = "data/1ss_model_data.rds")

# point to the `CASCADE_1.0.sif` topology file in the Zenodo dataset
edge_mat = emba::get_edges_from_topology_file(topology.file = "/media/disk/abmlog/CASCADE_1.0.sif")
edge_tbl = edge_mat %>% as_tibble()

node_names = colnames(res_ss_df)

# calculate in-degrees
in_degree = sapply(node_names, function(node) { edge_tbl %>% filter(target == node) %>% nrow() })

# calculate out-degress
out_degree = sapply(node_names, function(node) { edge_tbl %>% filter(source == node) %>% nrow() })

# degree distribution stats for CASCADE 1.0
dd_stats = dplyr::bind_cols(node = node_names, in_degree = in_degree, out_degree = out_degree)

saveRDS(dd_stats, file = "data/dd_stats.rds")
