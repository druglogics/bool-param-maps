# Embed the link-operator values of important nodes on the
# parameterization maps generated with unsupervised and supervised
# umap (using all the models of the CASCADE 1.0 dataset)
library(dplyr)
library(tibble)
library(ggplot2)

# get umap coordinates (see `param_ss_umap.R`)
lo_umap = readRDS(file = "data/lo_umap/lo_umap_20nn.rds")
lo_umap_ham = readRDS(file = "data/lo_umap/lo_umap_14nn_ham_0_5_min_dist.rds")
lo_sumap = readRDS(file = "data/lo_umap/lo_sumap_14nn_0_3_min_dist.rds")

# get link-operator data
lo_mat = readRDS(file = "data/lo_mat.rds") # see `get_lo_mat.R`

# The important nodes (2 last are the least important)
imp_nodes = c("MAPK14", "ERK_f", "MEK_f", "PTEN", "mTORC1_c", "CFLAR", "CYCS")

################
# Unsupervised #
################

# hamming metric, 14 neighbors, min_dist = 0.5
for (node in imp_nodes) {
  print(paste0("Node: ", node))

  lo_umap_ham %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_mat[,node] %>% as.factor()) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Unsupervised Hamming Map (14 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/imp_nodes_param_ss_maps/unsup_", node, "_ham.png"), dpi = "print", width = 7, height = 5)
}

# euclidean metric, 20 neighbors
for (node in imp_nodes) {
  print(paste0("Node: ", node))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_mat[,node] %>% as.factor()) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Unsupervised Map (20 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/imp_nodes_param_ss_maps/unsup_", node, ".png"), dpi = "print", width = 7, height = 5)
}

##############
# Supervised #
##############

# euclidean metric, 14 neighbors, min_dist = 0.3
for (node in imp_nodes) {
  print(paste0("Node: ", node))

  lo_sumap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_mat[,node] %>% as.factor()) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Supervised Map (20 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/imp_nodes_param_ss_maps/sup_", node, ".png"), dpi = "print", width = 7, height = 5)
}
