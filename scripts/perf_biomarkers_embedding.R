# Embed the link-operator values of some important (or not) nodes on the
# parameterization maps generated with umap (supervised and unsupervised)
library(dplyr)
library(tibble)
library(ggplot2)

# which nodes to check? (2 last are the least important)
nodes = c("ERK_f", "MAPK14", "MEK_f", "CTNNB1", "TCF7_f", "CYCS", "CFLAR")

# Read the link operator data (see `get_1ss_lo_data.R`)
lo_data = readRDS(file = "data/lo_data.rds")

# read the link-operator UMAP unsupervised result (14 neighbors)
lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_14nn.rds"))

for (node in nodes) {
  print(paste0("Node: ", node))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_data[,node] %>% as.factor()) %>%
    #slice_sample(n = 10000) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Unsupervised Map (14 Neighbors)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/nodes_lo_maps/unsup_", node, ".png"), dpi = "print", width = 7, height = 5)
}

# read the link-operator UMAP supervised result (14 neighbors) - MCC Classes
lo_sumap = readRDS(file = paste0("data/mcc_sumaps/class/sumap_14nn_0.5w_class.rds"))

for (node in nodes) {
  print(paste0("Node: ", node))

  lo_sumap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_data[,node] %>% as.factor()) %>%
    #slice_sample(n = 10000) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Supervised Map (14 Neighbors, MCC Classes)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/nodes_lo_maps/sup_", node, ".png"), dpi = "print", width = 7, height = 5)
}

# read the link-operator UMAP supervised result (14 neighbors) - MCC as continuous variable
lo_sumap_cont = readRDS(file = paste0("data/mcc_sumaps/sumap_14nn_0.5w.rds"))

for (node in nodes) {
  print(paste0("Node: ", node))

  lo_sumap_cont %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(node_lo = lo_data[,node] %>% as.factor()) %>%
    #slice_sample(n = 10000) %>%
    ggplot(aes(x = X, y = Y, colour = node_lo)) +
    geom_point(shape = '.') +
    scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
    guides(colour = guide_legend(title = paste0(node, " LO"), label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Supervised Map (14 Neighbors, MCC continuous)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/nodes_lo_maps/sup_cont_", node, ".png"), dpi = "print", width = 7, height = 5)
}
