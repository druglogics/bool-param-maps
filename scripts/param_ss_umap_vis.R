###############################################################
# Visualize parameterization of the whole CASCADE 1.0 dataset #
# Data results are from the `param_ss_umap.R` script          #
###############################################################
library(dplyr)
library(ggplot2)
library(scales)

# get number of ss per model
models_ss_stats = readRDS(file = "data/models_ss_stats.rds") # see `count_models_ss.R`
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()
model_number = models_ss_stats %>% pull(model_number)

################
# Unsupervised #
################

n_neighbors = c(8,14,20) # local to global

# visualize the results with `euclidean` distance metric
for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  # get umap coordinates
  lo_umap = readRDS(file = paste0("data/lo_umap/lo_umap_", i, "nn.rds"))

  # tidy data
  data = lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(ss_num = ss_num, model_num = model_number)

  # color by #fixpoints
  data %>%
    ggplot(aes(x = X, y = Y, color = ss_num)) +
    geom_point(shape = '.') +
    guides(colour = guide_legend(title = "#Fixpoints", label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Map (All CASCADE 1.0 models) - ", i, " Neighbors")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("img/all_models_maps/umap_", i, "nn.png"), dpi = "print", width = 7, height = 5)

  # color by model decimal representation
  data %>%
    ggplot(aes(x = X, y = Y, color = model_num)) +
    geom_point(shape = '.') +
    scale_color_distiller(palette = "Spectral", labels = scales::label_number_si(),
      guide = guide_colourbar(title = "Model Number")) +
    labs(title = paste0("Parameterization Map (All CASCADE 1.0 models) - ", i, " Neighbors")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("img/all_models_maps/umap_", i , "nn_model_num.png"), dpi = "print", width = 7, height = 5)
}

# visualize the results with `hamming` distance metric
lo_umap = readRDS(file = "data/lo_umap/lo_umap_14nn_ham_0_5_min_dist.rds")

# tidy data
data = lo_umap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(ss_num = ss_num, model_num = model_number)

# color by #fixpoints
data %>%
  ggplot(aes(x = X, y = Y, color = ss_num)) +
  geom_point(shape = '.') +
  guides(colour = guide_legend(title = "#Fixpoints", label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("Parameterization Map (All CASCADE 1.0 models) - ", i, " Neighbors")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("img/all_models_maps/ham_umap_", i, "nn.png"), dpi = "print", width = 7, height = 5)

# color by model decimal representation
data %>%
  ggplot(aes(x = X, y = Y, color = model_num)) +
  geom_point(shape = '.') +
  scale_color_distiller(palette = "Spectral", labels = scales::label_number_si(),
    guide = guide_colourbar(title = "Model Number")) +
  labs(title = paste0("Parameterization Map (All CASCADE 1.0 models) - ", i, " Neighbors")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("img/all_models_maps/ham_umap_", i , "nn_model_num.png"), dpi = "print", width = 7, height = 5)

##############
# Supervised #
##############

# visualize the results with `euclidean` distance metric
for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  # get umap coordinates
  lo_sumap = readRDS(file = paste0("data/lo_umap/lo_sumap_", i, "nn.rds"))

  # tidy data
  data = lo_sumap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(ss_num = ss_num)

  # color by #fixpoints
  data %>%
    ggplot(aes(x = X, y = Y, color = ss_num)) +
    geom_point(shape = '.') +
    guides(colour = guide_legend(title = "#Fixpoints", label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Parameterization Supervised Map (All CASCADE 1.0 models) - ", i, " Neighbors")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0("img/all_models_maps/sumap_", i, "nn.png"), dpi = "print", width = 7, height = 5)
}

# visualize the results with `euclidean` distance metric, 14 `n_neighbors` and `min_dist` = 0.3
lo_sumap = readRDS(file = "data/lo_umap/lo_sumap_14nn_0_3_min_dist.rds")

# tidy data
data = lo_sumap %>%
  `colnames<-` (c("X", "Y")) %>%
  tibble::as_tibble() %>%
  tibble::add_column(ss_num = ss_num)

# color by #fixpoints
data %>%
  ggplot(aes(x = X, y = Y, color = ss_num)) +
  geom_point(shape = '.') +
  guides(colour = guide_legend(title = "#Fixpoints", label.theme = element_text(size = 12),
    override.aes = list(shape = 19, size = 12))) +
  labs(title = paste0("Parameterization Supervised Map (All CASCADE 1.0 models)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0("img/all_models_maps/sumap_14nn_0_3_min_dist.png"), dpi = "print", width = 7, height = 5)
