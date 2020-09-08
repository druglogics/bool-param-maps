##########################
# UMAP visualize results #
##########################

library(dplyr)
library(ggplot2)

# get number of ss per model
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()

# visualize unsupervised dataset
lo_umap = readRDS(file = "data/lo_umap.rds")
model_samples_indexes = readRDS(file = "data/indexes.rds")

colnames(lo_umap) = c('X','Y')
data = lo_umap %>%
  as_tibble() %>%
  mutate(ss_num = ss_num[model_samples_indexes])

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = ss_num), size = 0.01) +
  theme_classic() +
  guides(colour = guide_legend(title = "#Fixpoints", override.aes = list(size = 12)))
ggsave(filename = "img/umap_unsupervised_s001.png", dpi = "print", width = 7, height = 5)

# visualize supervised dataset
lo_sumap = readRDS(file = "data/lo_sumap.rds")
model_samples_indexes_sup = readRDS(file = "data/indexes_sup.rds")

colnames(lo_sumap) = c('X','Y')
data = lo_sumap %>%
  as_tibble() %>%
  mutate(ss_num = ss_num[model_samples_indexes_sup])

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = ss_num), size = 0.01) +
  theme_classic() +
  guides(colour = guide_legend(title = "#Fixpoints", override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_s001.png", dpi = "print", width = 7, height = 5)
