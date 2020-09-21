library(emba)
library(dplyr)
library(tibble)
library(Ckmeans.1d.dp)
library(ggplot2)
library(ggpubr)

# read fitness data (see `get_1ss_fitness_data.R`)
fit_data = readRDS(file = "data/fit_data.rds")

# Fitness density figure
fit_data %>% as_tibble() %>%
  ggplot() +
  geom_density(aes(x = value), fill = "#4DAF4A", alpha = 0.7, adjust = 3, show.legend = FALSE) +
  xlim(c(0,1.05)) +
  xlab("Fitness to AGS steady state") +
  ggtitle("Fitness Density  (All link-operator, 1 stable state models)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/1ss_models_fit_density.png", dpi = "print", width = 7, height = 5)

# Fitness Maps

# check: model name order is correct (`lo_umap` comes from `lo_data`)
lo_data = readRDS(file = "data/lo_data.rds")
stopifnot(all(rownames(lo_data) == names(fit_data)))

n_neighbors = c(2,4,6,8,11,14)
for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_", i, "nn.rds"))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(fitness = fit_data) %>%
    #slice_sample(n = 10000) %>%
    ggplot(aes(x = X, y = Y, colour = fitness)) +
    geom_point(shape = '.') +
    scale_color_gradient2(low = "red", mid = "grey", high = "green",
      n.breaks = 6, midpoint = mean(fit_data)) +
    labs(title = paste0("Parameterization Map, Fitness colored - ", i, " Neighbours")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/fit_maps/", i, "nn.png"), dpi = "print", width = 7, height = 5)
}

# Fitness vs MCC correlation
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds")
mcc_data = mcc_res$models.mcc
res = Ckmeans.1d.dp(x = mcc_data, k = 4) # 4 MCC classes

# check model name order
stopifnot(all(names(fit_data) == names(mcc_data)))

# MCC vs fitness boxplot
my_comparisons = list(c(1,2), c(1,3), c(1,4), c(3,4))
bind_cols(fitness = fit_data, MCC = mcc_data, mcc_class = as.factor(res$cluster)) %>%
  #slice_sample(n = 10000) %>%
  ggplot(aes(x = mcc_class, y = fitness, fill = mcc_class)) +
    geom_boxplot() +
    guides(fill = guide_legend(title = "MCC Class")) +
    ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
    labs(title = "MCC vs Fitness (All 1 stable state Models)", x = "MCC Class", y = "Fitness") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/mcc_vs_fit.png", dpi = "print", width = 7, height = 5)

## ggpubr way
# bind_cols(fitness = fit_data, MCC = mcc_data, mcc_class = as.factor(res$cluster)) %>%
#   ggpubr::ggboxplot(x = 'mcc_class', y = 'fitness', fill = 'mcc_class') +
#   ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")

# scatter plot
# bind_cols(fitness = fit_data, MCC = mcc_data) %>%
#   slice_sample(n = 10000) %>%
#   ggplot(aes(x = fitness, y = MCC)) +
#   geom_point() +
#   geom_smooth(method = lm, formula = y ~ x) +
#   stat_cor(method = "kendall", label.y = 0.6) +
#   theme_classic(base_size = 14)
