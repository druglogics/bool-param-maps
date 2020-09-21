#############################################
# Embed the gitsbe-generated best-fit models
# in the UMAP 2D parameterization map
#############################################

library(dplyr)
library(tibble)
library(emba)
library(ggplot2)

lo_data = readRDS(file = "data/lo_data.rds")

# tar -xzvf data/cascade_1.0_ss_1000sim_fixpoints_hsa.tar.gz
lo_gitsbe_models = emba::get_link_operators_from_models_dir(models.dir = "data/cascade_1.0_ss_1000sim_fixpoints_hsa_20200921_101955/models")

df_args = c(lo_gitsbe_models, sep = "")
gitsbe_bin_models = do.call(paste, df_args)
gitsbe_num_models = paste0("network_", base::strtoi(gitsbe_bin_models, base = 2))

# from the total of 3000 models, 2566 have a unique link-operator parameterization
gitsbe_num_models = unique(gitsbe_num_models)

# data check: gitsbe model numerical representation is part of the whole set!
all(gitsbe_num_models %in% rownames(lo_data))

indexes = which(rownames(lo_data) %in% gitsbe_num_models)

is_gitsbe_model = sapply(rownames(lo_data), function(model_name) {
  ifelse(model_name %in% gitsbe_num_models, yes = 1, no = 0)
}, USE.NAMES = FALSE) %>% as.factor()

# data check
all(which(is_gitsbe_model == 1) == indexes)

# Plot
# set.seed(42)
# rnd_indexes = sample(x = 1:nrow(lo_umap), size = 100000)
# my_sizes  = ifelse(is_gitsbe_model[rnd_indexes] == 1, 2, 0.01)
# my_shapes = ifelse(is_gitsbe_model[rnd_indexes] == 1, 'x', '.')
# my_alphas = ifelse(is_gitsbe_model[rnd_indexes] == 1, 1, 0.1)

my_sizes  = ifelse(is_gitsbe_model == 1, 2, 0.01)  # larger size for crosses to show them up
my_shapes = ifelse(is_gitsbe_model == 1, 'x', '.') # 'x' is a cross, '.' is a pixel
my_alphas = ifelse(is_gitsbe_model == 1, 1, 0.01)

n_neighbors = c(2,4,6,8,11,14)

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_", i, "nn.rds"))
  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(is_gitsbe = is_gitsbe_model) %>%
    #slice(rnd_indexes) %>%
    ggplot(aes(x = X, y = Y, colour = is_gitsbe)) +
    geom_point(size = my_sizes, shape = my_shapes, alpha = my_alphas) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = paste0("Parameterization Map - ", i, " Neighbours")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/gitsbe_umaps/", i, "nn.png"), dpi = "print", width = 7, height = 5)
}
