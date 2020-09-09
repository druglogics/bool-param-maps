##########################
# UMAP visualize results #
##########################

library(dplyr)
library(tibble)
library(usefun)
library(ggplot2)

# get number of ss per model
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
ss_num = models_ss_stats %>% pull(ss_num) %>% as.factor()
model_number = models_ss_stats %>% pull(model_number)

# get important node indexes
lo_mat = readRDS(file = "data/lo_mat.rds")
mapk14_index = which(colnames(lo_mat) == "MAPK14")
erkf_index   = which(colnames(lo_mat) == "ERK_f")
mekf_index   = which(colnames(lo_mat) == "MEK_f")
pten_index   = which(colnames(lo_mat) == "PTEN")
mtorc1_index = which(colnames(lo_mat) == "mTORC1_c")
cflar_index  = which(colnames(lo_mat) == "CFLAR")
n_bits = ncol(lo_mat) # n_bits = 23

# UMAP supervised results
# 3 main clusters corresponding to the 0, 1 and 2 fixpoint models
lo_sumap = readRDS(file = "data/lo_sumap.rds")
colnames(lo_sumap) = c('X','Y')
# 33% of the models used
model_samples_indexes_sup = readRDS(file = "data/indexes_sup.rds")

# get the model numbers in decimal format
model_num = model_number[model_samples_indexes_sup]

print("Get MAPK14 link operators")
mapk14_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, mapk14_index, mapk14_index) %>% as.integer()
})

print("Get ERK_f link operators")
erkf_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, erkf_index, erkf_index) %>% as.integer()
})

print("Get MEK_f link operators")
mekf_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, mekf_index, mekf_index) %>% as.integer()
})

print("Get PTEN link operators")
pten_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, pten_index, pten_index) %>% as.integer()
})

print("Get mTORC1_c link operators")
mtorc1_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, mtorc1_index, mtorc1_index) %>% as.integer()
})

print("Get CFLAR link operators")
cflar_param = sapply(model_num, function(num) {
  bin_n = usefun::dec_to_bin(decimal_num = num, bits = n_bits)
  substr(bin_n, cflar_index, cflar_index) %>% as.integer()
})

# combine all data in a tibble
data = lo_sumap %>%
  as_tibble() %>%
  mutate(ss_num = ss_num[model_samples_indexes_sup]) %>%
  mutate(model_num = model_num) %>%
  add_column(mapk14_param = mapk14_param, erkf_param = erkf_param,
    mekf_param = mekf_param, pten_param = pten_param,
    mtorc1_param = mtorc1_param, cflar_param = cflar_param)

# plot UMAP maps
ggplot(data, mapping = aes(x = X, y = Y)) +
  geom_point(mapping = aes(color = factor(mapk14_param)), size = 0.01) +
  theme_classic() +
  scale_color_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "MAPK14 LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12))) #+
# stat_ellipse(mapping = aes(fill = factor(ss_num)), alpha = 0.15, geom = "polygon", type = "norm") +
#   guides(fill = guide_legend(title = "#Fixpoints", label.theme = element_text(size = 12), override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_MAPK14.png", dpi = "print", width = 7, height = 5)

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = factor(erkf_param)), size = 0.01) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "ERK_f LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_ERK_f.png", dpi = "print", width = 7, height = 5)

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = factor(mekf_param)), size = 0.01) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "MEK_f LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_MEK_f.png", dpi = "print", width = 7, height = 5)

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = factor(pten_param)), size = 0.01) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "PTEN LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_PTEN.png", dpi = "print", width = 7, height = 5)

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = factor(mtorc1_param)), size = 0.01) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "mTORC1_c LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_mTORC1_c.png", dpi = "print", width = 7, height = 5)

ggplot(data) +
  geom_point(aes(x = X, y = Y, color = factor(cflar_param)), size = 0.01) +
  theme_classic() +
  scale_colour_brewer(palette = "Set1", labels = c("AND-NOT", "OR-NOT")) +
  guides(colour = guide_legend(title = "CFLAR LO", label.theme = element_text(size = 12),
    override.aes = list(size = 12)))
ggsave(filename = "img/umap_supervised_CFLAR.png", dpi = "print", width = 7, height = 5)