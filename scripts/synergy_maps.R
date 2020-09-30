library(dplyr)
library(tibble)
library(emba)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)

# Read the link operator data (see `get_1ss_lo_data.R`)
lo_data = readRDS(file = "data/lo_data.rds")

# Read the model predictions (see `get_model_predictions.R` script)
pred_data = readRDS(file = "data/model_predictions.rds")

# order check
stopifnot(rownames(lo_data) == rownames(pred_data))

# observed synergies
obs_syn = emba::get_observed_synergies(file = "data/observed_synergies_cascade_1.0")

# Find the number of models that predicted all poassible observed synergy subsets
synergy_subset_stats = emba::get_synergy_subset_stats(
  model.predictions = pred_data %>% select(all_of(obs_syn)),
  synergies = obs_syn)
# remove the subsets where no model predicted them
synergy_subset_stats = synergy_subset_stats[synergy_subset_stats > 0]
saveRDS(object = synergy_subset_stats, file = "data/synergy_subset_stats.rds")

# Tidy models prediction synergy data
pred_data = pred_data %>%
  as_tibble() %>% # removes model names
  select(all_of(obs_syn)) %>% # keep only the observed synergies
  mutate(across(everything(), ~tidyr::replace_na(.x, -1))) # NA values as -1

#################
# Synergy Stats #
#################
stat_data = pred_data %>%
  purrr::map(function(col) { table(col) %>% as.integer() }) %>%
  bind_cols() %>%
  add_column(is_synergy_tally = c("NA", "no", "yes"))

stat_data %>%
  tidyr::pivot_longer(cols = all_of(obs_syn)) %>%
  ggplot(aes(x = name, y = value, fill = is_synergy_tally)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    geom_text(aes(label = paste0(as.character(round(value/1e3)), "K")),
      vjust = -0.5, position = position_dodge(0.9), size = 3.5) +
    scale_fill_manual(values = c("black", "red2", "green"), labels = c("NA", "NO", "YES")) +
    scale_y_continuous(labels = scales::label_number_si(accuracy = 0.1)) +
    guides(fill = guide_legend(title = "Synergy", label.theme = element_text(size = 12))) +
    labs(x = "Observed Synergies", y = "Number of Models") +
    theme_classic(base_size = 14)
ggsave(filename = "img/synergy_stats.png", dpi = "print", width = 7, height = 5)

################
# Synergy Maps #
################

# set number of neighbors
n_neighbors = c(2,8,11,14)

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_", i, "nn.rds"))

  for (drug_combo in obs_syn) {
    print(paste0("Drug Combination: ", drug_combo))

    is_synergy = pred_data %>% pull(drug_combo)
    #my_sizes  = ifelse(is_synergy == 1, 1.5, 0.01)  # larger size for crosses to show up
    #my_shapes = ifelse(is_synergy == 1, 'x', '.') # 'x' is a cross, '.' is a pixel

    #indx = sample(x = 1:nrow(lo_umap), size = 50000)

    lo_umap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(is_synergy = as.factor(is_synergy)) %>%
      #slice(indx) %>%
      ggplot(aes(x = X, y = Y, colour = is_synergy)) +
      geom_point(shape = '.') +
      #geom_point(shape = my_shapes[indx], size = my_sizes[indx]) +
      scale_colour_manual(values = c("black", "red2", "green"), labels = c("NA", "NO", "YES")) +
      guides(colour = guide_legend(title = "Synergy", label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12))) +
      labs(title = paste0(drug_combo, " Synergy Parameterization Map - ", i, " Neighbours")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(filename = paste0("img/synergy_maps/", drug_combo, "_", i, "nn.png"), dpi = "print", width = 7, height = 5)
  }
}

############################
# All synergies in one map #
############################

# from `synergy_subset_stats` we know that the only subsets predicted are:
# `PI-PD`, `PI-5Z`, `PD-AK`, `AK-5Z`, `PI-PD,PD-AK` and `PI-5Z,AK-5Z`
# Also, the models predicting `PD-AK` form a proper subset of the `PI-PD,PD-AK`
# models and the same for `PI-5Z` and `PI-5Z,AK-5Z`
which_syn_set = pred_data %>%
  mutate(syn_set = case_when(
    `PI-5Z` == 1 ~ "PI-5Z,AK-5Z",
    `AK-5Z` == 1 ~ "AK-5Z", # this does not include the models that have also `PI-5Z` == 1
    `PD-AK` == 1 ~ "PI-PD,PD-AK",
    `PI-PD` == 1 ~ "PI-PD", # this does not include the models that have also `PD-AK` == 1
    TRUE ~ "none")) %>%
  pull(syn_set) %>%
  forcats::as_factor()
which_syn_set = factor(which_syn_set, levels(which_syn_set)[c(1,3,2,4,5)])

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_", i, "nn.rds"))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(which_syn_set = which_syn_set) %>%
    ggplot(aes(x = X, y = Y, colour = which_syn_set)) +
    geom_point(shape = '.') +
    scale_colour_manual(values = c("black", "magenta", "green", "yellow", "steelblue")) +
    guides(colour = guide_legend(title = "Synergy Subsets", label.theme = element_text(size = 12),
      override.aes = list(shape = 19, size = 12))) +
    labs(title = paste0("Synergy Parameterization Map - ", i, " Neighbours")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/synergy_maps/all_syn_", i, "nn.png"), dpi = "print", width = 7, height = 5)
}
