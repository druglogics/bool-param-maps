##############################
# Gitsbe model Fitness Density
##############################

library(dplyr)
library(tibble)
library(emba)
library(usefun)
library(ggplot2)

lo_data = readRDS(file = "data/lo_data.rds")

# Zenodo: 10.5281/zenodo.3988424
ss_gitsbe_models = emba::get_stable_state_from_models_dir(models.dir = "/home/john/tmp/ags_paper_res/cascade_1.0/cascade_1.0_ss_1000sim_fixpoints_hsa_20200609_153545/models")

# read AGS `steady_state` vector (see `get_1ss_fitness_data.R`)
steady_state = readRDS(file = "data/ags_steady_state.rds")

# calculate models fitness to AGS steady state
gitsbe_models_fit = apply(ss_gitsbe_models[, names(steady_state)], 1,
  usefun::get_percentage_of_matches, steady_state)

gitsbe_models_fit %>% as_tibble() %>%
  ggplot() +
  geom_density(aes(x = value), fill = "#4DAF4A", alpha = 0.7, adjust = 2, show.legend = FALSE) +
  xlim(c(0,1)) +
  xlab("Fitness to AGS steady state") +
  ggtitle("Fitness Density of 3000 Gitbse Models") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/gitsbe_fit_density.png", dpi = "print", width = 7, height = 5)
