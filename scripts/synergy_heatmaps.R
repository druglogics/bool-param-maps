library(dplyr)
library(tibble)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# Read the model predictions (see `get_model_predictions.R` script)
pred_data = readRDS(file = "data/model_predictions.rds")

# observed synergies
obs_syn = emba::get_observed_synergies(file = "data/observed_synergies_cascade_1.0")

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

stopifnot(all(rownames(ss_data) == rownames(pred_data)))

# sample size of data
sample_size = 10000

# colors for 0 and 1's
my_colors = c("red", "lightyellow")
col_fun = circlize::colorRamp2(breaks = c(0, 1), colors = my_colors)

# Activity State Legend
activity_state_legend = ComplexHeatmap::Legend(nrow = 1, title = "Activity State",
  labels = c("Inhibited", "Active"), legend_gp = gpar(fill = c("red", "lightyellow")),
  direction = "horizontal")

for (drug_combo in obs_syn) {
  print(paste0("Drug combo: ", drug_combo))

  # get the model that predict the `drug_combo` as synergistic
  models = pred_data %>%
    filter(!!as.symbol(drug_combo) == 1) %>%
    rownames_to_column() %>%
    pull(rowname)

  # get the models stable states
  models_ss = ss_data[models,] %>% as.matrix()
  rownames(models_ss) = NULL

  # for reproducibility
  set.seed(42)

  if (drug_combo != "PD-AK") {
    mat = models_ss[sample(x = 1:nrow(models_ss), size = sample_size),]
  } else {
    mat = models_ss # only ~2000 PD-AK synergistic models, so no need for sampling
  }

  # make the heatmap and save it!
  png(file =  paste0("img/synergy_heatmaps/", drug_combo, "_ss_heat.png"),
      width = 7, height = 5, units = "in", res = 600)
  ss_heat = ComplexHeatmap::Heatmap(matrix = mat,
    col = col_fun, column_names_gp = gpar(fontsize = 6),
    use_raster = TRUE, raster_device = "png",
    column_title = paste0(drug_combo, " Synergistic Models Activity State Patterns"),
    show_heatmap_legend = FALSE, heatmap_legend_param = list(direction = "horizontal"))
  ComplexHeatmap::draw(ss_heat, annotation_legend_list = activity_state_legend,
                       annotation_legend_side = "bottom")
  dev.off()
}
