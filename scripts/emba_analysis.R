# In this analysis the data used includes only the models with 1 stable state
library(dplyr)
library(emba)

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

# Read the model predictions (see `get_model_predictions.R` script)
pred_data = readRDS(file = "data/model_predictions.rds")

# check the model names have the same order
stopifnot(all(rownames(ss_data) == rownames(pred_data)))

# Read the link operator data (see `get_1ss_lo_data.R`)
lo_data = readRDS(file = "data/lo_data.rds")

# check the model names have the same order
stopifnot(all(rownames(pred_data) == rownames(lo_data)))

# observed synergies
obs_syn = emba::get_observed_synergies(file = "data/observed_synergies_cascade_1.0")

# Performance Biomarker Analysis
# split models to a different number of MCC classes (the more classes,
# the more biomarkers will be found). 14 was the total number of unique MCC scores.
for (mcc_classes in 3:14) {
  print(mcc_classes)
  mcc_res = emba::biomarker_mcc_analysis(
    model.predictions = pred_data,
    models.stable.state = ss_data,
    models.link.operator = lo_data,
    observed.synergies = obs_syn,
    threshold = 0.6,
    num.of.mcc.classes = mcc_classes,
    penalty = 0.1)
  # in the repo we are going to keep only a few of them to save space :)
  saveRDS(object = mcc_res, file = paste0("data/mcc_", mcc_classes, "_res.rds"))
}

# Synergy Biomarker analysis
synergy_res = emba::biomarker_synergy_analysis(
  model.predictions = pred_data,
  models.stable.state = ss_data,
  models.link.operator = lo_data,
  observed.synergies = obs_syn,
  threshold = 0.6,
  calculate.subsets.stats = FALSE,
  penalty = 0.1)
saveRDS(object = synergy_res, file = "data/synergy_res.rds")
