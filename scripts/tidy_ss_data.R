library(dplyr)
library(emba)

# point to the `models` directory in the Zenodo dataset
data_dir = "/media/disk/abmlog/abmlog_cascade_1.0_models_fixpoints/models"
model_dirs = list.dirs(path = data_dir, recursive = FALSE)

ss_data = list()
index = 1
for (model_dir in model_dirs) {
  print(paste0("Reading directory: ", model_dir, " (", index, ")"))
  ss_data[[index]] = emba::get_stable_state_from_models_dir(models.dir = model_dir)
  index = index + 1
}

res = dplyr::bind_rows(ss_data)

saveRDS(res, file = "1ss_model_data.rds")