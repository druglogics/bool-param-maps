library(stringr)
library(dplyr)

# point to the `models` directory in the Zenodo dataset
data_dir = "/media/disk/abmlog/abmlog_cascade_1.0_models_fixpoints/models"
model_dirs = list.dirs(path = data_dir, recursive = FALSE)
model_index = 1

data_list = list()
data_index = 1
for (model_dir in model_dirs) {
  print(paste0("Reading directory: ", model_dir, " (", model_index, ")"))

  # filter only the gitsbe files
  gitsbe_files = stringr::str_subset(string =
      list.files(path = model_dir, full.names = TRUE), pattern = ".gitsbe")

  for (file in gitsbe_files) {
    # read first 6 lines, no more are needed
    lines = readLines(con = file, n = 6)
    model_number = stringr::str_remove(string =
        lines[stringr::str_which(lines, "network_")], pattern = "modelname: network_")
    ss_num = sum(stringr::str_count(string = lines, pattern = "stablestate"))
    data_list[[data_index]] = dplyr::bind_cols(model_number = as.integer(model_number), ss_num = ss_num)
    data_index = data_index + 1
  }

  model_index = model_index + 1
}

models_ss_stats = dplyr::bind_rows(data_list)
saveRDS(models_ss_stats, file = "data/models_ss_stats.rds")
