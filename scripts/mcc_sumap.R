# Use supervised umap on the 1 stable state model link operator data
# with their respective MCC score as response target data
library(dplyr)
library(tibble)
library(tibble)
library(uwot)
library(ggplot2)

# Read the link operator data (see `get_1ss_lo_data.R`)
lo_data = readRDS(file = "data/lo_data.rds")

# Read the mcc data (see `emba_analysis.R`)
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds")
models_mcc = mcc_res$models.mcc

# check data
stopifnot(all(rownames(lo_data) == names(models_mcc)))

# drop the model names and convert to matrix for faster simulation
models_mcc = unname(models_mcc)
rownames(lo_data) = NULL
lo_data = as.matrix(lo_data)

# Run supervised umap for different values of `n_neighbors` and `target_weight`
n_neighbors = c(2,8,14) # local to global
target_weights = c(0, 0.5, 1) # only data, equal data and target, only target

for (w in target_weights) {
  print(paste0("Target weight: ", w))
  for (i in n_neighbors) {
    print(paste0("#Neighbors: ", i))

    set.seed(42)
    mcc_sumap = uwot::umap(X = lo_data, y = models_mcc,
      n_threads = 4, n_neighbors = i, target_weight = w, verbose = TRUE)
    saveRDS(object = mcc_sumap, file = paste0("data/mcc_sumaps/sumap_", i, "nn_", w, "w.rds"))

    mcc_sumap %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(MCC = models_mcc) %>%
      ggplot(aes(x = X, y = Y, colour = MCC)) +
      geom_point(shape = '.') +
      scale_color_gradient2(low = "red", mid = "grey", high = "green",
        n.breaks = 6, midpoint = mean(models_mcc)) +
      labs(title = paste0("Supervised Parameterization Map (", i, " Neighbours, Target Weight: ", w, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(filename = paste0("img/mcc_sumaps/", i, "nn_", w, "w.png"), dpi = "print", width = 7, height = 5)
  }
}

# sample data and play
# n_samples = 10000
# indexes = sample(x = 1:nrow(lo_data), size = n_samples)
# mcc_sumap = uwot::umap(X = lo_data[indexes,], y = models_mcc[indexes], verbose = TRUE,
#   n_threads = 4, n_neighbors = 14, target_n_neighbors = 2, target_weight = 0.7)
