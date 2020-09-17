########################
# UMAP for 1 ss models #
########################

library(dplyr)
library(tibble)
library(ggplot2)
library(uwot)

lo_data = readRDS(file = "data/lo_data.rds")

for(n_neighbors in 2:20) {
  print(n_neighbors)
  set.seed(42)

  lo_umap = uwot::umap(X = lo_data, metric = "manhattan", n_threads = 8, n_neighbors = n_neighbors)

  # I have saved them all but they take a lot of space, so I will just keep
  # a small number of them in the repo, namely for n_neighbors = 2,4,6,8,11 and 14
  saveRDS(object = lo_umap, file = paste0("1ss_lo_umap/lo_umap_1ss_", n_neighbors, "nn.rds"))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    as_tibble() %>%
    ggplot(aes(x = X, y = Y)) + geom_point(size = 0.01) + theme_classic()

  # save the images
  ggsave(filename = paste0("img/1ss_umap/1ss_umap_unsup_", n_neighbors, ".png"), dpi = "print", width = 7, height = 5)
}
