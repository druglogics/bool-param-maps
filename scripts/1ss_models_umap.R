########################
# UMAP for 1 ss models #
########################

library(dplyr)
library(tibble)
library(ggplot2)
library(uwot)

lo_data = readRDS(file = "data/lo_data.rds")
rownames(lo_data) = NULL
lo_data = as.matrix(lo_data)

n_neighbors = 2:20

for(i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  set.seed(42)
  lo_umap = uwot::umap(X = lo_data, n_threads = 8, n_neighbors = i, verbose = TRUE)

  # The result files take a lot of space, so I will just keep a small number of them in the repo
  saveRDS(object = lo_umap, file = paste0("data/1ss_lo_umap/umap_1ss_", i, "nn.rds"))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    as_tibble() %>%
    ggplot(aes(x = X, y = Y)) +
      geom_point(shape = '.') +
      labs(title = paste0("Parameterization Map - ", i, " Neighbours")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

  # save the images
  ggsave(filename = paste0("img/1ss_umap/1ss_umap_unsup_", i, "nn.png"), dpi = "print", width = 7, height = 5)
}
