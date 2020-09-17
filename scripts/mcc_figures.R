library(emba)
library(Ckmeans.1d.dp)
library(dplyr)
library(tibble)
library(ggplot2)

# MCC Histogram
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds")
models.mcc = mcc_res$models.mcc
num.of.mcc.classes = 4
res = Ckmeans.1d.dp(x = models.mcc, k = num.of.mcc.classes)
models.cluster.ids = res$cluster
png(file = "img/mcc_hist.png", units = "in", width = 7, height = 5, res = 300)
emba::plot_mcc_classes_hist(models.mcc, models.cluster.ids, num.of.mcc.classes)
dev.off()

# MCC Maps

# data check
lo_data = readRDS(file = "data/lo_data.rds")
stopifnot(rownames(lo_data) == names(models.mcc))
# `lo_umap` comes from the `lo_data` (see `1ss_models_umap.R` script)
# so the model order is fine

n_neighbors = c(2,4,6,8,11,14)

for (i in n_neighbors) {
  print(paste0("#Neighbors: ", i))

  lo_umap = readRDS(file = paste0("data/1ss_lo_umap/lo_umap_1ss_", i, "nn.rds"))

  lo_umap %>%
    `colnames<-` (c("X", "Y")) %>%
    tibble::as_tibble() %>%
    tibble::add_column(MCC = models.mcc) %>%
    #slice_sample(n = 10000) %>%
    ggplot(aes(x = X, y = Y, colour = MCC)) +
      geom_point(shape = '.') +
      scale_color_gradient2(low = "red", mid = "grey", high = "green",
        n.breaks = 6, midpoint = mean(models.mcc)) +
      labs(title = paste0("Parameterization Map, MCC colored - ", i, " Neighbours")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("img/mcc_maps/", i, "nn.png"), dpi = "print", width = 7, height = 5)
}