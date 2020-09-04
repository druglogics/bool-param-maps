################
# PCA analysis # - Plotting did not work
################
library(ggbiplot)

lo_mat = readRDS(file = "data/lo_mat.rds")
models_lo_pca = prcomp(x = lo_mat) # `center = TRUE, scale. = TRUE`

# All variables have: sdev = 0.5, 23 PCs are found and I cannot plot them!
# The first variable/node had always rot = 0 and this generated some trouble
# with the plotting function
ggbiplot(models_lo_pca, ellipse = TRUE, circle = TRUE)
# to have colored groups add parameter: `groups = models_ss_stats %>% pull(ss_num) %>% as.factor()`