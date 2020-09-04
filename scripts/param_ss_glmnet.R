############################
# Multinomial glmnet LASSO #
############################

library(glmnet)
library(dplyr)

# read data
models_ss_stats = readRDS(file = "models_ss_stats.rds")
lo_mat = readRDS(file = "lo_mat.rds")
ss_num = models_ss_stats %>% pull(ss_num)

# fit glmnet multinomial lasso (a = 1) - all model dataset
fit_a1 = glmnet::glmnet(x = lo_mat, y = ss_num, family = "multinomial", alpha = 1, type.multinomial = "grouped", standardize = FALSE, intercept = FALSE, trace.it = 1)

# save result
saveRDS(object = fit_a1, file = "fit_a1.rds")

# fit glmnet multinomial lasso (a = 1) - perform Cross-validation on simple random samples of the dataset
set.seed(42)
cvfit_data = list()
for(i in 1:20) { # sample 20 x 100000 models and do CV
  print(i)
  model_samples_indexes = sample(x = 1:nrow(lo_mat), size = 100000)
  cvfit_data[[i]] = glmnet::cv.glmnet(x = lo_mat[model_samples_indexes, ], y = ss_num[model_samples_indexes], family = "multinomial", alpha = 1, type.multinomial = "grouped", standardize = FALSE, intercept = FALSE, trace.it = 1)
}

# save result
saveRDS(object = cvfit_data, file = "cvfit_data.rds")
