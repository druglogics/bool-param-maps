###############################################
# For all the models that have 1 stable state,
# get their fitness to the AGS steady state
###############################################

library(dplyr)
library(tibble)
library(stringr)
library(usefun)

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

# get the AGS steady state (literature-curated).
# See S4 Table in Flobak et al. (2015) for more info
steady_state_file = "data/steadystate"
lines = readLines(steady_state_file)
ss_str = unlist(strsplit(x = lines[8], split = "\t"))
ss_mat = stringr::str_split(string = ss_str, pattern = ":", simplify = TRUE)
colnames(ss_mat) = c("nodes", "states")
ss_tbl = ss_mat %>% as_tibble() %>% mutate_at(vars(states), as.integer)

steady_state = ss_tbl %>% pull(states)
names(steady_state) = ss_tbl %>% pull(nodes)

# save `steady_state` vector
saveRDS(object = steady_state, file = "data/ags_steady_state.rds")

# calculate models fitness to AGS steady state
fit_data = apply(ss_data[, names(steady_state)], 1,
  usefun::get_percentage_of_matches, steady_state)

# data check
stopifnot(all(names(fit_data) == rownames(ss_data)))

# save result
saveRDS(object = fit_data, file = "data/fit_data.rds")
