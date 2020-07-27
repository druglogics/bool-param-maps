library(gtools)
library(tibble)
library(dplyr)

bool_values = c(TRUE, FALSE)

# lv = c(T,T,F)
# apply_operator(lv, 'or')
# apply_operator(lv, 'and')
apply_operator = function(logical_vector, operator = c('and','or')) {
  stopifnot(operator %in% c('and','or'))
  res = logical_vector[1]

  if (length(logical_vector) == 1) return(res)

  for (index in seq.int(from = 2, to = length(logical_vector))) {
    if (operator == 'and')
      res = res & logical_vector[index]
    else
      res = res | logical_vector[index]
  }
  return(res)
}

# lv => logical vector (e.g. lv = c(T,F,F,T))
# num_act => number of activators
# num_inh => number of inhibitors

# (a OR b) AND NOT (c OR d)
and_not = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') &! apply_operator(lv[(num_act+1):(num_act+num_inh)], 'or')
}
# (a OR b) OR NOT (c OR d)
or_not = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') |! apply_operator(lv[(num_act+1):(num_act+num_inh)], 'or')
}
# (a OR b) AND (NOT c OR NOT d)
balance1 = function(lv, num_act, num_inh) {
  apply_operator(lv[1:num_act], 'or') & apply_operator(!lv[(num_act+1):(num_act+num_inh)], 'or')
}

# activators win on equality (more activators always win)
exp_act_win = function(lv, num_act, num_inh) {
  act_sum = sum(lv[1:num_act])
  if (act_sum == 0) return(FALSE)
  inh_sum = sum(lv[(num_act+1):(num_act+num_inh)])
  if (act_sum >= inh_sum) return(TRUE) else return(FALSE)
}

# inhibitors win on equality (more activators always win)
exp_inh_win = function(lv, num_act, num_inh) {
  act_sum = sum(lv[1:num_act])
  if (act_sum == 0) return(FALSE)
  inh_sum = sum(lv[(num_act+1):(num_act+num_inh)])
  if (act_sum > inh_sum) return(TRUE) else return(FALSE)
}


####################################
# 1 act + 1 inh
d3 = permutations(v = bool_values, n = 2, r = 2, repeats.allowed = TRUE)
and_not_res  = as.integer(apply(d3, 1, and_not, 1, 1))
or_not_res   = as.integer(apply(d3, 1, or_not, 1, 1))
balance1_res = as.integer(apply(d3, 1, balance1, 1, 1))
exp_act_res  = as.integer(apply(d3, 1, exp_act_win, 1, 1))
exp_inh_res  = as.integer(apply(d3, 1, exp_inh_win, 1, 1))

d3 = d3 %>%
  as_tibble() %>%
  mutate_all(as.integer) %>%
  add_column(and_not_res = and_not_res) %>%
  add_column(or_not_res = or_not_res) %>%
  add_column(balance1_res = balance1_res) %>%
  add_column(exp_act_res = exp_act_res) %>%
  add_column(exp_inh_res = exp_inh_res)

###################################
# 2 act + 1 inh
d3 = permutations(v = bool_values, n = 2, r = 3, repeats.allowed = TRUE)
and_not_res  = as.integer(apply(d3, 1, and_not, 2, 1))
or_not_res   = as.integer(apply(d3, 1, or_not, 2, 1))
balance1_res = as.integer(apply(d3, 1, balance1, 2, 1))
exp_act_res  = as.integer(apply(d3, 1, exp_act_win, 2, 1))
exp_inh_res  = as.integer(apply(d3, 1, exp_inh_win, 2, 1))

d3 = d3 %>%
  as_tibble() %>%
  mutate_all(as.integer) %>%
  add_column(and_not_res = and_not_res) %>%
  add_column(or_not_res = or_not_res) %>%
  add_column(balance1_res = balance1_res) %>%
  add_column(exp_act_res = exp_act_res) %>%
  add_column(exp_inh_res = exp_inh_res)

##################################
# 1 act + 2 inh
d3 = permutations(v = bool_values, n = 2, r = 3, repeats.allowed = TRUE)
and_not_res  = as.integer(apply(d3, 1, and_not, 1, 2))
or_not_res   = as.integer(apply(d3, 1, or_not, 1, 2))
balance1_res = as.integer(apply(d3, 1, balance1, 1, 2))
exp_act_res  = as.integer(apply(d3, 1, exp_act_win, 1, 2))
exp_inh_res  = as.integer(apply(d3, 1, exp_inh_win, 1, 2))

d3 = d3 %>%
  as_tibble() %>%
  mutate_all(as.integer) %>%
  add_column(and_not_res = and_not_res) %>%
  add_column(or_not_res = or_not_res) %>%
  add_column(balance1_res = balance1_res) %>%
  add_column(exp_act_res = exp_act_res) %>%
  add_column(exp_inh_res = exp_inh_res)

####################################
# 2 act + 2 inh
d3 = permutations(v = bool_values, n = 2, r = 4, repeats.allowed = TRUE)
and_not_res  = as.integer(apply(d3, 1, and_not, 2, 2))
or_not_res   = as.integer(apply(d3, 1, or_not, 2, 2))
balance1_res = as.integer(apply(d3, 1, balance1, 2, 2))
exp_act_res  = as.integer(apply(d3, 1, exp_act_win, 2, 2))
exp_inh_res  = as.integer(apply(d3, 1, exp_inh_win, 2, 2))

d3 = d3 %>%
  as_tibble() %>%
  mutate_all(as.integer) %>%
  add_column(and_not_res = and_not_res) %>%
  add_column(or_not_res = or_not_res) %>%
  add_column(balance1_res = balance1_res) %>%
  add_column(exp_act_res = exp_act_res) %>%
  add_column(exp_inh_res = exp_inh_res)


truth_density_andnot = 100 * sum(and_not_res)/length(and_not_res)
truth_density_ornot = 100 * sum(or_not_res)/length(or_not_res)
truth_density_balance1 = 100 * sum(balance1_res)/length(balance1_res)

############################################
# manual method (you have to write down the variables: the names of the regulators)
d3 = d3 %>%
  mutate(and_not = (a | b) &! (c)) %>%
  mutate(or_not = (a | b) |! (c)) %>%
  mutate(balance1 = (a | b) & (!c))
d3 = d3 %>%
  mutate(and_not = (a) &! (b | c)) %>%
  mutate(or_not = (a) |! (b | c)) %>%
  mutate(balance1 = (a) & (!b | !c))

d3 = d3 %>% mutate_all(as.numeric)
active_res = d3 %>%
  summarise_at(.vars = c('and_not', 'or_not', 'balance1'),
    .funs = c(function(x) {sum(x)/n()}))


d = permutations(v = bool_values, n = 2, r = 4, repeats.allowed = TRUE)
colnames(d) = c('a','b','c','d')
d = as_tibble(d)

d = d %>%
  mutate(and_not = (a | b) &! (c | d)) %>%
  mutate(or_not = (a | b) |! (c | d)) %>%
  mutate(new1 = (a | b) & (!c | !d))
active_res_2 = d %>%
  summarise_at(.vars = c('and_not', 'or_not', 'new1'),
    .funs = c(function(x) {sum(x)/n()}))

# TRUE = 1, FALSE = 0
d = d %>% mutate_all(as.numeric)





