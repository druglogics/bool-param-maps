---
title: "A study in boolean model parameterization"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 11 December, 2020"
description: "Investigations related to link operators mutations in boolean models"
url: 'https\://druglogics.github.io/bool-param-maps/'
github-repo: "druglogics/bool-param-maps"
bibliography: references.bib
link-citations: true
site: bookdown::bookdown_site
---

# Intro {-}

The purpose of this analysis is to make an investigation regarding the relation of **boolean model parameterization** (based on the [standardized equation form](https://druglogics.github.io/druglogics-doc/gitsbe-description.html#default-equation) inspired by [@Mendoza2006]) and various other model attributes, e.g. their drug combination prediction performance, their ability to predict synergies, their fitness to a specific cancerous cell line activity profile, as well as more general directives, such as the identification of **essential nodes that drive the change of dynamics** throughout the parameterization space.

The dataset used in this analysis is derived from the *CASCADE 1.0* topology [@cascade2020].
The cancer signaling network has a total of $77$ nodes, out of which $23$ have link operators in their respective boolean equation (both activating and inhibiting regulators).
Using the [abmlog](https://github.com/druglogics/abmlog) software, we generated all $2^{23} = 8388608$ possible link operator mutated models for the CASCADE 1.0 topology.

We use the relative new, *non-linear dimension reduction* method UMAP [@McInnes2018a] to place all the generated boolean models in a **boolean parameterization map** and most of our conclusions are based on observing patterns in the produced maps.

:::{.green-box}
The next picture gives a general overview:

```r
knitr::include_graphics(path = "img/abmlog_to_maps.png")
```

<img src="img/abmlog_to_maps.png" width="3121" />
:::

:::{.note}
The model dataset is available in Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4022783.svg)](https://doi.org/10.5281/zenodo.4022783).
:::

Libraries used in this analysis:

```r
library(xfun)
library(knitr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(stringr)
library(ggplot2)
library(DT)
library(usefun)
library(emba)
library(forcats)
library(scales)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(glmnet)
library(randomForest)
library(ranger)
library(uwot)
library(Ckmeans.1d.dp)
```

# CASCADE 1.0 - General {-}

## Network Properties {-}

:::{.blue-box}
In this section we demonstrate the **scale-free properties of the CASCADE 1.0 network**.
We show that both in- and out-degree distributions are asymptotically power-law.
:::

Use the script [get_distribution_stats.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/get_distribution_stats.R) to generate the degree distribution stats.
We load the results:


```r
dd_stats = readRDS(file = "data/dd_stats.rds")
```


```r
dd_stats %>% group_by(in_degree) %>% tally() %>%
  ggplot(aes(x = in_degree, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) + 
  geom_smooth(aes(color = "red"), se = FALSE, show.legend = FALSE) + 
  theme_classic() +
  labs(title = "In-Degree Distribution (CASCADE 1.0)", x = "In Degree", y = "Number of Nodes")

dd_stats %>% group_by(out_degree) %>% tally() %>%
  ggplot(aes(x = out_degree, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  geom_smooth(aes(color = "red"), se = FALSE, span = 0.58, show.legend = FALSE) + 
  theme_classic() +
  labs(title = "Out-Degree Distribution (CASCADE 1.0)", x = "Out Degree", y = "Number of Nodes")
```

<div class="figure">
<img src="index_files/figure-html/in-degree-fig-1.png" alt="Degree Distribution (CASCADE 1.0)" width="50%" /><img src="index_files/figure-html/in-degree-fig-2.png" alt="Degree Distribution (CASCADE 1.0)" width="50%" />
<p class="caption">(\#fig:in-degree-fig)Degree Distribution (CASCADE 1.0)</p>
</div>

## Model Stable State Statistics {-}

The `gitsbe` files of the model dataset include also the fixpoint attractors of each model (`.bnet` files have only the equations).
Thus we can find the *frequency distribution* of the number of fixpoints across all produced models (use the script [count_models_ss.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/count_models_ss.R)).
The model stable state (fixpoint) statistics are as follows:


```r
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")

models_ss_stats %>% group_by(ss_num) %>% tally() %>%
  ggplot(aes(x = ss_num, y = n, fill = as.factor(ss_num))) +
  geom_bar(stat = "identity", show.legend = FALSE) + 
  scale_y_continuous(labels = scales::label_number_si()) +
  geom_text(aes(label = n), vjust = -0.5) +
  geom_text(aes(label = paste0(100 * round(n/nrow(models_ss_stats), digits = 2), "%")), size = 10, vjust = c(2.5, 2.5, -2)) +
  theme_classic2() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Stable States Distribution", x = "Number of Stable States", y = "Number of models")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/ss-stats-all-1.png" alt="Stable States Distribution across all link-operator parameterized models (CASCADE 1.0)" width="2100" />
<p class="caption">(\#fig:ss-stats-all)Stable States Distribution across all link-operator parameterized models (CASCADE 1.0)</p>
</div>

:::{.green-box}
Less than $50\%$ of the total possible parameterized models have a single fixpoint attractor which corresponds to a single stable phenotype behavior.
:::

# CASCADE 1.0 - Parameterization vs #fixpoints {-}

:::{.blue-box}
In this section we identify the **key nodes** whose parameterization affects the *change of dynamics* of the CASCADE 1.0 network, i.e. are responsible for the **change in the number of fixpoint attractors (0,1 and 2)** across all link-operator mutated models.
:::

We will use several statistical methods, in each of the sub-sections below.

:::{.orange-box}
The training data is a **link-operator matrix**, where rows are models ($2^{23}$), columns/features/variables are link-operator nodes ($23$ in total) and the parameterization values correspond to $0$ (`AND-NOT`) or $1$ (`OR-NOT`). 
The ternary response for each model is a number denoting the number of fixpoints ($0,1$ or $2$).
:::

The matrix we can generate with the script [get_lo_mat.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/get_lo_mat.R) and the response is part of the previously generated data from the script [count_models_ss.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/count_models_ss.R).

## Multinomial LASSO {-}

Use the script [param_ss_glmnet.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_glmnet.R) to fit a **multinomial LASSO model** for the data [@Friedman2010].
We now simply load the result object:

```r
fit_a1 = readRDS(file = "data/fit_a1.rds")
plot(fit_a1, xvar = "dev", type.coef = "2norm")
plot(fit_a1, xvar = "lambda", type.coef = "2norm")
```

<div class="figure">
<img src="index_files/figure-html/param-ss-glmnet-1.png" alt="Euclidean Norm of glment coefficients vs lambda and deviance explained" width="50%" /><img src="index_files/figure-html/param-ss-glmnet-2.png" alt="Euclidean Norm of glment coefficients vs lambda and deviance explained" width="50%" />
<p class="caption">(\#fig:param-ss-glmnet)Euclidean Norm of glment coefficients vs lambda and deviance explained</p>
</div>

As we can see there is no $\lambda$ that could explain more than $44\%$ of the deviance and there are a lot of non-zero coefficients associated with smaller values of $\lambda$.
For example, choosing  $\lambda = 0.0142$ (tested prediction accuracy $\approx 0.72$ on a random subset of the data), we have the following coefficients, shown in a heatmap:

```r
# choose a lambda
lambda1 = fit_a1$lambda[32]

fit_a1_coef = coef(fit_a1, s = lambda1) %>%
  lapply(as.matrix) %>%
  Reduce(cbind, x = .) %>% t()
rownames(fit_a1_coef) = names(coef(fit_a1, s = lambda1)) # 0, 1 and 2 stable states

imp_nodes_colors = rep("black", length(colnames(fit_a1_coef)))
names(imp_nodes_colors) = colnames(fit_a1_coef)
imp_nodes_colors[names(imp_nodes_colors) %in% c("MEK_f", "PTEN", "MAPK14")] = "green4"

ComplexHeatmap::Heatmap(matrix = fit_a1_coef, name = "Coef", row_title = "Number of Fixpoints", 
  row_names_side = "right", row_dend_side = "right", row_title_side = "right",
  column_title = "Glmnet Coefficient Scores (λ = 0.0142)",
  column_dend_height = unit(20, "mm"), column_names_gp = gpar(col = imp_nodes_colors))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-glmnet-2-1.png" alt="Heatmap of coefficients of multinomial model (glmnet)" width="2100" />
<p class="caption">(\#fig:param-ss-glmnet-2)Heatmap of coefficients of multinomial model (glmnet)</p>
</div>

```r
# check accuracy
# lo_mat = readRDS(file = "data/lo_mat.rds")
# models_ss_stats = readRDS(file = "data/models_ss_stats.rds")
# ss_num = models_ss_stats %>% pull(ss_num)
# set.seed(42)
# model_indexes = sample(x = 1:nrow(lo_mat), size = 100000)
# pred = predict(object = fit_a1, newx = lo_mat[model_indexes, ],
#   type = "class", s = lambda1) %>% as.integer()
# acc = sum(pred == ss_num[model_indexes])/length(pred)
# print(acc) # should be ~ 0.72
```

:::{.green-box}
Even thought the *glmnet* classifier might not be accurate enough, we still find that the nodes `PTEN` and `MAPK14` are the **most important (larger coefficients)** for distinguishing between the models with different number of fixpoints.
Additional nodes (like `MEK_f` and `CTNNB1`) are likely to be important as well.
:::

If we cross-validate the regularization parameter $\lambda$ (using the same script, we randomly selected smaller model samples - $100000$, $20$ times in total), and choose the $\lambda_{1se}$ for each different run to get the coefficients, the results can be visualized as follows:

```r
cvfit_data = readRDS(file = "data/cvfit_data.rds")
cvfit_mat_list = lapply(cvfit_data, function(cvfit) {
  co = coef(cvfit) %>% 
    lapply(as.matrix) %>%
    Reduce(cbind, x = .) %>% 
    t()
  rownames(co) = names(coef(cvfit)) # 0,1 and 2 stable states
  return(co)
})

cvfit_mat = do.call(rbind, cvfit_mat_list)
```


```r
imp_nodes_colors = rep("black", length(colnames(cvfit_mat)))
names(imp_nodes_colors) = colnames(cvfit_mat)
imp_nodes_colors[names(imp_nodes_colors) %in% c("MEK_f", "PTEN", "MAPK14", "CTNNB1", "mTORC1_c")] = "green4"

ComplexHeatmap::Heatmap(matrix = cvfit_mat, name = "Coef", row_title = "Number of Fixpoints", 
  row_dend_side = "right", row_title_side = "left",
  column_title = "Glmnet Coefficient Scores", row_km = 3, row_km_repeats = 10, 
  show_row_names = FALSE, column_dend_height = unit(20, "mm"),
  column_names_gp = gpar(col = imp_nodes_colors),
  left_annotation = rowAnnotation(foo = anno_block(
    labels = c("2", "1", "0"), # with `show_row_names = TRUE` you can check this
    labels_gp = gpar(col = "black", fontsize = 12))))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-glmnet-cv-fig-1.png" alt="Heatmap of coefficients of multinomial model (glmnet - 20 CV models)" width="2100" />
<p class="caption">(\#fig:param-ss-glmnet-cv-fig)Heatmap of coefficients of multinomial model (glmnet - 20 CV models)</p>
</div>

:::{.green-box}
The top 5 most important nodes are seen in green in the above heatmap: `MAPK14`, `PTEN`, `CTNNB1`, `MEK_f` and `mTORC1_c`.
:::

## Random Forests {-}

We used the [param_ss_randf.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_randf.R) script to tune and train a random forest classifier on the dataset [@Liaw2002].
First we tuned the `mtry` parameter, the number of variables randomly selected at each tree split:

```r
mtry_data = readRDS(file = "data/mtry_data.rds")

mtry_data %>%
  ggplot(aes(x = mtry, y = OOBError, group = mtry)) +
  geom_boxplot(fill = "green4") +
  labs(title = "Tuning Random Forest mtry parameter") +
  theme_classic(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-rf-mtry-1.png" alt="Random Forest Tuning (mtry)" width="2100" />
<p class="caption">(\#fig:param-ss-rf-mtry)Random Forest Tuning (mtry)</p>
</div>

:::{.green-box}
A value between $14-18$ for `mtry` seems to minimize the Out-Of-Bag Error, so we choose $16$ for the rest of the analysis.
For the number of trees parameter, we stayed with the default value ($500$) as we observed that they were more than enough after a few test runs.
:::

Next, we randomly selected $100000$ models from the dataset - a total of $20$ times - to train the random forest classifier and find the **importance score of each node**.
We load the result data and tidy up a bit:

```r
rf_imp_data = readRDS(file = "data/rf_imp_data.rds")

# make a list of tibbles
tbl_list = lapply(rf_imp_data, function(mat) {
  nodes = rownames(mat)
  tibble::as_tibble(mat) %>% mutate(nodes = nodes)
})

# OneForAll
imp_res = dplyr::bind_rows(tbl_list)

# Get the importance stats
imp_stats = imp_res %>% 
  group_by(nodes) %>% 
  summarise(mean_acc  = mean(MeanDecreaseAccuracy), 
            sd_acc    = sd(MeanDecreaseAccuracy), 
            mean_gini = mean(MeanDecreaseGini), 
            sd_gini   = sd(MeanDecreaseGini), 
            .groups = 'keep') %>% 
  ungroup()
```

The importance scores for each node were the **mean decrease in accuracy and node impurity** (Gini Index).
We calculate the mean importance and standard deviation scores across all random samples for both measures of importance:

```r
# color first 5 nodes in x axis
imp_col = c(rep("green4", 5), rep("grey30", nrow(imp_stats)-5))
imp_stats %>%
  mutate(nodes = forcats::fct_reorder(nodes, desc(mean_acc))) %>%
  ggplot(aes(x = nodes, y = mean_acc, fill = mean_acc)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_gradient(low = "steelblue", high = "red") +
  geom_errorbar(aes(ymin=mean_acc-sd_acc, ymax=mean_acc+sd_acc), width = 0.2) +
  theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, colour = imp_col)) +
  labs(title = "Random Forest Variable Importance (Accuracy)", 
    x = "Nodes", y = "Mean Decrease Accuracy")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-rf-imp-fig-1.png" alt="Random Forest: Mean Decrease Accuracy per node" width="2100" />
<p class="caption">(\#fig:param-ss-rf-imp-fig)Random Forest: Mean Decrease Accuracy per node</p>
</div>


```r
imp_stats %>%
  mutate(nodes = forcats::fct_reorder(nodes, desc(mean_gini))) %>%
  ggplot(aes(x = nodes, y = mean_gini, fill = mean_gini)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_gradient(low = "steelblue", high = "red") +
  geom_errorbar(aes(ymin = mean_gini-sd_gini, ymax = mean_gini+sd_gini), width = 0.2) +
  theme_classic(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90, colour = imp_col)) +
  labs(title = "Random Forest Variable Importance (Gini Index)", 
    x = "Nodes", y = "Mean Decrease in Node Impurity")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-rf-imp-fig2-1.png" alt="Random Forest: Mean Decrease in Node Impurity (Gini Index)" width="2100" />
<p class="caption">(\#fig:param-ss-rf-imp-fig2)Random Forest: Mean Decrease in Node Impurity (Gini Index)</p>
</div>

:::{.green-box}
The top 5 important nodes by any of the two importance measures using random forests, include the nodes found as important with the LASSO method: `MAPK14`, `ERK_f`, `MEK_f`, `PTEN`, `mTORC1_c`.
:::

Same results we get when using a faster, more memory efficient and with multi-thread support, random forest R package, namely `ranger` [@Wright2017].
Use the script [param_ss_ranger.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_ranger.R) to reproduce the results:

```r
ranger_res = readRDS(file = "data/ranger_res.rds")
imp_res = tibble(nodes = names(ranger_res$variable.importance), 
  gini_index = ranger_res$variable.importance)

imp_res %>% 
  mutate(nodes = forcats::fct_reorder(nodes, desc(gini_index))) %>%
  ggplot(aes(x = nodes, y = gini_index, fill = gini_index)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_gradient(low = "steelblue", high = "red") +
    scale_y_continuous(labels = scales::label_number_si()) +
    theme_classic(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 90, colour = imp_col)) +
    labs(title = "Random Forest (Ranger) Variable Importance (Gini Index)", 
      x = "Nodes", y = "Mean Decrease in Node Impurity")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/param-ss-rf-imp-fig3-1.png" alt="Random Forest (ranger): Mean Decrease in Node Impurity (Gini Index)" width="2100" />
<p class="caption">(\#fig:param-ss-rf-imp-fig3)Random Forest (ranger): Mean Decrease in Node Impurity (Gini Index)</p>
</div>

## Parameterization Maps {-}

:::{.blue-box}
We use UMAP [@McInnes2018a] to **reduce the dimension of our dataset** from $23$ (number of nodes with link operators) to $2$ and visualize it, to see if there is any apparent *visual relation* between the model parameterization and number of fixpoints.
:::

We used the [param_ss_umap.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_umap.R) script to run the UMAP implementation offered by the `uwot` R package.
We tried various values for the `n_neighbors` parameter where larger values result in **more global views** of the dataset, while smaller values result in **more local data** being preserved.
Also, the *distance metric* between the model parameterization vectors was mostly set to the standard (*euclidean*), but we also tried the *hamming* distance which seemed appropriate because of the binary nature of the dataset.

We make the figures afterwards using the result UMAP data with the [param_ss_umap_vis.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_umap_vis.R) script.
See all the produced figures [here](https://github.com/druglogics/bool-param-maps/tree/master/img/all_models_maps).

### Unsupervised UMAP {-}

First we used UMAP in **unsupervised mode** (no *a priori* knowledge of the number of fixpoints per model provided or of any other information/label per model for that matter).
So, UMAP is given all binary numbers from $0$ ($23$ $0$'s) to $2^{23}-1$ ($23$ $1$'s) representing each possible link operator mutated model ($0$'s map to `AND-NOT`, $1$'s to `OR-NOT`) and places them in the 2D plane.

The following figures show us the 2D parameterization map of all CASCADE 1.0 models, **colored either by their decimal (base-10) number** (converted from the binary link-operator model representation) or **by their respective number of fixpoints**:

```r
knitr::include_graphics(path = "img/all_models_maps/umap_20nn_model_num.png")
knitr::include_graphics(path = "img/all_models_maps/umap_20nn.png")
```

<div class="figure">
<img src="img/all_models_maps/umap_20nn_model_num.png" alt="CASCADE 1.0 Model Parameterization Maps (euclidean distance)" width="50%" /><img src="img/all_models_maps/umap_20nn.png" alt="CASCADE 1.0 Model Parameterization Maps (euclidean distance)" width="50%" />
<p class="caption">(\#fig:umap-unsup-figs)CASCADE 1.0 Model Parameterization Maps (euclidean distance)</p>
</div>


```r
knitr::include_graphics(path = "img/all_models_maps/ham_umap_20nn_model_num.png")
knitr::include_graphics(path = "img/all_models_maps/ham_umap_20nn.png")
```

<div class="figure">
<img src="img/all_models_maps/ham_umap_20nn_model_num.png" alt="CASCADE 1.0 Model Parameterization Maps (hamming distance)" width="50%" /><img src="img/all_models_maps/ham_umap_20nn.png" alt="CASCADE 1.0 Model Parameterization Maps (hamming distance)" width="50%" />
<p class="caption">(\#fig:umap-unsup-figs-2)CASCADE 1.0 Model Parameterization Maps (hamming distance)</p>
</div>

:::{.green-box}
- Using the *hamming* distance metric the visualization of the dataset results in many more smaller clusters compared to the *euclidean* representation, which results in $8$ **neighborhoods/superclusters of similarly parameterized models**.
These seem to follow the numerical representation (i.e. models with close decimal numbering seem to be clustered together).
- Models with **similar parameterization seem to also have the same number of fixpoints** (i.e. there is some order in the placement of models that belong to the same fixpoint class - this phenomenon does not manifest chaotically - e.g. models in the same sub-cluster in the hamming parameterization map tend to have the same number of fixpoints).
- There is no distinct pattern that can match model parameterization with number of attractors. 
In other words, models with different number of fixpoints can manifest in no particular order whatsoever across the parameterization map/space.
:::

### Supervised UMAP {-}

Next, we used **UMAP in supervised mode** - i.e. the association between each model and the corresponding fixpoint group was given as input to UMAP:

```r
knitr::include_graphics(path = "img/all_models_maps/sumap_14nn_0_3_min_dist.png")
knitr::include_graphics(path = "img/all_models_maps/sumap_20nn.png")
```

<div class="figure">
<img src="img/all_models_maps/sumap_14nn_0_3_min_dist.png" alt="CASCADE 1.0 Model Parameterization Supervised Maps" width="50%" /><img src="img/all_models_maps/sumap_20nn.png" alt="CASCADE 1.0 Model Parameterization Supervised Maps" width="50%" />
<p class="caption">(\#fig:umap-sup-fig)CASCADE 1.0 Model Parameterization Supervised Maps</p>
</div>

:::{.green-box}
With **increased model complexity** (meaning dynamical complexity - i.e. **more attractors**), the subsequent fixpoint superclusters form sub-clusters that are very differently parameterized - i.e. they form **distinct families of models** and thus appear to be more *spread out* in the supervised parameterization map.

The simple argument for this is as follows: the larger distance between a model in the $0$ or $1$-fixpoint (supervised) supercluster is smaller than many of the model distances within the $2$-fixpoint sub-clusters.
:::

## Embedding Important Nodes in the Map {-}

Using random forest and the regularized LASSO method, we found important nodes whose parameterization affects the change of dynamics (number of fixpoints) in the CASCADE 1.0 signaling network.
Using UMAP we observed that closely parameterized models form clusters.

We will now color the UMAP parameterization maps according to the link-operator values of the top $5$ most important nodes found from the aforementioned methods as well as the $2$ least important node reported with random forests (use the [param_ss_umap_imp_nodes.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/param_ss_umap_imp_nodes.R) script and see all the produced figures [here](https://github.com/druglogics/bool-param-maps/tree/master/img/imp_nodes_param_ss_maps)).

The 3 most important nodes:

```r
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_MAPK14.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_ERK_f.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_MEK_f.png")
```

<img src="img/imp_nodes_param_ss_maps/unsup_MAPK14.png" width="33%" /><img src="img/imp_nodes_param_ss_maps/unsup_ERK_f.png" width="33%" /><img src="img/imp_nodes_param_ss_maps/unsup_MEK_f.png" width="33%" />

The next 2 most important nodes:

```r
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_mTORC1_c.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_PTEN.png")
```

<img src="img/imp_nodes_param_ss_maps/unsup_mTORC1_c.png" width="50%" /><img src="img/imp_nodes_param_ss_maps/unsup_PTEN.png" width="50%" />

`CFLAR` and `CYCS` were the **least important node** for assessing the number of fixpoints of a model from its parameterization:

```r
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_CFLAR.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/unsup_CYCS.png")
```

<img src="img/imp_nodes_param_ss_maps/unsup_CFLAR.png" width="50%" /><img src="img/imp_nodes_param_ss_maps/unsup_CYCS.png" width="50%" />

Same trend can be seen (and maybe a little more clearly) in the supervised corresponding maps:

```r
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/sup_MAPK14.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/sup_mTORC1_c.png")
knitr::include_graphics(path = "img/imp_nodes_param_ss_maps/sup_CFLAR.png")
```

<img src="img/imp_nodes_param_ss_maps/sup_MAPK14.png" width="33%" /><img src="img/imp_nodes_param_ss_maps/sup_mTORC1_c.png" width="33%" /><img src="img/imp_nodes_param_ss_maps/sup_CFLAR.png" width="33%" />

:::{.green-box}
We can see a **visual link** between node importance (related to #fixpoints) and link operator assignment: **the less important a node is, the more randomly distributed (chaotically) it's link-operator values are across the parameterization map**.

A collection of important nodes can be used to more accurately define **families of closely parameterized models** and as we've seen above this also translates to models belonging to the same fixpoint class.
:::

# CASCADE 1.0 Analysis - 1 ss models {-}

## Stable States Data {-}

:::{.note}
In this section, **only the boolean models that have 1 stable state** will be used in the analysis.
All variables of interest (stable state, link-operator parameterization, fitness to steady state, performance MCC score, etc.) will relate only to the 1 stable state models from now on.

To load the stable state data for the models that have **1 stable state** use the Zenodo dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4022783.svg)](https://doi.org/10.5281/zenodo.4022783) and the script [get_ss_data.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/get_ss_data.R)
:::

## Parameterization Maps {-}

In this section we present the results of using UMAP [@McInnes2018a] on the link-operator parameterization data of the CASCADE 1.0 models with 1 stable state.
We created several such *parameterization maps* by adjusting the *n_neighbors* parameter input (from $2$ to $20$), which is responsible for the **size of the local neighborhood** (in terms of number of neighboring sample points) used for the manifold approximation.
As the documentation says, larger values result in **more global views** of the manifold, while smaller values result in **more local data** being preserved.
To get these map images and the reduced dimension dataset, use the script [1ss_models_umap.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/1ss_models_umap.R) for more details.

:::{.blue-box}
Note that in all these mappings to the 2D space, **models that share similar link-operator parameterization will reside in the same area/cluster in the map**.
The larger the *n_neighbors* is, the more the smaller clusters merge into larger ones.
The images for $\ge 14$ *n_neighbors* are almost exactly the same.
:::

We present some of these maps below:

```r
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_2.png")
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_3.png")
```

<div class="figure">
<img src="img/1ss_umap/1ss_umap_unsup_2.png" alt="2D Parameterization map for 1 stable state models (2 and 3 neighbors)" width="50%" /><img src="img/1ss_umap/1ss_umap_unsup_3.png" alt="2D Parameterization map for 1 stable state models (2 and 3 neighbors)" width="50%" />
<p class="caption">(\#fig:param-maps-1ss-models-1)2D Parameterization map for 1 stable state models (2 and 3 neighbors)</p>
</div>


```r
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_5.png")
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_6.png")
```

<div class="figure">
<img src="img/1ss_umap/1ss_umap_unsup_5.png" alt="2D Parameterization map for 1 stable state models (5 and 6 neighbors)" width="50%" /><img src="img/1ss_umap/1ss_umap_unsup_6.png" alt="2D Parameterization map for 1 stable state models (5 and 6 neighbors)" width="50%" />
<p class="caption">(\#fig:param-maps-1ss-models-2)2D Parameterization map for 1 stable state models (5 and 6 neighbors)</p>
</div>


```r
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_8.png")
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_11.png")
```

<div class="figure">
<img src="img/1ss_umap/1ss_umap_unsup_8.png" alt="2D Parameterization map for 1 stable state models (8 and 11 neighbors)" width="50%" /><img src="img/1ss_umap/1ss_umap_unsup_11.png" alt="2D Parameterization map for 1 stable state models (8 and 11 neighbors)" width="50%" />
<p class="caption">(\#fig:param-maps-1ss-models-3)2D Parameterization map for 1 stable state models (8 and 11 neighbors)</p>
</div>


```r
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_12.png")
knitr::include_graphics(path = "img/1ss_umap/1ss_umap_unsup_15.png")
```

<div class="figure">
<img src="img/1ss_umap/1ss_umap_unsup_12.png" alt="2D Parameterization map for 1 stable state models (12 and 15 neighbors)" width="50%" /><img src="img/1ss_umap/1ss_umap_unsup_15.png" alt="2D Parameterization map for 1 stable state models (12 and 15 neighbors)" width="50%" />
<p class="caption">(\#fig:param-maps-1ss-models-4)2D Parameterization map for 1 stable state models (12 and 15 neighbors)</p>
</div>

:::{.green-box}
We observe the existence of **two large families (superclusters) of parameterization models**, especially for more *global* views of the dataset ($\text{n_neighbors}\ge 8$).

**Distinct smaller sub-clusters of closely parameterized models** also manifest for different number of neighbors, which suggests that there is some order in the parameterization of the 1 stable state models (as exemplified by the UMAP method) across multiple visualization scales.
:::

## Gitsbe Models on the Map {-}

[Gitsbe](https://druglogics.github.io/druglogics-doc/gitsbe.html) uses a genetic algorithm approach to produce boolean models that are fitted to activity-based, biomarker training data.

We used Gitsbe and tested the produced models performance (ensemble-wise drug combination predictions) against synergy data from [@Flobak2015] in another report (AGS paper Sup. Material).
The calibrated models performed very well in terms of both ROC and PR-AUC.

:::{.blue-box}
Here we want to check whether models produced by a method such as a genetic algorithm-based one **have similar parameterization** - i.e. they belong in the same neighborhood in the parameterization map.
:::

We will use models from $1000$ gitsbe simulations, calibrated to steady state (a total of $3000$ models, choosing the $3$ best-fit models from each simulation).
The results are provided in [this data file](https://github.com/druglogics/bool-param-maps/blob/master/data/cascade_1.0_ss_1000sim_fixpoints_hsa.tar.gz) and to reproduce them, follow the instructions [here](https://druglogics.github.io/ags-paper/reproduce-data-simulation-results.html), keeping the default configuration options for CASCADE 1.0 and changing only the number of simulations to $1000$).

All the Gitsbe models had a large fitness to the steady state AGS data (their stable states fitting almost exactly the states of the manually-curated 24 nodes), as it can be seen from the next figure (see [gitsbe_models_fit.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/gitsbe_models_fit.R)):

```r
knitr::include_graphics(path = "img/gitsbe_fit_density.png")
```

<div class="figure">
<img src="img/gitsbe_fit_density.png" alt="Gitsbe model fitness to AGS steady state" width="1050" />
<p class="caption">(\#fig:gitbse-fit-fig)Gitsbe model fitness to AGS steady state</p>
</div>

To generate the next figures (same map, same gitsbe models, different number of neighbors) use the [gitsbe_model_embedding.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/gitsbe_model_embedding.R) - all figures are available [here](https://github.com/druglogics/bool-param-maps/tree/master/img/gitsbe_umaps):


```r
knitr::include_graphics(path = "img/gitsbe_umaps/2nn.png")
knitr::include_graphics(path = "img/gitsbe_umaps/4nn.png")
```

<div class="figure">
<img src="img/gitsbe_umaps/2nn.png" alt="Gitsbe models in Parameterization map (2 and 4 neighbors)" width="50%" /><img src="img/gitsbe_umaps/4nn.png" alt="Gitsbe models in Parameterization map (2 and 4 neighbors)" width="50%" />
<p class="caption">(\#fig:gitsbe-maps-1)Gitsbe models in Parameterization map (2 and 4 neighbors)</p>
</div>


```r
knitr::include_graphics(path = "img/gitsbe_umaps/8nn.png")
knitr::include_graphics(path = "img/gitsbe_umaps/14nn.png")
```

<div class="figure">
<img src="img/gitsbe_umaps/8nn.png" alt="Gitsbe models in Parameterization map (6 and 8 neighbors)" width="50%" /><img src="img/gitsbe_umaps/14nn.png" alt="Gitsbe models in Parameterization map (6 and 8 neighbors)" width="50%" />
<p class="caption">(\#fig:gitsbe-maps-2)Gitsbe models in Parameterization map (6 and 8 neighbors)</p>
</div>

:::{.green-box}
Gitsbe-generated models that fit the biomarker steady state data for the AGS cell line have a **diverse structure that spans across the parameterization map** but nonetheless appear to gather in **smaller parameterization-specific sub-clusters** (better seen in the Figure with 14 neighbors which gives a more global view of the dataset).

Observing the distribution of the gitsbe models in the parameterization map, we see that most of them are being **placed at one of the two superclusters**.
:::

Of course, there are areas in the map that Gitsbe models do not cover, which may as well be high-performance model areas.
Since we have generated all possible link-operator models with CASCADE 1.0, we can proceed to generate a performance map (atop the parameterization one) and cross-check if the gitsbe models fall into high-performance areas or not.

## Performance Maps {-}

Every model generated via `abmlog` ([Model Stable State Statistics]) was tested against the synergy dataset of [@Flobak2015].
Among $21$ drug combinations, $4$ were found synergistic in that dataset, namely:
  

```r
obs_syn = emba::get_observed_synergies(file = "data/observed_synergies_cascade_1.0")
usefun::pretty_print_vector_values(obs_syn, vector.values.str = "synergies")
```

> 4 synergies: PI-PD, PI-5Z, PD-AK, AK-5Z

:::{.blue-box}
Using the [drabme](https://druglogics.github.io/druglogics-doc/drabme.html) software module, we tested every CASCADE 1.0 model against this dataset and got each model's predictions for each drug combination.
The *HSA* rule was used to define if a model is synergistic for a particular combination.
The results are available in the Zenodo dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4022783.svg)](https://doi.org/10.5281/zenodo.4022783), file `cascade_1.0_hsa_fixpoints.tar.gz`.
As previously said, we are going to use the 1 stable state model predictions only.
:::

We used the emba R package [@Zobolas2020] to perform a biomarker analysis on the 1 stable state models dataset and their predictions from `drabme` (see script [emba_analysis.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/emba_analysis.R)).
Part of the results from the emba analysis is the calculation of the **Matthews correlation coefficient (MCC) performance score** for each model.
We use these MCC model scores to draw the next figures (see [mcc_figures.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/mcc_figures.R) script)

Splitting all the 1 stable state models to $4$ MCC classes we can see that **a lot of models have performance close to random prediction or even worse**:


```r
knitr::include_graphics(path = "img/mcc_hist.png")
knitr::include_graphics(path = "img/mcc_hist2.png")
```

<div class="figure">
<img src="img/mcc_hist.png" alt="MCC Classes Histograms" width="50%" /><img src="img/mcc_hist2.png" alt="MCC Classes Histograms" width="50%" />
<p class="caption">(\#fig:mcc-histogram)MCC Classes Histograms</p>
</div>

:::{.green-box}
Most of the 1 stable state models have MCC performance close to random or worse, making it thus difficult for any algorithm to find the *best* performance models (i.e. the 3-4th MCC class models, which amounts to only $5\%$ of the total number of models with 1 stable state, given the standardized boolean equation parameterization inspired by [@Mendoza2006]).
:::

### Unsupervised UMAP {-}

If we draw the parameterization maps for different number of neighbors and **color the points/models according to their MCC score**, we get these images (to see all figures, check [here](https://github.com/druglogics/bool-param-maps/tree/master/img/mcc_maps)):


```r
knitr::include_graphics(path = "img/mcc_maps/2nn.png")
knitr::include_graphics(path = "img/mcc_maps/14nn.png")
```

<div class="figure">
<img src="img/mcc_maps/2nn.png" alt="MCC Parameterization map (2 and 14 neighbors)" width="50%" /><img src="img/mcc_maps/14nn.png" alt="MCC Parameterization map (2 and 14 neighbors)" width="50%" />
<p class="caption">(\#fig:mcc-maps-1)MCC Parameterization map (2 and 14 neighbors)</p>
</div>

:::{.green-box}
We observe that the two large parameterization superclusters (especially outlined by the figures with $2$ and $\ge 8$ neighbors) are closely related with the MCC performance metric.
Specifically, these **2 superclusters dichotomize the models performance landscape** into 2 areas, where only one of them has the majority of good performance models (i.e. those that have an MCC score $>0$).

Also, we observe that closely parameterized models tend to have same performance (**existence of smaller parameterization clusters of same performance models**).
:::

:::{.orange-box}
Comparing the corresponding **Gitsbe parameterization maps with the MCC maps**, we can clearly verify now that the Gitsbe models **might not always be the best performing ones** (in terms of MCC score) since they appear in both superclusters (and that's only a first-fine measure for determining performance based on structure) - but at least they **tend to be more in the higher performance supercluster**.

Reasons why the gitsbe models are not always on the high-performance supercluster can be justified by:

- The nature of the genetic algorithm, which is a *heuristic* method and as such produces *local* optima
- The training dataset does not provide enough *constraints* for the optimization problem, i.e. more steady state nodes are needed to be fitted in the models stable state which will result in models with more strict parameterization patterns
- A combination of the two above
:::

### Supervised UMAP {-}

:::{.blue-box}
We perform **supervised dimension reduction** using UMAP on the 1 stable state model dataset. 
The only difference with the previous performance maps is that we provide UMAP (as an input) the information of the models performance (the MCC score) and not just overlay it (via coloring) afterwards.

The main thing we want to answer here is **if UMAP can do a better job placing the models on the map**, given additional information about their performance.
:::

We tried UMAP with different *number of neighbors* (from local to a more global view of the dataset), *target weight* (weighting factor between data topology ($0$) and target/response topology ($1$) or somewhere in between) and setting the response either to the MCC score as a continuous, numeric value or to the MCC classes (discrete variable) as [shown above](#fig:mcc-histogram).
See script [mcc_sumap.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/mcc_sumap.R) for more details.

First, we show some maps, where the response **MCC was treated as a continuous variable**:

```r
knitr::include_graphics(path = "img/mcc_sumaps/4nn_0w.png")
knitr::include_graphics(path = "img/mcc_sumaps/4nn_0.5w.png")
knitr::include_graphics(path = "img/mcc_sumaps/4nn_1w.png")
```

<div class="figure">
<img src="img/mcc_sumaps/4nn_0w.png" alt="MCC Parameterization map (4 neighbors, MCC continuous)" width="33%" /><img src="img/mcc_sumaps/4nn_0.5w.png" alt="MCC Parameterization map (4 neighbors, MCC continuous)" width="33%" /><img src="img/mcc_sumaps/4nn_1w.png" alt="MCC Parameterization map (4 neighbors, MCC continuous)" width="33%" />
<p class="caption">(\#fig:mcc-sumaps-1)MCC Parameterization map (4 neighbors, MCC continuous)</p>
</div>


```r
knitr::include_graphics(path = "img/mcc_sumaps/14nn_0.5w.png")
knitr::include_graphics(path = "img/mcc_sumaps/14nn_1w.png")
```

<div class="figure">
<img src="img/mcc_sumaps/14nn_0.5w.png" alt="MCC Parameterization map (14 neighbors, MCC continuous)" width="50%" /><img src="img/mcc_sumaps/14nn_1w.png" alt="MCC Parameterization map (14 neighbors, MCC continuous)" width="50%" />
<p class="caption">(\#fig:mcc-sumaps-2)MCC Parameterization map (14 neighbors, MCC continuous)</p>
</div>


```r
knitr::include_graphics(path = "img/mcc_sumaps/20nn_0.5w.png")
knitr::include_graphics(path = "img/mcc_sumaps/20nn_1w.png")
```

<div class="figure">
<img src="img/mcc_sumaps/20nn_0.5w.png" alt="MCC Parameterization map (20 neighbors, MCC continuous)" width="50%" /><img src="img/mcc_sumaps/20nn_1w.png" alt="MCC Parameterization map (20 neighbors, MCC continuous)" width="50%" />
<p class="caption">(\#fig:mcc-sumaps-3)MCC Parameterization map (20 neighbors, MCC continuous)</p>
</div>

:::{.green-box}
Very artistic and knotty figures!
UMAP achieves some sort of distinction between the different MCC values, especially when the topology is heavily based on the response MCC value (*target weight* is $1$)
:::

Next, we show some maps, where the response **MCC was treated as a categorical variable**, denoting the MCC class the models belong to (higher the better):

```r
knitr::include_graphics(path = "img/mcc_sumaps/class/4nn_0.5w_class.png")
knitr::include_graphics(path = "img/mcc_sumaps/class/6nn_0.5w_class.png")
```

<div class="figure">
<img src="img/mcc_sumaps/class/4nn_0.5w_class.png" alt="MCC Parameterization map (4 and 6 neighbors, MCC Classes)" width="50%" /><img src="img/mcc_sumaps/class/6nn_0.5w_class.png" alt="MCC Parameterization map (4 and 6 neighbors, MCC Classes)" width="50%" />
<p class="caption">(\#fig:mcc-sumaps-4)MCC Parameterization map (4 and 6 neighbors, MCC Classes)</p>
</div>


```r
knitr::include_graphics(path = "img/mcc_sumaps/class/8nn_0.5w_class.png")
knitr::include_graphics(path = "img/mcc_sumaps/class/10nn_0.5w_class.png")
```

<div class="figure">
<img src="img/mcc_sumaps/class/8nn_0.5w_class.png" alt="MCC Parameterization map (8 and 10 neighbors, MCC Classes)" width="50%" /><img src="img/mcc_sumaps/class/10nn_0.5w_class.png" alt="MCC Parameterization map (8 and 10 neighbors, MCC Classes)" width="50%" />
<p class="caption">(\#fig:mcc-sumaps-5)MCC Parameterization map (8 and 10 neighbors, MCC Classes)</p>
</div>


```r
knitr::include_graphics(path = "img/mcc_sumaps/class/14nn_0w_class.png")
knitr::include_graphics(path = "img/mcc_sumaps/class/14nn_0.5w_class.png")
knitr::include_graphics(path = "img/mcc_sumaps/class/14nn_1w_class.png")
```

<div class="figure">
<img src="img/mcc_sumaps/class/14nn_0w_class.png" alt="MCC Parameterization map (14 neighbors, MCC Classes)" width="33%" /><img src="img/mcc_sumaps/class/14nn_0.5w_class.png" alt="MCC Parameterization map (14 neighbors, MCC Classes)" width="33%" /><img src="img/mcc_sumaps/class/14nn_1w_class.png" alt="MCC Parameterization map (14 neighbors, MCC Classes)" width="33%" />
<p class="caption">(\#fig:mcc-sumaps-6)MCC Parameterization map (14 neighbors, MCC Classes)</p>
</div>

:::{.green-box}
We observe that the higher the number of neighbors is ($\ge 10$ with a balanced *target weight* of $0.5$), the better UMAP classifies the models to distinct (super-) clusters representing the different MCC classes.
:::

## Performance vs Fitness {-}

:::{.blue-box}
We try to see if there is any correlation between **performance (MCC)** and **fitness to the AGS steady state**.

We study **two AGS steady states**: one with a total of $24$ nodes which is the result of manual literature curation (S4 Table in @Flobak2015) and serves as the *gold-standard*, and one with all $77$ nodes, which is the single (and only) fixpoint of the AGS boolean model in @Flobak2015.

Note that second, larger steady state reports the same activity values for the common $24$ nodes as the first one.

Use the [fit_figures.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/fit_figures.R) script to reproduce the figures for the literature-curated steady state and the [fit_vs_performance_ags_pub.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/fit_vs_performance_ags_pub.R) for the AGS model fixpoint/steady state.
:::

First, we check if the fitness density of all 1 stable state models covers the whole *fitness spectrum*:

```r
knitr::include_graphics(path = "img/1ss_models_fit_density.png")
knitr::include_graphics(path = "img/1ss_models_fit_density_pub.png")
```

<div class="figure">
<img src="img/1ss_models_fit_density.png" alt="All 1 stable state models fitness to AGS steady state(s)" width="50%" /><img src="img/1ss_models_fit_density_pub.png" alt="All 1 stable state models fitness to AGS steady state(s)" width="50%" />
<p class="caption">(\#fig:1ss-models-fit-density-fig)All 1 stable state models fitness to AGS steady state(s)</p>
</div>

The Pearson correlation between the two fitness vectors is high:

```r
fit_data = readRDS(file = "data/fit_data.rds")
fit_data_pub = readRDS(file = "data/fit_data_pub.rds")

# data check
stopifnot(all(names(fit_data) == names(fit_data_pub)))

# Pearson correlation
cor(x = fit_data, y = fit_data_pub, method = "pearson") %>% round(digits = 2)
```

[1] 0.94

:::{.green-box}
We observe that **most of the models have at least half of the nodes** in the same state as in the AGS steady state, no matter which steady state we fit the stable state data against.
This makes the **fitness score distribution skewed** towards the higher fitness values.

Of course, though the two steady states completely agree on the values of the $24$ nodes, it might be that the rest of the nodes have slightly different activity state values from the ones found in the fixpoint attractor of the AGS model in @Flobak2015.
If that is the case in reality, the fitness distribution might be more uniform.
:::

We follow the same classification scheme as [above](#fig:mcc-histogram), i.e. splitting the 1 stable models to $4$ MCC classes and comparing the fitness scores between these classes:

```r
knitr::include_graphics(path = "img/mcc_vs_fit.png")
knitr::include_graphics(path = "img/mcc_vs_fit_pub.png")
```

<div class="figure">
<img src="img/mcc_vs_fit.png" alt="MCC performance vs Fitness to AGS steady state(s) for all 1 stable state models. Left is fitness to the AGS-curated steady state (24 nodes), right is for the full, 77-node steady state." width="50%" /><img src="img/mcc_vs_fit_pub.png" alt="MCC performance vs Fitness to AGS steady state(s) for all 1 stable state models. Left is fitness to the AGS-curated steady state (24 nodes), right is for the full, 77-node steady state." width="50%" />
<p class="caption">(\#fig:fit-vs-perf-fig)MCC performance vs Fitness to AGS steady state(s) for all 1 stable state models. Left is fitness to the AGS-curated steady state (24 nodes), right is for the full, 77-node steady state.</p>
</div>

:::{.green-box}
The last MCC class has a **statistically significant higher median fitness** to the AGS steady state compared to all other lower MCC classes, even though the correlation between fitness and performance is not linear and most of the models have a fitness higher than $0.5$.
Note also that models that have significantly lower than $0.5$ fitness in each steady state case, belong to the first two, lower performance classes.

The last two figures points us to the fact that more constraints are needed for the fitness calculation of the Gitsbe genetic algorithm or any other for that matter (i.e. more literature-curated nodes in the AGS steady state - now only $24/77=31\%$ is included) in order to define more restrictive parameterized models that would allow a much more *uniform* fitness density spectrum (or at least **not skewed towards the higher fitness values**).
Such fitness spectrum would (hopefully) allow for more granularity in the corresponding performance behavior between the different MCC classes, and thus more distinctive correlation.
:::

## Fitness Maps {-}

If we draw the parameterization maps for different number of neighbors and **color the points/models according to their fitness to the AGS steady state (24 nodes)**, we get these images (see [fit_figures.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/fit_figures.R) script):


```r
knitr::include_graphics(path = "img/fit_maps/2nn.png")
knitr::include_graphics(path = "img/fit_maps/4nn.png")
```

<div class="figure">
<img src="img/fit_maps/2nn.png" alt="Fitness Parameterization map (2 and 4 neighbors)" width="50%" /><img src="img/fit_maps/4nn.png" alt="Fitness Parameterization map (2 and 4 neighbors)" width="50%" />
<p class="caption">(\#fig:fit-maps-1)Fitness Parameterization map (2 and 4 neighbors)</p>
</div>


```r
knitr::include_graphics(path = "img/fit_maps/11nn.png")
knitr::include_graphics(path = "img/fit_maps/14nn.png")
```

<div class="figure">
<img src="img/fit_maps/11nn.png" alt="Fitness Parameterization map (11 and 14 neighbors)" width="50%" /><img src="img/fit_maps/14nn.png" alt="Fitness Parameterization map (11 and 14 neighbors)" width="50%" />
<p class="caption">(\#fig:fit-maps-2)Fitness Parameterization map (11 and 14 neighbors)</p>
</div>

:::{.green-box}
Higher fitness models manifest in both superclusters, suggesting again the need for more literature-curated nodes in the training data (AGS steady state).
Also, closely parameterized models tend to have same fitness (**existence of smaller parameterization clusters of same fitness models**).

No apparent correlation can be observed between fitness and performance (MCC) maps.
:::

## Performance biomarkers {-}

:::{.blue-box}
We assess **important nodes (biomarkers)** whose *activity* and/or link-operator (*parameterization*) affects the models performance in terms of the achieved MCC score.
:::

### emba biomarkers {-}

We use the `emba` results with $4$ MCC classes (see script [emba_analysis.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/emba_analysis.R) and MCC classification histograms [above](#fig:mcc-histogram)).
What emba [@Zobolas2020] does is to group the models to $4$ MCC Classes and compare the average stable state activities and link operator values for each node between **every possible pair of MCC classes**.
Nodes that have **higher average differences** (either positively or negatively and between any comparison pair) are considered as **more important** and are thus annotated as *biomarkers* based on a given threshold.

In the next heatmap, a positive (resp. negative) state difference for a node denotes that its activity value was larger (resp. lower) in the higher MCC class of the corresponding comparison pair:

```r
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds")

# color some important nodes and lower the threshold
activity_biomarkers = emba::get_biomarkers(mcc_res$diff.state.mcc.mat, threshold = 0.48)

imp_nodes_colors = rep("black", ncol(mcc_res$diff.state.mcc.mat))
names(imp_nodes_colors) = colnames(mcc_res$diff.state.mcc.mat)
imp_nodes_colors[names(imp_nodes_colors) %in% activity_biomarkers$biomarkers.pos] = "green4"
imp_nodes_colors[names(imp_nodes_colors) %in% activity_biomarkers$biomarkers.neg] = "red4"

# define coloring function of state differences
col_fun = circlize::colorRamp2(c(min(mcc_res$diff.state.mcc.mat), 0, max(mcc_res$diff.state.mcc.mat)), 
  c("red", "white", "green"))

set.seed(42)
state_heat = ComplexHeatmap::Heatmap(matrix = mcc_res$diff.state.mcc.mat, 
  name = "State Difference", col = col_fun,
  row_title = "MCC Classes Pairs", row_names_side = "left", 
  row_title_side = "left", row_dend_side = "right",
  column_title = "Average State Differences between MCC Classes", 
  column_names_gp = gpar(fontsize = 6, col = imp_nodes_colors),
  heatmap_legend_param = list(direction = "horizontal"))

biomarkers_legend = Legend(title = "Activity State Biomarkers",
  labels = c("Active", "Inhibited"),
  legend_gp = gpar(fill = c("green4", "red4")))

draw(state_heat, annotation_legend_list = biomarkers_legend, 
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
  merge_legend = TRUE)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/emba-activity-biomarkers-1.png" alt="Heatmap of Average State Differences (emba, 4 MCC Classes)" width="672" />
<p class="caption">(\#fig:emba-activity-biomarkers)Heatmap of Average State Differences (emba, 4 MCC Classes)</p>
</div>

So, the following nodes are identified as *active* biomarkers (they tend to be active in the higher performance models):

```r
usefun::pretty_print_vector_values(activity_biomarkers$biomarkers.pos)
```

> 7 nodes: MAP3K7, GRAP2, MAP2K7, MAP2K3, TAB_f, ERK_f, NLK

and the following as *inhibited* biomarkers:

```r
usefun::pretty_print_vector_values(activity_biomarkers$biomarkers.neg)
```

> 1 node: MAPK14

Note that some nodes tend to have either no actual state differences between the MCC classes (so they are not important) or they take both negative and positive differences between the different MCC class pair comparisons (without surpassing the given threshold, so also not important for us here).


```r
# color some important nodes and lower the threshold
param_biomarkers = emba::get_biomarkers(mcc_res$diff.link.mcc.mat, threshold = 0.48)

imp_nodes_colors = rep("black", ncol(mcc_res$diff.link.mcc.mat))
names(imp_nodes_colors) = colnames(mcc_res$diff.link.mcc.mat)
imp_nodes_colors[names(imp_nodes_colors) %in% param_biomarkers$biomarkers.pos] = "green4"
imp_nodes_colors[names(imp_nodes_colors) %in% param_biomarkers$biomarkers.neg] = "red4"

# define coloring function of state differences
col_fun = circlize::colorRamp2(c(min(mcc_res$diff.link.mcc.mat), 0, max(mcc_res$diff.link.mcc.mat)), 
  c("red", "white", "green"))

set.seed(42)
state_heat = ComplexHeatmap::Heatmap(matrix = mcc_res$diff.link.mcc.mat, 
  name = "LO Difference", col = col_fun,
  row_title = "MCC Classes Pairs", row_names_side = "left", 
  row_title_side = "left", row_dend_side = "right",
  column_title = "Average Link Operator Differences between MCC Classes", 
  column_names_gp = gpar(fontsize = 11, col = imp_nodes_colors),
  heatmap_legend_param = list(direction = "horizontal"))

biomarkers_legend = Legend(title = "Link Operator Biomarkers",
  labels = c("OR-NOT (1)", "AND-NOT (0)"),
  legend_gp = gpar(fill = c("green4", "red4")))

draw(state_heat, annotation_legend_list = biomarkers_legend, 
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
  merge_legend = TRUE)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/emba-lo-biomarkers-1.png" alt="Heatmap of Average Link Operator Differences (emba, 4 MCC Classes)" width="672" />
<p class="caption">(\#fig:emba-lo-biomarkers)Heatmap of Average Link Operator Differences (emba, 4 MCC Classes)</p>
</div>

:::{.green-box}
We observe that `ERK_f` will have the `OR-NOT` link operator in most of the higher performance models and the `MAPK14` the `AND-NOT` link operator.
This of course relates to the fact that these nodes were found also as *active* and *inhibited* biomarkers respectively and that they have a very large observed agreement between stable state activity value and link operator parameterization (see [analysis here](https://druglogics.github.io/brf-bias/cascade-1-0-data-analysis.html#fig:ss-lo-agreement-prop)).

Interestingly, these two nodes (`ERK_f` and `MAPK14`) were 2 of the top most important nodes influencing the change of dynamics (number of attractors) in the link operator parameterization space of the CASCADE 1.0 network.
:::

### Random Forest biomarkers {-}

We use the `ranger` R package [@Wright2017] to find **important nodes/variables** that determine the difference in performance (MCC score) between the input models.
Both the stable state data as well the link operator parameterization data will be used as training sets (see the script [perf_biomarkers_ranger.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/perf_biomarkers_ranger.R)).


```r
ranger_ss_mcc_res = readRDS(file = "data/ss_mcc_ranger_res.rds")

imp_ss_rf = tibble(nodes = names(ranger_ss_mcc_res$variable.importance), 
  gini_index = ranger_ss_mcc_res$variable.importance)

# color first 12 nodes in x axis
imp_col = c(rep("green4", 14), rep("grey30", nrow(imp_ss_rf) - 14))

imp_ss_rf %>% 
  mutate(nodes = forcats::fct_reorder(nodes, desc(gini_index))) %>%
  ggplot(aes(x = nodes, y = gini_index, fill = gini_index)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_gradient(low = "steelblue", high = "red") +
    theme_classic(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 90, size = 7, colour = imp_col)) +
    labs(title = "RF Importance for MCC Classification (stable states)", 
      x = "Nodes", y = "Mean Decrease in Node Impurity")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/rf-mcc-state-imp-1.png" alt="Random Forest Important Nodes for Performance Classification (MCC) based on stable state data" width="2100" />
<p class="caption">(\#fig:rf-mcc-state-imp)Random Forest Important Nodes for Performance Classification (MCC) based on stable state data</p>
</div>

:::{.green-box}
We observe a lot of common important activity biomarkers (important nodes) between the random forest results and the `emba` results: `ERK_f`,`MAP3K7`,`MAP2K7`,`NLK`,`TAB_f` to name a few.
:::


```r
ranger_lo_mcc_res = readRDS(file = "data/lo_mcc_ranger_res.rds")

imp_lo_rf = tibble(nodes = names(ranger_lo_mcc_res$variable.importance), 
  gini_index = ranger_lo_mcc_res$variable.importance)

# color first 8 nodes in x axis
imp_col = c(rep("green4", 8), rep("grey30", nrow(imp_lo_rf) - 8))

imp_lo_rf %>% 
  mutate(nodes = forcats::fct_reorder(nodes, desc(gini_index))) %>%
  ggplot(aes(x = nodes, y = gini_index, fill = gini_index)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_gradient(low = "steelblue", high = "red") +
    theme_classic(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 90, size = 12, colour = imp_col)) +
    labs(title = "RF Importance for MCC Classification (parameterization)", 
      x = "Nodes", y = "Mean Decrease in Node Impurity")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/rf-mcc-param-imp-1.png" alt="Random Forest Important Nodes for Performance Classification (MCC) based on link operator data)" width="2100" />
<p class="caption">(\#fig:rf-mcc-param-imp)Random Forest Important Nodes for Performance Classification (MCC) based on link operator data)</p>
</div>

:::{.green-box}
`ERK_f` is a major link operator biomarker in both `emba` and random forest results
:::

### emba vs RF Results {-}

We check the correlation between `emba` (sum of absolute state differences for all comparison pairs per node, i.e. summing up the absolute values of the columns in the heatmaps above) and random forest (Gini Index) results:

```r
# compute total emba importance
imp_ss_emba = mcc_res$diff.state.mcc.mat %>% abs() %>% colSums()
imp_lo_emba = mcc_res$diff.link.mcc.mat %>% abs() %>% colSums()

# make sure it's in the same order
imp_ss_emba = imp_ss_emba[imp_ss_rf %>% pull(nodes)]
imp_lo_emba = imp_lo_emba[imp_lo_rf %>% pull(nodes)]

imp_ss_rf %>% add_column(emba_score = imp_ss_emba) %>%
  ggplot(aes(x = emba_score, y = gini_index)) +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm") + 
    ggpubr::stat_cor(method = "kendall", cor.coef.name = "tau") +
    theme_classic() +
    labs(title = "Random Forest vs emba variable importance results (Activity State)",
      x = "emba (Sum of absolute state differences)", y = "RF (Gini Index)")

imp_lo_rf %>% add_column(emba_score = imp_lo_emba) %>%
  ggplot(aes(x = emba_score, y = gini_index)) +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm") + 
    ggpubr::stat_cor(method = "kendall", cor.coef.name = "tau") +
    theme_classic() +
    labs(title = "Random Forest vs emba variable importance results (Parameterization)",
      x = "emba (Sum of absolute link operator differences)", y = "RF (Gini Index)")
```

<div class="figure">
<img src="index_files/figure-html/cor-emba-rf-1.png" alt="Correlation between emba and Random Forest Importance Results" width="50%" /><img src="index_files/figure-html/cor-emba-rf-2.png" alt="Correlation between emba and Random Forest Importance Results" width="50%" />
<p class="caption">(\#fig:cor-emba-rf)Correlation between emba and Random Forest Importance Results</p>
</div>

:::{.green-box}
- We use the non-parametric Kendall correlation method and it's rank-based correlation coefficient, $\tau$.
- There is **a correlation between the importance results from the two different methods**, even though the `emba` results are more *artificially-made* (like adding up absolute coefficient scores from linear models as a measure of importance) compared to the impurity metric which is a *built-in* feature of the tree-based algorithms such as random forests.
As a result for example, we observe a much greater granularity of the importance measurements in the activity state correlation figure.
- `emba` results are **more verbose** though in the sense that they allow *directionality* in their interpretation: we get which nodes must be activated or inhibited in the higher performance models or which should be parameterized with `OR-NOT` or `AND-NOT`, information which the random forests do not provide.
- **Larger correlation** is observed **for the activity state** compared to the parameterization results.
:::

### Embedding Link Operator Biomarkers in the Map {-}

:::{.blue-box}
Using as a base the **performance parameterization maps (supervised and unsupervised)** for the 1 stable state models, we color the points (models) according to the **link operator value** of some of the **performance biomarker** that we found in the previous analysis using `emba` and random forests.
:::

We give some examples of how the distribution of link-operator values looks like in the parameterization maps for the more important nodes (e.g. `ERK_f`, `MAPK14`) but also for the least important ones (e.g. `CYCS`, `CFLAR`).
See the script [perf_biomarkers_embedding.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/perf_biomarkers_embedding.R) for more details.
All the produced images by the script are accessible [here](https://github.com/druglogics/bool-param-maps/tree/master/img/nodes_lo_maps).

Using as a base the **unsupervised parameterization MAP with 14 neighbors**, we have:

```r
knitr::include_graphics(path = "img/nodes_lo_maps/unsup_MAPK14.png")
knitr::include_graphics(path = "img/nodes_lo_maps/unsup_ERK_f.png")
knitr::include_graphics(path = "img/nodes_lo_maps/unsup_CFLAR.png")
knitr::include_graphics(path = "img/nodes_lo_maps/unsup_CYCS.png")
```

<div class="figure">
<img src="img/nodes_lo_maps/unsup_MAPK14.png" alt="MCC Unsupervised Parameterization maps colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/unsup_ERK_f.png" alt="MCC Unsupervised Parameterization maps colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/unsup_CFLAR.png" alt="MCC Unsupervised Parameterization maps colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/unsup_CYCS.png" alt="MCC Unsupervised Parameterization maps colored by link-operator values of different performance biomarkers" width="50%" />
<p class="caption">(\#fig:nodes-lo-maps-1)MCC Unsupervised Parameterization maps colored by link-operator values of different performance biomarkers</p>
</div>

Using as a base the **supervised parameterization MAP with 14 neighbors (MCC as continuous variable)**, we have:

```r
knitr::include_graphics(path = "img/nodes_lo_maps/sup_cont_MAPK14.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_cont_ERK_f.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_cont_CFLAR.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_cont_CYCS.png")
```

<div class="figure">
<img src="img/nodes_lo_maps/sup_cont_MAPK14.png" alt="MCC Supervised Parameterization maps (MCC continuous) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_cont_ERK_f.png" alt="MCC Supervised Parameterization maps (MCC continuous) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_cont_CFLAR.png" alt="MCC Supervised Parameterization maps (MCC continuous) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_cont_CYCS.png" alt="MCC Supervised Parameterization maps (MCC continuous) colored by link-operator values of different performance biomarkers" width="50%" />
<p class="caption">(\#fig:nodes-lo-maps-2)MCC Supervised Parameterization maps (MCC continuous) colored by link-operator values of different performance biomarkers</p>
</div>

Using as a base the **supervised parameterization MAP with 14 neighbors (MCC as discrete class variable)**, we have:

```r
knitr::include_graphics(path = "img/nodes_lo_maps/sup_MAPK14.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_ERK_f.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_CFLAR.png")
knitr::include_graphics(path = "img/nodes_lo_maps/sup_CYCS.png")
```

<div class="figure">
<img src="img/nodes_lo_maps/sup_MAPK14.png" alt="MCC Supervised Parameterization maps (MCC discrete) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_ERK_f.png" alt="MCC Supervised Parameterization maps (MCC discrete) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_CFLAR.png" alt="MCC Supervised Parameterization maps (MCC discrete) colored by link-operator values of different performance biomarkers" width="50%" /><img src="img/nodes_lo_maps/sup_CYCS.png" alt="MCC Supervised Parameterization maps (MCC discrete) colored by link-operator values of different performance biomarkers" width="50%" />
<p class="caption">(\#fig:nodes-lo-maps-3)MCC Supervised Parameterization maps (MCC discrete) colored by link-operator values of different performance biomarkers</p>
</div>

:::{.green-box}
We observe that both in the unsupervised parameterization maps and the MCC-supervised ones (especially the continuous one), **the less important a node is, the more chaotically it's link operator value manifests across the map**.
:::

## Synergy Maps {-}

As stated in a [previous section](#performance-maps), the CASCADE 1.0 models produced by `abmlog` were tested for synergy against the drug combination dataset in @Flobak2015.
Using the *HSA* method to define a drug combination as synergistic or not (antagonistic), we first present some useful statistics (see script [synergy_maps.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/synergy_maps.R)):

```r
knitr::include_graphics(path = 'img/synergy_stats.png')
```

<div class="figure">
<img src="img/synergy_stats.png" alt="Synergy Statistics across all CASCADE 1.0 models with 1 stable state" width="1050" />
<p class="caption">(\#fig:syn-stats)Synergy Statistics across all CASCADE 1.0 models with 1 stable state</p>
</div>

:::{.orange-box}
For the 'NA' models, either one of the two single-drug perturbed models or the double-perturbed model, did not have any fixpoint attractors.
:::

We observe that **only a very few models** (compared to the total parameterization space) can predict the observed drug combinations as *synergistic*.
This can be visualized also in the next synergy maps:


```r
knitr::include_graphics(path = "img/synergy_maps/AK-5Z_2nn.png")
knitr::include_graphics(path = "img/synergy_maps/AK-5Z_14nn.png")
```

<div class="figure">
<img src="img/synergy_maps/AK-5Z_2nn.png" alt="AK-PZ Synergy Parameterization Maps" width="50%" /><img src="img/synergy_maps/AK-5Z_14nn.png" alt="AK-PZ Synergy Parameterization Maps" width="50%" />
<p class="caption">(\#fig:syn-maps-1)AK-PZ Synergy Parameterization Maps</p>
</div>


```r
knitr::include_graphics(path = "img/synergy_maps/PD-AK_2nn.png")
knitr::include_graphics(path = "img/synergy_maps/PD-AK_14nn.png")
```

<div class="figure">
<img src="img/synergy_maps/PD-AK_2nn.png" alt="PD-AK Synergy Parameterization Maps" width="50%" /><img src="img/synergy_maps/PD-AK_14nn.png" alt="PD-AK Synergy Parameterization Maps" width="50%" />
<p class="caption">(\#fig:syn-maps-2)PD-AK Synergy Parameterization Maps</p>
</div>


```r
knitr::include_graphics(path = "img/synergy_maps/PI-PD_2nn.png")
knitr::include_graphics(path = "img/synergy_maps/PI-PD_14nn.png")
```

<div class="figure">
<img src="img/synergy_maps/PI-PD_2nn.png" alt="PI-PD Synergy Parameterization Maps" width="50%" /><img src="img/synergy_maps/PI-PD_14nn.png" alt="PI-PD Synergy Parameterization Maps" width="50%" />
<p class="caption">(\#fig:syn-maps-3)PI-PD Synergy Parameterization Maps</p>
</div>


```r
knitr::include_graphics(path = "img/synergy_maps/PI-5Z_2nn.png")
knitr::include_graphics(path = "img/synergy_maps/PI-5Z_14nn.png")
```

<div class="figure">
<img src="img/synergy_maps/PI-5Z_2nn.png" alt="PI-PZ Synergy Parameterization Maps" width="50%" /><img src="img/synergy_maps/PI-5Z_14nn.png" alt="PI-PZ Synergy Parameterization Maps" width="50%" />
<p class="caption">(\#fig:syn-maps-4)PI-PZ Synergy Parameterization Maps</p>
</div>

:::{.green-box}
As we have seen, the UMAP results with higher number of neighbors cluster the similarly parameterized models better.
This makes it easier to spot the **small clusters** which include models that can predict each of the respective observed synergies as *synergistic*.

We note that all of these **synergistic sub-clusters** are part of the **high-performance supercluster** (right one in the images above).
:::

Interestingly, if we see consider **all possible observed synergy subsets** we see that there are models that can predict the double synergy set `PI-PD,PD-AK` and the `PI-5Z,AK-5Z`, but there is literally **no model that can predict three or four (all) of the observed synergies**:

```r
synergy_subset_stats = readRDS(file = "data/synergy_subset_stats.rds")
usefun::pretty_print_vector_names_and_values(synergy_subset_stats)
```

> : 2675360, PI-PD: 59616, PI-5Z: 31968, PD-AK: 1536, AK-5Z: 85248, PI-PD,PD-AK: 1536, PI-5Z,AK-5Z: 31968

This is another proof that this parameterization is *restrictive* since there is not even one possible parameterization which represents biological reality (i.e. treatment with the previously mentioned four drug combinations is synergistic). 

Also, looking at the numbers of models for each set of synergies, we see that the models predicting the `PD-AK` synergy form a proper subset of the `PI-PD,PD-AK` models and the same for `PI-5Z` and `PI-5Z,AK-5Z`.
Thus we include only the parent sets in the next synergy maps, where we have combined all synergies in one:


```r
knitr::include_graphics(path = "img/synergy_maps/all_syn_11nn.png")
knitr::include_graphics(path = "img/synergy_maps/all_syn_14nn.png")
```

<div class="figure">
<img src="img/synergy_maps/all_syn_11nn.png" alt="Combined Synergy Parameterization Maps" width="50%" /><img src="img/synergy_maps/all_syn_14nn.png" alt="Combined Synergy Parameterization Maps" width="50%" />
<p class="caption">(\#fig:syn-maps-5)Combined Synergy Parameterization Maps</p>
</div>

:::{.green-box}
We observe that there are **two synergy sub-clusters**, one primarily related to the `PD` drug and the other to the `5Z` drug.
Since the synergies belong to $2$ mutually exclusive (in terms of parameterization) clusters, this is a visual cue that there cannot be a model that predicts all 4 of these synergies.
:::

## Stable State Patterns in Synergistic models {-}

:::{.blue-box}
We want to identify if there are **activity patterns** in the models that predict the observed synergies.

To do this, we take a random sample of such synergistic models that predict each of the $4$ observed synergies and show the corresponding heatmaps of stable state activities.
More details on the [synergy_heatmaps.R](https://github.com/druglogics/bool-param-maps/blob/master/scripts/synergy_heatmaps.R) script.
:::


```r
knitr::include_graphics(path = "img/synergy_heatmaps/PD-AK_ss_heat.png")
```

<div class="figure">
<img src="img/synergy_heatmaps/PD-AK_ss_heat.png" alt="PD-AK stable state activity heatmap" width="2100" />
<p class="caption">(\#fig:synergy-heatmaps)PD-AK stable state activity heatmap</p>
</div>


```r
knitr::include_graphics(path = "img/synergy_heatmaps/PI-PD_ss_heat.png")
```

<div class="figure">
<img src="img/synergy_heatmaps/PI-PD_ss_heat.png" alt="PI-PD stable state activity heatmap" width="2100" />
<p class="caption">(\#fig:synergy-heatmaps-2)PI-PD stable state activity heatmap</p>
</div>

:::{.green-box}
- All the models that predict the `PD-AK` synergy show one **distinct activity state pattern**
- The `PI-PD` heatmap is more *heterogeneous* (though still patterns do exist).
The reason for this might be because the parameterization is more spread out across the map in the combined synergy [figure above](#fig:syn-maps-5) - i.e. some models are in the low-performance supercluster.
:::


```r
knitr::include_graphics(path = "img/synergy_heatmaps/AK-5Z_ss_heat.png")
```

<div class="figure">
<img src="img/synergy_heatmaps/AK-5Z_ss_heat.png" alt="AK-5Z stable state activity heatmaps" width="2100" />
<p class="caption">(\#fig:synergy-heatmaps-3)AK-5Z stable state activity heatmaps</p>
</div>


```r
knitr::include_graphics(path = "img/synergy_heatmaps/PI-5Z_ss_heat.png")
```

<div class="figure">
<img src="img/synergy_heatmaps/PI-5Z_ss_heat.png" alt="PI-5Z stable state activity heatmaps" width="2100" />
<p class="caption">(\#fig:synergy-heatmaps-4)PI-5Z stable state activity heatmaps</p>
</div>

:::{.green-box}
All the synergistic models in the `5Z` sub-cluster seem to follow the **same stable state activity patterns** (pretty much we get the same node names in the columns after the clustering in the heatmaps). 
There exist $2$ **main such patterns**, specified by the state vector (ON,{ON/OFF},{OFF/ON},ON) - where ON and OFF denote vectors of *active*/*inhibited* nodes that are clustered together.
:::

Note that in all $4$ heatmaps above, the `Prosurvival` node is always *active*, denoting proliferating activity patterns in the boolean models where synergies can manifest.

## Synergy Biomarkers {-}

:::{.blue-box}
We assess **important nodes (biomarkers)** whose *activity* and/or link-operator (*parameterization*) affects the manifestation of the $4$ observed synergies.
:::

Here, `emba` [@Zobolas2020] performs the classification to *synergistic* and *antagonistic* model groups per synergy observed (as exemplified by the [figure above](#fig:syn-stats) - the `NA` models are discarded) and compares the average stable state activities and link operator values of each node in the two groups.

In the next heatmap, a positive (resp. negative) state difference for a node denotes that its activity value was larger (resp. lower) in the corresponding synergistic model group:

```r
syn_res = readRDS(file = "data/synergy_res.rds")

# define coloring function of state differences
col_fun = circlize::colorRamp2(c(min(syn_res$diff.state.synergies.mat), 
  0, max(syn_res$diff.state.synergies.mat)), c("red", "white", "green"))

set.seed(42)
state_syn_heat = ComplexHeatmap::Heatmap(matrix = syn_res$diff.state.synergies.mat, 
  name = "State Difference", col = col_fun,
  row_title = "Synergistic Drug Combinations", row_names_side = "left", 
  row_title_side = "left", row_dend_side = "right",
  column_title = "Average State Differences (Synergistic vs Antagonistic)", 
  column_names_gp = gpar(fontsize = 5),
  heatmap_legend_param = list(direction = "horizontal"))

biomarkers_legend = Legend(title = "Activity State Biomarkers",
  labels = c("Active", "Inhibited"),
  legend_gp = gpar(fill = c("green4", "red4")))

draw(state_syn_heat, annotation_legend_list = biomarkers_legend, 
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
  merge_legend = TRUE)
```

<div class="figure">
<img src="index_files/figure-html/state-syn-heat-1.png" alt="Heatmap of Average State Differences between synergistic and antagonistic model groups for each observed synergy" width="672" />
<p class="caption">(\#fig:state-syn-heat)Heatmap of Average State Differences between synergistic and antagonistic model groups for each observed synergy</p>
</div>

:::{.green-box}
- The activity state biomarkers are clustered in two categories representing the two synergy sub-clusters we found previously, namely the `5Z` sub-cluster and the `PD` sub-cluster.
- The models predicting the `5Z` synergies seem to have more **active state biomarkers**, whereas the models predicting the `PD` synergies seem to have more **inhibited state biomarkers**.
- `MAPK14` seems to be a **definitive node distinguishing the two synergy subgroups**, since it's a strong inhibited biomarker for the `5Z` synergies, but a fairly unimportant node for the `PD` synergies.
- `ERK_f` is an active state biomarker (same observation was found in the *performance* biomarkers analysis)
:::

In the next heatmap, a positive (resp. negative) link operator difference for a node denotes that its link operator was mostly `OR-NOT` (resp. `AND-NOT`) in the corresponding synergistic model group:

```r
# define coloring function of state differences
col_fun = circlize::colorRamp2(c(min(syn_res$diff.link.synergies.mat), 
  0, max(syn_res$diff.link.synergies.mat)), c("red", "white", "green"))

set.seed(42)
lo_syn_heat = ComplexHeatmap::Heatmap(matrix = syn_res$diff.link.synergies.mat, 
  name = "Link Operator Difference", col = col_fun,
  row_title = "Synergistic Drug Combinations", row_names_side = "left", 
  row_title_side = "left", row_dend_side = "right",
  column_title = "Average Link Operator Differences (Synergistic vs Antagonistic)", 
  column_names_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(direction = "horizontal"))

biomarkers_legend = Legend(title = "Link Operator Biomarkers",
  labels = c("OR-NOT (1)", "AND-NOT (0)"),
  legend_gp = gpar(fill = c("green4", "red4")))

draw(lo_syn_heat, annotation_legend_list = biomarkers_legend, 
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
  merge_legend = TRUE)
```

<div class="figure">
<img src="index_files/figure-html/lo-syn-heat-1.png" alt="Heatmap of Average Link Operator Differences between synergistic and antagonistic model groups for each observed synergy" width="672" />
<p class="caption">(\#fig:lo-syn-heat)Heatmap of Average Link Operator Differences between synergistic and antagonistic model groups for each observed synergy</p>
</div>

:::{.green-box}
- `ERK_f` is a strong synergy OR-NOT biomarker for all synergies (same as in the *performance* link-operator biomarkers section and shown in this [figure](#fig:nodes-lo-maps-1))
- `MAPK14` is again the node that mostly distinguishes between the two synergy sub-clusters of `PD` and `5Z`. 
This also relates to the embedding of the `MAPK14` link-operator value as demonstrated by this [figure](#fig:nodes-lo-maps-1), where the `MAPK14` node has the `OR-NOT` link-operator in the `PD` sub-cluster and `AND-NOT` in the `5Z` sub-cluster.
:::

# R session info {-}


```r
xfun::session_info()
```

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

Locale:
  LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
  LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
  LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
  LC_PAPER=en_US.UTF-8       LC_NAME=C                 
  LC_ADDRESS=C               LC_TELEPHONE=C            
  LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

Package version:
  abind_1.4-5              assertthat_0.2.1         backports_1.1.10        
  base64enc_0.1.3          BH_1.72.0.3              bookdown_0.21           
  boot_1.3.25              broom_0.7.2              callr_3.5.1             
  car_3.0-10               carData_3.0-4            cellranger_1.1.0        
  circlize_0.4.10          Ckmeans.1d.dp_4.3.3      cli_2.1.0               
  clipr_0.7.1              clue_0.3-57              cluster_2.1.0           
  codetools_0.2-16         colorspace_1.4-1         compiler_3.6.3          
  ComplexHeatmap_2.2.0     conquer_1.0.2            corrplot_0.84           
  cowplot_1.1.0            cpp11_0.2.3              crayon_1.3.4            
  crosstalk_1.1.0.1        curl_4.3                 data.table_1.13.2       
  desc_1.2.0               digest_0.6.27            dplyr_1.0.2             
  dqrng_0.2.1              DT_0.16                  ellipsis_0.3.1          
  emba_0.1.8               evaluate_0.14            fansi_0.4.1             
  farver_2.0.3             FNN_1.1.3                forcats_0.5.0           
  foreach_1.5.1            foreign_0.8-75           gbRd_0.4-11             
  generics_0.0.2           GetoptLong_1.0.4         ggplot2_3.3.2           
  ggpubr_0.4.0             ggrepel_0.8.2            ggsci_2.9               
  ggsignif_0.6.0           glmnet_4.0-2             GlobalOptions_0.1.2     
  glue_1.4.2               graphics_3.6.3           grDevices_3.6.3         
  grid_3.6.3               gridExtra_2.3            gtable_0.3.0            
  haven_2.3.1              highr_0.8                hms_0.5.3               
  htmltools_0.5.0          htmlwidgets_1.5.2        igraph_1.2.6            
  irlba_2.3.3              isoband_0.2.2            iterators_1.0.13        
  jsonlite_1.7.1           knitr_1.30               labeling_0.4.2          
  later_1.1.0.1            lattice_0.20-41          lazyeval_0.2.2          
  lifecycle_0.2.0          lme4_1.1.25              magrittr_1.5            
  maptools_1.0.2           markdown_1.1             MASS_7.3.53             
  Matrix_1.2-18            MatrixModels_0.4.1       matrixStats_0.57.0      
  methods_3.6.3            mgcv_1.8.33              mime_0.9                
  minqa_1.2.4              munsell_0.5.0            nlme_3.1.149            
  nloptr_1.2.2.2           nnet_7.3.14              openxlsx_4.2.2          
  parallel_3.6.3           pbkrtest_0.4.8.6         pillar_1.4.6            
  pkgbuild_1.1.0           pkgconfig_2.0.3          pkgload_1.1.0           
  png_0.1-7                polynom_1.4.0            praise_1.0.0            
  prettyunits_1.1.1        processx_3.4.4           progress_1.2.2          
  promises_1.1.1           ps_1.4.0                 purrr_0.3.4             
  quantreg_5.74            R6_2.4.1                 randomForest_4.6-14     
  ranger_0.12.1            rbibutils_1.3            RColorBrewer_1.1-2      
  Rcpp_1.0.5               RcppAnnoy_0.0.16         RcppArmadillo_0.10.1.0.0
  RcppEigen_0.3.3.7.0      RcppProgress_0.4.2       Rdpack_2.0              
  readr_1.4.0              readxl_1.3.1             rematch_1.0.1           
  rio_0.5.16               rje_1.10.16              rjson_0.2.20            
  rlang_0.4.8              rmarkdown_2.5            rprojroot_1.3.2         
  RSpectra_0.16.0          rstatix_0.6.0            rstudioapi_0.11         
  scales_1.1.1             shape_1.4.5              sitmo_2.0.1             
  sp_1.4.4                 SparseM_1.78             splines_3.6.3           
  statmod_1.4.35           stats_3.6.3              stringi_1.5.3           
  stringr_1.4.0            survival_3.2-7           testthat_2.3.2          
  tibble_3.0.4             tidyr_1.1.2              tidyselect_1.1.0        
  tinytex_0.26             tools_3.6.3              usefun_0.4.8            
  utf8_1.1.4               utils_3.6.3              uwot_0.1.8              
  vctrs_0.3.4              viridisLite_0.3.0        visNetwork_2.0.9        
  withr_2.3.0              xfun_0.18                xml2_1.3.2              
  yaml_2.2.1               zip_2.1.1               
```

# References {-}
