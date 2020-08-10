---
title: "Balance Mutations in Logical Modeling"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 10 August, 2020"
description: "Investigations for Balance link operators paper"
url: 'https\://bblodfon.github.io/balance-paper/'
github-repo: "bblodfon/balance-paper"
bibliography: references.bib
link-citations: true
site: bookdown::bookdown_site
---

# Intro {-}

Several analyses/investigations relating to the balance logical operators paper.

# Input {-}

Loading libraries:

```r
library(xfun)
library(knitr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
```

# Analysis {-}

## Link Operators {-}

Let $f$ be a boolean function $f(x,y):\{0,1\}^n \rightarrow \{0,1\}$, with $m$ **activators** $x=\{x_i\}_{i=1}^{m}$ and $k$ **inhibitors** $y=\{y_i\}_{i=1}^{k}$, that is a total of $n=m+k$ regulators.
For a link operator to make sense, we have that $m,k \ge 1$ (at least one regulator in each category).
We next describe what each link operator looks like in the general form of $f$:

- `AND-NOT`: $$f(x,y) = \left(\bigvee_{i=1}^{m} x_i\right) \land \lnot \left(\bigvee_{i=1}^{k} y_i\right)$$
- `OR-NOT`: $$f(x,y) = \left(\bigvee_{i=1}^{m} x_i\right) \lor \lnot \left(\bigvee_{i=1}^{k} y_i\right)$$
- `BalanceOp1`: $$f(x,y) = \bigvee_{\forall (i,j)}^{m,k}(x_i\land \lnot y_j) = \left(\bigvee_{i=1}^{m} x_i\right) \land \left(\bigvee_{i=1}^{k} \lnot y_i\right)$$
- `exp_act_win`: $$f_{act-win}(x,y)=\begin{cases}
        1, & \text{for } \sum_{i=1}^{m} x_i \ge \sum_{i=1}^{k} y_i, \text{ with } \sum_{i=1}^{m} x_i \ne 0\\
        0, & \text{otherwise}
        \end{cases}$$
- `exp_inh_win`: $$f_{inh-win}(x,y)=\begin{cases}
        1, & \text{for } \sum_{i=1}^{m} x_i \gt \sum_{i=1}^{k} y_i, \text{ with } \sum_{i=1}^{m} x_i \ne 0\\
        0, & \text{otherwise}
        \end{cases}$$

Note that: $f_{inh-win}(x,y) = \lnot f_{act-win}(y,x)$.
I searched for an analytical formula for the two last functions (they get pertty big!).
See discussion on [math.stackexchange](https://math.stackexchange.com/questions/3767774/identify-boolean-function-that-satisfies-some-constrains/).

# Truth Density Data Analysis {-}

:::{.blue-box}
*Truth Density* of a boolean equation/expression, given it's equivalent truth table, is the **number of rows that the expression is active** divided to **the total number of rows** $(2^n)$.
:::

I created every possible truth table for up to $20$ variables (variables here means *regulators* for us) and calculated the `AND-NOT`, `OR-NOT`, `BalanceOp1`, `exp_act_win`, `exp_inh_win` results for every possible configuration of the number of activators and inhibitors that added up to the number of regulators.
Then, from the truth tables I calculated the **truth density** of each operator in each particular configuration.
See part of the data below:

```r
stats = readRDS(file = "stats.rds")
stats[1:5,1:6] %>% kable(caption = "Thuth Density Data", digits = 2)
```



Table: (\#tab:load-data)Thuth Density Data

| num_reg| num_act| num_inh| td_and_not| td_or_not| td_balance_op|
|-------:|-------:|-------:|----------:|---------:|-------------:|
|       2|       1|       1|       0.25|      0.75|          0.25|
|       3|       1|       2|       0.12|      0.62|          0.38|
|       3|       2|       1|       0.38|      0.88|          0.38|
|       4|       1|       3|       0.06|      0.56|          0.44|
|       4|       2|       2|       0.19|      0.81|          0.56|

:::{.orange-box}
Use the [fun.R](https://github.com/bblodfon/balance-paper/blob/master/fun.R) script to reproduce this data.
:::

Also, I have proved the exact formulas for the truth densities in the case of the `AND-NOT` and `OR-NOT` link operators (see [here](http://tiny.cc/link-proofs) for a proof sketch).
I write them here explicitly, as well as their long-term behaviour (for large $n$):

- `AND-NOT`: $$TD_{AND-NOT}=\frac{2^m-1}{2^n} \xrightarrow{n \text{ large}} \frac{1}{2^k}$$
- `OR-NOT`:  $$TD_{OR-NOT}=\frac{2^n-2^k}{2^n} \xrightarrow{n \text{ large}} 1-\frac{1}{2^m}$$

So, if there are a lot of regulators in a boolean equation and **only one of them is an inhibitor (activator)** the truth density of the `AND-NOT` (`OR-NOT`) equation case will converge to $1/2$.

We can use the data above to validate the formulas from the proof (up to $n=20$):

```r
# Validate AND-NOT Truth Density formula
formula_td_and_not = stats %>% 
  mutate(formula_td_and_not = (2^num_act - 1)/(2^num_reg)) %>%
  pull(formula_td_and_not)

all(stats %>% pull(td_and_not) == formula_td_and_not)
```

```
[1] TRUE
```

```r
# Validate OR-NOT Truth Density formula
formula_td_or_not = stats %>% 
  mutate(formula_td_or_not = (((2^num_act - 1) * (2^num_inh)) + 1)/(2^num_reg)) %>%
  pull(formula_td_or_not)

all(stats %>% pull(td_or_not) == formula_td_or_not)
```

```
[1] TRUE
```

Comparing the `AND-NOT` and `OR-NOT` truth densities across the number of regulators:

```r
# tidy up data
stats_and_or = pivot_longer(data = stats, cols = c(td_and_not, td_or_not), 
  names_to = "lo", values_to = "td") %>%
  select(num_reg, lo, td) %>%
  mutate(lo = replace(x = lo, list = lo == "td_and_not", values = "AND-NOT")) %>%
  mutate(lo = replace(x = lo, list = lo == "td_or_not", values = "OR-NOT")) %>%
  rename(`Link Operator` = lo)

ggboxplot(data = stats_and_or, x = "num_reg", y = "td", 
  color = "Link Operator", palette = "Set1",
  title = "AND-NOT vs OR-NOT Truth Densities", 
  xlab = "Number of regulators", ylab = "Truth Density") +
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/fig-and-or-not-1.png" alt="AND-NOT vs OR-NOT Truth Densities across all possible activators and inhibitors combinations up to 20 regulators" width="672" />
<p class="caption">(\#fig:fig-and-or-not)AND-NOT vs OR-NOT Truth Densities across all possible activators and inhibitors combinations up to 20 regulators</p>
</div>

:::{.green-box}
- **The more regulators** there are, the more likely it is that the `AND-NOT` link operator in the boolean equation will result in an **inhibited** target and that the `OR-NOT` link operator in an **active** target.
- For $n>6$, the points outside the boxplots (with a truth density of $\frac{1}{2}, \frac{1}{4}, 1-\frac{1}{4},\frac{1}{8},1-\frac{1}{8},...$) correspond to the **long-term behaviour** of the truth density formulas shown above, where there is **large imbalance between the number of activators and inhibitors**.
:::



```r
# ggboxplot(data = stats, x = "num_reg", y = "td_balance_op")
# ggboxplot(data = stats, x = "num_reg", y = "td_exp_act")
# ggboxplot(data = stats, x = "num_reg", y = "td_exp_inh")
```


# R session info {-}


```r
xfun::session_info()
```

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.5 LTS

Locale:
  LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
  LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
  LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
  LC_PAPER=en_US.UTF-8       LC_NAME=C                 
  LC_ADDRESS=C               LC_TELEPHONE=C            
  LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

Package version:
  abind_1.4-5         assertthat_0.2.1    backports_1.1.8    
  base64enc_0.1.3     BH_1.72.0.3         bookdown_0.20      
  boot_1.3.25         broom_0.5.6         callr_3.4.3        
  car_3.0-8           carData_3.0-4       cellranger_1.1.0   
  cli_2.0.2           clipr_0.7.0         colorspace_1.4-1   
  compiler_3.6.3      corrplot_0.84       cowplot_1.0.0      
  crayon_1.3.4        curl_4.3            data.table_1.12.8  
  desc_1.2.0          digest_0.6.25       dplyr_1.0.0        
  ellipsis_0.3.1      evaluate_0.14       fansi_0.4.1        
  farver_2.0.3        forcats_0.5.0       foreign_0.8-75     
  generics_0.0.2      ggplot2_3.3.2       ggpubr_0.4.0       
  ggrepel_0.8.2       ggsci_2.9           ggsignif_0.6.0     
  glue_1.4.1          graphics_3.6.3      grDevices_3.6.3    
  grid_3.6.3          gridExtra_2.3       gtable_0.3.0       
  haven_2.3.1         highr_0.8           hms_0.5.3          
  htmltools_0.5.0     isoband_0.2.2       jsonlite_1.7.0     
  knitr_1.29          labeling_0.3        lattice_0.20-41    
  lifecycle_0.2.0     lme4_1.1.23         magrittr_1.5       
  maptools_1.0.1      markdown_1.1        MASS_7.3.51.6      
  Matrix_1.2.18       MatrixModels_0.4.1  methods_3.6.3      
  mgcv_1.8.31         mime_0.9            minqa_1.2.4        
  munsell_0.5.0       nlme_3.1-148        nloptr_1.2.2.1     
  nnet_7.3.14         openxlsx_4.1.5      parallel_3.6.3     
  pbkrtest_0.4.8.6    pillar_1.4.4        pkgbuild_1.0.8     
  pkgconfig_2.0.3     pkgload_1.1.0       plyr_1.8.6         
  polynom_1.4.0       praise_1.0.0        prettyunits_1.1.1  
  processx_3.4.2      progress_1.2.2      ps_1.3.3           
  purrr_0.3.4         quantreg_5.55       R6_2.4.1           
  RColorBrewer_1.1.2  Rcpp_1.0.4.6        RcppEigen_0.3.3.7.0
  readr_1.3.1         readxl_1.3.1        rematch_1.0.1      
  reshape2_1.4.4      rio_0.5.16          rlang_0.4.6        
  rmarkdown_2.3       rprojroot_1.3.2     rstatix_0.6.0      
  rstudioapi_0.11     scales_1.1.1        sp_1.4.2           
  SparseM_1.78        splines_3.6.3       statmod_1.4.34     
  stats_3.6.3         stringi_1.4.6       stringr_1.4.0      
  testthat_2.3.2      tibble_3.0.1        tidyr_1.1.0        
  tidyselect_1.1.0    tinytex_0.24        tools_3.6.3        
  utf8_1.1.4          utils_3.6.3         vctrs_0.3.1        
  viridisLite_0.3.0   withr_2.2.0         xfun_0.15          
  yaml_2.2.1          zip_2.0.4          
```

# References {-}
