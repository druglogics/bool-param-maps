---
title: "Balance Mutations in Logical Modeling"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 24 August, 2020"
description: "Investigations for Balance link operators paper"
url: 'https\://bblodfon.github.io/balance-paper/'
github-repo: "bblodfon/balance-paper"
bibliography: references.bib
link-citations: true
site: bookdown::bookdown_site
---

# Intro {-}

Several analyses/investigations relating to the balance logical operators paper.

Loading libraries:

```r
library(xfun)
library(knitr)
library(dplyr)
library(tidyr)
library(corrplot)
library(latex2exp)
library(ggpubr)
library(DT)
```

# BBR Function Analysis {-}

## Balance Boolean Regulatory Functions (BBRs) {-}

:::{.green-box}
BBR => Balance Boolean Regulatory Functions
:::

Let $f$ be a boolean function $f(x,y):\{0,1\}^n \rightarrow \{0,1\}$, with $m$ **activators** $x=\{x_i\}_{i=1}^{m}$ and $k$ **inhibitors** $y=\{y_i\}_{i=1}^{k}$, that is a total of $n=m+k$ regulators.
A subset of BBRs have a non-DNF representation that puts the different category regulators in separate groups and a *link boolean operator* in between them.
As such, for a link operator to make sense, we have that $m,k \ge 1$ (at least one regulator in each category).
An example of such a function that has been used in the literature [@Mendoza2006] is the formula with the `AND-NOT` link operator: 

- `AND-NOT`: $$f(x,y) = \left(\bigvee_{i=1}^{m} x_i\right) \land \lnot \left(\bigvee_{i=1}^{k} y_i\right)$$

A **variant** of that one that shifts the balance in favor of the activators (as we will see the truth density significantly increases) is the function with the `OR-NOT` link operator:

- `OR-NOT`: $$f(x,y) = \left(\bigvee_{i=1}^{m} x_i\right) \lor \lnot \left(\bigvee_{i=1}^{k} y_i\right)$$

Another one of this type of functions is the next one:

- `BalanceOp1`: $$f(x,y) = \bigvee_{\forall (i,j)}^{m,k}(x_i\land \lnot y_j) = \left(\bigvee_{i=1}^{m} x_i\right) \land \left(\bigvee_{i=1}^{k} \lnot y_i\right)$$

Next, we introduce some **threshold functions** with *pseudo-Boolean constraints*:

- `exp_act_win`: $$f_{act-win}(x,y)=\begin{cases}
        1, & \text{for } \sum_{i=1}^{m} x_i \ge \sum_{i=1}^{k} y_i\\
        0, & \text{otherwise}
        \end{cases}$$
- `exp_inh_win`: $$f_{inh-win}(x,y)=\begin{cases}
        1, & \text{for } \sum_{i=1}^{m} x_i \gt \sum_{i=1}^{k} y_i\\
        0, & \text{otherwise}
        \end{cases}$$

Note that: $f_{inh-win}(x,y) = \lnot f_{act-win}(y,x)$.

I searched for an analytical formula for the two last functions (they get pretty big!).
More info and discussion about these 2 last formulas, see the [math.stackexchange question](https://math.stackexchange.com/questions/3767774/identify-boolean-function-that-satisfies-some-constrains/).

## Truth Density Data Analysis {-}

### Data {-}

:::{.blue-box}
*Truth Density (TD)* of a boolean equation/expression, given it's equivalent truth table, is the **number of rows that the expression is active** divided to **the total number of rows** $(2^n)$.
:::

I created every possible truth table for up to $20$ variables (variables here means *regulators* for us) and calculated the `AND-NOT`, `OR-NOT`, `BalanceOp1`, `exp_act_win`, `exp_inh_win` results for every possible configuration of the number of activators and inhibitors that added up to the number of regulators.
Then, from the truth tables I calculated the **truth density** of each operator in each particular configuration.
See part of the data below:

```r
stats = readRDS(file = "stats.rds")
stats[1:5,] %>% kable(caption = "Thuth Density Data", digits = 2)
```



Table: (\#tab:load-data)Thuth Density Data

| num_reg| num_act| num_inh| td_and_not| td_or_not| td_balance_op| td_exp_act| td_exp_inh|
|-------:|-------:|-------:|----------:|---------:|-------------:|----------:|----------:|
|       2|       1|       1|       0.25|      0.75|          0.25|       0.75|       0.25|
|       3|       1|       2|       0.12|      0.62|          0.38|       0.50|       0.12|
|       3|       2|       1|       0.38|      0.88|          0.38|       0.88|       0.50|
|       4|       1|       3|       0.06|      0.56|          0.44|       0.31|       0.06|
|       4|       2|       2|       0.19|      0.81|          0.56|       0.69|       0.31|

:::{.orange-box}
Use the [fun.R](https://github.com/bblodfon/balance-paper/blob/master/fun.R) script to reproduce this data.
:::

### Truth Density formulas {-}

Also, I have proved the exact formulas for the truth densities in the case of the `AND-NOT` and `OR-NOT` link operators (see [here](http://tiny.cc/link-proofs) for a proof sketch).
I write them here explicitly, as well as their long-term behaviour (for large $n$. number of regulators):

- `AND-NOT`: $$TD_{AND-NOT}=\frac{2^m-1}{2^n} \xrightarrow{n \text{ large}} \frac{1}{2^k}$$
- `OR-NOT`:  $$TD_{OR-NOT}=\frac{2^n-2^k}{2^n} \xrightarrow{n \text{ large}} 1-\frac{1}{2^m}$$

For large $n$, the $TD_{AND-NOT}$ depends only **on the number of inhibitors** while the $TD_{OR-NOT}$ depends only **on the number of activators**.

Also, again for large $n$, the extreme case of having a TD value equal to $0.5$ is a result of having **only one of the regulators being an inhibitor (activator)** of the `AND-NOT` (`OR-NOT`) equation.

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

### *AND-NOT* vs *OR-NOT* {-}

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
<img src="index_files/figure-html/fig-and-or-not-1.png" alt="AND-NOT vs OR-NOT Truth Densities across all possible activators and inhibitors combinations up to 20 regulators" width="2100" />
<p class="caption">(\#fig:fig-and-or-not)AND-NOT vs OR-NOT Truth Densities across all possible activators and inhibitors combinations up to 20 regulators</p>
</div>

:::{.green-box}
- **The more regulators** there are, the more likely it is that the `AND-NOT` link operator in the boolean equation will result in an **inhibited** target and that the `OR-NOT` link operator in an **active** target.
- For $n>6$, the points outside the boxplots (with a truth density of $\frac{1}{2}, \frac{1}{4}, 1-\frac{1}{4},\frac{1}{8},1-\frac{1}{8},...$) correspond to the **long-term behaviour** of the truth density formulas shown above, but where there is also **large imbalance between the number of activators and inhibitors**.
:::

We can also check the relation between TD and number of activators and inhibitors in each case.
The following two figures show us why **the number of inhibitors** are more decisive in the `AND-NOT` case:

```r
ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_inh", 
  y = "td_and_not", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Inhibitors", 
  title = "AND-NOT TD vs Number of Inhibitors") + 
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_act", 
  y = "td_and_not", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Activators", 
  title = "AND-NOT TD vs Number of Activators") + 
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure">
<img src="index_files/figure-html/and-not-reg-plot-1.png" alt="AND-NOT TD vs Number of Activators and Inhibitors" width="50%" /><img src="index_files/figure-html/and-not-reg-plot-2.png" alt="AND-NOT TD vs Number of Activators and Inhibitors" width="50%" />
<p class="caption">(\#fig:and-not-reg-plot)AND-NOT TD vs Number of Activators and Inhibitors</p>
</div>

In the `OR-NOT` case the number of activators is more important:

```r
ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_inh", 
  y = "td_or_not", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Inhibitors", 
  title = "OR-NOT TD vs Number of Inhibitors") + 
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_act", 
  y = "td_or_not", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Activators", 
  title = "OR-NOT TD vs Number of Activators") + 
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure">
<img src="index_files/figure-html/or-not-reg-plot-1.png" alt="OR-NOT TD vs Number of Activators and Inhibitors" width="50%" /><img src="index_files/figure-html/or-not-reg-plot-2.png" alt="OR-NOT TD vs Number of Activators and Inhibitors" width="50%" />
<p class="caption">(\#fig:or-not-reg-plot)OR-NOT TD vs Number of Activators and Inhibitors</p>
</div>

### *BalanceOp1* TD {-}

If we add the `BalanceOp1` formuls's TD results to the first plot we have:

```r
# tidy up data
stats_and_or_balance = pivot_longer(data = stats, cols = c(td_and_not, td_or_not, td_balance_op), 
  names_to = "lo", values_to = "td") %>%
  select(num_reg, lo, td) %>%
  mutate(lo = replace(x = lo, list = lo == "td_and_not", values = "AND-NOT")) %>%
  mutate(lo = replace(x = lo, list = lo == "td_or_not", values = "OR-NOT")) %>%
  mutate(lo = replace(x = lo, list = lo == "td_balance_op", values = "BalanceOp1")) %>%
  rename(`Link Operator` = lo)

ggboxplot(data = stats_and_or_balance, x = "num_reg", y = "td", 
  color = "Link Operator", palette = "Set1",
  title = "AND-NOT vs OR-NOT vs BalanceOp1 Truth Densities", 
  xlab = "Number of regulators", ylab = "Truth Density") +
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/fig-and-or-not-balanceop-1.png" alt="AND-NOT vs OR-NOT vs BalanceOp1 Truth Densities across all possible activators and inhibitors combinations up to 20 regulators" width="2100" />
<p class="caption">(\#fig:fig-and-or-not-balanceop)AND-NOT vs OR-NOT vs BalanceOp1 Truth Densities across all possible activators and inhibitors combinations up to 20 regulators</p>
</div>

:::{.green-box}
- The `BalanceOp1` TD values are closer to the TD values of the `OR-NOT` formula compared to the `AND-NOT` one.
- The `BalanceOp1` is less *biased* compared to the `OR-NOT` link operator, but still for large $n$ (regulators) it practically **makes the target activated**.
:::

As we can see in the following two figures, the `BalanceOp1` shows a more balanced dependency between the number of activators and inhibitors:

```r
ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_inh", 
  y = "td_balance_op", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Inhibitors", 
  title = "BalanceOp1 TD vs Number of Inhibitors") + 
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(data = stats %>% rename(`#Regulators` = num_reg), x = "num_act", 
  y = "td_balance_op", color = "#Regulators",
  ylab = "Truth Density", xlab = "Number of Activators", 
  title = "BalanceOp1 TD vs Number of Activators") + 
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure">
<img src="index_files/figure-html/balanceOp1-reg-plot-1.png" alt="BalanceOp1 TD vs Number of Activators and Inhibitors" width="50%" /><img src="index_files/figure-html/balanceOp1-reg-plot-2.png" alt="BalanceOp1 TD vs Number of Activators and Inhibitors" width="50%" />
<p class="caption">(\#fig:balanceOp1-reg-plot)BalanceOp1 TD vs Number of Activators and Inhibitors</p>
</div>

### Threshold Functions TD {-}

In contrast, if we check the truth density of the $f_{act-win}(x,y)$ and $f_{inh-win}(x,y)$ boolean functions we have:

```r
# tidy up data
stats_functions = pivot_longer(data = stats, cols = c(td_exp_act, td_exp_inh), 
  names_to = "fun", values_to = "td") %>%
  select(num_reg, fun, td) %>%
  mutate(fun = replace(x = fun, list = fun == "td_exp_act", values = "Activators Win")) %>%
  mutate(fun = replace(x = fun, list = fun == "td_exp_inh", values = "Inhibitors Win")) %>%
  rename(`Equation Formula` = fun)

ggboxplot(data = stats_functions, x = "num_reg", y = "td",
  color = "Equation Formula", palette = "lancet",
  title = TeX("Truth Densities of $f_{act-win}(x,y)$ and $f_{inh-win}(x,y)$"),
  xlab = "Number of regulators", ylab = "Truth Density") +
  theme(plot.title = element_text(hjust = 0.5))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/fig-two-bool-formulas-1.png" alt="Truth Desities of two robust boolean formulas across all possible activators and inhibitors combinations up to 20 regulators" width="2100" />
<p class="caption">(\#fig:fig-two-bool-formulas)Truth Desities of two robust boolean formulas across all possible activators and inhibitors combinations up to 20 regulators</p>
</div>

:::{.green-box}
- Both boolean functions have a **large variance of truth densities** irrespective of the number of regulators.
- The median values seem to converge to $0.5$ for both formulas.
- The median value of truth density for the $f_{act-win}(x,y)$ is always larger than the $f_{inh-win}(x,y)$ (as expected).
:::

### TD Data Distance {-}

We check how close are the truth density values of the different proposed BBRs, also compared to the **proportion of activators**, e.g. if a BBR has 1 activator and 5 inhibitors (resp. 5 activators and 1 inhibitor) I would expect my regulatory function's output to be statistically more inhibited (resp. activated).
We find the *euclidean distance* between the different truth density values and show them in a table and dendrogram format:


```r
act_prop = stats %>% mutate(act_prop = num_act/num_reg) %>% pull(act_prop)
td_and_not = stats %>% pull(td_and_not)
td_or_not = stats %>% pull(td_or_not)
td_balance_op = stats %>% pull(td_balance_op)
td_exp_act = stats %>% pull(td_exp_act)
td_exp_inh = stats %>% pull(td_exp_inh)

d = dist(rbind(act_prop, td_and_not, td_or_not, td_balance_op, td_exp_act, td_exp_inh))
```


```r
# color `act_prop` column
breaks = quantile(unname(as.matrix(d)[, "act_prop"]), probs = seq(.05, .95, .05), na.rm = TRUE)
col = round(seq(255, 40, length.out = length(breaks) + 1), 0) %>%
  {paste0("rgb(255,", ., ",", ., ")")} # red

caption.title = "Table 1: Euclidean Distances between vectors of truth density values (Symmetric)"
DT::datatable(data = d %>% as.matrix(), options = list(dom = "t", scrollX = TRUE),
  caption = htmltools::tags$caption(caption.title, style="color:#dd4814; font-size: 18px")) %>% 
  formatRound(1:6, digits = 3) %>%
  formatStyle(columns = c("act_prop"), backgroundColor = styleInterval(breaks, col))
```

<!--html_preserve--><div id="htmlwidget-c2a931752ae4f30d7598" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-c2a931752ae4f30d7598">{"x":{"filter":"none","caption":"<caption style=\"color:#dd4814; font-size: 18px\">Table 1: Euclidean Distances between vectors of truth density values (Symmetric)<\/caption>","data":[["act_prop","td_and_not","td_or_not","td_balance_op","td_exp_act","td_exp_inh"],[0,6.2357551189418,6.2357551189418,6.23071307451427,2.37324587183497,2.37324587183497],[6.2357551189418,0,11.5854262787016,10.842296589664,7.82178323744214,6.7696332482574],[6.2357551189418,11.5854262787016,0,2.49443825784938,6.7696332482574,7.82178323744214],[6.23071307451427,10.842296589664,2.49443825784938,0,7.12290730774761,7.92299559783636],[2.37324587183497,7.82178323744214,6.7696332482574,7.12290730774761,0,1.86374127113628],[2.37324587183497,6.7696332482574,7.82178323744214,7.92299559783636,1.86374127113628,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>act_prop<\/th>\n      <th>td_and_not<\/th>\n      <th>td_or_not<\/th>\n      <th>td_balance_op<\/th>\n      <th>td_exp_act<\/th>\n      <th>td_exp_inh<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","scrollX":true,"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":2,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":3,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":4,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":5,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":6,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"className":"dt-right","targets":[1,2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"rowCallback":"function(row, data) {\nvar value=data[1]; $(this.api().cell(row, 1).node()).css({'background-color':isNaN(parseFloat(value)) ? '' : value <= 0.5933 ? \"rgb(255,255,255)\" : value <= 1.1866 ? \"rgb(255,244,244)\" : value <= 1.7799 ? \"rgb(255,232,232)\" : value <= 2.3732 ? \"rgb(255,221,221)\" : value <= 2.3732 ? \"rgb(255,210,210)\" : value <= 2.3732 ? \"rgb(255,198,198)\" : value <= 2.3732 ? \"rgb(255,187,187)\" : value <= 2.3732 ? \"rgb(255,176,176)\" : value <= 3.3376 ? \"rgb(255,164,164)\" : value <= 4.302 ? \"rgb(255,153,153)\" : value <= 5.2663 ? \"rgb(255,142,142)\" : value <= 6.2307 ? \"rgb(255,131,131)\" : value <= 6.232 ? \"rgb(255,119,119)\" : value <= 6.2332 ? \"rgb(255,108,108)\" : value <= 6.2345 ? \"rgb(255,97,97)\" : value <= 6.2358 ? \"rgb(255,85,85)\" : value <= 6.2358 ? \"rgb(255,74,74)\" : value <= 6.2358 ? \"rgb(255,63,63)\" : value <= 6.2358 ? \"rgb(255,51,51)\" : \"rgb(255,40,40)\"});\n}"}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render","options.columnDefs.5.render","options.rowCallback"],"jsHooks":[]}</script><!--/html_preserve-->


```r
plot(hclust(dist(d)), main = "Distance Dendogram of Thruth Densities",
  ylab = "Euclidean Distance", sub = "BBR Truth Densities", xlab = "")
```

<img src="index_files/figure-html/dist-dendogram-1.png" width="672" style="display: block; margin: auto;" />

:::{.green-box}
- The **threshold functions** have truth densities values that are **closer to the proportion of activators** for a varying number of regulators, compared to the `AND-NOT` and `OR-NOT` formulas.
As such they represent more realistic candidates for regulatory functions from a statistical point of view.
- The TD values of `OR-NOT` and `BalanceOp1` are in general very close (as we've also seen in previous Figure)
- 
:::

### Correlation {-}

We will now check the *correlation* between each pair of operators/proposed functions, as well as the number of regulators, inhibitors and activators:

```r
M = cor(stats, method = "kendall")
res = cor.mtest(stats, method = "kendall")
corrplot(corr = M, type = "upper", p.mat = res$p, sig.level = c(.001, .01, .05), 
  pch.cex = 1, pch.col = "white", insig = "label_sig", tl.col = "black", tl.srt = 45)
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/cor-plot-1.png" alt="Correlation Matrix of Truth Densities and number of regulators" width="2100" />
<p class="caption">(\#fig:cor-plot)Correlation Matrix of Truth Densities and number of regulators</p>
</div>

:::{.green-box}
- The two functions results $f_{act-win}(x,y), f_{inh-win}(x,y)$ are highly correlated as expected
- Lower `AND-NOT` TD values highly correlate with *higher* number of inhibitors
- Higher `OR-NOT` TD values highly correlate with *higher* number of activators
:::

# CASCADE 1.0 Analysis {-}





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
  crayon_1.3.4        crosstalk_1.1.0.1   curl_4.3           
  data.table_1.12.8   desc_1.2.0          digest_0.6.25      
  dplyr_1.0.0         DT_0.14             ellipsis_0.3.1     
  evaluate_0.14       fansi_0.4.1         farver_2.0.3       
  forcats_0.5.0       foreign_0.8-75      generics_0.0.2     
  ggplot2_3.3.2       ggpubr_0.4.0        ggrepel_0.8.2      
  ggsci_2.9           ggsignif_0.6.0      glue_1.4.1         
  graphics_3.6.3      grDevices_3.6.3     grid_3.6.3         
  gridExtra_2.3       gtable_0.3.0        haven_2.3.1        
  highr_0.8           hms_0.5.3           htmltools_0.5.0    
  htmlwidgets_1.5.1   isoband_0.2.2       jsonlite_1.7.0     
  knitr_1.29          labeling_0.3        later_1.1.0.1      
  latex2exp_0.4.0     lattice_0.20-41     lazyeval_0.2.2     
  lifecycle_0.2.0     lme4_1.1.23         magrittr_1.5       
  maptools_1.0.1      markdown_1.1        MASS_7.3.51.6      
  Matrix_1.2.18       MatrixModels_0.4.1  methods_3.6.3      
  mgcv_1.8.31         mime_0.9            minqa_1.2.4        
  munsell_0.5.0       nlme_3.1-148        nloptr_1.2.2.1     
  nnet_7.3.14         openxlsx_4.1.5      parallel_3.6.3     
  pbkrtest_0.4.8.6    pillar_1.4.4        pkgbuild_1.0.8     
  pkgconfig_2.0.3     pkgload_1.1.0       plyr_1.8.6         
  polynom_1.4.0       praise_1.0.0        prettyunits_1.1.1  
  processx_3.4.2      progress_1.2.2      promises_1.1.1     
  ps_1.3.3            purrr_0.3.4         quantreg_5.55      
  R6_2.4.1            RColorBrewer_1.1.2  Rcpp_1.0.4.6       
  RcppEigen_0.3.3.7.0 readr_1.3.1         readxl_1.3.1       
  rematch_1.0.1       reshape2_1.4.4      rio_0.5.16         
  rlang_0.4.6         rmarkdown_2.3       rprojroot_1.3.2    
  rstatix_0.6.0       rstudioapi_0.11     scales_1.1.1       
  sp_1.4.2            SparseM_1.78        splines_3.6.3      
  statmod_1.4.34      stats_3.6.3         stringi_1.4.6      
  stringr_1.4.0       testthat_2.3.2      tibble_3.0.1       
  tidyr_1.1.0         tidyselect_1.1.0    tinytex_0.24       
  tools_3.6.3         utf8_1.1.4          utils_3.6.3        
  vctrs_0.3.1         viridisLite_0.3.0   withr_2.2.0        
  xfun_0.15           yaml_2.2.1          zip_2.0.4          
```

# References {-}
