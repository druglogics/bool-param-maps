---
title: "Balance Mutations in Logical Modeling"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 27 August, 2020"
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
library(tibble)
library(corrplot)
library(latex2exp)
library(ggpubr)
library(ggplot2)
library(DT)
library(usefun)
library(emba)
library(forcats)
library(scales)
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
stats = readRDS(file = "data/stats.rds")

DT::datatable(data = stats,
  caption = htmltools::tags$caption("Truth Density Data", style="color:#dd4814; font-size: 18px"),
  options = list(pageLength = 6, scrollX = TRUE, order = list(list(1, "asc")))) %>% 
  formatRound(4:8, digits = 2)
```

<!--html_preserve--><div id="htmlwidget-00dda29ef670179b0f78" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-00dda29ef670179b0f78">{"x":{"filter":"none","caption":"<caption style=\"color:#dd4814; font-size: 18px\">Truth Density Data<\/caption>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190"],[2,3,3,4,4,4,5,5,5,5,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20],[1,1,2,1,2,3,1,2,3,4,1,2,3,4,5,1,2,3,4,5,6,1,2,3,4,5,6,7,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,8,9,10,1,2,3,4,5,6,7,8,9,10,11,1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12,13,1,2,3,4,5,6,7,8,9,10,11,12,13,14,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],[1,2,1,3,2,1,4,3,2,1,5,4,3,2,1,6,5,4,3,2,1,7,6,5,4,3,2,1,8,7,6,5,4,3,2,1,9,8,7,6,5,4,3,2,1,10,9,8,7,6,5,4,3,2,1,11,10,9,8,7,6,5,4,3,2,1,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1],[0.25,0.125,0.375,0.0625,0.1875,0.4375,0.03125,0.09375,0.21875,0.46875,0.015625,0.046875,0.109375,0.234375,0.484375,0.0078125,0.0234375,0.0546875,0.1171875,0.2421875,0.4921875,0.00390625,0.01171875,0.02734375,0.05859375,0.12109375,0.24609375,0.49609375,0.001953125,0.005859375,0.013671875,0.029296875,0.060546875,0.123046875,0.248046875,0.498046875,0.0009765625,0.0029296875,0.0068359375,0.0146484375,0.0302734375,0.0615234375,0.1240234375,0.2490234375,0.4990234375,0.00048828125,0.00146484375,0.00341796875,0.00732421875,0.01513671875,0.03076171875,0.06201171875,0.12451171875,0.24951171875,0.49951171875,0.000244140625,0.000732421875,0.001708984375,0.003662109375,0.007568359375,0.015380859375,0.031005859375,0.062255859375,0.124755859375,0.249755859375,0.499755859375,0.0001220703125,0.0003662109375,0.0008544921875,0.0018310546875,0.0037841796875,0.0076904296875,0.0155029296875,0.0311279296875,0.0623779296875,0.1248779296875,0.2498779296875,0.4998779296875,6.103515625e-05,0.00018310546875,0.00042724609375,0.00091552734375,0.00189208984375,0.00384521484375,0.00775146484375,0.01556396484375,0.03118896484375,0.06243896484375,0.12493896484375,0.24993896484375,0.49993896484375,3.0517578125e-05,9.1552734375e-05,0.000213623046875,0.000457763671875,0.000946044921875,0.001922607421875,0.003875732421875,0.007781982421875,0.015594482421875,0.031219482421875,0.062469482421875,0.124969482421875,0.249969482421875,0.499969482421875,1.52587890625e-05,4.57763671875e-05,0.0001068115234375,0.0002288818359375,0.0004730224609375,0.0009613037109375,0.0019378662109375,0.0038909912109375,0.0077972412109375,0.0156097412109375,0.0312347412109375,0.0624847412109375,0.124984741210938,0.249984741210938,0.499984741210938,7.62939453125e-06,2.288818359375e-05,5.340576171875e-05,0.00011444091796875,0.00023651123046875,0.00048065185546875,0.00096893310546875,0.00194549560546875,0.00389862060546875,0.00780487060546875,0.0156173706054688,0.0312423706054688,0.0624923706054688,0.124992370605469,0.249992370605469,0.499992370605469,3.814697265625e-06,1.1444091796875e-05,2.6702880859375e-05,5.7220458984375e-05,0.000118255615234375,0.000240325927734375,0.000484466552734375,0.000972747802734375,0.00194931030273438,0.00390243530273438,0.00780868530273438,0.0156211853027344,0.0312461853027344,0.0624961853027344,0.124996185302734,0.249996185302734,0.499996185302734,1.9073486328125e-06,5.7220458984375e-06,1.33514404296875e-05,2.86102294921875e-05,5.91278076171875e-05,0.000120162963867188,0.000242233276367188,0.000486373901367188,0.000974655151367188,0.00195121765136719,0.00390434265136719,0.00781059265136719,0.0156230926513672,0.0312480926513672,0.0624980926513672,0.124998092651367,0.249998092651367,0.499998092651367,9.5367431640625e-07,2.86102294921875e-06,6.67572021484375e-06,1.43051147460938e-05,2.95639038085938e-05,6.00814819335938e-05,0.000121116638183594,0.000243186950683594,0.000487327575683594,0.000975608825683594,0.00195217132568359,0.00390529632568359,0.00781154632568359,0.0156240463256836,0.0312490463256836,0.0624990463256836,0.124999046325684,0.249999046325684,0.499999046325684],[0.75,0.625,0.875,0.5625,0.8125,0.9375,0.53125,0.78125,0.90625,0.96875,0.515625,0.765625,0.890625,0.953125,0.984375,0.5078125,0.7578125,0.8828125,0.9453125,0.9765625,0.9921875,0.50390625,0.75390625,0.87890625,0.94140625,0.97265625,0.98828125,0.99609375,0.501953125,0.751953125,0.876953125,0.939453125,0.970703125,0.986328125,0.994140625,0.998046875,0.5009765625,0.7509765625,0.8759765625,0.9384765625,0.9697265625,0.9853515625,0.9931640625,0.9970703125,0.9990234375,0.50048828125,0.75048828125,0.87548828125,0.93798828125,0.96923828125,0.98486328125,0.99267578125,0.99658203125,0.99853515625,0.99951171875,0.500244140625,0.750244140625,0.875244140625,0.937744140625,0.968994140625,0.984619140625,0.992431640625,0.996337890625,0.998291015625,0.999267578125,0.999755859375,0.5001220703125,0.7501220703125,0.8751220703125,0.9376220703125,0.9688720703125,0.9844970703125,0.9923095703125,0.9962158203125,0.9981689453125,0.9991455078125,0.9996337890625,0.9998779296875,0.50006103515625,0.75006103515625,0.87506103515625,0.93756103515625,0.96881103515625,0.98443603515625,0.99224853515625,0.99615478515625,0.99810791015625,0.99908447265625,0.99957275390625,0.99981689453125,0.99993896484375,0.500030517578125,0.750030517578125,0.875030517578125,0.937530517578125,0.968780517578125,0.984405517578125,0.992218017578125,0.996124267578125,0.998077392578125,0.999053955078125,0.999542236328125,0.999786376953125,0.999908447265625,0.999969482421875,0.500015258789062,0.750015258789062,0.875015258789062,0.937515258789062,0.968765258789062,0.984390258789062,0.992202758789062,0.996109008789062,0.998062133789062,0.999038696289062,0.999526977539062,0.999771118164062,0.999893188476562,0.999954223632812,0.999984741210938,0.500007629394531,0.750007629394531,0.875007629394531,0.937507629394531,0.968757629394531,0.984382629394531,0.992195129394531,0.996101379394531,0.998054504394531,0.999031066894531,0.999519348144531,0.999763488769531,0.999885559082031,0.999946594238281,0.999977111816406,0.999992370605469,0.500003814697266,0.750003814697266,0.875003814697266,0.937503814697266,0.968753814697266,0.984378814697266,0.992191314697266,0.996097564697266,0.998050689697266,0.999027252197266,0.999515533447266,0.999759674072266,0.999881744384766,0.999942779541016,0.999973297119141,0.999988555908203,0.999996185302734,0.500001907348633,0.750001907348633,0.875001907348633,0.937501907348633,0.968751907348633,0.984376907348633,0.992189407348633,0.996095657348633,0.998048782348633,0.999025344848633,0.999513626098633,0.999757766723633,0.999879837036133,0.999940872192383,0.999971389770508,0.99998664855957,0.999994277954102,0.999998092651367,0.500000953674316,0.750000953674316,0.875000953674316,0.937500953674316,0.968750953674316,0.984375953674316,0.992188453674316,0.996094703674316,0.998047828674316,0.999024391174316,0.999512672424316,0.999756813049316,0.999878883361816,0.999939918518066,0.999970436096191,0.999985694885254,0.999993324279785,0.999997138977051,0.999999046325684],[0.25,0.375,0.375,0.4375,0.5625,0.4375,0.46875,0.65625,0.65625,0.46875,0.484375,0.703125,0.765625,0.703125,0.484375,0.4921875,0.7265625,0.8203125,0.8203125,0.7265625,0.4921875,0.49609375,0.73828125,0.84765625,0.87890625,0.84765625,0.73828125,0.49609375,0.498046875,0.744140625,0.861328125,0.908203125,0.908203125,0.861328125,0.744140625,0.498046875,0.4990234375,0.7470703125,0.8681640625,0.9228515625,0.9384765625,0.9228515625,0.8681640625,0.7470703125,0.4990234375,0.49951171875,0.74853515625,0.87158203125,0.93017578125,0.95361328125,0.95361328125,0.93017578125,0.87158203125,0.74853515625,0.49951171875,0.499755859375,0.749267578125,0.873291015625,0.933837890625,0.961181640625,0.968994140625,0.961181640625,0.933837890625,0.873291015625,0.749267578125,0.499755859375,0.4998779296875,0.7496337890625,0.8741455078125,0.9356689453125,0.9649658203125,0.9766845703125,0.9766845703125,0.9649658203125,0.9356689453125,0.8741455078125,0.7496337890625,0.4998779296875,0.49993896484375,0.74981689453125,0.87457275390625,0.93658447265625,0.96685791015625,0.98052978515625,0.98443603515625,0.98052978515625,0.96685791015625,0.93658447265625,0.87457275390625,0.74981689453125,0.49993896484375,0.499969482421875,0.749908447265625,0.874786376953125,0.937042236328125,0.967803955078125,0.982452392578125,0.988311767578125,0.988311767578125,0.982452392578125,0.967803955078125,0.937042236328125,0.874786376953125,0.749908447265625,0.499969482421875,0.499984741210938,0.749954223632812,0.874893188476562,0.937271118164062,0.968276977539062,0.983413696289062,0.990249633789062,0.992202758789062,0.990249633789062,0.983413696289062,0.968276977539062,0.937271118164062,0.874893188476562,0.749954223632812,0.499984741210938,0.499992370605469,0.749977111816406,0.874946594238281,0.937385559082031,0.968513488769531,0.983894348144531,0.991218566894531,0.994148254394531,0.994148254394531,0.991218566894531,0.983894348144531,0.968513488769531,0.937385559082031,0.874946594238281,0.749977111816406,0.499992370605469,0.499996185302734,0.749988555908203,0.874973297119141,0.937442779541016,0.968631744384766,0.984134674072266,0.991703033447266,0.995121002197266,0.996097564697266,0.995121002197266,0.991703033447266,0.984134674072266,0.968631744384766,0.937442779541016,0.874973297119141,0.749988555908203,0.499996185302734,0.499998092651367,0.749994277954102,0.87498664855957,0.937471389770508,0.968690872192383,0.984254837036133,0.991945266723633,0.995607376098633,0.997072219848633,0.997072219848633,0.995607376098633,0.991945266723633,0.984254837036133,0.968690872192383,0.937471389770508,0.87498664855957,0.749994277954102,0.499998092651367,0.499999046325684,0.749997138977051,0.874993324279785,0.937485694885254,0.968720436096191,0.984314918518066,0.992066383361816,0.995850563049316,0.997559547424316,0.998047828674316,0.997559547424316,0.995850563049316,0.992066383361816,0.984314918518066,0.968720436096191,0.937485694885254,0.874993324279785,0.749997138977051,0.499999046325684],[0.75,0.5,0.875,0.3125,0.6875,0.9375,0.1875,0.5,0.8125,0.96875,0.109375,0.34375,0.65625,0.890625,0.984375,0.0625,0.2265625,0.5,0.7734375,0.9375,0.9921875,0.03515625,0.14453125,0.36328125,0.63671875,0.85546875,0.96484375,0.99609375,0.01953125,0.08984375,0.25390625,0.5,0.74609375,0.91015625,0.98046875,0.998046875,0.0107421875,0.0546875,0.171875,0.376953125,0.623046875,0.828125,0.9453125,0.9892578125,0.9990234375,0.005859375,0.03271484375,0.11328125,0.2744140625,0.5,0.7255859375,0.88671875,0.96728515625,0.994140625,0.99951171875,0.003173828125,0.019287109375,0.072998046875,0.19384765625,0.38720703125,0.61279296875,0.80615234375,0.927001953125,0.980712890625,0.996826171875,0.999755859375,0.001708984375,0.01123046875,0.046142578125,0.1334228515625,0.29052734375,0.5,0.70947265625,0.8665771484375,0.953857421875,0.98876953125,0.998291015625,0.9998779296875,0.00091552734375,0.0064697265625,0.0286865234375,0.08978271484375,0.21197509765625,0.395263671875,0.604736328125,0.78802490234375,0.91021728515625,0.9713134765625,0.9935302734375,0.99908447265625,0.99993896484375,0.00048828125,0.003692626953125,0.017578125,0.059234619140625,0.15087890625,0.303619384765625,0.5,0.696380615234375,0.84912109375,0.940765380859375,0.982421875,0.996307373046875,0.99951171875,0.999969482421875,0.0002593994140625,0.0020904541015625,0.0106353759765625,0.0384063720703125,0.105056762695312,0.227249145507812,0.401809692382812,0.598190307617188,0.772750854492188,0.894943237304688,0.961593627929688,0.989364624023438,0.997909545898438,0.999740600585938,0.999984741210938,0.0001373291015625,0.0011749267578125,0.0063629150390625,0.0245208740234375,0.0717315673828125,0.166152954101562,0.314529418945312,0.5,0.685470581054688,0.833847045898438,0.928268432617188,0.975479125976562,0.993637084960938,0.998825073242188,0.999862670898438,0.999992370605469,7.2479248046875e-05,0.0006561279296875,0.0037689208984375,0.01544189453125,0.048126220703125,0.118942260742188,0.240341186523438,0.407264709472656,0.592735290527344,0.759658813476562,0.881057739257812,0.951873779296875,0.98455810546875,0.996231079101562,0.999343872070312,0.999927520751953,0.999996185302734,3.814697265625e-05,0.000364303588867188,0.0022125244140625,0.00960540771484375,0.0317840576171875,0.0835342407226562,0.179641723632812,0.323802947998047,0.5,0.676197052001953,0.820358276367188,0.916465759277344,0.968215942382812,0.990394592285156,0.997787475585938,0.999635696411133,0.999961853027344,0.999998092651367,2.00271606445312e-05,0.000201225280761719,0.00128841400146484,0.00590896606445312,0.0206947326660156,0.0576591491699219,0.131587982177734,0.25172233581543,0.411901473999023,0.588098526000977,0.74827766418457,0.868412017822266,0.942340850830078,0.979305267333984,0.994091033935547,0.998711585998535,0.999798774719238,0.999979972839355,0.999999046325684],[0.25,0.125,0.5,0.0625,0.3125,0.6875,0.03125,0.1875,0.5,0.8125,0.015625,0.109375,0.34375,0.65625,0.890625,0.0078125,0.0625,0.2265625,0.5,0.7734375,0.9375,0.00390625,0.03515625,0.14453125,0.36328125,0.63671875,0.85546875,0.96484375,0.001953125,0.01953125,0.08984375,0.25390625,0.5,0.74609375,0.91015625,0.98046875,0.0009765625,0.0107421875,0.0546875,0.171875,0.376953125,0.623046875,0.828125,0.9453125,0.9892578125,0.00048828125,0.005859375,0.03271484375,0.11328125,0.2744140625,0.5,0.7255859375,0.88671875,0.96728515625,0.994140625,0.000244140625,0.003173828125,0.019287109375,0.072998046875,0.19384765625,0.38720703125,0.61279296875,0.80615234375,0.927001953125,0.980712890625,0.996826171875,0.0001220703125,0.001708984375,0.01123046875,0.046142578125,0.1334228515625,0.29052734375,0.5,0.70947265625,0.8665771484375,0.953857421875,0.98876953125,0.998291015625,6.103515625e-05,0.00091552734375,0.0064697265625,0.0286865234375,0.08978271484375,0.21197509765625,0.395263671875,0.604736328125,0.78802490234375,0.91021728515625,0.9713134765625,0.9935302734375,0.99908447265625,3.0517578125e-05,0.00048828125,0.003692626953125,0.017578125,0.059234619140625,0.15087890625,0.303619384765625,0.5,0.696380615234375,0.84912109375,0.940765380859375,0.982421875,0.996307373046875,0.99951171875,1.52587890625e-05,0.0002593994140625,0.0020904541015625,0.0106353759765625,0.0384063720703125,0.105056762695312,0.227249145507812,0.401809692382812,0.598190307617188,0.772750854492188,0.894943237304688,0.961593627929688,0.989364624023438,0.997909545898438,0.999740600585938,7.62939453125e-06,0.0001373291015625,0.0011749267578125,0.0063629150390625,0.0245208740234375,0.0717315673828125,0.166152954101562,0.314529418945312,0.5,0.685470581054688,0.833847045898438,0.928268432617188,0.975479125976562,0.993637084960938,0.998825073242188,0.999862670898438,3.814697265625e-06,7.2479248046875e-05,0.0006561279296875,0.0037689208984375,0.01544189453125,0.048126220703125,0.118942260742188,0.240341186523438,0.407264709472656,0.592735290527344,0.759658813476562,0.881057739257812,0.951873779296875,0.98455810546875,0.996231079101562,0.999343872070312,0.999927520751953,1.9073486328125e-06,3.814697265625e-05,0.000364303588867188,0.0022125244140625,0.00960540771484375,0.0317840576171875,0.0835342407226562,0.179641723632812,0.323802947998047,0.5,0.676197052001953,0.820358276367188,0.916465759277344,0.968215942382812,0.990394592285156,0.997787475585938,0.999635696411133,0.999961853027344,9.5367431640625e-07,2.00271606445312e-05,0.000201225280761719,0.00128841400146484,0.00590896606445312,0.0206947326660156,0.0576591491699219,0.131587982177734,0.25172233581543,0.411901473999023,0.588098526000977,0.74827766418457,0.868412017822266,0.942340850830078,0.979305267333984,0.994091033935547,0.998711585998535,0.999798774719238,0.999979972839355]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>num_reg<\/th>\n      <th>num_act<\/th>\n      <th>num_inh<\/th>\n      <th>td_and_not<\/th>\n      <th>td_or_not<\/th>\n      <th>td_balance_op<\/th>\n      <th>td_exp_act<\/th>\n      <th>td_exp_inh<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":6,"scrollX":true,"order":[[1,"asc"]],"columnDefs":[{"targets":4,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 2, 3, \",\", \".\"); }"},{"targets":5,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 2, 3, \",\", \".\"); }"},{"targets":6,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 2, 3, \",\", \".\"); }"},{"targets":7,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 2, 3, \",\", \".\"); }"},{"targets":8,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 2, 3, \",\", \".\"); }"},{"className":"dt-right","targets":[1,2,3,4,5,6,7,8]},{"orderable":false,"targets":0}],"autoWidth":false,"orderClasses":false,"lengthMenu":[6,10,25,50,100]}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render"],"jsHooks":[]}</script><!--/html_preserve-->

:::{.orange-box}
Use the [get_stats.R](https://github.com/bblodfon/balance-paper/blob/master/scripts/get_stats.R) script to reproduce this data.
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

### *AND-NOT* vs *OR-NOT* TD {-}

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

caption.title = "Euclidean Distances between vectors of truth density values (Symmetric)"
DT::datatable(data = d %>% as.matrix(), options = list(dom = "t", scrollX = TRUE),
  caption = htmltools::tags$caption(caption.title, style="color:#dd4814; font-size: 18px")) %>% 
  formatRound(1:6, digits = 3) %>%
  formatStyle(columns = c("act_prop"), backgroundColor = styleInterval(breaks, col))
```

<!--html_preserve--><div id="htmlwidget-5de7f6b5609ad95629f4" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-5de7f6b5609ad95629f4">{"x":{"filter":"none","caption":"<caption style=\"color:#dd4814; font-size: 18px\">Euclidean Distances between vectors of truth density values (Symmetric)<\/caption>","data":[["act_prop","td_and_not","td_or_not","td_balance_op","td_exp_act","td_exp_inh"],[0,6.2357551189418,6.2357551189418,6.23071307451427,2.37324587183497,2.37324587183497],[6.2357551189418,0,11.5854262787016,10.842296589664,7.82178323744214,6.7696332482574],[6.2357551189418,11.5854262787016,0,2.49443825784938,6.7696332482574,7.82178323744214],[6.23071307451427,10.842296589664,2.49443825784938,0,7.12290730774761,7.92299559783636],[2.37324587183497,7.82178323744214,6.7696332482574,7.12290730774761,0,1.86374127113628],[2.37324587183497,6.7696332482574,7.82178323744214,7.92299559783636,1.86374127113628,0]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>act_prop<\/th>\n      <th>td_and_not<\/th>\n      <th>td_or_not<\/th>\n      <th>td_balance_op<\/th>\n      <th>td_exp_act<\/th>\n      <th>td_exp_inh<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","scrollX":true,"columnDefs":[{"targets":1,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":2,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":3,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":4,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":5,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":6,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"className":"dt-right","targets":[1,2,3,4,5,6]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false,"rowCallback":"function(row, data) {\nvar value=data[1]; $(this.api().cell(row, 1).node()).css({'background-color':isNaN(parseFloat(value)) ? '' : value <= 0.5933 ? \"rgb(255,255,255)\" : value <= 1.1866 ? \"rgb(255,244,244)\" : value <= 1.7799 ? \"rgb(255,232,232)\" : value <= 2.3732 ? \"rgb(255,221,221)\" : value <= 2.3732 ? \"rgb(255,210,210)\" : value <= 2.3732 ? \"rgb(255,198,198)\" : value <= 2.3732 ? \"rgb(255,187,187)\" : value <= 2.3732 ? \"rgb(255,176,176)\" : value <= 3.3376 ? \"rgb(255,164,164)\" : value <= 4.302 ? \"rgb(255,153,153)\" : value <= 5.2663 ? \"rgb(255,142,142)\" : value <= 6.2307 ? \"rgb(255,131,131)\" : value <= 6.232 ? \"rgb(255,119,119)\" : value <= 6.2332 ? \"rgb(255,108,108)\" : value <= 6.2345 ? \"rgb(255,97,97)\" : value <= 6.2358 ? \"rgb(255,85,85)\" : value <= 6.2358 ? \"rgb(255,74,74)\" : value <= 6.2358 ? \"rgb(255,63,63)\" : value <= 6.2358 ? \"rgb(255,51,51)\" : \"rgb(255,40,40)\"});\n}"}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render","options.columnDefs.2.render","options.columnDefs.3.render","options.columnDefs.4.render","options.columnDefs.5.render","options.rowCallback"],"jsHooks":[]}</script><!--/html_preserve-->


```r
plot(hclust(dist(d)), main = "Distance Dendogram of Thruth Densities",
  ylab = "Euclidean Distance", sub = "BBR Truth Densities", xlab = "")
```

<img src="index_files/figure-html/dist-dendogram-1.png" width="672" style="display: block; margin: auto;" />

:::{.green-box}
- The **threshold functions** have truth densities values that are **closer to the proportion of activators** for a varying number of regulators, compared to the `AND-NOT` and `OR-NOT` formulas.
As such they represent more realistic candidates for regulatory functions from a statistical point of view.
- The TD values of `OR-NOT` and `BalanceOp1` are in general very close (as we've also seen in previous Figure)
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

## Scale-free Properties {-}

:::{.blue-box}
In this section we demonstrate the **scale-free properties of the CASCADE 1.0 network**.
We show that both in- and out-degree distributions are asymptotically power-law.
:::

Use the script [get_distribution_stats.R](https://github.com/bblodfon/balance-paper/blob/master/scripts/get_distribution_stats.R) to generate the degree distribution stats.
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
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/in-degree-fig-1.png" alt="In Degree Distribution (CASCADE 1.0)" width="2100" />
<p class="caption">(\#fig:in-degree-fig)In Degree Distribution (CASCADE 1.0)</p>
</div>


```r
dd_stats %>% group_by(out_degree) %>% tally() %>%
  ggplot(aes(x = out_degree, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  geom_smooth(aes(color = "red"), se = FALSE, span = 0.58, show.legend = FALSE) + 
  theme_classic() +
  labs(title = "Out-Degree Distribution (CASCADE 1.0)", x = "Out Degree", y = "Number of Nodes")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/out-degree-fig-1.png" alt="Out Degree Distribution (CASCADE 1.0)" width="2100" />
<p class="caption">(\#fig:out-degree-fig)Out Degree Distribution (CASCADE 1.0)</p>
</div>

## Model Stable State Statistics {-}

Using [abmlog](https://github.com/druglogics/abmlog) we generated all $2^{23} = 8388608$ possible link operator mutated models for the CASCADE 1.0 topology.
The models are stored in both `.gitsbe` and `.bnet` files in the Zenodo dataset [TOADD]. 
The `gitsbe` files include also the fixpoint attractors.
Thus we can find the *frequency distribution* of the number of fixpoints across all produced models (use the script [count_model_ss.R](https://github.com/bblodfon/balance-paper/blob/master/scripts/count_model_ss.R)).
The model stable state (fixpoint) statistics are as follows:


```r
models_ss_stats = readRDS(file = "data/models_ss_stats.rds")

models_ss_stats %>% group_by(ss_num) %>% tally() %>%
  ggplot(aes(x = ss_num, y = n, fill = as.factor(ss_num))) +
  geom_bar(stat = "identity", show.legend = FALSE) + 
  geom_text(aes(label = n), vjust = -0.5) +
  theme_classic2() +
  labs(title = "Stable States Distribution", x = "Number of Stable States", y = "Number of models")
```

## Stable States Data {-}

:::{.note}
To load the stable state data for the models that have **1 stable state** use the Zenodo dataset (**TODO**) and the script [tidy_ss_data.R](https://github.com/bblodfon/balance-p
aper/blob/master/scripts/tidy_ss_data.R)
:::

## Parameterization and Stable State Agreement {-}

We calculate the `node_stats` tibble object using the [get_node_stats.R](https://github.com/bblodfon/balance-paper/blob/master/scripts/get_node_stats.R) script.
This object includes the agreement statistics information for each node that has a link operator (i.e. it is targeted by both activators and inhibitors).
Only the **models that had 1 stable state** are considered in this analysis section.

Load the `node_stats`:

```r
node_stats = readRDS(file = "data/node_stats.rds")
```

:::{.note}
We are interested in two variables of interest:

- **Parameterization** of a link operator node: `AND-NOT` (0) vs `OR-NOT` (1)
- **Stable State** of a node: *inhibited* (0) vs *active* (1)

There exist are 4 different possibilities related to 2 cases:

1. `0-0`, `1-1` => parameterization and stable state match (e.g. node was parameterized with `AND-NOT` and it's state was inhibited or it had `OR-NOT` and its state was active)
2. `1-0`, `0-1` => parameterization and stable state differ (e.g. node had `OR-NOT` and its state was inhibited, or `AND-NOT` and it's state was active)
:::

In the next Figure we show the **total observed proportionate agreement** for each node, which is the number of models for which parameterization and stable state matched (case 1 above) divided by the total amount of models:

```r
node_stats %>% mutate(node = forcats::fct_reorder(node, desc(num_reg))) %>% 
  ggplot(aes(x = node, y = obs_prop_agreement, fill = as.factor(num_reg))) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels=scales::percent) +
    labs(title = "Agreement between Link Operator Parameterization and Stable State Activity", x = "Target Nodes with both activating and inhibiting regulators", y = "Observed Proportionate Agreement") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_brewer(guide = guide_legend(reverse=TRUE, title = "#Regulators"), palette = "Set1")
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/ss-lo-agreement-prop-1.png" alt="Parameterization and Stable State activity agreement" width="2100" />
<p class="caption">(\#fig:ss-lo-agreement-prop)Parameterization and Stable State activity agreement</p>
</div>

:::{.green-box}
The total barplot area covered (i.e. the **total agreement score** so to speak) is **77.7294779%**.

The above score means that the is a higher probability than chance to assign a node the `AND-NOT` (resp. `OR-NOT`) link operator in its respective boolean equation and that node at the same time having an inhibited (resp. activated) stable state of 0 (.resp 1) in any CASCADE 1.0 link operator parameterized model.
This suggests that the corresponding boolean formula used **is biased** and the previous analysis in this report showed that for larger networks this behaviour will become statistically more prevalent.

As such, even though the number of regulators are **less than 6**, we find that there is **strong agreement** between *link operator* and *stable state activity* across all the nodes that have both types of regulators (activators and inhibitors).
This agreement can be seen stronger for some nodes than others.
:::

Next, we calculate per node, the proportion of link operator assignments that retained their expected (i.e. keeping the same digit) stable state activity (e.g. the proportion of models corresponding to the cases `0-0`/(`0-0` + `0-1`) for the `AND-NOT` link operator - similar for `OR-NOT`):

```r
node_stats %>% 
  mutate(and_not_0ss_prop = and_not_0ss_agreement/(and_not_0ss_agreement + and_not_1ss_disagreement)) %>% 
  mutate(or_not_1ss_prop  = or_not_1ss_agreement/(or_not_1ss_agreement + or_not_0ss_disagreement)) %>%
  select(node, num_reg, and_not_0ss_prop, or_not_1ss_prop, active_prop) %>%
  rename(`AND-NOT` = and_not_0ss_prop, `OR-NOT` = or_not_1ss_prop) %>%
  mutate(node = forcats::fct_reorder(node, desc(num_reg))) %>%
  pivot_longer(cols = c(`AND-NOT`, `OR-NOT`)) %>%
  ggplot(aes(x = node, y = value, fill = name)) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_y_continuous(labels=scales::percent) +
    labs(title = "Link Operator Parameterization Agreement with Stable State Activity", 
      x = "Target Nodes with both activating and inhibiting regulators", 
      y = "Observed Proportionate Agreement") +
    theme_classic() + theme(axis.text.x = element_text(angle = 90)) + 
    scale_fill_brewer(guide = guide_legend(title = "Link Operator"), palette = "Set1") + 
    geom_line(aes(y = active_prop, color = active_prop), group = 1, size = 1.2) +
    scale_color_gradient(labels=scales::percent, low="grey", high="green", 
      name = "%Models:active node", limits = c(0,1)) + 
    theme(legend.title = element_text(size = 10))
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/ss-comp-agreement-props-1.png" alt="Parameterization and Stable State activity agreement 2" width="2100" />
<p class="caption">(\#fig:ss-comp-agreement-props)Parameterization and Stable State activity agreement 2</p>
</div>

:::{.green-box}
- Higher proportional activity for a node correlates with higher `OR-NOT`-activated state agreement.
- `LRP_f` has 4 activators and 1 inhibitor and from the previous statistical analysis with found that $TD_{AND-NOT,4+1}=0.469$, $TD_{OR-NOT,4+1}=0.969$, numbers which correspond really well with the proportionate agreement scores found across all the CASCADE 1.0 models.
- `TSC_f` has 1 activator and 4 inhibitors (which corresponds well to it's total inhibition profile in all the models).
- `TSC_f` and `mTORC2_c` are always found inhibited and thus the agreement with the `AND-NOT`-inhibited state is perfect and the `OR-NOT`-activated state agreement zero.
:::

In the above Figure, wherever there is less than **0.5 disagreement**, we can always explain it with the *activity* proportion value and the number of activators being more (or less resp.) than the number of inhibitors - see following table:


```r
caption.title = "Link Operator Statistics"
DT::datatable(data = node_stats %>% select(node, num_reg, num_act, num_inh), 
  caption = htmltools::tags$caption(caption.title, style="color:#dd4814; font-size: 18px"),
  options = list(order = list(list(2, "desc")))) %>% 
  formatRound(5:6, digits = 3)
```

<!--html_preserve--><div id="htmlwidget-9669e94a6a8fa026d628" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9669e94a6a8fa026d628">{"x":{"filter":"none","caption":"<caption style=\"color:#dd4814; font-size: 18px\">Link Operator Statistics<\/caption>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23"],["mTORC2_c","JNK_f","MAPK14","RTPK_f","MEK_f","SHC1","PTEN","SOS1","ERK_f","RAF_f","mTORC1_c","GAB_f","PDPK1","IKBKB","TSC_f","TP53","MDM2","CYCS","CFLAR","LRP_f","CTNNB1","TCF7_f","DKK_g"],[2,3,3,4,3,2,2,2,2,4,3,2,2,2,5,2,3,2,2,5,2,2,2],[1,2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,2,1,1,4,1,1,1],[1,1,1,2,1,1,1,1,1,3,1,1,1,1,4,1,1,1,1,1,1,1,1]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>node<\/th>\n      <th>num_reg<\/th>\n      <th>num_act<\/th>\n      <th>num_inh<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"order":[[2,"desc"]],"columnDefs":[{"targets":5,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"targets":6,"render":"function(data, type, row, meta) { return DTWidget.formatRound(data, 3, 3, \",\", \".\"); }"},{"className":"dt-right","targets":[2,3,4]},{"orderable":false,"targets":0}],"autoWidth":false,"orderClasses":false}},"evals":["options.columnDefs.0.render","options.columnDefs.1.render"],"jsHooks":[]}</script><!--/html_preserve-->





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
  base64enc_0.1.3     BH_1.72.0.3         bibtex_0.4.2.2     
  bookdown_0.20       boot_1.3.25         broom_0.5.6        
  callr_3.4.3         car_3.0-8           carData_3.0-4      
  cellranger_1.1.0    Ckmeans.1d.dp_4.3.2 cli_2.0.2          
  clipr_0.7.0         colorspace_1.4-1    compiler_3.6.3     
  corrplot_0.84       cowplot_1.0.0       crayon_1.3.4       
  crosstalk_1.1.0.1   curl_4.3            data.table_1.12.8  
  desc_1.2.0          digest_0.6.25       dplyr_1.0.0        
  DT_0.14             ellipsis_0.3.1      emba_0.1.7         
  evaluate_0.14       fansi_0.4.1         farver_2.0.3       
  forcats_0.5.0       foreign_0.8-75      gbRd_0.4-11        
  generics_0.0.2      ggplot2_3.3.2       ggpubr_0.4.0       
  ggrepel_0.8.2       ggsci_2.9           ggsignif_0.6.0     
  glue_1.4.1          graphics_3.6.3      grDevices_3.6.3    
  grid_3.6.3          gridExtra_2.3       gtable_0.3.0       
  haven_2.3.1         highr_0.8           hms_0.5.3          
  htmltools_0.5.0     htmlwidgets_1.5.1   igraph_1.2.5       
  isoband_0.2.2       jsonlite_1.7.0      knitr_1.29         
  labeling_0.3        later_1.1.0.1       latex2exp_0.4.0    
  lattice_0.20-41     lazyeval_0.2.2      lifecycle_0.2.0    
  lme4_1.1.23         magrittr_1.5        maptools_1.0.1     
  markdown_1.1        MASS_7.3.51.6       Matrix_1.2.18      
  MatrixModels_0.4.1  methods_3.6.3       mgcv_1.8.31        
  mime_0.9            minqa_1.2.4         munsell_0.5.0      
  nlme_3.1-148        nloptr_1.2.2.1      nnet_7.3.14        
  openxlsx_4.1.5      parallel_3.6.3      pbkrtest_0.4.8.6   
  pillar_1.4.4        pkgbuild_1.0.8      pkgconfig_2.0.3    
  pkgload_1.1.0       plyr_1.8.6          polynom_1.4.0      
  praise_1.0.0        prettyunits_1.1.1   processx_3.4.2     
  progress_1.2.2      promises_1.1.1      ps_1.3.3           
  purrr_0.3.4         quantreg_5.55       R6_2.4.1           
  RColorBrewer_1.1.2  Rcpp_1.0.4.6        RcppEigen_0.3.3.7.0
  Rdpack_1.0.0        readr_1.3.1         readxl_1.3.1       
  rematch_1.0.1       reshape2_1.4.4      rio_0.5.16         
  rje_1.10.16         rlang_0.4.6         rmarkdown_2.3      
  rprojroot_1.3.2     rstatix_0.6.0       rstudioapi_0.11    
  scales_1.1.1        sp_1.4.2            SparseM_1.78       
  splines_3.6.3       statmod_1.4.34      stats_3.6.3        
  stringi_1.4.6       stringr_1.4.0       testthat_2.3.2     
  tibble_3.0.1        tidyr_1.1.0         tidyselect_1.1.0   
  tinytex_0.24        tools_3.6.3         usefun_0.4.8       
  utf8_1.1.4          utils_3.6.3         vctrs_0.3.1        
  viridisLite_0.3.0   visNetwork_2.0.9    withr_2.2.0        
  xfun_0.15           yaml_2.2.1          zip_2.0.4          
```

# References {-}
