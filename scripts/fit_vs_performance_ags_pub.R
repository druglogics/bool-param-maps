#################################################################
# Correlation between fitness to AGS publication model          #
# steady state (Flobak et al. 2015) and models MCC performance  #
#################################################################
library(readr)
library(emba)
library(usefun)
library(dplyr)
library(tibble)
library(Ckmeans.1d.dp)
library(ggplot2)
library(ggpubr)

# Read the models stable state data (see `get_ss_data.R` script)
ss_data = readRDS(file = "data/ss_data.rds")

# get steady state from publication (see `flobak2015_ss.ipynb`)
ss_pub = readr::read_csv(file = 'data/steadystate_ags_pub.csv')
ss_pub = ss_pub %>% unlist()

# Booleanize the 4 multivalued nodes
ss_pub["Prosurvival"] = 1 # (3) range is {0,1,2,3}, `Antisurvival` is 0 so no need to change
ss_pub["CCND1"] = 1 # (2) range is {0,1,2}, `Caspase37` is 0 so no need to change
stopifnot(all(ss_pub %in% c(0,1)))

# Update the node names
stopifnot(ncol(ss_data) == length(ss_pub)) # 77 nodes each

# should be 57 mismatched names, got MANUAL WORK TO DO!
sum(!names(ss_pub) %in% colnames(ss_data))

# See also S1 Table in Supplementary Material of Flobak et al. (2015)
names(ss_pub)[names(ss_pub) == "AKT"] = "AKT_f"
names(ss_pub)[names(ss_pub) == "ASK1"] = "MAP3K5"
names(ss_pub)[names(ss_pub) == "Axin"] = "AXIN1"
names(ss_pub)[names(ss_pub) == "betacatenin"] = "CTNNB1"
names(ss_pub)[names(ss_pub) == "betaTrCP"] = "BTRC"
names(ss_pub)[names(ss_pub) == "Caspase37"] = "CASP3"
names(ss_pub)[names(ss_pub) == "Caspase8"] = "CASP8"
names(ss_pub)[names(ss_pub) == "Caspase9"] = "CASP9"
names(ss_pub)[names(ss_pub) == "CK1"] = "CK1_f"
names(ss_pub)[names(ss_pub) == "cMYC"] = "MYC"
names(ss_pub)[names(ss_pub) == "CytochromeC"] = "CYCS"
names(ss_pub)[names(ss_pub) == "DKK1"] = "DKK_f"
names(ss_pub)[names(ss_pub) == "DKK1gene"] = "DKK_g"
names(ss_pub)[names(ss_pub) == "Dvl"] = "DVL_f"
names(ss_pub)[names(ss_pub) == "Egr1"] = "EGR1"
names(ss_pub)[names(ss_pub) == "ERK"] = "ERK_f"
names(ss_pub)[names(ss_pub) == "FOXO"] = "FOXO_f"
names(ss_pub)[names(ss_pub) == "Fz"] = "FZD_f"
names(ss_pub)[names(ss_pub) == "GAB"] = "GAB_f"
names(ss_pub)[names(ss_pub) == "GSK3"] = "GSK3_f"
names(ss_pub)[names(ss_pub) == "IKKA"] = "CHUK"
names(ss_pub)[names(ss_pub) == "IKKB"] = "IKBKB"
names(ss_pub)[names(ss_pub) == "JNK"] = "JNK_f"
names(ss_pub)[names(ss_pub) == "LRP"] = "LRP_f"
names(ss_pub)[names(ss_pub) == "MDM2gene"] = "MDM2_g"
names(ss_pub)[names(ss_pub) == "MEK"] = "MEK_f"
names(ss_pub)[names(ss_pub) == "MEKK4"] = "MAP3K4"
names(ss_pub)[names(ss_pub) == "MKK3"] = "MAP2K3"
names(ss_pub)[names(ss_pub) == "MKK4"] = "MAP2K4"
names(ss_pub)[names(ss_pub) == "MKK7"] = "MAP2K7"
names(ss_pub)[names(ss_pub) == "MLK3"] = "MAP3K11"
names(ss_pub)[names(ss_pub) == "MMP"] = "MMP_f"
names(ss_pub)[names(ss_pub) == "MSK"] = "MSK_f"
names(ss_pub)[names(ss_pub) == "mTORC1"] = "mTORC1_c"
names(ss_pub)[names(ss_pub) == "mTORC2"] = "mTORC2_c"
names(ss_pub)[names(ss_pub) == "NFkB"] = "NFKB_f"
names(ss_pub)[names(ss_pub) == "p38alpha"] = "MAPK14"
names(ss_pub)[names(ss_pub) == "p53"] = "TP53"
names(ss_pub)[names(ss_pub) == "PDK1"] = "PDPK1"
names(ss_pub)[names(ss_pub) == "PI3K"] = "PIK3CA"
names(ss_pub)[names(ss_pub) == "pras40"] = "AKT1S1"
names(ss_pub)[names(ss_pub) == "PTENgene"] = "PTEN_g"
names(ss_pub)[names(ss_pub) == "Rac"] = "RAC_f"
names(ss_pub)[names(ss_pub) == "Raf"] = "RAF_f"
names(ss_pub)[names(ss_pub) == "Ras"] = "KRAS"
names(ss_pub)[names(ss_pub) == "Rheb"] = "RHEB"
names(ss_pub)[names(ss_pub) == "RSK"] = "RSK_f"
names(ss_pub)[names(ss_pub) == "RTPK"] = "RTPK_f"
names(ss_pub)[names(ss_pub) == "RTPKgene"] = "RTPK_g"
names(ss_pub)[names(ss_pub) == "S6K"] = "S6K_f"
names(ss_pub)[names(ss_pub) == "SFRP1gene"] = "SFRP1_g"
names(ss_pub)[names(ss_pub) == "SHP2"] = "PTPN11"
names(ss_pub)[names(ss_pub) == "SOS"] = "SOS1"
names(ss_pub)[names(ss_pub) == "TAB"] = "TAB_f"
names(ss_pub)[names(ss_pub) == "TAK1"] = "MAP3K7"
names(ss_pub)[names(ss_pub) == "TCF"] = "TCF7_f"
names(ss_pub)[names(ss_pub) == "TSC"] = "TSC_f"

# now node names should match!
stopifnot(all(names(ss_pub) %in% colnames(ss_data)))
stopifnot(all(colnames(ss_data) %in% names(ss_pub)))

# read the AGS "Gold Standard/Literature curated" `steady_state` vector (see `get_1ss_fitness_data.R`)
steady_state = readRDS(file = "data/ags_steady_state.rds")

# compare the common nodes (24) => all equal!!!
stopifnot(all(ss_pub[names(steady_state)] == steady_state))

# randomize ss_pub (just for testing)
# ss_pub_new = sample(c(0,1), length(ss_pub), replace = TRUE)
# names(ss_pub_new) = names(ss_pub)
# ss_pub = ss_pub_new

# calculate models fitness to AGS publication steady state (77 nodes)
fit_data_pub = apply(ss_data[, names(ss_pub)], 1, usefun::get_percentage_of_matches, ss_pub)

# data check
stopifnot(all(names(fit_data_pub) == rownames(ss_data)))

# save result
saveRDS(object = fit_data_pub, file = "data/fit_data_pub.rds")

# Fitness density figure
fit_data_pub %>% as_tibble() %>%
  ggplot() +
  geom_density(aes(x = value), fill = "#5ab4ac", alpha = 0.7, adjust = 2, show.legend = FALSE) +
  xlab("Fitness to the AGS model's steady state in Flobak et al. (2015)") +
  ggtitle("Fitness Density  (All link-operator, 1 stable state models)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/1ss_models_fit_density_pub.png", dpi = "print", width = 7, height = 5)

# Fitness vs MCC correlation
mcc_res = readRDS(file = "data/emba_mcc_res/mcc_4_res.rds")
mcc_data = mcc_res$models.mcc
res = Ckmeans.1d.dp(x = mcc_data, k = 4) # 4 MCC classes

# check model name order
stopifnot(all(names(fit_data_pub) == names(mcc_data)))

# MCC vs fitness boxplot
my_comparisons = list(c(1,2), c(1,3), c(1,4), c(3,4))
bind_cols(fitness = fit_data_pub, MCC = mcc_data, mcc_class = as.factor(res$cluster)) %>%
  ggplot(aes(x = mcc_class, y = fitness, fill = mcc_class)) +
  geom_boxplot() +
  guides(fill = guide_legend(title = "MCC Class")) +
  ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
  labs(title = "MCC vs Fitness (All 1 stable state Models)", x = "MCC Class", y = "Fitness") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = "img/mcc_vs_fit_pub.png", dpi = "print", width = 7, height = 5)
