##############################
# Project: Profiles of Child Mental Health Problems During the COVID-19 Pandemic in Jordan and Associations with Mothers’ Psychosocial Functioning
# Date: 15 April 2024
# Purpose: Code for analyses of the main manuscript
#############################
# Preamble
##########
options(scipen=999)
options(RGL_USE_NULL=TRUE)
packages <- c('dplyr', 'tidyr', 'ggplot2', 'kableExtra', 'gridExtra', 'reshape2', 'data.table',
              'cluster', 'factoextra', 'readr', 'kml','fBasics', 'nnet', 'mosaic', 'LongCART', 
              'rstatix', 'psych', 'qwraps2', 'gtsummary', 'forcats')
lapply(packages, require, character.only = TRUE)

#####################
# Data preparation
#####################
# load unidentified data
########################
dat <- read_csv('https://raw.githubusercontent.com/rali314/Child-Mental-Health-Problems-COVID-19-Pandemic-JO/main/publication_data_deid.csv',  na = c("", "NA", "."))

# Rename variables to make them more intuitive, create indices
df <- dat  %>% 
  dplyr::filter(!is.na(SDQ_Int_Base), !is.na(SDQ_Ext_Base)) %>% 
  mutate(LgHCC_W3_C = as.numeric(NA), LgHCC_W3_M = as.numeric(NA),
         GSE_W3 = as.numeric(NA), GSE_W4 = as.numeric(NA),
         PHQ_W1 = as.numeric(NA), PHQ_W2 = as.numeric(NA),
         ParStr_W1 = as.numeric(NA), ParStr_W2 = as.numeric(NA),
         RelatSp_W1 = as.numeric(NA), RelatSp_W2 = as.numeric(NA),
         LgHCC_W4_M = log10(HCC_W4_M),
         LgHCC_W4_C = log10(HCC_W4_C),
         d_CRPR_W1 = abs(CRPR_C_Bas - CRPR_W_Bas),
         d_CRPR_W2 = abs(CRPR_C_Jun20 - CRPR_W_Jun20),
         d_CRPR_W3 = abs(CRPR_C_Dec20 - CRPR_W_Dec20),
         d_CRPR_W4 = abs(CRPR_C_Jun21 - CRPR_W_Jun21),
         c_male = Gender,
         c_age = ChAge_Wave1,
         c_fborn = ifelse(C_Firstborn == 1, 1, 0),
         m_nat_jo = ifelse(M_Nat == 1, 1, 0),
         m_working = ifelse(M_Job == 1 | M_Job == 2, 1, 0),
         hh_size = Children + Adults,
         m_edu_uni = ifelse(M_School > 3, 1, 0),
         ) %>% .[-32,]

# Create vectors of relevant variable names to facilitate the analysis
control_vars <- c('c_male', 'c_age', 'c_fborn', 'm_nat_jo', 
                  'm_working', 'Income',
                  'm_edu_uni') #, 'M_jobloss_cov2', 'M_jobloss_cov3')
control_vars_plus <- c(control_vars, 'PHQ_4', 'ParStr_4')
control_vars2 <- c('c_male', 'c_age', 'c_fborn', 'm_nat_jo',
                  'm_edu_uni')
# Check outliers
##################
apply(df %>% dplyr::select(SDQ_Int_Base, SDQ_Int_Jun, SDQ_Int_Dec, SDQ_Int_Jun21, 
                           SDQ_Ext_Base, SDQ_Ext_Jun, SDQ_Ext_Dec, SDQ_Ext_Jun21, 
                           LgHCC_W1_C, LgHCC_W2_C, LgHCC_W3_C, LgHCC_W4_C, 
                           LgHCC_W1_M, LgHCC_W2_M, LgHCC_W3_M, LgHCC_W4_M, control_vars,
                           GSE_Mn, GSE_Mn_cov1, GSE_W3, GSE_W4, 
                           PHQ_W1, PHQ_W2, PHQ_Dec20, PHQ_Jun21,
                           ParStr_W1, ParStr_W2, ParStr_Dec20, ParStr_Jun21, 
                           RelatSp_W1, RelatSp_W2, RelatSp_Dec20, RelatSp_Jun21,
                           CRPR_W_Dec20, CRPR_C_Dec20, CRPR_W_Jun20, CRPR_C_Jun20, CRPR_W_Bas, CRPR_C_Bas, CRPR_W_Jun21, CRPR_C_Jun21,
                           d_CRPR_W1, d_CRPR_W2, d_CRPR_W3, d_CRPR_W4,
                           SDQ_TotDiff, SDQ_TotDiff_cov1, SDQ_TotDiff_cov2, SDQ_TotDiff_cov3), 
      MARGIN = 2, function(x){which(x < (mean(x) - 3*sd(x)) | x > (mean(x) + 3*sd(x)))})

# Fix outliers by replacing them with average values
df[14,'SDQ_Int_Base'] <- mean(df$SDQ_Int_Base)
df[16,'CRPR_W_Bas'] <- mean(df$CRPR_W_Bas)

# Descriptive stats of demographics
###################
demog_stats <- fBasics::basicStats(df %>%.[, control_vars] )
demog_stats_tbl <- as.data.frame(t(demog_stats)) %>%
  mutate(Observations = nobs - NAs,
         NA_p = NAs*100/nobs) %>%
  dplyr::select("NA_p", "Observations", "Mean", "Stdev", "Minimum", "Median", "Maximum") %>%
  kable("html", digits = 2, caption = "Demographic characteristics") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE)
demog_stats_tbl

sdq_w_temp <- df %>% dplyr::select(SDQ_Int_Base, SDQ_Int_Jun, SDQ_Int_Dec, SDQ_Int_Jun21, 
                SDQ_Ext_Base, SDQ_Ext_Jun, SDQ_Ext_Dec, SDQ_Ext_Jun21, 
                LgHCC_W1_C, LgHCC_W2_C, LgHCC_W3_C, LgHCC_W4_C, 
                LgHCC_W1_M, LgHCC_W2_M, LgHCC_W3_M, LgHCC_W4_M,
                # GSE_Mn, GSE_Mn_cov1, GSE_W3, GSE_W4, 
                PHQ_W1, PHQ_W2, PHQ_Dec20, PHQ_Jun21,
                ParStr_W1, ParStr_W2, ParStr_Dec20, ParStr_Jun21, 
                # RelatSp_W1, RelatSp_W2, RelatSp_Dec20, RelatSp_Jun21, 
                CRPR_W_Bas, CRPR_W_Jun20, CRPR_W_Dec20, CRPR_W_Jun21, CRPR_C_Bas, CRPR_C_Jun20, CRPR_C_Dec20, CRPR_C_Jun21,
                d_CRPR_W1, d_CRPR_W2, d_CRPR_W3, d_CRPR_W4,
                SDQ_TotDiff, SDQ_TotDiff_cov1, SDQ_TotDiff_cov2, SDQ_TotDiff_cov3,
                ID, control_vars,
                ) 

# Check for missing values
##########################
pMiss <- function(x){sum(is.na(x))/length(x)*100} # see % of missing values in each COL
na_tbl <- apply(sdq_w_temp,2,pMiss) %>% kable('html', digits = 2, caption = '% missing values') %>% 
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"), full_width=FALSE) 
na_tbl

# Reshape data to long format
###############################
sdq_l_temp <-  melt(setDT(sdq_w_temp), id.vars	= c('ID', 'c_male', 'c_age', 'c_fborn', 'm_nat_jo', 
                                                   'm_working', 'Income', 'm_edu_uni'), 
                            measure.vars = list(c(1:4), c(5:8), c(9:12), c(13:16), c(17:20), c(21:24),
                                                c(25:28), c(29:32),  c(33:36), c(37:40)),
                            value.name = c('SDQ_Int', 'SDQ_Ext', 'HCC_C', 'HCC_M', #'GSE', 
                                           'PHQ', 'ParStr', #'RelatSp',
                                           'CRPR_W', 'CRPR_C', 'd_CRPR', 'SDQ_TotDiff'),
                            variable.name = 'wave') %>%
  mutate(wave_int = as.integer(wave)) %>% 
  group_by(ID) %>%
  arrange(wave_int, .by_group = TRUE)


# Imputing missing values using LOCF method in the long data
sdq_l_temp[,c('SDQ_Int', 'SDQ_Ext')] <- nafill(sdq_l_temp[,c('SDQ_Int', 'SDQ_Ext')], 'locf')

# scaling variables
sdq_l <- sdq_l_temp %>% 
  mutate(s_SDQ_Int = as.numeric(scale(SDQ_Int, scale = TRUE)), # Adding standardized variables
         s_SDQ_Ext = as.numeric(scale(SDQ_Ext, scale = TRUE)))
summary(sdq_l)

# Reshape to wide again as will be used in the analysis
###############################
sdq_w <- sdq_l[,-20] %>% ungroup() %>% # had to drop wave_int here as it is redundant with wave
  tidyr::pivot_wider(names_from = wave, 
              values_from = c('SDQ_Int', 'SDQ_Ext', 's_SDQ_Int', 's_SDQ_Ext', 's_SDQ_Int', 
                              's_SDQ_Ext','HCC_C', 'HCC_M', #'GSE',
                              'PHQ', 'ParStr', #'RelatSp', 
                              'CRPR_W', 'CRPR_C', 'd_CRPR', 'SDQ_TotDiff'))

#########################
# Analysis: Longitudinal Clustering
#########################
# KML: Longitudinal k-means clustering 
#########################
# create clusterLongData objects
int_cld <- kml::clusterLongData(as.matrix(sdq_w[, c("SDQ_Int_1", "SDQ_Int_2", "SDQ_Int_3", "SDQ_Int_4")]), timeInData = 1:4)
ext_cld <- kml::clusterLongData(as.matrix(sdq_w[, c("SDQ_Ext_1", "SDQ_Ext_2", "SDQ_Ext_3", "SDQ_Ext_4")]), timeInData = 1:4)
sdq_cld <- kml::clusterLongData(as.matrix(sdq_w[, c("SDQ_TotDiff_1", "SDQ_TotDiff_2", "SDQ_TotDiff_3", "SDQ_TotDiff_4")]), timeInData = 1:4)

# Perform clustering for: 
# internalized behavioral problems, 
# externalized behavioral problems 
# and SDQ
set.seed(1)
kml::kml(int_cld, 1:6)
kml::kml(ext_cld, 1:6)
kml::kml(sdq_cld, 3:6)

# Plot and activate cluster selection chart
kml::choice(int_cld)
kml::choice(ext_cld)
kml::choice(sdq_cld)

kml::plot(ext_cld, 3,parTraj=parTRAJ(col="clusters"))
kml::plot(int_cld, 3,parTraj=parTRAJ(col="clusters"))
kml::plot(sdq_cld, 3,parTraj=parTRAJ(col="clusters"))

# Inspecting clustering criterias for 3,4,5 cluster choices
# Internalized behavioral problem clustering
int_cld@c1[[1]]@criterionValues
int_cld@c2[[1]]@criterionValues
int_cld@c3[[1]]@criterionValues
int_cld@c4[[1]]@criterionValues
int_cld@c5[[1]]@criterionValues

# Externalized behavioral problem clustering
ext_cld@c2[[1]]@criterionValues
ext_cld@c3[[1]]@criterionValues
ext_cld@c4[[1]]@criterionValues
ext_cld@c5[[1]]@criterionValues

# Extract clusters
int_cl_v <- getClusters(int_cld, nbCluster = 3)
ext_cl_v <- getClusters(ext_cld, nbCluster = 3)
sdq_cl_v <- getClusters(sdq_cld, nbCluster = 3)

# Combine the data with new clustering variables with the rest of the data
sdq_w_cl <- cbind(cl_int= int_cl_v, cl_ext= ext_cl_v, cl_sdq = sdq_cl_v,
                       sdq_w)
sdq_w_cl %>% group_by(cl_ext) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_int) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_sdq) %>% summarize(n()/50)

sdq_w_cl <- sdq_w_cl %>%  mutate(cl_ext=recode(cl_ext, 
                    `A`="B",
                    `B`="A"),
                    cl_sdq=recode(cl_sdq, 
                    `A`="B",
                    `B`="A")) 
sdq_w_cl$cl_ext <- fct_relevel(sdq_w_cl$cl_ext, "A")

# Inspect distribution in l=clusters
sdq_w_cl %>% group_by(cl_ext) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_int) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_sdq) %>% summarize(n()/50)

# Give meaningful names to clusters instead of A, B, C
sdq_w_cl <- sdq_w_cl %>% 
  mutate(cl_ext_str = case_when(cl_ext == 'A' ~ 'Profile 3: high and stable',
                                cl_ext == 'B' ~ 'Profile 2: mid and stable',
                                cl_ext == 'C' ~ 'Profile 1: low and decreasing'),
         cl_int_str = case_when(cl_int == 'A' ~ 'Profile 2: mid and decreasing',
                                cl_int == 'B' ~ 'Profile 3: high and increasing',
                                cl_int == 'C' ~ 'Profile 1: low and stable'),
         cl_sdq_str = case_when(cl_sdq == 'A' ~ 'Low',
                                cl_sdq == 'B' ~ 'Mid',
                                cl_sdq == 'C' ~ 'High'),
         cl_ext_str = factor(cl_ext_str, 
                levels = c('Profile 1: low and decreasing', 'Profile 2: mid and stable','Profile 3: high and stable'), 
                ordered = T),
         cl_int_str = factor(cl_int_str, 
                             levels = c('Profile 1: low and stable', 'Profile 2: mid and decreasing','Profile 3: high and increasing'), 
                             ordered = T),
         cl_sdq_str = factor(cl_sdq_str, 
                             levels = c('Low', 'Mid','High'), 
                             ordered = T)
         )
sdq_w_cl %>% group_by(cl_ext_str) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_int_str) %>% summarize(n()/50)
sdq_w_cl %>% group_by(cl_sdq_str) %>% summarize(n()/50)

(clust_freq_tbl <- sdq_w_cl %>% group_by(cl_ext, cl_int) %>% 
  summarize(count = n()) %>% 
  spread(., cl_int, count) %>%
  kable("html", digits = 2, caption = "Int and Ext behavioral problems clusters") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE))

# reshape to long
sdq_l_cl <-  melt(setDT(sdq_w_cl), id.vars	= c('ID', 'c_male', 'c_age', 'c_fborn', 'm_nat_jo', 
                                                   'm_working', 'Income',
                                                   'm_edu_uni', 'cl_int', 'cl_ext', 'cl_sdq', 'cl_int_str', 'cl_ext_str', 'cl_sdq_str'), 
                  measure.vars = list(c(12:15), c(16:19), c(20:23), c(24:27), c(28:31), c(32:35),
                                      c(36:39), c(40:43), c(44:47), c(48:51), c(52:55), c(56:59)),
                  value.name = c('SDQ_Int', 'SDQ_Ext', 's_SDQ_Int', 's_SDQ_Ext', 'HCC_C', 'HCC_M', #'GSE', 
                                 'PHQ', 'ParStr', #'SDQ_TotDiff_',
                                 'CRPR_W', 'CRPR_C', 'd_CRPR', 'SDQ_TotDiff'),
                  variable.name = 'wave') %>%
  mutate(wave_num = as.numeric(wave),
         wave_t = case_when(wave == 1 ~ "t1 (2019)",
                            wave == 2 ~ "t2 (June 2020)",
                            wave == 3 ~ "t3 (December 2020)",
                            wave == 4 ~ "t4 (June 2021)"))
dim(sdq_l_cl)

# Describing the clusters in terms of demographic variables
# Running equal means test (ANOVA) on demographic variables

# Male child gender (y/n)
aov(c_male ~ cl_int, data = sdq_w_cl) %>% summary()
aov(c_male ~ cl_ext, data = sdq_w_cl) %>% summary()

# Age
aov(c_age ~ cl_int, data = sdq_w_cl) %>% summary()
aov(c_age ~ cl_ext, data = sdq_w_cl) %>% summary()

# First born child (y/n)
aov(c_fborn ~ cl_int, data = sdq_w_cl) %>% summary()
aov(c_fborn ~ cl_ext, data = sdq_w_cl) %>% summary()

# Jordanian national (y/n)
aov(m_nat_jo ~ cl_int, data = sdq_w_cl) %>% summary()
aov(m_nat_jo ~ cl_ext, data = sdq_w_cl) %>% summary()

# Working mother (y/n)
aov(m_working ~ cl_int, data = sdq_w_cl) %>% summary()
aov(m_working ~ cl_ext, data = sdq_w_cl) %>% summary()

# Income (povery, low, mid-high)
aov(Income ~ cl_int, data = sdq_w_cl) %>% summary()
aov(Income ~ cl_ext, data = sdq_w_cl) %>% summary()

# Uni and higher mother education (y/n)
aov(m_edu_uni ~ cl_int, data = sdq_w_cl) %>% summary()
aov(m_edu_uni ~ cl_ext, data = sdq_w_cl) %>% summary()

############################
# Graphical inspection
############################
# Distribution across clusters (Figure 1)
########################
grid.arrange(
sdq_l_cl %>%
  ggplot() +
    geom_boxplot(aes(cl_int_str, SDQ_Int, linetype = wave_t)) +
    ylim(c(0.9,2.7)) +
    jtools::theme_apa() +
    scale_linetype_manual(name = "wave", values = c("solid", "dotted", "longdash", "twodash")) +
    scale_alpha_manual(name = "Community", values = c(1, 0.5)) +
    theme(legend.position="right", text = element_text(family = "Times New Roman"), 
        axis.text.x = element_text(face = "italic")) +    labs(x = "", y = "") +
    ggtitle('(a) Internalizing problems'),
sdq_l_cl %>%
  ggplot() +
    geom_boxplot(aes(cl_ext_str, SDQ_Ext, linetype = wave_t)) +
    ylim(c(0.9,2.7)) +
    scale_color_grey() + theme_classic() +
    jtools::theme_apa() +
    scale_linetype_manual(name = "wave", values = c("solid", "dotted", "longdash", "twodash")) +
    theme(legend.position="right", text = element_text(family = "Times New Roman"), 
          axis.text.x = element_text(face = "italic")) +
    labs(x = "", y = "") +
    ggtitle('(b) Externalizing problems'), 
ncol = 1, nrow = 2)



# Formal test of means of clusters
###################################
# Internalizing behavioral problems variable
lapply(1:4, function(x){aov(SDQ_Int ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()}) 
lapply(1:4, function(x){pairwise.t.test(sdq_l_cl[which(sdq_l_cl$wave == x)]$SDQ_Int, sdq_l_cl[which(sdq_l_cl$wave == x)]$cl_int, p.adjust.method = "bonf")})
aov(SDQ_Int ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == 1)]) %>% lsr::etaSquared()
aov(SDQ_Int ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == 2)]) %>% lsr::etaSquared()
aov(SDQ_Int ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == 3)]) %>% lsr::etaSquared()
aov(SDQ_Int ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == 4)]) %>% lsr::etaSquared()

# Externalizing behavioral problems variable
lapply(1:4, function(x){aov(SDQ_Ext ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})
lapply(1:4, function(x){pairwise.t.test(sdq_l_cl[which(sdq_l_cl$wave == x)]$SDQ_Ext, sdq_l_cl[which(sdq_l_cl$wave == x)]$cl_ext, p.adjust.method = "bonf")})
aov(SDQ_Ext ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == 1)]) %>% lsr::etaSquared()
aov(SDQ_Ext ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == 2)]) %>% lsr::etaSquared()
aov(SDQ_Ext ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == 3)]) %>% lsr::etaSquared()
aov(SDQ_Ext ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == 4)]) %>% lsr::etaSquared()

# Formal tests of equal means between cluster across the 4 rounds
####################################################
# Internalizing behavioral problem clustering
lapply(1:4, function(x){aov(CRPR_W ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})
lapply(1:4, function(x){aov(CRPR_C ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})
lapply(1:4, function(x){aov(d_CRPR ~ cl_int, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})

# Externalizing behavioral problem clustering
lapply(1:4, function(x){aov(CRPR_W ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})
lapply(1:4, function(x){aov(CRPR_C ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})
lapply(1:4, function(x){aov(d_CRPR ~ cl_ext, dat = sdq_l_cl[which(sdq_l_cl$wave == x)]) %>% summary()})

control_vars2 <- c('c_male', 'c_age', 'Income')

# Is the data changing over time?
#################################
# ANOVA tests to see if the means of the clustering variables change over time
aov(SDQ_Int ~ wave, data = sdq_l) %>% summary()
aov(SDQ_Ext ~ wave, data = sdq_l) %>% summary()

# ANOVA tests to see if the means of the parent variables change over time
aov(PHQ ~ wave, data = sdq_l) %>% summary()
aov(CRPR_C ~ wave, dat = sdq_l) %>% summary()
aov(CRPR_W ~ wave, dat = sdq_l) %>% summary()
aov(d_CRPR ~ wave, dat = sdq_l) %>% summary()

# Graphical representation
##########################
# Correlation matrix of variables
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)  
}

(cor_mat_sdq <- cor_mat(sdq_w[,9:16], 
                        method = "pearson",
                        alternative = "two.sided",
                        conf.level = 0.95) %>% 
    cor_mark_significant(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "")) %>%
    kable("html", digits = 2, caption = "Correlation matrix: SDQ variables") %>%
    kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive")))

cor_mat_scales <- cor_mat(sdq_w[,c(9:16,35, 36, 39:(NCOL(sdq_w)))], 
                          method = "pearson",
                          alternative = "two.sided",
                          conf.level = 0.95) %>% 
  cor_mark_significant(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "")) %>%
  kable("html", digits = 2, caption = "Correlation matrix: ") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive")) %>% print()

cor_mat_plt_sdq <- pairs(sdq_w[,9:16], lower.panel = panel.smooth, upper.panel = panel.cor, main = "Correlation between internalizing and externalizing behavior problems")


# Descriptive stats by profile/cluster
#######################################
(ext_cl_stats_tbl <- sdq_w_cl %>% 
  group_by(cl_ext_str) %>%
  summarize_at(control_vars_plus, list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never") %>%
  kable("html", digits = 2, caption = "Samply information by externalizing behavioral problems profiles") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE))

(int_cl_stats_tbl <- sdq_w_cl %>% 
  group_by(cl_int_str) %>%
  summarize_at(control_vars_plus, list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never") %>%
  kable("html", digits = 2, caption = "Samply information by internalizing behavioral problems profiles") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE))

(sdq_cl_stats_tbl <- sdq_w_cl %>% 
  group_by(cl_sdq_str) %>%
  summarize_at(control_vars_plus, list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never") %>%
  kable("html", digits = 2, caption = "Samply information by global SDQ scale profiles") %>%
  kable_styling(bootstrap_options=c("striped", "hover", "condensed", "responsive"),
                full_width=FALSE))

control_vars_names <- c('Child gender: male', 'Child age', 'Child first born', 'Nationality: Jo', 'Mother currently working', 'Income', 'Mother\'s education: uni and higher')

bind_rows(names = t(control_vars_names), sdq_w_cl %>% 
            group_by(cl_ext_str) %>%
            summarize_at(control_vars, list(mean_sd), digits = 3, denote_sd = "paren"))
int_avg  <- sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(control_vars, list(avg = mean), na.rm = T) %>% print()
int_sd  <- sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(control_vars, list(sd = sd), na.rm = T) %>% print()
bind_cols(tibble(Category = int_avg[,1]), 
          print(int_avg[, 2:8], "(", int_sd[, 2:8], ")", sep = ''))

# What predicts group membership?
# Multinomial logistic regression
##################################
# Exploring potential relationships between demographic characteristics and 
# children’s mental health trajectories by cluster membership. 
multinom_pivot_wider <- function(x) {
if (!inherits(x, "tbl_regression") || !inherits(x$model_obj, "multinom")) {
stop("`x=` must be class 'tbl_regression' summary of a `nnet::multinom()` model.")
}
# create tibble of results
df <- tibble::tibble(outcome_level = unique(x$table_body$groupname_col))
df$tbl <-
purrr::map(
df$outcome_level,
function(lvl) {
gtsummary::modify_table_body(
x,
~dplyr::filter(.x, .data$groupname_col %in% lvl) %>%
dplyr::ungroup() %>%
dplyr::select(-.data$groupname_col)
)
}
)
tbl_merge(df$tbl, tab_spanner = paste0("**", df$outcome_level, "**"))
}

# Tabulating model results for each clustering variable
######################################################
# Externalizing behavioral problem clustering
(multinom_ext <- multinom(cl_ext_str ~  c_male + c_age + Income, data = sdq_w_cl) %>% 
    tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE, label = multinom_label) %>%
    multinom_pivot_wider())

(multinom_ext_red <- multinom(cl_ext_str ~  c_male + c_age + Income, data = sdq_w_cl %>% dplyr::filter(cl_ext != 'A')) %>% 
    tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE, label = multinom_label) %>%
    multinom_pivot_wider())

# Internalizing behavioral problem clustering
multinom_label <- list(c_age ~ "Child age", c_male ~ "Child gender: male", c_male ~ "Child gender: male")
(multinom_int <- multinom(cl_int_str ~  c_male + c_age + Income, data = sdq_w_cl) %>% 
    tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE, label = multinom_label) %>%
    multinom_pivot_wider())

(multinom_int_red <- multinom(cl_int_str ~  c_male + c_age + Income, data = sdq_w_cl %>% dplyr::filter(cl_int != 'A')) %>% 
    tbl_regression(exponentiate = T, add_estimate_to_reference_rows = TRUE, label = multinom_label) %>%
    multinom_pivot_wider())


# Descriptive stats for response variables
##########################################
# Clustered by externalizing behavioral problem
###
# Parent health questionnaire response variable
sdq_w_cl %>% group_by(cl_ext_str) %>%
  summarize_at(c("PHQ_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

# Parent stress response variable
sdq_w_cl %>% group_by(cl_ext_str) %>%
  summarize_at(c("ParStr_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

# CRPR_W response variable
sdq_w_cl %>% group_by(cl_ext_str) %>%
  summarize_at(c("CRPR_W_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

# CRPR_C response variable
sdq_w_cl %>% group_by(cl_ext_str) %>%
  summarize_at(c("CRPR_C_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

# Clustered by internalizing behavioral problem
###
sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(c("PHQ_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(c("ParStr_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(c("CRPR_W_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")

sdq_w_cl %>% group_by(cl_int_str) %>%
  summarize_at(c("CRPR_C_4"), list(mean_sd), digits = 2, denote_sd = "paren", na_rm = T, show_n = "never")



