## ----setup, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------
# For ADNIMERGE, go to http://adni.loni.usc.edu/, https://adni.bitbucket.io/

library(Hmisc)
library(knitr)
library(tidyverse)
library(kableExtra)
library(gridExtra)
library(plotly)
library(nlme)
library(emmeans)
library(arsenal)
library(grid)
library(gridExtra)
library(mvtnorm)
library(mice)

options(digits=2)

theme_set(theme_bw())

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbbPalette <-
    c("#0072B2", "#D55E00", "#E69F00",
      "#009E73", "#F0E442", "#999999",
      "#000000", "#56B4E9", "#CC79A7")
scale_colour_discrete <-
    function(...) scale_colour_manual(..., values = cbbPalette)
scale_fill_discrete <-
    function(...) scale_fill_manual(..., values = cbbPalette)
scale_colour_discrete <-
    function(...) scale_colour_manual(..., values = cbbPalette)
scale_fill_discrete <-
    function(...) scale_fill_manual(..., values = cbbPalette)

theme_table <- function(..., levs=2){
  theme_minimal(...) + 
    theme(
      panel.grid = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_text(face='bold', color=cbbPalette[1:levs]),
      axis.title = element_blank())
}

center <- function(x) scale(x, scale = FALSE)










## ---------------------------------------------------------------------------------------------------------------
dd2 <- subset(ADNIMERGE::adnimerge, VISCODE=='bl' & !is.na(Hippocampus) & DX=='CN') %>%
  mutate(APOEe4 = ifelse(APOE4>0, 1, 0))


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
lm_fit1_gender <- lm(I(Hippocampus/ICV*100) ~ PTGENDER, data=dd2)


## ---------------------------------------------------------------------------------------------------------------
summary(lm_fit1_gender)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|t|)` = format.pval(round(`Pr(>|t|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
lm_fit2_gender <- lm(Hippocampus ~ ICV + PTGENDER, data=dd2)


## ---------------------------------------------------------------------------------------------------------------
summary(lm_fit2_gender)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|t|)` = format.pval(round(`Pr(>|t|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----hipp_gender_scatter, fig.height=4, fig.width=6-------------------------------------------------------------
ggplot(dd2, aes(y=Hippocampus, x=ICV, color=PTGENDER)) +
  geom_point(alpha=0.5) +
  geom_smooth(method='lm') +
  xlab(expression(paste("ICV (", m*m^{3},")"))) +
  ylab(expression(paste("Hippocampus (", m*m^{3},")"))) +
  theme(legend.position = c(0.2, 0.8))


## ----hipp_gender_box, fig.height=4, fig.width=4*2---------------------------------------------------------------
p1 <- ggplot(dd2, aes(y=Hippocampus, x=PTGENDER)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=100) +
  xlab('') + ylab(expression(paste("Hippocampus (", m*m^{3},")")))

p2 <- ggplot(dd2, aes(y=ICV, x=PTGENDER)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=15000) +
  xlab('') + ylab(expression(paste("ICV (", m*m^{3},")")))

p3 <- ggplot(dd2, aes(y=(Hippocampus/ICV*100), x=PTGENDER)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=0.0075) +
  xlab('') + ylab("Hippocampus (%ICV)")

p4 <- ggplot(dd2, 
  aes(y=Hippocampus-ICV*coef(lm_fit2_gender)['ICV'], x=PTGENDER)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', 
    dotsize=0.3, alpha=0.2, binwidth=100) +
  xlab('') + ylab(expression(paste("Hippocampus (", m*m^{3},", model-adjusted for ICV)")))

grid.arrange(p1,p2,p3,p4, nrow=1)


## ----hipp_adj_gender_scatter, fig.height=4, fig.width=6---------------------------------------------------------
ggplot(dd2, 
  aes(y=Hippocampus-ICV*coef(lm_fit2_gender)['ICV'], x=ICV, color=PTGENDER)) +
  geom_point(alpha=0.5) +
  geom_smooth(method='lm') +
  xlab(expression(paste("ICV (", m*m^{3},")"))) +
  ylab(expression(paste("Hippocampus (", m*m^{3},", model-adjusted for ICV)"))) +
  theme(legend.position = c(0.2, 0.8))


## ---------------------------------------------------------------------------------------------------------------
ggplot() + 
  geom_segment(aes(x=1.009, y=0, xend=2.076, yend=100),
    arrow = arrow(length = unit(0.03, "npc"), ends='both')) +
  geom_abline(slope = 100/(2.076-1.009), 
    intercept=-1.009*100/(2.076-1.009), linetype='dashed') +
  geom_point(aes(x=1.009, y=0)) +
  geom_point(aes(x=2.076, y=100)) +
  xlab('PiB SUVr') + ylab('Centiloid') +
  xlim(1, 3) + ylim(-5, 200)
   






## ----echo=FALSE-------------------------------------------------------------------------------------------------
pibids <- unique(subset(ADNIMERGE::adnimerge, !is.na(PIB))$RID)
av45ids <- unique(subset(ADNIMERGE::adnimerge, !is.na(AV45))$RID)
set.seed(20210506)
holdout <- sample(intersect(pibids, av45ids), size = 10)

dd <- ADNIMERGE::adnimerge %>%
  filter(!RID %in% holdout) %>%
  arrange(RID, EXAMDATE) %>%
  select(RID, DX, PIB, AV45) %>%
  rename(PiB = PIB, Florbetapir = AV45) %>%
  mutate(PiB = as.numeric(PiB), Florbetapir = as.numeric(Florbetapir)) %>%
  group_by(RID) %>%
  fill(DX, .direction = "downup") %>%
  filter(!is.na(PiB) | !is.na(Florbetapir)) %>%
  filter(!duplicated(RID)) %>%
  pivot_longer(cols=c('PiB', 'Florbetapir'), names_to='Tracer', values_to='SUVR') %>%
  mutate(Tracer = factor(Tracer, levels = c('PiB', 'Florbetapir'))) %>%
  filter(!is.na(SUVR))
  
ggplot(dd, aes(x=SUVR)) +
  stat_ecdf(geom = "step") +
  facet_grid(.~Tracer, scales='free_x')


## ---- echo=FALSE------------------------------------------------------------------------------------------------
t1 <- with(dd, table(DX, Tracer))
t2 <- round(with(dd, prop.table(table(DX, Tracer), margin = 2))*100, 1)
tt <- t1
tt[,1] <- paste(t1[,1], paste0("(", t2[,1], "%)"))
tt[,2] <- paste(t1[,2], paste0("(", t2[,2], "%)"))
kable(tt)


## ---------------------------------------------------------------------------------------------------------------
invproptab <- 1/with(dd, prop.table(table(DX, Tracer), margin = 2))
kable(invproptab)


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
# Record the sampling adjustment weights in the data
dd <- dd %>% mutate(
  wt = case_when(
    DX == 'CN' & Tracer == 'PiB' ~ invproptab['CN', 'PiB'],
    DX == 'MCI' & Tracer == 'PiB' ~ invproptab['MCI', 'PiB'],
    DX == 'Dementia' & Tracer == 'PiB' ~ invproptab['Dementia', 'PiB'],
    DX == 'CN' & Tracer == 'Florbetapir' ~ invproptab['CN', 'Florbetapir'],
    DX == 'MCI' & Tracer == 'Florbetapir' ~ invproptab['MCI', 'Florbetapir'],
    DX == 'Dementia' & Tracer == 'Florbetapir' ~ invproptab['Dementia', 'Florbetapir']
  ))


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
# Create adjusted ECDF functions (mapping SUVRs to Cumulative Probabilities)
# Hmisc::wtd.Ecdf returns a data.frame evaluating the ECDF at each observed value
PiB.ecdf.data <- with(subset(dd, Tracer == 'PiB'), 
  Hmisc::wtd.Ecdf(SUVR, weights=wt, normwt=TRUE))
Fbp.ecdf.data <- with(subset(dd, Tracer == 'Florbetapir'), 
  Hmisc::wtd.Ecdf(SUVR, weights=wt, normwt=TRUE))

# approxfun creates a function via linear interpolation 
# mapping SUVRs to cumulative probabilities (0 to 1 scale)
PiB.ecdf <- with(PiB.ecdf.data, approxfun(x, ecdf, rule=2))
Fbp.ecdf <- with(Fbp.ecdf.data, approxfun(x, ecdf, rule=2))

# Create adjusted **inverse** ECDF functions 
# mapping Cumulative Probabilities (0 to 1 scale) to SUVRs
Probs <- seq(0,1,by=0.01)
PiB.inv.ecdf.data <- as.numeric(with(subset(dd, Tracer == 'PiB'), 
  Hmisc::wtd.quantile(SUVR, weights=wt, normwt=TRUE, probs=Probs)))
Fbp.inv.ecdf.data <- as.numeric(with(subset(dd, Tracer == 'Florbetapir'), 
  Hmisc::wtd.quantile(SUVR, weights=wt, normwt=TRUE, probs=Probs)))

PiB.inv.ecdf <- approxfun(Probs, PiB.inv.ecdf.data, rule=2)
Fbp.inv.ecdf <- approxfun(Probs, Fbp.inv.ecdf.data, rule=2)


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
dd <- dd %>% mutate(
  `Adjusted cumulative probability` = case_when( # 
    Tracer == 'PiB' ~ PiB.ecdf(SUVR),
    Tracer == 'Florbetapir' ~ Fbp.ecdf(SUVR)),
  `Adjusted Z-score` = qnorm(`Adjusted cumulative probability`), # adjusted z-scores
  `Florbetapir to PiB adjusted SUVR` = case_when(
    Tracer == 'Florbetapir' ~ PiB.inv.ecdf(Fbp.ecdf(SUVR))),
  `PiB to Florbetapir adjusted SUVR` = case_when(
    Tracer == 'PiB' ~ Fbp.inv.ecdf(PiB.ecdf(SUVR))),
  CL = case_when(
    Tracer == 'PiB' ~ 100*(SUVR - 1.009)/1.067)) %>%
  arrange(Tracer, SUVR)


## ----weighted-ecdfs---------------------------------------------------------------------------------------------
ggplot(dd, aes(x=SUVR)) +
  stat_ecdf(geom = "step") +
  facet_grid(.~Tracer, scales='free_x') +
  geom_line(aes(x=SUVR, y=`Adjusted cumulative probability`), linetype='dashed') +
  ylab('Cumulative Probability')


## ---------------------------------------------------------------------------------------------------------------
dd %>% filter(Tracer == 'PiB') %>%
  ggplot(aes(x=CL, y=`Adjusted cumulative probability`)) +
  geom_line() +
  geom_point()


## ----pib-densities----------------------------------------------------------------------------------------------
dd %>% select(RID, DX, Tracer, SUVR, `Florbetapir to PiB adjusted SUVR`) %>%
  pivot_longer(c('SUVR', 'Florbetapir to PiB adjusted SUVR'), 
    names_to = 'Source', values_to = 'SUVR') %>%
  filter( (Tracer == 'PiB' & Source == 'SUVR') |
      (Tracer == 'Florbetapir' & Source == 'Florbetapir to PiB adjusted SUVR')) %>%
  mutate(Source = ifelse(Source == 'SUVR', 
    'Directly observed PiB SUVR', 'ECDF mapped Florbetapir to PiB SUVR')) %>%
ggplot(aes(x=SUVR, color=Source)) +
  geom_density(alpha=0.5)


## ----pib-densities-dx-------------------------------------------------------------------------------------------
dd %>% select(RID, DX, Tracer, SUVR, `Florbetapir to PiB adjusted SUVR`) %>%
  pivot_longer(c('SUVR', 'Florbetapir to PiB adjusted SUVR'), 
    names_to = 'Source', values_to = 'SUVR') %>%
  filter( (Tracer == 'PiB' & Source == 'SUVR') |
      (Tracer == 'Florbetapir' & Source == 'Florbetapir to PiB adjusted SUVR')) %>%
  mutate( Source = ifelse(Source == 'SUVR', 'Directly observed PiB SUVR', 'ECDF mapped Florbetapir to PiB SUVR')) %>%
  filter(!is.na(DX)) %>%
ggplot(aes(x=SUVR, color=Source)) +
  geom_density(alpha=0.5) +
  facet_grid(.~DX)


## ----z-score-densities------------------------------------------------------------------------------------------
ggplot() +
  geom_rug(data = subset(dd, Tracer == 'Florbetapir'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_rug(data = subset(dd, Tracer == 'PiB'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_density(data = subset(dd, Tracer == 'Florbetapir'), 
    aes(x=`Adjusted Z-score`, color = Tracer), alpha=0.5) +
  geom_density(data = subset(dd, Tracer == 'PiB'), aes(x=`Adjusted Z-score`, color = Tracer))




## ---------------------------------------------------------------------------------------------------------------
dd.validate <- full_join(
  ADNIMERGE::adnimerge %>% 
    arrange(RID, EXAMDATE) %>% 
    select(RID, EXAMDATE, PIB) %>%
    rename(PiB = PIB) %>%
    filter(RID %in% holdout) %>%
    filter(!is.na(PiB)) %>%
    filter(!duplicated(RID)),
  ADNIMERGE::adnimerge %>% 
    arrange(RID, EXAMDATE) %>% 
    select(RID, EXAMDATE, AV45) %>%
    rename(Florbetapir = AV45) %>%
    filter(RID %in% holdout) %>%
    filter(!is.na(Florbetapir)) %>%
    filter(!duplicated(RID)), by='RID', suffix = c('.PiB', '.Fbp')
  ) %>%
  mutate(
    `Years between scans` = as.numeric(EXAMDATE.Fbp - EXAMDATE.PiB)/365.25,
    `Royse et al linear map of PiB to Florbetapir` = 
      as.numeric(0.497 * PiB + 0.503),
    `Royse et al linear map of Florbetapir to PiB` = 
      as.numeric((Florbetapir - 0.503)/0.497),
    `Navitsky et al linear map of PiB to Florbetapir` = 
      as.numeric(0.536 * PiB + 0.502),
    `Navitsky et al linear map of Florbetapir to PiB` = 
      as.numeric((Florbetapir - 0.502)/0.536),
    `ECDF map of Florbetapir to PiB` = PiB.inv.ecdf(Fbp.ecdf(Florbetapir)),
    `ECDF map of PiB to Florbetapir` = Fbp.inv.ecdf(PiB.ecdf(PiB))
  )

fbb2pib_navitsky <- dd.validate %>%
  select(RID, `Years between scans`, PiB, 
    `Navitsky et al linear map of Florbetapir to PiB`,
    `ECDF map of Florbetapir to PiB`) %>%
  rename(
    `Navitsky et al linear map*` = `Navitsky et al linear map of Florbetapir to PiB`,
    `ECDF map` = `ECDF map of Florbetapir to PiB`) %>%
  pivot_longer(c('Navitsky et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated PiB SUVR')

fbb2pib_navitsky_RMSE <- fbb2pib_navitsky %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((PiB-`Estimated PiB SUVR`)^2))) %>%
  mutate(PiB = 2, `Estimated PiB SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(fbb2pib_navitsky, aes(x=PiB, y=`Estimated PiB SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual PiB SUVR') +
  geom_text(data=fbb2pib_navitsky_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 1.5, linetype='dashed')


## ---------------------------------------------------------------------------------------------------------------
pib2fbb_navitsky <- dd.validate %>%
  select(RID, `Years between scans`, Florbetapir, 
    `Navitsky et al linear map of PiB to Florbetapir`,
    `ECDF map of PiB to Florbetapir`) %>%
  rename(
    `Navitsky et al linear map*` = `Navitsky et al linear map of PiB to Florbetapir`,
    `ECDF map` = `ECDF map of PiB to Florbetapir`) %>%
  pivot_longer(c('Navitsky et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated Florbetapir SUVR')

pib2fbb_navitsky_RMSE <- pib2fbb_navitsky %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((Florbetapir-`Estimated Florbetapir SUVR`)^2))) %>%
  mutate(Florbetapir = 1.5, `Estimated Florbetapir SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(pib2fbb_navitsky, aes(x=Florbetapir, y=`Estimated Florbetapir SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual Florbetapir SUVR') +
  geom_text(data=pib2fbb_navitsky_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.11, linetype='dashed') +
  geom_vline(xintercept = 1.11, linetype='dashed')


## ---------------------------------------------------------------------------------------------------------------
fbb2pib_royse <- dd.validate %>%
  select(RID, `Years between scans`, PiB, 
    `Royse et al linear map of Florbetapir to PiB`,
    `ECDF map of Florbetapir to PiB`) %>%
  rename(
    `Royse et al linear map*` = `Royse et al linear map of Florbetapir to PiB`,
    `ECDF map` = `ECDF map of Florbetapir to PiB`) %>%
  pivot_longer(c('Royse et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated PiB SUVR')

fbb2pib_royse_RMSE <- fbb2pib_royse %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((PiB-`Estimated PiB SUVR`)^2))) %>%
  mutate(PiB = 2, `Estimated PiB SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(fbb2pib_royse, aes(x=PiB, y=`Estimated PiB SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual PiB SUVR') +
  geom_text(data=fbb2pib_royse_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.5, linetype='dashed') +
  geom_vline(xintercept = 1.5, linetype='dashed')


## ---------------------------------------------------------------------------------------------------------------
pib2fbb_royse <- dd.validate %>%
  select(RID, `Years between scans`, Florbetapir, 
    `Royse et al linear map of PiB to Florbetapir`,
    `ECDF map of PiB to Florbetapir`) %>%
  rename(
    `Royse et al linear map*` = `Royse et al linear map of PiB to Florbetapir`,
    `ECDF map` = `ECDF map of PiB to Florbetapir`) %>%
  pivot_longer(c('Royse et al linear map*', 'ECDF map'), 
    names_to = 'Method', values_to = 'Estimated Florbetapir SUVR')

pib2fbb_royse_RMSE <- pib2fbb_royse %>% group_by(Method) %>%
  summarise(RMSE = sqrt(mean((Florbetapir-`Estimated Florbetapir SUVR`)^2))) %>%
  mutate(Florbetapir = 1.5, `Estimated Florbetapir SUVR` = 0.75,
    text = paste0('RMSE=', round(RMSE, digits=3)))

ggplot(pib2fbb_royse, aes(x=Florbetapir, y=`Estimated Florbetapir SUVR`)) +
  geom_point(aes(color=`Years between scans`)) +
  geom_abline(intercept = 0, slope=1) +
  facet_grid(.~Method) +
  xlab('Actual Florbetapir SUVR') +
  geom_text(data=pib2fbb_royse_RMSE, aes(label=text)) +
  geom_hline(yintercept = 1.11, linetype='dashed') +
  geom_vline(xintercept = 1.11, linetype='dashed')






## ----echo=FALSE-------------------------------------------------------------------------------------------------
fit_adni <- lme(ADAS11 ~ PTGENDER + center(AGE) + I(Years.bl*12), 
  data=ADNIMERGE::adnimerge,
  random=~M|RID, subset = M<=24 & DX.bl=='AD', na.action=na.omit)
summary(fit_adni)


## ----generate_data1, echo=TRUE----------------------------------------------------------------------------------
# fixed effects parameters estimated from ADNI
Beta <- c(
   '(Intercept)'=20 , # mean ADAS at baseline
        'female'=0.1, # worse scores for females
         'age_c'=0.1, # worse change for older at baseline (age mean centered)
         'month'=0.4, # worse per month post baseline
  'month:active'=-0.05) # improvement per month with treatment

# standard deviations for random effects
sigma_random_intercept <- 6.0
sigma_random_slope <- 0.37
sigma_residual <- 3.3

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.40/18 # approx per month


## ---------------------------------------------------------------------------------------------------------------
# set seed so that simulation is reproducible
set.seed(20170701)

# simulate subject specific data
subjects <- data.frame(
  id = 1:(2*n),
  active = sample(c(rep(0,n), rep(1,n)), 2*n),
  female = sample(0:1, 2*n, replace=TRUE),
  age = rnorm(2*n, 75, 7.8), 
  censor = rexp(2*n,rate=attrition_rate), # attrition month
  ran.intercept = rnorm(2*n, sd=sigma_random_intercept),
  ran.slope     = rnorm(2*n, sd=sigma_random_slope)) %>%
  mutate(
    age_c = age-mean(age),
    Female = factor(female, levels = c(1, 0), labels = c('Female', 'Male')),
  )

# simulate data over time
trial <- right_join(subjects, 
  expand.grid(id = 1:(2*n), month=months)) %>%
  mutate(
    residual = rnorm(2*n*length(months), sd=sigma_residual),
    group = factor(active, 0:1, c('placebo', 'active')),
    missing = ifelse(month>censor, 1, 0)) %>%
  arrange(id, month) %>%
  select(id, month, everything())

# calculate the ADAS scores with random effects and residuals and 
# round to nearest digit in 0-70
trial <- mutate(trial,
  ADAS11 = (model.matrix(~ female+age_c+month+month:active, trial)[, names(Beta)] %*% Beta)[,1],
  ADAS11 = round(ADAS11 + ran.intercept + ran.slope*month + residual, digits = 0),
  ADAS11 = replace(ADAS11, ADAS11<0, 0),
  ADAS11 = replace(ADAS11, ADAS11>70, 70))

# filter out the missing observations
trial_obs <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_obs %>%
  select(id, month, female, age, age_c, active, group, ADAS11) %>% 
  mutate(month = paste0('ADAS11.m', month)) %>%
  pivot_wider(names_from = month, values_from = ADAS11)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS11.m0), 
  filter(trial_obs, month>0)) %>%
  mutate(ADAS11.ch = ADAS11 - ADAS11.m0)


## ---------------------------------------------------------------------------------------------------------------
head(select(trial_obs, -group, -missing)) %>% kable()


## ----results='asis'---------------------------------------------------------------------------------------------
tab <- tableby(group ~ Female + age + ADAS11, 
  data = subset(trial_obs, month == 0))

summary(tab, text=TRUE) %>% kable()


## ----spaghetti_plot---------------------------------------------------------------------------------------------
ggplot(trial_obs, aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1, 0.85), legend.background = element_rect(fill=NA))


## ---------------------------------------------------------------------------------------------------------------
summaryTable <- trial_obs %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS11),
    mean=mean(ADAS11),
    sd=sd(ADAS11),
    lower95 = Hmisc::smean.cl.normal(ADAS11)[['Lower']],
    upper95 = Hmisc::smean.cl.normal(ADAS11)[['Upper']],
    min=min(ADAS11),
    max=max(ADAS11))
as.data.frame(summaryTable) %>% kable() %>% collapse_rows(1)


## ----meanplot---------------------------------------------------------------------------------------------------
p <- ggplot(summaryTable, aes(x=month, y=mean, color=group)) +
  geom_line() +
  geom_errorbar(aes(min=lower95, max=upper95), position=position_dodge(0.2), width=0) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
countTab <- ggplot(summaryTable, aes(x=month, y=group, label=n)) + geom_text() + theme_table()
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))


## ---------------------------------------------------------------------------------------------------------------
chsummaryTable <- trial_mmrm %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS11.ch),
    mean=mean(ADAS11.ch),
    sd=sd(ADAS11.ch),
    lower95 = Hmisc::smean.cl.normal(ADAS11.ch)[['Lower']],
    upper95 = Hmisc::smean.cl.normal(ADAS11.ch)[['Upper']],
    min=min(ADAS11.ch),
    max=max(ADAS11.ch))
chsummaryTable <- bind_rows(
  filter(summaryTable, month == 0) %>% 
    mutate(mean=0, sd=0, lower95=0, upper95=0, min=0, max=0),
  chsummaryTable) %>%
  arrange(group, month)

as.data.frame(chsummaryTable) %>% kable() %>% collapse_rows(1)


## ----meanchplot-------------------------------------------------------------------------------------------------
p <- ggplot(chsummaryTable, aes(x=month, y=mean, color=group)) +
  geom_line() +
  geom_errorbar(aes(min=lower95, max=upper95), position=position_dodge(0.2), width=0) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS change') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))


## ---------------------------------------------------------------------------------------------------------------
m1 <- subset(chsummaryTable, group=='active' & month==18)[['mean']]
m2 <- subset(chsummaryTable, group=='placebo' & month==18)[['mean']]
n1 <- subset(chsummaryTable, group=='active' & month==18)[['n']]
n2 <- subset(chsummaryTable, group=='placebo' & month==18)[['n']]
sd1 <- subset(chsummaryTable, group=='active' & month==18)[['sd']]
sd2 <- subset(chsummaryTable, group=='placebo' & month==18)[['sd']]
s <- sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1+n2-2))
tt <- (m2-m1)/(s*sqrt(1/n2 + 1/n1)) 
DF <- n1+n2-2


## ----echo=FALSE-------------------------------------------------------------------------------------------------
print(t.test(ADAS11.ch ~ group, data = trial_mmrm, subset = month==18, 
  var.equal=TRUE), digits = 6)


## ---------------------------------------------------------------------------------------------------------------
x <- seq(-5, 5, by = 0.01)
dens <- data.frame(
	x        =  x,
	density  = dt(x, df = DF)
)
shadel <- rbind(
  c(x=-5, y=0),
  filter(dens, x <= -abs(tt)),
  c(x=-abs(tt),0))
shadeu <- rbind(
  c(x=abs(tt),0),
  filter(dens, x >= abs(tt)),
  c(x=5, y=0))
ggplot(dens, aes(x=x, y=density)) + geom_line() +
  geom_polygon(data = shadel, aes(x=x, y=density), fill=cbbPalette[1]) +
  geom_polygon(data = shadeu, aes(x=x, y=density), fill=cbbPalette[1]) +
  geom_vline(xintercept=qt(c(0.025, 0.975), df=DF), color='grey', linetype='dashed') +
  annotate('text', x=-1.96, y=0.4, label='x=-1.96') +
  annotate('text', x=1.96, y=0.4, label='x=1.96')


## ----ancovai, echo = FALSE, size = 'scriptsize'-----------------------------------------------------------------
summary(lm(ADAS11.ch ~ active + ADAS11.m0, data = trial_mmrm))


## ----ancovaii, echo = FALSE, size = 'scriptsize'----------------------------------------------------------------
summary(lm(ADAS11.ch ~ active*center(ADAS11.m0), data = trial_mmrm))


## ----ancovaii2cov, echo = FALSE, size = 'scriptsize'------------------------------------------------------------
summary(lm(ADAS11.ch ~ active*center(ADAS11.m0) + female + age_c, data = trial_mmrm))


## ----trial_stage1_plot, echo = FALSE----------------------------------------------------------------------------
ggplot(subset(trial, id %in% 1:4),
  aes(x=month, y=ADAS11, group = id, color = group)) +
  stat_smooth(method = 'lm') + geom_line() + geom_point() +
  facet_wrap(~id)


## ----trial_fit_stage1, eval = TRUE, echo = FALSE, size = 'scriptsize'-------------------------------------------
trial_stage1 <- as.data.frame(do.call(rbind, lapply(unique(trial$id),
  function(i){
    fit <- lm(ADAS11 ~ month,
      data = trial_obs, subset = id == i)
    c(id = i, beta = fit$coef, sigma = summary(fit)$sigma)
})))
trial_stage1 <- right_join(trial_stage1,
  filter(trial, month == 0) %>% 
    select(id, active, group, age_c, female))
head(trial_stage1) %>% kable()


## ----size = 'scriptsize'----------------------------------------------------------------------------------------
summary(lm(ADAS11 ~ month, data = trial_obs, subset = id == 1))


## ----trial_plot_stage2------------------------------------------------------------------------------------------
ggplot(trial_stage1,
  aes(x=group, y=beta.month, group = group, color = group)) +
  geom_boxplot(alpha = 0.25) + 
  ylab('ADAS11 change per month')


## ----trial_stage2_fit-------------------------------------------------------------------------------------------
summary(lm(beta.month ~ female + age_c + active, data = trial_stage1))


## ---------------------------------------------------------------------------------------------------------------
save(center, trial, trial_obs, trial_wide, trial_mmrm, 
  months, attrition_rate, n, 
  Beta, countTab, summaryTable, file='simulated-trial.Rdata')

