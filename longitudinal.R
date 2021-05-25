## ----setup, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------
# For ADNIMERGE, go to http://adni.loni.usc.edu/, https://adni.bitbucket.io/

library(Hmisc)
library(knitr)
library(kableExtra)
library(gridExtra)
library(plotly)
library(nlme)
library(emmeans)
library(multcomp)
library(arsenal)
library(grid)
library(gridExtra)
library(mvtnorm)
library(mice)
library(tidyverse)

options(digits=3)

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

load('simulated-trial.Rdata')
trial_obs <- trial_obs %>%
  arrange(id, month)










## ----trial_lme, size = 'tiny'-----------------------------------------------------------------------------------
fit_lme <- lme(ADAS11 ~ month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme)


## ----trial_lme_age, size = 'tiny'-------------------------------------------------------------------------------
fit_lme_cov <- lme(ADAS11 ~ age_c + female + month + month:active, data = trial_obs, random = ~month|id)
summary(fit_lme_cov)


## ----trial_lme_rcode, eval = FALSE, echo = TRUE-----------------------------------------------------------------
lme(ADAS11 ~ month + month:active, 
  data = trial_obs, random = ~month|id)

lme(ADAS11 ~ age_c + female + month + month:active, 
  data = trial_obs, random = ~month|id)


## ----trial_lme_profiles, echo = FALSE, size = 'scriptsize'------------------------------------------------------
em <- fit_lme_cov %>%
  ref_grid(at = list(
    month = unique(trial_obs$month),
    active = unique(trial_obs$active),
    female = 1, age.c = 0)) %>%
  emmeans(specs = c('active', 'month')) %>%
  as.data.frame()


## ----echo = FALSE, size = 'scriptsize'--------------------------------------------------------------------------
em %>% kable()


## ----echo = FALSE, size = 'scriptsize'--------------------------------------------------------------------------
p <- em %>%
  mutate(group = factor(active, levels = c(0,1), labels = c('Placebo', 'Active'))) %>%
  ggplot(aes(x = month, y = emmean, group = group)) +
  geom_line(aes(color=group)) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill=group), alpha=0.25) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))


## ----echo = FALSE, size = 'scriptsize'--------------------------------------------------------------------------
p <- em %>%
  mutate(group = factor(active, levels = c(0,1), labels = c('Placebo', 'Active'))) %>%
  ggplot(aes(x = month, y = emmean, group = group, color=group))+
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width=0, position=position_dodge(0.2)) +
  ylab('Mean ADAS (95% CI)') +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position=c(0.1, 0.75))
grid.draw(arrangeGrob(p,countTab,heights=c(3,1)))


## ----varPar, echo=FALSE-----------------------------------------------------------------------------------------
# simulate data with different variance parameters
varPar <- expand.grid(
  sigma_random_intercept = c(2, 10),
  sigma_random_slope = c(0.2, 0.75),
  sigma_residual = c(2, 8)
)

varPlot <- do.call(rbind, lapply(1:nrow(varPar), function(i){
  set.seed(20170714)
  subjects <- data.frame(
    id = 1:(2*n),
    active = sample(c(rep(0,n), rep(1,n)), 2*n),
    female = sample(0:1, 2*n, replace=TRUE),
    age = rnorm(2*n, 75, 7.8),
    censor = rexp(2*n,rate=attrition_rate),
    sigma_random_intercept = varPar[i, 'sigma_random_intercept'],
    sigma_random_slope = varPar[i, 'sigma_random_slope'],
    sigma_residual = varPar[i, 'sigma_residual'],
    ran.intercept = rnorm(2*n, sd=varPar[i, 'sigma_random_intercept']),
    ran.slope     = rnorm(2*n, sd=varPar[i, 'sigma_random_slope'])) %>%
    mutate(age_c = age - mean(age))

  trial <- right_join(subjects, 
    expand.grid(id = 1:(2*n), month=months)) %>%
    mutate(
      residual = rnorm(2*n*length(months), sd=varPar[i, 'sigma_residual']),
      group = factor(active, 0:1, c('placebo', 'active')),
      missing = ifelse(month>censor, 1, 0)) %>%
    arrange(id, month) %>%
    filter(!missing)
  trial$ADAS11 <- round(
    model.matrix(~ female+age_c+month+month:active, data = trial)[, names(Beta)] %*% 
    Beta +
    with(trial, ran.intercept + ran.slope*month + residual), 
    digits = 0
  )[,1]
  trial
}))


## ---------------------------------------------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_random_intercept==2 & sigma_random_slope==0.2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Residual SD =', sigma_residual)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))


## ---------------------------------------------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_slope==0.2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Random intercept SD =', sigma_random_intercept)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.6,0.8))


## ---------------------------------------------------------------------------------------------------------------
ggplot(filter(varPlot, sigma_residual==2 & sigma_random_intercept==2), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~paste('Random slope SD =', sigma_random_slope)) +
  ylim(0,70) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))


## ----trial_lme_apoe_int, size = 'tiny'--------------------------------------------------------------------------
fit_lme_int <- update(fit_lme, random = ~1|id)
summary(fit_lme_int)


## ----trial_lme_apoe_int_vs_slope, size = 'footnotesize', echo = TRUE--------------------------------------------
anova(fit_lme_int, fit_lme)


## ---------------------------------------------------------------------------------------------------------------
means <- expand.grid(female=1, age_c=0, month=0:18, active=0, id=1:3) %>%
  filter(id %in% c(1,2) | month %in% months) 
means <- mutate(means,
    ADAS11 = (model.matrix(~ female+age_c+month+month:active, means)[, names(Beta)] %*% Beta)[,1],
    ADAS11 = replace(ADAS11, ADAS11<0, 0),
    ADAS11 = replace(ADAS11, ADAS11>70, 70),
    ADAS11 = replace(ADAS11, id == 2, ADAS11 + month*month*0.1),
    ADAS11 = replace(ADAS11, id == 3 & month == 0, ADAS11 + 0),
    ADAS11 = replace(ADAS11, id == 3 & month == 6, ADAS11 + 10),
    ADAS11 = replace(ADAS11, id == 3 & month == 12, ADAS11 + 2),
    ADAS11 = replace(ADAS11, id == 3 & month == 18, ADAS11 + 25)) %>%
  filter(id %in% 1:3) %>%
  mutate(Mean = factor(id, labels=c('linear', 'quadratic', 'categorical')))

ggplot(means, aes(x=month, y=ADAS11, group=Mean, color=Mean)) + 
  geom_line() +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.15,0.7))


## ----varPar2, echo=FALSE----------------------------------------------------------------------------------------
# simulate data with different variance parameters
varPar <- expand.grid(
  variance = c('homogeneous', 'heterogeneous'),
  correlation = c('uncorrelated', 'correlated')
)

SD <- list(homogeneous = sqrt(rep(4, 4)), heterogeneous = (1:4)*2)
Cor <- list(uncorrelated = 0, correlated = 0.8)

varPlot <- do.call(rbind, lapply(1:nrow(varPar), function(i){
  set.seed(20170714)
  subjects <- data.frame(
    id = 1:(2*n),
    active = sample(c(rep(0,n), rep(1,n)), 2*n),
    female = sample(0:1, 2*n, replace=TRUE),
    age_c = rnorm(2*n, 0, 7.8),
    censor = rexp(2*n,rate=attrition_rate))
    
  vv <- diag(SD[[varPar[i,'variance']]])
  cc <- matrix(Cor[[varPar[i,'correlation']]], nrow=4, ncol=4)
  diag(cc) <- 1
  resids <- as.numeric(t(rmvnorm(nrow(subjects), mean=rep(0,4), sigma=vv%*%cc%*%vv)))

  trial <- right_join(subjects, 
    expand.grid(id = 1:(2*n), month=months)) %>%
    arrange(id, month) %>%
    mutate(residual = resids,
      group = factor(active, 0:1, c('placebo', 'active')),
      missing = ifelse(month>censor, 1, 0),
      variance = varPar[i,'variance'],
      correlation = varPar[i,'correlation']) %>%
    arrange(id, month) %>%
    filter(!missing)
  trial$ADAS11 <- round(
    model.matrix(~ female+age_c+month+month:active, data = trial)[, names(Beta)] %*% 
    Beta + trial$residual, 
    digits = 0
  )[,1]
  trial
}))


## ---------------------------------------------------------------------------------------------------------------
ggplot(filter(varPlot, correlation=='correlated'), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~variance) +
  ylim(0,50) +
  scale_x_continuous(breaks=months) +
  theme(legend.position=c(0.1,0.8))


## ---------------------------------------------------------------------------------------------------------------
ggplot(filter(varPlot, variance=='heterogeneous'), 
  aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'lm', size = 2) +
  facet_wrap(~correlation) +
  ylim(0,50) +
  scale_x_continuous(breaks=months) +
  theme(legend.position='none')


## ----echo=FALSE, eval=FALSE-------------------------------------------------------------------------------------
####################################################
### pilot estimates are from a model fit to ADNI ###
####################################################

# library(ADNIMERGE) # available for loni.usc.edu
# adni_ad <- filter(adnimerge, M<=24 & M!=18 & !is.na(ADAS11) & DX.bl=='AD') %>%
#   mutate(m = as.factor(M),
#      visNo = as.numeric(m))
#
# with(adni_ad, table(m, visNo))
# fit_adni <- gls(ADAS11 ~ PTGENDER + center(AGE) + m, data=adni_ad,
#   correlation = corSymm(form = ~ visNo | RID),
#   weights = varIdent(form = ~ 1 | m) )
# summary(fit_adni)


## ----echo=TRUE--------------------------------------------------------------------------------------------------
Beta <- c(
   '(Intercept)'= 19.8, # mean ADAS at baseline
        'female'=-0.51, # female perform better
         'age_c'= 0.04, # worse change for older at baseline (age mean centered)
            'm6'= 2.23, # worsening at month 6 in pbo
           'm12'= 4.46, # worsening at month 12 in pbo
           'm18'= 7.31, # worsening at month 18 in pbo
     'm6:active'=-0.20, # relative improvement at month 6 with treatment
    'm12:active'=-0.70, # relative improvement at month 12 with treatment
    'm18:active'=-1.75) # relative improvement at month 18 with treatment

# other design parameters
months <- c(0, 6, 12, 18)
n <- 200 # per group
attrition_rate <- 0.40/18 # approx per month

# var-cov parameters
SD <- 6.77                          # standard deviation scale parameter
vv <- diag(c(1, 1.2, 1.5, 1.8))          # heterogeneous variance weight matrix
cc <- matrix(0.75, nrow=4, ncol=4)   # correlation matrix
diag(cc) <- 1


## ---------------------------------------------------------------------------------------------------------------
# set seed so that simulation is reproducible
set.seed(20170714)

# simulate subject specific data
subjects <- data.frame(
  id = 1:(2*n),
  active = sample(c(rep(0,n), rep(1,n)), 2*n),
  female = sample(0:1, 2*n, replace=TRUE),
  age = rnorm(2*n, 75, 7.8), 
  censor = rexp(2*n,rate=attrition_rate)) %>%
  mutate(age_c = age-mean(age))
  
# simulate vector of correlated residuals
resids <- as.numeric(t(rmvnorm(nrow(subjects), mean=rep(0,nrow(vv)), sigma=SD^2*vv%*%cc%*%vv)))

# simulate data over time
trial <- right_join(subjects,
  expand.grid(id = 1:(2*n), month=months)) %>%
  arrange(id, month) %>%    ## WARNING: data must be properly sorted by subject and time 
  mutate(residual = resids, ## prior to appending residuals
    group = factor(active, 0:1, c('placebo', 'active')),
    missing = ifelse(month>censor, 1, 0),
    m = as.factor(month),
    visNo = as.numeric(m)) %>%
  arrange(id, month)

# create visit indicators
trial <- cbind(trial, model.matrix(~ -1+m, data = trial))

# calculate the ADAS scores with random effects and residuals and 
# round to nearest digit in 0-70
trial <- mutate(trial,
  ADAS11 = (model.matrix(~ female+age_c+m6+m12+m18+(m6+m12+m18):active, data = trial)[, names(Beta)] %*% Beta)[,1],
  ADAS11 = round(ADAS11 + residual, digits = 0),
  ADAS11 = replace(ADAS11, ADAS11<0, 0),
  ADAS11 = replace(ADAS11, ADAS11>70, 70))

# filter out the missing observations
trial_obs <- filter(trial, !missing)

# transfrom data from long to wide
trial_wide <- trial_obs %>%
  select(id, month, female, age, age_c, active, group, ADAS11) %>% 
  mutate(month = paste0('ADAS11.m', month)) %>%
  spread(month, ADAS11) %>%
  select(id:group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

# data for MMRM
trial_mmrm <- right_join(
  select(trial_wide, id, ADAS11.m0), 
  filter(trial_obs, month>0)) %>%
  mutate(ADAS11.ch = ADAS11 - ADAS11.m0,
    m = as.factor(month),
    visNo = as.numeric(m))


## ---------------------------------------------------------------------------------------------------------------
trial_mmrm %>% select(-group, -missing, -month, -visNo) %>% 
  select(id, m, m0, m6, m12, m18, everything()) %>%
  head() %>%
  kable()


## ----spaghetti_plot---------------------------------------------------------------------------------------------
ggplot(trial, aes(x=month, y=ADAS11, group=id, color=group)) + 
  geom_line(alpha=0.25) +
  geom_smooth(aes(group = NULL), method = 'loess', size = 2) +
  scale_x_continuous(breaks=months, lim=c(0,18)) +
  theme(legend.position=c(0.1, 0.85), legend.background = element_rect(fill=NA))


## ----echo = TRUE------------------------------------------------------------------------------------------------
# Symmetric correlation, heterogeneous variance
MMRMsymHet <- gls(ADAS11.ch ~ 
  -1+ADAS11.m0+female+age_c+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_mmrm, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# Compound Symmetric correlation, heterogeneous variance
MMRMcompSymHet <- gls(ADAS11.ch ~ 
  -1+ADAS11.m0+female+age_c+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_mmrm, correlation = corCompSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )


## ----echo=TRUE, eval=FALSE--------------------------------------------------------------------------------------
summary(MMRMsymHet)


## ---------------------------------------------------------------------------------------------------------------
x <- summary(MMRMsymHet)
print(summary(x$modelStruct), sigma = x$sigma)


## ---------------------------------------------------------------------------------------------------------------
cat('Coefficients:\n')
xtTab <- as.data.frame(x$tTable)
printCoefmat(x$tTable, eps=0.001, digits=3)
# cat("\nStandardized residuals:\n")
# print(x$residuals)
cat("\n")
cat("Residual standard error:", format(x$sigma),"\n")
# cat("Degrees of freedom:", dd[["N"]],"total;",dd[["N"]] - dd[["p"]],
#     "residual\n")


## ----echo = TRUE------------------------------------------------------------------------------------------------
# Symmetric correlation, heterogeneous variance
cLDAsymHet <- gls(ADAS11 ~ 
  -1+female+age_c+m0+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_obs, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )

# Compound Symmetric correlation, heterogeneous variance
cLDAcompSymHet <- gls(ADAS11 ~ 
  -1+female+age_c+m0+(m6+m12+m18)+(m6+m12+m18):active,
  data=trial_obs, correlation = corCompSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m) )


## ---------------------------------------------------------------------------------------------------------------
x <- summary(cLDAsymHet)
print(summary(x$modelStruct), sigma = x$sigma)


## ---------------------------------------------------------------------------------------------------------------
cat('Coefficients:\n')
xtTab <- as.data.frame(x$tTable)
printCoefmat(x$tTable, eps=0.001, digits=3)
# cat("\nStandardized residuals:\n")
# print(x$residuals)
cat("\n")
cat("Residual standard error:", format(x$sigma),"\n")
# cat("Degrees of freedom:", dd[["N"]],"total;",dd[["N"]] - dd[["p"]],
#     "residual\n")


## ----echo = TRUE------------------------------------------------------------------------------------------------
# Linear time
cLDAlin <- gls(ADAS11 ~ 
  female + age_c + month + month:active,
  data=trial_obs, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m))

# Quadratic time
cLDAquad <- gls(ADAS11 ~ 
  female + age_c + (month + I(month^2)) + (month + I(month^2)):active,
  data=trial_obs, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m))


## ----echo = FALSE-----------------------------------------------------------------------------------------------
# Natural cubic spline
b1 <- function(t){
  as.numeric(predict(splines::ns(trial_obs$month, df=2), t)[,1])
}
b2 <- function(t){
  as.numeric(predict(splines::ns(trial_obs$month, df=2), t)[,2])
}
cLDAncs <- gls(ADAS11 ~ 
  female + age_c + (b1(month) + b2(month)) + (b1(month) + b2(month)):active,
  data=trial_obs, correlation = corSymm(form = ~ visNo | id),
  weights = varIdent(form = ~ 1 | m))


## ----echo = FALSE, size = 'scriptsize'--------------------------------------------------------------------------
plotData0 <- filter(trial_obs, !duplicated(paste(month, active))) %>%
  arrange(active, month) %>%
  mutate(female = 1, age_c=0) %>%
  dplyr::select(-age, -censor, -residual, -ADAS11)

plotMatrix_cat <- model.matrix(~ -1 + female + age_c + m0 + (m6 + m12 + m18) + 
  (m6 + m12 + m18):active, 
  data = plotData0)

plotMatrix_quad <- model.matrix(~ female + age_c + (month + I(month^2)) + 
  (month + I(month^2)):active, 
  data = plotData0)

plotMatrix_ncs <- model.matrix(~ female + age_c + (b1(month) + b2(month)) + 
  (b1(month) + b2(month)):active,
  data = plotData0)

plotData <- bind_rows(
  bind_cols(plotData0, confint(glht(cLDAsymHet, linfct = plotMatrix_cat))$confint %>%
      as.data.frame()) %>%
    mutate(model = 'categorical time'),
  bind_cols(plotData0, confint(glht(cLDAquad, linfct = plotMatrix_quad))$confint %>%
      as.data.frame()) %>%
    mutate(model = 'quadratic time'),
  bind_cols(plotData0, confint(glht(cLDAncs, linfct = plotMatrix_ncs))$confint %>%
      as.data.frame()) %>%
    mutate(model = 'natural cubic splines'))


## ---------------------------------------------------------------------------------------------------------------
summaryTable <- trial_obs %>% 
  group_by(group, month) %>%
  summarise(
    n=length(ADAS11),
    mean=mean(ADAS11),
    sd=sd(ADAS11),
    lower95 = smean.cl.normal(ADAS11)[['Lower']],
    upper95 = smean.cl.normal(ADAS11)[['Upper']],
    min=min(ADAS11),
    max=max(ADAS11))
countTab <- ggplot(summaryTable, aes(x=month, y=group, label=n)) + geom_text() + theme_table()

ggplot(plotData, aes(x = month, y = Estimate))+
  geom_line(aes(color=group, linetype=model), position=position_dodge(width=2)) +
  geom_point(aes(color=group, shape=model), position=position_dodge(width=2)) +
  ylim(c(15,30)) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position='top')


## ---------------------------------------------------------------------------------------------------------------
contrastData0 <- filter(trial_obs, !duplicated(paste(month, active))) %>%
  arrange(active, month) %>%
  filter(active==1 & month>0) %>%
  mutate(female = 1, age_c=0) %>%
  dplyr::select(-age, -censor, -residual, -ADAS11)

contrastMatrix_cat <- as.data.frame(model.matrix(~ -1 + female + age_c + m0 + (m6 + m12 + m18) + (m6 + m12 + m18):active, 
  data = contrastData0)) %>%
  mutate(female=0, m0=0, m6=0, m12=0, m18=0)

contrastMatrix_quad <- as.data.frame(model.matrix(~ female + age_c + (month + I(month^2)) + (month + I(month^2)):active, 
  data = contrastData0)) %>%
  mutate('(Intercept)' = 0, female=0, month=0, 'I(month^2)' = 0)

contrastMatrix_ncs <- as.data.frame(model.matrix(~ female + age_c + (b1(month) + b2(month)) + (b1(month) + b2(month)):active, 
  data = contrastData0)) %>%
  mutate('(Intercept)' = 0, female=0, 'b1(month)'=0, 'b2(month)' = 0)

rownames(contrastMatrix_quad) <- rownames(contrastMatrix_cat) <- 
  rownames(contrastMatrix_ncs) <-paste0('m', months[-1])


## ---------------------------------------------------------------------------------------------------------------
plotData %>%
  filter(model == 'natural cubic splines') %>%
  ggplot(aes(x = month, y = Estimate))+
  geom_line(aes(color=group)) +
  geom_point(aes(color=group)) +
  ylim(c(15,30)) +
  scale_x_continuous(breaks=months) +
  ylab('Mean ADAS (95% CI)') +
  theme(legend.position='top')


## ----eval=FALSE, echo=TRUE--------------------------------------------------------------------------------------
gls(ADAS11 ~  female + age_c +   # bl covs
    ns(months, df=2) +        # natural spline for placebo (1 knot)
    ns(months, df=2):active + # natural spline for active (1 knot)
  correlation = corSymm(form = ~ visNo | id), # general correlation
  weights = varIdent(form = ~ 1 | visNo))     # heterogeneous variance


## ---------------------------------------------------------------------------------------------------------------
cts <- bind_rows(
  (glht(cLDAsymHet, linfct = as.matrix(contrastMatrix_cat)) %>% confint())$confint %>%
    as.data.frame() %>%
    mutate(model = 'categorical time'),
  (glht(cLDAquad, linfct = as.matrix(contrastMatrix_quad)) %>% confint())$confint %>%
    as.data.frame() %>%
    mutate(model = 'quadratic time'),
  (glht(cLDAncs, linfct = as.matrix(contrastMatrix_ncs)) %>% confint())$confint %>%
    as.data.frame() %>%
    mutate(model = 'natural cubic splines')) %>%
  mutate(month = rep(months[-1], 3))

ggplot(cts, aes(x=month, y=Estimate, color=model)) +
  geom_point(aes(shape=model), position=position_dodge(width=2)) +
  geom_errorbar(aes(x=month, ymax=upr, ymin=lwr), position=position_dodge(width=2))


## ---------------------------------------------------------------------------------------------------------------
summary(glht(cLDAsymHet, linfct = as.matrix(contrastMatrix_cat)))




## ---- fig.height=5, fig.width=5*2.2-----------------------------------------------------------------------------
options(digits=5)

hippvol <- ADNIMERGE::adnimerge %>%
  select(RID, COLPROT, VISCODE, Month.bl, SITE, DX.bl, AGE, PTGENDER, ICV, 
    FLDSTRENG, FSVERSION, Hippocampus) %>%
  filter(!is.na(Hippocampus)) %>%
  arrange(RID, Month.bl) %>%
  mutate(yrs = Month.bl/12) %>%
  group_by(RID) %>%
  fill(ICV, .direction = "downup") %>%
  mutate(ICV.bl = first(ICV))

ggplot(hippvol, aes(x=yrs, y=Hippocampus, color=SITE, group=RID)) +
  geom_line(alpha=0.5) +
  facet_wrap(.~DX.bl, scales = 'free_x') +
  scale_colour_hue()


## ---- echo=TRUE-------------------------------------------------------------------------------------------------
hipp_fit <- lme(Hippocampus ~ yrs*DX.bl + AGE + PTGENDER + ICV.bl,
  hippvol, random = ~yrs|RID)

hipp_fit_site <- lme(Hippocampus ~ yrs*DX.bl + AGE + PTGENDER + ICV.bl,
  hippvol, random = list(SITE = ~ yrs, RID = ~ yrs))

anova(hipp_fit_site, hipp_fit)


## ----size='tiny'------------------------------------------------------------------------------------------------
summary(hipp_fit_site)






## ----results = 'hide', echo = FALSE-----------------------------------------------------------------------------
# get default predictor matrix
ini_mi <- mice(trial_wide, maxit = 0, print = FALSE)
predictorMatrix <- ini_mi$predictorMatrix
# Don't want ADAS11.m12 predict by ADAS11.m18, etc.:
predictorMatrix['ADAS11.m12', 'ADAS11.m18'] <- 0
predictorMatrix['ADAS11.m6', 'ADAS11.m12'] <- 0
predictorMatrix['ADAS11.m6', 'ADAS11.m18'] <- 0


## ----results = 'hide', echo = TRUE------------------------------------------------------------------------------
trial_imp <- mice(trial_wide, predictorMatrix=predictorMatrix, seed = 20170714, maxit=100)


## ----echo = TRUE, size = 'scriptsize'---------------------------------------------------------------------------
print(head(trial_wide), digits = 2) # raw data with missing values:
print(head(complete(trial_imp)), digits = 2) # first complete version:


## ----echo = FALSE, size = 'scriptsize'--------------------------------------------------------------------------
fits_mi <- with(data=trial_imp, lm(ADAS11.m18~active*center(ADAS11.m0)))
summary(fits_mi)


## ----echo = FALSE-----------------------------------------------------------------------------------------------
summary(pool(fits_mi)) %>%
  remove_rownames() %>%
  column_to_rownames(var="term") %>%
  printCoefmat()


## ---------------------------------------------------------------------------------------------------------------
post <- trial_imp$post
k_tipping <- seq(0.5, 2, 0.25)
est_tipping <- lapply(k_tipping, function(tip){
  # increase imputed ADAS11.m18 in the active group by 
  # factor of k x MAR treatment estimate (3.1)
  # (nullify the imputed treatment effect to varying degrees)
  post["ADAS11.m18"] <- paste0("
  imp[[j]][data$active[!r[, j]] == 1, i] <- imp[[j]][data$active[!r[, j]] == 1, i] + ", 
  tip, 
  " * 4
  imp[[j]][data$active[!r[, j]] == 0, i] <- imp[[j]][data$active[!r[, j]] == 0, i]
  ")
  imp_k <- mice(trial_wide, post=post, predictorMatrix=predictorMatrix, seed = 20170714, maxit=100, print = FALSE)
  fit_k <- with(imp_k, lm(ADAS11.m18~active*center(ADAS11.m0)))
  bind_cols(tipping_factor = tip, summary(pool(fit_k))) %>%
    filter(term == 'active') %>%
    select(-term)
})

# inspect 4th case
tip_4 <- k_tipping[4]
post["ADAS11.m18"] <- paste0("
  imp[[j]][data$active[!r[, j]] == 1, i] <- imp[[j]][data$active[!r[, j]] == 1, i] + ", 
  tip_4, " * 4
  imp[[j]][data$active[!r[, j]] == 0, i] <- imp[[j]][data$active[!r[, j]] == 0, i]
")
imp_4 <- mice(trial_wide, post=post, predictorMatrix=predictorMatrix, seed = 20170714, maxit=100, print = FALSE)
fit_4 <- with(imp_4, lm(ADAS11.m18~active*center(ADAS11.m0)))


## ----echo=FALSE-------------------------------------------------------------------------------------------------
bind_rows(est_tipping) %>%
  printCoefmat(digit=4, eps=0.001)


## ---------------------------------------------------------------------------------------------------------------
imputed_ids <- subset(trial_wide, is.na(ADAS11.m18))[1:5,]$id
# filter(trial_wide, id %in% imputed_ids)


## ---------------------------------------------------------------------------------------------------------------
# Imputation assuming under MAR 
filter(complete(trial_imp), id %in% imputed_ids) %>%
  select(id, female, age, group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)


## ---------------------------------------------------------------------------------------------------------------
# Imputation under MNAR (k=0.5)
filter(mice::complete(imp_4), id %in% imputed_ids) %>%
  select(id, female, age, group, ADAS11.m0, ADAS11.m6, ADAS11.m12, ADAS11.m18)

