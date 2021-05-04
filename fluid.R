## ----setup, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------
# devtools::install_url('https://cran.rstudio.com/src/contrib/Archive/calibFit/calibFit_2.1.0.tar.gz')
# remotes::install_github('atrihub/SRS')
# For ADNIMERGE, go to http://adni.loni.usc.edu/, https://adni.bitbucket.io/

library(knitr)
library(tidyverse)
library(kableExtra)
library(gridExtra)
library(RefManageR)
library(plotly)
library(calibFit)
library(SRS)
library(mixtools)
library(pROC)
library(party)

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






## ----echo=FALSE, fig.align='center', out.width='57%'------------------------------------------------
knitr::include_graphics("./images/atri.png")


## ----echo=FALSE, fig.align='center', out.width='47%'------------------------------------------------
knitr::include_graphics("./images/actc_logo.png")


## ----generate_batch_data----------------------------------------------------------------------------
# simulated data with batch effects
set.seed(20200225)

batch_data <- 
  tibble(
    batch = 1:10,
    Sigma = rgamma(n=10, shape=360/10, scale=10), # variance for each batch
    Mean = rnorm(n=10, mean=850, sd=200)) %>% # Mean for each batch
  group_by(batch) %>%
  nest() %>%
  mutate(Biomarker = map(data, ~ rnorm(n=50, .$Mean, .$Sigma))) %>%
  unnest(Biomarker) %>%
  unnest(data) %>%
  ungroup() %>%
  arrange(batch) %>%
  mutate(
    id = 1:(10*50),
    batch = as.factor(batch),
    Biomarker = ifelse(Biomarker<0, 0, Biomarker))


## ----batch_data_plot--------------------------------------------------------------------------------
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, alpha=0.2)


## ----batch_data_summaries, results='asis', cache=FALSE----------------------------------------------
batch_data_sum <- batch_data %>% group_by(batch) %>%
  summarize(
    N=length(Biomarker),
    Mean=mean(Biomarker), 
    SD=sd(Biomarker),
    `SD/Mean = CV (%)`=SD/Mean*100)
batch_data_sum %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ---- echo=TRUE, results='markup'-------------------------------------------------------------------
anova(lm(Biomarker ~ batch, batch_data))


## ----batch_confounds--------------------------------------------------------------------------------
low_groups <- subset(batch_data_sum, Mean<median(batch_data_sum$Mean))$batch
batch_data <- batch_data %>%
  mutate(Group = ifelse(batch %in% low_groups, 'A', 'B'))
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(aes(color=Group, fill=Group), 
    binaxis='y', stackdir='center', dotsize=0.3, alpha=0.5)


## ----batch_randomized-------------------------------------------------------------------------------
batch_data$Group <- sample(batch_data$Group, size=nrow(batch_data))
ggplot(batch_data, aes(y=Biomarker, x=batch)) +
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(aes(color=Group, fill=Group), 
    binaxis='y', stackdir='center', dotsize=0.3, alpha=0.5)


## ----randomization----------------------------------------------------------------------------------
data(srs_data)
# head(srs_data)

p.func.greedy.if.possible <- function(overallImbalance, treatmentCounts, maxCounts)
{
    cant.go <- treatmentCounts > maxCounts
    if(all(cant.go)) 
      stop("Randomization impossible. Probably need another treatment group.")

    number.of.treatments <- length(overallImbalance)
    k <- which(overallImbalance == min(overallImbalance))
    p.vec <- rep(0, number.of.treatments)
    p.vec[k] <- 1
    p.vec/sum(p.vec)
    p.vec[cant.go] <- 0

    if(all(p.vec == 0)){ # try less greedy
      number.of.treatments <- length(overallImbalance)
      p.star <- 2/3
      k <- which(overallImbalance == min(overallImbalance))
      if (length(k) > 1) {
          k <- sample(k, 1)
      }
      p.vec <- rep((1 - p.star)/(number.of.treatments - 1), number.of.treatments)
      p.vec[k] <- p.star
      p.vec      
      p.vec[cant.go] <- 0
    }
    
    p.vec
}

get.counts <- function(object)
{
  expt <- object@expt
  treatment.names <- expt@treatment.names
  factor.names <- expt@factor.names
  factor.level.names <- expt@factor.level.names
  treatment.names <- expt@treatment.names
  state.matrix <- object@stateTable
  tr.ratios <- object@tr.ratios
  
  tr.assignments <- object@tr.assignments
  tr.assignments$Treatment <- factor(tr.assignments$Treatment, 
    levels = treatment.names)
  tr.assignments$Counts <- factor(tr.assignments$Counts, 
    levels = factor.level.names[[which(factor.names == "Counts")]])
  tr.assignments <- with(tr.assignments, table(Counts, Treatment)) * 
    as.numeric(factor.level.names[[which(factor.names == "Counts")]])
  colSums(tr.assignments)
}  

expt <- ClinicalExperiment(number.of.factors = 3,
  factor.names = c('Counts', 'Group', 'Age'),
  number.of.factor.levels = c(2, 5, 2),
  factor.level.names = 
    list(c(4, 5), 1:5, c('young', 'old')),
  number.of.treatments = 13,
  treatment.names = as.character(1:13))

g.func <- function(imbalances)
{
    factor.weights <- c (1, 100, 1)
    imbalances %*% factor.weights
}

r.obj <- new("cPocockSimonRandomizer", expt, as.integer(20130827), 
  g.func=g.func, p.func = p.func.greedy.if.possible, max.counts = 30)

for(i in 1:nrow(srs_data)){
  r.obj <- randomize(r.obj, as.character(srs_data[i, "ID"]), 
     as.character(srs_data[i, expt@factor.names]))
}


## ---- results='asis'--------------------------------------------------------------------------------
tr.assignments <- r.obj@tr.assignments %>%
  mutate(
    Treatment = factor(Treatment, levels = r.obj@expt@treatment.names),
    `Subject ID` = 1:nrow(r.obj@tr.assignments)
  ) %>%
  rename(
    Plate = Treatment,
    `Num. of samples` = Counts) %>%
  select(`Subject ID`, `Num. of samples`, Group, Age, Plate)

tr.assignments[1:10, ] %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ---------------------------------------------------------------------------------------------------
tab <- with(tr.assignments, table(Plate, Age)) %>%
  as.data.frame() %>%
  pivot_wider(names_from='Age', values_from='Freq') %>%
  t() 
tab <- rbind(tab, `Num. samples`=as.numeric(get.counts(r.obj)))

tab %>%
  kable(digits=2, format = 'html') %>%
  kable_styling(
    bootstrap_options=c('striped', 'condensed'),
    font_size=18, full_width=FALSE)


## ----calibFit_fits, out.width='100%', fig.height=4, fig.width=4*(2)---------------------------------
data(HPLC)
data(ELISA)

linmodel <- lm(Response~Concentration, data=HPLC)
# The predicted response
HPLC$Fitted <- fitted(linmodel)

p1 <- ggplot(HPLC, aes(x=Concentration, y=Response)) +
  geom_point() +
  geom_line(aes(y=Fitted)) +
	xlab("Concentration (ng/ml)") +
	ylab("Response") +
	ggtitle("HPLC with ordinary least squares fit")

fplmodel <- with(ELISA,
  calib.fit(Concentration, Response, type="log.fpl")
)
# The predicted response
ELISA$Fitted <- fplmodel@fitted.values

p2 <- ggplot(ELISA, aes(x=log(Concentration), y=Response)) +
  geom_point() +
  geom_line(aes(y=Fitted)) +
	xlab("log(Concentration (ng/ml))") +
	ylab("Response") +
	ggtitle("ELISA with 4 parameter logistic fit")

grid.arrange(p1,p2,nrow=1)


## ----calibFit_residuals, out.width='100%', fig.height=4, fig.width=4*(2)----------------------------
p1 <- ggplot(HPLC, aes(x=Concentration, y=Response-Fitted)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Concentration (ng/ml)") +
	ylab("Residuals") +
	ggtitle("HPLC with ordinary least squares fit")

p2 <- ggplot(ELISA, aes(x=log(Concentration), y=Response-Fitted)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("log(Concentration (ng/ml))") +
	ylab("Residuals") +
	ggtitle("ELISA with 4 parameter logistic fit")

grid.arrange(p1,p2,nrow=1)


## ---- calib_fit-------------------------------------------------------------------------------------
cal.fpl <- with(ELISA, calib.fit(Concentration,Response,type="log.fpl"))
cal.lin.pom <- with(HPLC, calib.fit(Concentration,Response,type="lin.pom"))
cal.fpl.pom <- with(ELISA, calib.fit(Concentration,Response,type="log.fpl.pom"))

linpom.fit <- cal.lin.pom@fitted.values
fplpom.fit <- cal.fpl.pom@fitted.values

sig.lin <- cal.lin.pom@sigma
sig.fpl <- cal.fpl.pom@sigma

theta.lin <- cal.lin.pom@theta
theta.fpl <- cal.fpl.pom@theta

linpom.res <- cal.lin.pom@residuals*(1/((linpom.fit^theta.lin)*sig.lin))
fplpom.res <- cal.fpl.pom@residuals*(1/((fplpom.fit^theta.fpl)*sig.fpl))


## ----calibFit_pom_residuals, out.width='100%', fig.height=4, fig.width=4*(2)------------------------
p1 <- ggplot(HPLC, aes(x=linpom.fit, y=linpom.res)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Fitted Values (LS-POM)") +
	ylab("Standardized Residuals") +
	ggtitle("HPLC with least squares POM")

p2 <- ggplot(ELISA, aes(x=fplpom.fit, y=fplpom.res)) +
  geom_point() +
  geom_hline(yintercept = 0) +
	xlab("Fitted Values (FPL-POM)") +
	ylab("Standardized Residuals") +
	ggtitle("ELISA with 4 parameter logistic POM")

grid.arrange(p1,p2,nrow=1)


## ----calib_hplc_pom, fig.height=4.5, fig.width=4.5*2------------------------------------------------
par(mfrow=c(1,2))
ciu <- fitted(linmodel) + summary(linmodel)$sigma*qt(.975,linmodel$df)
cil <- fitted(linmodel) - summary(linmodel)$sigma*qt(.975,linmodel$df)

plot(HPLC[, c('Concentration', 'Response')], main = "HPLC data fit without POM",col="blue",pch=16)
lines(HPLC$Concentration,fitted(linmodel),col="lightblue")
lines(HPLC$Concentration,ciu,col="lightblue",lty=2)
lines(HPLC$Concentration,cil,col="grey",lty=2)

plot(cal.lin.pom,print=FALSE,main = "HPLC data fit with POM",xlab = "Concentration", ylab = "Response")


## ----calib_elisa_pom, fig.height=4.5, fig.width=4.5*2-----------------------------------------------
par(mfrow=c(1,2))
#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl,print=FALSE,main = "ELISA fit without POM",xlab = "Concentration", ylab = "Response")

#par(mar = c(3.5,3.5,1.5,1.5))
plot(cal.fpl.pom,print=FALSE,main = "ELISA fit with POM",xlab = "Concentration", ylab = "Response")


## ----calibrated1, fig.height=4.5, fig.width=4.5, out.width='100%'-----------------------------------
calib.lin <- calib(cal.lin.pom, HPLC$Response)
plot(calib.lin, main="HPLC calribated with linear POM")


## ----calibrated2, fig.height=4.5, fig.width=4.5, out.width='100%'-----------------------------------
calib.fpl <- calib(cal.fpl.pom, ELISA$Response)
plot(calib.fpl, main="ELISA calribated with FPL POM")


## ----classification, fig.height=4, fig.width=6------------------------------------------------------
dd <- subset(ADNIMERGE::adnimerge, !is.na(ABETA)) %>%
  arrange(RID, EXAMDATE) %>%
  filter(!duplicated(RID)) %>%
  mutate(
    ABETA = as.numeric(gsub('<', '', gsub('>', '', ABETA))),
    TAU = as.numeric(gsub('<', '', gsub('>', '', TAU))),
    PTAU = as.numeric(gsub('<', '', gsub('>', '', PTAU)))
  )

ggplot(dd, aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#E69F00", "#D55E00")) +
  theme(legend.position = c(0.80,0.75))


## ----classification_no_mci--------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#D55E00"))


## ----ROC_abeta, fig.width=5, fig.height=5-----------------------------------------------------------
roc_abeta <- roc((DX=='Dementia') ~ ABETA, subset(dd, DX!='MCI'))
# ggroc(roc_abeta)
plot(roc_abeta, print.thres="best", print.thres.best.method="youden")


## ----ROC_abeta_tau, fig.width=5, fig.height=5-------------------------------------------------------
roc_tau <- roc((DX=='Dementia') ~ TAU, subset(dd, DX!='MCI'))
roc_tauabeta <- roc((DX=='Dementia') ~ I(TAU/ABETA), subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red',
  print.thres="best", print.thres.best.method="youden")
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta))),
       col=c("blue", "orange", "red"), lwd=2)



## ---- echo=FALSE, eval=FALSE------------------------------------------------------------------------
roc.test(roc_abeta, roc_tau)
roc.test(roc_abeta, roc_tauabeta)


## ----abeta_tau_scatter_youden, fig.height=2.5, fig.width=3*2----------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  geom_abline(intercept = 0, slope=0.394)


## ---------------------------------------------------------------------------------------------------
logistic.fit <- glm((DX=='Dementia') ~ scale(ABETA) + scale(TAU), 
  data = subset(dd, DX!='MCI'), family=binomial)
summary(logistic.fit)$coef %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|z|)` = format.pval(round(`Pr(>|z|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ----logistic_pred_prob-----------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=predict(logistic.fit, type='response'))) +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Pred. prob. of Dementia") +
  geom_abline(intercept = 0, slope=0.394)


## ----ratio_gradient---------------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=TAU/ABETA)) +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Tau/Abeta ratio") +
  geom_line(data=data.frame(ABETA=80:1700, TAU=80:1700), 
            aes(color=TAU/ABETA), linetype='dashed') +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1.5*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=2*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/2*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/3*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/4*(80:1700)), aes(color=TAU/ABETA)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=1/5*(80:1700)), aes(color=TAU/ABETA)) +
  coord_cartesian(xlim=c(200,1700), ylim=c(80,852))


## ----ratio_gradient_logistic------------------------------------------------------------------------
tau_fun <- function(abeta, p){
  # log(p/(1-p)) = y = a + b*abetaz + c*tauz
  # tauz = (y - a - b*abetaz)/c
  abetaz <- (abeta - attr(logistic.fit$model$`scale(ABETA)`, "scaled:center"))/
    attr(logistic.fit$model$`scale(ABETA)`, "scaled:scale")
  tauz <- (log(p/(1-p)) - logistic.fit$coef['(Intercept)']
    - logistic.fit$coef['scale(ABETA)']*abetaz)/
    logistic.fit$coef['scale(TAU)']
  tauz*attr(logistic.fit$model$`scale(TAU)`, "scaled:scale") +
    attr(logistic.fit$model$`scale(TAU)`, "scaled:center")
}
pd <- subset(dd, DX!='MCI') %>%
  mutate(prob = predict(logistic.fit, type='response'))
ggplot(pd, aes(x=ABETA, y=TAU, color=prob)) + 
  geom_point() +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  labs(color="Pred. prob. of Dementia") +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.01), prob=0.01)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.1), prob=0.1)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.25), prob=0.25)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.5), prob=0.5)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.75), prob=0.75)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.9), prob=0.9)) +
  geom_line(data=data.frame(ABETA=80:1700, TAU=tau_fun(80:1700, p=0.99), prob=0.99)) +
  coord_cartesian(xlim=c(200,1700), ylim=c(80,852))


## ----ROC_logistic, fig.width=5, fig.height=5--------------------------------------------------------
roc_logistic.fit <- roc((DX=='Dementia') ~ predict(logistic.fit, type='response'), subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red')
plot(roc_logistic.fit, add=TRUE, col='purple')
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta)), 'Logistic model'),
       col=c("blue", "orange", "red", "purple"), lwd=2)



## ---- echo=FALSE, eval=FALSE------------------------------------------------------------------------
roc.test(roc_abeta, roc_logistic.fit)


## ---------------------------------------------------------------------------------------------------
logistic.fit2 <- glm((DX=='Dementia') ~ scale(ABETA) + scale(TAU) + scale(I(AGE+Years.bl)) + as.factor(APOE4), family=binomial, data = subset(dd, DX!='MCI'))
summary(logistic.fit2)$coef %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var='Coefficient') %>%
  mutate(`Pr(>|z|)` = format.pval(round(`Pr(>|z|)`, digits = 3), eps = 0.001, digits=3)) %>%
  kable()


## ---- eval=FALSE------------------------------------------------------------------------------------
roc_logistic.fit2 <- roc((DX=='Dementia') ~ predict(logistic.fit2, type='response'), subset(dd, DX!='MCI'))
roc.test(roc_abeta, roc_logistic.fit2)
roc.test(roc_tauabeta, roc_logistic.fit2)
auc(roc_logistic.fit2); ci(roc_logistic.fit2)


## ----tree1, fig.height=4.5, fig.width=5*(2)---------------------------------------------------------
tree.fit <- ctree((DX=='Dementia') ~ ABETA + TAU, data = subset(dd, DX!='MCI'), 
  controls = ctree_control(maxdepth=2))
plot(tree.fit)


## ----tree2------------------------------------------------------------------------------------------
ggplot(subset(dd, DX!='MCI'), aes(x=ABETA, y=TAU)) + 
  geom_point(aes(color=DX)) +
  scale_color_manual(values=c("#0072B2", "#D55E00")) +
  geom_vline(xintercept = 886) +
  #geom_hline(yintercept = 210) +
  #geom_hline(yintercept = 440) +
  geom_segment(aes(x = 0, y = 210, xend = 886, yend = 210)) +
  geom_segment(aes(x = 886, y = 440, xend = 1700, yend = 440))


## ----ROC_rf, fig.width=5, fig.height=5--------------------------------------------------------------
tree.fit2 <- ctree((DX=='Dementia') ~ ABETA + TAU, data = subset(dd, DX!='MCI'))
rf.fit <- cforest((DX=='Dementia') ~ ABETA + TAU, data = subset(dd, DX!='MCI'))
roc_tree.fit <- roc((DX=='Dementia') ~ predict(tree.fit2, type='response')[,1], subset(dd, DX!='MCI'))
roc_rf.fit <- roc((DX=='Dementia') ~ predict(rf.fit, type='response')[,1], subset(dd, DX!='MCI'))

plot(roc_abeta, col='blue')
plot(roc_tau, add=TRUE, col='orange')
plot(roc_tauabeta, add=TRUE, col='red')
plot(roc_logistic.fit, add=TRUE, col='purple')
plot(roc_tree.fit, add=TRUE, col='grey')
plot(roc_rf.fit, add=TRUE, col='black')
legend("bottomright", 
  legend=c(expression(paste("A", beta)), "Tau", expression(paste("Tau/A", beta)), 
      'Logistic model', 'Binary Tree', 'Random Forest'),
    col=c("blue", "orange", "red", "purple", "grey", "black"), lwd=2)


## ---- echo=FALSE, eval=FALSE------------------------------------------------------------------------
roc.test(roc_abeta, roc_tree.fit)
roc.test(roc_abeta, roc_rf.fit)


## ----density_Abeta, fig.height=2.25, fig.width=3*(2)------------------------------------------------
ggplot(subset(dd, DX!='MCI' & ABETA<1700), aes(x=ABETA)) + 
  geom_histogram(aes(y=..density..), alpha=0.5) +
  geom_density() +
  geom_rug(aes(color=DX))


## ---- results='hide', cache=TRUE--------------------------------------------------------------------
#' Plot a Mixture Component
#' 
#' @param x Input data
#' @param mu Mean of component
#' @param sigma Standard deviation of component
#' @param lam Mixture weight of component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
mixmdl <- normalmixEM(subset(dd, DX!='MCI' & ABETA<1700)$ABETA, k = 2)
## number of iterations= 121

set.seed(1)
mvmixmdl <- mvnormalmixEM(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], k = 2)
## number of iterations= 85


## ----mixture_distribution_Abeta---------------------------------------------------------------------
data.frame(ABETA = mixmdl$x) %>%
  ggplot(aes(x=ABETA)) +
  geom_histogram(aes(y=..density..), alpha=0.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "#D55E00", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "#0072B2", lwd = 1.5) +
  ylab("Density")


## ---------------------------------------------------------------------------------------------------
post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior)) %>%
  arrange(x) %>% 
  rename(Abeta = x, `Prob. Abnormal` = comp.1, `Prob. Normal` = comp.2)
post.df %>% filter(Abeta > 1030 & Abeta < 1080) %>% kable()

# A reasonable cutoff might be:
# filter(post.df, `Prob. Abnormal` <= `Prob. Normal`)[1,]


## ----Bivariate_Density------------------------------------------------------------------------------
kd <- with(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], MASS::kde2d(ABETA, TAU, n = 50))
fig <- plot_ly(x = kd$x, y = kd$y, z = kd$z) %>% add_surface() %>% layout(
  showlegend = FALSE,
  margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0),
  scene = list(
      xaxis = list(title = "Abeta"),
      yaxis = list(title = "Tau"),
      zaxis = list(title = "Density")
    ))

htmlwidgets::saveWidget(fig, "bvdensity_csf_tau.html")


## ----bv_kernel_density------------------------------------------------------------------------------
kd <- with(subset(dd, DX!='MCI' & ABETA<1700)[, c('ABETA', 'TAU')], MASS::kde2d(ABETA, TAU, n = 50))
kdl <- expand.grid(i=1:50, j=1:50) %>%
  mutate(Abeta=NA, Tau=NA, density=NA)
for(r in 1:nrow(kdl)){
  kdl[r, 'Abeta'] <- kd$x[kdl$i[r]]
  kdl[r, 'Tau'] <- kd$y[kdl$j[r]]
  kdl[r, 'density'] <- kd$z[kdl$i[r], kdl$j[r]]
}

ggplot(kdl, aes(Abeta, Tau, z = density)) +
  geom_raster(aes(fill=density))+
  geom_contour(color='white')


## ----mvmix_post_prob, fig.width=5, fig.height=5, out.width='100%'-----------------------------------
post.df2 <- cbind(mvmixmdl$x, mvmixmdl$posterior) %>% 
  as.data.frame() %>% 
  rename(
    `Prob. Abnormal` = comp.1,
    `Prob. Normal` = comp.2
  )

ggplot(post.df2, aes(x=ABETA, y=TAU, color=`Prob. Abnormal`)) + 
  geom_point() +
  scale_colour_gradient(low="#0072B2", high="#D55E00") +
  geom_abline(intercept = 0, slope=0.394) +
  theme(legend.position=c(0.85, 0.82))


## ----mvmix_density, fig.width=5, fig.height=5, out.width='100%'-------------------------------------
plot(mvmixmdl, density = TRUE, alpha = c(0.01, 0.05, 0.10), 
  marginal = FALSE, whichplots=3, main2='', 
  xlab2='ABETA', ylab2='TAU')


## ----refs, echo=FALSE, results="asis"---------------------------------------------------------------
PrintBibliography(bib)

