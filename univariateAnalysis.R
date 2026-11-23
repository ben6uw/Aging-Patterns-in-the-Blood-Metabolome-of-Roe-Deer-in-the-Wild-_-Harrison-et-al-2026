
library(tidyverse)
library(data.table)
library(ggpubr)
library(GGally)
library(lmerTest)
library(pheatmap)
library(pals)
library(survival)
library(doParallel)
library(missForest)
library(caret)
library(RMTstat) # updated TW test
library(mediation)
library(segmented)

# analysis of roe deer data

# Jeff initially sent data for all >1yr old deer, and then sent the data for all deer that were sampled:  This is a message that accompanied the second send:
# Please find attached the updated dataset included the full set of samples (with samples from 1 year old individuals). As you will see, there is only life history / physiological data corresponding to 209 individuals, not 210.

# The reason for this is that the individual F338 has been captured twice in winter 2017 with blood samples (samples 2017-51-020 and 2017-51-048) collected each time (this normally does not happen). When sending the samples for the quantification of the metabolome, we made a mistake and we sent the blood samples corresponding to these two captures (they correspond to your samples 103 and 130). 

# As the best is probably to keep one sample per individual per year, we should remove one of these two samples. I thus only kept data for the sample collected during the second capture of the individuals (which corresponds to your sample 103). I hope it makes sense.


dir('~/Google Drive/')
getwd()
setwd('~/Google Drive/My Drive/Documents/RoeDeer')


rm(list=ls())
################################################
load('data/dataforUnivariateModels')


## use model selection to identify the most substantial sources of variation
stepLmList <- list()
segModList <- list()
p <- numeric() 

for(i in 1:length(mzs)) {
  tmp <- dplyr::select(d, c(mzs[i], age, body_mass, sex, stage, AM_PM, cohort_quality, field_year, date))
  colnames(tmp)[1] <- 'Y'
  stepLmList[[i]] <- step(lm(Y ~ (age + body_mass + sex)^2 + stage +  cohort_quality + field_year + date, tmp), 
                          direction = "both",
                          scope = list(lower = ~ age))
  
  segModList[[i]] <- tryCatch(
    segmented(stepLmList[[i]], seg.Z = ~age),
    error = function(e) {
      message("Skipping item ", i, " because of error: ", e$message)
      return(NA) # Assign NA to the results list and move on
    })
  p[[i]] <- tryCatch(
    davies.test(segModList[[i]], k=10)$p.value,
    error = function(e) {
      message("Skipping item ", i, " because of error: ", e$message)
      return(NA) # Assign NA to the results list and move on
    })
}

## screen models for: 
# those that either did not have a significant segmented regression, by P value:
segEval <- data.frame('mz'=mzs, 'davies_P'=p)
nonsegMZs <- which(segEval$davies_P>0.05)

## or whose segmentation was not identified (failure to identify breakpoint)
failedMZs <- which(sapply(segModList, function(x) is.na(x[1]))) # identifies any failed PC models

modelsToReplace <- union(nonsegMZs, failedMZs)
modelsToReplace <- modelsToReplace[order(modelsToReplace)]
modelsToReplace 
segModsToKeep <- c(1:length(mzs))[!c(1:length(mzs)) %in% modelsToReplace] # no models a significant segmented age term

## for any segmented age model that failed the above criteria: 
# repeat the stepwise linear regression, but allow the model to drop age, which is no longer evaluated for segmentation:

for(i in modelsToReplace) {
  tmp <- dplyr::select(d, c(mzs[i], age, body_mass, sex, stage, cohort_quality, field_year, date))
  colnames(tmp)[1] <- 'Y'
  segModList[[i]] <- step(lm(Y ~ (age + body_mass + sex)^2 + stage + cohort_quality + field_year + date, tmp), direction = "both") 
  }

names(segModList) <- mzs


save(segModList, d, mzs, segModsToKeep, file='univaraiteSegModList')



rm(list=ls())
#### incorporate random effect, so that we can avoid psudo replication and/or properly fit variation in the midst of potential wihin-deer main effects
load('univaraiteSegModList')


## incorporate random effect of deer
# fit the mixed model and test terms
mmList <- list()
singles <- character()

for(i in 1:length(mzs)) {
  tmp <- dplyr::select(d, c(mzs[i], age, body_mass, sex, stage, cohort_quality, field_year, deer, date))
  colnames(tmp)[1] <- 'Y'
  tryCatch( expr = {
    f0 <- formula(segModList[[i]]) # extract selected model 
    f1 <- update(f0, ~ . + (1 | deer)) # add random effect  
    mmList[[i]] <- lmer(f1, data=tmp) # fit mixed model
    singles[i] <- isSingular(mmList[[i]])
  },
  error = function(e) {
    return(NA)
  })
}

names(mmList) <- mzs

## if present, the loop above will fail when incorporating random effects into a segmented model, manually retain those as their segmented form
mmList[[segModsToKeep]] <- segModList[[segModsToKeep]]

aList <- lapply(mmList, anova) # ANOVA 
models <- sapply(mmList, function(x) as.character(formula(x))[3])

# Get marginal and conditional R2 (performance package)
library(performance)
r2s <- t(sapply(mmList, r2)) # Nakagawa's R2
head(r2s)
r2s <- bind_cols('R2_marginal'=unlist(r2s[ ,2]), 'R2_conditional'=unlist(r2s[ ,1]))
r2s$ICC <- r2s$R2_conditional - r2s$R2_marginal ## The contribution of random effects: In simple cases, the ICC corresponds to the difference between the conditional R2 and the marginal R2
head(r2s)

models <- sapply(mmList, function(x) as.character(formula(x))[3])

## add coefficients to the models?
models

coefs <- sapply(mmList, function(x) round(coef(summary(x))[ ,1], 3)) 
coefs <- bind_rows(coefs)
coefs # 'verbose', but complete

MZmodels <- data.frame(mz= mzs, model= models, r2s[ ,c(1,3)], 'singular'=singles)
rownames(MZmodels) = NULL
head(MZmodels)

table('singular'=MZmodels$singular, 'no ICC'=is.na(MZmodels$ICC)) # all singulars have no ICC (expect7ed), and vice versa, ...omit singular

MZmodels <- dplyr::select(MZmodels, -singular)
MZmodels <- bind_cols(MZmodels, coefs)
head(MZmodels)

# for the models that did not convege on a random effect of deer, remvoe the random effect 
MZmodels$model[which(is.na(MZmodels$ICC))] <- gsub('\\(1|deer)', '', MZmodels$model[which(is.na(MZmodels$ICC))]) 

MZmodels$model <- gsub('deer', 'ID', MZmodels$model) # rename

MZmodels$model[which(MZmodels$singular==T)] <- gsub('\\|', '', MZmodels$model[which(MZmodels$singular==T)]) 
MZmodels$model[which(MZmodels$singular==T)] <- gsub('\\+ \\+', '+ ', MZmodels$model[which(MZmodels$singular==T)]) 
## some manual repair will be needed

write.csv(MZmodels, quote=F, row.names = F, file='paper/Table S2 mz model selection.csv')



## get significant terms, and table of metabolites by term:
mmList

sList <- lapply(mmList, summary)
sList[[1]]

getCoefs <- function(x) bind_cols('term'=rownames(x$coefficients), x$coefficients)
m <- lapply(sList, getCoefs)
m <- lapply(m, function(x) {
  x$FDR <- p.adjust(x$`Pr(>|t|)`, 'fdr', n=length(mzs)) 
  return(x) })
names(m) <- mzs

o <- as.data.frame(bind_rows(m, .id = "metabolite"))
head(o)
o <- o[o$term!= '(Intercept)', ]
rownames(o) <- NULL
head(o)

rev(sort(table(o$term[o$FDR <=0.05])))


# manualy modifying the model table: make bold the significant terms
sub <- o[o$FDR<=0.05, ]
view(dplyr::select(sub, c(metabolite, term))) # only the terms sigificant for each pc 

o$metabolite[o$term == 'age:body_mass' & o$FDR <=0.05]
o$metabolite[o$term == 'age' & o$FDR <=0.05]



#################################
# figure panels:
#################################
## volcanoes
head(o)
tmp <- o[o$term %in% c('age', 'age:body_mass'), ]
head(tmp)

tmp$plotlabs <- tmp$metabolite
tmp$plotlabs[tmp$FDR>0.05] <- ''

tmp$term[tmp$term=="age:body_mass"] <- "age:body mass"

mzs
tmp$plotlabs[tmp$plotlabs=="S-adenosylhomocysteine"] <- 'SAH'
 
volc <- print(ggplot(tmp, aes(y=-log(FDR, 10), x=Estimate, label=plotlabs))+
  geom_hline(yintercept = -log(0.05, 10), color='red')+
  geom_point(size=1)+
  theme_bw(base_size = 14)+
  facet_wrap(~term, scales='free_x')+
  ggrepel::geom_text_repel(color='purple', size=4, max.overlaps = 3, 
                           nudge_y = 0.6)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14))+
  labs(x=expression(beta),
       y=bquote("-log"["10"]~"FDR")))



## example plots
head(o)

# age mzs:
a <- o[o$term == 'age' & o$FDR <=0.05, ]
a$metabolite[order(a$FDR)]

a1 <- ggplot(d, aes(y=`lysine`, x=age))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
  geom_point(size=0.5)+
  theme_bw(base_size = 14)+
  scale_x_sqrt(breaks=c(0,1,2,4,6,8,10,14))+
  labs(x='age (years)')

a2 <- ggplot(d, aes(y=`anserine`, x=age))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
  geom_point(size=0.5)+
  theme_bw(base_size = 14)+
  scale_x_sqrt(breaks=c(0,1,2,4,6,8,10,14))+
  labs(x='age (years)')

a3 <- ggplot(d, aes(y=`citrulline`, x=age))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
  geom_point(size=0.5)+
  theme_bw(base_size = 14)+
  scale_x_sqrt(breaks=c(0,1,2,4,6,8,10,14))+
  labs(x='age (years)')


ggarrange(volc, ggarrange(a1, a2, a3, ncol=1), widths=c(2,1))


## segmented age mz(s):
##########################################
### there were no metabolites with segmented age effects (that were an improvement over linear age effects).

o[o$metabolite=="N-N-dimethyl-arginine", ] # segmented age and body mass
summary(mmList[["N-N-dimethyl-arginine"]])

brkpt <- mmList[["N-N-dimethyl-arginine"]]$psi[2]

ggplot(d, aes(y=`N-N-dimethyl-arginine`, x=age))+
  geom_smooth(aes(group = age > brkpt ), method = "lm", se = TRUE, color='white') +
  geom_vline(xintercept = brkpt, linetype = "dashed", color = "grey40")+
  geom_point()+
  theme_minimal(base_size = 12)

d$residNDMA <- residuals(lm(`N-N-dimethyl-arginine` ~ body_mass + sex, d))

DMAplot <- print(ggplot(d, aes(y=residNDMA, x=age))+
  geom_smooth(aes(group = age > brkpt ), method = "lm", se = TRUE, color='white') +
  geom_vline(xintercept = brkpt, linetype = "dashed", color = "grey40")+
  geom_point(size=0.5)+
  theme_bw(base_size = 14)+
  labs(y='N-N-dimethylarginine',
       x='age (years)'))
######################################################################



# age by sex, by body mass interaction mzs:
######################################################################
o$metabolite[o$term == 'age:body_mass' & o$FDR <=0.05]
o$metabolite[o$term == 'age:sexM' & o$FDR <=0.05]

summary(mmList[["indole-3-acetamide"]]) # for indole-3-acetate, age and mass are both possitive, the interaction is negative, so, we anticpate less steep age by mass, or less steep mass effect amogn older deer:

table(cut_number(d$body_mass, 3))
d$sizeGrps <- cut_number(d$body_mass, 3)

table(cut(d$age, breaks =  c(0, 1, 3, 6, Inf)))
d$ageGrps <- as.factor(cut(d$age, breaks =  c(0, 1, 3, 6, Inf)))

ggplot(d, aes(y=`indole-3-acetamide`, x=body_mass))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
  geom_point(size=0.5)+
  theme_bw(base_size = 12)+
  facet_wrap(~ageGrps, scales='free_x', nrow=1)

ggplot(d, aes(y=`indole-3-acetamide`, x=age))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
  geom_point(size=0.5)+
  theme_bw(base_size = 12)+
  facet_wrap(~sizeGrps, scales='free_x', nrow=1)


table(cut(d$age, breaks =  c(0, 2, Inf)))
d$ageGrps <- as.factor(cut(d$age, breaks =  c(0, 2, Inf)))
levels(d$ageGrps) <- c('<2yrs old', '>2yrs old')

ageXmass <- print(ggplot(d, aes(y=`indole-3-acetamide`, x=body_mass))+
  geom_smooth(method = "lm", se = TRUE, color='white') +
    geom_point(size=0.5)+
  theme_bw(base_size = 14)+
  facet_wrap(~ageGrps, scales='free_x', nrow=1)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))+
  labs(x='body mass (kg)'))


ggarrange(volc + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)), 
          ggarrange(ggarrange(a1, a2, nrow=1), ageXmass, ncol=1), widths=c(1, 0.8))



save(mmList, d, mzs, file='univaraiteModels')
