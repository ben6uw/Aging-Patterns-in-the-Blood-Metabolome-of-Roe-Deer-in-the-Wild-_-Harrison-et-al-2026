library(ggplot2)
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
library(segmented)
library(RMTstat) # updated TW test

# analysis of roe deer data

# Jeff initially sent data for all >1yr old deer, and then sent the data for all deer that were sampled:  This is a message that accompanied the second send:
# Please find attached the upJulianDated dataset included the full set of samples (with samples from 1 year old individuals). As you will see, there is only life history / physiological data corresponding to 209 individuals, not 210.

# The reason for this is that the individual F338 has been captured twice in winter 2017 with blood samples (samples 2017-51-020 and 2017-51-048) collected each time (this normally does not happen). When sending the samples for the quantification of the metabolome, we made a mistake and we sent the blood samples corresponding to these two captures (they correspond to your samples 103 and 130). 

# As the best is probably to keep one sample per individual per year, we should remove one of these two samples. I thus only kept data for the sample collected during the second capture of the individuals (which corresponds to your sample 103). I hope it makes sense.

rm(list=ls())

dir('~/Google Drive/')
getwd()
setwd('~/Google Drive/My Drive/Documents/RoeDeer')

dir('data') # NOTE: the comprehensive data is tagged '-V2'.  The data in the previous set are all in V2, so just work with V2

mz <- read.csv("data/roeDeerAreaReport3-METABOLOME.csv", head=T, stringsAsFactors = F)
head(mz) # MacCoss lab data
str(mz)
table(rowSums(is.na(mz)))
colSums(is.na(mz))


#################################################################################
# characterize demographic data:
#################################################################################
dir('data')
h <- as.data.frame(readxl::read_xlsx('data/RoeDeer-LIFE HISTORY DATA-V6.xlsx'))
head(h) 
str(h)

table(h$pop) # 'pop' is invariant

head(h)
h$deer <- h$`Id_JMG-UniqueIndivual-ID` # deer id
length(unique(h$deer)) # 147 deer sampled, 84 are >1 year old

table(h$Sex)
h$Sex <- as.factor(h$Sex)
h$`Morning_Evening(Caputre)`
h$`Morning_Evening(Caputre)` <- as.factor(h$`Morning_Evening(Caputre)`)
colnames(h)[colnames(h)=="Morning_Evening(Caputre)"] <- 'AM_PM'
levels(h$AM_PM) <- c('evening', 'morning')
h$AM_PM <- factor(h$AM_PM, levels=c('morning', 'evening'))
colnames(h)[colnames(h)=="CohortQualityIndex-StandardiedJuvenieleMass-Kg"] <- 'CohortQuality'

h$YearOfSampling <- as.factor(h$YearOfSampling)
h$Cohort <- as.factor(h$Cohort)
h$Batch <- as.factor(h$Batch) # metabolome batch

table(h$Sex, h$AM_PM) # time of sampling balanced over sex
table(h$Cohort, h$AM_PM) # AM PM gets streched over cohort
table(h$Cohort, h$Sex) # as does sex

hist(h$Age, 20) 
table(h$Age)
h$stage <- as.factor(ifelse(h$Age <=1, 'juvenile', 'adult'))
h$stage 

p1 <- ggplot(h, aes(y=Age, x=Cohort, color=Sex))+
  geom_jitter(width=0.05, size=1)+
  ylab('age (years)')+
  theme_bw(base_size = 16)+
  ggtitle('age structure by cohort')+
  facet_wrap(~stage)+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=1))
p1

mean <- aggregate(Age ~ deer, h, FUN=function(x) mean(x, na.rm=T))
mean

h$orderedDeer <- factor(h$deer, levels=mean$deer[order(mean$Age)])
h$orderedDeer

p2 <- ggplot(h, aes(x=Age, y=orderedDeer, color=Sex, group=orderedDeer))+
  geom_line()+
  geom_point(size=1)+
  theme_classic(base_size = 16)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  ylab('deer (n=147)')+
  xlab('age (years)')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  ggtitle('some longitudinal sampling')
p2

ggarrange(p1, p2, nrow=1, common.legend = T, legend='right')


head(h)
n <- as.data.frame(table(h$deer))
colnames(n) <- c('deer', 'samples.per.deer')

h <- left_join(h, n) # number of samples per deer
h$samples.per.deer


table(table(h$deer, h$Age)) # should be a maximum of 1, this means that re-sampled deer were sampled at most, one time per year.
table(table(h$deer)) # number of deer samples, 1, 2 or 3 times
table(h$JulianDate) # all sampling occured between mid-January and mid Feb

h$JulianDate # the 'day of the (Jualian) year'; note sampling occured in a window of dates within Feb-Mar 
range(h$JulianDate)

head(h)

ggplot(h, aes(x=Age, y=Body_Mass, color=Cohort, group=deer))+
  geom_line()+
  geom_point(size=1)+
  theme_classic(base_size = 16)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  xlab('age (years)')+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
#  cohort is confounded with age, and thus partly with weight.

ggplot(h, aes(x=sqrt(Age), y=Body_Mass))+
  geom_point(size=1)+
  theme_classic(base_size = 16)+
  geom_smooth(linetype=0)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  facet_wrap(~ YearOfSampling)



## body mass and age distributions across study years, the pop looks well sampled acorss wt and age in each sampling year.

summary(lmer(Body_Mass ~ sqrt(Age) * YearOfSampling + JulianDate + (1|Cohort), h)) # controlling for cohort (which is confounded with age)
summary(lm(Body_Mass ~ sqrt(Age) * YearOfSampling + JulianDate, h))

## look at persistence of effect of juvanile mass (cohort quality) on mass with age:
tmp <- h

ggplot(tmp, aes(x=CohortQuality, y=Body_Mass))+
  geom_smooth(method='lm', alpha=0.5, linetype=0)+
  geom_jitter(width=0.01, size=1)+
  facet_wrap(~stage)+
  theme_bw(base_size = 14)+
  xlab('Cohort Quality (standardized juveniele mass)')+
  ggtitle('any infuence of cohort mass on later deer mass diminished with age\nnote: points are jittered for visibility')


## This population and sampling regime are complex, we can anticipate several issues:
# confounding (with age in particular) of: 
# 1. cohort and age
# 2. possible confounding of CBC measures with age
# 3. body size

# effects from the environmnent might be reflected in:
# 1. year of sampling 
# 2. cohort, including cohort quality
# 3. Julian date
# 4. body size


#####################################################################
## plot sample/study design:
ggpairs(h, columns=c('JulianDate', 'YearOfSampling', 'Run.order', 'Batch'), ggplot2::aes(colour = Batch, alpha=0.5)) + theme_bw() # there it is.  ...the order of sample collection is highly (completely?) colinear with batch and with run order.

# how to handle these confounds?
# the issue(s): effects from any of these variables can appear as effects of the others
# could use a model to test for effects of all in context of the others and look for the largest effect(s).  Include any singificant effects in a technical correction.  Technical in this sense includes both study design, LCMS/prep artifacts (Batch and run order) as well as some environmental variation (that manifesting in year of field sampling and by day-of-the-year).


#################################################################################
# merge metabolome and demographic data:
#################################################################################
head(h) 
table(h$Id_Metabolome, h$Batch)
mz$Sample
h$Id_Metabolome

mz <- left_join(mz, h)

mz$sampleType <- NA
c('sheep', 'pool', 'iroa', 'pb', 'db')
mz$sampleType[grepl('sheep', mz$Sample)] <- 'sheep'
mz$sampleType[grepl('iroa', mz$Sample)] <- 'iroa'
mz$sampleType[grepl('db', mz$Sample)] <- 'db'
mz$sampleType[grepl('pb', mz$Sample)] <- 'pb'
mz$sampleType[grepl('pool', mz$Sample)] <- 'pool'

unique(mz$Sample[is.na(mz$sampleType)]) # should all be study samples
mz$sampleType[is.na(mz$sampleType)] <- 'studySample'
mz$sampleType <- factor(mz$sampleType)
table(mz$sampleType) 

length(unique(h$deer))
x <- dplyr::select(h, c(deer, Sex))
table(x$Sex[!duplicated(x)])



######################################################################
### normalize metabolome data, assess run order and any batch effects
######################################################################
head(mz)
table(mz$Normalization.Method) # mzs were either normalized to global standards, or to an isotope
plot(mz$Normalized.Area, mz$TotalAreaRatioToGlobalStandards, pch=19, log='xy') # visualize the normalization effects

# use Normalized.Area:
ggplot(mz, aes(y=log(Normalized.Area), x=Run.order, color=sampleType)) +
  geom_point(size=0.5)+
  theme_classic(base_size = 7)+
  facet_wrap(~Molecule.Name, scales='free')+
  scale_color_manual(values=c(6:1))

# minor evidence of run-order effects, which may arise from co-linear batch effects:
## look for batch (and run) effects using the study samples, which will enable us to see batch-effects as well: requires merge, below
mzs <- unique(mz$Molecule.Name)
head(mz)

w <- pivot_wider(dplyr::select(mz, -c(TotalAreaRatioToGlobalStandards, Normalization.Method, Total.Area.Fragment)), names_from = Molecule.Name, values_from = Normalized.Area)

w <- w[w$sampleType=='studySample', ] # only retain study samples
w$sampleType <- droplevels(w$sampleType)

mzdat <- w[ ,mzs]
# normalization and scaling plots
par(mfrow=c(2,2))
boxplot(mzdat, ylab='raw m/z', cex=0.4)
boxplot(log(mzdat), ylab='log m/z', cex=0.4)
boxplot(t(log(mzdat)), xlab='samples', ylab='log m/z', cex=0.4)
boxplot(scale(t(log(mzdat)), scale=F), xlab='samples', ylab='sample-mean-centered log m/z', cex=0.4)

pairs(data.frame('run.order'=w$Run.order, 'batch'=w$Batch, prcomp(mzdat, scale=T)$x[ ,1:6]), pch=20, col=w$Batch, main='PCA of raw peak areas') # note: these are PCs of raw (not logged) data
 
pairs(data.frame('run.order'=w$Run.order, 'batch'=w$Batch, prcomp(log(mzdat), scale=T)$x[ ,1:6]), pch=20, col=w$Batch, main='PCA of log-peak areas') # logged, but not centered

pairs(data.frame('run.order'=w$Run.order, 'batch'=w$Batch, prcomp(t(scale(t(log(mzdat)), scale=F)), scale=T)$x[ ,1:6]), pch=20, col=w$Batch, main='PCA of mean-centered-log-peak areas')
  
### the mean-centered logged metabolite data is an improvement. But, while it appears to reduce batch, and maybe some run effects, it does not eliminate them, which is typical.


# accept the mean-centered loged mz data and evaluate any relationship between batch and the study design:
scat <- as.data.frame(t(scale(t(log(mzdat)), scale=F)))
w[ ,mzs] <- scat # accept log-scaled data

## remember about samples 130 and 103 (at the begining of the file):
colSums(is.na(w))
rowSums(is.na(w))
w[rowSums(is.na(w))>1, ]
w <- w[w$Sample!='130', ] # remove sample 130, this is also the sample with many (24) missing values


mzs
mzNameFix <- tolower(mzs)
mzNameFix <- gsub('n-', 'N-', mzNameFix)
mzNameFix <- gsub('-l-', '-L-', mzNameFix)
mzNameFix <- gsub('l-', 'L-', mzNameFix)
mzNameFix <- gsub('methyL', 'methyl', mzNameFix)
mzNameFix <- gsub('s-', 'S-', mzNameFix)
mzNameFix <- gsub('o-', 'O-', mzNameFix)
mzNameFix 
mzNameFix [mzNameFix=="N-dimethyl-arginine"] <- "N-N-dimethyl-arginine"

w <- w %>% 
  rename_at(vars(mzs), ~ mzNameFix)
colnames(w)

mzs <- mzNameFix

save(w, h, mzs, file='data/mergedData1')




rm(list=ls())
######################################################################
load('data/mergedData1')

# cosider confounding of batch and run order with study design, we may not have batch/run effects if instead the study covariates could better explain mz variation:

## review design:
colnames(w)

ggpairs(w, columns=c('JulianDate', 'YearOfSampling', 'Run.order', 'Batch', 'AM_PM', 'Age', 'Body_Mass', 'stage', 'Sex'), ggplot2::aes(colour = Batch, alpha=0.8)) + theme_bw(base_size = 10) + ggtitle('roe deer study') # there it is.  ...the order of sample collection is highly (completely?) colinear with batch and with run order.

# batch is confounded with:
# collection order
# run order
# Body_Mass*
# Age*
# stage*
# AM_PM*

# (*) these are mostly the last two batches, the last batch is small:
# test sensitivity of the analysis to the presence of the last batch

hist(w$Age, 30)
hist(w$Body_Mass, 20)
plot(Body_Mass ~ Age, w, pch=20)
w$YearOfSampling
w$AM_PM


# adjust data to reduce batch effects

# pre-transpose:
mzDat <- t(w[ ,mzs])

### try segmented model to catch non-linearity in mz ~ age
adj <- w[ ,mzs] 

seg = function(x) segmented(lm(Y ~ Age + Body_Mass + Sex + AM_PM + stage + CohortQuality + JulianDate + YearOfSampling, x), seg.Z = ~ Age) # function to fit segemnted model

for(i in 1:length(mzs)) {
  tmp <- bind_cols(w[ ,mzs[i]], dplyr::select(w, c(Age, Body_Mass, Sex, Batch, AM_PM, stage, CohortQuality, JulianDate, YearOfSampling)))
  colnames(tmp)[1] <- 'Y'
  
  if("Breakpoint estimate(s) outdistanced to allow finite estimates and st.errs" != tryCatch( seg(tmp), warning = function(w) {return(w$message)})[1]) { 
    segM <- seg(tmp)
    breakPoint <- segM$psi[ ,'Est.']
    M <- model.matrix( ~ Age + pmax(Age - breakPoint, 0) + Body_Mass + Sex + AM_PM + stage + CohortQuality + JulianDate + YearOfSampling, data=tmp)
    atmp <- sva::ComBat(mzDat, batch=tmp$Batch, mod = M, par.prior=T, prior.plots=F)
    adj[ ,mzs[i]] <- atmp[mzs[i], ] 
    
  } else {
    M <- model.matrix( ~ Age + Body_Mass + Sex + AM_PM + stage + CohortQuality + JulianDate + YearOfSampling, data=tmp)
    atmp <- sva::ComBat(mzDat, batch=tmp$Batch, mod = M, par.prior=T, prior.plots=F)
    adj[ ,mzs[i]] <- atmp[mzs[i], ] 
  }
}

par(mfrow=c(1,2))
plot(prcomp(w[ ,mzs], scale=T)) 
plot(prcomp(adj[ ,mzs], scale=T)) 


w[ ,mzs] <- adj # accept batch adjustment with segmented age model
d <- w

save(d, mzs, file="data/postBatch_segmentalAge_Data")




## did, the handling of batch effects sharpen any effects of run order?
###############################################################################
rm(list=ls())

load('data/mergedData1') # loads 'w' pre-batch data
load('data/postBatch_segmentalAge_Data') # loads 'd' post-batch data

pairs(data.frame('run.order'=w$Run.order, 'batch'=w$Batch, prcomp(w[ ,mzs], scale=T)$x[ ,1:6]), pch=20, col=w$Batch, main='pre-batch corrected')

pairs(data.frame('run.order'=d$Run.order, 'batch'=d$Batch, prcomp(d[ ,mzs], scale=T)$x[ ,1:6]), pch=20, col=d$Batch, main='post-batch corrected')
## clearly removing batch effects made clear run-order effects a prominent source of mz variation

## do run order effects obviously manfest in individual mzs?
l <- pivot_longer(d, mzs)

ggplot(l, aes(y=value, x=Run.order, color=Batch)) +
  geom_point(size=0.5)+
  theme_classic(base_size = 7)+
  facet_wrap(~name, scales='free')
## clearly they do

# adjusting for run order with an mz ~ run order model (perhaps a spline to account for non-linearity) at this point (post-batch) is ideal as its effect are more clear (fitable)


d <- as.data.frame(d) # convert to df from tibble

residList <- list()

for(i in 1:length(mzs)) {
  residList[[i]] <- residuals(lm(d[ ,mzs[i]] ~ Run.order, d)) }


names(residList) <- mzs
residMat <- do.call(bind_cols, residList)

r <- d
r[ ,mzs] <- residMat

## look at the fix:
l <- pivot_longer(r, mzs)

ggplot(l, aes(y=value, x=Run.order, color=Batch)) +
  geom_point(size=0.5)+
  theme_classic(base_size = 7)+
  facet_wrap(~name, scales='free')

## is some of the non-linearity of mz ~ run order due to field year or date?
ggplot(l, aes(y=value, x=Run.order, color=YearOfSampling)) +
  geom_point(size=0.5)+
  theme_classic(base_size = 7)+
  facet_wrap(~name, scales='free') # looks like there is variation by field year (e.g., valine, MTA, N-Ac-phe)


## accept the run order adjustment and proceed:
d <- r

save(d, mzs, file="data/postBatch_postRunOrder_Data")

#####################################################################





rm(list=ls())
######################################################################
# pca to evaluate sources of variation
######################################################################
load('data/postBatch_postRunOrder_Data')

pca <- prcomp(d[ ,mzs], scale=T)
propvar <- pca$sdev/sum(pca$sdev)
sigEigs <- 6 # will do real TW test in the main analysis
PCs <- paste0('PC', 1:sigEigs)
cumulativeVar <- print(sum(propvar[1:sigEigs])) # the top _ PCs explain _% of the metabolome variance

pca <- data.frame(d[ ,!colnames(d) %in% c(mzs)], pca$x[ ,PCs])
head(pca)

l <- as.data.frame(pivot_longer(pca, PCs, names_to = 'PC'))

# consider technical/environmental/biological effects:
# batch (had been removed via ComBat)
# run order
# year of sampling (? biological/environmental or technical?)
# biological confounds (AM/PM, Body_Mass, 'Cohort quality', stage) # CBC will be considered later

## test segmented age term:
i=1
car::Anova(segmented(lm(value ~ Age + AM_PM + Body_Mass + Sex + stage + CohortQuality + YearOfSampling + JulianDate + Batch + Run.order, l[l$PC==PCs[i], ]), seg.Z = ~ Age), type="III") # NOTE: U1.Age is the difference in slopes btw treshold, and psi.Age is any threshold that is identified

aList <- list()

for(i in 1:length(PCs)) {
  aList[[i]] <- car::Anova(segmented(lm(value ~ Age + AM_PM + Body_Mass + Sex + stage + CohortQuality + YearOfSampling + JulianDate + Batch + Run.order, l[l$PC==PCs[i], ]), seg.Z = ~ Age), type="III") }
names(aList) <- PCs


aList[[1]]

# make a Sum of Sq barplot to summarize the sources of variation.
SS <- as.data.frame(sapply(aList, function(x) x$`Sum Sq`))
SS$term <- rownames(aList[[1]])
SS <- SS[!SS$term %in% c('Residuals', '(Intercept)'), ]
SS <- pivot_longer(SS, cols=PCs, names_to = 'PC', values_to = 'SumOfSquares')
head(SS)
SS <- SS[SS$term!= 'psi1.Age', ] # remove the summary of breakpoint
SS$term <- as.factor(SS$term)
levels(SS$term)

SS$source <- ifelse(SS$term %in% c('Batch', 'Run.order'), 'technical', 'biological') # distinguish technical from biological variation
SS$source <- as.factor(SS$source)

head(SS)
largestEffect <- SS %>% group_by(term) %>% summarise(maxSS = max(SumOfSquares))

a1<- ggplot(SS, aes(y=SumOfSquares, x=PC, fill=term))+
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 12)+
  facet_wrap(~source)+
  xlab('')
a1 # batch effects at this point are potentially due to confounding with other effects

byterm <- ggplot(SS, aes(y=SumOfSquares, x=PC, fill=source))+
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 12)+
  facet_wrap(~term, nrow=2)+
  xlab('')

byterm

largestEffect[rev(order(largestEffect$maxSS)), ]

# example plots:
p1 <- ggplot(pca, aes(y=PC1, x=AM_PM, fill=AM_PM))+
  geom_violin(alpha=0.4)+
  geom_jitter(width=0.1, size=1)+
  theme_bw(base_size = 20)+
  theme(legend.position = 'none')

p2 <- ggplot(pca, aes(y=PC6, x=Age))+
  geom_point(size=1)+
  geom_smooth(method='lm', linetype=0)+
  theme_bw(base_size = 20)

p3 <- ggplot(pca, aes(y=PC1, x=Sex, fill=Sex))+
  geom_violin(alpha=0.4)+
  geom_jitter(width=0.1, size=1)+
  theme_bw(base_size = 20)+
  theme(legend.position = 'none')


ggarrange(p1, p2, p3, nrow=1)



################################################################
# adjust for AM_PM 
################################################################
# pre-transpose:
mzDat <- t(d[ ,mzs])

### try segmented model to catch non-linearity in mz ~ age
adj <- d[ ,mzs] 

seg = function(x) segmented(lm(Y ~ Age + Body_Mass + Sex + AM_PM + stage + CohortQuality + JulianDate + YearOfSampling, x), seg.Z = ~ Age) # function to fit segemnted model

for(i in 1:length(mzs)) {
  tmp <- bind_cols(d[ ,mzs[i]], dplyr::select(d, c(Age, Body_Mass, Sex, CohortQuality, JulianDate, YearOfSampling, AM_PM, stage)))
  colnames(tmp)[1] <- 'Y'
  
  if("Breakpoint estimate(s) outdistanced to allow finite estimates and st.errs" != tryCatch( seg(tmp), warning = function(z) {return(z$message)})[1]) { 
    segM <- seg(tmp)
    breakPoint <- segM$psi[ ,'Est.']
    M <- model.matrix( ~ Age + pmax(Age - breakPoint, 0) + Body_Mass + Sex + stage + CohortQuality + JulianDate + YearOfSampling, data=tmp)
    atmp <- sva::ComBat(mzDat, batch=tmp$AM_PM, mod = M, par.prior=T, prior.plots=F)
    adj[ ,mzs[i]] <- atmp[mzs[i], ] 
    
  } else {
    M <- model.matrix( ~ Age + Body_Mass + Sex + stage + CohortQuality + JulianDate + YearOfSampling, data=tmp)
    atmp <- sva::ComBat(mzDat, batch=tmp$AM_PM, mod = M, par.prior=T, prior.plots=F)
    adj[ ,mzs[i]] <- atmp[mzs[i], ] 
  }
}

d[ ,mzs] <- adj # accept adjustment


### rename terms to make more readily understood analysis/figures/tables:
################################################################
fixUs <- colnames(d)[!colnames(d) %in% mzs]
fixUs

fixed <- fixUs
fixed[fixed=="JulianDate"] <- 'date'
fixed[fixed=="CohortQuality"] <- 'cohort_quality'
fixed[fixed=="YearOfSampling"] <- 'field_year'
fixed[fixed=="Body_Mass"] <- 'body_mass'
fixed[fixed=="Sex"] <- 'sex'
fixed[fixed=="Age"] <- 'age'
d <- d %>% 
  rename_at(vars(fixUs), ~ fixed)
colnames(d)


save(d, mzs, file='data/postAMPMdata')
################################################################





rm(list=ls())
######################################################################
load('data/postAMPMdata')

colnames(d)

table(d$stage)
hist(d$age, 20)

demoVars <- c('age', 'body_mass', 'sex', 'stage', 'cohort_quality', 'date', 'field_year') # demographic variables

save(d, mzs, demoVars, file='data/dataforUnivariateModels')
