# roe deer clock:

library(data.table)
library(tidyverse)
library(caret)
library(gridExtra)
library(ggforce)
library(magrittr)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(smatr)
library(ggpubr)

rm(list=ls())

getwd()

setwd('~/Google Drive/My Drive/Documents/RoeDeer')
 # setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/RoeDeer')



# make sex-specific deer clocks, us them to depict the relative age/senescence of the opposite sex:
##########################################################################################################

# load data:
load('data/dataforUnivariateModels')



# load elastic net functions based on those written by Xiaqing Zhao
load('clockBuildingFunction') 

EN_LeaveOneDeerOut # set seed chuck (kinda complicated) to make the exact result the same each time

table(d$stage) # remove juveniles and compare clocks built with and without juvenile data
juvDat<- d
length(unique(juvDat$deer)) # 147 deer, 209 samples

d <- d[d$age>1, ] # omit juveniles 
length(unique(d$deer)) # 84 deer, 134 samples remain

ggplot(d, aes(x=age, fill=sex))+
  geom_histogram(bins=12)+
  theme_classic(base_size = 16)+
  facet_wrap(~sex, ncol=1)


# sex-specific clocks:

EN_LeaveOneDeerOut

mdat <- d[d$sex=='M', ]
fdat <- d[d$sex=='F', ]

## set seed
folds <- 5 # we will use 5-Fold CV 
seedsPerFold <- 3 # I got this from the error thrown when I did it wrong - I hope that halps! 
set.seed(77) # this is required for the call to 'sample()' below
seeds <- vector(mode = "list", length = folds + 1) #length is = folds + 1
for(k in 1:folds) seeds[[k]] <- sample.int(1000,  seedsPerFold) 
seeds[[folds + 1]] <- sample.int(1000, 1)

mClock <-train(
  as.matrix(mdat[ ,mzs]), 
  mdat$age,
  method = "glmnet", 
  trControl = trainControl(
    seed=seeds,
    method = "cv", 
    number = 5,
    verboseIter = F)
)

fClock <-train(
  as.matrix(fdat[ ,mzs]), 
  fdat$age,
  method = "glmnet", 
  trControl = trainControl(
    seed=seeds,
    method = "cv", 
    number = 5,
    verboseIter = F)
)

plot(predict(mClock, d[ ,mzs]) ~ d$age, col=d$sex)
plot(predict(fClock, d[ ,mzs]) ~ d$age, col=d$sex)


d$fClockPred <- predict(fClock, d[ ,mzs])
d$mClockPred <- predict(mClock, d[ ,mzs])


fLodo <- EN_LeaveOneDeerOut(as.matrix(fdat[ ,mzs]), fdat$age, fdat$deer) 
mLodo <- EN_LeaveOneDeerOut(as.matrix(mdat[ ,mzs]), mdat$age, mdat$deer) 

mLodo[[1]]
plot(mLodo[[1]] ~ mdat$age)
abline(0,1, col=2)

plot(fLodo[[1]] ~ fdat$age)
abline(0,1, col=2)

d$fClockPred[d$sex=='F'] <- fLodo[[1]] # the within-sex predictions now come from LODO 
d$mClockPred[d$sex=='M'] <- mLodo[[1]]

l <- pivot_longer(d, c(fClockPred, mClockPred), values_to = 'predAge', names_to = 'clock')

l$clock <- as.factor(l$clock)
levels(l$clock) <- c('female clock', 'male clock') 
levels(l$sex) <- c('female samples', 'male samples')

sexClockPlots <- ggplot(l, aes(predAge, x=age))+
  geom_line(aes(group=deer), alpha=0.5)+
  geom_smooth(method='lm', linetype = 0, alpha=0.3)+
  geom_point(size=1)+
  theme_bw(base_size = 14)+
  labs(y='predicted age (years)', x='age (years)')+
  theme(strip.text = element_text(size = 16),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_line(color = "grey96"))+
  facet_grid(~clock ~ sex, scales='free_x')+
  scale_x_continuous(breaks=c(2,4,6,8,10,12))


d$sexlabs <- d$sex
levels(d$sexlabs) <- c('female samples', 'male samples')

predictionCorPlot <- ggplot(d, aes(mClockPred, fClockPred, group=deer))+
  geom_line(alpha=0.3)+
  geom_abline(slope=1, linetype='dashed', alpha=0.5)+
  geom_point(size=1)+
  theme_bw(base_size = 14)+
  labs(y='female clock prediction (years)', x='male clock prediction (years)' )+
  theme(strip.text = element_text(size = 16),
        strip.background = element_blank(),
        legend.position = 'none')+
  facet_wrap(~ sexlabs, scales='free')


l$clockError <- l$predAge-l$age

errorPlot <- ggplot(l, aes(clockError, x=clock, fill=clock))+
  geom_abline(slope=0, linetype='dashed', alpha=0.5)+
  geom_violin(trim=F, alpha=0.3, color=0)+
  geom_jitter(size=1, width=0.05)+
  theme_bw(base_size = 14)+
  labs(x='', y='clock error (years)')+
  theme(strip.text = element_text(size = 16),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.grid = element_line(color = "grey90"))+
  facet_grid(~ sex, scales='free_x')


par(mfrow=c(2,2))
plot(d$fClockPred[d$sex=='F'] ~ d$age[d$sex=='F'])
plot(d$fClockPred[d$sex=='M'] ~ d$age[d$sex=='M'])
plot(d$mClockPred[d$sex=='F'] ~ d$age[d$sex=='F'])
plot(d$mClockPred[d$sex=='M'] ~ d$age[d$sex=='M'])

summary(lm(d$fClockPred[d$sex=='F'] ~ d$age[d$sex=='F'])) # female pred female
sqrt(mean((d$fClockPred[d$sex=='F'] - d$age[d$sex=='F'])^2)) # RMSE

summary(lm(d$fClockPred[d$sex=='M'] ~ d$age[d$sex=='M'])) # female pred male
sqrt(mean((d$fClockPred[d$sex=='M'] - d$age[d$sex=='M'])^2)) # RMSE

summary(lm(d$mClockPred[d$sex=='F'] ~ d$age[d$sex=='F'])) # male pred female
sqrt(mean((d$mClockPred[d$sex=='F'] - d$age[d$sex=='F'])^2)) # RMSE

summary(lm(d$mClockPred[d$sex=='M'] ~ d$age[d$sex=='M'])) # male pred male
sqrt(mean((d$mClockPred[d$sex=='M'] - d$age[d$sex=='M'])^2)) # RMSE

modList <- list(lm(d$fClockPred[d$sex=='F'] ~ d$age[d$sex=='F']), 
                lm(d$fClockPred[d$sex=='M'] ~ d$age[d$sex=='M']),
                lm(d$mClockPred[d$sex=='F'] ~ d$age[d$sex=='F']),
                lm(d$mClockPred[d$sex=='M'] ~ d$age[d$sex=='M'])
                )
names(modList) <- c('fclock pred female', 'fclock pred male', 'mclock pred female', 'mclock pred male')

sList <- lapply(modList, summary)
sList[['fclock pred female']]$coefficients
intercepts <- t(sapply(sList, function(x) return(x$coefficients[1, ])))
slopes <- t(sapply(sList, function(x) return(x$coefficients[2, ])))

intercepts
slopes

summary(lm(fClockPred ~ age * sex, d))
summary(lm(mClockPred ~ age * sex, d))



var.test(l$clockError[l$sex=='male samples']  ~ l$clock[l$sex=='male samples']) # ratio of variances = 1.06, p=0.821

var.test(l$clockError[l$sex=='female samples']  ~ l$clock[l$sex=='female samples']) # ratio of variances = 1.06, p=0.810

ggarrange(sexClockPlots,
ggarrange(predictionCorPlot, errorPlot, ncol=1, heights=c(1, 1.2)), nrow=1, widths=c(1.3,1))

