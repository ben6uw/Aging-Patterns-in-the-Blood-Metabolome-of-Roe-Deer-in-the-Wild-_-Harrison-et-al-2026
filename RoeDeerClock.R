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
dir("~")
 setwd('~/Google Drive/My Drive/Documents/RoeDeer')
 # setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/RoeDeer')

 dir()


rm(list=ls())
# deer clock:
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

mean(d$age) 

################################################################################################################################################
## Custom function to get predictions from leave one subject (deer) out, train by 5-fold CV on remaining data, then predicting on the held-out deer.  Repeated until all deer had been held out and their age(s) predicted

# [Note]: these elastic net functions are fully-reproducible ONLY under the CV regime and the tuning grid used in its train and trainControl steps.  The user will need to modify the seed-setting chunk (annotated in the code) to make their ML reproducible


lodo <- EN_LeaveOneDeerOut(as.matrix(d[ ,mzs]), d$age, d$deer) 

lodo.wjuveniles <- EN_LeaveOneDeerOut(as.matrix(juvDat[ ,mzs]), juvDat$age, juvDat$deer) 


#############################################################################
# output is a 2-element list:
# the first element is the age predictions for the deer, the second is the variable importance

par(mfrow=c(1,2))
plot(lodo[[1]] ~ d$age)
plot(lodo.wjuveniles[[1]] ~ juvDat$age)

d$lodo <- lodo[[1]]
juvDat$lodo <- lodo.wjuveniles[[1]]

save(juvDat, lodo.wjuveniles, d, lodo, file='clockTraining')
 



rm(list=ls())
#################################################
load('data/dataforUnivariateModels')
load('clockTraining')

ggplot(d, aes(lodo, x=age))+
  geom_point(size=1)+
  geom_line(aes(group=deer, alpha=0.5))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  geom_smooth(method='lm', linetype=0)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14))+
  ylab('prediction of left-out-deer (years)')+
  xlab('age (years)')

ggplot(d, aes(lodo, x=age))+
  geom_abline(linetype='dashed', color='grey40')+
  geom_smooth(method='lm', linetype=0)+
  geom_line(aes(group=deer, alpha=0.5))+
  geom_point(size=1)+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14))+
  ylab('prediction of left-out-deer (years)')+
  xlab('age (years)')

ggplot(juvDat, aes(lodo, x=age))+
  geom_point(size=1)+
  geom_line(aes(group=deer, alpha=0.5))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  geom_smooth(method='lm', linetype=0)+
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14))+
  ylab('prediction of left-out-deer (years)')+
  xlab('age (years)')

## calculate root mean squared error (RMSE):
sqrt(mean((d$age - d$lodo)^2)) 
sqrt(mean((juvDat$age - juvDat$lodo)^2))


summary(m <- lmerTest::lmer(lodo ~ age + (1|deer), d)) # P= 6.9x10^-10, ageB= 0.38
MuMIn::r.squaredGLMM(m) 
# R2m = the variance explained by the fixed effects = 0.30       
# R2c = the var exp. by the entire model = 0.68

summary(m <- lmerTest::lmer(lodo ~ age + (1|deer), juvDat)) # p <2x10^-16, slope=0.50 
MuMIn::r.squaredGLMM(m) # # R2m = the variance explained by the fixed effects = 0.50
# (R2c = the var exp. by the entire model) 0.68

## is the apparent accuracy of the all-deer clock sensitive to the sample size increase from the juvenile data?
## also, would the adult clock appear more accurate if it predcted ages of the juveniles and included those in the regression?

summary(m <- lmerTest::lmer(lodo ~ age + (1|deer), juvDat[juvDat$age>1, ])) # p 1.2e-08, slope=0.29
MuMIn::r.squaredGLMM(m) # # R2m = the variance explained by the fixed effects = 0.27
# (R2c = the var exp. by the entire model) 0.51

sqrt(mean((juvDat$age[juvDat$age>1] - juvDat$lodo[juvDat$age>1])^2))


## compare residuals over age from the clocks with and without juvenile data:

d$AA <- residuals(lm(lodo ~ age, d)) # age acceleration
d$dataSource <- 'adult'

juvDat$AA <- residuals(lm(lodo ~ age, juvDat)) # age acceleration
juvDat$dataSource <- 'all_deer'

comp <- bind_rows(dplyr::select(juvDat, deer, dataSource, age, lodo, AA), dplyr::select(d, deer, dataSource, age, lodo, AA))

head(comp)
l <- pivot_longer(comp, c(lodo, AA))
head(l)

l$dataSource <- factor(l$dataSource, levels=c('all_deer', 'adult'))
l$name <- as.factor(l$name)
levels(l$name) <- c('age accelteraion', 'predicted age')
l$name <- factor(l$name, levels=c('predicted age', 'age accelteraion'))

agePlot <- print(ggplot(subset(l, name=='predicted age'), aes(y=value, x=age, color=dataSource))+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text.x = element_text(size = 14))+
  geom_smooth(method='lm', linetype=0)+
  geom_point(size=1)+
  scale_x_continuous(breaks=c(1, 2,4,6,8,10,12,14))+
  ylab('prediction (years)')+
  xlab('age (years)')+
  facet_wrap(~ dataSource))
  
AAplot <- print(ggplot(subset(l, name=='age accelteraion'), aes(y=value, x=age, color=dataSource))+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text.x = element_text(size = 14))+
  geom_smooth(method='lm', se=F, linewidth = 0.8)+
  geom_smooth(linetype=0)+
  geom_point(size=1)+
  scale_x_continuous(breaks=c(1, 2,4,6,8,10,12,14))+
  ylab('age accelteraion (years)')+
  xlab('age (years)')+
  facet_grid(~dataSource))

ggarrange(agePlot, AAplot, ncol=1)




## does anything associate with age acceleration?

rm(list=ls())
######################################################################################
load('data/dataforUnivariateModels')
load('clockTraining')

library(lmerTest)

## look among demographic variables for assoc with AA:

d$AA <- residuals(lm(lodo ~ age, d)) # get age acccel.

colnames(d)
### include age to see if there may be an effect on AA ~ age x X, where X could be mass or sex.  Such an effect may occur if the rate of aging (slope) differed by weight or sex.
m <- lmer(AA ~ (age * body_mass * sex)^2 + cohort_quality + date + field_year + (1|deer), d) # I leave age in here to see if there are age-dependent interactions between AA and demographic vars.  
m <- print(lmerTest::step(m, direction = "both"))
m

ggplot(d, aes(AA, x=body_mass))+
  geom_hline(yintercept=0, linetype='dashed')+
  geom_point(size=1)+
  geom_line(aes(group=deer, alpha=0.5))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  geom_smooth(method='lm', linetype=0)+
  ylab('age accelteration (yrs)')

### I added lines to connect data from each deer.  If these associations are consistent within deer, then the global slope and the within deer slopes should agree to some extent.
## for AA ~ weight, there is not a convincing global effect, nor much agreement in within deer.

### I do not think we have a good case for any of the demographic variables (weight, sex, Cohort Quality and their interactions with age) assoc. with AA



rm(list=ls())

##################################################################################
# variable importance:
load('clockTraining')
d 

# the second element is the variable-importance, from varImp()), from each model
v <- lodo[[2]]
class(v[[1]])

v2 <- lapply(v, function(x) data.frame(importance = x$importance))
imps <- as.matrix(do.call(cbind, v2))
imps[1:4,1:4]
hist(colSums(imps!=0), 20) # these are the total number of features selected by the elastic nets
mean(colSums(imps!=0))
range(colSums(imps!=0))

imps <- imps[rev(order(rowMeans(imps))), ]
tmp <- rownames(imps)

colnames(imps) <- paste0('model_', 1:ncol(imps))
imps <- data.frame(mz=rownames(imps), imps)
imps$mz <- factor(imps$mz, levels=rev(tmp))

l <- pivot_longer(imps, -mz, values_to = 'importance', names_to='interation')

head(l)

ggplot(l, aes(y=importance, x=mz))+
  geom_jitter(width=0.1, size=0.2, alpha=0.5)+
  coord_flip()+
  theme_classic(base_size = 10)+
  xlab('')+
  ggtitle('variable importance among leave-deer-out builds')

ggplot(subset(l, mz %in% tmp[1:20]), aes(y=importance, x=mz))+
  geom_jitter(width=0.1, size=0.2, alpha=1)+
  coord_flip()+
  theme_classic(base_size = 16)+
  xlab('')


## where are the age-associated metabolites among these?
top20 <- tmp[1:20]

l <- pivot_longer(d, top20)
l$name <- factor(l$name, levels=top20)

p1 <- ggplot(l, aes(y=value, x=age))+
  geom_point(size=0.5)+
  scale_x_sqrt(breaks=c(2,4,6,8,10,12,14))+
  geom_smooth(alpha=0.2, lwd=0.5)+
  theme_classic()+
  facet_wrap(~name, scales='free')+
  ylab('normalized abundance')+
  xlab('age (years)')+
  theme(strip.text.x = element_text(size = 7), axis.text.x = element_text(size = 7))

p2 <- ggplot(d, aes(lodo, x=age))+
  geom_smooth(alpha=0.2, lwd=0.5)+
  geom_point(size=1)+
  scale_x_sqrt(breaks=c(2,4,6,8,10,12,14))+
  theme_classic(base_size = 16) +
  ylab('predicted age (years)')+
  xlab('age (years)')


ggarrange(p1, p2, widths=c(1, 0.7), common.legend = T, legend='right')


#### are clock-important metabolites enriched for age-assoc.?
load('uniPs') # these are the p-values from univariate mixed models of metabolite abundance 

meanImp <- data.frame(mz=imps$mz, meanImp = rowMeans(imps[ ,-c(1)]))
head(meanImp)
x <- left_join(meanImp, uniPs)

head(x)

dev.off()
plot(x$meanImp, -log(x$age), pch=20)
cor.test(x$meanImp, -log(x$age), method='kendall', alternative='greater')

top20

plot(x$meanImp, -log(x$age), pch=20)
cor.test(x$meanImp, -log(x$age2), method='kendall', alternative='greater')



