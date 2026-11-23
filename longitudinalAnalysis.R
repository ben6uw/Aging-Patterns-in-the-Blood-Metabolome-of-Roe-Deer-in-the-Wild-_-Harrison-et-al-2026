library(data.table)
library(tidyverse)
library(survival)
library(ggpubr)
library(survminer)
library(lmerTest)


# analysis of survival data from roe deer:

# email from Jeff, Jan 2025:

# With my colleague Jean-Michel (cc to this email) we have now finished to prepare the longevity data for the individuals for which with have information on their metabolome. These data are summarize in the columns R, S, T and U of the attached file. But let me give you a few explanations regarding these data.

# If the individuals are considered dead (i.e. ‘yes’ in the column R), there are three possible options depending on their causes of death: 

# i) 'Natural death' when the dead animal is retrieved in the reserve without any sign of accident. For these individuals, we know their exact lifespan (provided in the lifespan column, T).
# ii) 'Human caused' when the cause of death is hunting (in most cases). For these animals, we obviously have no information on their ’natural lifespan’ and the lifespan data can only be interpreted as the minimum age reached by the individuals. For instance, if we take the example of the individual F338, we entered ‘>7’in the lifespan data as the only biologically relevant information we have is that this females reached 7 years of age.
# iii) ‘Unknown’ when the individual was not observed within the last 3 years. For these individuals, we entered the age at last contact as the minimum age reached by these individuals.

# If the individuals are considered alive (i.e. ’NA’ in the column R), this is obviously much more straightforward. This means that we have evidence that these animals were alive within the last 3 years and we simply entered the age at the last observation (minimum lifespan) in the lifespan column. The date at last contact is also provided. As you will see, only 15 individuals from this dataset are still potentially alive.


## UPDATE April 8th 2025:
# Dear Ben,
# Please find attached the updated dataset. [RoeDeer-LIFE HISTORY DATA-V6.xlsx] As you will see, we have increased the resolution of the age and longevity data (for individuals who died of natural causes). We came back to the original files and refined the period of death by thirds of a year (for example, 8 years, 8.37 years or 8.67 years old). This works very well, except for one individual for whom capture and death are in the same period, but one month apart. For this individual (F351) we have therefore shifted the longevity by about a month. All the individuals now have a longevity later than the last capture. Do not hesitate if you have further questions,

# All best wishes, Jeff



rm(list=ls())

dir('~')
# setwd('~/brharrison2000@gmail.com - Google Drive/My Drive/Documents/RoeDeer')
setwd('~/Google Drive/My Drive/Documents/RoeDeer')
dir('data') # 

####################################################################################
s <- as.data.frame(readxl::read_xlsx('data/RoeDeer-LIFE HISTORY DATA-V6.xlsx'))
head(s)

# one fix to the data, re Jeff's email:
# sent email fo Jeff, about deer F387 "She was killed (hunted) at 13.5 years old (November 2018) so the correct value is > 13.5 (and not > 11.5)."
s[s$`Id_JMG-UniqueIndivual-ID`=='F387', ] # these data are corrected
###

table(table(s$idkit)) # 'idkit is unique to each sample, can use as unique identifier later.
table(table(s$`Id_JMG-UniqueIndivual-ID`)) # observations per deer

s$deer <- s$`Id_JMG-UniqueIndivual-ID`
s <- s %>% dplyr::rename('age'='Age', 'sex' = 'Sex', 'date' = 'JulianDate', 'body_mass' = 'Body_Mass', 'field_year'= 'YearOfSampling', 'cohort_quality' = "CohortQualityIndex-StandardiedJuvenieleMass-Kg")

table(s$sex, s$deer)
length(unique(s$deer)) # survival data for 147 deer
s$sex <- as.factor(s$sex)
table(s$pop) # pop is invariant
s$field_year <- as.factor(s$field_year)
levels(s$field_year)

colnames(s)
s$Cohort # year of birth
s$Cohort <- as.factor(s$Cohort)
levels(s$Cohort)

head(s)

table(s$CauseOfDeath)
# note: need 0/1 status for most survival models, so convert the cause of death 0-4 status to 0/1 for each cause:

## recruitment age (age at first blood sample) 
head(s)
youngestSample <- as.data.table(s)[ , .SD[which.min(age)], by = deer] # not necessarily P1
head(youngestSample)
recruitmentAge <- youngestSample$age
names(recruitmentAge) <- youngestSample$deer
s$recruitmentAge <- recruitmentAge[s$deer]

head(s)

LS <- as.numeric(gsub('>', '', s$Lifespan))
LS >=  s$age # all deer census ages are older than the sample data = good!
s$censusAge <-  as.numeric(gsub('>', '', s$Lifespan))

save(s, file='data/survivalDat1')



rm(list=ls())
##############################################################################
load('data/survivalDat1')

##############################################################################
### Temporary: for plotting to show the longitudinal blood collection data
tmp <- s

censusAge <- data.frame('deer'=tmp$deer, censusAge=tmp$censusAge) 
censusAge <- censusAge[!duplicated(censusAge), ]
head(censusAge)
orderedDeer <- censusAge$deer[order(censusAge$censusAge)]

tmp$deer <- factor(tmp$deer, levels=orderedDeer)
table(table(tmp$deer)) # events per deer (samples or death/censor)

orderedRecruites <- dplyr::select(tmp, c(deer, recruitmentAge))
orderedRecruites <- orderedRecruites[!duplicated(orderedRecruites), ]
orderedRecruites  <- orderedRecruites$deer[order(orderedRecruites$recruitmentAge)]

head(tmp)

l <- pivot_longer(tmp[ ,c('deer', 'CauseOfDeath', 'age', 'censusAge', 'recruitmentAge')], c('age', 'censusAge'), values_to = 'age', names_to = 'event')
head(l)
l <- l[!duplicated(l), ]

table(l$CauseOfDeath)
l$status <- 0
head(l)
l$status[l$event=='censusAge' & !is.na(l$CauseOfDeath)] <- 1
head(l)
l$event[l$event=='age'] <- 'sample'
l$event[l$event=='censusAge'] <- 'census'

l$fate <- l$CauseOfDeath
head(l)
l$fate[l$event=='sample'] <- 'sample'
table(l$fate)
l$fate[l$fate=='NA'] <- 'censored'
l$fate[l$fate=='Unknown'] <- 'censored'
l$fate[l$fate=='Human caused'] <- 'censored'
l$fate[l$fate=='Natural'] <- 'dead'

ggplot(l, aes(y=deer, x=age, group=deer, color=event))+
  geom_line(color='grey', size=0.5)+
  geom_point(size=1)+
  scale_color_manual(values=c('grey40', 'red'))+
  theme_classic(base_size = 20)+
  theme( axis.text.y=element_blank(), 
         axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks = round(1:20)) +
  ylab('147 deer')+
  xlab('age (years)')

ggplot(l, aes(y=deer, x=age, group=deer, color=fate))+
  geom_line(color='grey', linewidth=0.3)+
  geom_point(size=1)+
  theme_classic(base_size = 18)+
  theme( axis.text.y=element_blank(), 
         axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks = round(1:20)) +
  ylab('147 deer')+
  xlab('age (years)')+
  labs(color='')+
  scale_color_manual(values=c('orange', 'purple','grey'))

table(l$fate)

ggplot(subset(l, fate!='sample'), aes(x=age)) +
  geom_histogram(color="#e9ecef", alpha=0.9, position = 'identity', bins=6) +
  theme_bw(base_size = 16)+
  ylab('deer')+
  xlab('age (years)')+
  facet_wrap(~fate, ncol=1, scales='free')


l$CauseOfDeath[l$CauseOfDeath=='NA'] <- 'Unknown'
l$CauseOfDeath
l$CauseOfDeath <- factor(l$CauseOfDeath, levels=c('Natural', 'Human caused', 'Unknown'))

levels(l$CauseOfDeath) <- tolower(levels(l$CauseOfDeath))

ggplot(subset(l, fate!='sample'), aes(x=age)) +
  geom_histogram(color="#e9ecef", alpha=0.9, position = 'identity', bins=10) +
  theme_bw(base_size = 16)+
  ylab('deer')+
  xlab('age (years)')+
  facet_wrap(~CauseOfDeath, ncol=1, scales='free_y')



rm(list=ls())
#########################
# add metabolome data
load('data/survivalDat1')
load('data/dataforUnivariateModels') 


# match variable types:
inter <- as.character(intersect(colnames(s), colnames(d)))
inter
str(s[ ,inter])
str(d[ ,inter])

s$Id_Metabolome <- as.character(s$Id_Metabolome)
d$Date <- as.Date(d$Date, format="%m/%d/%y")
s$Date <- as.Date(s$Date)
s$Batch <- as.factor(s$Batch)
head(s)

d$Id_Metabolome <- as.character(d$Id_Metabolome)

dat <- full_join(s, d) 
head(dat)

colnames(dat)
dat$etime <- dat$age - dat$recruitmentAge # event time (etime) the time since recruitment (in years)

## get status data:
statusData <- dplyr::select(dat, c(deer, censusAge, CauseOfDeath, recruitmentAge))
head(statusData)

## see email from Jeff at beginning of code about the meaning of the causes or death:
statusData$status <- 0
statusData$status[statusData$CauseOfDeath %in% c('Natural')] <- 1 # only 'Natural' (found dead) will be status=1, the others, Human caused (right-censor), 'unknown' (lost contact for 3+yrs, with date of last contact) will be treated as 'lost of followup' and status age for them will be the last known contact, and NA = they have evidence that they are still alive.
statusData <- statusData[!duplicated(statusData), ] 
table(statusData$status, statusData$CauseOfDeath)
head(statusData)
colSums(is.na(statusData))
statusData <- statusData[rowSums(is.na(statusData))==0, ]
statusData$futime <- statusData$censusAge - statusData$recruitmentAge


# set the initial covariates for deer (these cannot change over time). Vars that change with time will be added below
d0 <- dat[ ,c('deer', 'recruitmentAge', 'sex', 'cohort_quality', 'CauseOfDeath')] 
d0 <- d0[!duplicated(d0), ] # this step has to give unique deer
colSums(is.na(d0)) # look at NA counts
d0 <- d0[rowSums(is.na(d0))==0, ]

d0 <- left_join(d0, dplyr::select(statusData, c(deer, futime)))

## add status:
statusData <- statusData[statusData$deer %in% d0$deer, ]
ldat <- tmerge(d0, statusData, deer, death = event(futime, status)) 

head(ldat)
ldat[1:10, ]
table(ldat$death)

## add alternative causes of death as status:
head(statusData)
statusData$human.caused <- ifelse(statusData$CauseOfDeath=='Human caused', 1, 0)
head(statusData)

ldat <- tmerge(ldat, statusData, deer, human.caused = event(futime, human.caused)) 
head(ldat)

hist(statusData$futime)
mean(statusData$futime)
median(statusData$futime)
range(statusData$futime)
table(statusData$status)
table(statusData$CauseOfDeath)

## add blood, weight and other time-dependent data:
colnames(dat)

ldat <- tmerge(ldat, dat, deer, body_mass = tdc(etime, body_mass)) 
ldat <- tmerge(ldat, dat, deer, AM_PM = tdc(etime, `Morning_Evening(Caputre)`)) 


## add mz data
varString <- character()
for(i in 1:length(mzs)){
  varString[i] <- paste0('`', mzs[i], '`', '=tdc(etime, ', '`',  mzs[i], '`),') }

cat(noquote(varString)) # prints covar input so that you can manually past into line below. Silly

cat(noquote(varString[40:72])) 


ldat <- tmerge(ldat, dat, id=deer, `1-methyladenosine`=tdc(etime, `1-methyladenosine`), `1-methyl-L-histidine`=tdc(etime, `1-methyl-L-histidine`), `1-methylnicotinamide`=tdc(etime, `1-methylnicotinamide`), `2-aminoisobutyrate`=tdc(etime, `2-aminoisobutyrate`), `3-methoxytyrosine`=tdc(etime, `3-methoxytyrosine`), `4-acetamidobutanoate`=tdc(etime, `4-acetamidobutanoate`), `adenosine`=tdc(etime, `adenosine`), `alanine`=tdc(etime, `alanine`), `anserine`=tdc(etime, `anserine`), `arginine`=tdc(etime, `arginine`), `aspartate`=tdc(etime, `aspartate`), `azelate`=tdc(etime, `azelate`), `biotin`=tdc(etime, `biotin`), `carnosine`=tdc(etime, `carnosine`), `choline`=tdc(etime, `choline`), `citrulline`=tdc(etime, `citrulline`), `cortisol`=tdc(etime, `cortisol`), `creatine`=tdc(etime, `creatine`), `creatinine`=tdc(etime, `creatinine`), `cystathionine`=tdc(etime, `cystathionine`), `cytidine`=tdc(etime, `cytidine`), `deoxycarnitine`=tdc(etime, `deoxycarnitine`), `deoxycytidine`=tdc(etime, `deoxycytidine`), `deoxyguanosine`=tdc(etime, `deoxyguanosine`), `glutamate`=tdc(etime, `glutamate`), `glutamine`=tdc(etime, `glutamine`), `glutarylcarnitine`=tdc(etime, `glutarylcarnitine`), `glycocholate`=tdc(etime, `glycocholate`), `guanosine`=tdc(etime, `guanosine`), `hippurate`=tdc(etime, `hippurate`), `histidine`=tdc(etime, `histidine`), `hypoxanthine`=tdc(etime, `hypoxanthine`), `indole-3-acetamide`=tdc(etime, `indole-3-acetamide`), `indole-3-acetate`=tdc(etime, `indole-3-acetate`), `inosine`=tdc(etime, `inosine`), `isoleucine`=tdc(etime, `isoleucine`), `kynurenine`=tdc(etime, `kynurenine`), `L-carnitine`=tdc(etime, `L-carnitine`), `leucine`=tdc(etime, `leucine`), `lumichrome`=tdc(etime, `lumichrome`), `lysine`=tdc(etime, `lysine`), `methionine`=tdc(etime, `methionine`), `methylthioadenosine`=tdc(etime, `methylthioadenosine`), `N-trimethyllysine`=tdc(etime, `N-trimethyllysine`), `N-N-dimethyl-arginine`=tdc(etime, `N-N-dimethyl-arginine`), `N-acetylleucine`=tdc(etime, `N-acetylleucine`), `N-acetylphenylalanine`=tdc(etime, `N-acetylphenylalanine`), `N-acetyltryptophan`=tdc(etime, `N-acetyltryptophan`), `nicotinamide`=tdc(etime, `nicotinamide`), `O-acetylcarnitine`=tdc(etime, `O-acetylcarnitine`), `ornithine`=tdc(etime, `ornithine`), `oxoproline`=tdc(etime, `oxoproline`), `pantothenate`=tdc(etime, `pantothenate`), `paraxanthine`=tdc(etime, `paraxanthine`), `phenylalanine`=tdc(etime, `phenylalanine`), `phenylethanolamine`=tdc(etime, `phenylethanolamine`), `phosphorylcholine`=tdc(etime, `phosphorylcholine`), `pipecolate`=tdc(etime, `pipecolate`), `proline`=tdc(etime, `proline`), `riboflavin`=tdc(etime, `riboflavin`), `S-adenosylhomocysteine`=tdc(etime, `S-adenosylhomocysteine`), `S-adenosylmethionine`=tdc(etime, `S-adenosylmethionine`), `sebacate`=tdc(etime, `sebacate`), `serine`=tdc(etime, `serine`), `spermidine`=tdc(etime, `spermidine`), `threonine`=tdc(etime, `threonine`), `thymine`=tdc(etime, `thymine`), `tryptophan`=tdc(etime, `tryptophan`), `tyrosine`=tdc(etime, `tyrosine`), `urocanate`=tdc(etime, `urocanate`), `valine`=tdc(etime, `valine`), `xanthosine`=tdc(etime, `xanthosine`))


save(ldat, mzs, d, dat, file='longitudinalData')




rm(list=ls())
######################################################################
load('longitudinalData')
load('clockTraining') # add age-acceleration from the clock (loads a new df: 'd'):


d$AA <- residuals(lm(lodo ~ age, d)) # get age acccel.

p1 <- print(ggplot(d, aes(lodo, x=age))+
        stat_smooth(method = "lm", size=0.8)+
        geom_line(aes(group=deer, alpha=0.5))+
        geom_point(size=1)+
        theme_bw(base_size = 12)+
        scale_x_sqrt(breaks=c(2,4,6,8,10,12,14,16)) +
        ylab('predicted age (years)')+
        xlab('age (years)') +
        theme(legend.position = 'none'))


p2 <- print(ggplot(d, aes(AA, x=age))+
              stat_smooth(method = "lm", size=0.8)+
              geom_line(aes(group=deer, alpha=0.5))+
              geom_point(size=1)+
              theme_bw(base_size = 12)+
              scale_x_sqrt(breaks=c(2,4,6,8,10,12,14,16)) +
              ylab('"age acceleration" (years)')+
              theme(legend.position = 'none'))


ggarrange(p1, p2)



#######################################################################################
# age acceleration and deer mortality

d$Id_Metabolome <- as.character(d$Id_Metabolome)
dat <- left_join(dat, dplyr::select(d, c(Id_Metabolome, AA)))

ldat <- tmerge(ldat, dat, deer, AA = tdc(etime, AA)) 

## use stepwise regression to identify the best model
x <- dplyr::select(ldat, c(tstart, tstop, death, AA, recruitmentAge, sex, body_mass))
rowSums(is.na(x)) # this will have 75 NAs for the fawn AA
x <- x[rowSums(is.na(x))==0, ]  # omit these
m <- step(coxph(Surv(tstart, tstop, death==1) ~ AA + recruitmentAge * sex * body_mass, na.action = na.exclude, x))

m


cox.zph(m) # test PH assumption

## scale the variables to get comparable HRs:
scat <- ldat
scat$sex <- scale(as.numeric(scat$sex))
scat$recruitmentAge <- scale(scat$recruitmentAge)
scat$AA <- scale(scat$AA)
scat$body_mass <- scale(scat$body_mass)

m <- coxph(Surv(tstart, tstop, death==1) ~ AA + recruitmentAge * sex , scat) # fit the selected model to the scaled data to get comparable hazards across terms, whos units could vary.  After scaling the units are hazard per SD 
summary(m)

coxM <- summary(m)
coxM <- data.frame(coxM$coefficients)
coxM$term <- rownames(coxM)
coxM

coxM$coef
coxM$LCI <- coxM$coef-1.96*coxM$se.coef.
coxM$UCI <- coxM$coef+1.96*coxM$se.coef.

# manual renaming:
termNames <- coxM$term 
termNames
termLabels <- termNames

termLabels <- gsub('AA', 'age accel.', termLabels)
termLabels <- gsub(':', ' x ', termLabels)
termLabels <- gsub('sex', 'sex (male)', termLabels)
termLabels <- gsub('recruitmentAge', 'age', termLabels)
termLabels <- tolower(termLabels)
termLabels

names(termLabels) <- termNames

termLabels[coxM$term]
coxM$termLabels <- termLabels[coxM$term]
coxM
as.character(coxM$termLabels[order(coxM$coef)])
coxM$termLabels <- factor(coxM$termLabels, levels=as.character(coxM$termLabels[order(coxM$coef)]))
                          
ggplot(coxM, aes(x=coef, y=termLabels))+
  geom_point()+
  geom_errorbar(aes(xmin=LCI, xmax=UCI, width=0.2))+
  theme_bw(base_size = 12)+
  xlab('Hazard (+/- 95% CI)')+
  ylab('')+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14))+
  geom_vline(xintercept=0, linetype='dashed')



## try depicting these effects:

library(survminer)
# a sex-dependent effect of recruitment age on mortality

# divide recruitment ages into bins of roughly equal sex-representation
ggplot(ldat, aes(x=recruitmentAge, fill=sex)) +
  geom_vline(xintercept=5, linetype='dashed')+
  geom_histogram(alpha=0.8, position = 'identity') +
  labs(fill="")+
  facet_wrap(~sex)+
  theme_bw(base_size = 12)

ldat$recruitmentAge

ldat$recruitmentBin <- ifelse(ldat$recruitmentAge>5, '>5y', '<5y')
table(ldat$recruitmentBin, ldat$sex)

s <- survfit(Surv(tstart, tstop, death==1) ~ sex, ldat)
ggsurvplot_facet(s,  facet.by = 'recruitmentBin', data = ldat, xlab='follow up up time (years)')

## males recruited younger are more likely to live than young females
ldat$recruitmentBin <- cut(ldat$recruitmentAge, 4)
table(ldat$recruitmentBin)
table(ldat$recruitmentBin, ldat$sex)

sM <- survfit(Surv(tstart, tstop, death==1) ~ recruitmentBin, ldat[ldat$sex=='M', ])
sF <- survfit(Surv(tstart, tstop, death==1) ~ recruitmentBin, ldat[ldat$sex=='F', ])

pM <- print(ggsurvplot(sM, risk.table = F, 
           legend = "right", 
           palette ="PuOr",
           ggtheme = theme_bw(), 
           title='male survival',
           xlab='follow up up time (years)'))

pF <- print(ggsurvplot(sF, risk.table = F, 
           legend = "right", 
           palette ="PuOr",
           ggtheme = theme_bw(), 
           title='female survival',
           xlab='follow up up time (years)'))


pM
pF



save(ldat, dat, mzs, statusData, file='data/longitudinalData')

#######################################################################################



rm(list=ls())
############################################################################
### analysis of metabolites and mortality
############################################################################

## two approaches to investigate metabolites or the metabolome that might assoc with mortality:
# univariate, looping through the mzs
# elastic net to pull a predictive set (will bias toward complementary mzs) : 'see enCox.R' - this failed to robustly identify mutlivaiate models of mortality from metabolites

load('data/longitudinalData')

## considering the potential for metabolites to associate with mortality in different contexts, run stepwise model selection:

# omit fawns
x <- dplyr::select(dat, c(age, field_year, deer))
x <- x[x$age>1, ]
table(paste(ldat$deer, ldat$field_year) %in% paste(x$deer, x$field_year)) # ugly but functional

ldat <- ldat[paste(ldat$deer, ldat$field_year) %in% paste(x$deer, x$field_year), ] # omit data from fawns

# scale
scat <- ldat
scat$sex <- scale(as.numeric(scat$sex))
scat$recruitmentAge <- scale(scat$recruitmentAge)
scat$body_mass <- scale(scat$body_mass)
scat[ ,mzs] <- scale(scat[ ,mzs])

# loop through the metabolites:
sList <- list()

for(i in 1:length(mzs)) {
    sList[[i]] <- summary(step(coxph(Surv(tstart, tstop, death==1) ~ ldat[ ,mzs[i]] * recruitmentAge + ldat[ ,mzs[i]] * body_mass + ldat[ ,mzs[i]] * sex, scat))) }

names(sList) <- mzs

sList[[1]]

### save cox models and bonferroni corrected p values 
coxList <- lapply(sList, function(x) data.frame('term'=rownames(x$coefficients), x$coefficients, 'BH'=p.adjust(x$coefficients[ ,5], 'BH', n=72)))
names(coxList) <- mzs

save(coxList, ldat, mzs, file='mz.coxModels.selected')



load('mz.coxModels.selected')
##################################################
coxList[[1]]

o <- bind_rows(coxList, .id = "metabolite") # 'output'
rownames(o) = NULL

m <- o[grepl('mzs', o$term), ] # examine metabolite-involved effects
head(m)
m[order(m$`Pr...z..`), ][1:10, ]


## summarize results with Hazard plot(s)
m$mzLabels <- m$metabolite
m$mzLabels[m$BH>0.05] <- ''

termNames <- print(m$term)

termLabels <- termNames
termLabels <- gsub("\\]", '', termLabels)
termLabels <- gsub("\\[", '', termLabels)
termLabels <- gsub('ldat, mzsi', 'metabolite', termLabels)
termLabels <- gsub(':', ' x ', termLabels)
termLabels <- gsub('sex', 'sex (male)', termLabels)
termLabels <- gsub('recruitmentAge', 'age', termLabels)
termLabels <- gsub('body_mass', 'Body Mass', termLabels)
termLabels <- tolower(termLabels)
termLabels

names(termLabels) <- termNames
m$termLabels <- termLabels

head(m)
m$coef
m$LCI <- m$coef-1.96*m$se.coef.
m$UCI <- m$coef+1.96*m$se.coef.

m$metabolite <- as.factor(m$metabolite)
m$metabolite
m$metabolite <- factor(m$metabolite, levels=rev(levels(m$metabolite)))

m$sig <- as.factor(ifelse(m$BH <0.05, paste0("\U03B1", '<0.05'), 'NS'))

ggplot(m, aes(x=coef, y=metabolite, color=sig))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbar(aes(xmin=LCI, xmax=UCI, width=0.2))+
  theme_classic(base_size = 12)+
  xlab('Hazard (+/- 95% CI)')+
  facet_wrap(~termLabels, nrow=1)+
  scale_color_manual(values=c('grey70','purple'))+
  ylab('')+
  theme(axis.text=element_text(size=10))+
  labs(color='')

ggplot(subset(m, m$sig!='NS'), aes(x=coef, y=metabolite))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbar(aes(xmin=LCI, xmax=UCI, width=0.2))+
  theme_bw(base_size = 14)+
  xlab('Hazard (+/- 95% CI)')+
  facet_wrap(~termLabels, nrow=1)+
  ylab('')+
  theme(axis.text=element_text(size=12))+
  labs(color='')

range(dat$date)

head(m)

ggplot(m, aes(y=-log(Pr...z.., 10), x=coef, label=mzLabels, color=BH<0.05))+
  geom_point(size=1)+
  ylab('-log(10) P')+
  theme_classic(base_size = 16)+
  ggrepel::geom_text_repel(max.overlaps = 10, size=3, nudge_x = 4)+
  scale_color_manual(values=c('grey', 'red'))+
  facet_wrap(~termLabels, nrow=1)+
  theme(legend.position = 'bottom')+
  xlab(expression(beta))+
  theme(axis.text = element_text(size = 10))+
  labs(color='Benjamini-Hochberg=0.05')


m[m$BH <=0.05, ] 

### try to depict these effects
coxList[["adenosine"]]
  
ldat$adenosine.grp <- cut_number(ldat$`adenosine`, 4)
ldat$BodyMass <- cut_number(ldat$body_mass, 2)

table('body mass'=ldat$BodyMass, ldat$adenosine.grp)

s <- survfit(Surv(tstart, tstop, death==1) ~ adenosine.grp, ldat)

ggsurvplot(s,  
           facet.by = 'BodyMass', 
           data = ldat, 
           palette = c("orange", "blue", "purple", 'red'),
           ggtheme = theme_bw(), 
           title='adenosine',
           xlab='follow up up time (years)')


mzs

coxList[["N-N-dimethyl-arginine"]]
ldat$dimethylarginine.grp <- cut_number(ldat$`N-N-dimethyl-arginine`, 3)
ldat$BodyMass <- cut_number(ldat$body_mass, 3)

table(ldat$BodyMass)

table('body mass'=ldat$BodyMass, ldat$dimethylarginine.grp)

s <- survfit(Surv(tstart, tstop, death==1) ~ dimethylarginine.grp, ldat)

ggsurvplot(s,  
           facet.by = 'BodyMass', 
           data = ldat, 
           palette = c("orange", "blue", "purple"),
           ggtheme = theme_bw(base_size = 12),
           xlab='follow up up time (years)')


ggplot(ldat, aes(`N-N-dimethyl-arginine`, x=dimethylarginine.grp, fill=dimethylarginine.grp, color=dimethylarginine.grp))+
  geom_violin(alpha=0.5)+
  geom_jitter(width=0.1, size=1)+
  scale_fill_manual(values=c("orange", "blue", "purple"))+
  scale_color_manual(values=c("orange", "blue", "purple"))+
  theme_bw(base_size = 14)+
  xlab('')+
  labs(fill='N-N-dimethyl-arginine', color='N-N-dimethyl-arginine')+
  guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))+
  theme(axis.text.x = element_text(angle = 40, vjust = 0.6, hjust=0.5))


coxList[["nicotinamide"]] # the three way interaction mz x sex, and mz by mass will likely be very difficult to depict in KM plots

sList[["nicotinamide"]]


aa <- aareg(Surv(tstart, tstop, death==1) ~ nicotinamide * recruitmentAge + nicotinamide * body_mass + nicotinamide * sex, scat)

aa$coefficient <- aa$coefficient[ ,-c(1)] # removes the intercept from the plot panels
colnames(aa$coefficient)

colnames(aa$coefficient) <- paste(1:ncol(aa$coefficient), colnames(aa$coefficient))
colnames(aa$coefficient) <- gsub('_', ' ', colnames(aa$coefficient)) 
colnames(aa$coefficient) <- gsub('\\:', ' x ', colnames(aa$coefficient)) 
colnames(aa$coefficient) <- gsub('recruitmentAge', 'age', colnames(aa$coefficient)) 
colnames(aa$coefficient) <- gsub('sex', 'sex (male)', colnames(aa$coefficient)) 

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

autoplot(aa) +
  theme_classic(base_size = 14)+
  geom_hline(yintercept=0, alpha=0.4, linetype='dashed')+
  theme(legend.position = 'none')+
  facet_wrap(~variable, scales='free', nrow=2)+
  xlab('follow-up time (years)')+
  ylab('additive hazard')



coxList[["hippurate"]]

ldat$hippurate.grp <- cut_number(ldat$`hippurate`, 3)

table(ldat$hippurate.grp, ldat$sex)

s <- survfit(Surv(tstart, tstop, death==1) ~ hippurate.grp, ldat)

ggsurvplot(s,  
           facet.by = 'sex', 
           data = ldat, 
           palette = c("orange", "blue", "purple"),
           ggtheme = theme_bw(14), 
           xlab='follow up up time (years)')


ggplot(ldat, aes(`hippurate`, x=hippurate.grp, fill=hippurate.grp, color=hippurate.grp))+
  geom_violin(alpha=0.5)+
  geom_jitter(width=0.1, size=1)+
  xlab('')+
  scale_fill_manual(values=c("orange", "blue", "purple"))+
  scale_color_manual(values=c("orange", "blue", "purple"))+
  theme_bw(base_size = 14)+
  labs(fill='hippurate', color='hippurate')+
  guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))+
  theme(axis.text.x = element_text(angle = 40, vjust = 0.6, hjust=0.5))



