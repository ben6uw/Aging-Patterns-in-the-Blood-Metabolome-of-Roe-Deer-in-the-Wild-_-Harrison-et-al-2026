
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
library(performance)


# analysis of roe deer data

# Jeff initially sent data for all >1yr old deer, and then sent the data for all deer that were sampled:  This is a message that accompanied the second send:
# Please find attached the updated dataset included the full set of samples (with samples from 1 year old individuals). As you will see, there is only life history / physiological data corresponding to 209 individuals, not 210.

# The reason for this is that the individual F338 has been captured twice in winter 2017 with blood samples (samples 2017-51-020 and 2017-51-048) collected each time (this normally does not happen). When sending the samples for the quantification of the metabolome, we made a mistake and we sent the blood samples corresponding to these two captures (they correspond to your samples 103 and 130). 

# As the best is probably to keep one sample per individual per year, we should remove one of these two samples. I thus only kept data for the sample collected during the second capture of the individuals (which corresponds to your sample 103). I hope it makes sense.


dir('~/Google Drive/')
getwd()
setwd('~/Google Drive/My Drive/Documents/RoeDeer')


rm(list=ls())
######################################################################
### PCA to evaluate sources of metabolome  variation
######################################################################
load('data/dataforUnivariateModels')

pca <- prcomp(d[ ,mzs], scale=T)
propvar <- pca$sdev/sum(pca$sdev)
pv <- data.frame(PC=paste0('PC', 1:length(mzs)), propvar)
pv$PC <- factor(pv$PC, levels=pv$PC)
  
## use the Tracy Widom distribution to evalaute how many eigenvals exceed expectations from random data:
eigs <- pca$sdev^2 # extract eigenvalues
sigEigs <- print(sum(RMTstat::ptw(eigs, lower.tail = F)<=0.05)) # the number of PCs which exceed those from random data.
pv$sig <- factor(ifelse(as.numeric(pv$PC) <= sigEigs, 'P<0.05', 'NS'), levels=c('P<0.05', 'NS'))

head(pv)
pv$PC <- as.factor(c(1:nrow(pv)))
pv$cumulativeVar <- cumsum(pv$propvar)


scree.all <- ggplot(pv[1:(sigEigs+4), ], aes(y=propvar, x=PC, fill=sig))+
  geom_col()+
  theme_classic(base_size = 14)+
  scale_fill_manual(values=c('purple', 'grey'))+
  labs(fill='')+
  ylab('proportion of variance')+
  ggtitle('metabolome PCA all 209 samples')

scree.all 

ggplot(pv[1:(sigEigs), ], aes(y=cumulativeVar, x=PC))+
  geom_col()+
  theme_classic(base_size = 14)+
  ylab('cumulative proportion of variance')+
  ggtitle('metabolome PCA all 209 samples')


PCs <- paste0('PC', 1:sigEigs)
propvar[1:sigEigs]
cumulativeVar <- print(sum(propvar[1:sigEigs])) # the top _ PCs explain _% of the metabolome variance 

pca <- data.frame(d[ ,!colnames(d) %in% c(mzs)], pca$x[ ,PCs])
head(pca)
l <- as.data.frame(pivot_longer(pca, PCs, names_to = 'PC'))

# consider technical/environmental/biological effects:
# batch (had been removed via Combat)
# year of sampling (? biological/environmental or technical?)
# biological confounds (AM/PM, body_mass, 'Cohort quality', stage) # CBC will be considered later
# [note]: AM/PM was removed


## use model selection to identify the most substantial sources of variation within each PC

stepLmList <- list()
segModList <- list()
pList <- list() 


for(i in 1:length(PCs)) {
  stepLmList[[i]] <- step(lm(value ~ (age + body_mass + sex)^2 + stage + cohort_quality + field_year + date, l[l$PC==PCs[i], ]), 
            direction = "both",
            scope = list(lower = ~ age))

  segModList[[i]] <- tryCatch(
    segmented(stepLmList[[i]], seg.Z = ~age),
    error = function(e) {
      message("Skipping item ", i, " because of error: ", e$message)
      return(NA) # Assign NA to the results list and move on
      })
    pList[[i]] <- tryCatch(
      davies.test(segModList[[i]], k=10)$p.value,
      error = function(e) {
        message("Skipping item ", i, " because of error: ", e$message)
        return(NA) # Assign NA to the results list and move on
        })
}


## screen models for: 
# those that either did not have a significant segmented regression, by delta BIC and/or P value:
segEval <- data.frame('PC'=PCs, 'davies_P'= unlist(pList))
nonSegModels <- print(which(segEval[ ,'davies_P'] > 0.05)) # Davies test P >0.05

## or whose segmentation was not identified (failure to identify breakpoint)
failedPCs <- print(which(sapply(segModList, function(x) is.na(x[1])))) # identifies any failed PC models, there don't have to be any

modelsToReplace <- union(nonSegModels, failedPCs) 
modelsToReplace <- print(modelsToReplace[order(modelsToReplace)]) 

## for any segmented age model that failed the above criteria: 
# repeat the stepwise linear regression, but allow the model to drop age, which is no longer evaluated for segmentation:

for(i in modelsToReplace) {
  segModList[[i]]  <- step(lm(value ~ (age + body_mass + sex)^2 + stage + cohort_quality + field_year + date, l[l$PC==PCs[i], ]), direction = "both") 
}
  
######## incorporate random effect, so that we can avoid pseudo replication
# incorporate random effect of deer
# fit the mixed model and test terms
mmList <- list()

for(i in 1:length(PCs)) {
  tmp <- dplyr::select(l[l$PC==PCs[i], ], c(value, age, body_mass, sex, stage, cohort_quality, field_year, date, deer))
  tryCatch( expr = {
    f0 <- formula(segModList[[i]]) # extract selected model 
    f1 <- update(f0, ~ . + (1 | deer)) # add random effect  
    mmList[[i]] <- lmer(f1, data=tmp) # fit mixed model
  },
  error = function(e) {
    return(NA)
  })
}

names(mmList) <- PCs
mmList[sapply(mmList, is.null)] <- segModList[sapply(mmList, is.null)]

aList <- lapply(mmList, anova) # ANOVA 
aList[[1]]


# Get marginal and conditional R2 (performance package)
r2s <- t(sapply(mmList, r2)) # Nakagawa's R2
head(r2s)
r2s <- bind_cols('R2_marginal'=unlist(r2s[ ,2]), 'R2_conditional'=unlist(r2s[ ,1]))
r2s$ICC <- r2s$R2_conditional - r2s$R2_marginal ## The contribution of random effects: In simple cases, the ICC corresponds to the difference between the conditional R2 and the marginal R2
head(r2s)

models <- sapply(mmList, function(x) as.character(formula(x))[3])
PCmodels <- data.frame(PC=PCs, propVar = pv$propvar[1:sigEigs], model=models, r2s)
rownames(PCmodels) = NULL
head(PCmodels)

write.csv(PCmodels, file='paper/Table S1 PC model selection.csv', row.names = F, quote=F)

## can I write script to identify the significant terms in each model to add to the table?:

# variance by term within PCs:
aList[[1]]

for(i in 1:length(aList)){
aList[[i]]$term <- as.character(rownames(aList[[i]])) 
aList[[i]]$FDR <- p.adjust(aList[[i]]$`Pr(>F)`, n=sigEigs, 'BH')
}

x <- as.data.frame(as.matrix(bind_rows(aList, .id = "PC")))
x$`Pr(>F)` <- as.numeric(x$`Pr(>F)`)
x$`Sum Sq` <- as.numeric(x$`Sum Sq`)
x$`Pr(>F)` <- as.numeric(x$`Pr(>F)`)
x$FDR <- as.numeric(x$FDR)
x <- x[!grepl('Resid', x$term), ]
rownames(x) <- NULL

head(x)
sort(table(x$term))
x$PC <- factor(x$PC, levels=PCs)

ggplot(x, aes(`Sum Sq`, x=PC, fill=term))+
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 12)+
  facet_wrap(~term, nrow=2)

ggplot(subset(x, PC %in% PCs[1:6]), aes(y=`Sum Sq`, x=term))+
  geom_bar(position="stack", stat="identity")+
  theme_bw(base_size = 14)+
  facet_wrap(~PC, nrow=2)+
  xlab('')+
  ylab('variance (sum of squares)')+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))

aList[[2]]

ggplot(pca, aes(PC2, x=field_year))+
  geom_violin()+
  theme_bw(base_size = 14)


# mark significant terms:
head(x)
x$sig <- 'ns'
x$sig[x$FDR < 0.05] <- '*'
x$sig[x$FDR < 0.005] <- '**'
x$sig[x$FDR < 0.0005] <- '***'

table(x$sig)
subset(x, PC %in% PCs[1:6])

## plot-friendly term names
pnames <- x$term
pnames <- gsub('body_', 'body ', pnames)
pnames <- gsub('\\:', ' x ', pnames)
x$pnames <- pnames

# Plot with asterisks
ggplot(subset(x, PC %in% PCs[1:6]), aes(y=`Sum Sq`, x=pnames))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = sig), vjust = -0.1, size = 5)+
  theme_bw(base_size = 14)+
  facet_wrap(~PC, nrow=2)+
  xlab('')+
  ylim(0, 400)+
  ylab('sum of squares')+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.background = element_blank())

ggplot(subset(x, PC %in% PCs[1:6]), aes(y=`Sum Sq`, x=pnames))+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label = sig), vjust = -0.1, size = 5)+
  theme_bw(base_size = 14)+
  facet_wrap(~PC, nrow=2)+
  xlab('')+
  ylim(0, 400)+
  ylab('sum of squares')+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1),
        strip.background = element_blank())

sort(table(x$term[x$FDR<0.05]))

head(x)
view(x[x$sig!='ns', c(1, 8)])

ggplot(pca, aes(PC3, x=age))+
  geom_point(size=1)+
  theme_bw(base_size = 14)

ggplot(pca, aes(PC1, x=sex))+
  geom_violin()+
  theme_bw(base_size = 14)


rm(list=ls())
######################################################################
### analyze adult (only) data
######################################################################

load('data/dataforUnivariateModels')

## how many fawn samples remain in the data:
table(d$stage)

d2 <- d[d$age>1, ] # omit fawns for this analysis
length(unique(d2$deer)) # 84 deer remain
nrow(d2) # among 134 samples

### repeat the PCA to see what influence fawn/maternal signal may have on hiding effects of age, body size, sex, etc
pca <- prcomp(d2[ ,mzs], scale=T)

propvar <- pca$sdev/sum(pca$sdev)
pv <- data.frame(PC=1:length(mzs), propvar)
head(pv)
pv$PC <- factor(pv$PC, levels=pv$PC)

pv$cumulativeVar <- cumsum(pv$propvar)

## use the Tracy Widom distribution to evalaute how many eigenvals exceed expectations from random data:
eigs <- pca$sdev^2 # extract eigenvalues
sigEigs <- print(sum(RMTstat::ptw(eigs, lower.tail = F)<=0.05)) # the number of PCs which exceed the 

pv$sig <- factor(ifelse(as.numeric(pv$PC) <= sigEigs, 'P<0.05', 'NS'), levels=c('P<0.05', 'NS'))


## compare the eigenvalues of PCs from the metabolome with (d) and without (d2) fawn data
eigComp <- data.frame(PC = pv$PC, eig_alldeer = prcomp(d[ ,mzs], scale=T)$sdev^2, eig_nofawns = prcomp(d2[ ,mzs], scale=T)$sdev^2)
eigComp$eig_alldeer
eigComp$eig_nofawns
eigComp$sig = pv$sig

eigCompPlot <- ggplot(subset(eigComp, eigComp$PC%in% 1:(sigEigs+4)), aes(x=eig_alldeer, y=eig_nofawns, color=sig))+
  geom_abline(slope=1, alpha=0.5)+
  geom_point(size=1.5)+
  xlab('eigenvalues all deer')+
  ylab('eigenvalues without fawns')+
  theme_bw(base_size = 14)+
  scale_color_manual(values=c('purple', 'grey'))+
  scale_x_log10(breaks=seq(from=1,to=13, by=2))+
  scale_y_log10(breaks=seq(from=1,to=13, by=2))+
  labs(color='')+ 
  theme(axis.text = element_text(size = 10))

eigCompPlot

cumVar <- data.frame('PC' = 1:19,
  allDeer = cumsum(prcomp(d[ ,mzs], scale=T)$sdev/sum(prcomp(d[ ,mzs], scale=T)$sdev))[1:19], 
  noFawn = cumsum(prcomp(d2[ ,mzs], scale=T)$sdev/sum(prcomp(d2[ ,mzs], scale=T)$sdev))[1:19])

head(cumVar)
l <- pivot_longer(cumVar, c('allDeer', 'noFawn'), values_to = 'cumulativeVar', names_to='data')
l$data <- ifelse(l$data=='allDeer', 'all deer', 'adult only')

ggplot(l, aes(y=cumulativeVar, x=PC, group=data, color=data))+
  geom_line(size=1.5)+
  theme_classic(base_size = 14)+
  labs(color='')+
  scale_x_continuous(breaks=seq(from=1, to=sigEigs, by=2))

propvar <- pca$sdev/sum(pca$sdev)
pv <- data.frame(PC=1:length(mzs), propvar)
head(pv)
pv$PC <- factor(pv$PC, levels=pv$PC)

pv$cumulativeVar <- cumsum(pv$propvar)


PCs <- paste0('PC', 1:sigEigs)
propvar[1:sigEigs]
cumulativeVar <- print(sum(propvar[1:sigEigs])) # the top _ PCs explain _% of the metabolome variance

