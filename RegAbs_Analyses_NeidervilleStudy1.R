## Hayoung Ahn & Erik Nook
## Regulation Abstraction analyses of Neiderville Study 1

##----------------Load packages-------------------------
if (!require(coin)) {install.packages("coin"); require(coin)} ## For permutation tests
if (!(require("gclus"))) install.packages("gclus"); require("gclus")  ## For creating scatterplot panels
if (!(require("corrplot"))) install.packages("corrplot"); require("corrplot")  ## For creating correlation panels
if (!require(MBESS)) install.packages("MBESS"); require(MBESS) ## For mediation   
if (!require(MASS)) install.packages("MASS"); require(MASS) ## For mediation   
if (!require(lme4)) install.packages("lme4"); require(lme4) ## For mixed effect model 
if (!require(lmerTest)) install.packages("lmerTest"); require(lmerTest) ## For mixed effect model 
if (!require(lattice)) install.packages("lattice"); require(lattice) ## for dotplot
if (!require(effects)) install.packages("effects"); require(effects) ## for allEffects plot
if (!require(psych)) install.packages("psych"); require(psych) ## for allEffects plot
if (!require(coefplot2)) install.packages("coefplot2",repos="http://www.math.mcmaster.ca/bolker/R",type="source"); require(coefplot2) ## for allEffects plot
if (!require(car)) {install.packages("car"); require(car)}
if (!require(stringi)) {install.packages("stringi"); require("stringi")}    
if (!require(stringr)) {install.packages("stringr"); require("stringr")}     
if (!require(reshape2)) {install.packages("reshape2"); require("reshape2")}      
if (!require(ez)) {install.packages("ez"); require("ez")}     
if (!require(effsize)) {install.packages("effsize"); require("effsize")}
if (!require(ggplot2)) {install.packages("ggplot2"); require("ggplot2")}
if (!require(Rmisc)) {install.packages("Rmisc"); require("Rmisc")}
if (!require(lm.beta)) {install.packages("lm.beta"); require("lm.beta")}
if (!require(lmSupport)) {install.packages("lmSupport"); require("lmSupport")}
if (!require(textclean)) {install.packages("textclean"); require("textclean")}
if (!require(tm)) {install.packages("tm"); require("tm")}
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load_gh('hrbrmstr/pluralize')
p_load(quanteda)

##----------------Prepare and attach dataframe---------------------------

## Load data
setwd('/Users/hahn/Dropbox/codingclub')
d <- read.csv('Neiderville_Study1_OASIS.csv')

## Produce important variables
d$CB <- ifelse(is.na(d$x), 2, 1)
N <- length(d$CB)
d$ID <- factor(1:N)

## Produce important factors
d$CB <- as.factor(d$CB)
d$Gender <- factor(d$Gender, levels = c(1,2), labels = c('M','F'))
d$Race <- as.factor(d$Race)
d$Income <- as.factor(d$Income)

## Cancel out Ps who did not indicate income level or race
d$Income[which(d$Income=='10')] <- NA
d$Race[which(d$Race=='11')] <- NA

## Unrandomize ratings
d$order1 <- gsub('R1.[0-9][0-9]','', d$DO.BL.ERegTaskCB1)
d$order2 <- gsub('R2.[0-9][0-9]','', d$DO.BL.ERegTaskCB2)

for (s in 1:N) {
  d$order[s] <- max(cbind(d$order1[s],d$order2[s]))
}

d$order <- gsub('  ',' ', d$order)
d$order <- gsub(' 2 ',' ', d$order)

rates <- matrix(rep(NA,60),nrow = 1)
ratecolnames <- c(paste("LookNeg0",1:9,sep=""), paste("LookNeg",10:20,sep=""), paste("ReapNeg0",1:9,sep=""), paste("ReapNeg",10:20,sep=""),paste("LookNeu0",1:9,sep=""), paste("LookNeu",10:20,sep=""))
colnames(rates) <- ratecolnames
rates <- rates[,order(colnames(rates))]

getvars <- names(d) %in% c(paste("R1.0",1:9,sep=""),paste("R1.",10:60,sep=""), paste("R2.0",1:9,sep=""),paste("R2.",10:60,sep=""))
ratecols <- d[getvars]

ratevec <- NULL
tnames <- NULL
for (s in 1:N) {
  c <- 1
  pics <- strsplit(d$order[s], " ")
  for (trial in 1:60) {
    picture <- pics[[1]][trial]
    if (substr(picture,3,3)=="1" || substr(picture,3,3) == "2") {
      if (substr(picture,3,3)=="1") {
        cond <- "LookNeg"} else {cond = "ReapNeg"}
    }  else cond <- "LookNeu"
    tnames[trial] <- paste(cond,substr(picture,5,6), sep = "")
    c <- c+7
  }
  
  if (substr(picture,1,1) == "1") {
    ratevec <- ratecols[s,1:60]
  } else {ratevec <- ratecols[s,61:120]}
  
  colnames(ratevec) <- tnames
  ratevec <- ratevec[,order(colnames(ratevec))]
  rates <- rbind(rates, ratevec)
}
rates <- rates[-1,]
d <- cbind(d,rates)

## Compute EReg Variables
getvars <- names(d) %in% c(paste("LookNeg0",1:9,sep=""), paste("LookNeg",10:20,sep=""))
d$LookNeg <- rowMeans(d[getvars], na.rm = T)

getvars <- names(d) %in% c(paste("ReapNeg0",1:9,sep=""), paste("ReapNeg",10:20,sep=""))
d$ReapNeg <- rowMeans(d[getvars], na.rm = T)

getvars <- names(d) %in% c(paste("LookNeu0",1:9,sep=""), paste("LookNeu",10:20,sep=""))
d$LookNeu <- rowMeans(d[getvars], na.rm = T)

d$ReapSuccess <- d$LookNeg - d$ReapNeg
d$Reactivity <- d$LookNeg - d$LookNeu

## Check timing
d$minTime <- NA
d$maxTime <- NA
d$meanTime <- NA
for (s in 1:N) {
  d$minTime[s] <- min(d[s,c("t_3",paste("t_3.",1:119,sep=""))],na.rm=T)
  d$maxTime[s] <- max(d[s,c("t_3",paste("t_3.",1:119,sep=""))],na.rm=T)
  d$meanTime[s] <- mean(t(d[s,c("t_3",paste("t_3.",1:119,sep=""))]),na.rm=T)
}
d$badtime <- ifelse(d$minTime < 20,1,0)
d$badtime <- as.factor(d$badtime)
d$badtime25 <- ifelse(d$minTime < 25,1,0)
d$badtime25 <- as.factor(d$badtime25)
d$badtime15 <- ifelse(d$minTime < 15,1,0)
d$badtime15 <- as.factor(d$badtime15)
d$badmeantime25 <- ifelse(d$meanTime < 25,1,0)
d$badmeantime25 <- as.factor(d$badmeantime25)
d$badmeantime28 <- ifelse(d$meanTime < 28,1,0)
d$badmeantime28 <- as.factor(d$badmeantime28)
d$badmeantime <- ifelse(d$meanTime < 20,1,0)
d$badmeantime <- as.factor(d$badmeantime)

## Produce LIWC regulation variables (difference between LookNeg and ReapNeg)
d$WC_reg <- d$WC_RNeg - d$WC_LNeg
d$pronoun_reg <- d$pronoun_RNeg - d$pronoun_LNeg
d$i_reg <- d$i_RNeg - d$i_LNeg
d$we_reg <- d$we_RNeg - d$we_LNeg
d$you_reg <- d$you_RNeg - d$you_LNeg
d$verb_reg <- d$verb_RNeg - d$verb_LNeg
d$past_reg <- d$past_RNeg - d$past_LNeg
d$present_reg <- d$present_RNeg - d$present_LNeg
d$future_reg <- d$future_RNeg - d$future_LNeg
d$social_reg <- d$social_RNeg - d$social_LNeg
d$affect_reg <- d$affect_RNeg - d$affect_LNeg
d$posemo_reg <- d$posemo_RNeg - d$posemo_LNeg
d$negemo_reg <- d$negemo_RNeg - d$negemo_LNeg
d$anx_reg <- d$anx_RNeg - d$anx_LNeg
d$anger_reg <- d$anger_RNeg - d$anger_LNeg
d$sad_reg <- d$sad_RNeg - d$sad_LNeg
d$cogmech_reg <- d$cogmech_RNeg - d$cogmech_LNeg
d$insight_reg <- d$insight_RNeg - d$insight_LNeg
d$cause_reg <- d$cause_RNeg - d$cause_LNeg
d$tentat_reg <- d$tentat_RNeg - d$tentat_LNeg
d$certain_reg <- d$certain_RNeg - d$certain_LNeg
d$percept_reg <- d$percept_RNeg - d$percept_LNeg
d$body_reg <- d$body_RNeg - d$body_LNeg
d$article_reg <- d$article_RNeg - d$article_LNeg
d$discrep_reg <- d$discrep_RNeg - d$discrep_LNeg
d$Sixltr_reg <- d$Sixltr_RNeg - d$Sixltr_LNeg

d$i_reg_raw <- d$i_RNeg_raw - d$i_LNeg_raw
d$present_reg_raw <- d$present_RNeg_raw - d$present_LNeg_raw
d$posemo_reg_raw <- d$posemo_RNeg_raw - d$posemo_LNeg_raw
d$negemo_reg_raw <- d$negemo_RNeg_raw - d$negemo_LNeg_raw
d$article_reg_raw <- d$article_RNeg_raw - d$article_LNeg_raw
d$discrep_reg_raw <- d$discrep_RNeg_raw - d$discrep_LNeg_raw
d$Sixltr_reg_raw <- d$Sixltr_RNeg_raw - d$Sixltr_LNeg_raw

## produce other datasets
d_nobadmeantime25 <- d[-which(d$badmeantime25==1),]
data <- d_nobadmeantime25[-which(d_nobadmeantime25$BadInstructions==1),]


## Produce Linguistic Psychological Distancing measure (computed at trial level, (see Mehl, Robbins, & Holleran, 2012))
LIWCdata <- read.csv('Neiderville_Study1_LinguisticData.csv')

LIWCdata$ID <- NULL
LIWCdata$Trial <- NULL
LIWCdata$Cond <- NULL
for (i in 1:dim(LIWCdata)[1]) { 
  s <- strsplit(paste(LIWCdata$Filename[i]), "_")  
  LIWCdata$ID[i] <- as.numeric(s[[1]][1])
  s2 <- str_extract_all(paste(s[[1]][2]),"[0-9]")[[1]]
  ifelse(length(s2) == 2, LIWCdata$Trial[i] <-  as.numeric(paste(s2[1],s2[2],sep="")), LIWCdata$Trial[i] <- as.numeric(s2))
  if (LIWCdata$Trial[i] < 41) {ifelse(LIWCdata$Trial[i] < 21, LIWCdata$Cond[i] <- 1, LIWCdata$Cond[i] <- 2)}
  else {LIWCdata$Cond[i] <- 3}
}

excludedSIDs <- setdiff(d$ID, data$ID)
for (i in 1:length(LIWCdata$ID)) {
  LIWCdata$exclude[i] <- sum(LIWCdata$ID[i]==excludedSIDs) > 0}
LIWCdata <- LIWCdata[-which(LIWCdata$exclude==1),]

#Z-score and reverse Z-score! 
LIWCdata$SixltrZ <- scale(LIWCdata$Sixltr)
LIWCdata$articleZ <- scale(LIWCdata$article)
LIWCdata$iZ <- -1*scale(LIWCdata$i)
LIWCdata$discrepZ <- -1*scale(LIWCdata$discrep)
LIWCdata$presentZ <- -1*scale(LIWCdata$present)

alpha_psychdist <- cronbach.alpha(cbind(LIWCdata$SixltrZ, LIWCdata$articleZ, LIWCdata$iZ, LIWCdata$discrepZ, LIWCdata$presentZ)); alpha_psychdist$alpha 

LIWCdata$psychdist <- rowMeans(cbind(LIWCdata$SixltrZ, LIWCdata$articleZ, LIWCdata$iZ, LIWCdata$discrepZ, LIWCdata$presentZ))
LIWCdata$psychdist_no6 <- rowMeans(cbind(LIWCdata$articleZ, LIWCdata$iZ, LIWCdata$discrepZ, LIWCdata$presentZ))

## Paste into large dataset
data$NoCondDiff <- NULL
for (i in data$ID) {
  data$psychdist_LNeg[which(data$ID==i)] <- mean(LIWCdata$psychdist[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==1))])
  data$psychdist_RNeg[which(data$ID==i)] <- mean(LIWCdata$psychdist[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==2))])
  data$psychdist_LNeu[which(data$ID==i)] <- mean(LIWCdata$psychdist[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==3))])
  data$psychdist_no6_LNeg[which(data$ID==i)] <- mean(LIWCdata$psychdist_no6[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==1))])
  data$psychdist_no6_RNeg[which(data$ID==i)] <- mean(LIWCdata$psychdist_no6[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==2))])
  data$psychdist_no6_LNeu[which(data$ID==i)] <- mean(LIWCdata$psychdist_no6[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==3))])
}

data$psychdist_reg <- data$psychdist_RNeg - data$psychdist_LNeg
data$psychdist_no6_reg <- data$psychdist_no6_RNeg - data$psychdist_no6_LNeg

detach(data)
attach(data)

## Produce Linguistic Categorial Model (LCM) Scores 
LIWCdata$LCM <- ((LIWCdata$DAV_count*1)+(LIWCdata$IAV_count*2)+(LIWCdata$SV_count*3)+(LIWCdata$Adjectives_count*4)+(LIWCdata$Nouns_count*5))/(LIWCdata$DAV_count+LIWCdata$IAV_count+LIWCdata$SV_count+LIWCdata$Adjectives_count+LIWCdata$Nouns_count)
min(LIWCdata$LCM, na.rm = TRUE)
max(LIWCdata$LCM, na.rm = TRUE)
sum(is.na(LIWCdata$LCM)) # missing values

## Paste into large dataset; Average across trials in each condition
data$NoCondDiff <- NULL
for (i in data$ID) {
  data$LCM_LNeg[which(data$ID==i)] <- mean(LIWCdata$LCM[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==1))], na.rm = TRUE)
  data$LCM_RNeg[which(data$ID==i)] <- mean(LIWCdata$LCM[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==2))], na.rm = TRUE)
  data$LCM_LNeu[which(data$ID==i)] <- mean(LIWCdata$LCM[intersect(which(LIWCdata$ID==i),which(LIWCdata$Cond==3))], na.rm = TRUE)
}

data$LCM_reg <- data$LCM_RNeg - data$LCM_LNeg

data$LCM_all <- rowMeans(cbind(data$LCM_LNeg, data$LCM_RNeg, data$LCM_LNeu))
data$psychdist_all = rowMeans(cbind(psychdist_LNeg, psychdist_RNeg, psychdist_LNeu))

detach(data)
attach(data)

##----------------Visualize Data------------------------------

## EReg DVs
x11() 
layout(matrix(1:6,nrow=3, byrow=T),heights = c(2,1,1))
par(oma=rep(2,4), mar = rep(2,4))
boxplot(ReapSuccess, main = "ReapSuccess", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(Reactivity, main = "Reactivity Boxplot", col = c("lightblue")); points(mean(Reactivity, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
hist(ReapSuccess, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(Reactivity, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
qqnorm(ReapSuccess, main = " ", col = "blue"); qqline(ReapSuccess, col = "blue")
qqnorm(Reactivity, main = " ", col = "blue"); qqline(Reactivity, col = "blue")
# boxplot(psychdist)

## psychdist
boxplot(psychdist_LNeg, main = "psychdist_LNeg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(psychdist_RNeg, main = "psychdist_RNeg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(psychdist_LNeu, main = "psychdist_LNeu", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(psychdist_reg, main = "psychdist_reg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
hist(psychdist_LNeg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(psychdist_RNeg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(psychdist_LNeu, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(psychdist_reg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
qqnorm(psychdist_LNeg, main = " ", col = "blue"); qqline(psychdist_LNeg, col = "blue")
qqnorm(psychdist_RNeg, main = " ", col = "blue"); qqline(psychdist_RNeg, col = "blue")
qqnorm(psychdist_LNeu, main = " ", col = "blue"); qqline(psychdist_LNeu, col = "blue")
qqnorm(psychdist_reg, main = " ", col = "blue"); qqline(psychdist_reg, col = "blue")

## LCM
boxplot(LCM_LNeg, main = "LCM_LNeg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(LCM_RNeg, main = "LCM_RNeg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(LCM_LNeu, main = "LCM_LNeu", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
boxplot(LCM_reg, main = "LCM_reg", col = c("lightblue")); points(mean(ReapSuccess, na.rm=T), col = c("blue"), lwd = 2, pch = 6)
hist(LCM_LNeg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(LCM_RNeg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(LCM_LNeu, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
hist(LCM_reg, col = "lightblue", ylab = "Frequency", xlab = " ", main = " ")
qqnorm(LCM_LNeg, main = " ", col = "blue"); qqline(LCM_LNeg, col = "blue")
qqnorm(LCM_RNeg, main = " ", col = "blue"); qqline(LCM_RNeg, col = "blue")
qqnorm(LCM_LNeu, main = " ", col = "blue"); qqline(LCM_LNeu, col = "blue")
qqnorm(LCM_reg, main = " ", col = "blue"); qqline(LCM_reg, col = "blue")

##----------------Analyses-------------------------

## Participant info
length(d$ID)  #Total who did task
sum(d$badmeantime25=="1")  #Number who didn't consistently follow timing Ix
sum(d$BadInstructions==1, na.rm=T)  #Number who didn't consistently follow timing Ix, but 1 P is in both columns, so actually 5!
      intersect(d$ID[which(d$badmeantime25=="1")], d$ID[which(d$BadInstructions==1)]) #60 is in both columns!
length(data$ID)  # Total included in analyses
excludedSIDs <- setdiff(d$ID, data$ID)

## Produce Long dataset for ANOVAs
d.long <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("LookNeg","ReapNeg","LookNeu"), variable.name = "Condition", value.name = "Rating")
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("WC_LNeg","WC_RNeg","WC_LNeu"), variable.name = "Condition", value.name = "WC"); d.long$WC <- d.long0$WC
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("pronoun_LNeg","pronoun_RNeg","pronoun_LNeu"), variable.name = "Condition", value.name = "pronoun"); d.long$pronoun <- d.long0$pronoun
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("i_LNeg","i_RNeg","i_LNeu"), variable.name = "Condition", value.name = "i"); d.long$i <- d.long0$i
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("we_LNeg","we_RNeg","we_LNeu"), variable.name = "Condition", value.name = "we"); d.long$we <- d.long0$we
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("you_LNeg","you_RNeg","you_LNeu"), variable.name = "Condition", value.name = "you"); d.long$you <- d.long0$you
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("verb_LNeg","verb_RNeg","verb_LNeu"), variable.name = "Condition", value.name = "verb"); d.long$verb <- d.long0$verb
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("past_LNeg","past_RNeg","past_LNeu"), variable.name = "Condition", value.name = "past"); d.long$past <- d.long0$past
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("present_LNeg","present_RNeg","present_LNeu"), variable.name = "Condition", value.name = "present"); d.long$present <- d.long0$present
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("future_LNeg","future_RNeg","future_LNeu"), variable.name = "Condition", value.name = "future"); d.long$future <- d.long0$future
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("negate_LNeg","negate_RNeg","negate_LNeu"), variable.name = "Condition", value.name = "negate"); d.long$negate <- d.long0$negate
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("social_LNeg","social_RNeg","social_LNeu"), variable.name = "Condition", value.name = "social"); d.long$social <- d.long0$social
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("affect_LNeg","affect_RNeg","affect_LNeu"), variable.name = "Condition", value.name = "affect"); d.long$affect <- d.long0$affect
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("posemo_LNeg","posemo_RNeg","posemo_LNeu"), variable.name = "Condition", value.name = "posemo"); d.long$posemo <- d.long0$posemo
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("negemo_LNeg","negemo_RNeg","negemo_LNeu"), variable.name = "Condition", value.name = "negemo"); d.long$negemo <- d.long0$negemo
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("anx_LNeg","anx_RNeg","anx_LNeu"), variable.name = "Condition", value.name = "anx"); d.long$anx <- d.long0$anx
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("anger_LNeg","anger_RNeg","anger_LNeu"), variable.name = "Condition", value.name = "anger"); d.long$anger <- d.long0$anger
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("sad_LNeg","sad_RNeg","sad_LNeu"), variable.name = "Condition", value.name = "sad"); d.long$sad <- d.long0$sad
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("cogmech_LNeg","cogmech_RNeg","cogmech_LNeu"), variable.name = "Condition", value.name = "cogmech"); d.long$cogmech <- d.long0$cogmech
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("insight_LNeg","insight_RNeg","insight_LNeu"), variable.name = "Condition", value.name = "insight"); d.long$insight <- d.long0$insight
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("cause_LNeg","cause_RNeg","cause_LNeu"), variable.name = "Condition", value.name = "cause"); d.long$cause <- d.long0$cause
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("tentat_LNeg","tentat_RNeg","tentat_LNeu"), variable.name = "Condition", value.name = "tentat"); d.long$tentat <- d.long0$tentat
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("certain_LNeg","certain_RNeg","certain_LNeu"), variable.name = "Condition", value.name = "certain"); d.long$certain <- d.long0$certain
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("percept_LNeg","percept_RNeg","percept_LNeu"), variable.name = "Condition", value.name = "percept"); d.long$percept <- d.long0$percept
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("bio_LNeg","bio_RNeg","bio_LNeu"), variable.name = "Condition", value.name = "bio"); d.long$bio <- d.long0$bio
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("body_LNeg","body_RNeg","body_LNeu"), variable.name = "Condition", value.name = "body"); d.long$body <- d.long0$body
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("ppron_LNeg","ppron_RNeg","ppron_LNeu"), variable.name = "Condition", value.name = "body"); d.long$ppron <- d.long0$ppron
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("article_LNeg","article_RNeg","article_LNeu"), variable.name = "Condition", value.name = "article"); d.long$article <- d.long0$article
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("discrep_LNeg","discrep_RNeg","discrep_LNeu"), variable.name = "Condition", value.name = "discrep"); d.long$discrep <- d.long0$discrep
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("Sixltr_LNeg","Sixltr_RNeg","Sixltr_LNeu"), variable.name = "Condition", value.name = "Sixltr"); d.long$Sixltr <- d.long0$Sixltr
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("psychdist_LNeg","psychdist_RNeg","psychdist_LNeu"), variable.name = "Condition", value.name = "psychdist"); d.long$psychdist <- d.long0$psychdist
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("psychdist_no6_LNeg","psychdist_no6_RNeg","psychdist_no6_LNeu"), variable.name = "Condition", value.name = "psychdist_no6"); d.long$psychdist_no6 <- d.long0$psychdist_no6
d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("LCM_LNeg","LCM_RNeg","LCM_LNeu"), variable.name = "Condition", value.name = "LCM"); d.long$LCM <- d.long0$LCM
# d.long0 <- melt(data, id.vars = c("ID", "CB"), measure.vars = c("BA_LNeg","BA_RNeg","BA_LNeu"), variable.name = "Condition", value.name = "BA"); d.long$BA <- d.long0$BA

# Participant Characteristics
mean(meanTime); sd(meanTime)
table(Gender)
sum(data$Gender=="M")/length(data$ID)
sum(data$Race==5, na.rm = T)/length(data$ID)
sum(is.na(data$Race)) #Didn't disclose
min(data$Age); max(data$Age)
mean(data$Age); sd(data$Age)
data$Income
sum(is.na(data$Income)) #Didn't disclose
mean(as.numeric(data$Income),na.rm = T)
table(CB)

## Q1: Does one’s language become more abstract when regulating?
# LCM analyses
LookNegdv <- LCM_LNeg
ReapNegdv <- LCM_RNeg
LookNeudv <- LCM_LNeu

M1 <- ezANOVA(d.long, dv = LCM, wid = ID, within = .(Condition), detailed = T); M1; M1$ANOVA$eta <- M1$ANOVA$SSn/(M1$ANOVA$SSn + M1$ANOVA$SSd); 
Lims <- conf.limits.ncf(F.value = M1$ANOVA$F[2], conf.level = 0.90, df.1 <- M1$ANOVA$DFn[2], df.2 <- M1$ANOVA$DFd[2]); Lower.lim <- Lims$Lower.Limit/(Lims$Lower.Limit + df.1 + df.2 + 1); Upper.lim <- Lims$Upper.Limit/(Lims$Upper.Limit + df.1 + df.2 + 1); paste('partial eta^2:',round(M1$ANOVA$eta[2],2)); paste('Lower 90% lim:', round(Lower.lim,2)); paste('Upper 90% lim:', round(Upper.lim,2))

matrix(round(c(mean(LookNegdv), mean(ReapNegdv), mean(LookNeudv), sd(LookNegdv), sd(ReapNegdv), sd(LookNeudv)),2), byrow=T, nrow=2, dimnames = list(c("mean","sd"),c("LookNeg","ReapNeg","LookNeu")))

t.test(ReapNegdv, LookNegdv, paired = T); 
effsize::cohen.d(ReapNegdv, LookNegdv, paired = T)
paste("Hedges g_av:", round((mean(ReapNegdv)-mean(LookNegdv))/sqrt(((sd(LookNegdv)^2+sd(ReapNegdv)^2)/2))*(1-(3/(4*(length(CB)-1)-1))),2))
t.test(ReapNegdv, LookNeudv, paired = T)
effsize::cohen.d(ReapNegdv, LookNeudv, paired = T)
paste("Hedges g_av:", round((mean(ReapNegdv)-mean(LookNeudv))/sqrt(((sd(LookNeudv)^2+sd(ReapNegdv)^2)/2))*(1-(3/(4*(length(CB)-1)-1))),2))

plot(c(mean(LookNegdv), mean(ReapNegdv), mean(LookNeudv)), ylim = c(3.2,3.8))

means <- summarySEwithin(d.long, measurevar="LCM", withinvar="Condition", idvar="ID")
ggplot(means, aes(x=Condition, y=LCM)) +
  geom_bar(position=position_dodge(.9), colour="black", stat="identity") +
  geom_errorbar(position=position_dodge(.9), width=.25, aes(ymin=LCM-ci, ymax=LCM+ci)) +
  # scale_fill_manual(values=c("#CCCCCC","#FFFFFF")) +
  ggtitle("LCM Scores by Condition") + theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw() +
  coord_cartesian(ylim = c(3,4))

## Q2: What is the relationship between how abstract and how distant one’s language is?

# LCM and psychdist analyses (comparing LCM average of 3 conditions with psychdist average)
lm(psychdist_all ~ LCM_all, data = data)
corrout <- corr.test(data$psychdist_all, data$LCM_all)
corrout$p

ggplot(data, aes(x=psychdist_all, y=LCM_all)) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste("Distancing")), y="Abstraction") + 
  ggtitle("Abstraction vs. Distancing Overall") + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  

# LCM and psychdist analyses (across each emotion regulation condition)

## Abstraction vs. Distancing for LookNeg
lm(psychdist_LNeg ~ LCM_LNeg, data = data)
corrout <- corr.test(data$psychdist_LNeg, data$LCM_LNeg)
corrout$p
corrout$r

ggplot(data, aes(x=psychdist_LNeg, y=LCM_LNeg)) + xlim(-0.8,1) + ylim(3,4.5) + 
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste("Distancing")), y="Abstraction") + 
  ggtitle("Abstraction vs. Distancing for Look Neg") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))

## Abstraction vs. Distancing for Reap Neg
lm(psychdist_RNeg ~ LCM_RNeg, data = data)
corrout <- corr.test(data$psychdist_RNeg, data$LCM_RNeg)
corrout$p
corrout$r

ggplot(data, aes(x=psychdist_RNeg, y=LCM_RNeg)) + xlim(-0.8,1) + ylim(3,4.5) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste("Distancing")), y="Abstraction") + 
  ggtitle("Abstraction vs. Distancing for Reap Neg") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))

## Abstraction vs. Distancing for LookNeu
lm(psychdist_LNeu ~ LCM_LNeu, data = data)
corrout <- corr.test(data$psychdist_LNeu, data$LCM_LNeu)
corrout$p
corrout$r

ggplot(data, aes(x=psychdist_LNeu, y=LCM_LNeu)) + xlim(-0.8,1) + ylim(3,4.5) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste("Distancing")), y="Abstraction") + 
  ggtitle("Abstraction vs. Distancing for Look Neu") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  

# Mixed effect model
M1 <- lmer(psychdist ~ Condition * LCM  + (1|ID), data = d.long)
summary(M1)
anova(M1)

# Steiger tests
corout <- corr.test(data[, c("LCM_LNeg","LCM_RNeg", "LCM_LNeu","psychdist_LNeg","psychdist_RNeg", "psychdist_LNeu")])
corout$r

r.test(n = length(data$ID), r12 = corout$r[4,1], r34 = corout$r[5,2], r13=corout$r[2,1], r23=corout$r[4,2],r14=corout$r[5,1],r24=corout$r[5,4]) #LNeg v RNeg    1=LCM_LNeg, 2 = psychdist_LNeg, 3 = LCM_RNeg, 4 = psychdist_RNeg
r.test(n = length(data$ID), r12 = corout$r[4,1], r34 = corout$r[6,3], r13=corout$r[3,1], r23=corout$r[4,3],r14=corout$r[6,1],r24=corout$r[6,4]) #LNeg v LNeu    1=LCM_LNeg, 2 = psychdist_LNeg, 3 = LCM_LNeu, 4 = psychdist_LNeu
r.test(n = length(data$ID), r12 = corout$r[5,2], r34 = corout$r[6,3], r13=corout$r[3,2], r23=corout$r[5,3],r14=corout$r[6,2],r24=corout$r[6,5]) #RNeg v LNeu    1=LCM_RNeg, 2 = psychdist_RNeg, 3 = LCM_LNeu, 4 = psychdist_LNeu

# Analyses for change in LCM and change in psychdist
data$psychdist_reg <- data$psychdist_RNeg - data$psychdist_LNeg
data$LCM_reg <- data$LCM_RNeg - data$LCM_LNeg

lm(psychdist_reg ~ LCM_reg, data = data)
corrout <- corr.test(data$psychdist_reg, data$LCM_reg)
corrout$p
corrout$r

ggplot(data, aes(x=psychdist_reg, y=LCM_reg)) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste(Delta,"Distancing")), y=expression(paste(Delta,"Abstraction"))) + 
  ggtitle("Change in Abstraction vs. Change in Distancing") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  

# Analyses for change in LCM and change in negative affect words
data$negemo_reg <- data$negemo_RNeg - data$negemo_LNeg
data$LCM_reg <- data$LCM_RNeg - data$LCM_LNeg

lm(negemo_reg ~ LCM_reg, data = data)
corrout <- corr.test(data$negemo_reg, data$LCM_reg)
corrout$p
corrout$r

ggplot(data, aes(x=negemo_reg, y=LCM_reg)) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste(Delta,"Negative Affect Words")), y=expression(paste(Delta,"Abstraction"))) + 
  ggtitle("Change in Abstraction vs. Negative Affect Words") +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  

# Analyses for change in LCM and change in positive affect words
data$posemo_reg <- data$posemo_RNeg - data$posemo_LNeg
data$LCM_reg <- data$LCM_RNeg - data$LCM_LNeg

lm(posemo_reg ~ LCM_reg, data = data)
corrout <- corr.test(data$posemo_reg, data$LCM_reg)
corrout$p
corrout$r

ggplot(data, aes(x=posemo_reg, y=LCM_reg)) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste(Delta,"Positive Affect Words")), y=expression(paste(Delta,"Abstraction"))) + 
  ggtitle("Change in Abstraction vs. Positive Affect Words") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  

# Reappraisal Success and Change in Abstraction
lm(ReapSuccess ~ LCM_reg, data = data)
corrout <- corr.test(data$ReapSuccess, data$LCM_reg)
corrout$p
corrout$r

ggplot(data, aes(x=LCM_reg, y=ReapSuccess)) +
  geom_point() +    # Use hollow circles
  geom_smooth(method=lm) +   # Add linear regression line 
  labs(x=expression(paste(Delta,"Abstraction")), y=expression(paste("Reappraisal Success"))) + 
  ggtitle("Reappraisal Success vs. Change in Abstraction") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=9, color = "black"), axis.text.y=element_text(size=9, color="black"))  
