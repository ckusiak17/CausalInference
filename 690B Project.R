setwd('D:/Academics/690B Fall 2017')
install.packages('mosaic')
library(mosaic)
library(ltmle)
library(SuperLearner)
load("obs_data.Rda")

#TMLE
SL.library <- c('SL.mean', 'SL.glm',
                'SL.bayesglm','SL.rpartPrune',
                "SL.polymars", 'SL.randomForest')
X <- subset(obs.data, select = -SleepTrouble)
names(X)
X1 <- X0 <- X
X1$PhysActive <- 1 #under exposure
X0$PhysActive <- 0 #under control

QbarSL <- SuperLearner(Y = obs.data$SleepTrouble,
                       X = X,
                       SL.library = SL.library,
                       family = 'binomial')
#Error message - 'SuperLearner' does not work with missing values
sum(is.na(obs.data))

#3534 missing values in the dataset

missing.obs <- c(NA, ncol(obs.data))
for(i in 1:ncol(obs.data)){
  missing.obs[i] <- sum(is.na(obs.data[,i]))
}

#missing values by covariate
cbind(names(obs.data),missing.obs)

#Option 1: use only complete cases
#n = 982
sum(complete.cases(obs.data))
obs.data.mod <- na.omit(obs.data)
dim(obs.data.mod)
obs.data.mod$SleepTrouble <- as.numeric(obs.data.mod$SleepTrouble)

obs.data.mod$SleepTrouble <- obs.data.mod$SleepTrouble - 1
obs.data.mod$PhysActive <- as.numeric(obs.data.mod$PhysActive)
obs.data.mod$PhysActive <- obs.data.mod$PhysActive - 1

str(obs.data.mod)

X <- subset(obs.data.mod, select = -SleepTrouble)
X1.mod <- X0.mod <- X
X1.mod$PhysActive <- 1 #under exposure
X0.mod$PhysActive <- 0 #under control

SL.library <- c('SL.mean', 'SL.glm',
                'SL.rpartPrune','SL.randomForest',
                'SL.bayesglm', 'SL.step')

QbarSL <- SuperLearner(Y = obs.data.mod$SleepTrouble,
                       X = X,
                       SL.library = SL.library,
                       family = 'binomial')

QbarSL

QbarAW <- predict(QbarSL, newdata = X)$pred
Qbar1W <- predict(QbarSL, newdata = X1.mod)$pred
Qbar0W <- predict(QbarSL, newdata = X0.mod)$pred

tail(data.frame(QbarAW,Qbar1W,Qbar0W))
mean(Qbar1W - Qbar0W)

#Estimating g0(A|W)
gHatSL <- SuperLearner(Y = obs.data.mod$PhysActive, 
                       X = subset(obs.data.mod,select = -c(PhysActive,SleepTrouble)),
                       SL.library = SL.library,
                       family = 'binomial')
gHatSL
gHatAW <- rep(NA,nrow(obs.data.mod))
gHat1W <- gHatSL$SL.predict
gHat0W <- 1 - gHat1W

summary(data.frame(gHat1W,gHat0W))

gHatAW[obs.data.mod$PhysActive==1]<- gHat1W[obs.data.mod$PhysActive==1]
gHatAW[obs.data.mod$PhysActive==0]<- gHat0W[obs.data.mod$PhysActive==0]

tail(data.frame(obs.data.mod$PhysActive, gHatAW, gHat1W, gHat0W))

#Creating the clever covariate
H.AW <- as.numeric(obs.data.mod$PhysActive==1)/gHat1W - as.numeric(obs.data.mod$PhysActive==0)/gHat0W
H.1W <- 1/gHat1W
H.0W <- -1/gHat0W

tail(data.frame(obs.data.mod$PhysActive,H.AW,H.1W,H.0W))

PsiHat.IPTW <- mean(H.AW*obs.data.mod$SleepTrouble)
PsiHat.IPTW

logitUpdate<- glm(obs.data.mod$SleepTrouble ~ -1 +offset(qlogis(QbarAW)) + H.AW, family='binomial')
epsilon <- logitUpdate$coef
epsilon


# obtain the targeted estimates
QbarAW.star<- plogis(qlogis(QbarAW)+ epsilon*H.AW)
Qbar1W.star<- plogis(qlogis(Qbar1W)+ epsilon*H.1W)
Qbar0W.star<- plogis(qlogis(Qbar0W)+ epsilon*H.0W)

QbarAW.star
PsiHat.TMLE <- mean(Qbar1W.star - Qbar0W.star)

CV.Sl.out <- CV.SuperLearner(Y = obs.data.mod$SleepTrouble,
                             X = X,
                             SL.library = SL.library,
                             family = 'binomial')
summary(CV.Sl.out)

#using the ltmle package
obs.data.mod2 <- obs.data.mod[,c(1,3:20,2)]

obs.data.mod$HomeOwn <- as.numeric(obs.data.mod$HomeOwn)
obs.data.mod$HealthGen <- as.numeric(obs.data.mod$HealthGen)
obs.data.mod$LittleInterest <- as.numeric(obs.data.mod$LittleInterest)
obs.data.mod$Depressed <- as.numeric(obs.data.mod$Depressed)
obs.data.mod$TVHrsDay <- as.numeric(obs.data.mod$TVHrsDay)
obs.data.mod$CompHrsDay <- as.numeric(obs.data.mod$CompHrsDay)
obs.data.mod$SmokeNow <- as.numeric(obs.data.mod$SmokeNow)
obs.data.mod$RegularMarij <- as.numeric(obs.data.mod$RegularMarij)
obs.data.mod$HardDrugs <- as.numeric(obs.data.mod$HardDrugs)

ltmle.SL<- ltmle(data=obs.data.mod2, Anodes='PhysActive',
                 Ynodes='SleepTrouble',
                 abar=list(1,0),
                 SL.library=SL.library)

summary(ltmle.SL)

