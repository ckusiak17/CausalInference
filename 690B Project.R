library(mosaic)
library(tmle)
library(SuperLearner)
load("obs_data_v3.Rda")

#TMLE handcoded
X <- subset(ObsData, select = -SleepTrouble)
X1 <- X0 <- X
X1$PhysActive <- 1 #under exposure
X0$PhysActive <- 0 #under control
Covariates <- subset(ObsData, select = -c(SleepTrouble, PhysActive))
SL.library <- c('SL.mean', 'SL.glm',
                'SL.rpartPrune','SL.randomForest',
                'SL.bayesglm', 'SL.step')

QbarSL <- SuperLearner(Y = ObsData$SleepTrouble,
                       X = X,
                       SL.library = SL.library,
                       family = 'binomial')

QbarSL

QbarAW <- predict(QbarSL, newdata = X)$pred
Qbar1W <- predict(QbarSL, newdata = X1)$pred
Qbar0W <- predict(QbarSL, newdata = X0)$pred

tail(data.frame(QbarAW,Qbar1W,Qbar0W))
PsiHat.SS = mean(Qbar1W - Qbar0W)
PsiHat.SS

#Estimating g0(A|W)
gHatSL <- SuperLearner(Y = ObsData$PhysActive, 
                       X = subset(ObsData,select = -c(PhysActive,SleepTrouble)),
                       SL.library = SL.library,
                       family = 'binomial')
gHatSL
gHatAW <- rep(NA,nrow(ObsData))
gHat1W <- gHatSL$SL.predict
gHat0W <- 1 - gHat1W

summary(data.frame(gHat1W,gHat0W))

gHatAW[ObsData$PhysActive==1]<- gHat1W[ObsData$PhysActive==1]
gHatAW[ObsData$PhysActive==0]<- gHat0W[ObsData$PhysActive==0]

tail(data.frame(ObsData$PhysActive, gHatAW, gHat1W, gHat0W))

#Creating the clever covariate
H.AW <- as.numeric(ObsData$PhysActive==1)/gHat1W - as.numeric(ObsData$PhysActive==0)/gHat0W
H.1W <- 1/gHat1W
H.0W <- -1/gHat0W

tail(data.frame(ObsData$PhysActive,H.AW,H.1W,H.0W))

PsiHat.IPTW <- mean(H.AW*ObsData$SleepTrouble)
PsiHat.IPTW

#Identifying epsilon
logitUpdate<- glm(ObsData$SleepTrouble ~ -1 +offset(qlogis(QbarAW)) + H.AW, family='binomial')
epsilon <- logitUpdate$coef
epsilon


# obtain the targeted estimates
QbarAW.star<- plogis(qlogis(QbarAW)+ epsilon*H.AW)
Qbar1W.star<- plogis(qlogis(Qbar1W)+ epsilon*H.1W)
Qbar0W.star<- plogis(qlogis(Qbar0W)+ epsilon*H.0W)

QbarAW.star
PsiHat.TMLE <- mean(Qbar1W.star - Qbar0W.star)
exp(c(PsiHat.SS, PsiHat.IPTW, PsiHat.TMLE))

#return point estimates, targeted estimates of Qbar0AW
estimates <- data.frame(cbind(PsiHat.SS=PsiHat.SS, PsiHat.IPTW, PsiHat.TMLE))
Qbar.star <- data.frame(cbind(QbarAW.star, Qbar1W.star, Qbar0W.star))
names(Qbar.star)<- c('QbarAW.star', 'Qbar1W.star', 'Qbar0W.star')
head(list(estimates=estimates, Qbar.star=Qbar.star, H.AW=H.AW))

#tmle package

tmle.SL <- tmle(Y = ObsData$SleepTrouble,
                A = ObsData$PhysActive,
                W = Covariates,
                Q.SL.library =  SL.library,
                g.SL.library = SL.library)
summary(tmle.SL)


#Appendix
#Comparing the performance of the SuperLearner algorithm against that of the candidate algorithms
CV.Sl.out <- CV.SuperLearner(Y = ObsData$SleepTrouble,
                             X = X,
                             V = 10,
                             SL.library = SL.library,
                             family = 'binomial')
summary(CV.Sl.out)
