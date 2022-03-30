

getwd() # check 
#setwd("C:\\Users\\Gerhard Tutz\\LRZ Sync+Share\\TuRaschThresholdM\\RPaperPsychom\\")


############ Fit fears data set: first fixed difficulty functions + plots
#                          later Splines + plots



### run first:

source("packages.R")
source("ProgramsFixed.R")
source("ProgramsFixedAlpha.R")
source("ProgramsBSplinesPure.R")
source("ProgramsFixedandBasis.R")


##### fears  data data set I x P matrix
dat <- read.table("fear", header=TRUE,sep= " ")
###################

I <- dim(dat)[1]
P <- dim(dat)[2]



lin <- "log"  
indicator<- "D"


### fitting with varying slopes logarithmic difficulty functions, discrete

ThrV <- ThreshModFixed(dat,commonslope="var" ,indicator= "D", der="der",lin =lin,start, startind="no")
ThrV

##### fitting with varying slopes logarithmic difficulty functions, continuous

ThrVC <- ThreshModFixed(dat,commonslope="var" ,indicator= "C", der="der",lin =lin,start, startind="no")
ThrVC



### fitting with fixed slopes logarithmic difficulty functions, discrete

Thr <- ThreshModFixed(dat,commonslope="com" ,indicator= "D", der="der",lin =lin,start, startind="no")
Thr

# restart until Thr$convergence=0
start <- Thr$par
Thr <- ThreshModFixed(dat,commonslope="com" ,indicator= "D", der="no",lin =lin,start=start, startind="yes")
Thr

### fitting with fixed slopes logarithmic difficulty functions, continuous

lin <- "log"  
ThrC <- ThreshModFixed(dat,commonslope="com" ,indicator= "C", der="der",lin =lin,start, startind="no")
ThrC

# restart until Thr$convergence=0
start <- ThrC$par
ThrC <- ThreshModFixed(dat,commonslope="com" ,indicator= "C", der="no",lin =lin,start=start, startind="yes")
ThrC




#### fits with discrimination parameter, continuous, commonslope: var (varying slopes)

commonslope<-"var"  
indicator<- "D" 
der<-"no"
lin <- "log"
start<-ThrV$par
start<-c(start,rep(1,3))
ThrCVar <- ThreshModFixedAlpha(dat,commonslope="var" ,indicator= "C", der="der",lin =lin,start, startind="no"
                            ,penalpha=0)

#####
ThrDVar <- ThreshModFixedAlpha(dat,commonslope="var" ,indicator= "D", der="der",lin =lin,start, startind="no"
                               ,penalpha=0)

#refit
start <- ThrDVar$par
ThrDVar <- ThreshModFixedAlpha(dat,commonslope="var" ,indicator= "D", der="no",lin =lin,start, startind="yes"
                               ,penalpha=0)
ThrDVar



#### Plots without discrimination
parmatrest <- ThrVC$parmatrix ### or other fitted models


### for data transformation if lin= inv
width<- max(dat)-min(dat)
norm<- 10  ### 10 means in (0.0833 0.91666)
mintr <- min(dat)- width/norm
maxtr <- max(dat)+ width/norm
datstd<- (dat -mintr)/(maxtr-mintr)
#####################################

###plot   delta fct 
#parmatrest <- ThrV$parmatrix ### 
min <- min(dat)
max <- max(dat)
ylims <- c(-5,2)
x<- seq(min,max,(max-min)/40)
pcdum <-c ('1','2','3','4','5','6')
#x<- seq(0,maxresp,1)
y<- 0*x
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*x
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+x)
if(lin == "inv"){xt <- (x -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
plot(x,delta,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Difficulty functions",xlab ="y", type="b",lwd=1.2,
     ylab="",pch =pcdum[1],cex = 1.2, ylim=ylims )

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*x
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+x)
  if(lin == "inv"){xt <- (x -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  lines(x,delta,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=1.2, type="b",pch =pcdum[i],cex = 1.2)
}

####PT function 
#parmatrest <- ThrV$parmatrix ### 
pcdum <-c ('1','2','3','4','5','6')

y <- seq(min,max,(max-min)/40)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
theta <- 0 
problargy <- pnorm(theta - delta,0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="P(Y>y)",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- pnorm(theta - delta,0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}

##### IC curves
#parmatrest <- Thr$parmatrix ### 
y=0 ### specified y-value

theta <- seq(-8,4,0.1)
numtheta <- length(theta)

if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
problargy <- pnorm(theta - delta,0,1)
plot(theta,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="IC function",xlab ="theta", type="b",
     lwd=2.0,ylim=c(0,1), pch =pcdum[1])

for (l in 2:I){
  if(lin == "lin") delta <- parmatrest[l,1]+ parmatrest[l,2]*y
  if(lin == "log") delta <- parmatrest[l,1]+ parmatrest[l,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[l,1]+ parmatrest[l,2]*log(xt/(1-xt))}
  problargy <- pnorm(theta - delta,0,1)
  
  lines(theta,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=2.0, type="b",pch =pcdum[l],cex = 1)
}

#### densities
theta <- 0  # specify theta

pcdum <-c ('1','2','3','4','5','6')
ylims <- c(0,0.5)
y <- seq(min,max,(max-min)/70)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}

problargy <- dnorm(theta - delta,0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Densities",xlab ="y", type="b",lwd=1.2,ylim=ylims,
     ylab="", pch =pcdum[1])

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- dnorm(theta - delta,0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}


###### plots with alpha

parmatrest <- ThrCVar$parmatrix ### for continuous
alpha<- ThrCVar$alpha

parmatrest <- ThrDVar$parmatrix ### for discrete
alpha<- ThrDVar$alpha

pcdum <-c ('1','2','3','4','5','6')

y <- seq(min,max,(max-min)/40)
if(lin == "lin")delta <- parmatrest[1,1]+ parmatrest[1,2]*y
if(lin == "log")delta <- parmatrest[1,1]+ parmatrest[1,2]*log(1+y)
if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
delta <- parmatrest[1,1]+ parmatrest[1,2]*log(xt/(1-xt))}
theta <- 0 

problargy <- pnorm(alpha[1]*(theta - delta),0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="P(Y>y)",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  if(lin == "lin")  delta <- parmatrest[i,1]+ parmatrest[i,2]*y
  if(lin == "log")  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(1+y)
  if(lin == "inv"){xt <- (y -mintr)/(maxtr-mintr)
  delta <- parmatrest[i,1]+ parmatrest[i,2]*log(xt/(1-xt))}
  problargy <- pnorm(alpha[1]*(theta - delta),0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}


####################################################
##### Posterior person parameters  fixed difficulty

parmatrest<-ThrVC$parmatrix
stdest <- ThrVC$stdmixt
grid <- seq(-3,5,0.1)  ##theta grid to compute estimates

PostFitV <- PosteriorEstimates(grid,dat,I,indicator,lin =lin,parmatrest,stdest) 

colSums(dat)  
plot(colSums(dat),PostFitV)
cor(colSums(dat),PostFitV)




####################################################
#################  Splines Fits
######################################### 

indicatorsingle  <- "C"
indicatorvec  <- indicatorsingle
for (l in 2:I )indicatorvec<- c(indicatorvec,indicatorsingle)

numknotsstart<-4 #( 6 produces 8 Bsplines)
ord <-4

ThrExts <- Threshbasis(dat,basis,indicatorvec,ord, numknotsstart,lambdpenin,lambdc=0,start, startind,der="der")


##################plots splines


itemmatrixest <- ThrExts$parmatrix
minresp <- min(dat)
maxresp <- max(dat)
dif<- maxresp-minresp
knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))


###plot splines - delta fct

#x<- seq(minresp,maxresp,(maxresp-minresp)/50)
x<- seq(1,7,0.15)
y<- 0*x
for (l in 1:length(x)){s <-splineDesign(knots,x[l], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
y[l]<- s%*%itemmatrixest[1,]
#y[l]<- s%*%par2
}
plot(x,y,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Difficulty functions",xlab ="y", type="b",lwd=1.2,
     ylim=c(-4,3), ylab="",pch =pcdum[1])
for(it in 2:I){
  for (l in 1:length(x)){s <-splineDesign(knots,x[l], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
  y[l]<- s%*%itemmatrixest[it,]}
  lines(x,y,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, xlab ="y", type="b",lwd=1.2,
        pch =pcdum[it])}
#######

####PT function  
pcdum <-c ('1','2','3','4','5','6')

#y <- seq(minresp,maxresp,(max-min)/50)
y<- seq(1,7,0.15)
s <-splineDesign(knots,y, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)

delta <-  s%*%itemmatrixest[1,]
#plot(y,delta)


theta <- 0 
problargy <- pnorm(theta - delta,0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="P(Y>y)",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  delta <-  s%*%itemmatrixest[i,]
  problargy <- pnorm(theta - delta,0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}


##### density continuous
theta <- 0 ### specify value

pcdum <-c ('1','2','3','4','5','6')

y <- seq(minresp,maxresp,(max-min)/50)
s <-splineDesign(knots,y, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)

delta <-  s%*%itemmatrixest[1,]

problargy <- dnorm(theta - delta,0,1)
plot(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Densities",xlab ="y", type="b",lwd=1.2,ylim=c(0,0.5),
     ylab="", pch =pcdum[1])

for (i in 2:I){
  delta <-  s%*%itemmatrixest[i,]
  problargy <- dnorm(theta - delta,0,1)
  
  lines(y,problargy,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lty=2,lwd=1.2, type="b",pch =pcdum[i],cex = 1)
}

##### density discrete
theta <- 0 ### specify value

pcdum <-c ('1','2','3','4','5','6')
limit<-30
y <- seq(0,limit,1)

s <-splineDesign(knots,y, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)

delta <-  s%*%itemmatrixest[1,]

problargy <- pnorm(theta - delta,0,1)
limit1<-limit+1
prob<- matrix(0,limit1,1)
prob[1,1]<-1-problargy[1]
for (l in 2:limit1) prob[l,1]<-problargy[l-1] -problargy[l]

#plot(y,problargy ,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Densities",xlab ="y", type="b",lwd=1.2,ylim=c(0,1),
#     ylab="", pch =pcdum[1])
plot(y,prob ,cex.axis=1.5,cex.lab=1.5,cex.main=1.5, main ="Densities",xlab ="y", type="b",lwd=1.6,
     ylab="", pch =pcdum[1],ylim=c(0,0.20), cex=1.3)
sum(prob)

for (i in 2:I){
  delta <-  s%*%itemmatrixest[i,]
  problargy <- pnorm(theta - delta,0,1)
  prob[1,1]<-1-problargy[1]
  for (l in 2:limit1) prob[l,1]<-problargy[l-1] -problargy[l]
  print(sum(prob))
  lines(y,prob,cex.axis=1.5,cex.lab=1.5,cex.main=1.5,lwd=1.2, type="b",pch =pcdum[i],cex = 1.3)
}
###### comparison links


##### Posterior person parameters  splines

parmatrest<-ThrExts$par
stdest <- ThrExts$stdmixt
grid <- seq(-3,5,0.1)  ##theta grid to compute estimates
stdest<-ThrExts$stdmixt

PostFit <- 
  PosteriorEstimatesSplines(grid,dat,I,indicatorvec,lin =lin,ord, numknotsstart,parmatrest,stdest=1) 
colSums(dat)  
plot(colSums(dat),PostFit)
cor(colSums(dat),PostFit) ### correlation with scores

cor(PostFit,PostFitV)

