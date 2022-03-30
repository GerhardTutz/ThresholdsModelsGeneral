


#### Programs fitting with B-splines -item specific curves penalized
####  here indicator is vector!!

#### also derivative discrete


### here with dummy for discrete
####  auch modifizierte programme 

#library("splines2")
library(splines)
library(isoreg)
#### 


###### programs

###  ### ############### B-Splines program (B splines only no additional fixed function)

#### General program for Fits
####
#Tr <-Threshbasis(dat,basis =basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind="no",der="no")
# former version BsplinesprogramsJune21 (R2)

######################################

Threshbasis <- function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der){
  
  #### datb I x P matrix responses
  ####  basis: splines for splines, lin for linear,
  #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
  #### ord:           order of B-Splines
  #### numknotsstart: initial number of basis functions
  #### lambdpenin:    penalizes violation of increase 
  #### lambdc:        penalizes variation of coefficients 0: varying coefficients
  
  #### mit stdmixt
  datb<- dat
  
  #Memlist <- list("PenM1" = I)
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  itemsums <- rowSums(datb)/P
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  if(basis =="splines"){
    dif<- maxresp-minresp
    knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
    
    ### determine number of knots and basis
    sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
    numknots <- dim(sdim)[2]
    
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
    
    for (i in 1:I){for (p in 1:P){
      bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
    
    #fact <-1
    #parp   <- seq(1,numknots,1)
    #par <- parp-itemsums[1]
    #for (l in 2:I) par <- c(par,parp-itemsums[l])/2
    #par <- c(par,1)
    
    fact <-1
    parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
    par <- parp
    for (l in 2:I) par <- c(par,parp)
    par <- c(par,1)
    }### end splines
  #length(par)
  
  if(basis =="lin"){
    numknots <- 2
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,datb[i,p])}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,1)}
    
    for (i in 1:I){for (p in 1:P){
      datdum <- datb[i,p]-1
      basdummy[i,p,]<- c(1,datdum)}}
    
    #fact <-1
    #parp   <- fact*seq(1,numknots,1)
    #par <- parp
    #for (l in 2:I) par <- c(par,parp)  
    #par <- c(par,1) 
    par1 <- -itemsums  ### 
    par <- c(par1[1],1) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],1)
    par <- c(par,1) ### stdmixt
  } # end lin
  
  
  
  if(basis =="log"){
    numknots <- 2
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,log(1+datb[i,p]))}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,1/(1+datb[i,p]))}
    
    for (i in 1:I)if (indicator[i] !="C"){{for (p in 1:P){
      datdum <- datb[i,p]-1
      if (datdum > -1) basdummy[i,p,]<- c(1,log(datdum+1))}}} ### discrete corrected
    
    #fact <-1
    #parp   <- fact*seq(1,numknots,1)
    #par <- parp
    #for (l in 2:I) par <- c(par,parp)   
    #par <- c(par,1) 
    
    par1 <- -itemsums/3  ### 
    par <- c(par1[1],1) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],2)
    par <- c(par,1) ### stdmixt 
  } # end log
  
  if(basis =="inv"){
    numknots <- 2
    bas <- array(0, dim=c(I,P,numknots))
    basderiv <- array(0, dim=c(I,P,numknots))
    basdummy <- array(0, dim=c(I,P,numknots)) 
    
    low <- 0.1  ## in interval (low, 1-low)
    
    datdum <- low -(1-2*low)*minresp/(maxresp-minresp)+(1-2*low)*datb/(maxresp-minresp)
    
    for (i in 1:I){for (p in 1:P)   bas[i,p,]<- c(1,log(datdum[i,p]/(1-datdum[i,p])))}
    for (i in 1:I){for (p in 1:P)   basderiv[i,p,]<- c(0,1/(datdum[i,p]*(1-datdum[i,p])))}
    
    for (i in 1:I){if (indicator[i] !="C"){for (p in 1:P){
      datdum <- datb[i,p]-1
      if (idatdum > 0) basdummy[i,p,]<- c(1,1/(datdum*(1-datdum)))}}} ### discrete not clear
    
    #fact <-1
    #parp   <- fact*seq(1,numknots,1)
    #par <- parp
    #for (l in 2:I) par <- c(par,parp)  
    #par <- c(par,1) 
    
    par1 <- -itemsums  ### 
    par <- c(par1[1],1) ### 1:slope difficulty
    for(i in 2:I) par <- c(par,par1[i],1)
    par <- c(par,1) ### stdmixt 
  } # end inv
  
  
  
  
  
  ##test
  #it <-3
  #val <--(par[it]+par[it+3]*datb[it,])
  #min(val)
  #max(val)
  #### 
  
  if(startind == "yes") par <- start 
  
  if(der!="der") {fittrw <- optim(par, LoglikIntExtBsplIN, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,
                  I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                  lower = -Inf, upper = Inf,control = list(), hessian = FALSE) }
  
  if(der=="der") {fittrw <- optim(par, LoglikIntExtBsplIN, gr = derLoglikBsplI,lambdpenin ,lambdc ,numknots=numknots,
                                  dat=datb,I=I, indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "BFGS",
                                 lower = -Inf, upper = Inf,control = list(), hessian = FALSE) }
  
  #LoglikIntExtBsplIN(par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
  #                     indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  ll <- length(par)
  itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,I)) ##I x numknots matrix
  stdmixt <- fittrw$par[ll]
  loglikorig<--fittrw$value
  
  derivfin <-0  ## dummy
  #derivative here
  
  derivfin <-derLoglikBsplI(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
                       indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  
  ## logorig
  loglikorig <- -LoglikIntExtBsplIN(fittrw$par,lambdpenin=0,lambdc=0,numknots=numknots,dat=datb,I=I,
                                    indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  ### positive stddev
  stdmixt<-abs(stdmixt) 
  fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
  ####
  AIC <-  -2*loglikorig +2*length(par)
  BIC <-  -2*loglikorig +2*log(length(par))
  
  #########
  newList <- list("parmatrix" = itemmatrixest,"stdmixt"=stdmixt,"par"=fittrw$par, "numbasis"=numknots, 
                  "Loglikwithpen"=-fittrw$value, "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                  "AIC"=AIC, "BIC"=BIC,"derivative"=derivfin)
  return(newList)
  
}
########################################################







##### Program for splines only
 
ThreshModSplines <- function(datb,commonslope,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind){
  
  #### datb I x P matrix responses
  ####  commonslope nicht benutzt
  #### indicator:     C continuous, discrete D1: 1,2,..   D0: 0,1,2 
  #### ord:           order of B-Splines
  #### numknotsstart: initial number of basis functions
  #### lambdpenin:    penalizes violation of increase 
  #### lambdc:        penalizes variation of coefficients 0: varying coefficients
  
  #### mit stdmixt
  
  
  I <- dim(datb)[1]  
  P <- dim(datb)[2]
  
  minresp <- min(datb)
  maxresp <- max(datb)
  dif<- maxresp-minresp
  knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
  
  ### determine number of knots and basis
  sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  bas <- array(0, dim=c(I,P,numknots))
  basderiv <- array(0, dim=c(I,P,numknots))
  basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
  
  for (i in 1:I){for (p in 1:P){
    bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  for (i in 1:I){for (p in 1:P){
    basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
  
  for (i in 1:I){for (p in 1:P){
    datdum <- datb[i,p]-1
    basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  fact <-1
  parp   <- fact*seq(1,numknots,1)-fact*(numknots/2)
  par <- parp
  for (l in 2:I) par <- c(par,parp)
  
  
  
  
  ###
  par <- c(par,1)
  if(startind == "yes") par <- start 
  
  fittrw <- optim(par, LoglikIntExtBsplIN, gr = NULL,lambdpenin ,lambdc ,numknots=numknots,dat=datb,I=I, 
                  indicator, bas=bas,basderiv=basderiv,basdummy=basdummy,method = "Nelder-Mead",
                  lower = -Inf, upper = Inf,control = list(), hessian = FALSE) 
  
  ll <- length(par)
  itemmatrixest<-  t(matrix(fittrw$par[1:ll-1],numknots,I)) ##I x numknots matrix
  stdmixt <- fittrw$par[ll]
  loglikorig <- -fittrw$value
  
  #dim(bas)
  #dim(basdummy)
  ### new l#oglik
  if (lambdpenin >0){
    upp<- numknots-1
    peninc<-0
    for (i in 1:I){
      for (l in 1:upp ){dif  <-itemmatrixest[i,l+1]- itemmatrixest[i,l]
      peninc <-peninc + dif *indicatt(-dif , 0)}}
    loglikorig <- -fittrw$value +  lambdpenin*peninc } 
  #### end new lik
  
  ### positive stddev
  stdmixt<-abs(stdmixt) 
  fittrw$par[length(fittrw$par)] <-abs(fittrw$par[length(fittrw$par)])  
  ####
  AIC <-  -2*loglikorig +2*length(par)
  BIC <-  -2*loglikorig +2*log(length(par))
  
  #########
  newList <- list("parmatrix" = itemmatrixest,"stdmixt"=stdmixt,"par"=fittrw$par, "numbasis"=numknots, 
                  "Loglikwithpen"=-fittrw$value, "convergence"=fittrw$convergence,"Loglik"=loglikorig,
                  "AIC"=AIC, "BIC"=BIC)
  return(newList)
  
}
########################################################


 


#### loglik mit stdmixt

LoglikIntExtBsplIN <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                              indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy){
  ##### with stdmixt
  ### item specific
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  
  ##trial
  #now <-Memlist$PenM1 
  
  P <- dim(dat)[2]
  
  numpar <- length(par) 
  sdnormal<- par[numpar]
  #paroriginal <- par[1:numpar-1]
  
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = sdnormal,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots)
    ind[p]<-log(int)
    #int<-abs(int)  # 
    
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
  
    } #  end p
    
    dimBspl<-length(par)-1
    parBspl <- par[1:dimBspl]
    itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
    
    
    ####penalty increase
    #lambdpenin <- 0
    if (lambdpenin >0){
      
      upp<- numknots-1
      peninc<-0
      for (i in 1:I){
        for (l in 1:upp ){dif  <-itemmatrix[i,l+1]- itemmatrix[i,l]
        peninc <-peninc + dif^2 *indicatt(-dif , 0)}}
      
      #pen <- dif%*%dif
      sum <- sum -  lambdpenin*peninc ### + weil negativ , jetzt - ????
      
    } #end if
    ############
    
    #### penalty common difficulties
    #lambdc <- 3
    if (lambdc >0){
      upp1<- I-1
      upp2<- numknots-1
      penc<-0
      for (i in 1:upp1){
        for (l in 1:upp2 ){dif  <-(itemmatrix[i+1,l+1]- itemmatrix[i+1,l])-(itemmatrix[i,l+1]- itemmatrix[i,l])
        penc <- penc +dif*dif}}
      
      ##### funktioniert: with matrix penc = pencc jetzt auskommentiert
      
      #upp1<- I-1
      #upp2<- numknots-1
      #counth <- 0
      #for (i in 1:upp1){
      #  for (l in 1:upp2 ){
      #  counth <- counth +1
      #  penM <- matrix(0,1,length(par))
      #  penM[(i-1)*numknots+l+1]<-1
      #  penM[(i-1)*numknots+l]<--1
      #  penM[(i)*numknots+l+1]<--1
      #  penM[(i)*numknots+l]<- 1
      #  if (counth <= 1) pennew <- penM
      # if(counth > 1) pennew <- rbind(pennew,penM)
      # }}
      #pencc <- t(pennew%*%par)%*% (pennew%*%par) 
      
      ####ende auskommentiert
      
      if (numknots == 2) {penc<-0
        for (i in 2:I) penc<- penc + (itemmatrix[i,2]-itemmatrix[i-1,2])^2
        
        ## alternative pennew=penc
        M<- matrix(0,I-1,I)
        d1 <- c(-1,1)
        M[1,]<- c(d1,rep(0,I-2))
        M[I-1,]<-c(rep(0,I-2),d1)
        I2 <- I-2
        for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
        pennew <- t(itemmatrix[,2])%*%t(M)%*%M%*%itemmatrix[,2]
        ###
        }
      
      
      sum <- sum -  lambdc*penc ### 
      #sum <- sum -  lambdc*pennew 
    } #end if
    
  
  sum<- -sum
  return(sum)}
###############




##################################
indicatt <-function(x,thresh){
  ret <-0
  if (x > thresh) ret <-1 
  return(ret)}

###########################
prodfctIntExtBsplINVec  <-function(theta,datitem=datitem,I=I,par=par,slope=slope,
                               indicator=indicator, bas,basderiv,basdummy,pers=pers,numknots=numknots){
  #### with std mixt only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  prod <-1    
 
    dimBspl<-length(par)-1
    parBspl <- par[1:dimBspl]
    itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
    
    for (i in 1:I){
      iteml <- itemmatrix[i,]
      item=i
      if(indicator[i]=="C"){
      prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                       item=item,pers=pers)} 
      if(indicator[i]=="D0"){
        prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                              item=item,pers=pers,indicator[i])}
      if(indicator[i]=="D1"){
        prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                              item=item,pers=pers,indicator[i])} ##idetisch mit D0
      }

  return(prod)}  
#########################



########################### here linear or log
deltalinNVdensBsplI <-function(iteml,slope , theta, bas,basderiv,basdummy=basdummy,item=item,pers=pers){
  #### computes density function continuous
  ###  iteml parameter for current item
  ### p is person
  #### 
  #### first item fixed
  
  #DFval <- matrix(0,length(y),1)
  #for (i in 1:length(y)){
  #####
  val <-  bas[item,pers,]%*%iteml
  valder <- basderiv[item,pers,]%*%iteml
  
  DFval<-   dnorm(slope*theta - val)*valder
  #}
  return(DFval)}
#########################


########################### here linear or log
deltalinNVdensDiscrBsplI <-function(obs,iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                    item=item,pers=pers,indicatorsingle){
  #### computes discrete density for count data
  ###  for observ 0,1,2... or 1,2,  
  ### ###  iteml parameter for current item
  
  
  if(indicatorsingle=="D0"){
    
    if(obs <= 0){val <-  bas[item,pers,]%*%iteml
    DFval<-   1- pnorm(slope*theta - val)}
    
    if(obs > 0){
      valprev <-  basdummy[item,pers,]%*%iteml
      val <-  bas[item,pers,]%*%iteml
      
      DFval<-   pnorm(slope*theta - valprev) - pnorm(slope*theta - val)}
  }#end indicator if
  
  if(indicatorsingle=="D1"){
    
    if(obs <= 1){val <-  bas[item,pers,]%*%iteml
    DFval<-   1- pnorm(slope*theta - val)}
    
    if(obs > 1){
      valprev <-  basdummy[item,pers,]%*%iteml
      val <-  bas[item,pers,]%*%iteml
      
      DFval<-   pnorm(slope*theta - valprev) - pnorm(slope*theta - val)}
  }#end indicator if
  
  
  return(DFval)}
#########################


################## derivative

derLoglikBsplI <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                          indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy){
  ### derivative, only continuous
  #parmatr <- matrix(par,I,2)
  # derivfin <-derLoglikBsplI(fittrw$par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=datb,I=I,
  #indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  
  slope <- 1 ### slope fixed  weight on theta!
  P <- dim(dat)[2]
  
  lpar <- length(par)
  der <- matrix(0,lpar,1)
  
  parmatr <- t(matrix(par[1:lpar-1],numknots,I))
  stdmixt <- par[lpar]
  
  # basis functions prodfctIntExtBsplINVec 
  
  for (p in 1:P){
    datitem <- dat[,p]
    
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                        numknots=numknots) 
    for (l in 1:I){
      for (j in 1:numknots){ 
        comp <- l
        pers<-p
        count <- (l-1)*numknots+j
        #der[l,1] <-der[l,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
        #                    indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp)/int } 
        der[count,1] <-der[count,1]+gauss.hermite(prodfctIntExtBsplINderiv1, mu = 0, sd = stdmixt,datitem=datitem,
                          numknots=numknots,I=I,par=par, slope=slope,indicator=indicator,order = 10, bas=bas,
                          basderiv=basderiv,basdummy=basdummy,pers=pers,comp=comp, j=j)/int } # end j 
    } #end l
  } # end p
  
  
  
  # stdmixt
  omit <-1  
  if(omit <= 0){
  for (p in 1:P){
    datitem <- dat[,p]
    pers<-p
    intm <-gauss.hermite(prodfctIntExtBsplINVec, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,slope=slope,
                         indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                         numknots=numknots) 
    
    smixt <- par[lpar]
    der[lpar,1] <-der[lpar,1]+
      gauss.hermite(prodfctIntExtBsplINderivmixt, mu = 0, sd = stdmixt,datitem=datitem,I=I,par=par,
                    slope=slope,indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers,
                    smixt=smixt)/intm
    
  } # end p
} #end omit
  #der[lpar,1]<- -der[lpar,1]
  
  
  #### now numeric !!!!!
  deltan<- 0.00005
  vect1<- par
  vect1[lpar] <- vect1[lpar]+deltan
  #loglikorig  <- LoglikN(par,dat,I,indicator, lin=lin)
  #loglikorignow <- LoglikN(vect1,dat,I,indicator, lin=lin)
  
  loglikorig<-LoglikIntExtBsplIN(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  loglikorignow<-LoglikIntExtBsplIN(vect1,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                                 indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy)
  derivnow <- (loglikorignow-loglikorig)/deltan
  der[lpar,1]<- -derivnow
  
  
  ##### end numeric
  
  
  
  #lambdc <- 3
  
  if (lambdpenin >0){
    numknots
    M<- matrix(0,numknots-1,numknots)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,numknots-2))
    numknots1 <- numknots-1
    if(numknots >2)for(l in 2:numknots1) M[l,]<- c(rep(0,l-1),d1,rep(0,numknots-l-1))
    dum <- diag(1,I)
    Mtot <- kronecker(dum,M)
    Mtot<- cbind(Mtot,rep(0,dim(Mtot)[1]))
    pen1 <- 2*t(Mtot)%*%(Mtot%*%par*indicattmatrix(-Mtot%*%par,0))
    
    der <- der - lambdpenin*pen1
    }
    ## end penin
   #dim(Mtot)
   
  if (lambdc >0){
    
    if (numknots > 2){
    upp1<- I-1
    upp2<- numknots-1
    counth <- 0
    for (i in 1:upp1){
      for (l in 1:upp2 ){
        counth <- counth +1
        penM <- matrix(0,1,length(par))
        penM[(i-1)*numknots+l+1]<-1
        penM[(i-1)*numknots+l]<--1
        penM[(i)*numknots+l+1]<--1
        penM[(i)*numknots+l]<- 1
        if (counth <= 1) pennew <- penM
        if(counth > 1) pennew <- rbind(pennew,penM)
      }}
    pencc <- t(pennew%*%par)%*% (pennew%*%par) 
    
    deriv1 <- 2*t(pennew)%*%pennew%*%par
    der <- der -  lambdc*deriv1
    }
      
    if (numknots == 2) {
    M<- matrix(0,I-1,I)
    d1 <- c(-1,1)
    M[1,]<- c(d1,rep(0,I-2))
    M[I-1,]<-c(rep(0,I-2),d1)
    I2 <- I-2
    for(l in 2:I2) M[l,]<- c(rep(0,l-1),d1,rep(0,I-l-1))
    
    deriv <- 2*t(M)%*%M%*%parmatr[,2]
    for (l in  1:I){der[l*2,1]<-der[l*2,1] -lambdc*deriv[l]}
   
  }} #end if
  
  
  der <- -der
  return(der)}
###############
indicattmatrix <-function(x,thresh){
  dim <- length(x)
  ret <-matrix(0,length(x),1)
  for (l in 1:dim)if (x[l] > thresh) ret[l] <-1 
  return(ret)}
 

###########################
prodfctIntExtBsplINderiv1 <-function(theta,datitem=datitem,numknots=numknots,I=I,par=par,slope=slope,
                                     indicator=indicator,order = order, bas,basderiv,basdummy,
                                     pers=pers,comp=comp, j=j){
  #### with std mixt only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  #####  j basis function
  
  prod <-1    
  set <- seq(1,I,1)
  
  
    dimBspl<-length(par)-1
    parBspl <- par[1:dimBspl]
    itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
    
    
    nonset <- set[-comp] 
    
    for (i in nonset){
      iteml <- itemmatrix[i,]
      item=i
      if(indicator[i]=="C"){prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,
                                                         basdummy=basdummy,item=item,pers=pers)} 
      if(indicator[i]=="D0"){prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,
                                                          basdummy=basdummy,item=item,pers=pers,indicator[i])}
      if(indicator[i]=="D1"){prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,
                                                           basdummy=basdummy,item=item,pers=pers,indicator[i])}
      }
    
    ## der f
    iteml <-itemmatrix[comp,]
    deltaf <-  bas[comp,pers,]%*%iteml
    prodderiv <- basderiv[comp,pers,]%*%iteml
    
    
    if(indicator[comp]=="C"){prod <- prod*(basderiv[comp,pers,j]*dnorm(theta-deltaf)+
                    (theta-deltaf)*dnorm(theta-deltaf)*bas[comp,pers,j]*prodderiv)}
    
    if(indicator[comp]=="D1"){ 
      if(datitem[comp] <= 1){
      val <-  bas[comp,pers,]%*%iteml
      prod <- prod*dnorm(slope*theta - val)*bas[comp,pers,j]}
    if(datitem[comp] > 1){
        val <-  bas[comp,pers,]%*%iteml
        valprev <-  basdummy[comp,pers,]%*%iteml
        prod <- prod*(-dnorm(slope*theta - valprev)*basdummy[comp,pers,j]+
                        dnorm(slope*theta - val)*bas[comp,pers,j])}
      }## end D1 
 
    if(indicator[comp]=="D0"){ 
      if(datitem[comp] <= 0){
        val <-  bas[comp,pers,]%*%iteml
        prod <- prod*dnorm(slope*theta - val)*bas[comp,pers,j]}
      if(datitem[comp] > 0){
        val <-  bas[comp,pers,]%*%iteml
        valprev <-  basdummy[comp,pers,]%*%iteml
        prod <- prod*(-dnorm(slope*theta - valprev)*basdummy[comp,pers,j]+
                        dnorm(slope*theta - val)*bas[comp,pers,j])}
    }## end D1 
    
   
  return(prod)}

#########################

prodfctIntExtBsplINderivmixt <-function(theta,datitem=datitem,numknots=numknots,I=I,par=par,slope=slope,
                                   indicator=indicator,order = order, bas,basderiv,basdummy,pers=pers,smixt=smixt){
  #### with std mixt only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  prod <-1    
  
    #low<- I+1
    #up <- I+numknots
    dimBspl<-length(par)-1
    parBspl <- par[1:dimBspl]
    itemmatrix<-  t(matrix(parBspl,numknots,I)) ##I x numknots matrix
    
   
    for (i in 1:I){
      iteml <- itemmatrix[i,]
      item=i
      if(indicator[i]=="C")prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, 
                                          bas,basderiv,basdummy=basdummy,item=item,pers=pers) 
      if(indicator[i]=="D0"){
        prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                              item=item,pers=pers,indicator[i])}
      if(indicator[i]=="D1"){
        prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                              item=item,pers=pers,indicator[i])} ##idetisch mit D0
      }
    
    prod <- prod*(theta^2-smixt^2)/smixt^3
    
   
  return(prod)}

#########################



##### Posterior Splines

#Threshbasis <- function(dat,basis,indicator,ord, numknotsstart,lambdpenin,lambdc,start, startind,der){

#PosteriorEstimates <-function(grid,dat,I,indicator,lin =lin,parmatrest,stdest)

#prodfctIntExtBsplINVec  <-function(theta,datitem=datitem,I=I,par=par,slope=slope,
#   indicator=indicator, bas,basderiv,basdummy,pers=pers,numknots=numknots)

PosteriorEstimatesSplines <-function(grid,dat,I,indicatorvec,lin =lin,ord, numknotsstart,parmatrest,stdest){ 
  
  ####computes centered posterior person parameters
  ######   parmatrest: estimated parameters
  
  
  numgrid <- length(grid)
  dens <- matrix(0,numgrid,1)
  P<-dim(dat)[2]
  esttheta <- matrix(0,P,1)
  
  ### P-splines
  datb<- dat
  I <- dim(datb)[1]  
  
  minresp <- min(datb)
  maxresp <- max(datb)
  
  
  dif<- maxresp-minresp
  knots <- seq(minresp-(ord-1)*dif/(numknotsstart-1),maxresp+(ord-1)*dif/(numknotsstart-1),dif/(numknotsstart-1))
  
  ### determine number of knots and basis
  sdim <-splineDesign(knots,datb[1,1], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)
  numknots <- dim(sdim)[2]
  
  bas <- array(0, dim=c(I,P,numknots))
  basderiv <- array(0, dim=c(I,P,numknots))
  basdummy <- array(0, dim=c(I,P,numknots))  #for discrete only
  
  for (i in 1:I){for (p in 1:P){
    bas[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  for (i in 1:I){for (p in 1:P){
    basderiv[i,p,]<- splineDesign(knots,datb[i,p], ord = ord, derivs=1, outer.ok = TRUE, sparse = FALSE)}}
  
  for (i in 1:I){for (p in 1:P){
    datdum <- datb[i,p]-1
    basdummy[i,p,]<- splineDesign(knots,datdum, ord = ord, derivs=0, outer.ok = TRUE, sparse = FALSE)}}
  
  # end splines
  
  ###
  
  for (p  in 1:P){obs<- dat[,p]
  for (l in 1:numgrid){
    theta <- grid[l]
    #dens[l,1]<-prodfct(theta,dat=obs,I=I,parmatr=parmatrest,slope=1,indicator=indicator, lin =lin)*dnorm(theta,0,stdest)
    dens[l,1]<-prodfctIntExtBsplINVec(theta,dat=obs,I=I,par=parmatrest,slope=1,
                                      indicator=indicatorvec,bas,basderiv,basdummy,pers=p,numknots=numknots)*dnorm(theta,0,stdest)
  }
  num<-  which.max(dens)
  esttheta[p]<- grid[num]
  }
  #### centering
  esttheta<-esttheta-mean(esttheta)
  
  return(esttheta)
}
###### end fct  


















###### older without stdmixt
#### loglik ohne stdmixt

LoglikIntExtBsplI <-function(par,lambdpenin=lambdpenin,lambdc=lambdc,numknots=numknots,dat=dat,I=I,
                             indicator=indicator, bas=bas,basderiv=basderiv,basdummy=basdummy){
  ### item specific
  #parmatr <- matrix(par,I,2)
  # par vector intercept items, length I+1 letzter is slope
  slope <- 1 ### slope fixed  weight on theta!
  #I <- dim(dat)[1]
  # lambdpenin penalizes decrease
  # lambdc penalizes variation across items
  P <- dim(dat)[2]
  sum <- 0
  ind <- rep(0,P)
  for (p in 1:P){
    datitem <- dat[,p]
    #slopeitem <- par[I+1]
    pers<-p
    int <-gauss.hermite(prodfctIntExtBsplI, mu = 0, sd = 1,datitem=datitem,numknots=numknots,I=I,par=par,slope=slope,
                        indicator=indicator,order = 10, bas=bas,basderiv=basderiv,basdummy=basdummy,pers=pers)
    ind[p]<-log(int)
    ####alternative integral
    #thetval<- seq(-5,5,0.05)
    #intalt<- 0
    #for (l in 1:length(thetval))intalt<- intalt+prodfctInt(thetval[l],datitem,I,par,slope)*dnorm(thetval[l])*0.05
    #int <- intalt
    #####  
    
    sum <- sum +log(int)
    
    itemmatrix <-  t(matrix(par,numknots,I)) ##I x numknots matrix
    
    ####penalty increase
    #lambdpenin <- 0
    if (lambdpenin >0){
      
      upp<- numknots-1
      peninc<-0
      for (i in 1:I){
        for (l in 1:upp ){dif  <-itemmatrix[i,l+1]- itemmatrix[i,l]
        peninc <-peninc + dif *indicatt(-dif , 0)}}
      
      #pen <- dif%*%dif
      sum <- sum +  lambdpenin*peninc ### + weil negativ
      
    } #end if
    ############
    
    #### penalty common difficulties
    #lambdc <- 3
    if (lambdc >0){
      upp1<- I-1
      upp2<- numknots-1
      penc<-0
      for (i in 1:upp1){
        for (l in 1:upp2 ){dif  <-(itemmatrix[i+1,l+1]- itemmatrix[i+1,l])-(itemmatrix[i,l+1]- itemmatrix[i,l])
        penc <- penc +dif*dif}}
      
      #pen <- dif%*%dif
      sum <- sum -  lambdc*penc ### 
      
    } #end if
    
  }
  sum<- -sum
  return(sum)}
###############

###########################
prodfctIntExtBsplI <-function(theta,datitem=datitem,numknots=numknots,I=I,par=par,slope=slope,
                              indicator=indicator, bas,basderiv,basdummy,pers=pers){
  #### only indicator C modified
  ### datitem vector length I for fixed person
  #####slope on y is 1
  
  ### Fehler eingebaut
  c <- dor
  ###
  
  prod <-1    
  if(indicator=="C"){
    #low<- I+1
    #up <- I+numknots
    #parBspl <- par[low:up]
    itemmatrix<-  t(matrix(par,numknots,I)) ##I x numknots matrix
    indc <- rep(0,I)
    for (i in 1:I){
      iteml <- itemmatrix[i,]
      item=i
      #prod <- prod*deltalinNVdensBsplI(datitem[i],iteml,slope , theta, bas,basderiv,item=item,pers=pers)}
      prod <- prod*deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,basdummy=basdummy,item=item,pers=pers) 
      indc[i]<-deltalinNVdensBsplI(iteml,slope , theta, bas,basderiv,basdummy=basdummy,item=item,pers=pers)}
  } #end if
  
  if(indicator=="D0"){
    itemmatrix<-  t(matrix(par,numknots,I)) ##I x numknots matrix
    indc <- rep(0,I)
    for (i in 1:I){
      iteml <- itemmatrix[i,]
      item=i
      #prod <- prod*deltalinNVdensBsplI(datitem[i],iteml,slope , theta, bas,basderiv,item=item,pers=pers)}
      prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                            item=item,pers=pers,indicator) 
      indc[i]<-deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                        item=item,pers=pers,indicator)}
  } #end if
  ####identisch mit D0
  if(indicator=="D1"){
    itemmatrix<-  t(matrix(par,numknots,I)) ##I x numknots matrix
    indc <- rep(0,I)
    for (i in 1:I){
      iteml <- itemmatrix[i,]
      item=i
      #prod <- prod*deltalinNVdensBsplI(datitem[i],iteml,slope , theta, bas,basderiv,item=item,pers=pers)}
      prod <- prod*deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                            item=item,pers=pers,indicator) 
      indc[i]<-deltalinNVdensDiscrBsplI(datitem[i],iteml,slope , theta, bas,basderiv,basdummy=basdummy,
                                        item=item,pers=pers,indicator)}
  } #end if
  
  if(indicator=="OG"){
    #for (i in 1:I){
    #  iteml <- c(par[i],1)
    #  prod <- prod*deltalinNVdensR(datitem[i],iteml,slope , theta)}
    
    for (i in 1:I){
      # Rasch
      #   DFval[i,1]<-   exp(y[i]*(slope*theta - val))/(1-exp(slope*theta - val))}
      # Ogive 
      prod<-   prod*((1-datitem[i])*pnorm(slope*theta - par[i])+datitem[i]*(1-pnorm(slope*theta -par[i])))}
    
  } #end if
  
  if(indicator=="R"){
    for (i in 1:I){
      exp1 <- exp(slope*theta - par[i])/(1+exp(slope*theta - par[i]))
      prod<-   prod*((1-datitem[i])*exp1+datitem[i]*(1-exp1))}
    
  } #end if
  return(prod)}  
#########################



