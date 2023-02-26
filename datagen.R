#UPLIFT MODEL USING HECK'S BETAS AND A REVERSE ENGINEERED INTERCEPT AND EQUATION 
#HECK'S BENZENE OUTDOOR DISTRIBUTION RESULTS IN A ZERO AVERAGE TREATMENT EFFECT
#BUT THE INDOOR BENZENE DISTRIBUTION RESULTS IN A POSITIVE AVERAGE TREATMENT EFFECT
#THE PROGRAM IS NOW USING THE INDOOR BENZENE DISTRIBUTION AND SEED (2)


#https://search.r-project.org/CRAN/refmans/osDesign/html/beta0.html
#https://data.library.virginia.edu/simulating-a-logistic-regression-model/
#https://stats.stackexchange.com/questions/439177/simulating-population-level-logistic-regression-model-with-pre-specified-prevale
#https://empirical-methods.com/logistic-regression-and-friends.html
#https://yury-zablotski.netlify.app/post/simple-logistic-regression/

setwd("C:/Users/User/Documents/R/Heck")   #include directory
getwd()

#install.packages("osDesign")
#install.packages("Hmisc")
#install.packages("rms")
library(osDesign)
library(Hmisc)
library(rms)





for (x in 1:200) {


  distr = 'outdoor' #Heck's distribution
  #distr = 'indoor'
  #distr = 'stove'  
  
  updated <- FALSE
  
  if (updated == FALSE){
    prevalence <- 0.00238899 #observed prevalence from Heck.xlsx N Sheet
    #prevalence <- 0.000009*70*0.60 #prevalence according to new Heck information 
    #0.000009 incidence of pediatric AML, 70 life expectancy, 0.60 overall survival rate
    HeckDF <- read.csv("HeckDF.csv")
    
  }else{
    HeckDF <- read.csv("HeckDF_updatedN.csv")
  }
  
  
  populationweighted <- TRUE
  
  #BEGIN HANEUSE MODIFIED CODE
  ########################################################################################
  #BEGIN HANEUSE MODIFIED CODE
  
  beta00 <-
    function(betaX,
             X,
             N,
             rhoY,
             expandX) #ADDED expandX
    {
      ##
      X <- expandCatX(X, expandX) #ADDED LINE
      #print("checking the expanded matrix")
      #print(X)
      #beta0eval is the name of the loss function
      value <- optimize(beta0eval,
                        interval=logit(rhoY) + c(-2,2),
                        betaX=betaX,
                        X=X,
                        N=N,
                        rhoY=rhoY,
                        tol=.Machine$double.eps^0.5)$minimum
      ##
      return(value)
    }
  
  logit <-
    function(p) log(p/(1-p))
  
  expit <-
    function(x) exp(x)/(1 + exp(x))
  
  #the first argument of this function, beta0, is the changing cell
  beta0eval <-
    function(beta0,
             betaX,
             X,
             N,
             rhoY)
    {
      ##
      #COMMENTED OUT THIS
      # designX <- X[,1]
      # for(i in 2:ncol(X))
      # {
      #   for(j in 1:max(X[,i]))
      #   {
      #     designX <- cbind(designX, as.numeric(X[,i] == j))
      #   }
      # }
      designX = as.matrix(X) #ADDED LINE
      ##
      #etaY  <- as.numeric(designX %*% c(beta0, betaX)) #ORIGINAL
      #etaY is the logit(p)=log(p/(1-p))
      etaY  <- as.numeric(designX %*% unlist(c(beta0, betaX))) #MODIFIED
      
      ##
      #to calculate the expected prevalence get the probabilities (expit(etaY)) 
      #--note that expit(logit(p)) = probs = e^logit(p)/(1+e^logit(p))
      #multiply the probabilities expit(etaY) by their corresponding weights (N/sum(N) and 
      #add them up to get the expected prevalence
      #minimize the difference between the expected prevalence and the observed prevalence
      value <- abs(sum(expit(etaY) * (N/sum(N))) - rhoY)
      return(value)
    }
  
  expandCatX <-
    function(X, expandX="all")
    {
      ##
      if(expandX == "all")
        colIndex <- 1:ncol(X)
      if(expandX != "all")
        colIndex <- c(1:ncol(X))[is.element(colnames(X), expandX)]
      
      ## Assumes a check has been performed to make sure the columns that are to be expanded
      ## adhere to the {0,1,2,...} coding convention
      ##
      value <- matrix(X[,1], ncol=1, dimnames=list(1:nrow(X), "Int"))
      ##
      if(ncol(X) > 1)
      {
        n.lev <- unlist(lapply(apply(X, 2, unique), FUN=length))
        for(i in 2:ncol(X))
        {
          if(!is.element(i, colIndex))
          {
            value <- cbind(value, X[,i])
            colnames(value)[ncol(value)] <- colnames(X)[i]
          }
          else
          {
            for(j in 1:(n.lev[i]-1))
            {
              value <- cbind(value, as.numeric(X[,i] == j))
              if(n.lev[i] == 2) colnames(value)[ncol(value)] <- colnames(X)[i]
              if(n.lev[i] > 2) colnames(value)[ncol(value)] <- paste(colnames(X)[i], ".", j, sep="")
            }
          }
        }
      }
      ##
      return(value)
    }
  
  #END OF NEW HANEUSE CODE
  ######################################################################################################
  #END OF NEW HANEUSE CODE
  
  #THIS IS HANEUSE'S EXAMPLE OF HOW TO USE HIS CODE
  # ##
  # data("Ohio")
  # 
  # ##
  # ## glm below is a multitarget version of logistic regression 
  # #with target using levels of cases and controls instead of single target = logit(p) = log(p/(1-p))
  # ## if prevalence = p=prob(Death), then logit(p) = log(p/(1-p)) = linearmodel(factor(Age) + Sex + Race)
  # ## factor(Age) constructs dummies based on age
  # #XM is the Design matrix of the logistic regression model (see below) using category levels (not dummies)
  # XM   <- cbind(Int=1, Ohio[,1:3]) #XM includes the intercept
  # fitM <- glm(cbind(Death, N-Death) ~ factor(Age) + Sex + Race, data=Ohio,
  #             family=binomial)
  # 
  # ## Overall prevalence in the observed data
  # ##
  # ## Note that the sum of the levels of cases and controls are here used to calculate this prevalence
  # sum(Ohio$Death)/sum(Ohio$N)
  # 
  # ## Intercept corresponding to the original vector of log-odds ratios
  # ## 
  # 
  # #beta0 is the intercept optimizing function (beta0 is the intercept)
  # #X=XM is the Design matrix for the logistic regression model  using category levels (not dummies)
  # #The first column should correspond to the intercept data (here, a column of ones, like in lm)
  # #For each exposure (Age, Sex, Race)
  # #the baseline group should be coded as 0, the first level as 1, and so on (no dummies used for XM)
  # #betaX=betaXm is the vector of betas=coefficients of logistic regression excluding the intercept. 
  # #These coefficients  are calculated using as input the dummies of XM (by factorization of Age)
  # #the beta0 calculated by the beta0() function
  # #turns out to be very similar to the original (intercept) output of fitM$coef
  # 
  # # fitM$coef
  # # beta0(betaX=fitM$coef[-1], X=XM, N=Ohio$N, rhoY=sum(Ohio$Death)/sum(Ohio$N))
  # 
  # 
  # # ## Reduction of Sex effect by 50%
  # # ##
  # # betaXm    <- fitM$coef[-1]
  # # betaXm[3] <- betaXm[3] * 0.5
  # # 
  # # #if you reduce the size of the beta of sex by 50%
  # # #the new intercept beta0 is much lower
  # # 
  # # beta0(betaX=betaXm, X=XM, N=Ohio$N, rhoY=sum(Ohio$Death)/sum(Ohio$N))
  # # 
  # # ## Doubling of Race effect
  # # ##
  # # betaXm    <- fitM$coef[-1]
  # # betaXm[4] <- betaXm[4] * 2
  # # 
  # # 
  # # beta0(betaX=betaXm, X=XM, N=Ohio$N, rhoY=sum(Ohio$Death)/sum(Ohio$N))
  # 
  # ##################################################################################################
  # 
  # #Using the modified beta00() instead of the original beta0()
  # cols <- colnames(x=XM) #includes Int
  # beta00(betaX=fitM$coef[-1], X=XM, N=Ohio$N, rhoY=sum(Ohio$Death)/sum(Ohio$N), expandX = cols)
  # beta0(betaX=fitM$coef[-1], X=XM, N=Ohio$N, rhoY=sum(Ohio$Death)/sum(Ohio$N)) #should be the same as with beta00
  # 
  # ##################################################################################################
  
  
  #NOW USING THE MODIFIED BETA00() FUNCTION TO GET THE INTERCEPT WITH HECK'S DATA
  
  #prepare the inputs for the BETA00() FUNCTION that gets our intercept
  #need: betaXM, XM, HeckDF$N, prevalence and cols
  
  drop <- c("N")
  XM = HeckDF[,!(names(HeckDF) %in% drop)] #without N, with the intercept Int
  
  drop <- c("benzene")
  XM_reduced = XM[,!(names(XM) %in% drop)] #without benzene and without N, with the intercept Int
  cols <- colnames(XM_reduced) #getting the names of the categorical features, including int, without benzene
  
  #prepare the betas in betaXM
  bet <- log(c(1.47,1.53,1.57,2.62,0.78,0.75,0.75,0.95,1.94)) #Heck's Odds Ratios transformed into betas
  nam <- c('socioeco.1','socioeco.2','socioeco.3','socioeco.4','race.1','race.2', 'birthplace', 'parity', 'benzene')
  betaXm <- as.data.frame(bet)
  betaXm <- as.data.frame(t(betaXm), n.col=9)
  colnames(betaXm) <- nam
  
  #Upload the 60 factuals and 60 counterfactuals into XM_uplift (total 120 rows)
  XM_uplift <- read.csv("XM_Uplift_continuous.csv") #with the intercept Int
  
  #simulation of benzene data in accordance with Heck's distribution and dividing by IQR as Heck does
  
  
  if (distr == 'stove'){
    
    ############################################################################USE THIS
    #stoves' benzene distribution as per literature (causes our_intercept to not converge)
    mult <- 1 #vs 1/1000 (1 ppm = 1000 ppb, so if mult=1 we are using ppb, as Heck does)
    IQR <- 1.197 #Heck's benzene distribution IQR
    len <- (XM_uplift$benzene) #120 benzene simulations
    minimum <- (0.7/IQR) * mult #min max of benzene divided by IQR of Heck's distribution
    maximum <- (11.5*1000/IQR) * mult 
    ave = (((11.5*1000-0.7)/2)/IQR) * mult
    desv = ((0.75*(11.5*1000-0.7) - 0.25*(11.5*1000-0.7))/IQR) * mult
    benzene_sim <- runif(len, min=minimum, max=maximum)
    benzene_sim <- rnorm(len, mean=ave, sd=desv) 
    benzene_sim[benzene_sim<minimum]<- minimum
    benzene_sim[benzene_sim>maximum]<- maximum
    our_split <- 0.903/IQR* mult
    
  }
  
  
  if (distr == 'indoor'){
    
    ############################################################################USE THIS
    #indoors' benzene distribution as per literature (gives an ATE above zero)
    mult <- 1 #vs 1/1000 (1 ppm = 1000 ppb, so if mult=1 we are using ppb, as Heck does)
    IQR <- 1.197 #Heck's benzene distribution IQR
    len <- (XM_uplift$benzene) #120 benzene simulations
    minimum <- (0.7/IQR) * mult #min max of benzene divided by IQR of Heck's distribution
    maximum <- (12/IQR) * mult 
    ave = (((12-0.7)/2)/IQR) * mult
    desv = ((0.75*(12-0.7) - 0.25*(12-0.7))/IQR) * mult
    benzene_sim <- runif(len, min=minimum, max=maximum)
    benzene_sim <- rnorm(len, mean=ave, sd=desv) 
    benzene_sim[benzene_sim<minimum]<- minimum
    benzene_sim[benzene_sim>maximum]<- maximum
    our_split <- 0.903/IQR* mult
    
  }
  
  if (distr == 'outdoor'){
    
    ###########################################################################USE THIS
    #Heck's benzene distribution (gives an ATE of zero)
    mult <- 1 #vs 1/1000 (1 ppm = 1000 ppb, so if mult=1 we are using ppb, as Heck does)
    IQR <- 1.197 #Heck's benzene distribution IQR
    len <- (XM_uplift$benzene) #120 benzene simulations
    minimum <- ( 0.151/IQR) * mult  #min max of benzene divided by IQR
    maximum <- (4.60/IQR) * mult
    quart_10 <- (0.410/IQR) * mult 
    quart_25 <- (0.591/IQR) * mult 
    quart_75 <- (1.788/IQR) * mult 
    quart_90 <- (2.574/IQR) * mult 
    ave <- (1.268/IQR) * mult 
    desv <- (0.830/IQR) * mult   
    benzene_sim <- runif(len, min=minimum, max=maximum)
    benzene_sim <- rnorm(len, mean=ave, sd=desv) 
    benzene_sim[benzene_sim<minimum]<- minimum
    benzene_sim[benzene_sim>maximum]<- maximum
    our_split <- 0.903/IQR* mult
  }
  
  
  
  #add the first 60 simulated continuous benzene to XM$benzene
  XM$benzene <- benzene_sim[1:60]
  
  #HAVING PREPARED THE INPUTS
  #USE BETA00() FUNCTION TO GET THE INTERCEPT WITH HECK'S DATA
  
  #just checking that this produces the correct dummies
  #XM_expanded <- expandCatX(XM, cols) s
  
  #Using the modified beta00() instead of the original beta0() to get the intercept we are looking for
  our_intercept = beta00(betaX=betaXm, X=XM, N=HeckDF$N, rhoY=prevalence, expandX=cols)
  
  
  ######################################################################################################
  #FACTUALS AND COUNTERFACTUALS TABLE
  
  len <- length(XM_uplift$benzene)/2 #60
  XM_uplift$benzene[1:len] <- minimum #first 60 are factuals
  XM_uplift$benzene[(len+1):(len*2)] <- maximum #second 60 are counterfactuals
  
  
  #calculate PROBaml and corresponding BINaml 
  #for the factuals and counterfactuals table
  
  #first get the our_logit (=the predicted "log of the Odds Ratios of AML")
  X <- expandCatX(XM_uplift, cols) #construct the dummy variables (one-hot-encoding the categorical levels)
  
  Xmatrix <- as.matrix(X)
  COEFvector = unlist(c(our_intercept, betaXm)) #coefficient vector including the intercept we discovered
  our_logit  <- as.vector(Xmatrix%*% COEFvector) #the logit(p)=log(p/(1-p))
  
  
  #now transform our_logit into probabilities
  #expit(logit(p)) = e^logit(p)/(1+e^logit(p)) = probabilities
  #if expit() breaks (because mult>=100) use plogis()
  #PROBaml <-  expit(our_logit)  #leads to error if mult >=100
  PROBaml <- plogis(our_logit)
  
  #now transform probability PROBaml to binary BINaml
  BINaml <- vector(length=length(PROBaml))
  BINaml <- rbinom(length(PROBaml),1,PROBaml) 
  
  #df will be the factuals and counterfactuals table
  df = as.data.frame(X)
  
  
  #add the probability of AML to the df
  df["PROBaml"] <- PROBaml #used for checking
  #add the corresponding simulated binary AML to the df 
  df["BINaml"] <- BINaml #used for checking
  df["outcome"] <- BINaml
  #df["outcome"] <- PROBaml 
  
  
  #transform min and max benzene to binary as per our_split 
  df[df$benzene >  our_split,]$benzene <- 9991
  df[df$benzene <= our_split,]$benzene <- 9990
  df[df$benzene == 9991,]$benzene <- 1
  df[df$benzene == 9990,]$benzene <- 0
  
  #having checked we can drop
  drop <- c("Int", "BINaml", "PROBaml")
  df <- df[,!(names(df) %in% drop)]
  dfwrite <- df
  
  ##################################################################################################
  
  #UPLIFT MODEL
  
  #To do a weighted average of treatment effect calculate the weights
  len <- length(df$outcome)/2
  df$N[1:len] <- HeckDF$N #first 60 are factuals
  df$N[(len+1):(len*2)] <- HeckDF$N #second 60 are counterfactuals
  
  
  #Causality (ATE)
  
  #UNCONDITIONAL
  
  
  #average treatment effect (ATE)
  len <- length(df$outcome)/2  #60
  ITE <- as.vector(matrix(0,nrow=len)) #individual treatment effect
  
  
  if(populationweighted == TRUE){
    
    for (n in 1:len)
    {
      ITE[n] <- (df$outcome[n+len]- df$outcome[n])*(df$N[n]/(.5*sum(df$N))) #individual treatment effect
    }
    ATE <-sum(ITE) #average treatment effect
  }else{
    for (n in 1:len)
    {
      ITE[n] <- df$outcome[n+len]- df$outcome[n] #individual treatment effect
    }
    ATE <- mean(ITE) #average treatment effect
  }
  
  df["ITE"] <- 999
  df$ITE[1:len] <- ITE[1:len]
  dfHOLD <- df
  
  ##################################################################################################
  
  
  
  #CONDITIONAL
  
  #race.X conditional average treatment effect
  len <- length(df[df$race.1 > 0,]$outcome)/2 # 20
  CITE <- as.vector(matrix(0,nrow=len)) #conditional individual treatment effect
  
  
  if(populationweighted == TRUE){
    for (n in 1:len)
    { 
      CITE[n] <- (df[df$race.1 >0,]$outcome[n+len]- df[df$race.1 > 0,]$outcome[n]) * (df[df$race.1 > 0,]$N[n])/(.5*sum(df[df$race.1 > 0,]$N)) ##race conditional individual treatment effect
    }
    CATE <-sum(CITE) #conditional average treatment effect
  }else{
    for (n in 1:len)
    {
      CITE[n] <- df[df$race.1 >0,]$outcome[n+len]- df[df$race.1 > 0,]$outcome[n] ##race conditional individual treatment effect
    }
    CATE <- mean(CITE) #conditional average treatment effect
  }
  
  
  #print("Race Conditional Average Treatment Effect") 
  #print(CATE) 
  
  print("Average Treatment Effect") 
  print(ATE) 
  
  
  
  if (x < 2)
  {
    dfaccum <- dfwrite
    ATEaccum <- ATE
  }
  else
  {
    dfaccum[(nrow(dfaccum) + 1):(nrow(dfaccum) + nrow(dfwrite)),] <- dfwrite
    ATEaccum = ATEaccum + ATE
  }
  
  print(x)
}

write.csv(dfaccum,"UpliftOutput.csv", row.names = TRUE)
print("Average Treatment Effect") 
print(ATEaccum/x) 


