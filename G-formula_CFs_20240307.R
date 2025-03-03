### 12/10/2023 Counterfactual scenario's

#set workingdirectory 
setwd("C:/Users/User/Documents/PhD/Project2_Treatment_patterns_subgroups/Scripts/2.R")
# load the survival library
library(survival)
library(dplyr)
library(magrittr)
library(doSNOW)
library(foreach)
library(data.table)
library(stats)
library(cfdecomp)
library(splines)
library(openxlsx)
#load source script with functions
source("Functions_v5.R")

# read in the data
setwd("C:/Users/User/Documents/PhD/Project2_Treatment_patterns_subgroups/Scripts/2.R/Data_from_scripts")
sample.dat <- read.csv("6.DF_long_20240307.csv")

#save minimum and maximum startage (will represent 16 and 45, but is standardized)
min.startage <- min(sample.dat$startage)
max.startage <- max(sample.dat$startage)

#! PICK SAMPLE
#sample.dat<-sample.dat[sample.dat$anopat %in% sample(unique(sample.dat$anopat),1000),]


#! SELECT ONLY ON MORE TURBULENT TREATMENT STEPS
#sample.dat <- sample.dat[sample.dat$startstep %in% 2:4,]
                         

#print variable list
names(sample.dat)

#set wd to new folder containing date and time.
setwd("C:/Users/User/Documents/PhD/Project2_Treatment_patterns_subgroups/Scripts/2.R")
new_folder_path <- file.path(paste0(getwd(), "/output"), paste0(format(Sys.time(), "%Y.%m.%d_%H.%M.%S"), "_CF"))
dir.create(new_folder_path)
setwd(new_folder_path)

sample.dat$startstep <- as.character(sample.dat$startstep)
sample.dat$startyear <- as.integer(sample.dat$startyear)
sample.dat$steptimer <- as.integer(sample.dat$steptimer)
sample.dat$steptimer0 <- as.integer(sample.dat$steptimer0)
sample.dat$steptimer1 <- as.integer(sample.dat$steptimer1)
sample.dat$steptimer2 <- as.integer(sample.dat$steptimer2)
sample.dat$steptimer3 <- as.integer(sample.dat$steptimer3)
sample.dat$steptimer4 <- as.integer(sample.dat$steptimer4)
sample.dat$l.steptimer <- as.integer(sample.dat$l.steptimer)
sample.dat$l.steptimer0 <- as.integer(sample.dat$l.steptimer0)
sample.dat$l.steptimer1 <- as.integer(sample.dat$l.steptimer1)
sample.dat$l.steptimer2 <- as.integer(sample.dat$l.steptimer2)
sample.dat$l.steptimer3 <- as.integer(sample.dat$l.steptimer3)
sample.dat$l.steptimer4 <- as.integer(sample.dat$l.steptimer4)
sample.dat$switch <- as.integer(sample.dat$switch)
sample.dat$l.switch <- as.integer(sample.dat$l.switch)
sample.dat$dur.switch <- as.integer(sample.dat$dur.switch)
sample.dat$l.dur.switch <- as.integer(sample.dat$l.dur.switch)
#startage should stay continuous as it is standardized.

#save in which column is which info stored
col.index.from <- NULL        
col.index.from[1] <- grep(paste0('^','step0','$'), colnames(sample.dat))
col.index.from[2] <- grep(paste0('^','step1','$'), colnames(sample.dat))
col.index.from[3] <- grep(paste0('^','step2','$'), colnames(sample.dat))
col.index.from[4] <- grep(paste0('^','step3','$'), colnames(sample.dat))
col.index.from[5] <- grep(paste0('^','step4','$'), colnames(sample.dat))
col.index.from[6] <- grep(paste0('^','switch','$'), colnames(sample.dat))
col.index.from[7] <- grep(paste0('^','dur.switch','$'), colnames(sample.dat))
col.index.from[8] <- grep(paste0('^','steptimer','$'), colnames(sample.dat))
col.index.from[9] <- grep(paste0('^','steptimer0','$'), colnames(sample.dat))
col.index.from[10] <- grep(paste0('^','steptimer1','$'), colnames(sample.dat))
col.index.from[11] <- grep(paste0('^','steptimer2','$'), colnames(sample.dat))
col.index.from[12] <- grep(paste0('^','steptimer3','$'), colnames(sample.dat))
col.index.from[13] <- grep(paste0('^','steptimer4','$'), colnames(sample.dat))
col.index.to <- NULL
col.index.to[1] <- grep(paste0('^','l.step0','$'), colnames(sample.dat))
col.index.to[2] <- grep(paste0('^','l.step1','$'), colnames(sample.dat))
col.index.to[3] <- grep(paste0('^','l.step2','$'), colnames(sample.dat))
col.index.to[4] <- grep(paste0('^','l.step3','$'), colnames(sample.dat))
col.index.to[5] <- grep(paste0('^','l.step4','$'), colnames(sample.dat))
col.index.to[6] <- grep(paste0('^','l.switch','$'), colnames(sample.dat))
col.index.to[7] <- grep(paste0('^','l.dur.switch','$'), colnames(sample.dat))
col.index.to[8] <- grep(paste0('^','l.steptimer','$'), colnames(sample.dat))
col.index.to[9] <- grep(paste0('^','l.steptimer0','$'), colnames(sample.dat))
col.index.to[10] <- grep(paste0('^','l.steptimer1','$'), colnames(sample.dat))
col.index.to[11] <- grep(paste0('^','l.steptimer2','$'), colnames(sample.dat))
col.index.to[12] <- grep(paste0('^','l.steptimer3','$'), colnames(sample.dat))
col.index.to[13] <- grep(paste0('^','l.steptimer4','$'), colnames(sample.dat))

#set variables
nid <- length(unique(sample.dat$anopat))
bssize <- 300 # = bootstrapping tim; small (e.g. 2) for trial purposes, then to 75/100 and for the final version to 500; to give stable confidence intervals.
mcsize <- 25  # = how many times Monte carlo is performed, for trial purposes 5, if everything works than 25

tunits <- 24  #number of time units, so in this case: 24 months(=30 days)


#the first value of each variable, just for trying out by hand
# b <- 1
# m <- 1
# t <- 2
# comorb <- "startage"

#Loop throught all comorbidities
comorb_list<-names(sample.dat)[7:16] #use [7:16] to include age, sex, and all comorbidities

for (comorb in comorb_list) {  

  print(paste(comorb, match(comorb, comorb_list), "/", length(comorb_list), Sys.time(), sep=" "))


  formula.outcome <- c(paste("outcome ~
  
                     # time constants #
  
                     splines::ns(startyear, df=4) +
                     startage +  # splines::ns(startage, df=4) 
                     startstep +     # coded as factor
                     sex +
  
                     diabetes + CVD + MHP +  arthritis +
  
                     hypothyroid + ATD + GORD + immunocompromised +
                     
                     # the above could be interacted with the steps
                     # but this will not have an effect on -MARGINAL- values over time
                     # so I have not done this yet
  
                     # time varying #
                     
                     splines::ns(time,df=5) +", 
                     
                     comorb, "*(l.steptimer +
  
                     # l.step1 +  # step 1 is ref
                     
                     # l.step1*l.steptimer + # interaction ref
  
                     l.step2 + l.step2*l.steptimer  +
  
                     l.step3 + l.step3*l.steptimer  +
  
                     l.step4 + l.step4*l.steptimer  +
  
                     l.step0 + l.step0*l.steptimer )"))
  
  formula.step0 <- as.formula(sub('outcome','step0',formula.outcome))
  formula.step2 <- as.formula(sub('outcome','step2',formula.outcome)) # Step 1 is now used as ref. category. Theoretically it does not matter which cat. is used for that, 
  formula.step3 <- as.formula(sub('outcome','step3',formula.outcome)) # but with bootstrapping the ref. cat. shouldn't be too small. Step 0 is empty at start, therefore step 1 seems more suitable.
  formula.step4 <- as.formula(sub('outcome','step4',formula.outcome))  

  # make arrays and matrices to put results of the Monte Carlo(mc) and Bootstrapping(bs) into
  step1toend.ex1.mc <- step1toend.ex2.mc <- matrix(rep(NA,5*1),ncol=5)
  step2toend.ex1.mc <- step2toend.ex2.mc <- matrix(rep(NA,5*1),ncol=5)
  step3toend.ex1.mc <- step3toend.ex2.mc <- matrix(rep(NA,5*1),ncol=5)
  step4toend.ex1.mc <- step4toend.ex2.mc <- matrix(rep(NA,5*1),ncol=5)
  colnames(step1toend.ex1.mc) <- colnames(step2toend.ex1.mc) <- 
    colnames(step3toend.ex1.mc) <- colnames(step4toend.ex1.mc) <-   
    colnames(step1toend.ex2.mc) <- colnames(step2toend.ex2.mc) <-  
    colnames(step3toend.ex2.mc) <- colnames(step4toend.ex2.mc) <- 
    c("step 0", "step 1", "step 2", "step 3", "step 4")
  # bootstrap for transition matrix
  step1toend.ex1.bs <- step1toend.ex2.bs <- matrix(rep(NA,5*bssize),ncol=5)
  step2toend.ex1.bs <- step2toend.ex2.bs <- matrix(rep(NA,5*bssize),ncol=5)
  step3toend.ex1.bs <- step3toend.ex2.bs <- matrix(rep(NA,5*bssize),ncol=5)
  step4toend.ex1.bs <- step4toend.ex2.bs <- matrix(rep(NA,5*bssize),ncol=5)
  colnames(step1toend.ex1.bs) <- colnames(step2toend.ex1.bs) <- 
    colnames(step3toend.ex1.bs) <- colnames(step4toend.ex1.bs) <- 
    colnames(step1toend.ex2.bs) <- colnames(step2toend.ex2.bs) <- 
    colnames(step3toend.ex2.bs) <- colnames(step4toend.ex2.bs) <- 
    c("step 0", "step 1", "step 2", "step 3", "step 4")
  
  # step proportions
  tsteps.ex1.mc <- tsteps.ex2.mc <- matrix(rep(NA,5*tunits),nrow=tunits, ncol=5)
  
  # create matrix for results of..
  ## - Duration till first switch (if any)
  dur.sw.ex1.mc <- dur.sw.ex2.mc <- matrix(rep(NA,2*1), ncol=2)
  dur.sw.ex1.bs <- dur.sw.ex2.bs <- matrix(rep(NA,2*bssize), ncol=2)
  colnames(dur.sw.ex1.mc) <- colnames(dur.sw.ex2.mc) <- c("mean", "sd")
  colnames(dur.sw.ex1.bs) <- colnames(dur.sw.ex2.bs) <- c("mean", "sd")
  ## - Number of switches
  num.sw.ex1.mc <- num.sw.ex2.mc <- matrix(rep(NA,2*1), ncol=2)
  num.sw.ex1.bs <- num.sw.ex2.bs <- matrix(rep(NA,2*bssize), ncol=2)
  colnames(num.sw.ex1.mc) <- colnames(num.sw.ex2.mc) <- c("mean", "sd")
  colnames(num.sw.ex1.bs) <- colnames(num.sw.ex2.bs) <- c("mean", "sd")
  ## - Switched yes or no
  switching.yn.ex1.mc <- switching.yn.ex2.mc <- matrix(rep(NA,2*1), ncol=2)
  switching.yn.ex1.bs <- switching.yn.ex2.bs <- matrix(rep(NA,4*bssize), ncol=4)
  colnames(switching.yn.ex1.mc) <- colnames(switching.yn.ex2.mc) <- c("switch", "no.switch")
  colnames(switching.yn.ex1.bs) <- colnames(switching.yn.ex2.bs) <- c("switch", "no.switch", "25th.switch", "75th.switch")
  
  ## let's start the actual G-formula: ##
  for (b in 1:bssize) {
    
    #register time #1 to calculate time to run
    t1 <- Sys.time() 
    #print(paste(b, "b", t1, sep=" "))
    
    # sample individuals (with replacement) from ex.dat
    ex.sample <- cluster.resample(sample.dat, "anopat")
    
    # (re)fit models to ex.sample
    # bsf = bootstrap fit
       # step2 and step3 might not converge as the sample is very very small, in some bootstrap iterations
  
      fit.bsf.step0 <- glm(formula.step0, family=binomial, data=ex.sample)
    fit.bsf.step2 <- glm(formula.step2, family=binomial, data=ex.sample,subset=ex.sample$step0==0)
    fit.bsf.step3 <- glm(formula.step3, family=binomial, data=ex.sample,subset=ex.sample$step0==0 & ex.sample$step2==0)
    fit.bsf.step4 <- glm(formula.step4, family=binomial, data=ex.sample,subset=ex.sample$step0==0 & ex.sample$step2==0 & ex.sample$step3==0)

    t2 <- Sys.time()
    #print(paste(b, "model fitting", t2 - t1, sep=" "))
    
    #Open cluster with multiple cores for doPar
    #Never use all 9, but only 3 to 5. Using 4 decreased the time to run from ~1.20 mins to <10 sec.
    cl <- makeCluster(4, type="SOCK")
    registerDoSNOW(cl)
    
    #Monte-Carlo simulation loop
    MC_XX = foreach(m=1:mcsize) %dopar% {
      
      # take individuals at time 1
      # and discard the other observations
      # ex.sample 1 represents the natural course (hopefully)
      ex.sample.1 <- ex.sample.2 <- sample.dat[sample.dat$time==1,]
      ex.sample.1$id <- ex.sample.2$id <- 1:nid
      
      # now let's create a counterfactual dataset where we intervene on the data in some way
      # we will have cured the comorbidity, so NO ONE HAS THE COMORBIDITY
      if (comorb == "startage") {
        ex.sample.1[comorb] <- min.startage
      } else {
      ex.sample.1[comorb] <- 0 
      }
      
      # start a loop that moves through the follow-up time units
      # this part of the g-formula tries to reproduce the data for the hypothetical scenario
      ex.sample.1.temp <- ex.sample.1
      
      #save step proportions at t1
      tsteps.ex1.mc[1,] <- c(sum(ex.sample.1.temp$step0),sum(ex.sample.1.temp$step1),sum(ex.sample.1.temp$step2),
                              sum(ex.sample.1.temp$step3),sum(ex.sample.1.temp$step4))/nrow(ex.sample.1.temp)
      # looks like the empirical data, so good

      for(t in 2:tunits) {
        
        # make a copy of the previous row and then update time
        ex.sample.1.temp$time <- ex.sample.1.temp$time+1
        
        # lag values
        ex.sample.1.temp[,col.index.to] <- ex.sample.1.temp[,col.index.from]
        
        # simulate time varying variables           #!Warnings: did not converge / fitted probabilities 0 or 1 occur
        # treatment step                              can happen with bootstrapping and sampling, is not problematic in this case.
        x <- multinomial.simulation.5var(fit.bsf.step0,fit.bsf.step2,fit.bsf.step3,
                                         fit.bsf.step4, ex.sample.1.temp)
        ex.sample.1.temp$step0 <- x[,1]
        ex.sample.1.temp$step2 <- x[,2]
        ex.sample.1.temp$step3 <- x[,3]
        ex.sample.1.temp$step4 <- x[,4]
        ex.sample.1.temp$step1 <- x[,5]     # notice th at the 'other' category here is 1, hence it gets the last column of simulated values
        rm(x)
        
        # save some info on switches
        ex.sample.1.temp$switch <- ifelse(    ex.sample.1.temp$step0!=ex.sample.1.temp$l.step0 |
                                               ex.sample.1.temp$step1!=ex.sample.1.temp$l.step1 |
                                               ex.sample.1.temp$step2!=ex.sample.1.temp$l.step2 |
                                               ex.sample.1.temp$step3!=ex.sample.1.temp$l.step3 |
                                               ex.sample.1.temp$step4!=ex.sample.1.temp$l.step4,
                                               ex.sample.1.temp$l.switch+1, ex.sample.1.temp$l.switch)
        ex.sample.1.temp$dur.switch <- ifelse(ex.sample.1.temp$switch==1 & is.na(ex.sample.1.temp$dur.switch), 
                                              ex.sample.1.temp$time, ex.sample.1.temp$dur.switch)
        
        #set steptimer
        ex.sample.1.temp$steptimer <- ifelse(ex.sample.1.temp$switch!=ex.sample.1.temp$l.switch,1,ex.sample.1.temp$l.steptimer+1)
        
        #save step proportions at t
        tsteps.ex1.mc[t,] <- c(sum(ex.sample.1.temp$step0),sum(ex.sample.1.temp$step1),sum(ex.sample.1.temp$step2),
                            sum(ex.sample.1.temp$step3),sum(ex.sample.1.temp$step4))/nrow(ex.sample.1.temp)
        
        # join this newly created data with the old data
        ex.sample.1 <- rbind(ex.sample.1,ex.sample.1.temp)
        
      }
      
      # now let's create a counterfactual dataset where we intervene on the data in some way
      # we will have cured the comorbidity, so EVERYONE HAS THE COMORBIDITY
      
      if (comorb == "startage") {
        ex.sample.2[comorb] <- max.startage
      } else {
        ex.sample.2[comorb] <- 1
      }
      
      # start a loop that moves through the follow-up time units
      # this part of the g-formula tries to reproduce the data for the hypothetical scenario
      ex.sample.2.temp <- ex.sample.2
      
      #save step proportions at t1
      tsteps.ex2.mc[1,] <- c(sum(ex.sample.2.temp$step0),sum(ex.sample.2.temp$step1),sum(ex.sample.2.temp$step2),
                             sum(ex.sample.2.temp$step3),sum(ex.sample.2.temp$step4))/nrow(ex.sample.2.temp)
      # looks like the empirical data, so good
      
      for(t in 2:tunits) {
        
        # make a copy of the previous row and then update time
        ex.sample.2.temp$time <- ex.sample.2.temp$time+1
        
        # lag values
        ex.sample.2.temp[,col.index.to] <- ex.sample.2.temp[,col.index.from]
        
        # simulate time varying variables           ##!Warnings: did not converge / fitted probabilities 0 or 1 occur - not problematic here
        # treatment step
        x <- multinomial.simulation.5var(fit.bsf.step0,fit.bsf.step2,fit.bsf.step3,
                                         fit.bsf.step4,ex.sample.2.temp)
        ex.sample.2.temp$step0 <- x[,1]
        ex.sample.2.temp$step2 <- x[,2]
        ex.sample.2.temp$step3 <- x[,3]
        ex.sample.2.temp$step4 <- x[,4]
        ex.sample.2.temp$step1 <- x[,5]     # notice th at the 'other' category here is 1, hence it gets the last column of simulated values
        rm(x)
        
        # save some info on switches
        # this is now also not updated for people who were censored; but that shouldn't happen in the real data either anyway
        ex.sample.2.temp$switch <- ifelse(    ex.sample.2.temp$step0!=ex.sample.2.temp$l.step0 |
                                                ex.sample.2.temp$step1!=ex.sample.2.temp$l.step1 |
                                                ex.sample.2.temp$step2!=ex.sample.2.temp$l.step2 |
                                                ex.sample.2.temp$step3!=ex.sample.2.temp$l.step3 |
                                                ex.sample.2.temp$step4!=ex.sample.2.temp$l.step4,
                                              ex.sample.2.temp$l.switch+1, ex.sample.2.temp$l.switch)
        ex.sample.2.temp$dur.switch <- ifelse(ex.sample.2.temp$switch==1 & is.na(ex.sample.2.temp$dur.switch), 
                                              ex.sample.2.temp$time, ex.sample.2.temp$dur.switch)
        
        #set steptimer
        ex.sample.2.temp$steptimer <- ifelse(ex.sample.2.temp$switch!=ex.sample.2.temp$l.switch,1,ex.sample.2.temp$l.steptimer+1)
        
        #save step proportions at t
        tsteps.ex2.mc[t,] <- c(sum(ex.sample.2.temp$step0),sum(ex.sample.2.temp$step1),sum(ex.sample.2.temp$step2),
                               sum(ex.sample.2.temp$step3),sum(ex.sample.2.temp$step4))/nrow(ex.sample.2.temp)

        # join this newly created data with the old data
        ex.sample.2 <- rbind(ex.sample.2,ex.sample.2.temp)
        
      }
      
      # save dummy information into categorical column, for transition matrix info
       ex.sample.1$steptvar <- ifelse(ex.sample.1$step0==1,0,
                                   ifelse(ex.sample.1$step1==1,1,
                                          ifelse(ex.sample.1$step2==1,2,
                                                 ifelse(ex.sample.1$step3==1,3,
                                                        ifelse(ex.sample.1$step4==1,4,
                                                                             NA)))))
      ex.sample.1$steptvar <- as.factor(ex.sample.1$steptvar)

      ex.sample.2$steptvar <- ifelse(ex.sample.2$step0==1,0,
                                     ifelse(ex.sample.2$step1==1,1,
                                            ifelse(ex.sample.2$step2==1,2,
                                                   ifelse(ex.sample.2$step3==1,3,
                                                          ifelse(ex.sample.2$step4==1,4,
                                                                        NA)))))
      ex.sample.2$steptvar <- as.factor(ex.sample.2$steptvar)
      
      #same for lagged columns
      ex.sample.1$l.steptvar <- ifelse(ex.sample.1$l.step0==1,0,
                                     ifelse(ex.sample.1$l.step1==1,1,
                                            ifelse(ex.sample.1$l.step2==1,2,
                                                   ifelse(ex.sample.1$l.step3==1,3,
                                                          ifelse(ex.sample.1$l.step4==1,4,
                                                                               NA)))))
      ex.sample.1$l.steptvar <- as.factor(ex.sample.1$l.steptvar)

      ex.sample.2$l.steptvar <- ifelse(ex.sample.2$l.step0==1,0,
                                       ifelse(ex.sample.2$l.step1==1,1,
                                              ifelse(ex.sample.2$l.step2==1,2,
                                                     ifelse(ex.sample.2$l.step3==1,3,
                                                            ifelse(ex.sample.2$l.step4==1,4,
                                                                          NA)))))
      ex.sample.2$l.steptvar <- as.factor(ex.sample.2$l.steptvar)
      
      # step at t=1 vs step at t=tunits
      step1toend.ex1.mc[1,] <- table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==1])/sum(table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==1]))
      step2toend.ex1.mc[1,] <- table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==2])/sum(table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==2]))
      step3toend.ex1.mc[1,] <- table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==3])/sum(table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==3]))
      step4toend.ex1.mc[1,] <- table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==4])/sum(table(ex.sample.1$steptvar[ex.sample.1$time==tunits & ex.sample.1$step==4]))

      step1toend.ex2.mc[1,] <- table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==1])/sum(table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==1]))
      step2toend.ex2.mc[1,] <- table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==2])/sum(table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==2]))
      step3toend.ex2.mc[1,] <- table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==3])/sum(table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==3]))
      step4toend.ex2.mc[1,] <- table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==4])/sum(table(ex.sample.2$steptvar[ex.sample.2$time==tunits & ex.sample.2$step==4]))

      #save time of first switch (excluding no switches (is na))
      dur.sw.ex1.mc[1,1]<-mean(ex.sample.1$dur.switch[ex.sample.1$time==tunits], na.rm=TRUE)
      dur.sw.ex1.mc[1,2]<-sd(ex.sample.1$dur.switch[ex.sample.1$time==tunits], na.rm=TRUE)

      dur.sw.ex2.mc[1,1]<-mean(ex.sample.2$dur.switch[ex.sample.2$time==tunits], na.rm=TRUE)
      dur.sw.ex2.mc[1,2]<-sd(ex.sample.2$dur.switch[ex.sample.2$time==tunits], na.rm=TRUE)

      #save number of switches
      num.sw.ex1.mc[1,1]<-mean(ex.sample.1$switch[ex.sample.1$time==tunits  & ex.sample.1$switch!=0], na.rm=TRUE)
      num.sw.ex1.mc[1,2]<-sd(ex.sample.1$switch[ex.sample.1$time==tunits & ex.sample.1$switch!=0], na.rm=TRUE)

      num.sw.ex2.mc[1,1]<-mean(ex.sample.2$switch[ex.sample.2$time==tunits & ex.sample.2$switch!=0], na.rm=TRUE)
      num.sw.ex2.mc[1,2]<-sd(ex.sample.2$switch[ex.sample.2$time==tunits & ex.sample.2$switch!=0], na.rm=TRUE)

      #Save number of patients switching vs. not switching
      switching.yn.ex1.mc[1,] <- table(ex.sample.1$switch[ex.sample.1$time==24]==0)     
      switching.yn.ex2.mc[1,] <- table(ex.sample.2$switch[ex.sample.2$time==24]==0)     
      
      return(list(
        step1toend.ex1.mc,
        step2toend.ex1.mc,
        step3toend.ex1.mc,
        step4toend.ex1.mc,
        dur.sw.ex1.mc,
        num.sw.ex1.mc,
        switching.yn.ex1.mc,
        step1toend.ex2.mc,
        step2toend.ex2.mc,
        step3toend.ex2.mc,
        step4toend.ex2.mc,
        dur.sw.ex2.mc,
        num.sw.ex2.mc,
        switching.yn.ex2.mc,
        tsteps.ex1.mc,
        tsteps.ex2.mc
      ))
      
      }
    
    #close cluster
    stopCluster(cl)
    
    #Create tables from Multicore processing results
    step1toend.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[1]]}))
    step2toend.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[2]]}))
    step3toend.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[3]]}))
    step4toend.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[4]]}))
    dur.sw.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[5]]}))
    num.sw.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[6]]}))
    switching.yn.ex1.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[7]]}))
    step1toend.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[8]]}))
    step2toend.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[9]]}))
    step3toend.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[10]]}))
    step4toend.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[11]]}))
    dur.sw.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[12]]}))
    num.sw.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[13]]}))
    switching.yn.ex2.mc1<-do.call(rbind, lapply(1:mcsize, function(x){MC_XX[[x]][[14]]}))
    
    # Take mean over tables with treatment step distribution per time-point
    # Ex1
    tsteps.ex1.mc1 <- MC_XX[[1]][[15]]                            #set the first table
    for (i in 2:mcsize) { 
      tsteps.ex1.mc1 <- tsteps.ex1.mc1 + MC_XX[[i]][[15]]         #add the others (sum)
    }
    tsteps.ex1.mc1 <- tsteps.ex1.mc1/mcsize                       #divide by mcsize, so proportions will add up to 1
    # Ex 2
    tsteps.ex2.mc1 <- MC_XX[[1]][[16]]                            #set the first table
    for (i in 2:mcsize) { 
      tsteps.ex2.mc1 <- tsteps.ex2.mc1 + MC_XX[[i]][[16]]         #add the others (sum)
    }
    tsteps.ex2.mc1 <- tsteps.ex2.mc1/mcsize                       #divide by mcsize, so proportions will add up to 1
    #add columnames
    colnames(tsteps.ex1.mc1) <- colnames(tsteps.ex2.mc1) <- c("step 0", "step 1", "step 2", "step 3", "step 4")
    #Save tables 
    saveRDS(tsteps.ex1.mc1, paste0("temp.tsteps.ex1.b", b, comorb, ".rds"))
    saveRDS(tsteps.ex2.mc1, paste0("temp.tsteps.ex2.b", b, comorb, ".rds"))
    #remove
    rm(tsteps.ex1.mc1)
    rm(tsteps.ex2.mc1)
    
    #Merge results from all MC loops
    step1toend.ex1.bs[b,] <- apply(step1toend.ex1.mc1, 2, mean)
    step2toend.ex1.bs[b,] <- apply(step2toend.ex1.mc1, 2, mean)
    step3toend.ex1.bs[b,] <- apply(step3toend.ex1.mc1, 2, mean)
    step4toend.ex1.bs[b,] <- apply(step4toend.ex1.mc1, 2, mean)

    step1toend.ex2.bs[b,] <- apply(step1toend.ex2.mc1, 2, mean)
    step2toend.ex2.bs[b,] <- apply(step2toend.ex2.mc1, 2, mean)
    step3toend.ex2.bs[b,] <- apply(step3toend.ex2.mc1, 2, mean)
    step4toend.ex2.bs[b,] <- apply(step4toend.ex2.mc1, 2, mean)

    dur.sw.ex1.bs[b,] <- apply(dur.sw.ex1.mc1, 2, mean)
    dur.sw.ex2.bs[b,] <- apply(dur.sw.ex2.mc1, 2, mean)
    
    num.sw.ex1.bs[b,] <- apply(num.sw.ex1.mc1, 2, mean)
    num.sw.ex2.bs[b,] <- apply(num.sw.ex2.mc1, 2, mean)
    
    switching.yn.ex1.bs[b, 1:2] <- apply(switching.yn.ex1.mc1, 2, mean)
    switching.yn.ex1.bs[b, 3] <- quantile(switching.yn.ex1.mc1[,1], 0.25)
    switching.yn.ex1.bs[b, 4] <- quantile(switching.yn.ex1.mc1[,1], 0.75)
    
    switching.yn.ex2.bs[b, 1:2] <- apply(switching.yn.ex2.mc1, 2, mean)
    switching.yn.ex2.bs[b, 3] <- quantile(switching.yn.ex2.mc1[,1], 0.25)
    switching.yn.ex2.bs[b, 4] <- quantile(switching.yn.ex2.mc1[,1], 0.75)
    
    t5 <- Sys.time()
    print(paste(b, "total b loop", t5 - t1, sep=" "))
  }
  
  # step-to-end.bs omzetten in matrix
  # create matrix
    ex1.matrix <- matrix(data = NA, nrow=4, ncol=5, 
                         dimnames = list(c("step1", "step2", "step3", "step4"),
                                         c("step0", "step1", "step2", "step3", "step4"))) 
    ex2.matrix <- matrix(data = NA, nrow=4, ncol=5, 
                         dimnames = list(c("step1", "step2", "step3", "step4"),
                                         c("step0", "step1", "step2", "step3", "step4"))) 
  # fill with means of bootstrapping results
  ex1.matrix[1,] <- apply(step1toend.ex1.bs, 2, mean)
  ex1.matrix[2,] <- apply(step2toend.ex1.bs, 2, mean)
  ex1.matrix[3,] <- apply(step3toend.ex1.bs, 2, mean)
  ex1.matrix[4,] <- apply(step4toend.ex1.bs, 2, mean)
  ex1.m <- round(ex1.matrix*100, 1)
  ex2.matrix[1,] <- apply(step1toend.ex2.bs, 2, mean)
  ex2.matrix[2,] <- apply(step2toend.ex2.bs, 2, mean)
  ex2.matrix[3,] <- apply(step3toend.ex2.bs, 2, mean)
  ex2.matrix[4,] <- apply(step4toend.ex2.bs, 2, mean)
  ex2.m <- round(ex2.matrix*100, 1)

#################################################################################
  # duration till switch
  # take mean
  dur.sw.ex1 <- apply(dur.sw.ex1.bs, 2, mean, na.rm=TRUE)
  dur.sw.ex1.PI <- apply(dur.sw.ex1.bs, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)
  dur.sw.ex2 <- apply(dur.sw.ex2.bs, 2, mean, na.rm=TRUE)
  dur.sw.ex2.PI <- apply(dur.sw.ex2.bs, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)

  #number of switches
  #mean
  num.sw.ex1 <- apply(num.sw.ex1.bs, 2, mean, na.rm=TRUE)
  num.sw.ex1.PI <- apply(num.sw.ex1.bs, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)
  num.sw.ex2 <- apply(num.sw.ex2.bs, 2, mean, na.rm=TRUE)
  num.sw.ex2.PI <- apply(num.sw.ex2.bs, 2, quantile, prob=c(0.025, 0.975), na.rm=TRUE)

  #switching yes no    - there should be no Inf/NAs here.
  #mean
  sw.yn.ex1 <- apply(switching.yn.ex1.bs, 2, mean)
  sw.yn.ex1.PI <- apply(switching.yn.ex1.bs, 2, quantile, prob=c(0.025, 0.975))
  sw.yn.ex2 <- apply(switching.yn.ex2.bs, 2, mean)
  sw.yn.ex2.PI <- apply(switching.yn.ex2.bs, 2, quantile, prob=c(0.025, 0.975))
  
  
  #to obtain p-values (ex1 = no disease; ex2 = everyone disease)
  num.sign <- pval(num.sw.ex1.bs, num.sw.ex2.bs)
  dur.sign <- pval(dur.sw.ex1.bs, dur.sw.ex2.bs)
  sw.sign  <- pval(switching.yn.ex1.bs, switching.yn.ex2.bs)

  
  #Read in proportions per time point document1
  trsteps.ex1 <- readRDS(paste0("temp.tsteps.ex1.b", b, comorb, ".rds"))
  trsteps.ex2 <- readRDS(paste0("temp.tsteps.ex2.b", b, comorb, ".rds"))
  #for each bootstrap loop, read in the next document, and sum it to the previous ones
  for(i in 2:bssize){
    trsteps.ex1 <- trsteps.ex1 + readRDS(paste0("temp.tsteps.ex1.b", i, comorb, ".rds"))
    trsteps.ex2 <- trsteps.ex2 + readRDS(paste0("temp.tsteps.ex2.b", i, comorb, ".rds"))
  }
  #divide it by bssize, and translate to percentages with 2 decimals
  trsteps.ex1 <- round((trsteps.ex1 / bssize)*100, 2)
  trsteps.ex2 <- round((trsteps.ex2 / bssize)*100, 2)
  #if you want: double check, should be approximately add up to 100
  # print(rowSums(trsteps.ex1))
  # print(rowSums(trsteps.ex2))
  #delete previous documents from WD
  for(i in 1:bssize){
    unlink(paste0("temp.tsteps.ex1.b", i, comorb, ".rds"))
    unlink(paste0("temp.tsteps.ex2.b", i, comorb, ".rds"))
  }

  
  # Save results
  
  #create lists of df's
  df_list1 <- list(ex1.m,
                   as.data.frame(dur.sw.ex1),
                   dur.sw.ex1.PI,
                   as.data.frame(num.sw.ex1),
                   num.sw.ex1.PI,
                   as.data.frame(sw.yn.ex1),
                   sw.yn.ex1.PI,
                   trsteps.ex1)
  df_list2 <- list(ex2.m,
                   as.data.frame(dur.sw.ex2),
                   dur.sw.ex2.PI,
                   as.data.frame(num.sw.ex2),
                   num.sw.ex2.PI,
                   as.data.frame(sw.yn.ex2),
                   sw.yn.ex2.PI,
                   trsteps.ex2)
  Pval <- list(as.data.frame(unlist(dur.sign)), 
               as.data.frame(unlist(num.sign)),
               as.data.frame(unlist(sw.sign)))
  
  #open workbook and add sheets
  wb <- createWorkbook()
  addWorksheet(wb, paste0("no_",comorb))
  addWorksheet(wb, paste(comorb))
  addWorksheet(wb, paste("significance"))
  
  #write data into sheets
  curr_row <- 1
  for(i in seq_along(df_list1)) {
    writeData(wb, paste0("no_", comorb), names(df_list1)[i], startCol = 1, startRow = curr_row, rowNames = TRUE)
    writeData(wb, paste0("no_", comorb), df_list1[[i]], startCol = 1, startRow = curr_row+1, rowNames = TRUE)
    curr_row <- curr_row + nrow(df_list1[[i]]) + 2
  }
  
  curr_row <- 1
  for(i in seq_along(df_list2)) {
    writeData(wb, paste(comorb), names(df_list2)[i], startCol = 1, startRow = curr_row, rowNames = TRUE)
    writeData(wb, paste(comorb), df_list2[[i]], startCol = 1, startRow = curr_row+1, rowNames = TRUE)
    curr_row <- curr_row + nrow(df_list2[[i]]) + 2
  }
  
  curr_row <- 1
  for(i in seq_along(Pval)) {
    writeData(wb, paste("significance"), names(Pval)[i], startCol = 1, startRow = curr_row, rowNames = TRUE)
    writeData(wb, paste("significance"), Pval[[i]], startCol = 1, startRow = curr_row+1, rowNames = TRUE)
    curr_row <- curr_row + nrow(Pval[[i]]) + 2
  }
  
  #save
  saveWorkbook(wb, paste0(comorb,"_", format(Sys.Date(), "%Y%m%d"), ".xlsx"), returnValue = TRUE)
  
  write.csv(dur.sw.ex1.bs, file=paste0(comorb, ".dur.sw.ex1.bs.csv"), row.names=FALSE)
  write.csv(dur.sw.ex2.bs, file=paste0(comorb, ".dur.sw.ex2.bs.csv"), row.names=FALSE)
  write.csv(num.sw.ex1.bs, file=paste0(comorb, ".num.sw.ex1.bs.csv"), row.names=FALSE)
  write.csv(num.sw.ex2.bs, file=paste0(comorb, ".num.sw.ex2.bs.csv"), row.names=FALSE)
  write.csv(switching.yn.ex1.bs, file=paste0(comorb, ".switching.yn.ex1.bs.csv"), row.names=FALSE)
  write.csv(switching.yn.ex2.bs, file=paste0(comorb, ".switching.yn.ex2.bs.csv"), row.names=FALSE)
  
}

print("THE END")
