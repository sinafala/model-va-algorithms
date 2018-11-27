# Model VA Algorithm Functions
# Sam Clark
# 2018-10-30


###### <<Function>>
storeSCI <- function (SCIlist) {
  
  # converts a list of p(s|c) SCI vectors into a labeled matrix
  
  # SCIlist: a list of 5-element SCI vectors, each element betweten 0.0 and 1.0
  # return a labeled cause x p(s|c) matrix
  
  # create matrix from list
  SCI.mat <- matrix(unlist(SCIlist), ncol = 5, byrow = TRUE)
  # label columns
  colnames(SCI.mat) <- c("s1","s2","s3","s4","s5")
  # label rows
  # initialize row labels vector
  sci.rownames <- "cause1"
  # loop through causes
  for (i in 2:nrow(SCI.mat)) {
    # append each cause to row labels vector
    sci.rownames <- c(sci.rownames,paste("cause",i,sep=""))
  }
  rownames(SCI.mat) <- sci.rownames
  
  # return matrix of SCI
  return(SCI.mat)
  
}


###### <<Function>>
simulateDeaths <- function (simDeaths.num,csmf,psc.mat) {
  
  # simulates a set of VA deaths using p(s|c) SCI and CSMFs provided
  
  # simDeaths.num is the number of daeths to simulate
  # csmf is a vector of fractions, one for each cause
  # psc.mat is matrix of conditional probabilities p(s|c)
  # return 
  #   * simDeaths: a labeled deaths x cause,symptoms matrix of simulated deaths
  
  # initialize vector simulated deaaths matrix
  simDeaths <- c()
  #loop over deaths
  for (i in 1:simDeaths.num) {
    # draw a cause, 1 or 2
    cause <- which(rmultinom(1,1,csmf)==1)
    # select he correct set of conditional probabilities for drawing symptoms
    psc <- psc.mat[cause,]
    # draw symptoms conditional on this cause
    # initialize symptoms vector
    symptoms <- c()
    # loop over symptoms
    for (s in 1:5) {
      symptoms <- c(symptoms,rbinom(1,1,psc[s]))
    }
    # add this death to the deaths matrix
    simDeaths <- rbind(simDeaths,c(cause,symptoms))
  }
  # name the columns
  # initialize column names vector
  death.colnames <- c("cause","s1")
  # add symptom labels
  for (s in 2:5) {
    death.colnames <- c(death.colnames,paste("s",s,sep=""))
  }
  colnames(simDeaths) <- death.colnames
  # initialize row names vector
  death.rownames <- "death1"
  # add row labels
  for (i in 2:nrow(simDeaths)) {
    death.rownames <- c(death.rownames,paste("death",i,sep=""))
  }
  rownames(simDeaths) <- death.rownames
  
  # return simulated deaths
  return(simDeaths)
  
}


###### <<Function>>
calcCSMF <- function (deaths) {
  
  # calculate CSMFs for a set of deaths
  
  # deaths is a deaths x cause,symptoms matrix of deaths
  # retuns a causes-element vector of CSMFs
  
  # identify unique causes
  causes <- sort(unique(deaths[,1]))
  # initialize vector of CSMFs
  csmfs <- c()
  # loop over causes 
  for (c in 1:length(causes)) {
    # calculate the fraction of deaths with cause c and add it to the matrix of CSMFs
    csmfs <- c(csmfs,nrow(deaths[which(deaths[,1]==c),])/nrow(deaths))
  }
  # make csmfs result into a matrix
  csmfs <- matrix(csmfs,nrow=1)
  # label columns
  csmfs.colnames <- "cause1"
  for (i in 2:length(causes)) {
    csmfs.colnames <- c(csmfs.colnames,paste("cause",i,sep=""))
  }
  colnames(csmfs) <- csmfs.colnames
  # label rows - subsitute() returns the name of its input
  rownames (csmfs) <- c(substitute(deaths))
  
  # return CSMFs
  return(csmfs)
  
}


###### <<Function>>
calcSCI <- function(deaths) {
  
  # calculate SCI from a set of deaths with causes and symptoms
  
  # deaths is a deaths x cause,symptoms matrix of deaths
  # return a cause x symptoms matrix of conditional probabilities p(s|c)
  
  # identify unique causes
  causes <- sort(unique(deaths[,1]))
  # initialize list of vectors of SCI
  sci <- list()
  # loop over causes
  for (c in 1:length(causes)) {
    # calculate the fraction of deaths with each symptom for this cause
    # add the result vector to the list of SCI
    sci[[c]] <- colSums(deaths[which(deaths[,1]==c),2:6])/nrow(deaths[which(deaths[,1]==c),])
  }
  # convert the list to a matrix
  sci <- matrix(unlist(sci), ncol = 5, byrow = TRUE)
  colnames(sci) <- c("s1","s2","s3","s4","s5")
  # initialize row labels vector
  sci.rownames <- "cause1"
  # loop over causes and create row labels vector
  for (i in 2:length(causes)) {
    # append each cause to row labels vector
    sci.rownames <- c(sci.rownames,paste("cause",i,sep=""))
  }
  rownames(sci) <- sci.rownames
  
  # return cause x symptoms matrix of SCI
  return(sci)
  
}


###### <<Function>>
storeVADeaths <- function (deathList) {
  
  # converts a list of VA death vectors into a labeled matrix
  
  # deathList: a list of 5-element VA symptom vectors, each element 0 or 1
  # return a labeled death x symptoms matrix
  
  # create a matrix with the deaths & signs/symptoms: rows for deaths, columns for signs/symptoms
  deaths.mat <- matrix(unlist(deathList), ncol = 5, byrow = TRUE)
  # label columns
  colnames <- c("s1","s2","s3","s4","s5")
  colnames(deaths.mat) <- colnames
  # label rows
  # initialize row labels vector
  death.rownames <- "death1"
  # loop over deaths
  for (i in 2:nrow(deaths.mat)) {
    # append each death top row labels vector
    death.rownames <- c(death.rownames,paste("death",i,sep=""))
  }
  rownames(deaths.mat) <- death.rownames
  
  # return labeled matrix of ddeaths
  return(deaths.mat)
  
}


###### <<Function>>
calcNBC <- function(SCI,SCI.c,deaths) {
  
  # implements the NBC VA algorithm
  
  # SCI is a cause x symptoms matrix of SCI conditional probabilities p(s|c)
  # SCI.c is p(c) 
  # deaths is a deaths x symptoms matrix of VA deaths containing binary symptom values
  # returns a list containing
  #   * deaths x cause,(largest)cause-probability matrix of NBC-identified causes
  #   * deaths x cause-probabilities matrix of all cause-probabilities
  #   * CSMFs calculated by summing 'fractional deaths' for each cause,  
  #   * CSMFs calculated by calculating the fraction of deaths assigned to each cause
  #   * when: date and time when finished
  #   * deaths: the input deaths
  #   * sci: the input SCI
  
  # NBC.probs are the NBC probabilities for each death: deaths x causes
  NBC.probs <- array(0,dim=c(nrow(deaths),nrow(SCI)))
  
  # NBC.causes are causes for maximum NBC probabilities
  NBC.causes <- data.frame(
    cause=as.character(rep("",nrow(deaths))),
    probability=rep(0,nrow(deaths)),
    stringsAsFactors = FALSE
  )
  rownames(NBC.causes) <- rownames(deaths)
  
  # calculate the NBC conditional probabilities of causes given symptoms vectors, p(c|s)
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # loop over causes and calculate NBC products
    # initialize a vector of probabilities for each cause
    prob.causes <- c()
    for (c in 1:nrow(SCI)) {
      # loop over symptoms and update NBC products
      # initialize the NBC product for this cause
      nbc.product <- 1
      # loop over symptoms
      for (s in 1:ncol(SCI)) {
        # update the NBC product for this symptom
        nbc.product <- nbc.product * SCI[c,s]^deaths[d,s] * (1-SCI[c,s])^(1-deaths[d,s])
      }
      # multiple NBC product by p(c) and add to the results by cause
      prob.causes <- c(prob.causes,SCI.c[c] * nbc.product)
    }
    # now that we have NBC products for each cause, loop over causes again and 
    #   calculate NBC probability for each cause and this death's set of symptoms, 
    #   e.g. normalize the NBC products
    # loop over causes
    for (c in 1:nrow(SCI)) {
      NBC.probs[d,c] <- prob.causes[c] / sum(prob.causes)
    }
    # for each death identify the cause with the largest probability
    # take the first element from the which() statement because sometimes all values 
    #   of NBC.probs[d,] are equal
    NBC.causes[d,1] <- paste("cause",which(NBC.probs[d,]==max(NBC.probs[d,]))[1],sep="")
    NBC.causes[d,2] <- NBC.probs[d,which(NBC.probs[d,]==max(NBC.probs[d,]))[1]]
  }
  
  # store all NBC cause-specific probabilities for each death in a labeled matrix
  NBC.summaries <- NBC.probs
  rownames(NBC.summaries) <- rownames(deaths)
  # initialize vector of column labels, cause-probabilities
  causes.NBC <- c("cause1.probability")
  # loop over causes
  for (c in 2:nrow(SCI)) {
    # append each cause-probability to the column labels
    causes.NBC <- c(causes.NBC,paste("cause",c,".probability",sep=""))
  }
  colnames(NBC.summaries) <- causes.NBC
  
  # calculate NBC CSMFs by summing 'fractional deaths', i.e. cause-probabilities for each cause
  csmf.fracDeath <- colSums(NBC.summaries)/nrow(NBC.summaries)
  # make csmf.fracDeath a matrix
  csmf.fracDeath <- matrix(csmf.fracDeath,nrow=1)
  # label columns with cause names
  colnames(csmf.fracDeath) <- rownames(SCI)
  # label rows - subsitute() returns the name of its input
  rownames(csmf.fracDeath) <- c(paste(substitute(deaths),".fracDeath",sep=""))
  
  # calcualte CSMF by summing up classified deaths for each cause
  # initialize vector of CSMFs
  csmf.topCause <- c()
  # loop over causes
  for (c in 1:nrow(SCI)) {
    # calculate and append fraction of deaths assigned to each cause
    csmf.topCause <- c(csmf.topCause,length(which(NBC.causes[,1]==paste("cause",c,sep="")))/nrow(deaths))
  }
  # make csmf.topCause a matrix
  csmf.topCause <- matrix(csmf.topCause,nrow=1)
  # label columns
  colnames(csmf.topCause) <- rownames(SCI)
  # label rows - subsitute() returns the name of its input
  rownames(csmf.topCause) <- c(paste(substitute(deaths),".topCause",sep=""))
  
  # create list of outputs to return
  ret.list = list(
    # deaths x cause,(largest)cause-probability matrix of NBC-identified causes
    top.cause = NBC.causes,
    # deaths x cause-probabilities matrix of all cause-probabilities
    cause.probabilities = NBC.summaries,
    # CSMFs calculated by summing 'fractional deaths' for each cause,  
    #   i.e. cause-probabilities for each cause 
    csmf.fracDeath = csmf.fracDeath,
    # CSMFs calculated by calculating the fraction of deaths assigned to each cause
    csmf.topCause = csmf.topCause,
    # datetime when saved
    when = format(Sys.time(),"%a %b %d %X %Y"),
    # deaths
    deaths = deaths,
    # SCI
    sci = SCI
  )
  
  # return list containing various outputs
  return(ret.list)
  
}


###### <<Function>>
calcInterVA <- function(SCI,SCI.c,deaths) {
  
  # implements the InterVA VA algorithm
  
  # SCI is a cause x symptoms matrix of SCI conditional probabilities p(s|c)
  # SCI.c is p(c) 
  # deaths is a deaths x symptoms matrix of VA deaths containing binary symptom values
  # returns a list containing
  #   * top.cause: deaths x cause,(largest)cause-propensity matrix of InterVA-identified causes
  #   * cause.probabilities: deaths x cause-propensities matrix of all cause-propensities
  #   * csmf.fracDeath: CSMFs calculated by summing 'fractional deaths' for each cause,  
  #   * csmf.topCause: CSMFs calculated by calculating the fraction of deaths assigned to each cause
  #   * when: date and time when finished
  #   * deaths: the input deaths
  #   * sci: the input SCI
  
  # InterVA.props are the InterVA propensities for each death: deaths x causes
  InterVA.props <- array(0,dim=c(nrow(deaths),nrow(SCI)))
  
  # InterVA.causes are causes for maximum InterVA propensities
  InterVA.causes <- data.frame(
    cause=as.character(rep("",nrow(deaths))),
    propensity=rep(0,nrow(deaths)),
    stringsAsFactors = FALSE
  )
  rownames(InterVA.causes) <- rownames(deaths)
  
  # calculate the InterVA propensities for each cause given symptoms vectors, p(c|s)
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # initialize the propensities
    InterVA.props[d,] <- SCI.c
    # loop over symptoms
    for (s in 1:ncol(SCI)) {
      # if the symptom exists
      if (deaths[d,s]==1) {
        # loop over causes
        for (c in 1:nrow(SCI)) {
          # update InterVA.propsensity for the current cause
          InterVA.props[d,c] <- InterVA.props[d,c]*SCI[c,s]
        }
        # normalize the propensities
        InterVA.props[d,] <- InterVA.props[d,]/sum(InterVA.props[d,])
        # # the following is an InterVA4.02 bug that we can turn on by uncommenting these lines
        # # truncate the propensities at 0.00001
        # for (c in 1:nrow(SCI)) {
        #   if (InterVA.props[d,c]<0.00001) {InterVA.props[d,c] <- 0}
        # }
      }
    }
    # for each death identify the cause with the largest propensity greater than 0.4,
    #   or if none, label as 'indeterminate'
    # take the first element from the which() statement because sometimes all values 
    #   of InterVA.props[d,] are equal
    InterVA.causes[d,1] <- ifelse(InterVA.props[d,which(InterVA.props[d,]==max(InterVA.props[d,]))[1]]>=0.4,
                                  paste("cause",which(InterVA.props[d,]==max(InterVA.props[d,]))[1],sep=""),"indeterminate")
    InterVA.causes[d,2] <- ifelse(InterVA.props[d,which(InterVA.props[d,]==max(InterVA.props[d,]))[1]]>=0.4,
                                  InterVA.props[d,which(InterVA.props[d,]==max(InterVA.props[d,]))[1]],0)
  }
  
  # store all InterVA cause-specific propensities for each death in a labeled matrix
  InterVA.summaries <- InterVA.props
  rownames(InterVA.summaries) <- rownames(deaths)
  # initialize vector of column labels, cause-probabilities
  causes.InterVA <- c("cause1.propensity")
  # loop over causes
  for (c in 2:nrow(SCI)) {
    # append each cause-probability to the column labels
    causes.InterVA <- c(causes.InterVA,paste("cause",c,".propensity",sep=""))
  }
  colnames(InterVA.summaries) <- causes.InterVA
  
  # calculate InterVA CSMFs by summing 'fractional deaths', i.e. cause-probabilities for each cause
  csmf.fracDeath <- colSums(InterVA.summaries)/nrow(InterVA.summaries)
  # make csmf.fracDeath a matrix
  csmf.fracDeath <- matrix(csmf.fracDeath,nrow=1)
  # label columns with cause names
  colnames(csmf.fracDeath) <- rownames(SCI)
  # label rows - subsitute() returns the name of its input
  rownames(csmf.fracDeath) <- c(paste(substitute(deaths),".fracDeath",sep=""))
  
  # calcualte CSMF by summing up classified deaths for each cause
  # initialize vector of CSMFs
  csmf.topCause <- c()
  # loop over causes
  for (c in 1:nrow(SCI)) {
    # calculate and append fraction of deaths assigned to each cause
    csmf.topCause <- c(csmf.topCause,length(which(InterVA.causes[,1]==paste("cause",c,sep="")))/nrow(deaths))
  }
  # make csmf.topCause a matrix
  csmf.topCause <- matrix(csmf.topCause,nrow=1)
  # label columns
  colnames(csmf.topCause) <- rownames(SCI)
  # label rows - subsitute() returns the name of its input
  rownames(csmf.topCause) <- c(paste(substitute(deaths),".topCause",sep=""))
  
  # create list of outputs to return
  ret.list = list(
    # deaths x cause,(largest)cause-probability matrix of InterVA-identified causes
    top.cause = InterVA.causes,
    # deaths x cause-probabilities matrix of all cause-probabilities
    cause.propensities = InterVA.summaries,
    # CSMFs calculated by summing 'fractional deaths' for each cause,  
    #   i.e. cause-probabilities for each cause 
    csmf.fracDeath = csmf.fracDeath,
    # CSMFs calculated by calculating the fraction of deaths assigned to each cause
    csmf.topCause = csmf.topCause,
    # datetime when saved
    when = format(Sys.time(),"%a %b %d %X %Y"),
    # deaths
    deaths = deaths,
    # SCI
    sci = SCI
  )
  
  # return list containing various outputs
  return(ret.list)
  
}


###### <<Function>>
calcInSilicoVA <- function(SCI,SCI.c,deaths,samples=500,alpha=0.05) {
  
  # implements the InSilicoVA VA algorithm
  
  # SCI is a cause x symptoms matrix of SCI conditional probabilities p(s|c)
  # SCI.c is p(c) 
  # deaths is a deaths x symptoms matrix of VA deaths containing binary symptom values
  # samples is the number of samples to take from the joint distribution in the estimation step,
  #   default to 500
  # alpha is the Dirichlet alpha parameter, default to 0.05
  # returns a list containing
  #   * top.cause: deaths x cause,(largest)cause-probability matrix of InterVA-identified causes
  #   * cause.probabilities: deaths x cause-probabilities matrix of all cause-probabilities
  #   * csmf.median: CSMFs calculated as the median of the CSMF distributions
  #   * csmf.fivenum: five number summary of the CSMF distributions:
  #   * time: time taken to execute
  #   * when: date and time when finished
  #   * pi: array of individual probabilities of being assigned each cause for 
  #   * Y: array of causes assigned for each sample: sample x causes x deaths
  #   * C: cause-specific mortality fractions for each sample: sample x causes
  #   * deaths: the input deaths
  #   * sci: the input SCI
  
  # record start time
  start.time <- Sys.time()
  
  # InSilicoVA.pi is array of individual probabilities of being assigned each cause for 
  #   each sample: sample x causes x deaths
  InSilicoVA.pi <- array(0,dim=c(samples,nrow(SCI),nrow(deaths))) 
  # InSilicoVA.Y is array of causes assigned for each sample: sample x causes x deaths
  InSilicoVA.Y <- array(0,dim=c(samples,nrow(SCI),nrow(deaths)))  
  # InSilicoVA.C stores cause-specific mortality fractions for each sample: sample x causes
  InSilicoVA.C <- array(0,dim=c(samples,nrow(SCI))) 
  # C is vector of cause-specific mortality fractions in each sample, initialize to SCI.c
  C <- SCI.c
  # Y is array of causes assigned in each sample: causes by deaths
  Y <- array(0,dim=c(nrow(SCI),nrow(deaths))) 
  # pi is an array of probabilities of being assigned each cause in each sample: causes by deaths
  pi <- array(0,dim=c(nrow(SCI),nrow(deaths)))  
  # m is vector of Dirichlet parameters
  m <- c(rep(0,nrow(SCI))) 
  
  # loop over the sample, in each pass draw one value for the cause assignments and the CSMFs
  for (samp in 1:samples) {
    # loop over deaths, initialize pi to first guess of cause-specific mortality fractions
    for (i in 1:nrow(deaths)) {
      pi[,i] <- C
    }
    # for each death, update the probability of being assigned each cause conditional on current C
    # loop over deaths
    for (d in 1:nrow(deaths)) {
      # loop over symptoms
      for (s in 1:ncol(SCI)) {	
        # loop over causes		
        for (c in 1:nrow(SCI)) {
          # naive Bayes classifier used to update probabilities of being assigned each cause
          # multiplies the current fraction of deaths due to that cause by the naive Bayes liklihood 
          # 	of observing the recorded signs/symptoms for each death, conditional on the cause
          pi[c,d] <- pi[c,d] * SCI[c,s]^deaths[d,s] * (1-SCI[c,s])^(1-deaths[d,s])
        }		
      }	
    }
    # draw new causes for each death
    # loop over deaths
    for (d in 1:nrow(deaths)) {
      # draw a cause from a categorical distribution parameterized with the updated probabilities of 
      #	dying from each cause (death-specific, above)
      Y[,d] <- rmultinom(1,1,pi[,d])
      # for each death, record the normalized (sum to 1.0) death-specific probabilities of 
      #   dying from each cause (the area under the full joint distribution adds to 1.0,
      #   not the area under the p(c|s) curve for a specific c and s)
      # loop over causes
      for (c in 1:nrow(SCI)) {
        InSilicoVA.pi[samp,c,d] <- pi[c,d]/sum(pi[,d])
      }
    }
    # record the new cause assignments
    InSilicoVA.Y[samp,,] <- Y
    # create the updated Dirichlet parameter vector
    for (c in 1:nrow(SCI)) {
      # for each cause of death, sum alpha and the total number of deaths just assigned to that cause
      m[c] <- alpha + sum(Y[c,])
    }
    # draw a new set of cause-specific mortality fractions from a Dirichlet using the updated 
    #   parameter vector
    C <- rdirichlet(1,m)	
    # record the updated cause-specific mortality fractions
    InSilicoVA.C[samp,] <- C
  }
  
  # store the median and mean cause-assignment probabilities
  # initialize summaries with one row per death and columns for causes
  InSilicoVA.summaries <- array(0,dim=c(nrow(deaths),(nrow(SCI)*2)))
  rownames(InSilicoVA.summaries) <- rownames(deaths)
  # initilize column names
  sums.colnames <- c()
  # loop over causes
  for (c in 1:nrow(SCI)) {
    # make a column name for median and mean associated with this cause
    tmp1 <- paste("cause",c,".median",sep="")
    tmp2 <- paste("cause",c,".mean",sep="")
    # append the column names 
    sums.colnames <- c(sums.colnames,tmp1,tmp2)
  }
  colnames(InSilicoVA.summaries) <- sums.colnames
  # calculate and store the median and mean cause-assignment probabilities for each death
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # loop over causes
    for (c in 1:nrow(SCI)) {
      # calculate and store median and mean cause-assignment probabilities
      InSilicoVA.summaries[d,(c*2-1)] <- round(summary(InSilicoVA.pi[,c,d])[3],4)
      InSilicoVA.summaries[d,(c*2)] <- round(summary(InSilicoVA.pi[,c,d])[4],4)
    }
  }
  
  # InSilicoVA.causes are causes with largest median cause-assignment probability
  InSilicoVA.causes <- data.frame(
    cause=as.character(rep("",nrow(deaths))),
    median.probability=rep(0,nrow(deaths)),
    stringsAsFactors = FALSE
  )
  rownames(InSilicoVA.causes) <- rownames(deaths)
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # loop over causes to calculate median cause-assignment probabilities
    # initialize medians
    medians <- rep(0,nrow(SCI))
    # loop over causes
    for (c in 1:nrow(SCI)) {
      # calculate and store the median cause-assignment probability for each cause
      medians[c] <- summary(InSilicoVA.pi[,c,d])[3]  
    }
    # for each death identify the cause with the largest probability and store both
    #   the cause label and the probability
    # take the first element from the which() statement because sometimes all values 
    #   of medians are equal
    InSilicoVA.causes[d,1] <- paste("cause",which(medians==max(medians))[1],sep="")
    InSilicoVA.causes[d,2] <- medians[which(medians==max(medians))[1]]
  }
  
  # define the CSMFs using median of the CSMF distributions
  # matrix of cause-specific medians, transpose of apply output to create 1x5 matrix
  csmf.median <- t(matrix(apply(InSilicoVA.C,2,median,nrow=1)))
  # make column names
  # initialize cause labels vector
  causes <- "cause1"
  # loop over causes
  for (c in 2:nrow(SCI)) {
    # append each cause to the column labels
    causes <- c(causes,paste("cause",c,sep=""))
  }
  colnames(csmf.median) <- causes
  # label rows - subsitute() returns the name of its input
  rownames(csmf.median) <- c(paste(substitute(deaths),".median",sep=""))
  
  # calculate the five number summary of the CSMF distributions
  csmf.fivenum <- apply(InSilicoVA.C,2,fivenum)
  colnames(csmf.fivenum) <- causes
  rownames(csmf.fivenum) <- c("min","Q1","median","Q3","max")
  
  # record end time
  end.time <- Sys.time()
  
  # create list of outputs to return
  ret.list = list(
    # deaths x cause,(largest)cause-probability matrix of InterVA-identified causes
    top.cause = InSilicoVA.causes,
    # deaths x cause-probabilities matrix of all cause-probabilities
    cause.probabilities = InSilicoVA.summaries,
    # CSMFs calculated as the median of the CSMF distributions
    csmf.median = csmf.median,
    # five number summary of the CSMF distributions:
    #   smallest, 1st quartile, median, 3rd quartile, largest
    csmf.fivenum = csmf.fivenum,
    # time taken to execute
    time = end.time - start.time,
    # datetime when saved
    when = format(Sys.time(),"%a %b %d %X %Y"),
    # pi is array of individual probabilities of being assigned each cause for 
    #   each sample: sample x causes x deaths
    pi = InSilicoVA.pi,
    # Y is array of causes assigned for each sample: sample x causes x deaths
    Y = InSilicoVA.Y,
    # C stores cause-specific mortality fractions for each sample: sample x causes
    C = InSilicoVA.C,
    # deaths
    deaths = deaths,
    # SCI
    sci = SCI
  )
  
  # return list containing various outputs
  return(ret.list)
  
}


###### <<Function>>
plotInSilicoVA <- function(results,path,conv=0,pi=0,csmfH=0,csmfB=0) {
  
  # plotting functions for InSilicoVA results
  
  # results is an InSilicoVA results list
  # path points to where you want PDF files save, relative to the working directory
  #   must be in form './figs/' including last slash
  # if conv=1, create convergence plots
  # if pi=1, create histograms of pi
  # if csmfH=1, create CSMF histograms
  # if csmfB=1, create CSMF boxplots
  
  # create a list of cause names
  causes <- c("Cause 1")
  for (c in 2:nrow(results$sci)) {
    causes <- c(causes,paste("Cause",c))
  }
  # create a set of colors, one for each cause
  colors.cause <- rainbow(nrow(results$sci),alpha=0.5)
  
  # grab the deaths from the results
  deaths <- results$deaths
  
  # assess sample convergence: plot the probability of being assigned
  #   to each cause by sample
  if (conv==1) {
    # loop over deaths
    for (d in 1:nrow(deaths)) {
      # loop over causes
      for (c in 1:nrow(results$sci)) {
        title <- paste("death ",d,","," cause ",c,sep="")
        file.name <- paste(path,"convergence-death",d,",","cause",c,".pdf",sep="")
        pdf(file=file.name)
        # plot pi by sample
        plot(results$pi[,c,d],ylim=c(0,1),main=title,xlab="Sample",ylab="Probability")
        dev.off()
      }
    }
  }
  
  # a single histogram for each death with all cause-assignment probability distributions
  if (pi==1) {
    # loop over deaths
    for (d in 1:nrow(deaths)) {
      pdf(paste(path,"cause-distributions-death",d,".pdf",sep=""))
      hist(results$pi[,1,d],col=colors.cause[1],xlim=c(0,1),main=paste("Death",d),xlab="Probability")
      for (c in 2:nrow(results$sci)) {
        hist(results$pi[,c,d],col=colors.cause[c],add=T)
      }
      legend("topright",causes,fill=colors.cause)
      dev.off()
    }
  }
  
  # plot CSMFs on one plot
  if (csmfH==1) {
    pdf(paste(path,"CSMFs.pdf",sep=""))
    hist(results$C[,1],col=colors.cause[1],main="CSMFs",xlab="Probability",xlim=c(0,1))
    for (c in 2:nrow(results$sci)) {
      hist(results$C[,c],col=colors.cause[c],add=T)
    }
    legend("topright",causes,fill=colors.cause)
    dev.off()
  }
  
  # boxplots of the CSMFs
  if (csmfB==1) {
    pdf(paste(path,"CSMF-boxplots.pdf",sep=""))
    boxplot(results$C,names=causes,main="Summary of CSMF Distributions")
    dev.off()
  }
  
}


###### <<Function>>
calcTariff <- function(GSdeaths,deaths,normalize,resample) {
  
  # implements the Tariff VA algorithm
  
  # GSdeaths is deaths by cause,sypmtoms matrix of gold standard deaths containing 
  #   reference cause and binary symptoms values
  # deaths is a deaths x symptoms matrix of VA deaths containing binary symptom values
  # normalize is boolean indicating whether or not to normalize the tariff sums
  # resample is a boolean indicating whether or not to resample the gold standard deaths 
  # returns a list containing
  #   * top.cause: deaths x cause,tariff-Q matrix of Tariff-identified causes
  #   * cause.tariffScores: deaths x cause-Q, tariff rank, matrix with all cause-Qs
  #   * csmf: CSMFs calculated from identified causes
  #   * when: datetime when saved
  #   * deaths: input deaths 
  #   * tariffs: tariff values calculated from gold standard deaths, i.e. SCI for Tariff
  
  # calculate the Tariffs from the gold standard deaths
  
  # identify the causes in the gold standard deaths
  causes <- sort(unique(GSdeaths[,1]))
  # store the number of causes
  num.causes <- length(causes)
  # store the number of symptoms
  num.symptoms <- ncol(GSdeaths)-1
  # create a vector of cause labels
  cause.labels <- "cause1"
  for (i in 2:num.causes) {
    # create and append the cause label
    cause.labels <- c(cause.labels,paste("cause",i,sep=""))
  }
  # create vector of symptoms labels
  symptom.labels <- "s1"
  for (s in 2:num.symptoms) {
    # create and append the symptom label
    symptom.labels <- c(symptom.labels,paste("s",s,sep=""))
  }
  
  # for each cause, calculate the fraction of deaths with each symptom, the tariff x
  # initilize the list of (raw) tariff x vectors
  x <- list()
  # loop over causes
  for (c in 1:num.causes) {
    # calculate the fraction of gold standard deaths that endorse each symptom, a by-symptom vector
    # append that vector to the list
    if (normalize) {
      x[[c]] <- colSums(GSdeaths[which(GSdeaths[,1]==c),2:6])/nrow(SCIdeaths[which(SCIdeaths[,1]==c),])
    } else {
      x[[c]] <- colSums(GSdeaths[which(GSdeaths[,1]==c),2:6])
    }
  }
  # create a cause x symptom matrix from the raw tariff list
  x.mat <- matrix(unlist(x), ncol = 5, byrow = TRUE)
  
  # create the tariff values by standardizing the (raw) tariff x values using the median 
  #   and interquartile range *across causes*
  # create an empty T matrix
  T <- x.mat-x.mat
  # calculate the medians of x across causes, i.e. by column
  medians <- apply(x.mat,2,median)
  # calculate the interquartile ranges of x across causes
  iqrs <- apply(x.mat,2,IQR)
  # loop over causes
  for (c in 1:num.causes) {
    # calculate and store the median-IQR standardized tariff values
    T[c,] <- (x.mat[c,] - medians) / iqrs
  }
  colnames(T) <- symptom.labels
  rownames(T) <- cause.labels
  
  # create the gold standard tariff scores for standardizing Q values
  # resample the gold standard deaths with replacement to create a uniform distribution by cause
  pi.bootstrap <- rdirichlet(1, rep(1, length(causes)))
  # n.bootstrap <- round(pi.bootstrap * dim(GSdeaths)[1])
  n.bootstrap <- table(sample(causes,nrow(GSdeaths),replace=TRUE,prob=pi.bootstrap))
  # ensure that each cause is sampled at least once
  for (c in 1:num.causes) {
    if (is.na(n.bootstrap[c])) {
      n.bootstrap[c] <- 1
    }
  }
  # create array for bootstrapped gold standard deaths
  GSdeaths.bootstrap <-  array(0,dim=c(sum(n.bootstrap), ncol(GSdeaths)))
  
  # sample the gold standard deaths with replacement by cause with the number for each cause
  #   determined by the dirichlet sample drawn just above
  # initialize the index for the array of bootstrapped gold standard deaths to 0, 
  #   none stored yet
  nsum <- 0
  # loop over causes
  for(c in 1:length(causes)){
    # store the number of deaths to sample for this cause
    nn <- n.bootstrap[c]
    # store the sampled deaths for this cause in the bootstrapped gold standard deaths array
    GSdeaths.bootstrap[(nsum+1):(nsum+nn),] <- 
      GSdeaths[sample(which(GSdeaths[,1]==c),nn,replace=TRUE),]
    # increment the index for the array of bootstrapped gold standard deaths by the number of 
    #   deaths just sampled for this cause
    nsum <- nsum + nn
  }
  # ensure that the bootstrapped data set has at most the original number of deaths
  GSdeaths.bootstrap <- GSdeaths.bootstrap[1:nrow(GSdeaths),]
  
  # pick the gold standard deaths to use when calculating the gold standard tariff scores, just below
  #   if 'resample' is TRUE, then use the resampled (bootstrapped) gold standard deaths
  if (resample) {
    GSdeaths.S <- GSdeaths.bootstrap
  } else {
    GSdeaths.S <- GSdeaths 
  }
  
  # create the gold standard tariff scores for each cause and each death for 
  #   standardizing Q values, later
  # initilize a deaths x causes array of S values
  S.gs <-  array(0,dim=c(nrow(GSdeaths.S),num.causes))  
  # loop over gold standard deaths
  for (d in 1:nrow(GSdeaths.S)) {
    # loop over causes
    for (c in 1:num.causes) {
      # loop over symptoms
      for (s in 1:num.symptoms) {
        # calculate and store the tariff score for this death and cause
        S.gs[d,c] <- S.gs[d,c] + T[c,s] * GSdeaths.S[d,(s+1)] 
      }
    }
  }
  
  # implement cause-assignment using tariffs
  
  # calculate tariff scores for each cause for each death in the data (not gold standard)
  # initialize a deaths x cause arrage of S
  S <-  array(0,dim=c(nrow(deaths),num.causes))  
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # loop over causes
    for (c in 1:num.causes) {
      # loop over symptoms
      for (s in 1:num.symptoms) {
        # calcualte and store the tariff score for this death and cause
        S[d,c] <- S[d,c] + T[c,s] * deaths[d,s] 
      }
    }
  }
  
  # calculate tariff score rank in the gold standard distribution of tariff scores by cause
  # create an empty Q matrix
  Q <- S-S
  # loop over deaths
  for (d in 1:nrow(deaths)) {
    # loop over causes
    for (c in 1:num.causes) {
      # rank of the tariff score for this cause in the distribution of gold standard tariff 
      #   scores for this cause add the tariff score for this death to the gold standard 
      #   distribuion and calculate a normalized rank score
      Q[d,c] <- tail(rank(c(S.gs[which(GSdeaths[,1]==c),c],S[d,c]))
                     /length(c(S.gs[which(GSdeaths[,1]==c),c],S[d,c])),n=1)
    }
  }
  
  # for each death identify causes with largest tariff scores
  # create an empty data frame to hold results
  Tariff.causes <- data.frame(
    cause=as.character(rep("",nrow(deaths))),
    Q=rep(0,nrow(deaths)),
    stringsAsFactors = FALSE
  )
  rownames(Tariff.causes) <- rownames(deaths)
  # loop over deaths 
  for (d in 1:nrow(deaths)) {
    # identify and store the cause and Q value for the cause with the largest Q for this death
    Tariff.causes[d,1] <- paste("cause",which(Q[d,]==max(Q[d,]))[1],sep="")
    Tariff.causes[d,2] <- Q[d,which(Q[d,]==max(Q[d,]))[1]]
  }
  
  # store all Tariff results
  Tariff.summaries <- Q
  rownames(Tariff.summaries) <- rownames(deaths)
  # create a vector of tariff Q labels
  causes.tariff <- c("cause1.Q")
  for (c in 2:num.causes) {
    # create and append each tariff Q label
    causes.tariff <- c(causes.tariff,paste("cause",c,".Q",sep=""))
  }
  colnames(Tariff.summaries) <- causes.tariff
  
  # calcualte and store tariff CSMFs
  # initialize vector of CSMFs
  csmf <- c()
  # loop over causes
  for (c in 1:num.causes) {
    # calculate the fraction of deaths assigned to this cause and append to the vector
    csmf <- c(csmf,length(which(Tariff.causes[,1]==paste("cause",c,sep="")))/nrow(deaths))
  }
  # create a matrix from the vector
  csmf <- matrix(csmf,nrow=1)
  colnames(csmf) <- cause.labels
  # label rows - subsitute() returns the name of its input
  rownames(csmf) <- c(substitute(deaths))
  
  # create list of outputs to return
  ret.list = list(
    # deaths x cause,tariff-Q matrix of Tariff-identified causes
    top.cause = Tariff.causes,
    # deaths x cause-Q, tariff rank, matrix with all cause-Qs
    cause.tariffScores = Tariff.summaries,
    # CSMFs calculated from identified causes
    csmf = csmf,
    # datetime when saved
    when = format(Sys.time(),"%a %b %d %X %Y"),
    # deaths
    deaths = deaths,
    # SCI
    tariffs = T
  )
  
  # return list containing various outputs
  return(ret.list)
  
}


###### <<Function>>
calcConcordance <- function(topCauses,refDeaths) {
  
  # calculate concordance among real cause and assigned causes
  
  # topCauses is a named list of top assigned causes, e.g.
  # topCauses <- list(
  #   NBC =  NBC.results$top.cause,
  #   InterVA = InterVA.results$top.cause,
  #   InSilicoVA = InSilicoVA.results$top.cause,
  #   Tariff = Tariff.results$top.cause
  # )
  # refDeaths is a set of deaths with a reference cause, e.g.
  # refCSMF <- testDeaths.CSMF
  # returns a list containing
  #   * discordant, matrix of discordant cause counts
  #   * percent.correct, matrix of percent count of causes that are the same 
  
  # creat a dataframe with the cause assignments from all the inputs
  # first, strip off the cause-assignment probabilities, etc.
  # loop over the inputs
  for (i in 1:length(topCauses)) {
    # store just the cause label for this input
    tmp <- topCauses[[i]][,1]
    # bring the death names along as labels for the causes
    names(tmp) <- rownames(topCauses[[i]])
    # store the assigned causes back into the 'topCauses' list
    topCauses[[i]] <- tmp
  }
  # second, merge the inputs into a dataframe
  # merge on the row names that identify the deaths, the ensures that deaths are 
  #   matched across inputs
  # merge the first two inputs' assigned causes and store in a dataframe called 'cause.df'
  cause.df <- merge(topCauses[[1]],topCauses[[2]],by="row.names")
  # assign the merged death names to rownames of the dataframe
  rownames(cause.df) <- cause.df$Row.names
  # remove the merged death names column from the dataframe
  cause.df <- cause.df[,2:ncol(cause.df)]
  # name the columns of the dataframe using the names from teh 'topCauses' list
  colnames(cause.df) <- c(names(topCauses)[1:2])
  # if there are more than two inputs, repeat the previous few steps for the remainder of the inputs
  if (length(topCauses)>2) {
    for (i in 3:length(topCauses)) {
      column.names <- c(colnames(cause.df),names(topCauses)[i])
      cause.df <- merge(cause.df,topCauses[[i]],by="row.names")
      rownames(cause.df) <- cause.df$Row.names
      cause.df <- cause.df[,2:ncol(cause.df)]
      colnames(cause.df) <- column.names
    }
  }
  
  # identify all the cause names
  # initialize a vector of all the cause names using the names from the first input
  cause.names <- as.character(cause.df[,1])
  # loop over the inputs, columns in 'cause.df'
  for (c in 2:ncol(cause.df)) {
    # append the cause names for this input
    cause.names <- c(cause.names,as.character(cause.df[,c]))
  }
  # find the set of unique cause names
  cause.names <- unique(cause.names)
  # store number of unique cause names
  causes <- length(cause.names)
  # store the levels of the cause names 
  cause.levels <- levels(as.factor(cause.names))
  # recode test death causes from numeric labels to character, i.e. 'cause1' or 'cause2'
  # initialize new cause labels vector
  refDeaths.causes <- rep("",nrow(refDeaths))
  # create recode string
  recode.str <- "1='cause1'"
  # loop over causes
  for (c in 2:causes) {
    recode.str <- paste(recode.str,";",c,"='cause",c,"'",sep="")
  }
  # recode!
  suppressWarnings(refDeaths.causes <- recode(refDeaths[,1],recode.str))
  
  # add real causes to Results.causes
  # store the existing column names
  column.names <- colnames(cause.df)
  # merge the reference death causes 
  cause.df <- merge(refDeaths.causes,cause.df,by="row.names")
  # name the deaths using the merged death names
  rownames(cause.df) <- cause.df[,1]
  # remove the merged death names column
  cause.df <- cause.df[,2:ncol(cause.df)]
  # name the columns
  colnames(cause.df) <- c("Reference",column.names)
  # ensure that all columns of 'cause.df' have the same levels
  for (c in 1:ncol(cause.df)) {
    levels(cause.df[,c]) <- cause.levels
  }
  
  # count discordant causes for all combinations of inputs and refernece cause
  # number of inputs and reference
  comparitors <- ncol(cause.df)
  # initialize array of discordant counts
  discordant <- array(0,dim=c(comparitors,comparitors))  
  # loop over comparitors
  for (ref in 1:comparitors) {
    # loop over comparitors
    for (comp in 1:comparitors) {
      # count dicordant causes for current combination of comparitors
      discordant[ref,comp] <- length(which(cause.df[,ref]!=cause.df[,comp]))
    }
  }
  # name the rows and columns of discordant with the names of the comparitors
  rownames(discordant) <- colnames(cause.df)
  colnames(discordant) <- colnames(cause.df)
  
  # create the return list
  ret <- list(
    # the discordant matrix
    discordant = discordant,
    # the percent correct matrix
    percent.correct = 100-100*discordant/nrow(VAdeaths)
  )
  
  return(ret)
  
}


###### <<Function>>
calcCSMFaccuracies <- function(CSMFs,refCSMF) {
  
  # calculate CSMF accuracies, etc.
  
  # CSMFs is a named list of variable number algorithm result objects, e.g.
  # CSMFs <- list(
  #   NBC.frac = NBC.results$csmf.fracDeath,
  #   NBC.top = NBC.results$csmf.topCause,
  #   InterVA.frac = InterVA.results$csmf.fracDeath,
  #   InterVA.top = InterVA.results$csmf.topCause,
  #   InSilicoVA = InSilicoVA.results$csmf.median,
  #   Tariff = Tariff.results$csmf
  # )
  # refCSMF is a reference set of CSMFs produced by calcCSMF(), e.g.
  # refCSMF <- testDeaths.CSMF
  # returns a list containing
  #   * CSMF.errors, raw CSMF errors compared to the reference
  #   * CSMF.tae, total absolute errors
  #   * CSMF.percErr, percent proportional errors
  #   * CSMF.percErr.tae, total absolute percent proportional errors
  #   * CSMF.accuracies, CSMF Accuracies
  
  # wrap up CSMFs nicely
  CSMF.mat <- rbind(refCSMF,CSMFs[[1]])
  for (c in 2:length(CSMFs)) {
    CSMF.mat <- rbind(CSMF.mat,CSMFs[[c]])
  }
  rownames(CSMF.mat) <- c("reference",names(CSMFs))
  CSMF.mat
  
  # Error in CSMF, subtract first (reference) row from all rows
  CSMF.errors <- sweep(CSMF.mat,2,CSMF.mat[1,])
  # total absolute error in CSMFs
  CSMF.tae <- rowSums(abs(CSMF.errors))
  
  # proportional error
  CSMF.percErr <- 100*sweep(CSMF.errors,2,CSMF.mat[1,],FUN="/")
  # total absolute proportional error
  CSMF.percErr.tae <- rowSums(abs(CSMF.percErr))
  
  # CSMF accuracy
  CSMF.accuracies <- c()
  for (i in 1:nrow(CSMF.mat)) {
    CSMF.accuracies <- 
      c(CSMF.accuracies, 1 - sum(abs(CSMF.mat[i,] - CSMF.mat[1,])) / (2 * (1 - min(CSMF.mat[1,]))))
  }
  names(CSMF.accuracies) <- rownames(CSMF.mat)
  CSMF.accuracies
  
  ret <- list(
    CSMF.errors = CSMF.errors,
    CSMF.tae = CSMF.tae,
    CSMF.percErr = CSMF.percErr,
    CSMF.percErr.tae = CSMF.percErr.tae,
    CSMF.accuracies = CSMF.accuracies
  )
  
  return(ret)
  
}


###### <<Function>>
crossValidate <- function(SCI.denovo,csmf,sim.deaths,samp.frac,insilico.reps=500,insilico.alpha=0.05,mods=NULL) {
  
  # run a cross validation set
  #   * simulate deaths
  #   * draw a sample of the simulated deaths
  #   * use the sample to create SCI
  #   * run all the algorithms using the SCI from the sample to assign
  #     casues to the deaths that were not sampled, the 'test deaths'
  #   * compare the assigned causes from each algorithm to the real 
  #     causes in the test deaths
  
  # SCI.denovo is a list of <symptom>-element SCI vectors, one for each cause
  # csmf is a <cause>-element vector of CSMFs
  # sim.deaths is number of deaths to simulate
  # samp.frac is a decimal number representing the fraction of deaths to sample
  # insilico.reps is the number of sampling repetitions for InSilicoVA
  # insilico.alpha is InSilicoVA's alpha parameter
  
  # returns a list of
  #   * simDeaths, a matrix of the simulated deaths 
  #   * result objects for each algorithm, NBC.results, etc.
  #   * concordance, a concordance object comparing the results
  #   * CSMF.accuracies, a CSMF accuracies object comparing CSMF results
  
  # store SCI defined above
  SCI <- storeSCI(SCI.denovo)
  # simulate a bunch of deaths using the SCI and csmf defined above
  simDeaths <- simulateDeaths(sim.deaths,csmf,SCI)
  # draw sample for to prepare SCI
  sample <- sample(nrow(simDeaths),samp.frac*nrow(simDeaths))
  # select the sampled deaths to be used to create SCI
  SCIdeaths <- simDeaths[sample,]
  # select the inverse of the sample as deaths for testing
  testDeaths <- simDeaths[(-sample),]
  # calcualte teh CSMF of the SCI and test deaths
  testDeaths.CSMF <- calcCSMF(testDeaths)
  # calcualte SCI from the SCI deaths
  valSCI <- calcSCI(SCIdeaths)
  
  # if 'mods', create modified SCI
  if (!is.null(mods)) {
    # scale the whole SCI matrix
    SCI.less <- expit(logit(valSCI) * mods[1])
    # initialize SCI.modified with SCI.less
    SCI.modified <- SCI.less
    
    # reduce the influence of symptom 2 on logit scale
    SCI.modified[,2] <- expit(logit(SCI.less[,2]) + mods[2])
    # increase the influence of symptom 3 by 20% on logit scale
    SCI.modified[,3] <- expit(logit(SCI.less[,3]) + mods[3])
    # increase the influence of symptom 5 by 75% on logit scale
    SCI.modified[,5] <- expit(logit(SCI.less[,5]) + mods[4])
    # make all symptoms less likely for cause 3
    SCI.modified[3,] <- expit(logit(SCI.less[3,]) + mods[5])
    
    valSCI <- SCI.modified
    
    # for Tariff
    valSCIdeaths <- simulateDeaths(length(sample),csmf,valSCI)
    
  }
  
  
  # create a starting p(c) for all the methods, a uniform distribution
  prob.c <- rep(1/nrow(SCIdeaths.SCI),nrow(SCIdeaths.SCI)) 
  VAdeaths <- testDeaths[,2:6]
  # run NBC
  NBC.results <- calcNBC(valSCI,prob.c,VAdeaths)
  # run InterVA
  InterVA.results <- calcInterVA(valSCI,prob.c,VAdeaths)
  # run InSilicoVA - it take a lot longer!
  InSilicoVA.results <- calcInSilicoVA(valSCI,prob.c,VAdeaths,insilico.reps,insilico.alpha)
  # run Tariff
  Tariff.results <- calcTariff(SCIdeaths,VAdeaths,TRUE,TRUE)
  
  # concordance
  topCauses <- list(
    NBC =  NBC.results$top.cause,
    InterVA = InterVA.results$top.cause,
    InSilicoVA = InSilicoVA.results$top.cause,
    Tariff = Tariff.results$top.cause
  )
  # calculate and have a look at the concordance matrices
  concordances <- calcConcordance(topCauses,testDeaths)
  
  # CSMF accuracies
  # create named list of CSMFs; names will be used throughout, so must be meaningful
  CSMFs <- list(
    NBC.frac = NBC.results$csmf.fracDeath,
    NBC.top = NBC.results$csmf.topCause,
    InterVA.frac = InterVA.results$csmf.fracDeath,
    InterVA.top = InterVA.results$csmf.topCause,
    InSilicoVA = InSilicoVA.results$csmf.median,
    Tariff = Tariff.results$csmf
  )
  # calculate CSMF accuracies
  CSMF.accuracies <- calcCSMFaccuracies(CSMFs,testDeaths.CSMF)
  
  # make a return list
  rep <- list(
    simDeaths = simDeaths,
    NBC.results = NBC.results,
    InterVA.results = InterVA.results,
    InSilicoVA.results = InSilicoVA.results,
    Tariff.results = Tariff.results,
    concordance.list = concordances,
    CSMFaccuracies.list = CSMF.accuracies
  )
  
  return(rep)
  
}


###### <<Function>>
logit <- function(x) {log(x/(1-x))}


###### <<Function>>
expit <- function(x) {exp(x)/(1+exp(x))}


