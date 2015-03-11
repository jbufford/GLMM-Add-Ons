########################## Make AIC Table for GLMMs ###########################
                            ## Jennifer Bufford ##
                          ## jlbufford@gmail.com ##
                           ## February 16, 2015 ##

###############################################################################

#Input: A GLMM, optionally a previous model to compare it to or a table to
#         add the results to

#Requires: AICcmodavg and either lme4, nlme or glmmADMB

#Output: Table of AIC(c) values for comparing a suite of a priori models

#make.tab makes a table of loglik, AIC(c), delta AIC(c), p-values to compare models
#May not work properly with some correlation structures in nlme (e.g. VarPower)


##### Make AIC Table Function #################################################


make.tab <- function(new.mod, prev.mod=NULL, m.tab=NA){


  ##### Calculate K #####

  K <- calc.AICc(new.mod, calc.K=T)

  data.n <- calc.AICc(new.mod, calc.n=T)


  ##### Calculate P-Value #####

  if(!is.null(prev.mod)) {

    if(class(new.mod) %in% c("lmerMod", "glmerMod", "glmmadmb")) {
      pv <- anova(new.mod, prev.mod)$Pr[2]
    }

    if(class(new.mod)=="lme") {pv <- anova(new.mod, prev.mod)$p[2]}

  } else {pv <- NA}


  ##### Create Table #####

  if(((data.n/K) < 40) | "AICc" %in% names(m.tab)) {

    mtab <- data.frame("Model" = as.character(formula(new.mod))[3],
                       "LogLik" = logLik(new.mod), "AICc" = calc.AICc(new.mod),
                       "dAICc"=NA, "PValue" = pv)
    print("Use AICc")

  } else {

    mtab <- data.frame("Model" = as.character(formula(new.mod))[3],
                       "LogLik" = logLik(new.mod), "AIC" = AIC(new.mod),
                       "dAIC"=NA, "PValue" = pv)
    print("Use AIC")
  }

  mtab$Model <- as.character(mtab$Model)


  ##### Merge Tables, Calc dAIC(c) #####

  if(length(m.tab)>1){

    mtab <- rbind(m.tab, mtab)

    mtab[,4] <- mtab[,3] - min(mtab[,3])

  }

  return(mtab)

}





##### Calc AICc ###############################################################


calc.AICc <- function(mod, calc.K=F, calc.n=F) {

  ###### Calculate K ######

  if(class(mod)=="lmerMod") {

    library(lme4)
    library(AICcmodavg)

    K <- length(fixef(mod)) + length(ranef(mod)) + 1 #for residual rand var

    data.n <- nrow(mod@frame)

    print(paste("K =", K))
    print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
  }

  if(class(mod)=="glmerMod") {

    library(lme4)
    library(AICcmodavg)

    K <- length(fixef(mod)) + length(ranef(mod)) #no resid var est (go figure)

    data.n <- nrow(mod@frame)

    if(calc.K){
      print(paste("K =", K))
      print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
    }
  }

  if(class(mod)=="lme") {

    library(nlme)
    library(AICcmodavg)

    K <- (length(fixef(mod)) + length(ranef(mod)) +
            length(attributes(mod$modelStruct$varStruct)$groupNames) +
            #             length(attributes(mod$modelStruct$varStruct)$groupNames) +
            length(coef(mod$modelStruct$corStruct)))
    #not sure how robust correl structure is
    # Still not working with varPower

    data.n <- nrow(mod$data)

    if(calc.K){
      print(paste("K =", K))
      print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
    }
  }

  if(class(mod)=='glmmadmb') {

    library(glmmADMB)

    K <- (length(fixef(mod)) + length(ranef(mod)) + length(mod$pz) +
            length(mod$alpha))

    data.n <- nrow(mod$frame)
  }

  if(calc.K){return(K)}
  if(calc.n){return(data.n)}


  ## Calc AICc ##

  return(AIC(mod) + ((2*K*(K+1))/(data.n-K-1)))

}





##### Test Candidate Models ###################################################


test.mods <- function(response, code.mods,tab.mods=NULL,rand=NULL, full.mod, dat=NULL){

  cand.list <- list()

  if(!is.null(tab.mods)){

    mod.tab <- data.frame('Model'=tab.mods[1],
                          'Code.Model'=as.character(formula(full.mod))[3],
                          'logLik'=logLik(full.mod), 'AICc'=calc.AICc(full.mod))

    mod.tab$Model <- as.character(mod.tab$Model)
    mod.tab$Code.Model <- as.character(mod.tab$Code.Model)

    cand.list[[1]] <- full.mod
    print(as.character(formula(full.mod))[3])


    for (i in 2:length(code.mods)){

      if(class(full.mod)=="lme"){

        cand.char <- paste(response, "~", code.mods[i])
        print(cand.char)

        cand.mod <- update(full.mod, fixed=cand.char)

      } else {

        cand.char <- paste(response, "~", code.mods[i], "+", rand)
        cand.form <- formula(cand.char)
        print(cand.char)

        if(!is.null(dat)){
          cand.mod <- update(full.mod, formula=cand.form, data=dat) }
        else {
          cand.mod <- update(full.mod, formula=cand.form)
        }

      }

      cand.list[[i]] <- cand.mod

      mod.tab[i,] <- c(as.character(tab.mods[i]), cand.char,
                       logLik(cand.mod), calc.AICc(cand.mod))

    }

    mod.tab$logLik <- as.numeric(mod.tab$logLik)
    mod.tab$AICc <- as.numeric(mod.tab$AICc)
    mod.tab$dAICc <- mod.tab$AICc - min(mod.tab$AICc)

  } else {

    mod.tab <- data.frame('Model'=code.mods[1], 'logLik'=logLik(full.mod),
                          'AICc'=calc.AICc(full.mod))

    mod.tab$Model <- as.character(mod.tab$Model)

    cand.list[[1]] <- full.mod
    print(code.mods[1])


    for (i in 2:length(code.mods)){

      if(class(full.mod)=="lme"){

        cand.char <- paste(response, "~", code.mods[i])
        print(cand.char)

        cand.mod <- update(full.mod, fixed=cand.char)

      } else {

        cand.char <- paste(response, "~", code.mods[i], "+", rand)
        cand.form <- formula(cand.char)
        print(cand.char)

        cand.mod <- update(full.mod, formula=cand.form)

      }


      cand.list[[i]] <- cand.mod

      mod.tab[i,] <- c(cand.char, logLik(cand.mod), calc.AICc(cand.mod))

    }

    mod.tab$logLik <- as.numeric(mod.tab$logLik)
    mod.tab$AICc <- as.numeric(mod.tab$AICc)
    mod.tab$dAICc <- mod.tab$AICc - min(mod.tab$AICc)
  }

  return(list(cand.list, mod.tab))

}
###############################################################################


##### Make AIC Table Function #################################################


make.tab <- function(new.mod, prev.mod=NULL, m.tab=NA){


  ##### Calculate K #####

  K <- calc.AICc(new.mod, calc.K=T)

  data.n <- calc.AICc(new.mod, calc.n=T)


  ##### Calculate P-Value #####

  if(!is.null(prev.mod)) {

    if(class(new.mod) %in% c("lmerMod", "glmerMod", "glmmadmb")) {
      pv <- anova(new.mod, prev.mod)$Pr[2]
    }

    if(class(new.mod)=="lme") {pv <- anova(new.mod, prev.mod)$p[2]}

  } else {pv <- NA}


  ##### Create Table #####

  if(((data.n/K) < 40) | "AICc" %in% names(m.tab)) {

    mtab <- data.frame("Model" = as.character(formula(new.mod))[3],
                       "LogLik" = logLik(new.mod), "AICc" = calc.AICc(new.mod),
                       "dAICc"=NA, "PValue" = pv)
    print("Use AICc")

  } else {

    mtab <- data.frame("Model" = as.character(formula(new.mod))[3],
                       "LogLik" = logLik(new.mod), "AIC" = AIC(new.mod),
                       "dAIC"=NA, "PValue" = pv)
    print("Use AIC")
  }

  mtab$Model <- as.character(mtab$Model)


  ##### Merge Tables, Calc dAIC(c) #####

  if(length(m.tab)>1){

    mtab <- rbind(m.tab, mtab)

    mtab[,4] <- mtab[,3] - min(mtab[,3])

  }

  return(mtab)

}





##### Calc AICc ###############################################################


calc.AICc <- function(mod, calc.K=F, calc.n=F) {

  ###### Calculate K ######

  if(class(mod)=="lmerMod") {

    library(lme4)
    library(AICcmodavg)

    K <- length(fixef(mod)) + length(ranef(mod)) + 1 #for residual rand var

    data.n <- nrow(mod@frame)

    print(paste("K =", K))
    print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
  }

  if(class(mod)=="glmerMod") {

    library(lme4)
    library(AICcmodavg)

    K <- length(fixef(mod)) + length(ranef(mod)) #no resid var est (go figure)

    data.n <- nrow(mod@frame)

    if(calc.K){
      print(paste("K =", K))
      print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
    }
  }

  if(class(mod)=="lme") {

    library(nlme)
    library(AICcmodavg)

    K <- (length(fixef(mod)) + length(ranef(mod)) +
            length(attributes(mod$modelStruct$varStruct)$groupNames) +
            #             length(attributes(mod$modelStruct$varStruct)$groupNames) +
            length(coef(mod$modelStruct$corStruct)))
    #not sure how robust correl structure is
    # Still not working with varPower

    data.n <- nrow(mod$data)

    if(calc.K){
      print(paste("K =", K))
      print(paste("AICcmodavg K =", AICc(mod, return.K=T)))
    }
  }

  if(class(mod)=='glmmadmb') {

    library(glmmADMB)

    K <- (length(fixef(mod)) + length(ranef(mod)) + length(mod$pz) +
            length(mod$alpha))

    data.n <- nrow(mod$frame)
  }

  if(calc.K){return(K)}
  if(calc.n){return(data.n)}


  ## Calc AICc ##

  return(AIC(mod) + ((2*K*(K+1))/(data.n-K-1)))

}





##### Test Candidate Models ###################################################


test.mods <- function(response, code.mods,tab.mods=NULL,rand=NULL, full.mod, dat=NULL){

  cand.list <- list()

  if(!is.null(tab.mods)){

    mod.tab <- data.frame('Model'=tab.mods[1],
                          'Code.Model'=as.character(formula(full.mod))[3],
                          'logLik'=logLik(full.mod), 'AICc'=calc.AICc(full.mod))

    mod.tab$Model <- as.character(mod.tab$Model)
    mod.tab$Code.Model <- as.character(mod.tab$Code.Model)

    cand.list[[1]] <- full.mod
    print(as.character(formula(full.mod))[3])


    for (i in 2:length(code.mods)){

      if(class(full.mod)=="lme"){

        cand.char <- paste(response, "~", code.mods[i])
        print(cand.char)

        cand.mod <- update(full.mod, fixed=cand.char)

      } else {

        cand.char <- paste(response, "~", code.mods[i], "+", rand)
        cand.form <- formula(cand.char)
        print(cand.char)

        if(!is.null(dat)){
          cand.mod <- update(full.mod, formula=cand.form, data=dat) }
        else {
          cand.mod <- update(full.mod, formula=cand.form)
        }

      }

      cand.list[[i]] <- cand.mod

      mod.tab[i,] <- c(as.character(tab.mods[i]), cand.char,
                       logLik(cand.mod), calc.AICc(cand.mod))

    }

    mod.tab$logLik <- as.numeric(mod.tab$logLik)
    mod.tab$AICc <- as.numeric(mod.tab$AICc)
    mod.tab$dAICc <- mod.tab$AICc - min(mod.tab$AICc)

  } else {

    mod.tab <- data.frame('Model'=code.mods[1], 'logLik'=logLik(full.mod),
                          'AICc'=calc.AICc(full.mod))

    mod.tab$Model <- as.character(mod.tab$Model)

    cand.list[[1]] <- full.mod
    print(code.mods[1])


    for (i in 2:length(code.mods)){

      if(class(full.mod)=="lme"){

        cand.char <- paste(response, "~", code.mods[i])
        print(cand.char)

        cand.mod <- update(full.mod, fixed=cand.char)

      } else {

        cand.char <- paste(response, "~", code.mods[i], "+", rand)
        cand.form <- formula(cand.char)
        print(cand.char)

        cand.mod <- update(full.mod, formula=cand.form)

      }


      cand.list[[i]] <- cand.mod

      mod.tab[i,] <- c(cand.char, logLik(cand.mod), calc.AICc(cand.mod))

    }

    mod.tab$logLik <- as.numeric(mod.tab$logLik)
    mod.tab$AICc <- as.numeric(mod.tab$AICc)
    mod.tab$dAICc <- mod.tab$AICc - min(mod.tab$AICc)
  }

  return(list(cand.list, mod.tab))

}