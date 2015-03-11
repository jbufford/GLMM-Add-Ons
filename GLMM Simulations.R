############################## GLMM Simulations ###############################
                            ## Jennifer Bufford ##
                          ## jlbufford@gmail.com ##
                           ## February 16, 2015 ##

###############################################################################

#Input: An lmer model, the data used to make the model, optional specifiers

#Requires: Either lme4, nlme or glmmADMB

#Output: Graphs of simulated coefficients and actual coefficients
#        Optional: returns simulated data and models as a list

#Takes an existing model and data and simulates data based on the model
#   Then creates models for each simulated data set and calculates coefficients
#   Compares simulated coefficients to either the initial model or to
#     another test model (e.g. a model that does or does not account for a term)

###############################################################################

## Note - this function works well for nlme w/ blocked random effects and het var
##    I haven't fully tested it on other possible models, but it should work for
##    lme4, glmmadmb and nlme models w/ one or more random effects and optionally with
##      an offset term, for nlme models w/ 1 random effect, 2 nested random effects, or
##      blocked random effects, and possibly for gls models as well
##    It does not incorporate weights or non-gaussian models (at the moment)

sim.glmm <- function(M, dat, nsim=100, test.mod=NA, respvar=NA, extra=NA, to.return=F) {

  library(plyr, quietly=T)

  os <- NA
  het.var <- NA


  ###### Extract Model Formula, Terms ######

  if(class(M)=="lmerMod" | class(M)=="glmerMod") {

    library(lme4)

    if(is.na(respvar)) { respvar <- names(M@frame)[1] }
    randvar <- names(M@flist)
    fixvar <- attributes(terms(M@frame))$term.labels
    fixvar <- fixvar[!(fixvar%in% randvar)]
    Mterms <- names(M@frame)[2:length(names(M@frame))]
    if("(weights)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if(sum(grepl("offset", Mterms))>0) {
      os <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}
    if(class(M)=="glmerMod") { fam <- M@resp$family[[1]] } else {fam <- 'gaussian'}

    rf <- ldply(ranef(M), function(x){
      x <- data.frame(x)
      x$cat <- ifelse(grepl("/", row.names(x)),
                      gsub("[[:print:]]*/", '', row.names(x)), row.names(x))
      x$coef <- x[,1]
      return(x[,c('coef','cat')])
    })
    names(rf)[names(rf)==".id"] <- 'var'

    resid.dev <- ifelse(is.na(M@devcomp$cmp['sigmaREML']), M@devcomp$cmp['sigmaML'],
                        M@devcomp$cmp['sigmaREML'])
  }

  if(class(M)=="lme") {

    library(nlme)

    if(is.na(respvar)) { respvar <- as.character(attributes(M$terms)$variables[2]) }

    fixvar <- attributes(M$terms)$term.labels

    #If there are fixed effects
    if(length(attributes(M$terms)$term.labels) > 0) {
      Mterms <- c(as.character(attributes(M$terms)$variables
                               [3:length(attributes(M$terms)$variables)]),
                  attributes(M$modelStruct$reStruct)$names, extra)
      if(sum(grepl("offset", Mterms))>0) {
        os <- gsub('offset\\(|\\)',"",grep("offset", Mterms, value=T))
        Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}

      randvar <- if(is.na(extra[1])){attributes(M$modelStruct$reStruct)$names} else {
        extra}

    } else { #if there are no fixed effects except intercept
      randvar <- if(is.na(extra[1])){ attributes(M$modelStruct$reStruct)$names } else {
        extra}
      Mterms <- randvar
    }

    resid.dev <- M$sigma

    if('varStruct' %in% names(M$modelStruct)){
      het.var <- attributes(M$modelStruct$varStruct)$weights
    }

    if(is.list(ranef(M)) & length(randvar)>1) {
      rf <- data.frame(t(ranef(M)[1,]), 'var'=NA, 'cat'=NA)
      names(rf)[1] <- 'coef'
      for(i in randvar){
        n <- which(substr(names(ranef(M)), 1, nchar(i)) %in% i)
        rf$var[n] <- substr(names(ranef(M)), 1, nchar(i))[n]
        rf$cat[n] <- substr(names(ranef(M)), nchar(i)+1, 20)[n]
      }
    } else {
      if(is.list(ranef(M))){
        rf <- ldply(ranef(M), function(x){
          x$cat <- ifelse(grepl("/", row.names(x)),
                          gsub("[[:print:]]*/", '', row.names(x)), row.names(x))
          x$coef <- x[,1]
          return(x[,c('coef','cat')])
        })
        names(rf)[names(rf)==".id"] <- 'var'
      } else {
        rf <- data.frame('coef'=ranef(M)[,1], 'var'=randvar, 'cat'=row.names(ranef(M)))
      }
    }

    fam <- 'gaussian'
  }

  if(class(M)=='gls') {

    library(nlme)

    fixvar <- names(M$parAssign)[!names(M$parAssign) %in% "(Intercept)"]
    Mterms <- names(M$parAssign)[-grep(names(M$parAssign), pattern=":")][-1]
    randvar <- NA
  }

  if(class(M)=='glmmadmb') {

    library(glmmADMB)

    if(is.na(respvar)) { respvar <- names(M$frame)[1] }
    randvar <- names(M$S)
    fixvar <- attributes(M$terms)$term.labels
    Mterms <- c(names(M$frame)[2:length(M$frame)], randvar)
    #     if("(weights)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if(sum(grepl("offset", Mterms))>0) {
      os <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}
    dat$Resid <- as.vector(M$residuals) #pearson's residuals
    #(not sure why but resid() isn't working)

    rf <- ldply(ranef(M), function(x){
      x <- data.frame(x)
      x$cat <- ifelse(grepl("/", row.names(x)),
                      gsub("[[:print:]]*/", '', row.names(x)), row.names(x))
      x$coef <- x[,1]
      return(x[,c('coef','cat')])
    })
    names(rf)[names(rf)==".id"] <- 'var'

    resid.dev <- ifelse(is.na(M@devcomp$cmp['sigmaREML']), M@devcomp$cmp['sigmaML'],
                        M@devcomp$cmp['sigmaREML'])

    fam <- M$family
  }


  ##### Get Coefficients #####

  cf <- data.frame('Coef'=fixef(M), 'Code'=names(fixef(M)), 'Term'=names(fixef(M)),
                   'Cat'=NA, 'Term2'=NA, 'Cat2'=NA, stringsAsFactors = F)

  for (j in which(!cf$Code %in% c('(Intercept)', Mterms))){
    if(grepl(":", cf$Code[j])){
      cf$Term2[j] <- strsplit(cf$Code[j], split=":")[[1]][2]
    }

    for(k in 1:nchar(cf$Code[j])){
      if(substr(cf$Code[j], 1,k) %in% Mterms){
        cf$Term[j] <- substr(cf$Code[j],1,k)
        cf$Cat[j] <- substr(cf$Code[j],(k+1),nchar(cf$Code[j]))
        if(is.na(cf$Term2[j])|cf$Term2[j] %in% Mterms) {break} else {
          for(l in 1:nchar(cf$Term2[j])){
            if(substr(cf$Term2[j], 1,l) %in% Mterms){
              cf$Cat2[j] <- substr(cf$Term2[j],(l+1),nchar(cf$Term2[j]))
              cf$Term2[j] <- substr(cf$Term2[j],1,l)
              break
            }
          }
        }
      }
    }

    if(grepl(":", cf$Cat[j])){
      cf$Cat[j] <- strsplit(cf$Cat[j], split=":")[[1]][1]
    }
  }


  ##### Simulate Response Var #####

  sim.dat <- list()

  for (i in 1:nsim){

    dsim <- dat

    if(fam=="gaussian"){
      sim.resp <- cf[cf$Term=="(Intercept)", 'Coef'] +
        rnorm(nrow(dsim), 0, resid.dev + ifelse(!is.na(het.var[1]), het.var, 0))
    } else {

      if(fam=="Gamma"|fam=="gamma"){
        sim.resp <- cf[cf$Term=="(Intercept)", 'Coef'] +
          rgamma(nrow(dsim), 1, resid.dev^2)
      }

      if(fam=="binomial"){
        sim.resp <- cf[cf$Term=="(Intercept)", 'Coef'] +
          rbinom(nrow(dsim), 0, resid.dev + ifelse(!is.na(het.var[1]), het.var, 0))
      }

      if(fam=="poisson"){
        sim.resp <- cf[cf$Term=="(Intercept)", 'Coef'] +
          rpoi(nrow(dsim), 0, resid.dev + ifelse(!is.na(het.var[1]), het.var, 0))
      }


    }

    for(j in 2:nrow(cf)){

      if(is.na(cf[j,'Term2'])){
        if(is.na(cf[j,'Cat'])){
          sim.resp <- sim.resp + cf[j, 'Coef']*dat[,cf[j,'Term']]
        } else {
          sim.resp <- sim.resp +cf[j,'Coef']*ifelse(dat[,cf[j,'Term']]==cf[j,'Cat'],1,0)
        }
      } else {
        if(is.na(cf[j,'Cat'])){
          if(is.na(cf[j,'Cat2'])){
            sim.resp <- sim.resp + cf[j, 'Coef']*dat[,cf[j,'Term']]*dat[,cf[j,'Term2']]
          } else {
            sim.resp <- sim.resp + cf[j, 'Coef']*dat[,cf[j,'Term']]*
              ifelse(dat[,cf[j,'Term2']]==cf[j,'Cat2'],1,0)
          }
        } else {
          if(is.na(cf[j,'Cat2'])){
            sim.resp <- sim.resp + cf[j, 'Coef']*
              ifelse(dat[,cf[j,'Term']]==cf[j,'Cat'],1,0)*dat[,cf[j,'Term2']]
          } else {
            sim.resp <- sim.resp + cf[j, 'Coef']*
              ifelse(dat[,cf[j,'Term']]==cf[j,'Cat'],1,0)*
              ifelse(dat[,cf[j,'Term2']]==cf[j,'Cat2'],1,0)
          }
        }
      }
    }

    if(class(M) != 'gls'){
      for(k in unique(rf$var)){
        rfk <- rf[rf$var==k,]
        sim.resp <- sim.resp + rfk[match(dat[,k], rfk$cat),'coef']
      }
    }

    if(!is.na(os[1])){
      sim.resp <- sim.resp + dat[,os]
    }

    #     if(fam != 'gaussian'){
    #       if(fam=="Gamma"|fam=="gamma"){
    #         sim.resp <- -1/sim.resp
    #       }
    #
    #       if(fam=="binomial"){
    #         sim.resp <- exp(sim.resp)/(exp(sim.resp)+1)
    #
    #       }
    #       if(fam=="poisson"){
    #         sim.resp <- exp(sim.resp)
    #       }
    #     }

    dsim[,respvar] <- sim.resp

    sim.dat[[i]] <- dsim
  }

  ##### Make New Models #####

  if(!is.na(test.mod)){
    sim.mod <- ldply(sim.dat, function(x){
      mod <- update(test.mod, data=x)
      return(fixef(mod))
    })
  } else {
    sim.mod <- ldply(sim.dat, function(x){
      mod <- update(M, data=x)
      return(fixef(mod))
    })
  }

  a_ply(sim.mod, .margins = 2, function(x){
    hist(x[,1], xlab=names(x), main=NULL)
    abline(v=mean(x[,1]), lty=1, col='blue')
    abline(v=mean(x[,1])+2*sd(x[,1]), lty=1, col='blue')
    abline(v=mean(x[,1])-2*sd(x[,1]), lty=1, col='blue')
    abline(v=mean(x[,1])+1*sd(x[,1]), lty=1, col='blue')
    abline(v=mean(x[,1])-1*sd(x[,1]), lty=1, col='blue')
    abline(v=cf[cf$Code==names(x),'Coef'], lty=2, col='black')
    legend(x = 'topright', legend=c('Simulated Mean', 'True Coef'), lty=c(1,2),
           col=c('blue','black'))
  })

  if(to.return){
    return(list(sim.mod, sim.dat))
  }
}
