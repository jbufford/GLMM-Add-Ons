############################ Model Check for GLMM #############################
                            ## Jennifer Bufford ##
                           ## jlbufford@gmail.com ##
                            ## February 16, 2015 ##

###############################################################################

#Input: A GLMM model, the data used to make the model, optional specifiers
#       incl an optional name for the (optional) pdf file

#Requires: car, influence.ME, InfluenceJB, ggplot2, HLMdiag, coefplot, lattice, stringr
#           and either lme4, nlme or glmmADMB
#          Optional: influence.ME, InfluenceJB, RMarkdown

#Output: Dotchart and qqplot of raw data
#        Graphs of all terms vs the response variable (raw data)
#        Graphs of resids vs fits, resids as qqnorm
#        Resids vs all model terms
#        Levene's test results for categorical model terms
#        Graphs of mean vs abs(resids) and mean vs sd of resids
#        Test for overdispersion in poisson and binomial models
#        Graph of leverage by obs (for lme4 models with less than 3 rand var only)
#        Caterpillar plot of random effects in lme4 or glmmADMB
#        Graph of coefficients and standard errors (in lme4 only)
#        Graphs of studentized residuals, partial regression plots for an equivalent GLM
#   `    Graphs of Cook's D and dfbetas by obs and/or min.unit (optional, lme4 only)
#        Graphs saved as a pdf or all output saved as a markdown document (optional)
#        Returns as a list influence.ME objects by obs and sp

#model.check runs a series of model checks for mixed models (see output above)
#levenes performs a levene's test of the residuals by grouping factors
#Bootstrap function bootstraps data/model to generate 95% CI
#   Re-runs if model was unstable and didn't work
#   Optionally, re-samples if data did not include a min number of reps in every level of #       a variable (although that doesn't seem to help)

###############################################################################





##### Model Check #############################################################


model.check <- function(M, dat, min.unit, make.pdf=F, make.markdown=F, name="Model",
                        infl=F, infl.obs=F, do.lm=F, off=T, respvar=NA, extra=NULL,
                        to.files=""){

  if(to.files!="" & substr(to.files,nchar(to.files),nchar(to.files)) != "/"){
    to.files <- paste(to.files, "/", sep="")
  }

  jb <- F

  if(is.null(min.unit)){
    warning('Please specify a minimum unit for diagnostic plots')
    break
  }

  library(ggplot2, quietly=T)
  if(infl) {library(influence.ME, quietly=T)}
  library(car, quietly=T)
  library(HLMdiag, quietly=T)
  library(coefplot, quietly=T)
  library(lattice, quietly=T)

  print(summary(M))


  ###### Create Markdown PDF ######

  if (make.markdown) {

    library(rmarkdown, quietly=T)

    if(grepl('/', name)){
      folder <- paste(strsplit(name, "/")[[1]][1:(length(strsplit(name, "/")[[1]])-1)],
                      collapse="/")
      wd <- paste(getwd(), folder, sep="/")
    } else {wd <- getwd()}

    name2 <- strsplit(name, "/")[[1]][length(strsplit(name, "/")[[1]])]
    name2 <- gsub(" ", "_", name2)

    render(paste(to.files, "Model_Checks_JB.Rmd", sep=""),
           output_file=paste(name2, "_Model_Checks_", Sys.Date(), ".pdf", sep=""),
           output_dir=wd)

    if(infl|infl.obs) {return(infldat)} else {return()}
  }


  ##### Prep Data #####

  dat$Fit <- fitted(M)
  dat$Resid <- resid(M) #for lme4 = deviance resids, for glmmadmb = pearson
  dat$MU <- dat[,min.unit]


  ###### Extract Model Formula, Terms ######

  if(class(M)=="lmerMod" | class(M)=="glmerMod") {

    library(lme4)

    if(is.na(respvar)) { respvar <- names(M@frame)[1] }
    randvar <- names(M@flist)
    fixvar <- attributes(terms(M@frame))$term.labels
    fixvar <- fixvar[!(fixvar%in% randvar)]
    Mterms <- names(M@frame)[2:length(names(M@frame))]
    if("(weights)" %in% Mterms) {
      jb <- T
      Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if(sum(grepl("offset", Mterms))>0) {
      jb <- T
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}
    if(class(M)=="glmerMod") { fam <- M@resp$family[[1]] }
  }

  if(class(M)=="lme") {

    library(nlme)

    if(is.na(respvar)) { respvar <- as.character(attributes(M$terms)$variables[2]) }

    fixvar <- attributes(M$terms)$term.labels

    #If there are fixed effects
    if(length(attributes(M$terms)$term.labels) > 0) {
      Mterms <- c(as.character(attributes(M$terms)$variables[3:
                                            length(attributes(M$terms)$variables)]),
                attributes(M$modelStruct$reStruct)$names, extra)
      if(sum(grepl("offset", Mterms))>0) {
        jb <- T
        Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}

      randvar <- c(attributes(M$modelStruct$reStruct)$names, extra)

    } else { #if there are no fixed effects except intercept
      randvar <- c(attributes(M$modelStruct$reStruct)$names, extra)
      Mterms <- randvar
    }
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
      jb <- T
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}
    dat$Resid <- as.vector(M$residuals) #pearson's residuals
    #(not sure why but resid() isn't working)
    fam <- M$family
  }


  ###### Create PDF ######

  if (make.pdf) {

    plot(dat$Fit, dat$Resid)

    dev.off()

    pdf(file=paste(name, " Model Checks ", Sys.Date(), ".pdf", sep=""),
        width=16, height=8)
  }


  ###### Plot Raw Data ######

  dotchart(dat[, respvar], xlab = respvar, main= paste("Cleveland Dotplot of", respvar))

  qqnorm(dat[, respvar], main = "Normal Q-Q Plot of Response Variable")


  ##### Plot Raw Terms #####

  for (i in Mterms){

    if(class(dat[,i])=="character") { dat[,i] <- factor(dat[,i]) }

    plot(x = dat[,i], y = dat[,respvar], ylab = respvar, xlab = i,
         main = paste('Plot of Raw Data'))
  }


  ###### Plot Residuals for Normality ######

  plot(fitted(M),dat$Resid, main = 'Plot of Fits vs Resids', xlab='Fitted',
       ylab='Residuals')
  qqnorm(dat$Resid, main = 'Normal Q-Q Plot of Residuals')

#   if(class(M)=="lme") { print(qqnorm(M, ~resid(.)|MU)) }


  ###### Plot Residuals for Homoscedacticity ######

  for (i in Mterms) {

    plot(dat[,i], dat$Resid, ylab = "Residuals", xlab = i,
         main = 'Plot of Residuals to Check for Homoscedasticity')
    abline(h=0)
  }


  ###### Levene's Test ######

  for (i in Mterms) {

    print(i)
    if(!is.factor(dat[,i])){
      warning("Not a Factor, Levene's Test Not Applicable")
      next
    }

    if(length(levels(dat[,i]))<2) {warning ("Single-level Factor"); next}
    print(leveneTest(dat$Resid, dat[,i]))
  }


  ###### Check for Var related to Mean ######

  MR <- aggregate(dat$Resid, by=list("U"=dat$MU), FUN=mean)
  names(MR)[2] <- "Mean"

  SDR <- aggregate(dat$Resid,by=list("U"=dat$MU), FUN=sd)
  names(SDR)[2] <- "SD"

  MSDR <- merge(MR, SDR)

  MSDR[is.na(MSDR$SD), "SD"] <- 0

  print(
    ggplot(MSDR, aes(x = Mean, y = SD)) +
      scale_y_continuous(expand=c(0,0.005)) + theme_bw() +
      geom_smooth(aes(group=1)) +
      geom_text(aes(label=U)) + xlab('Mean of Residuals') +ylab('Std Dev of Residuals')+
      annotate("text", x = mean(MSDR$Mean), y = max(MSDR$SD)*0.95,
               label=paste("Spearman:", round(cor(MSDR$Mean, MSDR$SD,
                              use="complete.obs", method="spearman"),2)))
    )


  ###### Check for Var related to Fit ######

  print(
    ggplot(dat, aes(x = Fit, y = abs(Resid))) +
      scale_y_continuous(expand=c(0,0.005)) + theme_bw() +
      geom_point() + geom_smooth(aes(group=1)) +
      annotate("text", x = mean(dat$Fit), y = max(abs(dat$Resid))*0.95,
               label=paste("Spearman:", round(cor(dat$Fit, abs(dat$Resid),
                                  use="complete.obs", method="spearman"),2)))
    )


  ##### Check for Overdispersion #####

  if((class(M)=="glmerMod" | class(M)=="glmmadmb")) {

    if(fam %in% c("poisson", 'binomial')) {

      rdf <- nrow(model.frame(M)) -
        (sum(sapply(VarCorr(M),FUN=function(x){nrow(x)*(nrow(x)+1)/2}))
         + length(fixef(M)))
      rp <- residuals(M, type="pearson")
      pval <- pchisq(sum(rp^2), df=rdf, lower.tail=FALSE)

      print('Test for Overdispersion:')
      print(c(chisq=sum(rp^2),ratio=sum(rp^2)/rdf,rdf=rdf,p=pval))
    }
    ## modified by Orou Gaoue from: http://glmm.wikidot.com/faq
    ## number of variance parameters in an n-by-n variance-covariance matrix
  }


  ###### Leverage Plots ######

  if(length(randvar) < 3 & class(M)=="lmerMod"){

    dat <- cbind(dat, leverage(M, level=1))

    print(
        ggplot(dat, aes(x = MU, y = overall)) +
          scale_y_continuous(name="Overall Leverage",expand=c(0,0.01)) +
          theme_bw() + ggtitle("Leverage Plot for LMM") +
          geom_point(size=4)
      )
    }


  ##### Plot Random Effects #####

  if(class(M)=="lmerMod" | class(M)=="glmerMod") {
    print(dotplot(ranef(M, condVar=T))) }

  if(class(M)=='glmmadmb') {
    re <- ranef(M, condVar=T)[[1]]
    re <- re[order(re),]
    print(dotplot(re))}


  ##### Plot Coefficients #####

  if((class(M)=="lmerMod" | class(M)=="glmerMod") & length(fixvar)>0) {
    print(coefplot(M)) }


  ###### Create Linear Model on Fixed Effects ######

  if(length(fixvar) > 0 & do.lm) {

    modl <- paste(respvar, " ~ ", fixvar[1])

  if(length(fixvar) > 1){

    for(i in 2:length(fixvar)){
      modl <- paste(modl, "+ ", fixvar[i])
    }
  }

  if(class(M)=='lmerMod' | class(M)=="glmerMod") {

    if("(weights)" %in% names(M@frame)) {

      dat$Wts <- dat[,as.character(M@call$weights)]

      if("(offset)" %in% names(M@frame)) {

        lmm <- lm(modl, weights=Wts, offset=LM, dat)

      } else {   lmm <- lm(modl, weights=Wts, dat)  }

    } else {

      if("(offset)" %in% names(M@frame)) {

        lmm <- lm(modl, offset=LM, dat)

      } else {  lmm <- lm(modl, dat) }
    }

  } else { lmm <- lm(modl, dat) }

  print("Linear model for studentized residuals, added-variable plots:")
  print(summary(lmm))


  ###### Studentized Residuals (and Others) ######

  infIndexPlot(lmm, main="Diagnostic Plots of LM on Fixed Effects Only")


  ###### Added-Variable (Partial Regression) Plots ######

  avPlots(lmm, id.method="mahal", id.n=3,
          main="Added-Variable plots of LM on Fixed Effects Only")

  }


  ###### Influential Observations ######

  if(infl & !(class(M)=='lmerMod' | class(M)=="glmerMod")) {

    warning(paste("Influence.ME is not implemented for", class(M)))
    infl <- F
  }

  if(infl){

    if(!'influenceJB' %in% ls() & jb){
      source(paste(to.files, "InfluenceJB.R", sep=""))
    }
    #If not already loaded, will load InfluenceJB.R to calc infl for model w/ weights


    if(infl.obs) {

      Infl <- if(jb) {influenceJB(model = M, obs=TRUE)} else {
        influence(model = M, obs=TRUE)}
      plot(Infl, which="cook", sort=T, cutoff=(4/(nrow(dat)-length(c(fixvar,randvar))
                                                  -1)), main="Cook's D for LMM")
      plot(Infl, which="dfbetas", cutoff=2/sqrt((nrow(dat)-length(c(fixvar,randvar))
                                                 -1)), main="Dfbetas for LMM")
    }


    if (min.unit %in% Mterms) {
      Infl.mu <- if(jb) {influenceJB(model = M, group = min.unit)} else {
        influence(model = M, group = min.unit)
      }
      plot(Infl.mu, which="cook", sort=T,
           cutoff=(4/(length(unique(dat$MU))-length(c(fixvar,randvar))-1)),
           main = paste("Cook's D by", min.unit))
    }

    if(make.pdf & off) { dev.off() }

    if(infl.obs) { return(list(Infl, Infl.mu)) } else { return(Infl.mu) }
  }


  ###### Close pdf ######

  if(make.pdf & off) { dev.off() }

}





##### Levene's Tests ##########################################################


levenes <- function(M, dat, extra=NULL) {

  library(car)
  library(stringr)


  ###### Extract Model Formula, Terms ######

  dat$Resid <- resid(M)

  if(class(M)=="lmerMod" | class(M)=="glmerMod") {

    Mterms <- names(M@frame)[2:length(names(M@frame))]
    if("(weights)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if("(offset)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(offset)")]}
  }

  if(class(M)=="lme") {

    #If there are fixed effects
    if(length(attributes(M$terms)$term.labels) > 0) {
      Mterms <- c(as.character(attributes(M$terms)$variables
                               [3:length(attributes(M$terms)$variables)]),
                  attributes(M$modelStruct$reStruct)$names)
    } else { #if there are no fixed effects except intercept
      Mterms <- attributes(M$modelStruct$reStruct)$names
    }
  }

  if(class(M)=='gls') {

    Mterms <- names(M$parAssign)[-grep(names(M$parAssign), pattern=":")][-1]
  }

  if(class(M)=='glmmadmb') {

    Mterms <- c(names(M$frame)[2:length(M$frame)], names(M$S))
    #     if("(weights)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    #     if("(offset)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(offset)")]}
    if(is.na(dat$Resid[1])) {print("manual"); dat$Resid <- as.vector(M$residuals)}
  }


  ##### Calculate Levene's Test #####

  print("Residals")

  for (i in Mterms) {

    print(i)
    if(!is.factor(dat[,i])){
      warning("Not a Factor, Levene's Test Not Applicable")
      next
    }

    if(length(levels(dat[,i]))<2) {warning("Single-level Factor"); next}
    print(leveneTest(dat$Resid, dat[,i]))
  }

  if(length(extra)>0) {
    for(i in extra) {

      print(i)
      print(leveneTest(dat$Resid, dat[,i]))
    }
  }
}



##### Bootstrap Function ######################################################


library(boot)

bootm <- function(dat, i, M, min.unit=NA, min.reps=0) {

  if(class(M)=="lmerMod" | class(M)=="glmerMod") {library(lme4)}
  if(class(M)=="lme") {library(nlme)} #to avoid recursive errors

  dats <- dat[i,]

  if(!is.na(min.unit) & min.reps>0){
    if(sum(xtabs( ~ dats[, min.unit])< min.reps)>0){
      boot.samp <- sample(1:nrow(dat), nrow(dat), replace=T)
      warning("Error - Re-sampling")
      bootm(dat, boot.samp, M, min.unit, min.reps)
    }
  }
  mod <- try(update(M, data=dats), T)

  if(class(mod)=="try-error") {
    boot.samp <- sample(1:nrow(dat), nrow(dat), replace=T)
    warning("Error - Re-running")
    bootm(dat, boot.samp, M, min.unit, min.reps)

  } else {fixef(mod)}
}