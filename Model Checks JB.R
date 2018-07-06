############################ Model Check for GLMM #############################
                            ## Jennifer Bufford ##
                           ## jlbufford@gmail.com ##
                             ## May 24, 2017 ##

###############################################################################

#Input: A GLMM model, the data used to make the model, optional specifiers
#       incl an optional name for the (optional) pdf file

#Requires: car, coefplot, ggplot2, lattice
#           and either lme4, nlme or glmmADMB
#          Optional: influence.ME, InfluenceJB, rmarkdown, HLMDiag, rsquaredglmm.R

#Output: Model summary and R2
#        Dotchart and qqplot of raw data
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

###############################################################################





##### Model Check #############################################################


model.check <- function(M, dat, min.unit=NA, make.pdf=F, make.markdown=F, name="Model",
                        infl=F, infl.obs=F, do.lm=F, off=T, respvar=NA, extra=NULL,
                        to.files="", parallel=F, ncpus=NULL){

  if(to.files!="" & substr(to.files,nchar(to.files),nchar(to.files)) != "/"){
    to.files <- paste(to.files, "/", sep="")
  }

  jb <- F

  library(ggplot2, quietly=T)
  if(infl|infl.obs) {library(influence.ME, quietly=T); infl <- T}
  library(car, quietly=T)
  library(coefplot, quietly=T)
  library(lattice, quietly=T)
  #Depending on call, may load lme4, nlme, glmmADMB, rmarkdown, HLMdiag or rsquaredglmm.R

  #Attractive theme for ggplot plots
  theme_jb <- function (base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(axis.text = element_text(size = 14), axis.title=element_text(size=16),
            axis.ticks=element_line(colour="black"),
            axis.line.x=element_line(color='black'),
            axis.line.y=element_line(color='black'),
            legend.key = element_rect(colour = "black"), panel.spacing.x = unit(3, 'mm'),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.border = element_blank(), panel.grid = element_blank(),
            strip.text = element_text(size=16), strip.text.y = element_text(angle=0),
            strip.background = element_rect(fill=NA, colour=NA, size = 0.2))
  }

  
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

    if(infl) {return(infldat)} else {return()}
  }


  ###### Extract Model Formula, Terms ######

  if('glm' %in% class(M)){

    if(missing(dat)) {dat <- M$data}

    fixvar <- attributes(M$terms)$term.labels
    Mterms <- fixvar[!grepl(':', fixvar)]
    randvar <- NA
    if(is.na(respvar)) { respvar <- names(M$model[1]) }
    fam <- M$family$family

    do.lm <- F

    if(is.na(min.unit)){return('Please specify a minimum unit for graphing.')}
  }

  if(any(class(M) %in% c("lmerMod","glmerMod"))) {

    library(lme4)

    #Retrieve data
    if(missing(dat)) {dat <- M@frame}

    #Set response variable
    if(is.na(respvar)) { respvar <- names(M@frame)[1] }

    #Set random variable(s) and min.unit
    randvar <- names(M@flist)
    if(any(grepl(':', randvar))) {randvar <- sub(':[[:print:]]*', '', randvar)} 
    #get lowest level of nested random effects

    if(is.na(min.unit) & length(randvar)==1) {min.unit <- randvar}
    if(is.na(min.unit)) {
      min.unit <- apply(dat[,randvar], MARGIN=2, FUN=function(x) {length(unique(x))})
      if(any(min.unit<nrow(dat))){
        min.unit <- names(min.unit[min.unit<nrow(dat)]) #prevents ID from being min.unit
        min.unit <- min.unit[1]
      } else {min.unit <- randvar[1]}
    }

    #Set fixed variables (excludes interactions)
    fixvar <- attributes(terms(M@frame))$term.labels
    fixvar <- fixvar[!(fixvar%in% randvar)]

    #Set list of terms (includes interactions)
    Mterms <- names(M@frame)[2:length(names(M@frame))]

    #Deal with weights or offsets
    if("(weights)" %in% Mterms) {
      jb <- T
      Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if(sum(grepl("offset", Mterms))>0) {
      jb <- T
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}
    if("glmerMod" %in% class(M)) { fam <- M@resp$family[[1]] }
  }

  if("lme" %in% class(M)) {

    library(nlme)

    if(missing(dat)) {dat <- M$data}

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

    if(is.na(min.unit)) {min.unit <- randvar[1]}
  }

  if('gls' %in% class(M)) {

    library(nlme)

    if(missing(dat)) {dat <- M$data}

    fixvar <- names(M$parAssign)[!names(M$parAssign) %in% "(Intercept)"]
    Mterms <- names(M$parAssign)[-grep(names(M$parAssign), pattern=":")][-1]
    randvar <- NA
    if(is.na(min.unit)){return('Please specify a minimum unit for graphing.')}
  }

  if('glmmadmb' %in% class(M)) {

    library(glmmADMB)

    if(missing(dat)) {
      dat <- M$frame
      warning('Data should be provided as an argument for glmmadmb models')}


    #Set response variable
    if(is.na(respvar)) { respvar <- names(M$frame)[1] }

    #Set random variable(s) and min unit
    randvar <- names(M$S)
    if(is.na(min.unit)) {min.unit <- randvar[1]}

    #Set fixed variables (excludes interactions)
    fixvar <- attributes(M$terms)$term.labels

    #Set list of terms (includes interactions)
    Mterms <- c(names(M$frame)[2:length(M$frame)], randvar)

    #Deal with weights or offsets
    #     if("(weights)" %in% Mterms) {Mterms <- Mterms[!(Mterms %in% "(weights)")]}
    if(sum(grepl("offset", Mterms))>0) {
      jb <- T
      Mterms <- gsub('offset\\(|\\)',"",Mterms[!(Mterms %in% "(offset)")])}

    #Set family
    fam <- M$family
  }

  
  ###### Print Summary, R-Squared, Wald Test ######
  
  print(summary(M))
  
  if(!any(c('glm','glmmadmb') %in% class(M))){
    
    try(source(paste(to.files, 'rsquaredglmm.R', sep='')))
    
    cat('\nR-Squared Values:\n')
    try(print(r.squared(M)))
    
    if(length(fixvar) > 0){
      
      cat('\nWald Chi-Square Test:\n')
      print(Anova(M))
    }
  }
  
  if('glm' %in% class(M)){
    cat('\nPseudo R-Squared (explained deviance):\n')
    print(1-M$deviance/M$null.deviance)
    
    cat('F-Test with Dispersion Estimate Based on Pearson Residuals:\n')
    print(Anova(M))
  }
  
  
  ##### Prep Data #####

  dat$Fit <- fitted(M)
  dat$Resid <- if('glmmadmb' %in% class(M)) {as.vector(M$residuals)} else {resid(M)}
  #for lme4 = deviance resids, for glmmadmb = pearson
  #not sure why but resid() isn't working for glmmadmb
  dat$MU <- dat[,min.unit]
  dat$RespVar <- dat[,respvar]


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

    if(any(class(dat[,i])=="character")) { dat[,i] <- factor(dat[,i]) }

    plot(x = dat[,i], y = dat[,respvar], ylab = respvar, xlab = i,
         main = paste('Plot of Raw Data:', i))
  }


  ###### Plot Residuals for Normality ######

  plot(fitted(M),dat$Resid, main = 'Plot of Fits vs Resids', xlab='Fitted',
       ylab='Residuals')
  qqnorm(dat$Resid, main = 'Normal Q-Q Plot of Residuals')

#   if("lme" %in% class(M)) { print(qqnorm(M, ~resid(.)|MU)) }


  ###### Plot Residuals for Homoscedacticity ######

  for (i in Mterms) {

    plot(dat[,i], dat$Resid, ylab = "Residuals", xlab = i,
         main = 'Plot of Residuals to Check for Homoscedasticity')
    abline(h=0)
  }


  ##### Plot Fits vs. Observed #####

  print(
  ggplot(dat, aes(x=RespVar, y = Fit)) + geom_point(size=3) +
    geom_abline(aes(intercept=0, slope=1), linetype=2) +
    xlab(paste("Observed", respvar)) + ylab(paste("Predicted", respvar)) +
    scale_y_continuous(expand=c(0.005,0)) + scale_x_continuous(expand=c(0.005,0))+
    theme_jb()
  )


  ###### Levene's Test ######

  for (i in Mterms) {

    cat('\n',i,'\n')
    if(!is.factor(dat[,i])){
      cat("Not a Factor, Levene's Test Not Applicable\n")
      next
    }

    if(length(levels(dat[,i]))<2) {cat("Single-level Factor\n"); next}
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
      scale_y_continuous(expand=c(0.005,0)) + theme_jb() +
      geom_smooth(aes(group=1)) +
      geom_text(aes(label=U)) + xlab('Mean of Residuals') +ylab('Std Dev of Residuals')+
      annotate("text", x = mean(MSDR$Mean), y = max(MSDR$SD)*0.95,
               label=paste("Spearman:", round(cor(MSDR$Mean, MSDR$SD,
                              use="complete.obs", method="spearman"),2)))
    )


  ###### Check for Var related to Fit ######

  print(
    ggplot(dat, aes(x = Fit, y = abs(Resid))) +
      scale_y_continuous(expand=c(0.005,0)) + theme_jb() +
      geom_point() + geom_smooth(aes(group=1)) +
      annotate("text", x = mean(dat$Fit), y = max(abs(dat$Resid))*0.95,
               label=paste("Spearman:", round(cor(dat$Fit, abs(dat$Resid),
                                  use="complete.obs", method="spearman"),2)))
    )


  ##### Check for Overdispersion #####

  if(any(class(M) %in% c("glmerMod","glmmadmb",'glm'))) {

    if(fam %in% c("poisson", 'binomial')) {

      if('glm' %in% class(M)){
        cat('\nTest for Overdispersion:\n')
        print(M$deviance/M$df.residual)
        cat('\n(should be close to 1)\n')
      }

      if(any(class(M) %in% c('glmerMod','glmmadmb'))){
        rdf <- nrow(model.frame(M)) -
          (sum(sapply(VarCorr(M),FUN=function(x){nrow(x)*(nrow(x)+1)/2}))
           + length(fixef(M)))
        rp <- residuals(M, type="pearson")
        pval <- pchisq(sum(rp^2), df=rdf, lower.tail=FALSE)

        cat('\nTest for Overdispersion:\n')
        print(c(chisq=sum(rp^2),ratio=sum(rp^2)/rdf,rdf=rdf,p=pval))
      }
    }
    ## modified by Orou Gaoue from: http://glmm.wikidot.com/faq
    ## number of variance parameters in an n-by-n variance-covariance matrix
  }


  ###### Leverage Plots ######

  if(length(randvar) < 3 & "lmerMod" %in% class(M)){

    try(library(HLMdiag, quietly=T))

    if(class(try(leverage(M, level=1)))!='try-error'){
      dat <- cbind(dat, leverage(M, level=1))

      print(
        ggplot(dat, aes(x = MU, y = overall)) +
          scale_y_continuous(name="Overall Leverage",expand=c(0,0.01)) +
          theme_jb() + ggtitle("Leverage Plot for LMM") +
          geom_point(size=4)
      )
    }
  }


  ##### Plot Random Effects #####

  if(any(class(M) %in% c("lmerMod","glmerMod"))) {
    print(dotplot(ranef(M, condVar=T))) }

  if('glmmadmb' %in% class(M)) {
    re <- ranef(M, condVar=T)[[1]]
    re <- re[order(re),]
    print(dotplot(re))}


  ##### Plot Coefficients #####

  if(any(class(M) %in% c("lmerMod", "glmerMod", 'glm','glmmadmb')) & length(fixvar)>0) {
    print(coefplot(M) + theme_jb())
  }


  ###### Create Linear Model on Fixed Effects ######

  if(length(fixvar) > 0 & do.lm) {

    modl <- paste(respvar, " ~ ", fixvar[1])

  if(length(fixvar) > 1){

    for(i in 2:length(fixvar)){
      modl <- paste(modl, "+ ", fixvar[i])
    }
  }

  if(any(class(M) %in% c('lmerMod',"glmerMod"))) {

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

  cat("\nLinear model for studentized residuals, added-variable plots:\n")
  print(summary(lmm))


  ###### Studentized Residuals (and Others) ######

  infIndexPlot(lmm, main="Diagnostic Plots of LM on Fixed Effects Only")


  ###### Added-Variable (Partial Regression) Plots ######

  avPlots(lmm, id.method="mahal", id.n=3,
          main="Added-Variable plots of LM on Fixed Effects Only")

  par(mfrow = c(1,1))

  }


  ##### Additional Plots for (G)LM ######

  if(any(class(M) %in% c('glm','lm'))){

    plot(M)

    infIndexPlot(M)

    avPlots(M, id.method="mahal", id.n=3)

    par(mfrow = c(1,1))
  }


  ###### Influential Observations ######

  if(infl & !any(class(M) %in% c('lmerMod',"glmerMod"))) {

    warning(paste("Influence.ME is not implemented for", class(M)))
    infl <- F
  }

  if(infl){

    if(!'influenceJB' %in% ls() & (jb|parallel)){
      source(paste(to.files, "InfluenceJB.R", sep=""))
    }
    #If not already loaded, will load InfluenceJB.R to calc infl for model w/ weights

    if (min.unit %in% Mterms) {
      Infl.mu <- if(jb|parallel) {
        influenceJB(model=M, group=min.unit, parallel=parallel, ncpus=ncpus)
      } else { influence(model = M, group = min.unit) }
      plot(Infl.mu, which="cook", sort=T,
           cutoff=(4/(length(unique(dat$MU))-length(c(fixvar,randvar))-1)),
           main = paste("Cook's D by", min.unit))
    } else {
      cat('\nMinimum unit must be a model term to calculate influence\n')
      Infl.mu <- NA
    }

    if(infl.obs) {

      Infl <- if(jb|parallel) {
        influenceJB(model=M, obs=T, parallel=parallel, ncpus=ncpus)
      }
      else { influence(model = M, obs=TRUE) }

      plot(Infl, which="cook", sort=T, cutoff=(4/(nrow(dat)-length(c(fixvar,randvar))-1)),
           main="Cook's D by Observation")

      plot(Infl, which="dfbetas", cutoff=2/sqrt((nrow(dat)-length(c(fixvar,randvar)) -1)),
           sort=T, to.sort=colnames(M@pp$X)[2], main="Dfbetas")

      if(any(cooks.distance(Infl) > 0.1)){
        cat("\nThe following observations have a Cook's D > 0.1:\n")
        print(dat[cooks.distance(Infl) > 0.1,])

        M2 <- update(M, .~., data=dat[cooks.distance(Infl) < 0.1,])
        print(coefplot(M2) + theme_jb() +
                ggtitle('Coefficients of Model without Influential Points'))
      }

    }


    if(infl.obs) { return(list(Infl, Infl.mu)) } else { return(Infl.mu) }
  }


  ###### Close pdf ######

  if(make.pdf & off) { dev.off() }

}
