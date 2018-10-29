########################### Graphing (G)LMM Results ###########################
                            ## Jennifer Bufford ##
                      ## jennifer.bufford@lincoln.ac.nz ##
                            ## October 24, 2018 ##

###############################################################################

#Input: a data.frame or list with CI, a vector of code names, a vector of graph names
#         (code and graph names must be ordered to match)
#       a list of data.frames w/ CI, a vector of names for the data.frames

#Requires: boot, ggplot2, grid

#Output: data.frame with variable names for graphing, coefs, and 95% CI

#Bootstrap function bootstraps data/model to generate 95% CI
#   Re-runs if model was unstable and didn't work
#   Optionally, re-samples if data did not include a min number of reps in every level of
#       a variable (although that doesn't seem to help)
#getci takes a boot object or a list from boot.ci or a data.frame from as.mcmc
#         creates a data.frame of variable names (for graphing) and coef and 95% CI
#merge.ci takes a list of data.frames with CI and a vector of names
#         creates a single data.frame for plotting facetted graphs
#plot.coef plots the coefficients and 95% CI from one or multiple models

###############################################################################



##### Bootstrap Function ######################################################


bootm <- function(dat, i, M, min.unit=NA, min.reps=0) {

  library(boot)

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



##### Extract CI from Bootstrap or MCMC #######################################


getci <- function(dat, CodeN, GraphN, mod, boot.type='perc'){

  if(class(dat)=='boot'){

    library(boot)

    dat.bt <- dat
    dat <- list()

    for(i in 1:length(dat.bt$t0)) {
      dat[[i]] <- boot.ci(dat.bt, index=i, type=boot.type)
      names(dat)[i] <- names(dat.bt$t0)[i]
    }
  }


  if(class(dat)=="list") {

    for(i in names(dat)){

      x <- data.frame("Coef"=dat[[i]]$t0,
                      "LowCI"=unlist(dat[[i]][4])[length(unlist(dat[[i]][4]))-1],
                      "HiCI"=unlist(dat[[i]][4])[length(unlist(dat[[i]][4]))])
      x$Sig <- ifelse((x$HiCI<0 | x$LowCI>0), "sig", 'NS')

      if(!missing(CodeN) & !missing(GraphN)){
        x$CName <- if(any(CodeN==i)){
          as.character(GraphN[which(CodeN==i)])} else {i}

      } else { x$CName <- i }

      if(i==names(dat)[1]) {cidat <- x} else {cidat <- rbind(cidat, x)}
    }

    if(!missing(CodeN) & !missing(GraphN)){
      cidat$CName <- factor(cidat$CName, order=T,
                         levels=unique(c(GraphN[GraphN %in% cidat$CName], cidat$CName)))

    }

    return(cidat)
  }

  if(class(dat)=='data.frame') {

    dat$Sig <- ifelse((dat$HiCI<0 | dat$LowCI>0), "sig", 'NS')

    if(!missing(CodeN) & !missing(GraphN)){
      for (i in dat$CMName) {

        dat[dat$CMName==i, "CName"] <- if(any(CodeN==i)){
          as.character(GraphN[which(CodeN==i)])} else {i}
      }

      dat$CName <-factor(dat$CName, order=T,
                         levels=unique(c(GraphN[GraphN %in% dat$CName], dat$CName)))
    } else {dat$CName <- dat$CMName}

    return(dat)
  }

    if(class(dat)=='matrix') { #to process output from confint.merMod

      if(missing(mod)) {stop('Model required to process output from confint.merMod')}

      dat <- data.frame(dat)
      dat$CMName <- rownames(dat)
      dat <- dat[dat$CMName %in% names(fixef(mod)),]
      dat$Coef <- fixef(mod)
      names(dat)[1:2] <- c('LowCI','HiCI')

      dat$Sig <- ifelse((dat$HiCI<0 | dat$LowCI>0), "sig", 'NS')

      if(!missing(CodeN) & !missing(GraphN)){
        for (i in dat$CMName) {

          dat[dat$CMName==i, "CName"] <- if(any(CodeN==i)) {
            as.character(GraphN[which(CodeN==i)])} else {i}
        }

        dat$CName <- factor(dat$CName, order=T,
                            levels=unique(c(GraphN[GraphN %in% dat$CName], dat$CName)))

      } else {dat$CName <- dat$CMName}

    return(dat)
  }
}



##### Create Single Data.frame with CI ########################################


merge.ci <- function(ci.list, ci.names, GraphN){

  for(i in 1:length(ci.list)){

    if(!missing(ci.names)){ ci.list[[i]]$Var <- ci.names[i] } else {

        if(length(names(ci.list))>0){
          ci.list[[i]]$Var <- names(ci.list)[i] } else {ci.list[[i]]$Var <- i}
        }

    if(i==1) {ci.merg <- ci.list[[1]]} else {ci.merg <- rbind(ci.merg, ci.list[[i]])}
  }

  #Format
  row.names(ci.merg) <- 1:nrow(ci.merg)

  if(!missing(ci.names)){
    ci.merg$Var <- factor(ci.merg$Var, levels=ci.names, order=T)
  }

  if(!missing(GraphN)){
    ci.merg$CName <- factor(ci.merg$CName, order=T,
                            levels=unique(c(GraphN[GraphN %in% ci.merg$CName],
                                            as.character(ci.merg$CName))))
  }

  return(ci.merg)
}



##### Plot Coefs ##############################################################


plot.coef <- function(dat, get.ci=T, boot.type='perc', CodeN, GraphN, mod, ci.names,
                      facet.scale = 'free', intercept = T, return.plot=F, ...){

  library(ggplot2)
  library(grid)

  #Attractive theme for ggplot
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

  if(class(dat)=="list"){

    datl <- list()

    if(get.ci){

      if(!missing(mod) & class(dat[[1]])=='matrix'){
        if(length(mod)<length(dat)) {stop('Must provide a model for each confint object listed in the same order as the confint objects')}
      }

      for(i in 1:length(dat)){

        if(!missing(mod)) {
          datl[[i]] <- getci(dat[[i]], CodeN, GraphN, mod=mod[[i]], boot.type)
        }
        else { datl[[i]] <- getci(dat[[i]], CodeN, GraphN, mod, boot.type) }
      }
    }

      toplot <- merge.ci(datl, ci.names, GraphN)

    } else {

      if(get.ci){ toplot <- getci(dat, CodeN, GraphN, mod, boot.type) } else {toplot<-dat}
    }

  if(!intercept) {
    toplot <- toplot[!toplot$CName %in% c('Intercept','(Intercept)','intercept'),]}

  SigShape <- sort(unique(ifelse(toplot$Sig=='NS', 1, 16)))

  coef.plot <- ggplot(data=toplot, aes(x=Coef, y=CName, shape=Sig)) +
    scale_color_manual(guide=F, name="Significance") +
    xlab("Effect Size") + ylab("") + geom_vline(aes(xintercept=0)) +
    scale_shape_manual(guide=F, values=SigShape) + theme_jb()

  if(class(dat)=="list"){

    coef.plot <- coef.plot + geom_point(size=4) + facet_grid(.~Var, scales =facet.scale) +
      geom_errorbarh(aes(xmin=LowCI, xmax=HiCI, height=0.15))

    if(return.plot) {return(coef.plot)} else {print(coef.plot)}

  } else {

    coef.plot <- coef.plot + geom_point(size=4) +
      geom_errorbarh(aes(xmin=LowCI, xmax=HiCI, height=0.15))

    if(return.plot) {return(coef.plot)} else {print(coef.plot)}

  }
}


