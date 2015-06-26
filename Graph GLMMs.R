########################### Graphing (G)LMM Results ###########################
                            ## Jennifer Bufford ##
                      ## jennifer.bufford@lincoln.ac.nz ##
                             ## June 2, 2015 ##

###############################################################################

#Input: a data.frame or list with CI, a vector of code names, a vector of graph names
#         (code and graph names must be ordered to match)
#       a list of data.frames w/ CI, a vector of names for the data.frames

#Requires: boot, ggplot2, grid

#Output: data.frame with variable names for graphing, coefs, and 95% CI

#Bootstrap function bootstraps data/model to generate 95% CI
#   Re-runs if model was unstable and didn't work
#   Optionally, re-samples if data did not include a min number of reps in every level of #       a variable (although that doesn't seem to help)
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


getci <- function(dat, CodeN = NA, GraphN = NA){

  if(class(dat)=='boot'){

    library(boot)

    dat.bt <- dat
    dat <- list()

    for(i in 1:length(dat.bt$t0)) {
      dat[[i]] <- boot.ci(dat.bt, index=i, type=c("norm", "basic"))
      names(dat)[i] <- names(dat.bt$t0)[i]
    }
  }


  if(class(dat)=="list") {

    for(i in names(dat)){

      x <- data.frame("Coef"=dat[[i]]$t0, "LowCI"=dat[[i]]$normal[2],
                      "HiCI"=dat[[i]]$normal[3])
      x$Sig <- ifelse((x$HiCI<0 | x$LowCI>0), "sig", 'NS')

      if(!is.na(CodeN[1]) & !is.na(GraphN[1])){
        x$CName <- as.character(GraphN[which(CodeN==i)])
      } else { x$CName <- i }

      if(i==names(dat)[1]) {cidat <- x} else {cidat <- rbind(cidat, x)}
    }

    if(!is.na(CodeN[1]) & !is.na(GraphN[1])){
      cidat$CName <- factor(cidat$CName,
                            levels=unique(GraphN[GraphN %in% cidat$CName]),order=T)
    }

    return(cidat)
  }

  if(class(dat)=='data.frame') {

    dat$Sig <- ifelse((dat$HiCI<0 | dat$LowCI>0), "sig", 'NS')

    if(!is.na(CodeN[1]) & !is.na(GraphN[1])){
      for (i in dat$CMName) {

        dat[dat$CMName==i, "CName"] <- as.character(GraphN[which(CodeN==i)])
      }

      dat$CName <-factor(dat$CName, levels=unique(GraphN[GraphN %in% dat$CName]),
                         order=T)
    }

    return(dat)
  }
}



##### Create Single Data.frame with CI ########################################


merge.ci <- function(ci.list, ci.names=NA, GraphN=NA){

  for(i in 1:length(ci.list)){

    if(!is.na(ci.names[1])){ ci.list[[i]]$Var <- ci.names[i] } else {

        if(length(names(ci.list))>0){
          ci.list[[i]]$Var <- names(ci.list)[i] } else {ci.list[[i]]$Var <- i}
        }

    if(i==1) {ci.merg <- ci.list[[1]]} else {ci.merg <- rbind(ci.merg, ci.list[[i]])}
  }

  #Format
  row.names(ci.merg) <- 1:nrow(ci.merg)

  if(!is.na(ci.names[1])){
    ci.merg$Var <- factor(ci.merg$Var, levels=ci.names, order=T)
  }

  if(!is.na(GraphN[1])){
    ci.merg$CName <- factor(ci.merg$CName,
                            levels=unique(GraphN[GraphN %in% ci.merg$CName]),
                            order=T)
  }

  return(ci.merg)
}



##### Plot Coefs ##############################################################


plot.coef <- function(dat, get.ci=T, CodeN=NA, GraphN=NA, ci.names=NA,
                      facet.scale = 'free', intercept = T, ...){

  library(ggplot2)
  library(grid)

  if(class(dat)=="list" & class(dat[[1]])=="list"){

    datl <- list()

    if(get.ci){
      for(i in 1:length(dat)){
        datl[[i]] <- getci(dat[[i]], CodeN, GraphN)
      }
    }

    toplot <- merge.ci(datl, ci.names, GraphN)

  } else {

    if(get.ci){ toplot <- getci(dat, CodeN, GraphN) } else { toplot <- dat }
  }

  if(!intercept) {
    toplot <- toplot[!toplot$CName %in% c('Intercept','(Intercept)','intercept'),]}

  coef.plot <- ggplot(data=toplot, aes(x=Coef, y=CName, shape=Sig)) +
    scale_color_manual(guide=F, name="Significance") +
    xlab("Effect Size") + ylab("") + theme_bw() + geom_vline(aes(x=0)) +
    scale_shape_manual(guide=F, values=c(1,16)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 14), axis.title = element_text(size = 16),
          strip.text = element_text(size=16),
          strip.background = element_rect(fill="white", color="White"),
          panel.margin.x = unit(3, 'mm'), ...)

  if(class(dat)=="list" & class(dat[[1]])=="list"){

    print(
      coef.plot + geom_point(size=4) +
        geom_errorbarh(aes(xmin=LowCI, xmax=HiCI, height=0.15)) +
        facet_grid(.~Var, scales = facet.scale)
    )
  } else {

    print(
      coef.plot + geom_point(size=4) +
        geom_errorbarh(aes(xmin=LowCI, xmax=HiCI, height=0.15))
    )
  }

}


