influenceJB <- function (model, group = NULL, select = NULL, obs = FALSE, gf = "single",
          count = FALSE, delete = TRUE, parallel=F, ncpus=NULL, cl=NULL, stopCL=T, ...)
{
  if (is.null(group) & !obs) {
    stop("Please specify either the 'group' parameter, or specify 'obs=TRUE'")
  }
  if (!is.null(group) & obs) {
    stop("Either specify the 'group' parameter, or specify 'obs=TRUE', but not both.")
  }
  data.adapted <- model.frame(model)
  #### Change Made Here ####
  if("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(weights)"] <-
      as.character(model@call$weights)}

  if("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(offset)"] <-
      as.character(model@call$offset)}

  if(sum(grepl("offset", names(data.adapted)))>0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <-
      gsub('offset\\(|\\)',"",names(data.adapted)[grep("offset", names(data.adapted))])}

  original.no.estex <- which(substr(names(lme4::fixef(model)), 1, 6) != "estex.")
  n.pred <- length(lme4::fixef(model)[original.no.estex])

  if (!obs) {
    grouping.names <- grouping.levels(model, group)
    n.groups <- length(grouping.names)
  }

  if (obs) { n.obs <- nrow(data.adapted) }

  or.fixed <- matrix(ncol = n.pred, nrow = 1, data =lme4::fixef(model)[original.no.estex])
  dimnames(or.fixed) <- list(NULL, names(lme4::fixef(model))[original.no.estex])
  or.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model)[original.no.estex])
  dimnames(or.se) <- list(NULL, names(lme4::fixef(model))[original.no.estex])
  or.vcov <- as.matrix(vcov(model)[original.no.estex, original.no.estex])
  dimnames(or.vcov) <- list(names(lme4::fixef(model)[original.no.estex]),
                            names(lme4::fixef(model)[original.no.estex]))
  or.test <- coef(summary(model))[original.no.estex, 3]

  if (!obs) {

    if (is.null(select)) {
      alt.fixed <- matrix(ncol = n.pred, nrow = n.groups, data = NA)
      dimnames(alt.fixed) <-
        list(grouping.names, names(lme4::fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.groups, data = NA)
      dimnames(alt.se) <- list(grouping.names,
                               names(lme4::fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.groups, data = NA)
      dimnames(alt.test) <- list(grouping.names,
                                 names(lme4::fixef(model))[original.no.estex])

      if(parallel){

        library(doParallel)

        if(is.null(cl) & !is.null(ncpus)){
          cl <- makeCluster(ncpus)
        }

        if(!is.null(cl)){
          registerDoParallel(cl)
        } else {registerDoParallel()}

        model.updated <-
          foreach (i = 1:n.groups, .inorder=F, .export=
                     c('exclude.influenceJB','lmer','glmer','glmerControl')) %dopar% {
                       #This is where a new model is created#
                       exclude.influenceJB(model, group, grouping.names[i], gf = gf,
                                           delete = delete)
                     }

        for(i in 1:n.groups){
          altered.no.estex <- which(substr(names(lme4::fixef(model.updated[[i]])),1,6)
                                    !="estex.")
          alt.fixed[i, ] <- as.matrix(lme4::fixef(model.updated[[i]])[altered.no.estex])
          alt.se[i, ] <- as.matrix(se.fixef(model.updated[[i]])[altered.no.estex])
          alt.vcov[[i]] <- as.matrix(vcov(model.updated[[i]])[altered.no.estex,
                                                         altered.no.estex])
          alt.test[i,] <- as.matrix(coef(summary(model.updated[[i]]))[,3]
                                    [altered.no.estex])

        }

      } else {
        for (i in 1:n.groups) {
          if (count == TRUE) {
            print(n.groups + 1 - i)
          }
          #This is where a new model is created#
          model.updated <- exclude.influenceJB(model, group,
                                               grouping.names[i], gf = gf, delete =delete)
          altered.no.estex<-which(substr(names(lme4::fixef(model.updated)),1,6)!="estex.")
          alt.fixed[i, ] <- as.matrix(lme4::fixef(model.updated)[altered.no.estex])
          alt.se[i, ] <- as.matrix(se.fixef(model.updated)[altered.no.estex])
          alt.vcov[[i]] <- as.matrix(vcov(model.updated)[altered.no.estex,
                                                         altered.no.estex])
          alt.test[i, ] <- as.matrix(coef(summary(model.updated))[,3][altered.no.estex])
        }
      }
    }

    if (!is.null(select)) {
      model.updated <- exclude.influenceJB(model, group, select, gf = gf, delete = delete)
      altered.no.estex <- which(substr(names(lme4::fixef(model.updated)),1,6) != "estex.")
      alt.fixed <- matrix(ncol = n.pred, nrow = 1,
                          data = lme4::fixef(model.updated)[altered.no.estex])
      dimnames(alt.fixed) <- list("Altered model",
                                  names(lme4::fixef(model.updated))[altered.no.estex])
      alt.se <- matrix(ncol=n.pred, nrow=1,data=se.fixef(model.updated)[altered.no.estex])
      dimnames(alt.se) <- list("Altered model",
                               names(lme4::fixef(model.updated))[altered.no.estex])
      alt.vcov <- list()
      alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
      dimnames(alt.vcov[[1]]) <- list(names(lme4::fixef(model.updated)[altered.no.estex]),
                                      names(lme4::fixef(model.updated)[altered.no.estex]))
      alt.test <- matrix(ncol = n.pred, nrow=1, data=coef(summary(model.updated))[,3]
                         [altered.no.estex])
      dimnames(alt.test) <- list("Altered model",
                                 names(lme4::fixef(model.updated))[altered.no.estex])
    }
  }

  if (obs) {
    if (is.null(select)) {
      alt.fixed <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.fixed) <- list(1:n.obs, names(lme4::fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.se) <- list(1:n.obs, names(lme4::fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.test) <- list(1:n.obs, names(lme4::fixef(model))[original.no.estex])

      if(parallel){

        library(doParallel)

        if(is.null(cl) & !is.null(ncpus)){
          cl <- makeCluster(ncpus)
        }

        if(!is.null(cl)){
          registerDoParallel(cl)
        } else {registerDoParallel()}


        model.updated<-
          foreach(i=1:n.obs, .inorder=F,
                  .export=c('exclude.influenceJB','lmer','glmer','glmerControl')) %dopar%{
                    exclude.influenceJB(model, obs = i)
                  }

        for (i in 1:n.obs){
          altered.no.estex <- which(substr(names(lme4::fixef(model.updated[[i]])),1,6)
                                    !="estex.")
          alt.fixed[i, ] <- as.matrix(lme4::fixef(model.updated[[i]])[altered.no.estex])
          alt.se[i, ] <- as.matrix(se.fixef(model.updated[[i]])[altered.no.estex])
          alt.vcov[[i]] <-as.matrix(vcov(model.updated[[i]])[altered.no.estex,
                                                        altered.no.estex])
          alt.test[i, ] <- as.matrix(coef(summary(model.updated[[i]]))[,3]
                                     [altered.no.estex])
        }
      } else {
        for (i in 1:n.obs) {
          if (count == TRUE) {
            print(n.obs + 1 - i)
          }
          model.updated <- exclude.influenceJB(model, obs = i)
          altered.no.estex<-which(substr(names(lme4::fixef(model.updated)),1,6)!="estex.")
          alt.fixed[i, ] <- as.matrix(lme4::fixef(model.updated)[altered.no.estex])
          alt.se[i, ] <- as.matrix(se.fixef(model.updated)[altered.no.estex])
          alt.vcov[[i]] <-as.matrix(vcov(model.updated)[altered.no.estex,
                                                        altered.no.estex])
          alt.test[i, ] <- as.matrix(coef(summary(model.updated))[, 3][altered.no.estex])
        }
      }
    }
    if (!is.null(select)) {
      model.updated <- exclude.influenceJB(model, obs = select)
      altered.no.estex <- which(substr(names(lme4::fixef(model.updated)),1,6) != "estex.")
      alt.fixed <- matrix(ncol = n.pred, nrow = 1,
                          data = lme4::fixef(model.updated)[altered.no.estex])
      dimnames(alt.fixed) <- list("Altered model",
                                  names(lme4::fixef(model.updated))[altered.no.estex])
      alt.se <-matrix(ncol=n.pred, nrow=1, data=se.fixef(model.updated)[altered.no.estex])
      dimnames(alt.se) <- list("Altered model",
                               names(lme4::fixef(model.updated))[altered.no.estex])
      alt.vcov <- list()
      alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, altered.no.estex])
      dimnames(alt.vcov[[1]]) <- list(names(lme4::fixef(model.updated)[altered.no.estex]),
                                      names(lme4::fixef(model.updated)[altered.no.estex]))
      alt.test <- matrix(ncol = n.pred, nrow = 1, data = coef(summary(model.updated))[,3]
                         [altered.no.estex])
      dimnames(alt.test) <- list("Altered model",
                                 names(lme4::fixef(model.updated))[altered.no.estex])
    }
  }

if(is.null(cl) & parallel){ stopImplicitCluster() } else if(stopCL & parallel) {
  stopCluster(cl)}

estex <- list(or.fixed = or.fixed, or.se = or.se, or.vcov = or.vcov, or.test = or.test,
                alt.fixed=alt.fixed, alt.se=alt.se, alt.vcov=alt.vcov, alt.test =alt.test)
  class(estex) <- "estex"
  return(estex)
}


exclude.influenceJB <- function (model, grouping = NULL, level = NULL, obs = NULL,
                                 gf = "single", delete = TRUE)
{
  data.adapted <- model.frame(model)
  ##### Change Made Here #####
  if("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(weights)"] <-
      as.character(model@call$weights)}

  if("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(offset)"] <-
      as.character(model@call$offset)}

  if(sum(grepl("offset", names(data.adapted)))>0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <-
      gsub('offset\\(|\\)',"",names(data.adapted)[grep("offset", names(data.adapted))])}

  added.variables <- character()

  if (!is.null(obs)) {
    if (!is.null(grouping) | !is.null(level)) {
      warning("Specification of the 'obs' parameter overrules specification of the 'grouping' and 'level' parameters.")
    }
    data.adapted <- data.adapted[-obs, ]
    model.updated <- update(model, data = data.adapted)
    return(model.updated)
  }

  if (delete == TRUE) {
    group.var <- which(names(data.adapted) == grouping)
    for (i in 1:length(level)) {
      data.adapted <- subset(data.adapted, data.adapted[, group.var] != level[i])
    }
    model.updated <- update(model, data = data.adapted)
    return(model.updated)
  }

  if (names(data.adapted)[2] != "intercept.alt") {
    data.adapted$intercept.alt <- ifelse(model@flist[, grouping] == level[1], 0, 1)
    data.adapted[, ncol(data.adapted) + 1] <- ifelse(model@flist[,grouping]==level[1],1,0)
    added.variables <- make.names(paste("estex.", as.character(level[1]), sep = ""))
    colnames(data.adapted)[ncol(data.adapted)] <- added.variables

    if (length(level) > 1) {
      for (i in 2:length(level)) {
        data.adapted$intercept.alt[model@flist[, grouping] == level[i]] <- 0
        data.adapted[, ncol(data.adapted) + 1] <-
          ifelse(model@flist[,grouping] == level[i], 1, 0)
        added.variables <- append(added.variables, values=
                                    make.names(paste("estex.",as.character(level[i]),
                                                     sep="")))
        colnames(data.adapted)[ncol(data.adapted)] <-
          added.variables[length(added.variables)]
      }
    }

    if (gf == "single") {
      grnr <- which(names(ranef(model)) == grouping)
      if (length(names(ranef(model)[[grnr]])) == 1) {
        model.updated <- update(model,formula=as.formula(
          paste(". ~ 0 + intercept.alt +", paste(added.variables,collapse="+"),"+ .",
                "- (1 |", grouping, ") + (0 + intercept.alt |", grouping, ")")),
          data = data.adapted)
      }

      if (length(names(ranef(model)[[grnr]])) > 1) {
        model.updated <- update(model,formula=as.formula(
          paste(". ~ 0 + intercept.alt + ", paste(added.variables,collapse="+"),"+ .",
                paste(" - (", paste(names(ranef(model)[[grnr]])[-1],collapse = "+"),
                      "|", grouping, ")"), " + (0 + intercept.alt +",
                paste(names(ranef(model)[[grnr]])[-1], collapse = "+"), "|", grouping,
                ")")), data = data.adapted)
      }
    }

    if (gf == "all") {
      delete.gf <- vector()

      for (i in 1:length(ranef(model))) {
        if (length(names(ranef(model)[[i]])) > 1) {
          delete.gf[i] <- paste("- (", paste(names(ranef(model)[[i]][-1]), collapse="+"),
                                "|", names(ranef(model))[i], ")")
        }

        if (length(names(ranef(model)[[i]])) == 1) {
          delete.gf[i] <- paste("- ( 1 |", names(ranef(model))[i], ")")
        }
      }

      delete.gf <- paste(delete.gf, collapse = " ")
      new.gf <- vector()

      for (i in 1:length(ranef(model))) {
        if (length(names(ranef(model)[[i]])) > 1) {
          new.gf[i] <- paste("+ (0 + intercept.alt +",
                             paste(names(ranef(model)[[i]][-1]), collapse = "+"),
                             "|", names(ranef(model))[i], ")")
        }

        if (length(names(ranef(model)[[i]])) == 1) {
          new.gf[i] <- paste("+ (0 + intercept.alt |", names(ranef(model))[i], ")")
        }
      }

      new.gf <- paste(new.gf, collapse = " ")
      model.updated <- update(model, formula =as.formula(
        paste(". ~ 0 + intercept.alt + ", paste(added.variables, collapse = "+"), "+ . ",
              delete.gf, new.gf)), data = data.adapted)
    }
  }

  if (names(data.adapted)[2] == "intercept.alt") {
    for (i in 1:length(level)) {
      data.adapted$intercept.alt[model@flist[, grouping] == level[i]] <- 0
      data.adapted[, ncol(data.adapted)+1] <- ifelse(model@flist[,grouping]==level[i],1,0)
      added.variables <- append(added.variables, values = make.names(
        paste("estex.", as.character(level[i]), sep = "")))
      colnames(data.adapted)[ncol(data.adapted)]<-added.variables[length(added.variables)]
    }
    model.updated <- update(model, formula = as.formula(
      paste(". ~ 0 + intercept.alt + ", paste(added.variables, collapse = "+"), "+ .")),
                            data = data.adapted)
  }

  return(model.updated)
}

