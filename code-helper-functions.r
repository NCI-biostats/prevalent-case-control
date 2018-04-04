#####################################
# Marlena Maziarz                   #
# March 29, 2018                    #
# marlena.maziarz@nih.gov (current) #
# marlenam@uw.edu (permanent)       #
#####################################


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set.program.control(...)
#  this function generates a list of parameters and settings that are needed for fitting the model
#  the list is passed to all functions.
#
# data should have at least the following columns:
#   case.status    0 = control, 1 = incident case, 2 = prevalent case
#   backward.time  numeric ranging from >= 0
#   all covariates listed in logistic.vars and survival.vars
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.program.control <- function(data,            # one entry per individual
                                logistic.vars,   # list of names, as in the dataset
                                survival.vars = NULL,   # list of names, as in the dataset
                                const.bh = F,
                                param.vec.fitting.start.values = NULL,
                                xi.a = NULL,
                                print.ok = F){

    # data shouldn't contain any NA's
    if(sum(is.na(data)) > 0){
        print('The data cannot have any missing values')
    }

    ix.0 <- data$case.status == 0
    ix.1 <- data$case.status == 1
    ix.2 <- data$case.status == 2

    n0 <- sum(ix.0)
    n1 <- sum(ix.1)
    n2 <- sum(ix.2)

    # column indices of the covariates in the dataset, order should correspond to the parameter vectors below
    logistic.x.col.ix <- which(colnames(data) %in% logistic.vars)
    survival.x.col.ix <- which(colnames(data) %in% survival.vars)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # !!! CAUTION - THIS PART NOT GENERAL !!!
    # used by the asymptotic variance functions
    # ie. assumes that covariates passed to the logistic and survival submodels are the same
    x.col.ix <- logistic.x.col.ix
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    param.vec.logistic <- paste(colnames(data)[logistic.x.col.ix], '.cc.inc.cc.prev', sep = '')
    param.vec.survival <- paste(colnames(data)[survival.x.col.ix], '.cc.surv.gamma', sep = '')

    if(const.bh){ # edit the covariate names so that they can be "grepped" by the fitting functions
        param.vec.fitting.names <-  c('lb0.cc.inc', 'nu0.cc.prev', param.vec.logistic, 'k2.scale.baseline.haz', param.vec.survival)
        param.vec.fitting.names.nice <- c('alpha', 'nu', logistic.vars, 'k2.scale', survival.vars)
        n.bh.params <- 1
    }else{ # Weibull
        param.vec.fitting.names <-  c('lb0.cc.inc', 'nu0.cc.prev', param.vec.logistic, 'k1.shape.baseline.haz', 'k2.scale.baseline.haz', param.vec.survival)
        param.vec.fitting.names.nice <- c('alpha', 'nu', logistic.vars, 'k1.shape', 'k2.scale', survival.vars)
        n.bh.params <- 2
    }
    if(is.null(param.vec.fitting.start.values)){
        param.vec.fitting.start.values <- rep.int(0.1, length(param.vec.fitting.names))
    }

    if(n2 == 0){ # no survival part and no nu*
        ix <- c(1,3:(3+length(logistic.x.col.ix)-1))
        param.vec.fitting.start.values <- param.vec.fitting.start.values[ix]
        param.vec.fitting.names        <- param.vec.fitting.names[ix]
        param.vec.fitting.names.nice   <- param.vec.fitting.names.nice[ix]
    }

    if(is.null(xi.a)){
        ix <- data$case.status == 2
        xi.a <- max(data$backward.time[ix])*1.01
        rm(ix)
    }

    # update the data control object
    data.control      <- list(n0 = n0,
                              n1 = n1,
                              n2 = n2,
                              ix.0 = ix.0,
                              ix.1 = ix.1,
                              ix.2 = ix.2,
                              N = n0 + n1 + n2,
                              logistic.vars = logistic.vars,
                              survival.vars = survival.vars,
                              logistic.x.col.ix = logistic.x.col.ix,
                              survival.x.col.ix = survival.x.col.ix,
                              x.col.ix = x.col.ix, # used by the asy variance
                              xi.a = xi.a,
                              param.vec.fitting.start.values = param.vec.fitting.start.values,
                              param.vec.fitting.names = param.vec.fitting.names,
                              param.vec.fitting.names.nice = param.vec.fitting.names.nice,
                              print.ok = print.ok,
                              const.bh = const.bh,
                              n.bh.params = n.bh.params)

    return(data.control)
}







###########################################################
# GLM.PREV.CC
#
# data should have at least the following columns:
#   case.status    0 = control, 1 = incident case, 2 = prevalent case
#   backward.time  numeric ranging from >= 0
#   all covariates listed in logistic.vars and survival.vars
# param.vec.fitting.start.values expects start values for: alpha, nu, betas (logistic params), baseline hazard params, gammas (survival parameters)

glm.prev.cc <- function(data,
                        logistic.vars,
                        survival.vars,
                        const.bh = F,   # const.bh = constant baseline hazard, 1 parameter, if FALSE = Weibull, 2 parameters
                        se.type = NULL, # 'asy', 'jack', 'boot', or NULL (no SEs)
                        param.vec.fitting.start.values = NULL,
                        optim.method = 'L-BFGS-B', #method = 'BFGS',
                        xi.a = NULL,
                        print.ok = F, # not currently used, can be used to print intermediate results
                        boot = 500){  # this is only used if se.type == 'boot'

    # if starting values not provided, use logistic model values for the logistic part, and best-guess constants for the survival model parameters
    if(is.null(param.vec.fitting.start.values)){
        formula.inc <- as.formula(paste('case.status ~ ', paste(logistic.vars, collapse = ' + ')))
        m.i         <- glm(formula.inc, data = data, subset = case.status < 2, family = binomial())
        ests.inc    <- summary(m.i)$coef[-1, 1]
        param.vec.fitting.start.values = c(-2, -4, ests.inc, rep(1, ifelse(const.bh, 1, 2)), rep(-0.1, length(survival.vars)))
    }

    pc <- set.program.control(data = data,
                              logistic.vars = logistic.vars,
                              survival.vars = survival.vars,
                              const.bh = const.bh,
                              param.vec.fitting.start.values = param.vec.fitting.start.values,
                              xi.a = xi.a,
                              print.ok = print.ok)

    # performs minimization of the log likelihood, if not converging, try tweaking the start values of this function
    optimx.out <- try(optimx(pc$param.vec.fitting.start.values, fn = constrained.likelihood, hessian = T,
                             method = optim.method,
                             lower = c(rep(-Inf, length(pc$logistic.vars) + 2), rep(0.1, pc$n.bh.params), rep(-Inf, length(pc$survival.vars))), # just constrain the constant/Weibull baseline hazard parameters
                             upper = rep(Inf, length(pc$logistic.vars) + 2 + pc$n.bh.params + length(pc$survival.vars)),
                             control = list(fnscale = 1),
                             data = data,
                             pc = pc)) # BFGS, Nelder-Mead, CG

    if(class(optimx.out) == 'try-error'){ # optimization failed
        print('optimx did not converge')
        out <- NULL
    }else if(optimx.out$convcode == 0){ # successfully converged
        ests <- matrix(coef(optimx.out), ncol = 1)
        rownames(ests) <- pc$param.vec.fitting.names.nice
        loglik <- -constrained.likelihood(ests, data, pc, print.loglik = T)
        out <- list(ests = ests, loglik = loglik)

        # now get the SEs of the estimates
        if(is.null(se.type)){ # no SEs are to be estimated
            out <- c(out, list(ses = rep(NA, length(out$ests)))) # ests, loglik, ses = NA
        }else if(se.type == 'jack'){ # jackknife
            m.p.jack <- glm.prev.cc.jack(data = data, logistic.vars = logistic.vars, survival.vars = survival.vars, const.bh = const.bh, param.vec.fitting.start.values =  as.numeric(ests),
                                         optim.method = optim.method, xi.a = xi.a)
            out <- c(out, m.p.jack) # ests, loglik, ses, i.succ, conv.prop
        }else if(se.type == 'asy'){ # asymptotic
            inner    <- get.avar(data = data, pc = pc, param.ests = ests)
            hess <- matrix(unlist(attr(optimx.out, "details")[optim.method ,"nhatend"]), ncol = length(pc$param.vec.fitting.start.values))
            hess.inv <- solve(hess)
            asy.var <- hess.inv %*% inner %*% hess.inv
            out <- c(out, list(ses = sqrt(diag(asy.var)))) # ests, loglik, ses
        }else if(se.type == 'boot'){
            m.p.boot<- glm.prev.cc.boot(data = data, logistic.vars = logistic.vars, survival.vars = survival.vars, const.bh = const.bh, param.vec.fitting.start.values =  as.numeric(ests),
                                        optim.method = optim.method, xi.a = xi.a, boot = boot)
            out <- c(out, m.p.boot) # ests, loglik, ses, quant.025, quant.975, conv.prop
        }else{
            print('No such se.type option, use one of: NULL, "jack", "asy" or "boot"', q = F)
            out <- c(out, list(ses = rep(NA, length(out$ests))))
        }
    }else{ # optimization did not converge
        out <- NULL
    }
    return(out)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in:
# param.vec is a vector of parameters (with 'greppable' names)
# data = CCH dataset
# pc = list of program control variables
#
# out:
# - log likelihood
constrained.likelihood <- function(param.vec, data, pc, print.loglik = F){

    if(pc$n1 > 0 & pc$n2 > 0){
        lc1 <- -sum(log(1
                        + exp(beta.x(param.vec, data, pc, ix.type = 'all', param.type = 'incident'))
                        + exp(beta.x(param.vec, data, pc, ix.type = 'all', param.type = 'prevalent') + log(get.mu.x(param.vec, data, pc, ix.type = 'all')))))

        lc2 <- sum(beta.x(param.vec, data, pc, ix.type = 'incident', param.type = 'incident'))
        lc3 <- sum(beta.x(param.vec, data, pc, ix.type = 'prevalent', param.type = 'prevalent') + log(get.surv.a(param.vec, data, pc, ix.type = 'prevalent')))
        log.lik <- lc1 + lc2 + lc3

    }else if(pc$n1 > 0 & pc$n2 == 0){
        lc1 <- -sum(log(1 + exp(beta.x(param.vec, data, pc, ix.type = 'all', param.type = 'incident'))))
        lc2 <- sum(beta.x(param.vec, data, pc, ix.type = 'incident', param.type = 'incident'))
        log.lik <- lc1 + lc2
    }else{ # pc$n1 == 0 & pc$n2 > 0
        lc1 <- -sum(log(1 + exp(beta.x(param.vec, data, pc, ix.type = 'all', param.type = 'prevalent') + log(get.mu.x(param.vec, data, pc, ix.type = 'all')))))
        lc3 <- sum(beta.x(param.vec, data, pc, ix.type = 'prevalent', param.type = 'prevalent') + log(get.surv.a(param.vec, data, pc, ix.type = 'prevalent')))
        log.lik <- lc1 + lc3
    }
    if(print.loglik){
        print(paste('log likelihood:', round(log.lik, 2)))
    }
    return(-log.lik)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# out: T/F index of individuals of interest, length = dim(data)[1]
get.ix <- function(ix.type, data){
    if(ix.type == 'all'){
        ix <- matrix(rep(T, dim(data)[1]), ncol = 1)
    }else if(ix.type == 'incident'){
        ix <- matrix(data$case.status == 1, ncol = 1)
    }else if(ix.type == 'prevalent'){
        ix <- matrix(data$case.status == 2, ncol = 1)
    }else{
        return(NULL)
    }
    return(ix)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# XB for the logistic model
# the set of individuals is specified by ix.type (all, incident, prevalent) and the parameter type by param.type (prevalent or incident)
# out:
# matrix XB of length sum(ix) by 1
beta.x <- function(param.vec, data, pc, ix.type = NULL, param.type = NULL){

    ix <- get.ix(ix.type, data)

    if(param.type == 'incident'){           # alpha.star, betas
        logistic.params <- get.param.val('cc.inc', param.vec, pc)
    }else if(param.type == 'prevalent'){    # nu.star, betas
        logistic.params <- get.param.val('cc.prev', param.vec, pc)
    }else{
        return(NULL)
    }

    x.design.mx.logistic <- as.matrix(cbind(rep.int(1, sum(ix)), data[ix, pc$logistic.x.col.ix])) # sum(ix) x p_l
    # as above, possibly faster:
    # x.design.mx.logistic <- matrix(1, nrow = sum(ix), ncol = (length(pc$logistic.x.col.ix) + 1))
    # x.design.mx.logistic[, -1] <- data[ix, pc$logistic.x.col.ix]
    return(x.design.mx.logistic %*% logistic.params)                                              # sum(ix) x 1
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# in:
# pattern   = string grep pattern to search for the parameters of interest
# param.vec = vector of parameter values
# pc        = program control list (should contain names of the parameters listed in the same order as in param.vec)
#
# out:
# values of the relevant parameters
get.param.val <- function(pattern, param.vec, pc){
    param.val <- param.vec[grep(pattern, pc$param.vec.fitting.names)]

    if(length(param.val) == 1){
        return(param.val)
    }else{
        return(matrix(as.numeric(param.val), ncol = 1))
    }
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# out:
# (1/k2)^k1 * exp(x*gamma)
get.theta <- function(param.vec, data, pc, ix.type){
    ix     <- get.ix(ix.type, data)

    kappas <- get.param.val('baseline.haz', param.vec, pc)                # kappas: shape, scale (Weibull), scale only (const)
    gammas <- get.param.val('gamma', param.vec, pc)                       # gammas

    x.mx   <- as.matrix(data[ix, pc$survival.x.col.ix])                    # x covariates in ix (sum(ix) x p_s)

    if(pc$const.bh){
        theta  <- (1/kappas) * exp(x.mx %*% gammas)     # n2 x 1
    }else{
        theta  <- (1/kappas[2])^kappas[1] * exp(x.mx %*% gammas)     # n2 x 1
    }

    return(theta)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# out: mu(x) [length = sum(ix)]
get.mu.x <- function(param.vec, data, pc, ix.type){
    theta       <- get.theta(param.vec, data, pc, ix.type) # (1/k2)^k1 * exp(x*gamma)
    if(pc$const.bh){
        rescaling   <- 1/theta
        upper.limit <- pc$xi.a

        mu.x.part   <- matrix(NA, nrow = length(theta), ncol = 1)
        for(i in 1:length(theta)){
            mu.x.part[i] <- pgamma(upper.limit, shape = 1, rate = theta[i], lower.tail = TRUE)
        }
        mu.x        <- mu.x.part * rescaling

    }else{ # weibull
        k1          <- get.param.val('baseline.haz', param.vec, pc)[1]

        rescaling   <- gamma(1/k1)/(k1 * theta^(1/k1))
        upper.limit <- pc$xi.a^k1

        mu.x.part   <- matrix(NA, nrow = length(theta), ncol = 1)
        for(i in 1:length(theta)){
            mu.x.part[i] <- pgamma(upper.limit, shape = 1/k1, rate = theta[i], lower.tail = TRUE)
        }
        mu.x        <- mu.x.part * rescaling
    }
    return(mu.x)
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# out: S(a|x)
get.surv.a <- function(param.vec, data, pc, ix.type){
    ix     <- get.ix(ix.type, data)
    theta  <- get.theta(param.vec, data, pc, ix.type)  # (1/k2)^k1 * exp(x*gamma)
    if(pc$const.bh){
        k1     <- 1 # not estimated in this case
    }else{
        k1     <- get.param.val('baseline.haz', param.vec, pc)[1]
    }
    a      <- data$backward.time[ix]

    return(exp(-a^k1 * theta))
}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANALYTIC VARIANCE FUNCTIONS
# !!!!!!!!!!!!!!!!!!                 CAUTION                     !!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!! CURRENTLY THIS ONLY WORKS FOR 2 COVARIATES  !!!!!!!!!!!!!!!!
# !!!!!! IN THE LOGISTIC AND SURVIVAL SUBMODELS WITH A WEIBULL BASELINE HAZARD !!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ASYMPTOTIC VARIANCE ASSUMING CONSTANT HAZARD
# E(dl/dtheta dl/dtheta')
# where dl/dtheta = first derivative of the likelihood with respect to theta, with theta = (alpha, nu, betas, gammas)
# INPUT: data, parameter vector, pc (program control list)
# RETURNS: a square matrix E(dl/dtheta dl/dtheta') of size n parameters in theta x n parameters in theta

get.avar <- function(data, pc, param.ests = NULL){
    if(pc$n1 > 0 & pc$n2 > 0){
        asy.var <- get.avar.all(data, pc, param.ests)
    }else if(pc$n1 > 0 & pc$n2 == 0){
        asy.var <- get.avar.no.prev(data, pc, param.ests)
    }else if(pc$n1 == 0 & pc$n2 > 0){
        asy.var <- get.avar.no.inc(data, pc, param.ests)
    }
    return(asy.var)
}



get.avar.all <- function(data, pc, param.ests = NULL){
    if(pc$const.bh){
        out <- get.avar.all.const(data, pc, param.ests)
    }else{
        out <- get.avar.all.weibull(data, pc, param.ests)
    }
    return(out)
}


get.avar.all.const <- function(data, pc, param.ests = NULL){

    al <- param.ests[1]
    nu <- param.ests[2]
    betas  <- matrix(param.ests[3:4], ncol = 1)
    k2.scale <- matrix(param.ests[5], ncol = 1)
    gammas <- matrix(param.ests[6:7], ncol = 1)

    # INNER
    ix0  <- data[, 'case.status'] == 0
    ix1  <- data[, 'case.status'] == 1
    ix2  <- data[, 'case.status'] == 2
    x.mx <- as.matrix(data[, pc$x.col.ix])                       # vector length n0+n1+n2
    a    <- data[, 'backward.time']                              # vector length n, n0 and n1 = 0, n2 = a|x

    w1   <- exp(al + x.mx%*%betas)                               # vector length n0+n1+n2

    mu.x <- get.mu.x(param.ests, data, pc, ix.type = 'all')
    w2   <- exp(nu + x.mx%*%betas + log(mu.x))                   # k2 is free

    eta <- 1 + w1 + w2                                           # vector length n0+n1+n2

    d.dk2.scale.log.mu.x  <- jacobian(func = num.deriv.kappa.log.mu.x,  x = k2.scale, data = data, param.ests = param.ests, pc = pc)
    d.dk2.scale.log.s.a.x <- jacobian(func = num.deriv.kappa.log.s.a.x, x = k2.scale, data = data, param.ests = param.ests, pc = pc)

    d.dgammas.log.mu.x  <- jacobian(func = num.deriv.gamma.log.mu.x,  x = gammas, data = data, param.ests = param.ests, pc = pc)

    dl.dtheta <- matrix(0, nrow = pc$N, ncol = 7) # const.bh = 7      # matrix n0+n1+n2 x n params

    # alpha
    dl.dtheta[,1] <- -w1/eta

    # nu
    dl.dtheta[,2] <- -w2/eta

    # betas
    dl.dtheta[,3:4] <- v2mx(-(w1+w2)/eta) * x.mx + x.mx * v2mx(vec = (ix1|ix2))  # k2 free

    # k2.scale
    dl.dtheta[,5] <- -w2/c(eta) * d.dk2.scale.log.mu.x + d.dk2.scale.log.s.a.x  # k2 free

    # gammas
    dl.dtheta[,6:7]   <- v2mx(-w2/eta) * d.dgammas.log.mu.x + v2mx(-(a/k2.scale)*exp(x.mx%*%gammas)) * x.mx   # ix2 is not necessary, a is 0 for ix0 and ix1

    # add the pieces separately
    inner <- pc$n0*cov(dl.dtheta[ix0,]) + pc$n1*cov(dl.dtheta[ix1,]) + pc$n2*cov(dl.dtheta[ix2,]) # CORRECT

    return(inner = inner)
}





get.avar.all.weibull <- function(data, pc, param.ests = NULL){

    al <- param.ests[1]
    nu <- param.ests[2]
    betas  <- matrix(param.ests[3:4], ncol = 1)
    kappas <- matrix(param.ests[5:6], ncol = 1)
    gammas <- matrix(param.ests[7:8], ncol = 1)

    # INNER
    ix0  <- data[, 'case.status'] == 0
    ix1  <- data[, 'case.status'] == 1
    ix2  <- data[, 'case.status'] == 2
    x.mx <- as.matrix(data[, pc$x.col.ix])                       # vector length n0+n1+n2
    a    <- data[, 'backward.time']                              # vector length n, n0 and n1 = 0, n2 = a|x

    w1   <- exp(al + x.mx%*%betas)                               # vector length n0+n1+n2

    mu.x <- get.mu.x(param.ests, data, pc, ix.type = 'all')
    w2   <- exp(nu + x.mx%*%betas + log(mu.x))                   # k1, k2 are free

    eta <- 1 + w1 + w2                                           # vector length n0+n1+n2

    d.dkappas.log.mu.x  <- jacobian(func = num.deriv.kappa.log.mu.x,  x = kappas, data = data, param.ests = param.ests, pc = pc)
    d.dkappas.log.s.a.x <- jacobian(func = num.deriv.kappa.log.s.a.x, x = kappas, data = data, param.ests = param.ests, pc = pc)

    d.dgammas.log.mu.x  <- jacobian(func = num.deriv.gamma.log.mu.x,  x = gammas, data = data, param.ests = param.ests, pc = pc)

    dl.dtheta <- matrix(0, nrow = pc$N, ncol = 8)                # matrix n0+n1+n2 x n params

    # alpha
    dl.dtheta[,1] <- -w1/eta

    # nu
    dl.dtheta[,2] <- -w2/eta

    # betas
    dl.dtheta[,3:4] <- v2mx(-(w1+w2)/eta) * x.mx + x.mx * v2mx(vec = (ix1|ix2))  # k1, k2 free

    # kappas
    dl.dtheta[,5:6] <- v2mx(-w2/eta) * d.dkappas.log.mu.x + d.dkappas.log.s.a.x  # k1, k2 free

    # gammas
    dl.dtheta[,7:8]   <- v2mx(-w2/eta) * d.dgammas.log.mu.x + v2mx(-(a/kappas[2])^kappas[1]*exp(x.mx%*%gammas)) * x.mx   # ix2 is not necessary, a is 0 for ix0 and ix1

    # add the pieces separately
    inner <- pc$n0*cov(dl.dtheta[ix0,]) + pc$n1*cov(dl.dtheta[ix1,]) + pc$n2*cov(dl.dtheta[ix2,]) # CORRECT

    return(inner = inner)
}



# same for const and weibull
get.avar.no.prev <- function(data, pc, param.ests = NULL){

    # if no prev, then no nu, and no kappas or gammas
    al <- param.ests[1] # + log(pc$cch.incident.sample.size/pc$cch.control.sample.size)
    betas  <- matrix(param.ests[2:3], ncol = 1)

    # INNER
    ix0  <- data[, 'case.status'] == 0
    ix1  <- data[, 'case.status'] == 1
    x.mx <- as.matrix(data[, pc$x.col.ix])             # vector length n0+n1+n2

    # notation as in Dr. X's writeup
    w1  <- exp(al + x.mx%*%betas)                      # vector length n0+n1+n2
    eta <- 1 + w1                                      # vector length n0+n1+n2

    dl.dtheta <- matrix(0, nrow = pc$N, ncol = 3)      # matrix n0+n1+n2 x 4

    # alpha
    dl.dtheta[,1]   <- -w1/eta

    # betas
    dl.dtheta[,2:3] <- v2mx(-w1/eta) * x.mx + x.mx * v2mx(vec = ix1)

    # add the pieces separately
    inner <- pc$n0*cov(dl.dtheta[ix0,]) + pc$n1*cov(dl.dtheta[ix1,]) # CORRECT

    return(inner = inner)
}



get.avar.no.inc <- function(data, pc, param.ests = NULL){
    if(pc$const.bh){
        out <- get.avar.no.inc.const(data, pc, param.ests)
    }else{
        out <- get.avar.no.inc.weibull(data, pc, param.ests)
    }
    return(out)
}

get.avar.no.inc.const <- function(data, pc, param.ests = NULL){

    # if no incident, no alpha
    nu <- param.ests[1] # + log(pc$cch.prevalent.sample.size/pc$cch.control.sample.size)
    betas  <- matrix(param.ests[2:3], ncol = 1)
    k2.scale <- matrix(param.ests[4], ncol = 1)
    gammas <- matrix(param.ests[5:6], ncol = 1)

    # INNER
    ix0  <- data[, 'case.status'] == 0
    ix2  <- data[, 'case.status'] == 2
    x.mx <- as.matrix(data[, pc$x.col.ix])             # vector length n0+n1+n2
    a    <- data[, 'backward.time'] # vector length n, n0 and n1 = 0, n2 = a|x

    mu.x <- get.mu.x(param.ests, data, pc, ix.type = 'all')
    w2   <- exp(nu + x.mx%*%betas + log(mu.x))

    eta  <- 1 + w2                         # vector length n0+n1+n2

    d.dk2.scale.log.mu.x  <- jacobian(func = num.deriv.kappa.log.mu.x, x = k2.scale, data = data, param.ests = param.ests, pc = pc)
    d.dk2.scale.log.s.a.x <- jacobian(func = num.deriv.kappa.log.s.a.x, x = k2.scale, data = data, param.ests = param.ests, pc = pc)

    d.dgammas.log.mu.x  <- jacobian(func = num.deriv.gamma.log.mu.x, x = gammas, data = data, param.ests = param.ests, pc = pc)


    dl.dtheta <- matrix(0, nrow = pc$N, ncol = 7) # matrix n0+n1+n2 x n params

    # nu
    dl.dtheta[,1]   <- -w2/eta

    # betas
    dl.dtheta[,2:3] <- v2mx(-w2/eta) * x.mx + x.mx * v2mx(vec = ix2)

    # k2.scale
    dl.dtheta[,4] <- v2mx(-w2/eta) * d.dk2.scale.log.mu.x + d.dk2.scale.log.s.a.x

    # gammas

    dl.dtheta[,5:6] <- v2mx(-w2/eta) * d.dgammas.log.mu.x + v2mx(-(a/k2.scale)*exp(x.mx%*%gammas)) * x.mx  # ix2 is not necessary, a is 0 for ix0 and ix1

    # add the pieces separately
    inner <- pc$n0*cov(dl.dtheta[ix0,]) + pc$n2*cov(dl.dtheta[ix2,]) # CORRECT

    return(inner = inner)
}



get.avar.no.inc.weibull <- function(data, pc, param.ests = NULL){

    # if no incident, no alpha
    nu <- param.ests[1] # + log(pc$cch.prevalent.sample.size/pc$cch.control.sample.size)
    betas  <- matrix(param.ests[2:3], ncol = 1)
    kappas <- matrix(param.ests[4:5], ncol = 1)
    gammas <- matrix(param.ests[6:7], ncol = 1)

    # INNER
    ix0  <- data[, 'case.status'] == 0
    ix2  <- data[, 'case.status'] == 2
    x.mx <- as.matrix(data[, pc$x.col.ix])             # vector length n0+n1+n2
    a    <- data[, 'backward.time'] # vector length n, n0 and n1 = 0, n2 = a|x

    mu.x <- get.mu.x(param.ests, data, pc, ix.type = 'all')
    w2   <- exp(nu + x.mx%*%betas + log(mu.x))

    eta  <- 1 + w2                         # vector length n0+n1+n2

    d.dkappas.log.mu.x  <- jacobian(func = num.deriv.kappa.log.mu.x, x = kappas, data = data, param.ests = param.ests, pc = pc)
    d.dkappas.log.s.a.x <- jacobian(func = num.deriv.kappa.log.s.a.x, x = kappas, data = data, param.ests = param.ests, pc = pc)

    d.dgammas.log.mu.x  <- jacobian(func = num.deriv.gamma.log.mu.x, x = gammas, data = data, param.ests = param.ests, pc = pc)


    dl.dtheta <- matrix(0, nrow = pc$N, ncol = 7) # matrix n0+n1+n2 x n params

    # nu
    dl.dtheta[,1]   <- -w2/eta

    # betas
    dl.dtheta[,2:3] <- v2mx(-w2/eta) * x.mx + x.mx * v2mx(vec = ix2)

    # kappas
    dl.dtheta[,4:5] <- v2mx(-w2/eta) * d.dkappas.log.mu.x + d.dkappas.log.s.a.x

    # gammas
    # dl.dtheta[,6:7] <- x.mx*v2mx(w2/eta) - x.mx * v2mx(a*exp(x.mx%*%gammas))   # assumes k1 = k2 = 1, ix2 is not necessary, a is 0 for ix0 and ix1

    dl.dtheta[,6:7] <- v2mx(-w2/eta) * d.dgammas.log.mu.x + v2mx(-(a/kappas[2])^kappas[1]*exp(x.mx%*%gammas)) * x.mx  # ix2 is not necessary, a is 0 for ix0 and ix1

    # add the pieces separately
    inner <- pc$n0*cov(dl.dtheta[ix0,]) + pc$n2*cov(dl.dtheta[ix2,]) # CORRECT

    return(inner = inner)
}





# for kappas (Weibull) or k2.scale (const)
# only calculated for prevalent, the rest are 0 because backward time is 0 for all except prevalent cases
num.deriv.kappa.log.s.a.x <- function(x, data, param.ests, pc){# }, ix.type = 'prevalent'){
    gammas <- get.param.val('gamma', param.ests, pc)
    x.mx   <- as.matrix(data[, pc$x.col.ix])
    a      <- data$backward.time

    if(pc$const.bh){
        kappa.part <- (a/x)
    }else{
        kappa.part <- (a/x[2])^x[1]
    }
    log.s.a.x <- -kappa.part * exp(x.mx %*% gammas)

    return(log.s.a.x)
}


num.deriv.kappa.log.mu.x <- function(x, data, param.ests, pc){
    if(pc$const.bh){
        out <- num.deriv.kappa.log.mu.x.const(x, data, param.ests, pc)
    }else{
        out <- num.deriv.kappa.log.mu.x.weibull(x, data, param.ests, pc)
    }
    return(out)
}


# for k2.scale for const
num.deriv.kappa.log.mu.x.const <- function(x, data, param.ests, pc){
    gammas <- get.param.val('gamma', param.ests, pc)
    x.mx   <- as.matrix(data[, pc$x.col.ix])
    theta  <- exp(x.mx %*% gammas)/c(x)

    rescaling <- 1/theta
    upper.limit <- pc$xi.a

    mu.x.part <- matrix(NA, nrow = length(theta), ncol = 1)
    for(i in 1:length(theta)){
        mu.x.part[i] <- pgamma(upper.limit, shape = 1, rate = theta[i], lower.tail = TRUE)
    }
    log.mu.x <- log(mu.x.part * rescaling)
    return(log.mu.x)
}



# for kappas (Weibull)
num.deriv.kappa.log.mu.x.weibull <- function(x, data, param.ests, pc){
    gammas <- get.param.val('gamma', param.ests, pc)
    x.mx   <- as.matrix(data[, pc$x.col.ix])
    theta  <- (1/x[2])^x[1] * exp(x.mx %*% gammas)

    rescaling <- gamma(1/x[1])/(x[1] * theta^(1/x[1]))
    upper.limit <- pc$xi.a^x[1]

    mu.x.part <- matrix(NA, nrow = length(theta), ncol = 1)
    for(i in 1:length(theta)){
        mu.x.part[i] <- pgamma(upper.limit, shape = 1/x[1], rate = theta[i], lower.tail = TRUE)
    }
    log.mu.x <- log(mu.x.part * rescaling)
    return(log.mu.x)
}



num.deriv.gamma.log.mu.x <- function(x, data, param.ests, pc){
    if(pc$const.bh){
        out <- num.deriv.gamma.log.mu.x.const(x, data, param.ests, pc)
    }else{
        out <- num.deriv.gamma.log.mu.x.weibull(x, data, param.ests, pc)
    }
    return(out)
}



num.deriv.gamma.log.mu.x.const <- function(x, data, param.ests, pc){
    k2.scale <- get.param.val('baseline.haz', param.ests, pc)

    x.mx   <- as.matrix(data[, pc$x.col.ix])
    theta  <- exp(x.mx[,1] + x.mx[,2] * c(x))/c(k2.scale)

    rescaling <- 1/theta
    upper.limit <- pc$xi.a

    mu.x.part <- matrix(NA, nrow = length(theta), ncol = 1)
    for(i in 1:length(theta)){
        mu.x.part[i] <- pgamma(upper.limit, shape = 1, rate = theta[i], lower.tail = TRUE)
    }
    log.mu.x <- log(mu.x.part * rescaling)
    return(log.mu.x)
}




num.deriv.gamma.log.mu.x.weibull <- function(x, data, param.ests, pc){
    kappas <- get.param.val('baseline.haz', param.ests, pc)
    k1     <- kappas[1]

    x.mx   <- as.matrix(data[, pc$x.col.ix])
    theta  <- (1/kappas[2])^k1 * exp(x.mx[,1] * x[1] + x.mx[,2] * x[2])

    rescaling <- gamma(1/k1)/(k1 * theta^(1/k1))
    upper.limit <- pc$xi.a^k1

    mu.x.part <- matrix(NA, nrow = length(theta), ncol = 1)
    for(i in 1:length(theta)){
        mu.x.part[i] <- pgamma(upper.limit, shape = 1/k1, rate = theta[i], lower.tail = TRUE)
    }
    log.mu.x <- log(mu.x.part * rescaling)
    return(log.mu.x)
}




# vector to matrix
v2mx <- function(vec, n.cols = 2){
    return(matrix(rep(vec, n.cols), ncol = n.cols, byrow = F))
}










# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# data should have at least the following columns:
#   case.status    0 = control, 1 = incident case, 2 = prevalent case
#   backward.time  numeric ranging from >= 0
#   all covariates listed in logistic.vars and survival.vars
#   param.vec.fitting.start.values expects start values for: alpha, nu, betas (logistic params), baseline hazard parameters, gammas (survival parameters)
glm.prev.cc.boot <- function(data,
                             logistic.vars,
                             survival.vars,
                             const.bh = F,
                             param.vec.fitting.start.values = NULL,
                             optim.method = 'L-BFGS-B',
                             xi.a = NULL,
                             print.ok = F,
                             boot = 500){

    if(is.null(param.vec.fitting.start.values)){
        param.vec.fitting.start.values <- glm.prev.cc(data = data, logistic.vars = logistic.vars, survival.vars = survival.vars, const.bh = const.bh)$ests
    }

    pc <- set.program.control(data = data,
                              logistic.vars = logistic.vars,
                              survival.vars = survival.vars,
                              const.bh = const.bh,
                              param.vec.fitting.start.values = param.vec.fitting.start.values,
                              xi.a = xi.a,
                              print.ok = print.ok)

    est.mx.boot <- matrix(NA, nrow = length(pc$param.vec.fitting.names), ncol = boot)
    rownames(est.mx.boot) <- pc$param.vec.fitting.names.nice

    i.iter <- 0
    i.succ <- 0
    while(i.succ < boot & i.iter <= (boot*3)){

        data.boot <- get.boot.sample(data, pc)

        optimx.out <- try(optimx(pc$param.vec.fitting.start.values,
                                 fn = constrained.likelihood,
                                 method = optim.method,
                                 lower = c(rep(-Inf, length(pc$logistic.vars) + 2), rep(0.1, pc$n.bh.params), rep(-Inf, length(pc$survival.vars))), # just constrain the constant/Weibull baseline hazard parameters
                                 upper = rep(Inf, length(pc$logistic.vars) + 2 + pc$n.bh.params + length(pc$survival.vars)),
                                 control = list(fnscale = 1),
                                 data = data.boot,
                                 pc = pc))

        if(class(optimx.out) == 'try-error'){
            # do nothing, move on. this iteration will be counted as unsuccessful.
        }else if(optimx.out$convcode == 0){ # successful iteration
            i.succ <- i.succ + 1
            est.mx.boot[,i.succ] <- coef(optimx.out)
        }
        i.iter <- i.iter + 1
    }
    ses       <- apply(est.mx.boot, 1, sd)
    quant.025 <- apply(est.mx.boot, 1, quantile, 0.025)
    quant.975 <- apply(est.mx.boot, 1, quantile, 0.975)

    out <- list(ses = ses, quant.025 = quant.025, quant.975 = quant.975, conv.prop = i.succ/i.iter)
    return(out)
}



glm.inc.only.boot <- function(data,
                              logistic.vars,
                              boot = 500){

    pc <- set.program.control(data = data, logistic.vars = logistic.vars)


    est.mx.boot <- matrix(NA, nrow = length(logistic.vars), ncol = boot)
    rownames(est.mx.boot) <- logistic.vars

    curr.model <- as.formula(paste('case.status ~ ', paste(logistic.vars, collapse = ' + ')))
    for(i in 1:boot){

        data.boot <- get.boot.sample(data, pc)

        m.i <- glm(curr.model, data = data.boot, subset = case.status < 2, family = binomial())
        est.mx.boot[,i] <- summary(m.i)$coef[-1, 1]
    }
    ses <- apply(est.mx.boot, 1, sd)

    out <- list(est.mx.boot = est.mx.boot, ses = ses)
    return(out)
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
get.boot.sample <- function(data, pc){
    full.ix <- 1:pc$N
    ix.0.boot <- sample(x = full.ix[pc$ix.0], size = pc$n0, replace = T)
    ix.1.boot <- sample(x = full.ix[pc$ix.1], size = pc$n1, replace = T)
    ix.2.boot <- sample(x = full.ix[pc$ix.2], size = pc$n2, replace = T)

    data.boot <- rbind(data[ix.0.boot,], data[ix.1.boot,], data[ix.2.boot,])
    return(data.boot)
}











######################################################################
# JACKKNIFE
######################################################################


# data should have at least the following columns:
#   case.status    0 = control, 1 = incident case, 2 = prevalent case
#   backward.time  numeric ranging from >= 0
#   all covariates listed in logistic.vars and survival.vars
#   param.vec.fitting.start.values expects start values for: alpha, nu, betas (logistic params), baseline hazard parameters, gammas (survival parameters)
glm.prev.cc.jack <- function(data,
                             logistic.vars,
                             survival.vars,
                             const.bh = F,
                             param.vec.fitting.start.values = NULL,
                             optim.method = 'L-BFGS-B',
                             xi.a = NULL,
                             print.ok = F){


    if(is.null(param.vec.fitting.start.values)){
        param.vec.fitting.start.values <- glm.prev.cc(data = data, logistic.vars = logistic.vars, survival.vars = survival.vars, const.bh = const.bh)$ests
    }

    pc <- set.program.control(data = data,
                              logistic.vars = logistic.vars,
                              survival.vars = survival.vars,
                              const.bh = const.bh,
                              param.vec.fitting.start.values = param.vec.fitting.start.values,
                              xi.a = xi.a,
                              print.ok = print.ok)

    est.mx.jack <- matrix(NA, nrow = length(pc$param.vec.fitting.names), ncol = pc$N)
    rownames(est.mx.jack) <- pc$param.vec.fitting.names.nice

    i.iter <- 0
    i.succ <- 0
    for(i in 1:pc$N){

        data.jack <- data[-i, ]

        optimx.out <- try(optimx(pc$param.vec.fitting.start.values,
                                 fn = constrained.likelihood,
                                 method = optim.method,
                                 lower = c(rep(-Inf, length(pc$logistic.vars) + 2), rep(0.1, pc$n.bh.params), rep(-Inf, length(pc$survival.vars))), # just constrain the constant/Weibull baseline hazard parameters
                                 upper = rep(Inf, length(pc$logistic.vars) + 2 + pc$n.bh.params + length(pc$survival.vars)),
                                 control = list(fnscale = 1),
                                 data = data.jack,
                                 pc = pc))

        if(class(optimx.out) == 'try-error'){
            # do nothing, move on. this iteration will be counted as unsuccessful.
        }else if(optimx.out$convcode == 0){ # successful iteration
            i.succ <- i.succ + 1
            est.mx.jack[,i.succ] <- coef(optimx.out)
        }
        i.iter <- i.iter + 1
    }


    print(paste(i.succ, 'of', pc$N, 'attampted jackknife iterations were successful.'))
    N <- i.succ

    est.mx.jack <- est.mx.jack[, 1:N]
    est.mx.full <- matrix(rep(param.vec.fitting.start.values, N), nrow = length(pc$param.vec.fitting.names), ncol = N)
    t.pseudo.mx <- N*est.mx.full - (N-1)*est.mx.jack

    t.pseudo.mn <- rowMeans(t.pseudo.mx)
    t.pseudo.mn.mx <- matrix(rep(t.pseudo.mn), nrow = length(pc$param.vec.fitting.names), ncol = N)
    ses <- sqrt(rowSums((t.pseudo.mx - t.pseudo.mn.mx)^2)/(N*(N-1)))

    out <- list(ses = ses, i.succ = i.succ, conv.prop = i.succ/i.iter)
    return(out)
}




