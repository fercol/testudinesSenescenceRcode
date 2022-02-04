# ============================= CODE METADATA ================================ #
# AUTHOR: Fernando Colchero
# SCRIPT TYPE: Functions
# DATE CREATED: 2022-02-01
# DESCRIPTION: Functions for demo code for reproducibility of manuscript 
#              Da Silva et al.
# ============================== START CODE ================================== ## ======================================= #
# ==== MAIN BAYESIAN PGLS FUNCTIONS: ====
# ======================================= #
# Main function:
RunBayesPGLS <- function(formula, data, Sigma, estLambda = TRUE, niter = 30000,
                         burnin = 10001, thinning = 10, nsim, ncpus) {
  
  # Function to correct residuals for phylogeny:
  Dfun <- function(Cmat) {
    iCmat <- solve(Cmat, tol = .Machine$double.eps)
    svdCmat <- La.svd(iCmat)
    D <- svdCmat$u %*% diag(sqrt(svdCmat$d)) %*% t(svdCmat$v)
    return(t(D))
  }
  
  # Start time counter:
  Start <- Sys.time()
  
  # Extract response:
  n <- nrow(data)
  chForm <- as.character(formula)
  y <- data[, chForm[2]]
  
  # Create design matrix:
  X <- model.matrix(as.formula(paste("~", chForm[3])), data = data)
  p <- ncol(X)
  parnames <- colnames(X)
  
  # Correlation between predictors:
  if (p == 1) {
    corPred <- NA
  } else if (ncol(X) == 2 & any(grepl("Intercept", parnames))) {
    corPred <- NA
  } else {
    corPred <- cor(X[, -grep("Intercept", parnames)])
    pp <- ifelse(any(grepl("Intercept", parnames)), p - 1, p)
    if (is.matrix(corPred)) {
      for (ii in 1:(pp-1)) {
        for (jj in ii:pp) {
          corPred[jj, ii] <- NA
        }
      }
      corPred <- corPred[-pp, -1]
    }
  }
  
  # Run sequence to find jumps for lambda:
  if (estLambda) {
    cat("Running sequence to calculate lambda jump sd... ")
    outjumps <- RunMCMC(sim = 1, y = y, X = X, Sigma = Sigma, n = n, p = p, 
                        niter = burnin, burnin = burnin, thinning = 1,
                        estLambda = TRUE, updateJumps = TRUE)
    jumpLam <- outjumps$jump
    cat("Done.\n")
  } else {
    jumpLam <- NULL
  }
  
  # run parallel estimation:
  cat("Running MCMC...\n ")
  sfInit(parallel = TRUE, cpus = ncpus)
  
  # Upload paramDemo:
  suppressMessages(sfSource("02code/TestuSenesfunctions.R"))
  
  # libraries:
  sfLibrary(mvtnorm)
  
  # Run parallel function:
  outparal <- sfClusterApplyLB(1:ncpus, RunMCMC, y = y, X = X, Sigma = Sigma,
                               n = n, p = p, niter = niter, burnin = burnin, 
                               thinning = thinning, estLambda = estLambda, 
                               updateJumps = FALSE, jumpLam = jumpLam)
  
  # Stop application:
  sfStop()
  cat("Done.\n")
  
  # rename runs:
  cat("Calculating summary statistics... ")
  names(outparal) <- sprintf("run%s", 1:nsim)
  
  # ========================== #
  # ==== EXTRACT RESULTS: ====
  # ========================== #
  # index for included:
  keep <- seq(burnin, niter, thinning)
  nkeep <- length(keep)
  
  # Parameter, likelihood and posterior tables:
  parsmat <- outparal[[1]]$pars[keep, ]
  lpmat <- outparal[[1]]$likepost[keep, ]
  
  for (isim in 1:nsim) {
    parsmat <- rbind(parsmat, outparal[[isim]]$pars[keep, ])
    lpmat <- rbind(lpmat, outparal[[isim]]$likepost[keep, ])
  }
  
  # Convergence statistics:
  PSRF <- CalcPSRF(outparal, keep = keep, nsim = nsim)
  
  # extract coefficients:
  coefs <- cbind(Mean = apply(parsmat, 2, mean), SD = apply(parsmat, 2, sd),
                 Lower = apply(parsmat, 2, quantile, 0.025), 
                 Upper = apply(parsmat, 2, quantile, 0.975))
  if (any(grepl("Intercept", parnames))) {
    idpvals <- 2:p
  } else {
    idpvals <- 1:p
  }
  pvals <- sapply(idpvals, function(pp) {
    pn <- pnorm(0, mean = coefs[pp, "Mean"], sd = coefs[pp, "SD"])
    if (pn > 0.5) {
      pn <- 2 * (1 - pn)
    } else {
      pn <- 2 * pn
    }
    if (pn < 0.0001) {
      pn <- "< 0.001"
    } else {
      pn <- format(round(pn, 3), scientific = FALSE)
    }
    return(pn)
  })
  
  allpvals <- rep(NA, nrow(coefs))
  allpvals[idpvals] <- pvals
  coefs <- data.frame(coefs, zeroCoverage = allpvals, 
                      Rhat = PSRF[, "Rhat"])
  
  # =================================== #
  # ==== EXTRACT PPOINT ESTIMATES: ====
  # =================================== #
  # Residual analysis:
  betahat <- coefs[1:p, 1]
  sighat <- coefs["sigma", 1]
  if (estLambda) {
    lambdahat <- coefs["lambda", 1]
  } else {
    lambdahat <- 0
  }
  SigHat <- Sigma * lambdahat
  diagSig <- rep(1, nrow(Sigma))
  diag(SigHat) <- diagSig
  detSigHat <- det(SigHat)
  SigInvHat <- solve(SigHat)
  betahat <- solve(t(X) %*% SigInvHat %*% X) %*% t(X) %*% SigInvHat %*% y
  muhat <- X %*% betahat
  reshat <- y - muhat
  Vhat <- Dfun(SigHat)
  phreshat <- Vhat %*% reshat
  fitted <- data.frame(mu = muhat, resid = phreshat)
  
  # =================================== #
  # ==== INFLUENTIAL OBSERVATIONS: ====
  # =================================== #
  # leverages:
  hhat <- diag(X %*% solve(t(X) %*% SigInvHat %*% X) %*% t(X) %*% SigInvHat)
  names(hhat) <- rownames(SigInvHat)
  
  # Estimated sigma:
  shat <- c((t(y - muhat) %*% SigInvHat %*% (y - muhat)) / (n - p))
  
  # standardized residuals:
  rhat <- reshat / sqrt(shat * (1 - hhat))
  
  # Studentized residuals:
  ti <- rhat * sqrt((n - p - 1) / (n - p - rhat^2))
  
  # Find potential outliers:
  alpha <- 0.05
  pv <- 2 * (1 - pt(abs(ti), df = n - p -1))
  
  # Potential influential obs:
  idinfl <- which(hhat > 2 * p / n | pv < alpha)
  if (length(idinfl) > 0) {
    inflObs <- cbind(leverages = hhat[idinfl], outliers = pv[idinfl])
  } else {
    inflObs <- NA
  }
  
  # DIC:
  likehat <- multiNorm(x = y, mean = muhat, invSig = 1 / sighat * SigInvHat, 
                       detSig = detSigHat)
  
  likeMean <- mean(lpmat[, "Likelihood"])
  pDIC <- 2 * (likehat - likeMean)
  DIC <- c(likeHat = likehat, likeMean = likeMean, 
           pDIC = pDIC, DIC = - 2 * likehat + 2 * pDIC)
  
  # Formula:
  form <- sprintf("%s ~ %s", chForm[2], chForm[3])
  cat("Done.\n")
  # ======================== #
  # ==== GATHER OUTPUT: ====
  # ======================== #
  finout <- list(coefficients = coefs, DIC = DIC, pars = parsmat, 
                 likepost = lpmat, runs = outparal, fitted = fitted,
                 settings = c(niter = niter, burnin = burnin, 
                              thinning = thinning, nsim = nsim),
                 form = form, estLambda = estLambda, corPred = corPred,
                 potInflObs = inflObs)
  class(finout) <- "BayesPGLS"
  End <- Sys.time()
  compTime <- End - Start
  timeUnits <- units(compTime)
  cat(sprintf("Total computing time %s %s\n", round(as.numeric(compTime), 2), 
              timeUnits))
  return(finout)
}

# Plotting:
plot.BayesPGLS <- function(x, plot.type = "traces") {
  op <- par(no.readonly = TRUE)
  if (plot.type == "traces") {
    p <- ncol(x$pars)
    parnames <- colnames(x$pars)
    nsim <- length(x$runs)
    niter <- x$settings["niter"]
    thinning <- x$settings["thinning"]
    plincl <- seq(1, niter, thinning)
    xlim <- range(plincl)
    if (nsim <= 9) {
      cols <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
                '#FFFF33', '#A65628', '#F781BF', '#999999')[1:nsim]
    } else {
      cols <- grey((nsim:1 + 1) /  nsim)
    }
    par(mfrow = c(ceiling(p / 2), 2))
    for (pp in 1:p) {
      ylim <- range(sapply(1:nsim, function(ss) {
        yran <- range(x$runs[[ss]]$pars[, pp])
        return(yran)
      }))
      plot(xlim, ylim, col = NA, xlab = "", ylab = "", main = parnames[pp])
      for (ss in 1:nsim) {
        lines(plincl, x$runs[[ss]]$pars[plincl, pp], col = cols[ss])
      }
      if (!parnames[pp] %in% c("sigma", "lambda")) {
        abline(h = 0, col = 'grey40')
      }
    }
  } else if (plot.type == "density") {
    p <- ncol(x$pars)
    parnames <- colnames(x$pars)
    nsim <- length(x$runs)
    niter <- x$settings["niter"]
    burnin <- x$settings["burnin"]
    thinning <- x$settings["thinning"]
    plincl <- seq(burnin, niter, thinning)
    xlim <- range(plincl)
    if (nsim <= 9) {
      cols <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
                '#FFFF33', '#A65628', '#F781BF', '#999999')[1:nsim]
    } else {
      cols <- grey((nsim:1 + 1) /  nsim)
    }
    par(mfrow = c(ceiling(p / 2), 2))
    for (pp in 1:p) {
      densx <- density(x$pars[, pp])
      idcis <- which(densx$x >= x$coefficients[pp, "Lower"] & 
                       densx$x <= x$coefficients[pp, "Upper"])
      densxci <- list(x = densx$x[idcis], y = densx$y[idcis])
      ncis <- length(idcis)
      idmean <- which(abs(densx$x - x$coefficients[pp, "Mean"]) ==
                        min(abs(densx$x - x$coefficients[pp, "Mean"])))[1]
      ymean <- c(0, densx$y[idmean])
      xmean <- rep(densx$x[idmean], 2)
      plot(range(densx$x), c(0, max(densx$y)), 
           col = NA, xlab = "", ylab = "", main = parnames[pp])
      polygon(c(densxci$x, rev(densxci$x)), c(rep(0, ncis), rev(densxci$y)),
              col = adjustcolor('dark red', alpha.f = 0.25), border = NA)
      lines(xmean, ymean, col = 'white')
      lines(densx$x, densx$y, col = 'dark red')
      if (!parnames[pp] %in% c("sigma", "lambda")) {
        abline(v = 0, col = 'grey40')
      }
    }
  } else if (plot.type == "diagnostics") {
    par(mfrow = c(2, 1))
    # Residuals vs fitted values:
    plot(x$fitted$mu, x$fitted$resid, xlab = "Fitted values", 
         ylab = "Residuals")
    abline(h = 0)
    
    # QQ-norm plots on residuals:
    res <- x$fitted$resid / sd(x$fitted$resid)
    qqnorm(res, xlab = "Theoretical quantiles", ylab = "Residuals")
    qqline(res)
    # lines(c(-10, 10), c(-10, 10), lty = 2)
  }
  par(op)
}

# Summary and print:
summary.BayesPGLS <- function(x, ...) {
  args <- list(...)
  if ("digist" %in% names(args)) {
    digits <- args$digits
  } else {
    digits <- 3
  }
  cat("Model:\n")
  cat(x$form)
  
  if (!is.na(x$corPred[1])) {
    cat("\n\nCorrelations between predictors:\n")
    print(x$corPred, digits = digits)
  }
  
  if (!is.na(x$potInflObs[1])) {
    cat("\n\nPotential influential obs.:\n")
    print(x$potInflObs, digits = digits)
  }
  
  cat("\n\nCoefficients: \n")
  print(x$coefficients, digits = digits)
  
  cat(sprintf("\nDIC = %s\n", round(x$DIC["DIC"], digits = digits)))
}

print.BayesPGLS <- function(x) {
  digits <- 3
  cat("Model:\n")
  cat(as.character(x$form))
  
  cat("\n\nCoefficients: \n")
  print(x$coefficients, digits = digits)
  
  cat(sprintf("\nDIC = %s\n", round(x$DIC["DIC"], digits = digits)))
}

# Function to identify and evaluate potential influential obs.:
TestInflObs <- function(formula, responseData, predictorData, predictors, 
                        estLambda = TRUE, niter = 30000, burnin = 10001, 
                        thinning = 10, nsim, ncpus, KLthresh = 0.8) {
  
  # Iteratively remove infl. obs and check fits:
  inflAn <- list()
  
  # Run full regression to identify influential obs.:
  regrDat <- PrepRegrData(responseData = responseData, 
                          predictorData = predictorData, 
                          predictors = predictors, phyloAll = phyloAll)
  
  inflAn$All <- RunBayesPGLS(formula = formula, 
                             data = regrDat$data, Sigma = regrDat$Sigma, 
                             niter = niter, burnin = burnin, 
                             thinning = thinning, nsim = nsim, ncpus = ncpus)
  
  # Extract species:
  if (!is.na(inflAn$All$potInflObs[1])) {
    sumInfl <- inflAn$All$potInflObs
    infSps <- gsub("_", " ", rownames(sumInfl))
    ninfsp <- length(infSps)
    for (sp in infSps) {
      regrDat <- PrepRegrData(responseData = responseData, predictorData = predictorData, 
                              predictors = predictors, exclSps = sp, 
                              phyloAll = phyloAll)
      
      inflAn[[sp]] <- RunBayesPGLS(formula = formula, 
                                   data = regrDat$data, Sigma = regrDat$Sigma, 
                                   nsim = 4, ncpus = 4)
    }
    
    # Parameter names:
    parNames <- rownames(inflAn[[1]]$coefficients)
    np <- length(parNames)
    
    # Calculate Kullback-Leibler discrepancies when excluding observations:
    klmat <- matrix(NA, ninfsp, np, dimnames = list(infSps, parNames))
    for (sp in infSps) {
      for (pname in parNames) {
        kl <- CalcKLc(m1 = inflAn$All$coefficients[pname, "Mean"], 
                      sd1 = inflAn$All$coefficients[pname, "SD"],
                      m2 = inflAn[[sp]]$coefficients[pname, "Mean"],
                      sd2 = inflAn[[sp]]$coefficients[pname, "SD"])
        klmat[sp, pname] <- kl["mqKl"]
      }
    }
    
    # Identify influential obs:
    idpot <- which(c(klmat) > KLthresh)
    if (length(idpot) > 0) {
      infcol <- ceiling(idpot / ninfsp)
      infrow <- idpot - (infcol - 1) * ninfsp
      potInf <- data.frame(species = infSps[infrow], coeff = parNames[infcol],
                           KL = klmat[idpot], stringsAsFactors = FALSE)
      
    } else {
      potInf <- NA
      cat("\nNo potential influential obs. found.\n")
    }
  } else {
    klmat <- NA
    potInf <- NA
    cat("\nNo potential influential obs. found.\n")
  }
  infres <- list(regrs = inflAn, kl = klmat, poteInflObs = potInf)
  class(infres) <- "potInflObs"
  return(infres)
}

summary.potInflObs <- function(x, digits = 2) {
  if (!is.data.frame(x$poteInflObs)) {
    cat("\nNo influential observations found.\n")
    
    cat("\nKullback-Leibler discrepancies table:\n")
    print(x$kl, digits = digits)
    
  } else {
    cat("\n Potential influential obs:\n")
    potInf <- x$poteInflObs
    for (xi in 1:nrow(potInf)) {
      sp <- potInf$species[xi]
      param <- potInf$coeff[xi]
      kl <- round(potInf$KL[xi], 2)
      comp <- rbind(All = x$regrs$All$coefficients[param, 1:5],
                    Excl = x$regrs[[sp]]$coefficients[param, 1:5])
      cat(sprintf("\nSps: %s\nCoef: %s\nK-L = %s\n", sp, param, kl))
      print(comp, digits = digits)
    }
  }
}
print.potInflObs <- function(x) {
  summary(x)
}

# ========================== #
# ==== MCMC FUNCTIONS: =====
# ========================== #
# Multivariate normal density:
multiNorm <- function(x, mean, invSig, detSig) {
  dens <- log(sqrt(detSig)) + 
    c(- 1 / 2 * t(x - mean) %*% invSig %*% (x - mean))
  return(dens)
}


# Truncated normal:
rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- qnorm(ru, mean, sd)
  return(rx)
}

dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# MCMC for phylogenetic regression:
RunMCMC <- function(sim, y, X, Sigma, n, p, niter, burnin, thinning,
                    estLambda = TRUE, updateJumps = TRUE, jumpLam = NULL) {
  # Starting parameters:
  sigNow <- 0.05
  if (estLambda) {
    lambdaNow <- 0.5
  } else {
    lambdaNow <- 0
    updateJumps <- FALSE
  }
  
  # Starting parameters:
  betaNow <- rep(0, p)
  names(betaNow) <- colnames(X)
  muNow <- c(X %*% betaNow)
  sigNow <- 0.05
  SigNow <- Sigma * lambdaNow
  diagSig <- rep(1, n)
  diag(SigNow) <- diagSig
  SigInvNow <- solve(SigNow)
  detSigNow <- det(SigNow)
  
  # Priors:
  betaPriorM <- rep(0, p)
  betaPriorV <- rep(100, p)
  lamPriorM <- 0.5
  lamPriorSD <- 0.5
  s1 <- 0.01
  s2 <- 0.01
  
  # Calculate initial likelyhood, prior and posterior:
  likeNow <- multiNorm(x = y, mean = muNow, invSig = 1 / sigNow * SigInvNow, 
                       detSig = detSigNow)
  postNow <- likeNow + dtnorm(x = lambdaNow, mean = lamPriorM, sd = lamPriorSD, 
                              lower = 0, upper = 1, log = TRUE)
  
  # Output matrices:
  outbeta <- matrix(NA, niter, p)
  colnames(outbeta) <- colnames(X)
  outvars <- matrix(NA, niter, 2)
  colnames(outvars) <- c("sigma", "lambda")
  outpost <- matrix(NA, niter, 2)
  colnames(outpost) <- c("Likelihood", "Posterior")
  
  # Jump update for lambda:
  if (updateJumps) {
    jumpLam <- 0.01
    updTarg <- 0.25
    updIters <- 50
    updInd <- rep(0, updIters)
    jumpVec <- jumpLam
    updSeq <- seq(updIters, niter, updIters)
    updMax <- max(updSeq)
    ucnt <- 0
  }
  
  # Run MCMC:
  for (iter in 1:niter) {
    # Jump update counter:
    if (updateJumps) {
      ucnt <- ucnt + 1
    }
    
    # ================== #
    # 1. SAMPLING BETA:
    # ================== #
    # 1.a. Sample beta's (Direct sampling through conjugate distrs.):
    v <- (t(X) %*% SigInvNow %*% y) / sigNow + betaPriorM / betaPriorV
    V <- solve((t(X) %*% SigInvNow %*% X) / sigNow + 1 / betaPriorV)
    betaNow <- t(rmvnorm(1, V %*% v, V))
    muNow <- c(X %*% betaNow)
    
    # =================== #
    # 2. SAMPLING SIGMA:
    # =================== #
    # 2.a. Sample sigma (Direct sampling; inverse gamma):
    u1 <- s1 + n / 2
    u2 <- s2 + .5 * (t(muNow - y) %*% SigInvNow %*% (muNow - y))
    sigNow <- 1 / rgamma(1, u1, u2)
    
    # update likelihood and posterior:
    likeNow <- multiNorm(x = y, mean = muNow, invSig = 1 / sigNow * SigInvNow, 
                         detSig = detSigNow)
    postNow <- likeNow + dtnorm(x = lambdaNow, mean = lamPriorM, 
                                sd = lamPriorSD, lower = 0, upper = 1, 
                                log = TRUE)
    
    # ==================== #
    # 3. SAMPLING LAMBDA:
    # ==================== #
    # 3.a. Sample lambda (Metropolis-Hastings):
    if (estLambda) {
      lambdaNew <- rtnorm(1, lambdaNow, jumpLam, lower = 0, upper = 1)
      SigNew <- Sigma * lambdaNew
      diag(SigNew) <- diagSig
      SigInvNew <- solve(SigNew)
      detSigNew <- det(SigNew)
      
      likeNew <- multiNorm(x = y, mean = muNow, invSig = 1 / sigNow * SigInvNew, 
                           detSig = detSigNew)
      postNew <- likeNew + dtnorm(x = lambdaNew, mean = lamPriorM, 
                                  sd = lamPriorSD, lower = 0, upper = 1, 
                                  log = TRUE)
      
      HastRatio <- dtnorm(x = lambdaNow, mean = lambdaNew, sd = jumpLam, 
                          lower = 0, upper = 1, log = TRUE) -
        dtnorm(x = lambdaNew, mean = lambdaNow, sd = jumpLam, lower = 0, 
               upper = 1, log = TRUE)
      r <- exp(postNew - postNow + HastRatio)
      z <- runif(1)
      if(r > z) {
        lambdaNow <- lambdaNew
        SigNow <- SigNew
        SigInvNow <- SigInvNew
        detSigNow <- detSigNew
        likeNow <- likeNew
        postNow <- postNew
        if (updateJumps & iter <= niter) {
          updInd[ucnt] <- 1
          # updInd[iter] <- 1
        }
      }      
    }
    
    # ======================== #
    # FILL-IN OUTPUT MATRICES: 
    # ======================== #
    outbeta[iter, ] <- betaNow
    outvars[iter, ] <- c(sigNow, lambdaNow)
    outpost[iter, ] <- c(likeNow, postNow)
    
    # ======================= #
    # UPDATE JUMP FOR LAMBDA:
    # ======================= #
    if (updateJumps) {
      if (iter %in% updSeq) {
        updRate <- sum(updInd) / updIters
        if (updRate == 0) updRate <- 1e-2
        jumpLam <- jumpLam * updRate / updTarg
        jumpVec <- c(jumpVec, jumpLam)
        updInd <- updInd * 0
        ucnt <- 0
        if (iter == updMax) {
          njumps <- length(jumpVec)
          idmean <- floor(njumps / 2):njumps
          jumpLam <- mean(jumpVec[idmean])
        }
      }
    }
  }
  if (updateJumps) {
    resList <- list(jump = jumpLam)
  } else {
    if (estLambda) {
      parmat <- cbind(outbeta, sigma = outvars)
    } else {
      parmat <- cbind(outbeta, sigma = outvars[, "sigma"])
    }
    resList <- list(pars = parmat, likepost = outpost,
                    jump = jumpLam, settings = c(niter = niter, burnin = burnin,
                                                 thinning = thinning))
  }
  return(resList)
}

# Function to calculate convergence statistics 
# based on Gelman et al. (2014).
CalcPSRF <- function(object, keep, nsim) {
  nthin <- length(keep)
  Means <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$pars[keep, ], 2, mean)
  }))
  Vars <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$pars[, ], 2, var)
  }))
  meanall <- apply(Means, 2, mean)
  B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
  W <- 1 / nsim * apply(Vars, 2, sum)
  Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
  Rhat <- sqrt(Varpl / W)
  Rhat[Varpl==0] <- 1
  conv <- cbind(B, W, Varpl, Rhat)
  rownames(conv) <- colnames(Means)
  return(conv)
}

# Function to calculate Kullback-Leibler discrepancies between parameter 
# posterior densities:
CalcKLc <- function(m1, sd1, m2, sd2) {
  mv <- c(m1, m2)
  sdv <- c(sd1, sd2)
  q1 <- qnorm(c(0.001, 0.999), mean = m1, sd = sd1)
  q2 <- qnorm(c(0.001, 0.999), mean = m2, sd = sd2)
  parRan <- range(c(q1, q2))
  parVec <- seq(parRan[1], parRan[2], length = 1000)
  dp <- parVec[2] - parVec[1]
  parDens <- sapply(1:2, function(pp) 
    dnorm(x = parVec, mean = mv[pp], sd = sdv[pp]))
  p1dens <- parDens[, 1]
  p2dens <- parDens[, 2]
  idp <- which(p1dens > 0 & p2dens > 0)
  kld1 <- sum(p1dens[idp] * log(p1dens[idp] / p2dens[idp]) * dp)
  kld2 <- sum(p2dens[idp] * log(p2dens[idp] / p1dens[idp]) * dp)
  qKlc1 <- (1 + (1 - exp(-2 * kld1)^(1 / 2))) / 2
  qKlc2 <- (1 + (1 - exp(-2 * kld2)^(1 / 2))) / 2
  mqKl <- (qKlc1 + qKlc2) / 2
  outList <- c(kl12 = kld1, kl21 = kld2, qkl12 = qKlc1, 
               qkl21 = qKlc2, mqKl = mqKl)
  return(outList)
}

# ==================================================== #
# ==== FUNCTIONS TO EXTRACT SUMMARY STATS. FROM BASTA:
# ==================================================== #
# Function to calculate ageing rates from BaSTA outputs in parallel:
CalcParalAR <- function(sim) {
  idseq <- (iseq[sim] + 1):iseq[sim + 1]
  nseq <- length(idseq)
  if (model == "EX") {
    arseq <- matrix(0, nseq, nsxl, 
                    dimnames = list(NULL, c("AR50", "AR20", "AR10", "AR05")))
  } else {
    arseq <- t(sapply(idseq, function(ith) {
      theta <- theMat[ith, ]
      Sxi <- CalcSurv(theta, x = x, model = model, shape = shape)
      ari <- rep(0, nsxl)
      names(ari) <- c("AR50", "AR20", "AR10", "AR05")
      for (ss in 1:nsxl) {
        idsx <- which(abs(Sxi - sxlevs[ss]) ==  min(abs(Sxi - sxlevs[ss])))[1]
        tempar <- CalcAgeingRate(theta, x[idsx], model = model, 
                                 shape = shape)
        ari[ss] <- tempar[1, "AR"]
      }
      return(ari)
    }))
  }
  return(arseq)
}


# ================================================ #
# ==== FUNCTION TO PREPARE DATA FOR ANALYSIS: ====
# ================================================ #
PrepRegrData <- function(responseData, predictorData, exclSps = NA, 
                         predictors = "ALL", phyloAll) {
  
  # tip labels (species names):
  physpAll <- as.character(phyloAll$tip.label)
  
  # Choose species to include in analysis:
  if (!is.na(exclSps[1])) {
    idincl <- which(!responseData$species %in% exclSps)
    responseData <- responseData[idincl, ]
  }
  
  # Covariates to be included:
  if (predictors[1] == "ALL") {
    predictors <- colnames(predictorData)[-1]
  }
  
  # Find if there are NAs in response data:
  naid <- sort(unique(unlist(apply(responseData[, -1], 2, function(xx) {
    id <- which(is.na(xx))
  }))))
  
  if (length(naid) > 0) {
    responseData <- responseData[-naid, ]
  }
  
  # Extract species from x:
  species <- responseData$species
  rownames(responseData) <- species
  
  # Design matrix:
  idsps <- which(predictorData$species %in% species)
  temp <- predictorData[idsps, predictors]
  rownames(temp) <- predictorData$species[idsps]
  tempPredDat <- temp[species, ]
  
  # Eliminate records with NAs:
  naid <- sort(unique(unlist(apply(tempPredDat, 2, function(xx) {
    id <- which(is.na(xx))
  }))))
  
  if (length(naid) > 0) {
    tempPredDat <- tempPredDat[-naid, ]
    responseData <- responseData[-naid, ]
  }
  
  # Extract species from x:
  spsub <- gsub(" ", "_", responseData$species)
  
  # Load BaSTA
  phylo <- drop.tip(phyloAll, which(!physpAll %in% spsub))
  
  # Match species to tree:
  physp <- as.character(phylo$tip.label)
  
  # Species orders following phylogeny:
  species <- gsub("_", " ", physp)
  nsp <- length(species)
  responseData <- responseData[species, ]
  tempPredDat <- tempPredDat[species, ]
  
  # Create data matrix:
  ndat <- data.frame(species = physp, responseData[, -1], tempPredDat)
  colnames(ndat) <- c("species", colnames(responseData)[-1], 
                      colnames(tempPredDat))
  
  # Extract Sigma matrix:
  pglsdat <- comparative.data(phy = phylo, data = ndat, vcv = T, 
                              names.col = species)
  Sigma <- pglsdat$vcv
  
  # Extract Sigma matrix:
  Sigma <- Sigma[1:nrow(Sigma), 1:ncol(Sigma)]
  Sigma <- Sigma / Sigma[1]
  nrow(Sigma)
  
  # Return list:
  retlist <- list(data = ndat, pglsDat = pglsdat, Sigma = Sigma, 
                  excluded = exclSps)
  return(retlist)
}

