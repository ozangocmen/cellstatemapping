#     This file is part of FAiR, a program to conduct Factor Analysis in R
#     Copyright 2008 Benjamin King Goodrich

#     https://github.com/cran/FAiR/blob/master/R/FAiR.R

#     FAiR is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     FAiR is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with FAiR.  If not, see <http://www.gnu.org/licenses/>.

## Factanal() estimates all factor analysis models

Factanal <-
  function(manifest, restrictions, scores = "none", seeds = 12345, 
           lower = sqrt(.Machine$double.eps), analytic = TRUE, 
           reject = TRUE, NelderMead = TRUE, impatient = FALSE, ...) {
   
     ## Arguments
    #  manifest: an object that inherits from manifest class; see make_manifest
    #  restrictions: an object that inherits from restrictions class; see make_restrictions
    #  scores: character indicating whether / how to calculate factor scores
    #  seeds: PRNG seeds for the unif.seed and int.seed arguments of genoud()
    #  lower: lower bound on uniquenesses in factanal()-style EFA models and the lower bound
    #         on eigenvalues when determining whether a matrix is computationally posdef
    #  analytic: use analytic gradients, etc.?
    #  reject: reject starting values that fail one or more constraints?
    # NelderMead: use method = "Nelder-Mead" in a call to optim() to polish the solution?
    # impatient: logical indicating whether to skip the slow parts
    #    ...: arguments that get passed to genoud()
    
    ## Preliminaries
    
    if(!is(manifest, "manifest")) {
      stop("'manifest' must inherit from class 'manifest'")
    }
    if(!is(restrictions, "restrictions")) {
      stop("'restrictions' must inherit from class 'restrictions'")
    }
    if(lower < 0) {
      stop("'lower' must be nonnegative")
    }
    kall <- match.call()
    
    S <- cormat(manifest)
    factors  <- restrictions@factors
    scores   <-  match.arg(scores, c("none", "regression", "Bartlett", "Thurstone",
                                     "Ledermann", "Anderson-Rubin", "McDonald",
                                     "Krinjen", "Takeuchi", "Harman"))
    
    if(analytic) {
      safe <- FAiR_analytic_safe(restrictions)
      if(!safe) {
        analytic <- FALSE
        warning("'analytic' coerced to FALSE because analytic gradients",
                " are not possible in this case")
      }
    }
    
    ## copy restrictions to give it new memory
    restrictions_copy <- FAiR_copy_restrictions(restrictions)
    
    # shortcut when the user wants the factanal() behavior
    if(is(restrictions_copy, "restrictions.factanal") && impatient) {
      opt <- factanal(covmat = cormat(manifest), scores = "none",
                      factors = factors[1],
                      rotation = "none", control = list(lower = lower))
      opt$par <- opt$uniquenesses
      FAobject <- create_FAobject(restrictions_copy, manifest, opt, 
                                  kall, scores, lower, analytic)
      FAobject@seeds <- matrix(NA_integer_, ncol = 2,
                               dimnames = list("extraction", 
                                               c("unif.seed", "int.seed")))
      return(FAobject)
    } # else proceed to use genoud()
    
    ## Prepare to call genoud() via the well-known model.frame() trick
    mc <- match.call(expand.dots = TRUE)
    mc[[1]] <- as.name("genoud")
    mc[names(formals(Factanal))] <- NULL
    
    # Arguments for genoud() that are *logically required* by Factanal()
    mc$nvars    <- restrictions_copy@nvars
    mc$max      <- FALSE
    mc$hessian  <- FALSE
    mc$lexical  <- ifelse(is(restrictions_copy, "restrictions.factanal"), FALSE, TRUE)
    mc$Domains  <- restrictions_copy@Domains
    mc$default.domains <- NULL
    mc$data.type.int   <- FALSE
    if(any(restrictions_copy@Domains[,1] == restrictions_copy@Domains[,2])) {
      mc$boundary.enforcement <- 2
    }
    if(restrictions@model == "SEFA") mc$fn <- function(par) { # lexical fit function
      fitS4(par, restrictions_copy, manifest, lower, TRUE)
    }
    else mc$fn <- function(par) { # lexical fit function
      fitS4(par, restrictions_copy, manifest, lower, FALSE)
    }
    
    if(is(restrictions_copy, "restrictions.factanal")) mc$BFGSfn <- NULL
    else mc$BFGSfn <- function(par, helper = NA) { # scalar fit function (continuous)
      bfgs_fitS4(par, restrictions_copy, manifest, helper, lower)
    }
    mc$BFGShelp <- function(initial, done = FALSE) { # helper for BFGSfn and gr
      bfgs_helpS4(initial, restrictions_copy, manifest, done, lower)
    }
    if(restrictions_copy@discrepancy == "MLE") {
      mc$gr <- function(par, helper) { # gradient for ML
        gr_fitS4(par, restrictions_copy, manifest, helper, lower)
      }
    }
    else if(FAiR_is.QD(restrictions_copy)) {
      mc$gr <- function(par, helper) { # gradient for QD
        FAiR_gradient_QD(par, restrictions_copy, manifest, 
                         helper, lower)
      }
    }
    else          mc$gr <- NULL # numeric gradient
    
    if(!analytic) mc$gr <- NULL # overwrite whatever was previous
    
    # Workaround to get replicatability from genoud()
    if(any(!is.null(mc$unif.seed) | !is.null(mc$int.seed))) {
      warning("Use the seeds argument to Factanal() instead of the unif.seed",
              "and int.seed arguments to genoud(). Using 12345 as the seed.")
    }
    if(!is.null(seeds)) {
      mc$unif.seed <- seeds[1]
      mc$int.seed  <- if(length(seeds) == 1) seeds else seeds[2]
    }
    else { # use genoud() defaults
      mc$unif.seed <- 812821
      mc$int.seed  <- 53058
    }
    
    # "Default" arguments for genoud() that can be superceded if explicitly specified
    if(is.null(mc$boundary.enforcement)) mc$boundary.enforcement <- 1
    if(is.null(mc$pop.size))             mc$pop.size <- 1000
    if(is.null(mc$MemoryMatrix))         mc$MemoryMatrix <- FALSE
    if(is.null(mc$print.level))          mc$print.level <- 1
    if(is.null(mc$P9mix))                mc$P9mix <- 1
    if(is.null(mc$BFGSburnin))           mc$BFGSburnin <- -1
    if(is.null(mc$max.generations))      mc$max.generations <- 1000
    if(is.null(mc$project.path))         mc$project.path <- paste(tempfile(), 
                                                                  "Factanal.txt", sep = "")
    
    # Deal with starting values
    if(FAiR_is.FA(SV <- eval(mc$starting.values))) { # use old result to start
      par <- c(coef(SV), log(SV@scale))
      if(is(restrictions_copy, "restrictions.factanal")) {
        par <- uniquenesses(SV)
      }
      else if(is(restrictions_copy, "restrictions.independent")) {
        par <- log(SV@scale)[restrictions_copy@free]
      }
      else if(is(restrictions_copy, "restrictions.orthonormal")) {
        Phi <- diag(factors[1])
        par <- c(Phi[lower.tri(Phi)], par)
        par <- par[restrictions_copy@free]
      }
      else if(is(restrictions_copy, "restrictions.2ndorder")) {
        if(is(SV, "FA.2ndorder")) {
          Xi <- cormat(SV, level = 2)
          Delta <- loadings(SV, level = 2)
          par <- c(Xi[lower.tri(Xi)], Delta, par)
        }
        else {
          Xi <- diag(factors[2])
          par <- c(Xi[lower.tri(Xi)], rep(0, prod(factors)), par)
        }
        par <- par[restrictions_copy@free]
      }
      else if(is(restrictions_copy, "restrictions.general")) {
        if(is(SV, "FA.general")) {
          par <- c(loadings(SV, level = 2), par)
        }
        else    par <- c(rep(0, factors[1]), par)
        par  <- par[restrictions_copy@free]
      }
      else if(is(restrictions_copy, "restrictions.1storder")) {
        Phi <- cormat(SV)
        par <- c(Phi[lower.tri(Phi)], par)
        par <- par[restrictions_copy@free]
      }
      else if(length(SV@optimization$extraction$par) == mc$nvars) {
        par <- SV@optimization$extraction$par
      }
      else {
        stop("it does not seem possible to extract starting values from ",
             "the object passed to 'starting.values'\n",
             "Please supply alternate starting values or leave ",
             "'starting.values' unspecified")
      }
      mc$starting.values <- par
    }
    else if(is(restrictions_copy, "restrictions.independent")) { # independence model
      if(is.null(SV)) {
        mc$starting.values <- log(manifest@sds)[restrictions_copy@free]
      }
    }
    
    if(is.null(mc$starting.values)) { # usual case of no starting values
      cat("Generating good starting values, have patience ...\n")
      flush.console()
      pop.size <- eval(mc$pop.size)
      if(impatient) {
        efa <- factanal(covmat = S, factors = factors[1],
                        rotation = "none")
        start <- 1 - efa$uniquenesses
      }
      else { # Use approximate PACE to get initial communality estimates
        start <- FAiR_PACE_by_RGENOUD(S, factors[1], seeds = seeds)
        start <- as.matrix(start)
      }
      
      # Call method for creating start values for all parameters
      if(!is.null(seeds)) set.seed(mc$unif.seed)
      mc$starting.values <- create_start(pop.size, start, restrictions_copy,
                                         manifest, reject)
    }
    else if(any(is.na(eval(mc$starting.values)))) mc$starting.values <- NULL
    else if(length(eval(mc$starting.values)) == ncol(S)) { # starting communalities
      pop.size <- eval(mc$pop.size)
      # Call method for creating start values for all parameters
      start <- as.matrix(eval(mc$starting.values))
      if(ncol(start) > nrow(start)) start <- t(start)
      if(!is.null(seeds)) set.seed(mc$unif.seed)
      mc$starting.values <- create_start(pop.size, start, restrictions_copy,
                                         manifest, reject)
    }
    else if(length(eval(mc$starting.values)) == mc$nvars) {
      if(any(eval(mc$starting.values) < restrictions_copy@Domains[,1])) {
        stop("some starting values are below their lower bounds")
      }
      if(any(eval(mc$starting.values) > restrictions_copy@Domains[,2])) {
        stop("some starting values are above their upper bounds")
      }
    }
    else if(!is.matrix(eval(mc$starting.values))) {
      stop("'starting.values' must either be an object of class FA, a numeric",
           " vector of length ", mc$nvars, " or a numeric matrix")
    }
    
    ## Estimate model
    opt <- eval(mc) # does all the real estimation work
    
    if(NelderMead && !is(restrictions_copy, "restrictions.factanal") &&
       !is(restrictions_copy, "restrictions.independent")) {
      help  <- bfgs_helpS4(opt$par, restrictions_copy, manifest, FALSE, lower)
      foo <- function(par) {
        bfgs_fitS4(par, restrictions_copy, manifest, help, lower)
      }
      par <- opt$par
      if(restrictions_copy@model == "SEFA") par[help$squashed] <- 0
      optNM <- optim(par = par, fn = foo, method = "Nelder-Mead")
      if(optNM$value < opt$value[help$marker]) {
        opt$par <- optNM$par
        opt$value <- fitS4(opt$par, restrictions_copy, 
                           manifest, lower, TRUE)
        opt$counts <- optNM$counts
        opt$convergence <- optNM$convergence
        opt$message <- optNM$message
        if(optNM$convergence == 10) {
          warning("Nelder-Mead solution is probably on a boundary")
          print(paste("Convergence code:", opt$convergence, 
                      "Convergence message:", opt$message))
          cat("\n")
        }
      }
      else if(optNM$value == opt$value[help$marker]) {
        cat("Nelder-Mead resulted in no improvement; ",
            "convergence presumably achieved\n")
      }
      else cat("Nelder-Mead optimizer went backwards; thus ignored\n")
      
      if(restrictions_copy@model == "SEFA") par[help$squashed] <- 0
    }
    
    # Call method to postprocess opt and bake object of correct class
    FAobject <- create_FAobject(restrictions_copy, manifest, opt, 
                                kall, scores, lower, analytic)
    FAobject@seeds <- matrix(c(mc$unif.seed, mc$int.seed), ncol = 2,
                             dimnames = list("extraction", c("unif.seed", "int.seed")))
    return(FAobject)
  }
save.image()
