fma.glmnet <-  function( x, y, focus, family, nlambda, nfolds, grouped,
                         penalty.factor = NULL, force.nlambda,
                         type.measure, rule, solnptol, qp.scale ){

    #### Note: Only jackknife is supported in this version.
    ####       x must not contain the intercept column.
    #### SJ: Define some constants
    n <- nrow(x) #### Sample size
    nbeta <- ncol(x) + 1 #### Number of coefficients, including intercept
    res <- list() # Initiate returned list
    if(is.null(penalty.factor) == TRUE){
        penalty.factor <- rep(1.0, ncol(x))
    }

    ########## Functions needed for weight estimations ##########
    #### SJ: Only useful if focus = "beta"
    if( focus == "beta" ){

        ####  SJ: fun_prob computes the cross-validation MSE
        fun_prob = function(w, B, x, y, family){

            mu = matrix(NA, nrow = nrow(x), ncol = 1)
            for(i in 1 : nrow(x)){

                beta_bar = B[, , i] %*% matrix(w, ncol = 1)

                if( family == "binomial" ){
                    mu[i] = plogis( sum( c(beta_bar) * c(1, x[i,]) ) )
                } else if( family == "poisson" ){
                    mu[i] = exp( sum( c(beta_bar) * c(1, x[i,]) ) )
                }

            }
            return(mean((y - mu) ^ 2))

        }

        ####  Sum of the weights
        fun_eq = function(w, B, x, y, family){
            z = sum(w)
            return(z)
        }


        #### VZ: New function needed for Singleton and hybrid averaging, generally a bit
        #### VZ: messier than the one above, since regressors are not the same in all models
        if("Singleton" %in% rule){
            fun_prob_sing <- function(w, B, x, y, family){
                mu <- matrix(0, nrow = nrow(x))

                for(i in 1:nrow(x)){   # Looping through folds

                    if(family == "binomial"){
                        mu[i] <- c(w[1]) * sum( c(B[, 1, i]) * c(1, x[i, 1]))
                        for(j in 2:(nbeta - 1) ){  # Looping through individual Singletons
                            mu[i] <-  mu[i] + c(w[j]) * sum( c(B[, j, i]) * c(1, x[i, j]))
                        }
                        mu[i] <- plogis(mu[i])
                    }

                    if(family == "poisson"){
                        mu[i] <- c(w[1]) * sum( c(B[, 1, i]) * c(1, x[i, 1]))
                        for(j in 2:(nbeta - 1) ){  # Looping through individual Singletons
                            mu[i] <-  mu[i] + c(w[j]) * sum( c(B[, j, i]) * c(1, x[i, j]))
                        }
                        mu[i] <- exp(mu[i])
                    }
                }
                return(mean((y - mu)^2))
            }
        }


        if("Hybrid" %in% rule){
            fun_prob_hybrid <- function(w, B, x, y, family, comb){
                mu <- matrix(0, nrow = nrow(x))

                for(i in 1:nrow(x)){   # Looping through folds

                    #### VZ: Since the betas are not all of the same dimension, I could
                    #### VZ: not find a better way of storing them than in a list

                    if(family == "binomial"){
                        for(j in 1:ncol(comb) ){  # Looping through individual models
                            sel <- which(comb[, j] %in% 1)
                            mu[i] <-  mu[i] + c(w[j] %*% cbind(1, x)[i, sel] %*% B[[i]][[j]])
                        }
                        mu[i] <- plogis(mu[i])
                    }

                    if(family == "poisson"){
                        for(j in 1:ncol(comb) ){  # Looping through individual models
                            sel <- which(comb[, j] %in% 1)
                            mu[i] <-  mu[i] + c(w[j] %*% cbind(1, x)[i, sel] %*% B[[i]][[j]])
                        }
                        mu[i] <- exp(mu[i])
                    }
                }
                return(mean((y - mu)^2))
            }

            fun_eq_hybrid = function(w, B, x, y, family, comb){
                z = sum(w)
                return(z)
            }
        }

    }

    #### ---------------------------------------------------- ####
    ####  SJ: Full model estimates
    if( family == "binomial" ){
        fulcoef <- glm.fit( x = cbind(1, x), y = y, family = binomial() )$coefficients
    } else if( family == "poisson" ){
        fulcoef <- glm.fit( x = cbind(1, x), y = y, family = poisson() )$coefficients
    }
    res$fullcoef <- fulcoef
    #### ---------------------------------------------------- ####

    #############################################################
    #### SJ: Run glmnet() to get the lambda sequence. (Any better way?)
    #### Note: lambda sequence is decreasing.
    fit_alldata = glmnet(x, y, family = family, nfolds = nfolds,
                         grouped = grouped, nlambda = nlambda,
                         type.measure = type.measure,
                         penalty.factor = penalty.factor)
    lambda_raw <- fit_alldata$lambda
    nlambda_raw <- length(lambda_raw)
    if(force.nlambda == TRUE & nlambda_raw != nlambda){
        #### I think glmnet use equal space on the log scale.
        fit_alldata = glmnet(x, y, family = family, nfolds = nfolds,
                             grouped = grouped,
                             lambda = exp(seq(max(log(lambda_raw)), min(log(lambda_raw)), length.out = nlambda)),
                             type.measure = type.measure,
                             penalty.factor = penalty.factor)
        lambda_raw <- fit_alldata$lambda
        nlambda_raw <- length(lambda_raw)
    }
    rm(nlambda)

    #### ---------------------------------------------------- ####
    #### SJ: We need to create the "correct" lambda sequence from the lambda sequence in glmnet()
    ####     This step depends on the quadrature rule that we want to use.
    ####     If we use "raw" (the discrete weight unit-simplex) or "Riemann" for integral, we do not need to introduce new lambda values.
    ####     At this stage, we do not include the full model (lambda = 0) into the lambda sequence
    lambda <- list()
    if("raw" %in% rule){
        lambda[["raw"]] <- lambda_raw
    }
    if("Riemann" %in% rule){
        lambda[["Riemann"]] <- lambda_raw
    }
    if( "Simpson(1/3)" %in% rule ){

        lambda_temp <- lambda_raw
        ####  Add the mid points to the lambda sequence
        l_seq = sort(c(lambda_temp, 0))
        S = 2 : length(l_seq)
        l_mid = c((l_seq[S] + l_seq[S - 1])/2)
        lambda_temp = numeric(length(l_mid) + length(l_seq))
        even = c(seq(from = 2, to = length(lambda_temp) - 2, by = 2), length(lambda_temp) - 1)
        lambda_temp[even] = l_mid
        lambda_temp[-even] = l_seq
        lambda_temp <- lambda_temp[-1] ####  We temporarily remove the full model.
        lambda[["Simpson(1/3)"]] <- sort(lambda_temp, decreasing = TRUE) ####  glmnet() works with decreasing sequence.

    }
    if( "Simpson(3/8)" %in% rule ){

        #### SJ: Not finished yet!!!!
        lambda_temp <- c( 0.0, sort(lambda_raw, decreasing = FALSE) )
        ####  Add the mid points to the lambda sequence
        lambdamat <- cbind( lambda_temp[-length(lambda_temp)], NA, NA, lambda_temp[-1] )
        lambdamat[, 2] <- (2.0 * lambdamat[,1] + lambdamat[,4]) / 3.0
        lambdamat[, 3] <- (lambdamat[,1] + 2.0 * lambdamat[,4]) / 3.0
        lambda_temp <- c(lambda_raw, lambdamat[,2], lambdamat[,3]) ####  We do not need the full model for cv.glmnet()
        lambda[["Simpson(3/8)"]] <- sort(lambda_temp, decreasing = TRUE) ####  glmnet() works with decreasing sequence.

    }

    #### SJ: If Simpson is considered, its lambda sequence always nests "raw".
    ####     So we only need to figure out what the "longest" lambda sequence is.
    if("Simpson(3/8)" %in% rule){
        lambda_long <- lambda[["Simpson(3/8)"]]
    } else if("Simpson(1/3)" %in% rule){
        lambda_long <- lambda[["Simpson(1/3)"]]
    } else if("raw" %in% rule){
        lambda_long <- lambda[["raw"]]
    } else if("Riemann" %in% rule){
        lambda_long <- lambda[["Riemann"]]
    }
    #### ---------------------------------------------------- ####


    #### ---------------------------------------------------- ####
    #### SJ: Using the above lambda sequence, we can start the CV step.
    if( focus == "beta" ) {

        # B: Array for beta
        B <- array(0, dim = c(ncol(x) + 1, length(lambda_long) + 1, nfolds))
        B_raw <- array(0, dim = c(ncol(x) + 1, length(lambda_raw) + 1, nfolds))

        #### VZ: Adding corresponding storage for hybrid and singletons
        if("Hybrid" %in% rule){
            betas_lasso <-  as.matrix(coef(fit_alldata, s = lambda_raw))  # Coefficients for hybrid lambda
            u <-  betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
            hybrid_comb <- unique(u, MARGIN = 2)     # Filtering unique model forms
            #### SJ: Check the full model. If it is already added, we do nothing.
            ####     If not added, we add the full model.
            if(all(colSums(hybrid_comb) != nbeta)){
                #### SJ: Full model is absent, we add the full model
                hybrid_comb <- cbind(hybrid_comb, 1)
            } else {
                #### SJ: If the full model is present, we make sure it is not added more than once
                full.index <- which(colSums(hybrid_comb) == nbeta)
                if(length(full.index) >= 2){
                    hybrid_comb <- hybrid_comb[, -full.index[-1]]
                }
            }
            B_hybrid = list()
        }
        if("Singleton" %in% rule){
            B_sing = array(0, dim = c(2, (nbeta - 1), nfolds))
        }

        ####  CV part. This part takes a long time, especially Hybrid and singleton.
        for(i in 1 : nfolds){

            x_i <- x[-i, , drop = FALSE]
            y_i <- y[-i]

            ####  glmnet()
            cvfit_long <- glmnet(x = x_i, y = y_i, family = family, lambda = lambda_long,
                                 penalty.factor = penalty.factor)
            #fit <-  glmnet(x = x_i, y = y_i, family = family, lambda = lambda)
            #fit_raw <- glmnet(x = x_i, y = y_i, family = family, lambda = lambda_raw)
            ####  full model and others.
            if( family == "binomial" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    B_hybrid[[i]] <- list()

                    for( j in 1:ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i,
                                           family = binomial())             # Fitting model
                        B_hybrid[[i]][[j]] <- coef(sel_fit)                 # Extracting coefficients
                    }
                }

                if( "Singleton" %in% rule ){
                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = binomial())
                        B_sing[, l, i] <-  matrix(coef(sing_fit))             # Extracting coefficients
                    }
                }
            } else if( family == "poisson" ){

                ####  Full model.
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    B_hybrid[[i]] <- list()

                    for( j in 1 : ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i,
                                           family = poisson())             # Fitting model
                        B_hybrid[[i]][[j]] <- coef(sel_fit)                 # Extracting coefficients

                    }
                }

                if( "Singleton" %in% rule ){

                    for( l in 1 : (nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        B_sing[, l, i] <-  matrix(coef(sing_fit))             # Extracting coefficients

                    }

                }

            }

            #### SJ: We add the full model as the last column, as lambda is a decreasing sequence
            B[, , i] = cbind( as.matrix(coef(cvfit_long)),
                              cvfit$coefficients )
            B_raw[, , i] = cbind( as.matrix(coef(cvfit_long, s = lambda_raw)),
                                  cvfit$coefficients )

        }

        #### SJ: In order to estimate the weights, we add 0 to the lambda sequence here.
        ####     Need to make sure that lambda remains a decreasing sequence!
        fit_long <- glmnet(x = x, y = y, family = family, lambda = lambda_long,
                           penalty.factor = penalty.factor)
        if("raw" %in% rule){

            coef_raw <- cbind(as.matrix(coef(fit_long, s = lambda[["raw"]])),
                              fulcoef)
            lambda[["raw"]] <- c(lambda[["raw"]], 0)
            init_w = rep( 1 / length(lambda[["raw"]]), length(lambda[["raw"]]) )
            nlopt_raw <- try( solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                    LB = rep(0, length(lambda[["raw"]])),
                                    UB = rep(1, length(lambda[["raw"]])),
                                    B = B_raw, x = x, y = y, family = family,
                                    control = list(trace = FALSE, tol = solnptol)),
                              TRUE )
            if( inherits(nlopt_raw, "try-error") == FALSE ){

                res[["raw"]] <- list()
                res[["raw"]]$w <- nlopt_raw$pars
                res[["raw"]]$convergence <- nlopt_raw$convergence #### 0 = converged
                res[["raw"]]$avecoef <- coef_raw %*% matrix(nlopt_raw$pars, ncol = 1)
                res[["raw"]]$elapsed <- nlopt_raw$elapsed
                res[["raw"]]$lambda <- lambda[["raw"]]

            } else {
                res[["raw"]] <- nlopt_raw
            }

        }
        if("Riemann" %in% rule){
            ####  SJ: Not finished yet
            lambda[["Riemann"]] <- c(lambda[["Riemann"]], 0)
        }
        if("Simpson(1/3)" %in% rule){

            coef_simp <- cbind(as.matrix(coef(fit_long, s = lambda[["Simpson(1/3)"]])),
                               fulcoef)
            lambda[["Simpson(1/3)"]] <- c(lambda[["Simpson(1/3)"]], 0)
            init_w = rep( 1 / length(lambda[["Simpson(1/3)"]]), length(lambda[["Simpson(1/3)"]]) )
            nlopt_simp <- try( solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                     LB = rep(0, length(lambda[["Simpson(1/3)"]])),
                                     UB = rep(1, length(lambda[["Simpson(1/3)"]])),
                                     B = B, x = x, y = y, family = family,
                                     control = list(trace = FALSE, tol = solnptol)),
                               TRUE )
            if( inherits(nlopt_simp, "try-error") == FALSE ){

                res[["Simpson(1/3)"]] <- list()
                res[["Simpson(1/3)"]]$v <- nlopt_simp$pars
                res[["Simpson(1/3)"]]$convergence <- nlopt_simp$convergence #### 0 = converged
                res[["Simpson(1/3)"]]$avecoef <- coef_simp %*% matrix(nlopt_simp$pars, ncol = 1)
                res[["Simpson(1/3)"]]$elapsed <- nlopt_simp$elapsed
                lambda_raw_zero <- c(lambda_raw, 0.0)
                scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
                for( i in 1 : nlambda_raw ){
                    scale <- c(scale,
                               2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                               (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
                }
                scale[length(scale)] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
                res[["Simpson(1/3)"]]$w <- nlopt_simp$pars / scale
                res[["Simpson(1/3)"]]$lambda <- lambda[["Simpson(1/3)"]]

            } else {
                res[["Simpson(1/3)"]] <- nlopt_simp
            }

        }
        if("Simpson(3/8)" %in% rule){
            ####  SJ: Not finished yet
            lambda[["Simpson(3/8)"]] <- c(lambda[["Simpson(3/8)"]], 0)
        }
        if("Singleton" %in% rule){

            init_w_sing = matrix( 1 / (nbeta - 1), nrow = (nbeta - 1) )
            nlopt_sing <- try( solnp(pars = init_w_sing, fun = fun_prob_sing, eqfun = fun_eq,
                                     eqB = c(1),
                                     LB = rep(0, (nbeta - 1)), UB = rep(1, (nbeta - 1)),
                                     B = B_sing, x = x, y = y, family = family,
                                     control = list(trace = FALSE, tol = solnptol)),
                               TRUE)

            if( inherits(nlopt_sing, "try-error") == FALSE){
                res[["Singleton"]] <- list()
                res[["Singleton"]]$w <- nlopt_sing$pars
                res[["Singleton"]]$convergence <- nlopt_sing$convergence
                res[["Singleton"]]$elapsed <- nlopt_sing$elapsed
                coef_sing <- matrix(0, nbeta, nbeta - 1)
                for( l in 1 : (nbeta - 1) ){
                    if(family == "binomial"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = binomial())
                    } else if(family == "poisson"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = poisson())
                    }
                    coef_sing[c(1, l + 1), l] <-  coef(sing_fit)
                }
                res[["Singleton"]]$avecoef <- coef_sing %*% matrix(nlopt_sing$pars, ncol = 1)
            }else{
                res[["Singleton"]] <- nlopt_sing
            }
        }
        if("Hybrid" %in% rule){

            init_w_hybrid = matrix( 1 / ncol(hybrid_comb), nrow = ncol(hybrid_comb) )
            nlopt_hybrid <- try( solnp(pars = init_w_hybrid, fun = fun_prob_hybrid,
                                       eqfun = fun_eq_hybrid, eqB = c(1),
                                       LB = rep(0, ncol(hybrid_comb)),
                                       UB = rep(1, ncol(hybrid_comb)),
                                       B = B_hybrid, x = x, y = y,
                                       family = family, comb = hybrid_comb,
                                       control = list(trace = FALSE, tol = solnptol)),
                                 TRUE)

            if( inherits(nlopt_hybrid, "try-error") == FALSE){
                res[["Hybrid"]] <- list()
                res[["Hybrid"]]$w <- nlopt_hybrid$pars
                res[["Hybrid"]]$convergence <- nlopt_hybrid$convergence
                res[["Hybrid"]]$elapsed <- nlopt_hybrid$elapsed
                res[["Hybrid"]]$models <- hybrid_comb
                coef_hybrid <- matrix(0, nbeta, ncol(hybrid_comb))
                for(j in 1 : ncol(hybrid_comb)){

                    sel <- which(hybrid_comb[, j] %in% 1)
                    if(family == "binomial"){
                        sel_fit <- glm.fit(cbind(1, x)[, sel], y, family = binomial())
                    } else if(family == "poisson"){
                        sel_fit <- glm.fit(cbind(1, x)[, sel], y, family = poisson())
                    }
                    coef_hybrid[sel, j] <- coef(sel_fit)

                }
            }else{
                res[["Hybrid"]] <- nlopt_hybrid
            }
        }

    } else if( focus == "mu" ){

        # Fitting the LASSO model
        cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                              grouped = grouped, keep = TRUE, lambda = lambda_long,
                              type.measure = type.measure,
                              penalty.factor = penalty.factor), TRUE)

        #### VZ: Finding unique model forms
        betas_lasso = as.matrix(coef(cvglm, s = lambda_raw))  # Coefficients for hybrid lambda
        u = betas_lasso != 0                    # Logical, TRUE if coefficient is not 0
        hybrid_comb = unique(u, MARGIN = 2)     # Filtering unique model forms
        #### SJ: Check the full model. If it is already added, we do nothing.
        ####     If not added, we add the full model.
        if(all(colSums(hybrid_comb) != nbeta)){
            #### SJ: Full model is absent, we add the full model
            hybrid_comb <- cbind(hybrid_comb, 1)
        } else {
            #### SJ: If the full model is present, we make sure it is not added more than once
            full.index <- which(colSums(hybrid_comb) == nbeta)
            if(length(full.index) >= 2){
                hybrid_comb <- hybrid_comb[, -full.index[-1]]
            }
        }

        ####  SJ: fit full model (returning eta)
        ####  VZ: Fit hybrid and Singletons
        fulmod_eta <- rep( NA, nrow(x) )
        hybrid_eta <- matrix(NA, nrow = nrow(x), ncol = ncol(hybrid_comb) )
        singleton_eta <- matrix(0, nrow = nrow(x), ncol = (nbeta - 1) )
        for( i in 1 : nrow(x) ){

            x_i <- x[-i,, drop = FALSE]
            y_i <- y[-i]
            if( family == "binomial" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )

                #### VZ: Hybrid CV
                if( "Hybrid" %in% rule ){
                    for( j in 1:ncol(hybrid_comb) ){                         # Looping through model forms
                        sel <-  which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <-  glm.fit(cbind(1, x_i)[, sel], y_i,
                                            family = binomial())               # Fitting model
                        b_sel <-  matrix(coef(sel_fit))                      # Extracting coefficients
                        hybrid_eta[i, j] <- cbind(1, x)[i, sel] %*% b_sel    # Predicting# Computing eta
                    }
                }

                #### VZ: Singleton CV
                if( "Singleton" %in% rule ){
                    for( l in 1:(nbeta - 1) ){                              # Looping throug Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,          # Fitting Singletons
                                            family = binomial())
                        b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
                        singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting
                    }
                }

            } else if( family == "poisson" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )

                #### VZ: Hybrid CV
                if( "Hybrid" %in% rule ){
                    for( j in 1:ncol(hybrid_comb) ){                       # Looping through model forms
                        sel <- which(hybrid_comb[, j] %in% 1)                # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel], y_i,
                                           family = poisson())               # Fitting model
                        b_sel <-  matrix(coef(sel_fit))                      # Extracting coefficients
                        hybrid_eta[i, j] <- cbind(1, x)[i, sel] %*% b_sel    # Predicting# Computing eta

                    }
                }

                #### VZ: Singleton CV
                if( "Singleton" %in% rule ){
                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
                        singleton_eta[i, l] <- cbind(1, x[i, l]) %*% b_sing   # Predicting

                    }

                }
            }
            fulmod_eta[i] <- sum(cvfit$coefficients * c(1.0, x[i,]) )

        }

        #### SJ: need to make sure fit.preval is the linear predictor scale
        ####     which is true in version 4.0-2
        if("raw" %in% rule){

            lambda_index <- match(lambda[["raw"]], lambda_long)
            lambda[["raw"]] <- c(lambda[["raw"]], 0)
            nlambda_raw_0 <- length(lambda[["raw"]])
            cv_eta_raw <- cbind(cvglm$fit.preval[, lambda_index],
                                fulmod_eta)
            if( family == "binomial" ){
                mu_cv_raw = plogis( cv_eta_raw )
            } else if( family == "poisson" ){
                mu_cv_raw = exp( cv_eta_raw )
            }

            #### The symmetric matrix in the quadratic programming
            qmat_raw <- crossprod( matrix(y, nrow = n, ncol = nlambda_raw_0) - mu_cv_raw )
            if(qp.scale == "max"){
                qmat_raw <- qmat_raw / max(qmat_raw)
            } else if(qp.scale == "1/n"){
                qmat_raw <- qmat_raw / n
            }

            ####  Quadratic programming
            qrsolve_raw <- try( kernlab::ipop( c = rep(0.0, nlambda_raw_0),
                                               H = qmat_raw,
                                               A = matrix(1.0, nrow = 1, ncol = nlambda_raw_0),
                                               b = 1.0,
                                               l = rep(0.0, nlambda_raw_0),
                                               u = rep(1.0, nlambda_raw_0),
                                               r = 0.0),
                                TRUE)

            if(inherits(qrsolve_raw, "try-error") == FALSE){
                res[["raw"]] <- list()
                res[["raw"]]$w <- qrsolve_raw@primal
                res[["raw"]]$convergence <- qrsolve_raw@how
                res[["raw"]]$pos_def <- all(eigen(qmat_raw)$values > 0)
            }

        }
        if("Riemann" %in% rule){

        }
        if("Simpson(1/3)" %in% rule){

            lambda_index <- match(lambda[["Simpson(1/3)"]], lambda_long)
            lambda[["Simpson(1/3)"]] <- c(lambda[["Simpson(1/3)"]], 0)
            nlambda_simp_0 <- length(lambda[["Simpson(1/3)"]])
            cv_eta_simp <- cbind(cvglm$fit.preval[, lambda_index],
                                fulmod_eta)
            if( family == "binomial" ){
                mu_cv_simp = plogis( cv_eta_simp )
            } else if( family == "poisson" ){
                mu_cv_simp = exp( cv_eta_simp )
            }

            #### The symmetric matrix in the quadratic programming
            qmat_simp <- crossprod( matrix(y, nrow = n, ncol = nlambda_simp_0) - mu_cv_simp )
            if(qp.scale == "max"){
                qmat_simp <- qmat_simp / max(qmat_simp)
            } else if(qp.scale == "1/n"){
                qmat_simp <- qmat_simp / n
            }

            ####  Quadratic programming
            qrsolve_simp <- try( kernlab::ipop( c = rep(0.0, nlambda_simp_0),
                                               H = qmat_simp,
                                               A = matrix(1.0, nrow = 1, ncol = nlambda_simp_0),
                                               b = 1.0,
                                               l = rep(0.0, nlambda_simp_0),
                                               u = rep(1.0, nlambda_simp_0),
                                               r = 0.0),
                                TRUE)

            res[["Simpson(1/3)"]] <- list()
            if(inherits(qrsolve_simp, "try-error") == FALSE){

                res[["Simpson(1/3)"]]$v <- qrsolve_simp@primal
                lambda_raw_zero <- c(lambda_raw, 0.0)
                scale <- (lambda_raw_zero[1] - lambda_raw_zero[2]) / 6
                for( i in 1 : nlambda_raw ){
                    scale <- c(scale,
                               2 * (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 3,
                               (lambda_raw_zero[i] - lambda_raw_zero[i + 2]) / 6 )
                }
                scale[nlambda] <- (lambda_raw_zero[i] - lambda_raw_zero[i + 1]) / 6
                res[["Simpson(1/3)"]]$w <- qrsolve_simp@primal / scale
                res[["Simpson(1/3)"]]$convergence <- qrsolve_simp@how
                res[["Simpson(1/3)"]]$pos_def <- all(eigen(qmat_simp)$values > 0)

            }

        }
        #### VZ: Adding QP for hybrid averaging
        if( "Hybrid" %in% rule ){

            if( family == "binomial" ){
                mu_cv_hybrid = plogis( hybrid_eta )
            } else if( family == "poisson" ){
                mu_cv_hybrid = exp( hybrid_eta )
            }

            # Symmetric matrix for QP
            qmat_hybrid <- crossprod( matrix(y, nrow = n, ncol = ncol(hybrid_eta)) - mu_cv_hybrid)
            if(qp.scale == "max"){
                qmat_hybrid <- qmat_hybrid / max(qmat_hybrid)
            } else if(qp.scale == "1/n"){
                qmat_hybrid <- qmat_hybrid / n
            }

            #Quadratic programmint
            qrsolve_hybrid <- try(kernlab::ipop( c = rep(0.0, ncol(hybrid_eta)),
                                                 H = 2.0 * qmat_hybrid,
                                                 A = matrix(1.0, ncol = ncol(hybrid_eta)),
                                                 b = 1.0,
                                                 l = rep(0.0, ncol(hybrid_eta)),
                                                 u = rep(1.0, ncol(hybrid_eta)),
                                                 r = 0.0),
                                  TRUE)

            #### VZ: Adding results to function output
            res[["Hybrid"]] <- list()
            if( inherits(qrsolve_hybrid, "try-error") == FALSE){
                res[["Hybrid"]]$w <- qrsolve_hybrid@primal
                res[["Hybrid"]]$convergence <- qrsolve_hybrid@how
                res[["Hybrid"]]$models <- hybrid_comb
                res[["Hybrid"]]$pos_def <- all(eigen(qmat_hybrid)$values > 0)
            }

        }

        #### VZ: Adding QP for Singleton averaging
        if( "Singleton" %in% rule){

            if( family == "binomial" ){
                mu_cv_singleton = plogis( singleton_eta )
            } else if( family == "poisson" ){
                mu_cv_singleton = exp( singleton_eta )
            }

            # Symmetric matrix for QP
            qmat_singleton <- crossprod( matrix(y, nrow = n, ncol = ncol(singleton_eta)) - mu_cv_singleton)
            if(qp.scale == "max"){
                qmat_singleton <- qmat_singleton / max(qmat_singleton)
            } else if(qp.scale == "1/n"){
                qmat_singleton <- qmat_singleton / n
            }

            #Quadratic programming
            qrsolve_singleton <- try(kernlab::ipop( c = rep(0.0, ncol(singleton_eta)),
                                                    H = 2.0 * qmat_singleton,
                                                    A = matrix(1.0, ncol = ncol(singleton_eta)),
                                                    b = 1.0,
                                                    l = rep(0.0, ncol(singleton_eta)),
                                                    u = rep(1.0, ncol(singleton_eta)),
                                                    r = 0.0),
                                     TRUE)

            #### VZ: Adding results to function output
            res[["Singleton"]] <- list()
            if( inherits(qrsolve_singleton, "try-error") == FALSE){
                res[["Singleton"]]$w <- qrsolve_singleton@primal
                res[["Singleton"]]$convergence <- qrsolve_singleton@how
                res[["Singleton"]]$pos_def <- all(eigen(qmat_singleton)$values > 0)
            }

        }

    }
    #### ---------------------------------------------------- ####

    res

}


