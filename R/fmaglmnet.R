fma.glmnet <-  function( x, y, focus, family, nlambda = 100, nfolds = NULL, grouped,
                         penalty.factor = NULL, force.nlambda, singleton.intercept.only = FALSE,
                         type.measure, rule, solnptol, qp.scale, glmnetfit = NULL ){

    ####       x must not contain the intercept column.
    #### SJ: Define some constants
    n <- nrow(x) #### Sample size
    nbeta <- ncol(x) + 1 #### Number of coefficients, including intercept
    res <- list() # Initiate returned list
    if(is.null(nfolds) == TRUE){
        nfolds = n
    }
    #### Create folds for CV
    folds <- caret::createFolds(1 : n, k = nfolds)

    ########## Functions needed for weight estimations ##########
    #### SJ: Only useful if focus = "beta"
    if( focus == "eta" ){

        ####  SJ: fun_prob computes the cross-validation MSE
        fun_prob = function(w, Eta, y, family){

            ave.eta <- Eta %*% matrix(w, ncol = 1)
            if( family == "binomial" ){
                mu = plogis(ave.eta)
            } else if( family == "poisson" ){
                mu = exp(ave.eta)
            }
            return(mean((y - mu) ^ 2))

        }

        ####  Sum of the weights
        fun_eq = function(w, Eta, y, family){
            return(sum(w))
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
    if(is.null(penalty.factor) == TRUE){
        penalty.factor <- rep(1.0, ncol(x))
    } else if(penalty.factor == "adaptiveLasso"){
        penalty.factor <- 1.0 / abs(fulcoef[-1])
    }
    res$penalty.factor <- penalty.factor
    #### ---------------------------------------------------- ####

    #############################################################
    #### SJ: Run glmnet() to get the lambda sequence. (Any better way?)
    #### Note: lambda sequence is decreasing.
    if(is.null(glmnetfit) == TRUE){
        glmnetfit = glmnet(x, y, family = family, nfolds = nfolds,
                             grouped = grouped, nlambda = nlambda,
                             type.measure = type.measure,
                             penalty.factor = penalty.factor)
    }
    lambda <- glmnetfit$lambda
    if(force.nlambda == TRUE & length(lambda) != nlambda){
        #### I think glmnet use equal space on the log scale.
        glmnetfit = glmnet(x, y, family = family, nfolds = nfolds,
                           grouped = grouped,
                           lambda = exp(seq(max(log(lambda)), min(log(lambda)), length.out = nlambda)),
                           type.measure = type.measure,
                           penalty.factor = penalty.factor)
        lambda <- glmnetfit$lambda
    }  else {
        lambda <- glmnetfit$lambda
    }
    rm(nlambda)

    #### ---------------------------------------------------- ####
    #### SJ: Using the above lambda sequence, we can start the CV step.
    if( focus == "eta" ) {

        #### VZ: Adding corresponding storage for hybrid and singletons
        if("Hybrid" %in% rule){

            betas_lasso <-  as.matrix(coef(glmnetfit, s = lambda))  # Coefficients for hybrid lambda
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
            eta_hybrid = matrix(0, n, ncol(hybrid_comb))

        }
        if("Singleton" %in% rule){
            if(singleton.intercept.only == TRUE){
                eta_sing = matrix(0, n, nbeta)
            } else {
                eta_sing = matrix(0, n, nbeta - 1)
            }
        }

        ####  CV part. This part takes a long time, especially Hybrid and singleton.
        fuleta <- rep(NA, n)
        for(i in 1 : nfolds){

            x_i <- x[-folds[[i]], , drop = FALSE]
            y_i <- y[-folds[[i]]]

            ####  full model and others.
            if( family == "binomial" ){

                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    for( j in 1:ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel, drop = F], y_i,
                                           family = binomial())             # Fitting model
                        eta_hybrid[folds[[i]], j] <- cbind(1.0, x[folds[[i]], , drop = F])[, sel, drop = F] %*% matrix(coef(sel_fit), ncol = 1)
                    }

                }

                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = binomial())
                        eta_sing[folds[[i]], l] <- cbind(1.0, x[folds[[i]], l, drop = F]) %*% matrix(coef(sing_fit), ncol = 1)

                    }
                    if(singleton.intercept.only == TRUE){
                        eta_sing[folds[[i]], nbeta] = qlogis(mean(y_i))
                    }

                }

            } else if( family == "poisson" ){

                ####  Full model.
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = poisson() )

                #### VZ: Adding hybrid and singleton CV
                if( "Hybrid" %in% rule ){

                    for( j in 1 : ncol(hybrid_comb) ){                      # Looping through model forms

                        sel <- which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <- glm.fit(cbind(1, x_i)[, sel, drop = F], y_i,
                                           family = poisson())             # Fitting model
                        eta_hybrid[folds[[i]], j] <- cbind(1.0, x[folds[[i]], , drop = F])[, sel, drop = F] %*% matrix(coef(sel_fit), ncol = 1)

                    }

                }

                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons

                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        eta_sing[folds[[i]], l] <- cbind(1.0, x[folds[[i]], l, drop = F]) %*% matrix(coef(sing_fit), ncol = 1)

                    }
                    if(singleton.intercept.only == TRUE){
                        eta_sing[folds[[i]], nbeta] = log(mean(y_i))
                    }

                }

            }

            fuleta[folds[[i]]] <- cbind(1.0, x[folds[[i]], , drop = F]) %*% matrix(cvfit$coefficients, ncol = 1)

        }
        ####  SJ: The essence of averaging beta is just to average the linear predictor.
        ####      Hence, we can use cv.glmnet to help us to get the CV linear predictor.
        cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                              grouped = grouped, keep = TRUE, lambda = lambda,
                              type.measure = type.measure,
                              penalty.factor = penalty.factor), TRUE)
        cv.eta <- cbind(cvglm$fit.preval, fuleta)

        #### VZ: Adds cvglmnet object to output
        res[["cvglmnet"]] <- cvglm

        #### SJ: In order to estimate the weights, we add 0 to the lambda sequence here.
        ####     Need to make sure that lambda remains a decreasing sequence!
        if("IMSE" %in% rule){

            res[["IMSE"]] <- NULL
            coef_raw = cbind(as.matrix(coef(cvglm, s = lambda)),
                             fulcoef) # We get out the coefficients first.
            lambda <- c(lambda, 0)
            init_w = rep( 1 / length(lambda), length(lambda) )
            nlopt_raw <- try( Rsolnp::solnp(pars = init_w, fun = fun_prob, eqfun = fun_eq, eqB = c(1),
                                            LB = rep(0, length(lambda)),
                                            UB = rep(1, length(lambda)),
                                            Eta = cv.eta, y = y, family = family,
                                            control = list(trace = FALSE, tol = solnptol)),
                              TRUE )
            ####  If nonlinear programming works, we can get the v function according to the approximation method.
            if( inherits(nlopt_raw, "try-error") == FALSE ){

                res[["IMSE"]] <- list()
                res[["IMSE"]]$w <- nlopt_raw$pars
                res[["IMSE"]]$convergence <- nlopt_raw$convergence #### 0 = converged
                res[["IMSE"]]$avecoef <- coef_raw %*% matrix(nlopt_raw$pars, ncol = 1)
                res[["IMSE"]]$elapsed <- nlopt_raw$elapsed
                res[["IMSE"]]$lambda <- lambda

            }

        }


        if("Singleton" %in% rule){

            if(singleton.intercept.only == TRUE){
                init_w_sing = matrix( 1 / nbeta, nrow = nbeta )
                nlopt_sing <- try( Rsolnp::solnp(pars = init_w_sing, fun = fun_prob, eqfun = fun_eq,
                                                 eqB = c(1),
                                                 LB = rep(0, nbeta), UB = rep(1, nbeta),
                                                 Eta = eta_sing, y = y, family = family,
                                                 control = list(trace = FALSE, tol = solnptol)),
                                   TRUE)
            } else {
                init_w_sing = matrix( 1 / (nbeta - 1), nrow = (nbeta - 1) )
                nlopt_sing <- try( Rsolnp::solnp(pars = init_w_sing, fun = fun_prob, eqfun = fun_eq,
                                                 eqB = c(1),
                                                 LB = rep(0, (nbeta - 1)), UB = rep(1, (nbeta - 1)),
                                                 Eta = eta_sing, y = y, family = family,
                                                 control = list(trace = FALSE, tol = solnptol)),
                                   TRUE)
            }

            res[["Singleton"]] <- list()
            if( inherits(nlopt_sing, "try-error") == FALSE){

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
                if(singleton.intercept.only == TRUE){
                    coef_sing = cbind(coef_sing, rep(0, nbeta))
                    coef_sing[1, nbeta] = log(mean(y))
                }
                res[["Singleton"]]$avecoef <- coef_sing %*% matrix(nlopt_sing$pars, ncol = 1)
            }

        }
        if("Hybrid" %in% rule){

            init_w_hybrid = matrix( 1 / ncol(hybrid_comb), nrow = ncol(hybrid_comb) )
            nlopt_hybrid <- try( Rsolnp::solnp(pars = init_w_hybrid, fun = fun_prob,
                                               eqfun = fun_eq, eqB = c(1),
                                               LB = rep(0, ncol(hybrid_comb)),
                                               UB = rep(1, ncol(hybrid_comb)),
                                               Eta = eta_hybrid, y = y, family = family,
                                               control = list(trace = FALSE, tol = solnptol)),
                                 TRUE)

            res[["Hybrid"]] <- list()
            if( inherits(nlopt_hybrid, "try-error") == FALSE){

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
                res[["Hybrid"]]$avecoef <- coef_hybrid %*% matrix(nlopt_hybrid$pars, ncol = 1)

            }

        }

    } else if( focus == "mu" ){

        # Fitting the LASSO model
        cvglm = try(cv.glmnet(x, y, family = family, nfolds = nfolds,
                              grouped = grouped, keep = TRUE, lambda = lambda,
                              type.measure = type.measure,
                              penalty.factor = penalty.factor), TRUE)
        #### VZ: Adds cvglmnet object to output
        res[["cvglmnet"]] <- cvglm

        #### VZ: Finding unique model forms
        betas_lasso = as.matrix(coef(cvglm, s = lambda))  # Coefficients for hybrid lambda
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
        if(singleton.intercept.only == TRUE){
            singleton_eta <- matrix(0, nrow = nrow(x), ncol = nbeta )
        } else {
            singleton_eta <- matrix(0, nrow = nrow(x), ncol = (nbeta - 1) )
        }

        ####  SJ: It is possible to write the following loop in Rcpp to gain speed.
        for( i in 1 : nfolds ){

            x_i <- x[-folds[[i]], , drop = FALSE]
            y_i <- y[-folds[[i]]]

            if( family == "binomial" ){
                cvfit <- glm.fit( x = cbind(1, x_i), y = y_i, family = binomial() )

                #### VZ: Hybrid CV
                if( "Hybrid" %in% rule ){
                    for( j in 1:ncol(hybrid_comb) ){                         # Looping through model forms
                        sel <-  which(hybrid_comb[, j] %in% 1)               # Indices of regressors to use
                        sel_fit <-  glm.fit(cbind(1, x_i)[, sel, drop = F], y_i,
                                            family = binomial())               # Fitting model
                        b_sel <-  matrix(coef(sel_fit), ncol = 1)                      # Extracting coefficients
                        hybrid_eta[folds[[i]], j] <- cbind(1, x)[folds[[i]], sel] %*% b_sel    # Predicting# Computing eta
                    }
                }

                #### VZ: Singleton CV
                if( "Singleton" %in% rule ){
                    for( l in 1:(nbeta - 1) ){                              # Looping throug Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,          # Fitting Singletons
                                            family = binomial())
                        b_sing <-  matrix(coef(sing_fit))                     # Extracting coefficients
                        singleton_eta[folds[[i]], l] <- cbind(1.0, x[folds[[i]], l, drop = F]) %*% b_sing   # Predicting
                    }
                    if(singleton.intercept.only == TRUE){
                        singleton_eta[folds[[i]], nbeta] = qlogis(mean(y_i))
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
                        hybrid_eta[folds[[i]], j] <- cbind(1, x)[folds[[i]], sel] %*% b_sel    # Predicting# Computing eta

                    }
                }

                #### VZ: Singleton CV
                if( "Singleton" %in% rule ){

                    for( l in 1:(nbeta - 1) ){                              # Looping through Singletons
                        sing_fit <- glm.fit(cbind(1, x_i[ , l]), y_i,         # Fitting Singletons
                                            family = poisson())
                        b_sing <-  matrix(coef(sing_fit), ncol = 1)                     # Extracting coefficients
                        singleton_eta[folds[[i]], l] <- cbind(1, x[folds[[i]], l]) %*% b_sing   # Predicting

                    }
                    if(singleton.intercept.only == TRUE){
                        singleton_eta[folds[[i]], nbeta] = log(mean(y_i))
                    }

                }
            }

            fulmod_eta[folds[[i]]]  <- cbind(1.0, x[folds[[i]], , drop = F]) %*% matrix(cvfit$coefficients, ncol = 1)

        }

        #### SJ: need to make sure fit.preval is the linear predictor scale
        ####     which is true in version 4.0-2
        all.coef <- as.matrix(coef(cvglm, s = lambda))
        if("IMSE" %in% rule){

            lambda <- c(lambda, 0)
            cv_eta_raw <- cbind(cvglm$fit.preval,
                                fulmod_eta)
            if( family == "binomial" ){
                mu_cv_raw = plogis( cv_eta_raw )
            } else if( family == "poisson" ){
                mu_cv_raw = exp( cv_eta_raw )
            }

            #### The symmetric matrix in the quadratic programming
            qmat_raw <- crossprod( matrix(y, nrow = n, ncol = length(lambda)) - mu_cv_raw )
            if(qp.scale == "max"){
                qmat_raw <- qmat_raw / max(qmat_raw)
            } else if(qp.scale == "1/n"){
                qmat_raw <- qmat_raw / n
            }

            ####  Quadratic programming
            qrsolve_raw <- try( kernlab::ipop( c = rep(0.0, length(lambda)),
                                               H = qmat_raw,
                                               A = matrix(1.0, nrow = 1, ncol = length(lambda)),
                                               b = 1.0,
                                               l = rep(0.0, length(lambda)),
                                               u = rep(1.0, length(lambda)),
                                               r = 0.0),
                                TRUE)

            res[["IMSE"]] <- NULL
            if(inherits(qrsolve_raw, "try-error") == FALSE){
                res[["IMSE"]] <- list()
                res[["IMSE"]]$w <- qrsolve_raw@primal
                res[["IMSE"]]$convergence <- qrsolve_raw@how
                res[["IMSE"]]$pos_def <- all(eigen(qmat_raw)$values > 0)
                res[["IMSE"]]$lambda <- lambda
                res[["IMSE"]]$coef <- all.coef # excluding the full model
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
            res[["Hybrid"]] <- NULL
            if( inherits(qrsolve_hybrid, "try-error") == FALSE){

                res[["Hybrid"]] <- list()
                res[["Hybrid"]]$w <- qrsolve_hybrid@primal
                res[["Hybrid"]]$convergence <- qrsolve_hybrid@how
                res[["Hybrid"]]$models <- hybrid_comb
                res[["Hybrid"]]$pos_def <- all(eigen(qmat_hybrid)$values > 0)
                res[["Hybrid"]]$coef <- matrix(0, nbeta, ncol(hybrid_comb))
                for( j in 1 : ncol(hybrid_comb) ){
                    sel <-  which(hybrid_comb[, j] %in% 1)
                    if(family == "binomial"){
                        sel_fit <-  glm.fit(cbind(1, x)[, sel], y, family = binomial())
                    } else if(family == "poisson"){
                        sel_fit <-  glm.fit(cbind(1, x)[, sel], y, family = poisson())
                    }
                    res[["Hybrid"]]$coef[sel, j] <- coef(sel_fit)
                }

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
                res[["Singleton"]]$coef <- matrix(0, nbeta, nbeta - 1)
                for( l in 1 : (nbeta - 1) ){
                    if(family == "binomial"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = binomial())
                    } else if(family == "poisson"){
                        sing_fit <- glm.fit(cbind(1, x[ , l]), y, family = poisson())
                    }
                    res[["Singleton"]]$coef[c(1, l + 1), l] <- coef(sing_fit) ##  add +1 for the intercept.
                }
                if(singleton.intercept.only == TRUE){

                    res[["Singleton"]]$coef <- cbind(res[["Singleton"]]$coef,
                                                     matrix(0, nbeta, 1))
                    if(family == "binomial"){
                        res[["Singleton"]]$coef[1, nbeta] = qlogis(mean(y))
                    } else if(family == "poisson"){
                        res[["Singleton"]]$coef[1, nbeta] = log(mean(y))
                    }

                }

            }

        }

    }
    #### ---------------------------------------------------- ####

    res

}


