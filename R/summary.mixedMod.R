summary.mixedMod <- function(object, print = TRUE, ...)
{
    ## fixed effects
    nef <- sum(object$argsloglik$idef)
    FE <- matrix(NA, nef, 3)
    FE[, 1] <- object$btot[1:nef]

    ## random effects
    nea <- sum(object$argsloglik$idea)
    nvc <- nea * (nea + 1) / 2
    RE <- matrix(0, nea, nea)
    RE[upper.tri(RE, diag = TRUE)] <- object$btot[nef + 1:nvc]
    colnames(RE) <- names(object$btot)[nef + 1:nvc]

    ## measurement error
    nerr <- ifelse(object$argsloglik$typelong == 1, 1, 0)
    sigma <- NA
    if(nerr > 0) sigma <- object$btot[nef + nvc + 1]
    names(sigma) <- "error"

    ## if converged, add SE and p-value
    if(object$istop == 1)
    {
        Vprm <- vcov(object)
        se <- sqrt(diag(Vprm)[1:nef])
        FE[, 2] <- se

        FE[, 3] <- 1 - pchisq((FE[, 1] / FE[, 2])^2, 1)
    }
    colnames(FE) <- c("coef", "se", "p-value")
    rownames(FE) <- names(object$btot)[1:nef]

    ## print
    if(print)
    {
        cat("Mixed model estimated with JMBordo::mixedMod \n")
        cat("\n")
        
        conv <- as.numeric(object$istop == 1)
        if(conv)
            cat("Model converged! \n")
        else
            cat("Be careful, model did not converge properly! \n")
        cat("Convergence criteria: parameters=", signif(object$convcrit[1],2), "\n")
        cat("                    : likelihood=", signif(object$convcrit[2],2), "\n")
        cat("                    : second derivatives=", signif(object$convcrit[3],2), "\n")
        cat("Number of Marquardt-Levenberg iterations: ", object$niter, "\n")
        cat("\n")

        cat("Maximum log-likelihood:", round(object$loglik, 2), "\n")
        cat("\n")

        cat("Fixed effects:\n")
        print(FE)

        if(nea > 0)
        {
            cat("Matrix of variance-covariance of the random effects:\n")
            print(RE)
        }

        if(nerr > 0)
        {
            cat("Residual standard error: ", sigma, "\n")
        }
            
    }

    
    return(invisible(list(fixed = FE, random = RE, error = sigma)))
}


summary.mixedMod <- function(object, print = TRUE, ...)
{
    if(print) cat("Mixed model estimated with JMBordo::mixedMod \n \n")
    
    allprm <- summary.JMBordo(object, print = as.numeric(print == TRUE))

    ## fixed effects
    nef <- sum(object$argsloglik$idef)
    FE <- NA
    if(nef > 0) FE <- allprm[1:nef, ]

    ## random effects
    nea <- sum(object$argsloglik$idea)
    nvc <- nea * (nea + 1) / 2
    RE <- NA
    if(nea > 0)
    {
        RE <- matrix(0, nea, nea)
        RE[upper.tri(RE, diag = TRUE)] <- object$btot[nef + 1:nvc] # b ou btot???
        RE <- t(RE)
        RE[upper.tri(RE, diag = TRUE)] <- object$btot[nef + 1:nvc] # b ou btot???
        colnames(RE) <- names(object$btot)[nef + 1:nvc]
    }
    
    ## measurement error
    nerr <- ifelse(object$argsloglik$typelong == 1, 1, 0)
    sigma <- NA
    if(nerr > 0) sigma <- object$btot[nef + nvc + 1]
    names(sigma) <- "error"

    if(print)
    {
        cat("Fixed effects:\n")
        print(FE)
        
        if(nea > 0)
        {
            cat("\nMatrix of variance-covariance of the random effects:\n")
            print(RE)
        }

        if(nerr > 0)
        {
            cat("\nResidual standard error: ", sigma, "\n")
        }
    }
    
    return(invisible(list(fixed = FE, random = RE, error = sigma)))
}


