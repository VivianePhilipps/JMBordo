summary.JMBordo <- function(object, print = 3, ...)
{
    ## table with coef, se and p-value
    npm <- length(object$b)
    res <- matrix(NA, npm, 3)
    
    ## all parameters
    res[, 1] <- coef(object)

    ## if converged, add SE and p-value
    if(object$istop == 1)
    {
        Vprm <- vcov(object)
        se <- sqrt(diag(Vprm))
        res[, 2] <- se

        res[, 3] <- 1 - pchisq((res[, 1] / res[, 2])^2, 1)
    }
    colnames(res) <- c("coef", "se", "p-value")
    rownames(res) <- names(object$b)

    ## print
    if(print %in% c(1, 3))
    {
        conv <- as.numeric(object$istop == 1)
        if(conv)
            cat("Model converged! \n")
        else
            cat("Be careful, model did not converge properly! \n")
        cat("Convergence criteria: parameters =", signif(object$convcrit[1], 2), "\n")
        cat("                    : likelihood =", signif(object$convcrit[2], 2), "\n")
        cat("                    : second derivatives =", signif(object$convcrit[3], 2), "\n")
        cat("Number of Marquardt-Levenberg iterations:", object$niter, "\n")
        cat("\n")

        cat("Log-likelihood:", round(object$loglik, 2), "\n")
        cat("AIC:", round(object$AIC, 2), "\n")
        cat("BIC:", round(object$BIC, 2), "\n")
        cat("\n")
    }
    if(print %in% c(2, 3))
    {
        cat("Estimated parameters: \n")
        print(res)
    }

    
    return(invisible(res))
}


vcov.JMBordo <- function(object, which = "estimated", ...)
{
    npm <- length(object$b)
    Vprm <- matrix(NA, npm, npm)
    Vprm[upper.tri(Vprm, diag = TRUE)] <- object$v
    Vprm <- t(Vprm)
    Vprm[upper.tri(Vprm, diag = TRUE)] <- object$v

    res <- NA
    if(which == "estimated")
        res <- Vprm
    else
    {
        V <- matrix(0, length(object$btot), length(object$btot))
        V[which(object$argsloglik$fix == 0), which(object$argsloglik$fix == 0)] <- Vprm

        if(which == "all")
            res <- V[which(object$argsloglik$fix < 2), which(object$argsloglik$fix < 2)]
        else if(which == "extended")
            res <- V           
    }
    
    return(res)
}


coef.JMBordo <- function(object, which = "estimated", ...)
{
    res <- NA
    if(which == "estimated")
        res <- object$btot[which(object$argsloglik$fix == 0)]
    else if(which == "all")
        res <- object$btot[which(object$argsloglik$fix < 2)]
    else if(which == "extended")
        res <- object$btot
    
    return(res)
}
