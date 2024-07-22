
print.JMestim <- function(x, ...)
{
    cat("Joint model estimated with JMBordo::JMestim \n")

    cat("\n")
    conv <- as.numeric(x$istop == 1)
    if(conv)
        cat("Model converged! \n")
    else
        cat("Be careful, model did not converge properly! \n")
    cat("\n")

    cat("Convergence code:", x$istop, "\n")
    cat("Number of iterations:", x$niter, "\n")
    cat("Log-likelihood:", x$loglik, "\n")
    cat("\n")
    
    nlong <- x$argsloglik$nlong
    nsurv <- x$argsloglik$nsurv

    sumnpm1 <- 0
    if(nsurv > 0)
    {
        for(ke in 1:nsurv)
        {
            if(nsurv == 1)
                cat(paste("Parameters of survival submodel: \n"))
            else
                cat(paste("Parameters of survival submodel ", ke, ":\n", sep = ""))
            print(x$btot[sumnpm1 + 1:x$argsloglik$npmsurv[ke]])
            cat("\n")

            sumnpm1 <- sumnpm1 + x$argsloglik$npmsurv[ke]
        }        
    }

    sumnpm2 <- sum(x$argsloglik$npmsurv)
    if(nlong > 0)
    {
        for(k in 1:nlong)
        {
            if(nlong == 1)
                cat(paste("Parameters of longitudinal submodel: \n"))
            else
                cat(paste("Parameters of longitudinal submodel ", k, ":\n", sep = ""))
            print(x$btot[sumnpm2 + 1:x$argsloglik$npmlong[k]])
            cat("\n")

            sumnpm2 <- sumnpm2 + x$argsloglik$npmlong[k]
        }
    }

    if(x$ncovRE > 0)
    {
        cat("Correlation parameters:\n")
        print(x$btot[(sum(x$argsloglik$npmsurv) + sum(x$argsloglik$npmlong) + 1):length(x$btot)])
    }

    return(invisible(x))
}




print.mixedMod <- function(x, ...)
{    
    cat("Mixed model estimated with JMBordo::mixedMod \n")

    cat("\n")
    conv <- as.numeric(x$istop == 1)
    if(conv)
        cat("Model converged! \n")
    else
        cat("Be careful, model did not converge properly! \n")
    cat("\n")
    cat("Convergence code:", x$istop, "\n")
    cat("Number of iterations:", x$niter, "\n")
    cat("Log-likelihood:", x$loglik, "\n")
    cat("\n")

    cat(paste("Parameters: \n"))
    print(x$btot)
    cat("\n")

    return(invisible(x))
}


print.phMod <- function(x, ...)
{    
    cat("Proportional hazard model estimated with JMBordo::phMod \n")

    cat("\n")
    conv <- as.numeric(x$istop == 1)
    if(conv)
        cat("Model converged! \n")
    else
        cat("Be careful, model did not converge properly! \n")
    cat("\n")
    cat("Convergence code:", x$istop, "\n")
    cat("Number of iterations:", x$niter, "\n")
    cat("Log-likelihood:", x$loglik, "\n")
    cat("\n")

    cat(paste("Parameters: \n"))
    print(x$btot)
    cat("\n")

    return(invisible(x))
}
