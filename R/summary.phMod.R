summary.phMod <- function(object, print = TRUE, ...)
{
    if(print) cat("Proportional hazard model estimated with JMBordo::phMod \n \n")
    
    allprm <- summary.JMBordo(object, print = as.numeric(print == TRUE))

    ## baseline risk
    nrisq <- 0
    if(object$argsloglik$hazard == 1) nrisq <- 2
    if(object$argsloglik$hazard == 2) nrisq <- object$argsloglik$nhazardnodes + 2
    baz <- NA
    if(nrisq > 0) baz <- allprm[1:nrisq, ]

    ## covariates
    nes <- sum(object$argsloglik$idsurv)
    covariates <- NA
    if(nes > 0) covariates <- allprm[nrisq + 1:nes, ]

    if(print)
    {
        cat("Baseline hazard: \n")
        print(baz)

        if(nes > 0)
        {
            cat("\nCovariates: \n")
            print(covariates)
        }
    }

    return(invisible(list(baseline = baz, covariates = covariates)))
}
