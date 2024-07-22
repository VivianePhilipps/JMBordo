
## association through random effects

RE <- function(x, nonlinear = NA)
{
    mc <- match.call()
    modelname <- as.character(mc[2])

    allNL <- c(NA, "polynomial", "logistic")
    if(!(nonlinear %in% allNL)) stop("Unknown nonlinear type")

    type <- match(tolower(nonlinear), allNL)
    
    npm <- NA
    
    btot <- NA

    res <- list(asso = "RE", name = modelname, type = type, npm = npm, btot = btot)
    class(res) <- "association"

    return(res)
}


## association through current level

CL <- function(x, nonlinear = NA)
{
    mc <- match.call()
    modelname <- as.character(mc[2])

    allNL <- c(NA, "polynomial", "logistic")
    if(!(nonlinear %in% allNL)) stop("Unknown nonlinear type")

    type <- match(tolower(nonlinear), allNL)
    
    npm <- 1 * as.numeric(type == 1) + 3 * as.numeric(type == 2) + 3 * as.numeric(type == 3)
    
    btot <- rep(0, npm)
    
    res <- list(asso = "CL", name = modelname, type = type, npm = npm, btot = btot)
    class(res) <- "association"

    return(res)
}


## association through current slope

CS <- function(x, nonlinear = NA)
{
    mc <- match.call()
    modelname <- as.character(mc[2])

    allNL <- c(NA, "polynomial", "logistic")
    if(!(nonlinear %in% allNL)) stop("Unknown nonlinear type")

    type <- match(tolower(nonlinear), allNL)
    
    npm <- 1 * as.numeric(type == 1) + 3 * as.numeric(type == 2) + 3 * as.numeric(type == 3)
    
    btot <- rep(0, npm)

    res <- list(asso = "CS", name = modelname, type = type, npm = npm, btot = btot)
    class(res) <- "association"

    return(res)
}

