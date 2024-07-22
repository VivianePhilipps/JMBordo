## survival : two-sided formula (left hand = survBordo)
## hazard : "Weibull"
## data : data.frame
## control : list

phMod <- function(survival, subject, hazard, data, control = list(), ...)
{
    ## convert to data.frame
    data <- as.data.frame(data)
    
    ## remove NA
    subjectform <- ~-1
    if(!missing(subject)) subjectform <- as.formula(paste("~", subject))
    survleft <- as.formula(paste("~", paste(survival[2], collapse = "+")))
    survright <- survival[-2]
    withoutNA <- removeNA(list(subjectform, survleft, survright), data = data)
    newdata <- withoutNA$newdata

    ## outcomes and number of subjects
    survnames <- eval(survival[[2]])
    Tstart <- 0
    lefttrunc <- 0
    if(survnames[1] != "")
    {
        Tstart <- newdata[, survnames[1]]
        lefttrunc <- 1
    }
    Tstop <- newdata[, survnames[2]]
    Event <- newdata[, survnames[3]]
    if(!missing(subject))
    {
        idni <- rle(newdata[, subject])
        nTi <- as.numeric(idni$lengths)
        ns <- length(nTi)
    }
    else
    {
        nTi <- rep(1, length(Event))
        ns <- length(Event)
    }
    
    ## covariates
    survform <- as.formula(paste("~", survival[3], "-1"))
    covariates <- createX0(list(survform), data = newdata)
    X <- covariates$X0
    idsurv <- covariates$idform[[1]]

    ## baseline hazard
    if(all(hazard == "splines")) hazard <- Msplines(df = 5) # M-splines with minimal number of nodes
    if(is.character(hazard))
    {
        hazard <- tolower(hazard)
        if(!(hazard %in% c("weibull", "splines"))) stop("hazard!")
        hazard <- switch(hazard, "weibull" = 1, "splines" = 2)

        if(hazard == 1) # Weibull
        {
            nrisq <- 2
            nhazardnodes <- 0
            hazardnodes <- 0
            binit_risq <- c(0.5, 1)
        }
    }
    else # splines(...) specified
    {
        if(any(is.na(hazard$nodes)))
        {
            tfornodes <- Tstop
            if(lefttrunc == 1) tfornodes <- c(Tstart, Tstop)
            
            hazardcall <- hazard$call
            hazardcall$x <- tfornodes
            hazard <- eval(hazardcall)
        }
        
        nhazardnodes <- hazard$nnodes
        nrisq <- hazard$npm
        hazardnodes <- hazard$nodes
        hazard <- 2
        binit_risq <- rep(0.1, nrisq)
    }
    
    ## parameters
    nes <- sum(idsurv)
    npm <- nrisq + nes

    ## control
    ctrl <- list(init = NULL, posfix = NULL, nQMC = 1, maxiter = 50, epsa = 0.0001, epsb = 0.0001, epsd = 0.0001, nproc = 1, clustertype = NULL, .packages = NULL, verbose = FALSE)
    control <- c(control, list(...))
    ctrl[names(control)] <- control

    
    ## initial values
    if(!is.null(ctrl$init))
    {
        binit <- ctrl$init
        if(length(binit) != npm) stop(paste("init should be of length", npm, "and it is of length", length(binit)))
    }
    else
    {
        binit <- c(binit_risq, rep(0, nes))
    }
    
    ## estimated and unestimated parameters
    fix <- rep(0, npm)
    if(!is.null(ctrl$posfix))
    {
        if(any(!(ctrl$posfix %in% 1:npm))) stop("posfix!")
        fix[ctrl$posfix] <- 1
    }
    bfix <- binit[which(fix  == 1)]
    b <- binit[which(fix == 0)]

    ## Quasi Monte-Carlo sequence is not needed
    seqQMC <- 0
    nQMC <- 1

    ## optmization
    if(ctrl$maxiter == 0)
    {
        deb <- proc.time()
 
        vrais <- loglik(b = b, bfix = bfix, fix = fix, nlong = 0, nsurv = 1, npmlong = 0, npmsurv = npm, Y = 0, Xlong = 0, ni = 0, nvlong = 0, idef = 0, idea = 0, ido = 0, typelong = 0, Tstart = Tstart, Tstop = Tstop, Event = Event, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = X, nTi = nTi, nvsurv = ncol(X), idsurv = idsurv, asso = 0, nlinesasso = 0, XlongT = 0, XlongTdt = 0, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = nQMC, indexRE = 0, nRE = 0, ns = ns)
        fin <- proc.time()
        
        res <- list(istop = 2, niter = 0, loglik = vrais, b = b, v = rep(NA,length(b)), convcrit = rep(NA,3), time = (fin-deb)[3], nproc = 1)
        
    }
    else
    {
        res <- marqLevAlg::mla(b = b, m = FALSE, fn = loglik, minimize = FALSE,
                               clustertype = ctrl$clustertype, .packages = ctrl$.packages,
                               epsa = ctrl$epsa, epsb = ctrl$epsb, epsd = ctrl$epsd,
                               digits = 8, print.info = ctrl$verbose, blinding = FALSE,
                               multipleTry = 25, file = "",
                               nproc = ctrl$nproc, maxiter = ctrl$maxiter,
                               bfix = bfix, fix = fix, nlong = 0, nsurv = 1, npmlong = 0, npmsurv = npm, Y = 0, Xlong = 0, ni = 0, nvlong = 0, idef = 0, idea = 0, ido = 0, typelong = 0, Tstart = Tstart, Tstop = Tstop, Event = Event, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = X, nTi = nTi, nvsurv = ncol(X), idsurv = idsurv, asso = 0, nlinesasso = 0, XlongT = 0, XlongTdt = 0, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = nQMC, indexRE = 0, nRE = 0, ns = ns)
        
        res <- list(istop = res$istop, niter = res$ni, loglik = res$fn.value,
                    b = res$b, v = res$v, convcrit = c(res$ca, res$cb, res$rdm),
                    time = res$time, nproc = ctrl$nproc)
    }
    

    ## btot
    btot <- rep(NA, length(fix))
    btot[which(fix == 0)] <- res$b
    btot[which(fix == 1)] <- bfix

    argsloglik <- list(bfix = bfix, fix = fix, nlong = 0, nsurv = 1, npmlong = 0, npmsurv = npm, ni = 0, nvlong = 0, idef = 0, idea = 0, ido = 0, typelong = 0, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = X, nTi = nTi, nvsurv = ncol(X), idsurv = idsurv, asso = 0, nlinesasso = 0, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = nQMC, indexRE = 0, nRE = 0, ns = ns)

    ## args
    if(missing(subject)) subject <- NA
    if(nhazardnodes == 0) hazardnodes <- NA
    args <- list(subjectform = subjectform, survleft = survleft,
                 survright = withoutNA$terms[[3]], 
                 subject = subject, survform = survform,
                 hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes)

    ## Tname and Xnames
    T0name <- survnames[1]
    Tname <- survnames[2]
    Dname <- survnames[3]
    Xnames <- setdiff(withoutNA$X, survnames)
    
    res <- c(res, args = list(args), argsloglik = list(argsloglik), list(T0name = T0name, Tname = Tname, Dname = Dname, Xnames = Xnames, btot = btot, estim = btot, AIC = 2 * (length(res$b) - res$loglik), BIC = length(res$b) * log(ns) - 2 * res$loglik, lefttrunc = lefttrunc))
    class(res) <- c("phMod", "JMBordo")
    
    return(res)
}


survBordo <- function(time, time2, event)
{
    mc <- match.call()

    if(length(mc) == 4)
    {
        T0 <- as.character(match.call()[2])
        T <- as.character(match.call()[3])
        D <- as.character(match.call()[4])
    }
    else
    {
        if(length(mc) == 3)
        {
            T0 <- ""
            T <- as.character(match.call()[2])
            D <- as.character(match.call()[3])
        }
        else
        {
            stop("Specify at least time and event")
        }
    }
    
    return(c(T0, T, D))
}
