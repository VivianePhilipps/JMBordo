## fixed  : two-sided formula
## random : right-sided formula
## subject : character
## family : "Gaussian" or "Poisson"
## data : data.frame
## offset : character
## control : list

mixedMod <- function(fixed, random, subject, family = "Gaussian", data, offset = NULL, control = list(), ...)
{
    ## convert to data.frame
    data <- as.data.frame(data)

    ## checks
    family <- tolower(family)
    if(!(family %in% c("gaussian", "poisson"))) stop("family!")

    ## remove NA
    subjectform <- as.formula(paste("~", subject))
    offsetform <- ~1
    if(!is.null(offset)) offsetform <- as.formula(paste("~", offset))
    fixedleft <- fixed[-3]
    fixedright <- fixed[-2]
    withoutNA <- removeNA(list(subjectform, fixedleft, fixedright, random, offsetform), data = data)
    newdata <- withoutNA$newdata

    ## sort by subject
    initialorder <- order(newdata[, subject])
    newdata <- newdata[order(newdata[, subject]), ]

    ## outcome and number of repeated measures
    ni <- as.numeric(table(newdata[, subject]))
    ns <- length(ni)
    Y <- newdata[, as.character(fixed[2])]
    typelong <- ifelse(family == "gaussian", 1, 2)

    ## covariates
    covariates <- createX0(list(fixed, random, offsetform), data = newdata)
    X <- covariates$X0
    idef <- covariates$idform[[1]]
    idea <- covariates$idform[[2]]
    ido <- covariates$idform[[3]]

    ## parameters
    nef <- sum(idef)
    nea <- sum(idea)
    nvc <- nea * (nea + 1) / 2
    nerr <- ifelse(family == "gaussian", 1, 0)
    npm <- nef + nvc + nerr
    
    ## control
    ctrl <- list(init = NULL, posfix = NULL, nQMC = 2000, maxiter = 50, epsa = 0.0001, epsb = 0.0001, epsd = 0.0001, nproc = 1, clustertype = NULL, .packages = NULL, verbose = FALSE)
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
        binit <- c(rep(0, nef), diag(nea)[upper.tri(diag(nea), diag = TRUE)], rep(1, nerr))
    } # a revoir si idiag possible !**

    ## cholesky
    if(nea > 0)
    {
        varcov <- matrix(0, nea, nea)
        varcov[upper.tri(varcov, diag = TRUE)] <- binit[nef + 1:nvc]

        ch <- chol(varcov)

        binit[nef + 1:nvc] <- ch[upper.tri(ch, diag = TRUE)]
    }

    ## estimated and unestimated parameters
    fix <- rep(0, npm)
    if(!is.null(ctrl$posfix))
    {
        if(any(!(ctrl$posfix %in% 1:npm))) stop("posfix!")
        fix[ctrl$posfix] <- 1
    }
    bfix <- binit[which(fix == 1)]
    b <- binit[which(fix == 0)]

    ## Quasi Monte-Carlo sequence
    if(nea > 0)
    {
        seqQMC <- randtoolbox::sobol(n = ctrl$nQMC, dim = nea, normal = TRUE, scrambling = 1)
        nQMC <- ctrl$nQMC
    }
    else
    {
        seqQMC <- 0
        nQMC <- 1
    }
    
    if(ctrl$maxiter == 0)
    {
        deb <- proc.time()
        
        vrais <- loglik(b = b, bfix = bfix, fix = fix, nlong = 1, nsurv = 0, npmlong = npm, npmsurv = 0, Y = Y, Xlong = X, ni = ni, nvlong = ncol(X), idef = idef, idea = idea, ido = ido, typelong = typelong, Tstart = 0, Tstop = 0, Event = 0, hazard = 0, hazardnodes = 0, nhazardnodes = 0, Xsurv = 0, nTi = 0, nvsurv = 0, idsurv = 0, asso = 0, nlinesasso = 0, XlongT = 0, XlongTdt = 0, lefttrunc = 0, seqQMC = seqQMC, nQMC = nQMC, indexRE = c(nef+1:nvc)-1, nRE = nea, ns = ns)
        fin <- proc.time()
        
        res <- list(istop = 2, niter = 0, loglik = vrais, b = b, v = rep(NA, length(b) * (length(b) + 1) / 2), convcrit = rep(NA,3), time = (fin-deb)[3], nproc = 1)
        
    }
    else
    {
        res <- marqLevAlg::mla(b = b, m = FALSE, fn = loglik, minimize = FALSE,
                               clustertype = ctrl$clustertype, .packages = ctrl$.packages,
                               epsa = ctrl$epsa, epsb = ctrl$epsb, epsd = ctrl$epsd,
                               digits = 8, print.info = ctrl$verbose, blinding = FALSE,
                               multipleTry = 25, file = "",
                               nproc = ctrl$nproc, maxiter = ctrl$maxiter,
                               bfix = bfix, fix = fix, nlong = 1, nsurv = 0, npmlong = npm, npmsurv = 0, Y = Y, Xlong = X, ni = ni, nvlong = ncol(X), idef = idef, idea = idea, ido = ido, typelong = typelong, Tstart = 0, Tstop = 0, Event = 0, hazard = 0, hazardnodes = 0, nhazardnodes = 0, Xsurv = 0, nTi = 0, nvsurv = 0, idsurv = 0, asso = 0, nlinesasso = 0, XlongT = 0, XlongTdt = 0, lefttrunc = 0, seqQMC = seqQMC, nQMC = nQMC, indexRE = c(nef+1:nvc)-1, nRE = nea, ns = ns)
        
        res <- list(istop = res$istop, niter = res$ni, loglik = res$fn.value,
                    b = res$b, v = res$v, convcrit = c(res$ca, res$cb, res$rdm),
                    time = res$time, nproc = ctrl$nproc)
    }

    ## btot
    btot <- rep(NA, length(fix))
    btot[which(fix == 0)] <- res$b
    btot[which(fix == 1)] <- bfix
    
    ## estim (btot with varcov instead of cholesky)
    estim <- btot
    if(nea > 0)
    {
        ch <- matrix(0, nea, nea)
        ch[upper.tri(ch, diag = TRUE)] <- btot[nef + 1:nvc]

        varcov <- t(ch) %*% ch

        estim[nef + 1:nvc] <- varcov[upper.tri(varcov, diag = TRUE)]
    }
    
    argsloglik <- list(bfix = bfix, fix = fix, nlong = 1, nsurv = 0, npmlong = npm, npmsurv = 0, ni = ni, nvlong = ncol(X), idef = idef, idea = idea, ido = ido, typelong = typelong, hazard = 0, hazardnodes = 0, nhazardnodes = 0, Xsurv = 0, nTi = 0, nvsurv = 0, idsurv = 0, asso = 0, nlinesasso = 0, lefttrunc = 0, seqQMC = seqQMC, nQMC = nQMC, indexRE = c(nef+1:nvc)-1, nRE = nea, ns = ns)

    ## args
    args <- list(subjectform = subjectform, fixedleft = fixedleft,
                 fixedright = withoutNA$terms[[3]], random = withoutNA$terms[[4]],
                 offsetform = offsetform, subject = subject, fixed = fixed,
                 typelong = typelong)

    ## Yname and Xnames
    Yname <- as.character(fixed[2])
    Xnames <- setdiff(withoutNA$X, c(Yname, subject))
    
    res <- c(res, args = list(args), argsloglik = list(argsloglik), list(Yname = Yname, Xnames = Xnames, btot = btot, estim = estim, AIC = 2 * (length(res$b) - res$loglik), BIC = length(res$b) * log(ns) - 2 * res$loglik))
    class(res) <- c("mixedMod", "JMBordo")
    
    return(res)
}
