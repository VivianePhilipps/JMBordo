JMestim <- function(longitudinal, survival, association, var.time, correlation=TRUE, dataLong, dataSurv, control=list(), ...)
{
##### Inputs #####
    
    ## check longitudinal
    if(missing(longitudinal))
    {
        nlong <- 0
        modelLongname <- NA
    }
    else
    {
        if(inherits(longitudinal, "mixedMod"))
            longitudinal <- list(longitudinal)
        else
            if(!is.list(longitudinal)) stop("long!")
        
        if(any(!sapply(longitudinal, function(x) inherits(x, "mixedMod")))) stop("longitudinal should only contain mixedMod objects")
        nlong <- length(longitudinal)

        modelLongname <- all.vars(match.call()[[2]])
        
        ## conversion to data frame
        dataLong <- as.data.frame(dataLong)
    }

    ## check survival
    if(missing(survival))
    {
        nsurv <- 0
        modelSurvname <- NA
    }
    else
    {
        if(inherits(survival, "phMod"))
            survival <- list(survival)
        else
            if(!is.list(survival)) stop("surv!")
        
        if(any(!sapply(survival, function(x) inherits(x, "phMod")))) stop("survival should only contain phMod objects")
        nsurv <- length(survival)

        modelSurvname <- all.vars(match.call()[[3]])

        ## conversion to data frame
        dataSurv <- as.data.frame(dataSurv)

        ## specify subject
        if(nlong > 0)
        {
            if(any(is.na(sapply(survival, function(x) x$args$subject)))) stop("subject is missing in survival models")
        }
    }


    ## Data for longitudinal part
    if(nlong > 0)
    {
        allDataLong <- mergeDataLong(longitudinal, dataLong)

        Y <- allDataLong$Y
        Xlong <- allDataLong$Xlong
        idef <- allDataLong$idef
        idea <- allDataLong$idea
        ido <- allDataLong$ido
        nmes <- allDataLong$nmes
        idLong <- nmes[, 1]
        ni <- as.matrix(nmes[, -1])
        nvlong <- ncol(Xlong)
        typelong <- sapply(longitudinal, function(x){x$args$typelong})
    }
    else
    {
        ## no longitudinal part
        Y <- 0
        Xlong <- 0
        idef <- matrix(0, 1, 1)
        idea <- matrix(0, 1, 1)
        ido <- 0
        ni <- 0
        nvlong <- 0
        typelong <- ""
        idLong <- NULL
    }
    
    ## number of parameters for longitudinal
    nef <- apply(idef, 1, sum)
    neatot <- sum(idea)
    nRE <- apply(idea, 1, sum)
    nvc <- sapply(nRE, function(x){x * (x + 1) / 2})
    ncovRE <- neatot * (neatot + 1) / 2 - sum(nvc)
    ncovREtot <- ncovRE
    if(correlation != TRUE) ncovRE <- 0
    nerr <- ifelse(typelong == 1, 1, 0)
    npmlong <- nef + nRE * (nRE + 1) / 2 + nerr
    npmlongtot <- npmlong
    ##!** npmlongtot et npmlong seront differents si idiag
    

    ## Data for survival part
    if(nsurv > 0)
    {
        allDataSurv <- mergeDataSurv(survival, dataSurv)

        T0 <- allDataSurv$T0
        T <- allDataSurv$T
        D <- allDataSurv$D
        rangeTsurv <- range(c(T0, T))
        Xsurv <- allDataSurv$Xsurv
        idsurv <- allDataSurv$idsurv
        nmes <- allDataSurv$nmes
        nTi <- as.matrix(nmes[, -1])
        idSurv <- nmes[, 1]
        nvsurv <- ncol(Xsurv)
        hazard <- sapply(survival, function(x){x$args$hazard})
        nhazardnodes <- sapply(survival, function(x){x$args$nhazardnodes})
        hazardnodes <- na.omit(unlist(sapply(survival, function(x){x$args$hazardnodes})))
        lefttrunc <- 0
        if(any(T0 > 0)) lefttrunc <- 1

        if(nlong > 0)
        {
            if(length(setdiff(idLong, idSurv)) | length(setdiff(idSurv, idLong))) stop("Not the same subjects in Long and Surv data")
        }
    }
    else
    {
        ## no survival part
        T0 <- 0
        T <- 0
        D <- 0
        rangeTsurv <- NA
        Xsurv <- 0
        nTi <- matrix(0, 1, 1)
        idsurv <- matrix(0, 1, 1)
        nvsurv <- 0
        hazard <- 0
        nhazardnodes <- 0
        hazardnodes <- 0
        lefttrunc <- 0
        idSurv <- NULL
    }
    ## check if all subjects have Y and T measures
    ## -> ?
    ns <- max(length(idLong), length(idSurv))

    ## check association
    if(missing(association))
    {
        nasso <- rep(0, nsurv)
        nassotot <- rep(0, nsurv)
        asso <- matrix(0, 1, 1)
    }
    else
    {
        if((nlong == 0) | (nsurv == 0)) stop("association is only needed if longitudinal and survival are specified")
        
        if(inherits(association, "formula"))
            association <- list(association)
        else
            if(!is.list(association)) stop("asso!")
        
        if(any(!sapply(association, function(x) inherits(x, "formula")))) stop("association should contain two-sided formulas")
        if(length(association) != nsurv) stop("association should contain one formula per survival model")

        npmbytype <- c(1, 3, 3) # 1 prm for linear asso, 3 for poly, 3 for logistic

        asso <- NULL # one line per association term
        jasso <- 0
        nasso <- rep(0, nsurv)
        for(ke in 1:nsurv)
        {
            aform <- association[[ke]]
            kmodelSurv <- which(modelSurvname == as.character(aform[[2]])) # which event
            if(!length(kmodelSurv)) stop(paste("unknown model", as.character(aform[[2]]), "in the association"))

            termsaform <- terms(aform)                            
            assoke <- labels(termsaform)
            orderAsso <- attr(termsaform, "order")
            
            for(k in 1:length(assoke))
            {
                if(orderAsso[k] == 1) #not an interaction term
                {
                    isAsso <- try(inherits(eval(str2lang(assoke[k])), "association"), silent = TRUE)
                    if(isAsso == TRUE) # it is an association (RE, CL, etc)
                    {
                        elementAsso <- eval(str2lang(assoke[k]))
                        jasso <- jasso + 1
                        kmodelLong <- which(modelLongname == elementAsso$name) # which longitudinal model is involved

                        ## fill the 'asso' table
                        asso <- rbind(asso, c(match(elementAsso$asso, c("RE", "CL", "CS")),
                                              kmodelLong,
                                              kmodelSurv,
                                              elementAsso$type,
                                              0))
                        ## number of parameters for the association term
                        if(asso[jasso, 1] == 1)
                            nasso[kmodelSurv] <- nasso[kmodelSurv] + nRE[kmodelLong] * npmbytype[elementAsso$type]
                        else
                            nasso[kmodelSurv] <- nasso[kmodelSurv] + elementAsso$npm
                    }
                    else # not an association
                    {
                        ## check that it is a covariate used in the survival model
                        if(!(assoke[k] %in% colnames(dataSurv)))
                            stop(paste("Unknown term in association:", assoke[k]))
                        else
                        {
                            elementmodelmatrix <- colnames(model.matrix(paste("~-1+", assoke[k]), data = dataSurv))
                            if(!all(assoke[k] %in% colnames(Xsurv)))
                                stop(paste("Variable", assoke[k], "should also appear in the survival model"))
                        }
                    }
                }
                else if(orderAsso[k] == 2) # interaction
                {
                    ## split the two elements involved in the interaction term
                    splitted <- strsplit(assoke[k], split = ":", fixed = TRUE)

                    ## which of this two element is an association
                    isAsso1 <- try(inherits(eval(str2lang(splitted[1])), "association"), silent = TRUE)
                    isAsso2 <- try(inherits(eval(str2lang(splitted[2])), "association"), silent = TRUE)
                    kasso <- which(c(isAsso1, isAsso2) == TRUE)
                    elementAsso <- eval(str2lang(assoke[kasso]))

                    ## no interaction between two associations and no interaction betwwen two covariates
                    if(all(c(isAsso1, isAsso2) == TRUE) | !any(c(isAsso1, isAsso2) == TRUE)) stop(paste("Term", assoke[k]), "not allowed in the association")

                    ## which longitudinal model is involved
                    kmodelLong <- which(modelLongname == elementAsso$name) 

                    ## covariate involved in the interaction
                    elementmodelmatrix <- colnames(model.matrix(paste("~", splitted[which(c(isAsso1, isAsso2) == FALSE)]), data = dataSurv))[-1]

                    for(l in 1:length(elementmodelmatrix))
                    {
                        jasso <- jasso + 1
                        
                        ## fill the asso table
                        asso <- rbind(asso, c(match(elementAsso$asso, c("RE", "CL", "CS")),
                                              kmodelLong,
                                              kmodelSurv,
                                              elementAsso$type,
                                              which(colnames(Xsurv) == elementmodelmatrix[l])))
                        ## number of parameters for the association term
                        if(asso[jasso, 1] == 1)
                            nasso[kmodelSurv] <- nasso[kmodelSurv] + nRE[kmodelLong] * npmbytype[elementAsso$type]
                        else
                            nasso[kmodelSurv] <- nasso[kmodelSurv] + elementAsso$npm
                    }
                }
                else
                    stop("interactions of order > 2 are not allowed in the association")

                
            } # end k
        } # end ke

        ## order asso according to event, marker, interation
        asso <- asso[order(asso[,3], asso[,2], asso[,5]),, drop = FALSE]
    }
    nlinesasso <- nrow(asso)
       
    ## number of parameters for survival
    nrisq <- rep(2, nsurv)
    nrisq[which(hazard == 2)] <- nhazardnodes[which(hazard == 2)] + 2
    nes <- apply(idsurv, 1, sum)
    npmsurv <- npmsurvtot <- nrisq + nes + nasso
    nassotot <- nasso # a supprimer !**
    npmsurvtot <- nrisq + nes + nassotot

    ## number of parameters
    npm <- sum(npmsurv) + sum(npmlong) + ncovRE # length of binit
    npmtot <- sum(npmsurvtot) + sum(npmlongtot) + ncovREtot # length of btot

    ## control
    ctrl <- list(init = NULL, posfix = NULL, nQMC = 2000, algo = "mla", maxiter = 50, epsa = 0.0001, epsb = 0.0001, epsd = 0.0001, nproc = 1, clustertype = NULL, .packages = NULL, verbose = FALSE, hessian = TRUE)
    control <- c(control, list(...))
    ctrl[names(control)] <- control
    unkn <- setdiff(names(control), names(ctrl))
    if(length(unkn)) stop(paste("Unknown options:", paste(unkn, collapse=" ")))

    ## inital values
    binit <- rep(NA, npm)
    btot <- rep(NA, npmtot)
    fix <- rep(0, npmtot)
    indexRE <- rep(NA, neatot * (neatot + 1) / 2)
    if(nsurv > 0)
    {
        before <- 0
        beforetot <- 0
        for(ke in 1:nsurv)
        {
            ## baselinerisk and time fixed covariates
            binit[before + 1:length(survival[[ke]]$estim)] <- survival[[ke]]$estim
            btot[beforetot + 1:length(survival[[ke]]$estim)] <- binit[before + 1:length(survival[[ke]]$estim)]
            before <- before + length(survival[[ke]]$estim)
            beforetot <- beforetot + length(survival[[ke]]$estim)

            ## estimated association parameters
            if(nasso[ke] > 0)
            {
                binit[before + 1:nasso[ke]] <- 0
                before <- before + nasso[ke]
                
                btot[beforetot + 1:nasso[ke]] <- 0
                beforetot <- beforetot + nasso[ke]
            }
            
            ## ## all association parameters
            ## if(nassotot[ke] > 0)
            ## {
            ##     aform <- association[[ke]]
            ##     termsaform <- terms(aform)
            ##     assoke <- eval(attr(termsaform, "variables")[-2])

            ##     getbtot <- function(x, nRE, modelLongname)
            ##     {
            ##         res <- getElement(x, "btot")
            ##         if(is.na(res)) res <- rep(0, nRE[which(modelLongname == x$name)])
            ##         return(res)
            ##     }

            ##     btot[beforetot + 1:nassotot[ke]] <- unlist(sapply(assoke, getbtot, nRE = nRE, modelLongname = modelLongname))
            ##     fix[beforetot + 1:nassotot[ke]] <- rep(0, nassotot[ke])

            ##     #btot[beforetot + 1:nassotot[ke]] <- unlist(sapply(assoke, function(x){getElement(x, "btot")}))
            ##     beforetot <- beforetot + nassotot[ke]
            ## }
        }
    }
    
    if(nlong > 0)
    {
        before <- sum(npmsurv)
        beforetot <- sum(npmsurvtot)
        for(k in 1:nlong)
        {
            ## parameters of long model k
            binit[before + 1:length(longitudinal[[k]]$estim)] <- longitudinal[[k]]$estim
            btot[beforetot + 1:length(longitudinal[[k]]$estim)] <- longitudinal[[k]]$estim
            before <- before + length(longitudinal[[k]]$estim)
            beforetot <- beforetot + length(longitudinal[[k]]$estim)
        }

        ## covariance of RE
        if(nlong > 1)
        {
            if(correlation == TRUE)
            {
                binit[before + 1:ncovRE] <- 0
            }
            else
            {
                fix[beforetot + 1:ncovREtot] <- 2
            }
            btot[beforetot + 1:ncovREtot] <- 0

            before <- before + ncovRE
            beforetot <- beforetot + ncovRE
        }
    }

    if(any(is.na(binit))) {print(binit); stop("binit!")} # a enlever !**

    ## user specified initial values
    if(!is.null(ctrl$init))
    {
        if(length(ctrl$init) != npm) stop(paste("init should be of length", npm, "instead of", length(ctrl$init)))

        binit <- ctrl$init
        btot[which(fix == 0)] <- ctrl$init
    }

    ## user specified posfix
    if(!is.null(ctrl$posfix))
    {
        if(any(!(ctrl$posfix %in% c(1:npm)))) stop(paste("posfix should contain integers from 1 to",npm))

        ## NB : posfix are indices in binit
        ## transform index in binit to index in btot
        jest <- which(fix == 0)
        jfix <- jest[ctrl$posfix]       
        fix[jfix] <- 1
    }
    
    ## index of random effects and cholesky transformation
    indexRE <- matrix(0, neatot, neatot)
    if(neatot > 0)
    {
        beforek <- sum(npmsurvtot)
        beforecov <- sum(npmsurvtot) + sum(npmlongtot)
        up <- upper.tri(indexRE, diag = TRUE)

        sumnRE <- 0
        for(k in 1:nlong)
        {
            nvc <- 0
            ## RE of model k
            if(nRE[k] > 0)
            {
                nvc <- nRE[k] * (nRE[k] + 1) / 2

                posk <- matrix(FALSE, neatot, neatot)
                posk[sumnRE + 1:nRE[k], sumnRE + 1:nRE[k]] <- TRUE

                indexRE[up & posk] <- beforek + nef[k] + 1:nvc
            }
            beforek <- beforek + nef[k] + nvc + nerr[k]

            ## covariance between RE of model k and all previous models' RE (1,...,k-1)
            if((nRE[k] > 0) & (sumnRE > 0))
            {
                indexRE[1:sumnRE, sumnRE + 1:nRE[k]] <- beforecov + 1:(sumnRE * nRE[k])

                beforecov <- beforecov + sumnRE * nRE[k]
            }

            sumnRE <- sumnRE + nRE[k]
        }
        
        indexRE <- as.numeric(indexRE[up])

        ## cholesky
        varcov <- matrix(0, neatot, neatot)
        varcov[upper.tri(varcov, diag = TRUE)] <- btot[indexRE]
        ch <- chol(varcov)
        btot[indexRE] <- ch[upper.tri(ch, diag = TRUE)]
    }

    ## estimated and unestimated parameters
    b <- btot[which(fix == 0)]
    bfix <- btot[which(fix > 0)]

    ## check var.time
    #if(any(assoLevel > 0) | any(assoSlope > 0))
    if(nlong > 0)    
    {
        if(missing(var.time)) stop("var.time is missing")
        if(length(var.time) == 1) var.time <- rep(var.time, nlong)

        if(any(!(var.time %in% colnames(dataLong)))) stop("var.time should be in dataLong")
        rangeTlong <- range(dataLong[, var.time]) # attention : toutes les donnees de dataLong ne sont pas forcement utilisees (on enleve les NA) !**
    }
    else
    {
        var.time <- NA
        rangeTlong <- NA
    }


    ## XlongT (needed for any asso CL or asso CS for Poisson outcome)
    needXlongT <- any(asso[,1] > 1) #any(assoLevel > 0) | (any(assoSlope[, which(typelong == 2)] > 0))
    sumnTi <- apply(nTi, 2, sum)
    if(needXlongT)
    {
        XlongT <- createMatrixGK(longitudinal, survival, var.time, dataLong, dataSurv, sumnTi, nvlong, deriv = FALSE)
    }
    else
    {
        ##XlongT <- matrix(0, ns * nsurv * nlong * (16 + 15 * lefttrunc), nvlong)
        XlongT <- matrix(0, sum(sumnTi) * nlong * (16 + 15 * lefttrunc), nvlong)
    }
    
    ## XlongTdt (for asso CS)
    if(any(asso[,1] == 3)) #any(assoSlope > 0))
    {
        XlongTdt <- createMatrixGK(longitudinal, survival, var.time, dataLong, dataSurv, sumnTi, nvlong, deriv = TRUE)
    }
    else
    {
        ##XlongTdt <- matrix(0, ns * nsurv * nlong * (16 + 15 * lefttrunc), nvlong)
        XlongTdt <- matrix(0, sum(sumnTi) * nlong * (16 + 15 * lefttrunc), nvlong)
    }

    ## QMC
    if(neatot > 0)
    {
        nQMC <- ctrl$nQMC[1]
        seqQMC <- randtoolbox::sobol(n = nQMC, dim = neatot, normal = TRUE, scrambling = 1)
    }
    else
    {
        nQMC <- 1
        seqQMC <- 0
    }
##### inputs OK #####

##browser()
    ## loglik arguments (without data)
    argsloglik <- list(b = b, bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv, npmlong = npmlongtot, npmsurv = npmsurvtot, Y = 0, Xlong = 0, ni = 0, nvlong = nvlong, idef = idef, idea = idea, ido = ido, typelong = typelong, Tstart = 0, Tstop = 0, Event = 0, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = 0, nTi = 0, nvsurv = nvsurv, idsurv = idsurv, asso = asso, nlinesasso = nlinesasso, XlongT = 0, XlongTdt = 0, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = ctrl$nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)
    
#browser()
##### Optimization #####
    deb <- proc.time()
    
    ## first run
    if(ctrl$maxiter[1] == 0)
    {
        ## compute log-likelihood
        vrais <- loglik(b = b, bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv, npmlong = npmlongtot, npmsurv = npmsurvtot, Y = Y, Xlong = Xlong, ni = ni, nvlong = nvlong, idef = idef, idea = idea, ido = ido, typelong = typelong, Tstart = T0, Tstop = T, Event = D, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = Xsurv, nTi = nTi, nvsurv = nvsurv, idsurv = idsurv, asso = asso, nlinesasso = nlinesasso, XlongT = XlongT, XlongTdt = XlongTdt, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)

        res <- list(istop = 2, niter = 0, loglik = vrais, b = b, v = rep(NA, length(b) * (length(b) + 1) / 2), convcrit = rep(NA, 3), time = NA, nproc = 1)       
    }
    else
    {
        if(ctrl$algo == "mla")
        {
            res <- marqLevAlg::mla(b = b, m = FALSE, fn = loglik, minimize = FALSE,
                                   clustertype = ctrl$clustertype, .packages = ctrl$.packages,
                                   epsa = ctrl$epsa, epsb = ctrl$epsb, epsd = ctrl$epsd,
                                   digits = 8, print.info = ctrl$verbose, blinding = FALSE,
                                   multipleTry = 25, file = "",
                                   nproc = ctrl$nproc, maxiter = ctrl$maxiter[1],
                                   bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv,
                                   npmlong = npmlongtot, npmsurv = npmsurvtot, Y = Y, Xlong = Xlong,
                                   ni = ni, nvlong = nvlong, idef = idef, idea = idea, ido = ido,
                                   typelong = typelong, Tstart = T0, Tstop = T, Event = D,
                                   hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = Xsurv, nTi = nTi, nvsurv = nvsurv, idsurv = idsurv,
                                   asso = asso, nlinesasso = nlinesasso, XlongT = XlongT,
                                   XlongTdt = XlongTdt, lefttrunc = lefttrunc, seqQMC = seqQMC,
                                   nQMC = nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)
            
            res <- list(istop = res$istop, niter = res$ni, counts = NA,
                        loglik = res$fn.value, b = res$b, v = res$v,
                        convcrit = c(res$ca, res$cb, res$rdm), time = res$time,
                        nproc = ctrl$nproc)
        }
        else
        {
            starts <- proc.time()
            res <- stats::optim(par =  b, fn = loglik, method = ctrl$algo, control = list(fnscale = -1), hessian = ctrl$hessian,
                                bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv,
                                npmlong = npmlongtot, npmsurv = npmsurvtot, Y = Y, Xlong = Xlong,
                                ni = ni, nvlong = nvlong, idef = idef, idea = idea, ido = ido,
                                typelong = typelong, Tstart = T0, Tstop = T, Event = D,
                                hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = Xsurv, nTi = nTi, nvsurv = nvsurv, idsurv = idsurv,
                                asso = asso, nlinesasso = nlinesasso, XlongT = XlongT,
                                XlongTdt = XlongTdt, lefttrunc = lefttrunc, seqQMC = seqQMC,
                                nQMC = nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)
            stops <- proc.time()

            if(ctrl$hessian)
                v <- solve(-res$hessian)
            else
                v <- matrix(NA, length(b), length(b))
            
            res <- list(istop = res$convergence + 1, niter = NA, counts = res$counts,
                        loglik = res$value, b = res$par, v = v[upper.tri(v, diag = TRUE)],
                        convcrit = rep(NA, 3), time = (stops - starts)[3], nproc = 1)
        }
        
    }

    ## second run
    if(length(ctrl$nQMC) == 2)
    {
        niter1 <- res$niter
        
        if(ctrl$nproc > 1)
        {
            clustpar <- parallel::makeCluster(ctrl$nproc)
            doParallel::registerDoParallel(clustpar)
        }

        ## compute Hessienne
        nQMC <- ctrl$nQMC[2]
        seqQMC <- randtoolbox::sobol(n = ctrl$nQMC[2], dim = neatot, normal = TRUE, scrambling = 1)

        derivees <- marqLevAlg::deriva(funcpa=loglik, nproc=ctrl$nproc, .packages = ctrl$.packages,  b = res$b, bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv, npmlong = npmlongtot, npmsurv = npmsurvtot, Y = Y, Xlong = Xlong, ni = ni, nvlong = nvlong, idef = idef, idea = idea, ido = ido, typelong = typelong, Tstart = T0, Tstop = T, Event = D, hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = Xsurv, nTi = nTi, nvsurv = nvsurv, idsurv = idsurv, asso = asso, nlinesasso = nlinesasso, XlongT = XlongT, XlongTdt = XlongTdt, lefttrunc = lefttrunc, seqQMC = seqQMC, nQMC = nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)
            
        if(ctrl$nproc > 1)  parallel::stopCluster(clustpar)
            
        H <- matrix(0, npm, npm)
        H[upper.tri(H, diag = TRUE)] <- derivees$v[1:(npm * (npm + 1) / 2)]
        H <- t(H)
        H[upper.tri(H, diag = TRUE)] <- derivees$v[1:(npm * (npm + 1) / 2)]
        
        Vprm <- try(solve(H), silent = TRUE)
        SEprm <- try(sqrt(diag(Vprm)), silent = TRUE)
        ##tryCatch(sqrt(diag(Vprm)), error = function(e) e, warning = function(w) w) # -> de class simpleError ou simpleWarning
        ## if it fails, optimize again
        if(!inherits(Vprm, "try-error") & all(is.finite(SEprm)))
        {
            res$v <- Vprm[upper.tri(Vprm, diag = TRUE)]
            niter2 <- NA
            res$istop <- 5
        }
        else
        {
            if(ctrl$maxiter[2] > 0)
            {
                res <- marqLevAlg::mla(b = res$b, m = FALSE, fn = loglik, minimize = FALSE,
                                       clustertype = ctrl$clustertype, .packages = ctrl$.packages,
                                       epsa = ctrl$epsa, epsb = ctrl$epsb, epsd = ctrl$epsd,
                                       digits = 8, print.info = ctrl$verbose, blinding = FALSE,
                                       multipleTry = 25, file = "",
                                       nproc = ctrl$nproc, maxiter = ctrl$maxiter[2],
                                       bfix = bfix, fix = fix, nlong = nlong, nsurv = nsurv,
                                       npmlong = npmlongtot, npmsurv = npmsurvtot, Y = Y, Xlong = Xlong,
                                       ni = ni, nvlong = nvlong, idef = idef, idea = idea, ido = ido,
                                       typelong = typelong, Tstart = T0, Tstop = T, Event = D,
                                       hazard = hazard, hazardnodes = hazardnodes, nhazardnodes = nhazardnodes, Xsurv = Xsurv, nTi = nTi, nvsurv = nvsurv, idsurv = idsurv,
                                       asso = asso, nlinesasso = nlinesasso, XlongT = XlongT,
                                       XlongTdt = XlongTdt, lefttrunc = lefttrunc, seqQMC = seqQMC,
                                       nQMC = nQMC, indexRE = indexRE - 1, nRE = nRE, ns = ns)
                niter2 <- res$niter
            }
            else
            {
                res$v <- rep(NA, length(b) * (length(b) + 1) / 2)
                niter2 <- 0
            }
        }
        
        res$niter <- c(niter1, niter2)
    }
    
    fin <- proc.time()
    res$time <- (fin - deb)[3]
##### optimization OK #####

    
##### Outputs #####
    
    ## btot
    btot <- rep(NA, length(fix))
    btot[which(fix == 0)] <- res$b
    btot[which(fix > 0)] <- bfix

    ## estim (btot with varcov instead of cholesky)
    estim <- btot
    if(neatot > 0)
    {
        ch <- matrix(0, neatot, neatot)
        ch[upper.tri(ch, diag = TRUE)] <- btot[indexRE]
 
        varcov <- t(ch) %*% ch

        estim[indexRE] <- varcov[upper.tri(varcov, diag = TRUE)]
    }
    ## a faire : renvoyer xlevels et formula/predvars
    ## renvoyer les models long et surv??
    

    res <- c(res, list(estim = estim, btot = btot, fix = fix, indexRE = indexRE, posfix = ctrl$posfix, modelLongname = modelLongname, modelSurvname = modelSurvname, asso = asso, var.time = var.time, correlation = correlation, ncovRE = ncovRE, AIC = 2 * (length(res$b) - res$loglik), BIC = length(res$b) * log(ns) - 2 * res$loglik, rangeTsurv = rangeTsurv, rangeTlong = rangeTlong), argsloglik = list(argsloglik))
    
##### outputs OK #####

    class(res) <- c("JMestim", "JMBordo")
    
    return(res)
}
