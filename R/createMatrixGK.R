## create matrix named XlongT (if deriv is FALSE) or XlongTdt (if deriv is TRUE) in loglik

createMatrixGK <- function(longitudinal, survival, var.time, dataLong, dataSurv, ns, nvlong, deriv = FALSE)
{
    nlong <- length(longitudinal)
    nsurv <- length(survival)
    
    lefttrunc <- 0
    if(survival[[1]]$T0name != "") lefttrunc <- 1
    np <- 16 + 15 * lefttrunc
    #res <- matrix(0, ns * nsurv * nlong * np, nvlong)
    #id <- rep(NA, ns * nsurv * nlong * np)
    res <- matrix(NA, sum(ns) * nlong * np, nvlong)
    id <- rep(NA, sum(ns) * nlong * np) # ns = nb de lignes pour chaque evt (>= nb sujets)

    avtjj <- 0
    for(ke in 1:nsurv)
    {
        ## id, T0, T, D without NA
        withoutNASurv <- removeNA(list(survival[[ke]]$args$subjectform,
                                       survival[[ke]]$args$survleft), data=dataSurv)
        newdataSurv <- withoutNASurv$newdata

        ## GK integration points on [-1, 1]
        ptsGK <- c(0,
                   0.207784955007898467600689403773245,
                   0.405845151377397166906606412076961,
                   0.586087235467691130294144838258730,
                   0.741531185599394439863864773280788,
                   0.864864423359769072789712788640926,
                   0.949107912342758524526189684047851,
                   0.991455371120812639206854697526329,
                   -0.207784955007898467600689403773245,      
                   -0.405845151377397166906606412076961,
                   -0.586087235467691130294144838258730,
                   -0.741531185599394439863864773280788,
                   -0.864864423359769072789712788640926,
                   -0.949107912342758524526189684047851,
                   -0.991455371120812639206854697526329)
        
        ## id and T0
        if(!lefttrunc)
        {
            ds0 <- matrix(0, nrow = 0, ncol = 5)
        }
        else
        {
            ds0 <- newdataSurv[, c(survival[[ke]]$args$subject, survival[[ke]]$T0name)]
            ds0 <- data.frame(id = rep(newdataSurv[, survival[[ke]]$args$subject], each = 15),
                              T0 = rep(newdataSurv[, survival[[ke]]$T0name], each = 15),
                              tgk = rep(ptsGK, ns[ke]))
            ds0$tgk <- (ds0$tgk + 1) / 2 * ds0$T0
            ds0$forsorting <- 3
            ds0$fakeid <- rep(1:ns[ke], each = 15)
        }
        colnames(ds0) <- c("id", "T", "ttttt", "forsorting", "fakeid")

        ## id and T
        ds <- newdataSurv[, c(survival[[ke]]$args$subject, survival[[ke]]$Tname)]
        ds <- data.frame(id = rep(newdataSurv[, survival[[ke]]$args$subject], each = 16),
                          T = rep(newdataSurv[, survival[[ke]]$Tname], each = 16),
                          tgk = rep(c(1, ptsGK), ns[ke]))
        ds$tgk <- (ds$tgk + 1) / 2 * ds$T
        ds$forsorting <- rep(c(1, rep(2, 15)), ns[ke])
        ds$fakeid <- rep(1:ns[ke], each = 16)
        colnames(ds) <- c("id", "T", "ttttt", "forsorting", "fakeid")

        ## data with GK integration points
        dsurv <- rbind(ds, ds0)

        ## add covariates from long models
        for(k in 1:nlong)
        {
            ## indices in res corresponding to event ke and outcome k
            ##jj <- (ke - 1) * ns * np * nlong + (k - 1) * ns * np + 1:(ns * np)
            jj <- avtjj + 1:(ns[ke] * np)

            ## id, X without NA           
            withoutNALong <- removeNA(list(longitudinal[[k]]$args$subjectform,
                                           longitudinal[[k]]$args$fixedright,
                                           longitudinal[[k]]$args$random,
                                           longitudinal[[k]]$args$offsetform), data=dataLong)
            newdataLong <- withoutNALong$newdata

            ## remove var.time
            dl <- newdataLong[, setdiff(colnames(newdataLong), var.time[k]), drop = FALSE]
            colnames(dl)[which(colnames(dl) == longitudinal[[k]]$args$subject)] <- "id"
            dl <- dl[!duplicated(dl$id), , drop = FALSE]

            ## merge X and T
            dlongT <- merge(dsurv, dl, by = "id", all.x = TRUE) # not in the right order even if sort = FALSE
            dlongT <- dlongT[order(as.numeric(dlongT$forsorting == 3), dlongT$fakeid, dlongT$forsorting),]
            id[jj] <- dlongT[, "id"]

            ## create XlongT
            fixed <- as.character(longitudinal[[k]]$args$fixedright)[2]
            fixed <- gsub(paste("\\b", var.time[k], "\\b", sep = ""), "ttttt", fixed)
            fixed <- formula(paste("~", fixed))

            random <- as.character(longitudinal[[k]]$args$random)[2]
            random <- gsub(paste("\\b", var.time[k], "\\b", sep = ""), "ttttt", random)
            random <- formula(paste("~", random))

            covariates <- createX0(list(fixed,
                                        random,
                                        longitudinal[[k]]$args$offsetform), data=dlongT)
            Xtmp <- covariates$X0

            if(deriv == FALSE)
            {
                ## columns of Xtmp in the same order as Xlong
                if(k == 1)
                {
                    res[jj, 1:ncol(Xtmp)] <- Xtmp
                    xnames <- colnames(Xtmp)
                    plusnames <- NULL
                }
                else
                {
                    plusnames <- setdiff(colnames(Xtmp), xnames)
                    Xlongk <- matrix(0, nrow(Xtmp), length(xnames) + length(plusnames))
                    colnames(Xlongk) <- c(xnames, plusnames)
                    Xlongk[, sapply(colnames(Xtmp), match, colnames(Xlongk))] <- Xtmp

                    res[jj, 1:ncol(Xlongk)] <- Xlongk
                    xnames <- colnames(Xlongk)
                }
            }

            ## derive Xlongk with respect to time "ttttt"
            if(deriv == TRUE)
            {
                h <- 0.0001 #?? !**
                h <- sapply(dlongT$ttttt, function(x){max(1E-7, 1E-4 * abs(x))})

                ## t + h
                dlongTplus <- dlongT
                dlongTplus$ttttt <- dlongTplus$ttttt + h
                covariates <- createX0(list(fixed,
                                            random,
                                            longitudinal[[k]]$args$offsetform), data=dlongTplus)
                Xtmpplus <- covariates$X0
                
                ## t - h
                dlongTmoins <- dlongT
                dlongTmoins$ttttt <- dlongTmoins$ttttt - h                
                covariates <- createX0(list(fixed,
                                            random,
                                            longitudinal[[k]]$args$offsetform), data=dlongTmoins)
                Xtmpmoins <- covariates$X0

                ## dt
                Xtmp <- (Xtmpplus - Xtmpmoins) / (2 * h)
                
                ## columns of Xtmpdt in the same order as Xlong
                if(k == 1)
                {
                    res[jj, 1:ncol(Xtmp)] <- Xtmp
                    xnames <- colnames(Xtmp)
                    plusnames <- NULL
                }
                else
                {
                    plusnames <- setdiff(colnames(Xtmp), xnames)
                    Xlongk <- matrix(0, nrow(Xtmp), length(xnames) + length(plusnames))
                    colnames(Xlongk) <- c(xnames, plusnames)
                    Xlongk[, sapply(colnames(Xtmp), match, colnames(Xlongk))] <- Xtmp
                    
                    res[jj, 1:ncol(Xlongk)] <- Xlongk
                    xnames <- colnames(Xlongk)
                }           
            }

            avtjj <- avtjj + ns[ke] * np
        } # end k
       
    } # end ke

    res <- res[order(id), , drop = FALSE]

    return(res)
}

        
