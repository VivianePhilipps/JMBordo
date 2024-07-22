mergeDataSurv <- function(survival, dataSurv)
{
    ## dataSurv without missing values
    withoutNASurv <- removeNA(list(survival[[1]]$args$subjectform,
                                   survival[[1]]$args$survleft,
                                   survival[[1]]$args$survright), data=dataSurv)
    newdataSurv <- withoutNASurv$newdata
    
    ## subjects
    idni <- rle(newdataSurv[, survival[[1]]$args$subject]) 
    nmes <- data.frame(id = idni$values, ni = idni$lengths)

    ## outcome
    if(survival[[1]]$T0name == "")
    {
        T0 <- rep(0, nrow(newdataSurv))
    }
    else
    {
        T0 <- newdataSurv[, survival[[1]]$T0name]
    }
    T <- newdataSurv[, survival[[1]]$Tname]
    D <- newdataSurv[, survival[[1]]$Dname]

    ## covariates
    covariates <- createX0(list(survival[[1]]$args$survform), data=newdataSurv)
    Xsurv <- covariates$X0
    xnames <- colnames(Xsurv)
    idsurv <- matrix(covariates$idform[[1]], nrow = 1)

    nsurv <- length(survival)
    ## the same for all survival models
    if(nsurv > 1)
    {
        dat <- data.frame(id = newdataSurv[, survival[[1]]$args$subject],
                          T0 = T0, T = T, D = D, cause = 1, Xsurv)
        xnames <- colnames(Xsurv)
        
        
        for(ke in 2:nsurv)
        {
            rm(withoutNASurv, newdataSurv, covariates, T0, T, D)

            withoutNASurv <- removeNA(list(survival[[ke]]$args$subjectform,
                                           survival[[ke]]$args$survleft,
                                           survival[[ke]]$args$survright), data=dataSurv)
            newdataSurv <- withoutNASurv$newdata
            
            ## subjects
            idni <- rle(newdataSurv[, survival[[1]]$args$subject]) 
            nmes <- merge(nmes, data.frame(id = idni$values, n = idni$lengths), all = TRUE)
            nmes[which(is.na(nmes[, 1 + ke])), 1 + ke] <- 0

            ## outcomes
            if(survival[[ke]]$T0name == "")
            {
                T0 <- rep(0, nrow(newdataSurv))
            }
            else
            {
                T0 <- newdataSurv[, survival[[ke]]$T0name]
            }
            T <- newdataSurv[, survival[[ke]]$Tname]
            D <- newdataSurv[, survival[[ke]]$Dname]

            ## covariates
            covariates <- createX0(list(survival[[ke]]$args$survform), data=newdataSurv)
            Xtmp <- covariates$X0

            ## columns of Xtmp that are not yet in Xsurv
            plusnames <- setdiff(colnames(Xtmp), xnames)#;browser()
###
            ## ## union of subjects and covariates
            ## #dat <- merge(dat,
            ## #             data.frame(id = newdataSurv[, survival[[ke]]$args$subject], as.data.frame(matrix(-1, nrow = nrow(dat), ncol = length(plusnames), dimnames = list(1:nrow(dat), plusnames)))),
            ## #             all = TRUE)
            ## dat <- cbind(dat,
            ##              data.frame(id = newdataSurv[, survival[[ke]]$args$subject], as.data.frame(matrix(-1, nrow = nrow(dat), ncol = length(plusnames), dimnames = list(1:nrow(dat), plusnames)))))
            ## datk <- merge(data.frame(id = newdataSurv[, survival[[ke]]$args$subject],
            ##                          T0 = T0, T = T, D = D, cause = ke, Xtmp),
            ##               dat[, "id", drop=FALSE],
            ##               all = TRUE)
            ## datk <- datk[, colnames(dat)] # same order as dat
            ## ## revoir le cas ou on a moins de X dans phMod2 que dans phMod1 !**
###
            fivenames <- c("id", "T0", "T", "D", "cause")
            datk <- matrix(0, nrow(Xtmp), 5 + length(xnames) + length(plusnames))
            colnames(datk) <- c(fivenames, xnames, plusnames)
            datk[, sapply(c(fivenames, colnames(Xtmp)), match, colnames(datk))] <- cbind(newdataSurv[, survival[[ke]]$args$subject], T0, T, D, ke, Xtmp)
            
            if(length(plusnames))
            {
                dat <- cbind(dat, matrix(0, nrow(dat), length(plusnames)))
                colnames(dat) <- colnames(datk)
                dat <- rbind(dat, datk)
                xnames <- c(xnames, plusnames)
            }
            else
            {
                dat <- rbind(dat, datk)
            }

            ## indicator
            idsurvk <- rep(0, length(xnames))
            idsurvk[sapply(colnames(Xtmp), match, xnames)] <- covariates$idform[[1]]
            idsurv <- t(apply(matrix(idsurv, nrow = ke - 1), 1, function(x, n){c(x, rep(0, n))}, n = length(plusnames)))
            idsurv <- rbind(idsurv, idsurvk)

            rm(Xtmp, datk)
        }
        
        dat <- dat[order(dat$id, dat$cause),]
        T0 <- dat$T0
        T <- dat$T
        D <- dat$D
        Xsurv <- dat[, -c(1, 2, 3, 4, 5), drop = FALSE]
        
        ## subjects
        ##idni <- rle(dat[, "id"]) # ou compter slt event pas NA??
        ##nmes <- data.frame(id = idni$values, ni = idni$lengths)
        nmes[which(is.na(nmes[, 2])), 2] <- 0
    }

    return(list(T0 = T0, T = T, D = D, Xsurv = Xsurv, idsurv = idsurv, nmes=nmes))
}
