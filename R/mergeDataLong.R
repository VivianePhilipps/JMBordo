
mergeDataLong <- function(longitudinal, dataLong)
{
    nlong <- length(longitudinal)
    
    ## dataLong without missing values
    withoutNALong <- removeNA(list(longitudinal[[1]]$args$subjectform,
                                   longitudinal[[1]]$args$fixedleft,
                                   longitudinal[[1]]$args$fixedright,
                                   longitudinal[[1]]$args$random,
                                   longitudinal[[1]]$args$offsetform), data=dataLong)
    newdataLong <- withoutNALong$newdata
    
    ## subjects
    idni <- rle(newdataLong[, longitudinal[[1]]$args$subject])
    nmes <- data.frame(id = idni$values, ni = idni$lengths)

    ## outcome 
    Y <- newdataLong[, longitudinal[[1]]$Yname]

    ## covariates
    covariates <- createX0(list(longitudinal[[1]]$args$fixed,
                               longitudinal[[1]]$args$random,
                               longitudinal[[1]]$args$offsetform), data=newdataLong)
    Xlong <- covariates$X0
    xnames <- colnames(Xlong)
    idef <- matrix(covariates$idform[[1]], nrow = 1)
    idea <- matrix(covariates$idform[[2]], nrow = 1)
    ido <- matrix(covariates$idform[[3]], nrow = 1)
    ## -> a faire : controler si on a le meme nombre de covariables que dans modele univarie !**

    
    ## the same for all longitudinal models
    if(nlong > 1)
    {
        id <- newdataLong[, longitudinal[[1]]$args$subject]
        
        for(k in 2:nlong)
        {
            rm(withoutNALong, covariates, idni)
            
            withoutNALong <- removeNA(list(longitudinal[[k]]$args$subjectform,
                                           longitudinal[[k]]$args$fixedleft,
                                           longitudinal[[k]]$args$fixedright,
                                           longitudinal[[k]]$args$random,
                                           longitudinal[[k]]$args$offsetform), data=dataLong)
            newdataLong <- withoutNALong$newdata
            
            ## subjects
            idni <- rle(newdataLong[, longitudinal[[k]]$args$subject])
            nmesk <- data.frame(id = idni$values, ni = idni$lengths)
            nmes <- merge(nmes, nmesk, by = c("id"), suffixes = c("", k))
            nmes[which(is.na(nmes[, paste("ni", k, sep = "")])), paste("ni", k, sep = "")] <- 0
            id <- c(id, newdataLong[, longitudinal[[k]]$args$subject])
            
            ## outcome 
            Y <- c(Y, newdataLong[, longitudinal[[k]]$Yname])
            
            ## covariates
            covariates <- createX0(list(longitudinal[[k]]$args$fixed,
                                        longitudinal[[k]]$args$random,
                                        longitudinal[[k]]$args$offsetform), data=newdataLong)
            Xtmp <- covariates$X0

            ## columns of Xtmp in the same order as Xlong
            plusnames <- setdiff(colnames(Xtmp), xnames)
            Xlongk <- matrix(0, nrow(Xtmp), length(xnames) + length(plusnames))
            colnames(Xlongk) <- c(xnames, plusnames)
            Xlongk[, sapply(colnames(Xtmp), match, colnames(Xlongk))] <- Xtmp

            if(length(plusnames)) # if Xtmp contains covariates that are not yet in Xlong
            {
                XLtmp <- cbind(Xlong, matrix(0, nrow(Xlong), length(plusnames)))
                Xlong <- rbind(XLtmp, Xlongk)
                xnames <- colnames(Xlong)
            }
            else # all covariates are already in Xlong
            {
                Xlong <- rbind(Xlong, Xlongk)
            }

            ## indicators (fixed effects, random effects, offset)
            idefk <- rep(0, length(xnames))
            idefk[sapply(colnames(Xtmp), match, colnames(Xlongk))] <- covariates$idform[[1]]
            ideak <- rep(0, length(xnames))
            ideak[sapply(colnames(Xtmp), match, colnames(Xlongk))] <- covariates$idform[[2]]
            idok <- rep(0, length(xnames))
            idok[sapply(colnames(Xtmp), match, colnames(Xlongk))] <- covariates$idform[[3]]
            
            idef <- t(apply(matrix(idef, nrow = k - 1), 1, function(x, n){c(x, rep(0, n))}, n = length(plusnames)))
            idef <- rbind(idef, idefk)
            idea <- t(apply(matrix(idea, nrow = k - 1), 1, function(x, n){c(x, rep(0, n))}, n = length(plusnames)))
            idea <- rbind(idea, ideak)
            ido <- t(apply(matrix(ido, nrow = k - 1), 1, function(x, n){c(x, rep(0, n))}, n = length(plusnames)))
            ido <- rbind(ido, idok)
        } # end k

        Y <- Y[order(id)]
        nmes <- nmes[order(nmes[,1]),]
        Xlong <- Xlong[order(id), , drop = FALSE]
    }

    return(list(Y = Y, Xlong = Xlong, idef = idef, idea = idea, ido = ido, nmes = nmes))
}
