## log-vraisemblance individuelle
loglik_i <- function(b, bfix, fix, nlong, nsurv, npmlong, npmsurv, Y, Xlong, ni, nvlong, idef, idea, ido, typelong, Tstart, Tstop, Event, hazard, Xsurv, nvsurv, idsurv, asso, nlinesasso, XlongT, XlongTdt, lefttrunc, seqQMC, nQMC, indexRE, nRE)
{
    .Call(C_loglik_indiv, as.double(b), as.double(bfix), as.integer(fix), as.integer(nlong), as.integer(nsurv), as.integer(npmlong), as.integer(npmsurv), as.double(Y), as.double(Xlong), as.integer(ni), as.integer(nvlong), as.integer(idef), as.integer(idea), as.integer(ido), as.integer(typelong), as.double(Tstart), as.double(Tstop), as.integer(Event), as.integer(hazard), as.double(Xsurv), as.integer(nvsurv), as.integer(idsurv), as.integer(asso), as.integer(nlinesasso), as.double(XlongT), as.double(XlongTdt), as.integer(lefttrunc), as.double(seqQMC), as.integer(nQMC), as.integer(indexRE), as.integer(nRE))
}

## log-vraisemblance totale
loglik <- function(b, bfix, fix, nlong, nsurv, npmlong, npmsurv, Y, Xlong, ni, nvlong, idef, idea, ido, typelong, Tstart, Tstop, Event, hazard, hazardnodes, nhazardnodes, Xsurv, nTi, nvsurv, idsurv, asso, nlinesasso, XlongT, XlongTdt, lefttrunc, seqQMC, nQMC, indexRE, nRE, ns)
{
    res <- 0
    np <- 16 + 15 * lefttrunc
    avtlong <- 0
    avtsurv <- 0
    avtlongT <- 0
    for(i in 1:ns)
    {
        if(nlong > 0)
        {
            if(is.matrix(ni))
                n_i <- ni[i,]
            else
                n_i <- ni[i]
            Y_i <- Y[avtlong + 1:sum(n_i)]
            Xlong_i <- Xlong[avtlong + 1:sum(n_i),]
        }
        else
        {
            n_i <- 0
            Y_i <- 0
            Xlong_i <- 0
        }

        if(nsurv > 0)
        {
            if(is.matrix(nTi))
                nT_i <- nTi[i,]
            else
                nT_i <- nTi[i]
            Tstart_i <- Tstart[avtsurv + 1:sum(nT_i)]
            Tstop_i <- Tstop[avtsurv + 1:sum(nT_i)]
            Event_i <- Event[avtsurv + 1:sum(nT_i)]
            Xsurv_i <- Xsurv[avtsurv + 1:sum(nT_i), ]
        }
        else
        {
            Tstart_i <- 0
            Tstop_i <- 0
            Event_i <- 0
            Xsurv_i <- 0
            nT_i <- 0
        }
        
        if((nsurv > 0) & (nlong > 0))
        {
            XlongT_i <- XlongT[avtlongT + 1:(np * nsurv * nlong),]
            XlongTdt_i <- XlongTdt[avtlongT + 1:(np * nsurv * nlong),]
        }
        else
        {
            XlongT_i <- 0
            XlongTdt_i <- 0
        }

        tmp <- .Call(C_loglik_indiv, as.double(b), as.double(bfix), as.integer(fix), as.integer(nlong), as.integer(nsurv), as.integer(npmlong), as.integer(npmsurv), as.double(Y_i), as.double(t(Xlong_i)), as.integer(n_i), as.integer(nvlong), as.integer(t(idef)), as.integer(t(idea)), as.integer(t(ido)), as.integer(typelong), as.double(Tstart_i), as.double(Tstop_i), as.integer(Event_i), as.integer(hazard), as.double(hazardnodes), as.integer(nhazardnodes), as.double(t(Xsurv_i)), as.integer(nT_i), as.integer(nvsurv), as.integer(t(idsurv)), as.integer(t(asso)), as.integer(nlinesasso),  as.double(t(XlongT_i)), as.double(t(XlongTdt_i)), as.integer(lefttrunc), as.double(t(seqQMC)), as.integer(nQMC), as.integer(indexRE), as.integer(nRE))
 
        
        res <- res + tmp
        avtlong <- avtlong + sum(n_i)
        avtsurv <- avtsurv + sum(nT_i)
        avtlongT <- avtlongT + np * nsurv * nlong
    }

    return(res)
}
