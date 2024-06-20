#' A testing method for inference of microbial networks
#' 
#' This function generates a p-value and q-value for each (linear, nonlinear, and omnibus) test of a pair of taxa.
#' 
#' @param otu.tab the \code{n.obs} by \code{n.otu} matrix of read counts. 
#' @param clustered.data a logical variable indicating whether the samples are clustered. The default is FALSE. 
#' @param CTRL the phylogeneic tree. The default is NULL.
#' @param fdr.nominal the nominal FDR level. The default is 0.1. 
#' @param n.perm.max the maximum number of permutation replicates. The default is NULL, in which case a maximum of \code{n.otu} * \code{n.rej.stop} * (1/\code{fdr.nominal}) are used, where \code{n.rej.stop} is set to 20. 
#' @param seed a single-value integer seed for the random process of drawing permutation replicates. 
#' The seed is user supplied or internally generated. The default is 123.
#' @return a list consisting of 
#'   \item{p.linear}{a \code{n.otu} by \code{n.otu} matrix of p-values for the linear test}
#'   \item{q.linear}{a matrix of q-values for the linear test}
#'   \item{p.nonlinear}{a matrix of p-values for the nonlinear test}
#'   \item{q.nonlinear}{a matrix of q-values for the nonlinear test}
#'   \item{p.omni}{a matrix of p-values for the omnibus test}
#'   \item{q.omni}{a matrix of q-values for the omnibus test}
#'   \item{which.pmin}{a \code{n.otu} by \code{n.otu} matrix of 0, 1, and 2 values, where 0 and 1 indicate that the nonlinear and linear tests, 
#' respectively, achieved the minimum p-value between the two, and 2 indicates that both tests yielded similar p-values.}
#' @keywords microbiome, network
#' @author Yi-Juan Hu <yijuan.hu@emory.edu>
#' @importFrom permute shuffleSet
#' @importFrom stats p.adjust
#' @import matrixStats
#' @export
#' @references Su C, He M, Van Doren VE, Kelley CF, Hu YJ (2024). A general testing method for inference of microbial networks with compositional data.
#'   XXX.
#' @examples
#' TestNet.res <- TestNet(otu.tab = sim.otu.tab)
TestNet <- function(otu.tab, clustered.data = FALSE, CTRL = NULL, fdr.nominal = 0.1, n.perm.max = NULL, seed = 123) {
    
    n <- nrow(otu.tab)
    J1 <- ncol(otu.tab)
    absence <- 1*(otu.tab==0)
    
    #----------------------
    # observed statistics
    #----------------------
    
    P <- (otu.tab + 0.5)/rowSums(otu.tab + 0.5)
    logP <- log(P)
    clrP <- logP - rowMeans(logP)
    impute <- colMeans(ifelse(absence, clrP, NA), na.rm=TRUE)
    clrP_mod <- t(ifelse(t(absence), impute, t(clrP))) # de-correlated
    
    cov_clrP_l1 <- crossprod(t(t(clrP_mod) - colMeans(clrP_mod)) ) / (n-1)
    cov_clrP_l0 <- matrix(NA, nrow=J1, ncol=J1)
    for (k in 1:(J1-1)) for (l in (k+1):J1) cov_clrP_l0[l,k] <- dcov::dcov2d(clrP_mod[,k], clrP_mod[,l]) 
    
    otu.l1 = cov_clrP_l1[lower.tri(cov_clrP_l1)] # 'otu' = 'pair' = 'test'
    otu.l0 = cov_clrP_l0[lower.tri(cov_clrP_l0)]
    
    #------------------------
    # permutation statiatics
    #------------------------
    
    n.otu = J1*(J1-1)/2
    inv.J1 = 1/J1
    correction = sqrt(J1)/sqrt(J1-1)
    
    n.rej.stop = 20
    n.perm.block = 100
    n.perm.minP = 100   
    tol.eq = 10^-6
    if (is.null(n.perm.max)) n.perm.max = n.otu * n.rej.stop * (1/fdr.nominal)
    
    n.otu.smallp = n.otu
    otu.smallp = 1:n.otu
    otu.smallp.2d = combn(seq_len(J1), 2)

    perm_absence = matrix(0, nrow=n, ncol=J1)
    perm_clrP = matrix(0, nrow=n, ncol=J1)
    otu.l1.perm = matrix(NA, nrow=n.otu, ncol=n.perm.minP)
    otu.l0.perm = matrix(NA, nrow=n.otu, ncol=n.perm.minP)
    Aset.l1 = rep(TRUE, n.otu)
    Aset.l0 = rep(TRUE, n.otu)
    Aset.omni = rep(TRUE, n.otu)
    p.otu.l1 = rep(NA, n.otu)
    p.otu.l0 = rep(NA, n.otu)
    p.otu.omni = rep(NA, n.otu)
    n.otu.l1 = rep(0, n.otu)
    n.otu.l0 = rep(0, n.otu)
    n.otu.l1.left = rep(0, n.otu)
    n.otu.l1.right = rep(0, n.otu)
    which.pmin.otu = rep(NA, n.otu)
    n.perm.completed = 0
    otu.tests.stopped = FALSE
    
    set.seed(seed)
    
    if (clustered.data) {
        ttl <- J1*1000
        perm.ttl <- t(shuffleSet(n, ttl, CTRL))
    } else {
        ttl <- J1*1000
        perm.ttl <- t(shuffleSet(n, ttl))
    }
    
    #------------------------------------------------
    # permutation
    #------------------------------------------------
    
    for (i.sim in 1:n.perm.max) {
        
        if (i.sim == 10001) {
            n.perm.block = 1000
        }
        
        perm <- perm.ttl[,sample(1:ttl, J1, replace=TRUE)]
        for (j in 1:J1) {
            perm_absence[,j] <- absence[perm[,j],j]
            perm_clrP[,j] <- clrP[perm[,j],j]
        }
            
        perm_clrP_tmp <- sweep(perm_clrP, 1, rowMeans(perm_clrP)) * correction
        perm_impute <- colMeans(ifelse(perm_absence, perm_clrP_tmp, NA), na.rm=TRUE)
        perm_clrP_mod <- t(ifelse(t(perm_absence), perm_impute, t(perm_clrP_tmp))) # de-correlated
        
        n.perm.completed = n.perm.completed + 1
        
        #------------------------------------------------
        # reduce OTUs to speed-up computation
        #------------------------------------------------
        
        if (i.sim >= n.perm.minP+1 & i.sim %% n.perm.block == 1) { # 1001, 2001, 3001, ...
            
            w.otu.smallp = which(apply(rbind(Aset.l1, Aset.l0, Aset.omni), 2, any))
            
            cat("number of tests do not meet early stopping criterion:", length(w.otu.smallp), "\n")
            
            if (length(w.otu.smallp) < n.otu.smallp) {
            
                otu.smallp = otu.smallp[w.otu.smallp]
                otu.smallp.2d = otu.smallp.2d[,w.otu.smallp,drop=FALSE]
                n.otu.smallp = length(otu.smallp)
                inv.n.otu.smallp = 1/n.otu.smallp
                
                n.otu.l1.left = n.otu.l1.left[w.otu.smallp]
                n.otu.l1.right = n.otu.l1.right[w.otu.smallp]
                n.otu.l1 = n.otu.l1[w.otu.smallp]
                n.otu.l0 = n.otu.l0[w.otu.smallp]
                Aset.l1 = Aset.l1[w.otu.smallp]
                Aset.l0 = Aset.l0[w.otu.smallp]
                Aset.omni = Aset.omni[w.otu.smallp]
            }
            
            otu.l1.perm.smallp = otu.l1.perm[w.otu.smallp,1:(i.sim-1),drop=FALSE]
            otu.l1.perm = matrix(NA, nrow=n.otu.smallp, ncol=i.sim-1+n.perm.block)
            otu.l1.perm[,1:(i.sim-1)] = otu.l1.perm.smallp
            rm(otu.l1.perm.smallp)
            
            otu.l0.perm.smallp = otu.l0.perm[w.otu.smallp,1:(i.sim-1),drop=FALSE]
            otu.l0.perm = matrix(NA, nrow=n.otu.smallp, ncol=i.sim-1+n.perm.block)
            otu.l0.perm[,1:(i.sim-1)] = otu.l0.perm.smallp
            rm(otu.l0.perm.smallp)
        }
        
        #------------------------------------------------
        # permutation statistics, p-values
        #------------------------------------------------

        for (k in unique(otu.smallp.2d[1,])) {  
            w_k <- which(otu.smallp.2d[1,]==k)
            l_seq <- otu.smallp.2d[2, w_k]
            len_l_seq <- length(l_seq)
            
            perm_clrP_modc_k <- perm_clrP_mod[,k] - mean(perm_clrP_mod[,k]) 
            perm_clrP_modc_l <- t(t(perm_clrP_mod[,l_seq,drop=FALSE]) - colMeans(perm_clrP_mod[,l_seq,drop=FALSE])) 
            otu.l1.perm[w_k, i.sim] <- crossprod(perm_clrP_modc_k, perm_clrP_modc_l) / (n-1)
            
            for (l in 1:len_l_seq) {
                otu.l0.perm[w_k[l], i.sim] <- dcov::dcov2d(perm_clrP_mod[,k], perm_clrP_mod[,l_seq[l]])
            }
        }
        
        # p-value
        diff.l1 <- otu.l1.perm[,i.sim] - c(otu.l1[otu.smallp])
        n.otu.l1.left <- n.otu.l1.left + (diff.l1 < -tol.eq) + 0.5*(abs(diff.l1) <= tol.eq)
        n.otu.l1.right <- n.otu.l1.right + (diff.l1 > tol.eq) + 0.5*(abs(diff.l1) <= tol.eq)
        n.otu.l1 <- pmin(n.otu.l1.left, n.otu.l1.right)*2

        diff.l0 <- otu.l0.perm[,i.sim] - c(otu.l0[otu.smallp])
        n.otu.l0 = n.otu.l0 + (diff.l0 >= tol.eq) + 0.5*(abs(diff.l0) < tol.eq)

        #--------------------------------------------------
        # check if the sequential procedure can be stopped
        #--------------------------------------------------
        
        if ((i.sim >= n.perm.minP & i.sim %% n.perm.block == 0) | (i.sim == n.perm.max)) {
            
            cat("permutations:", i.sim, "\n")
            
            inv.n.perm.completed = 1/n.perm.completed
            inv.n.perm.completed.1 = 1/(n.perm.completed+1)

            if (any(Aset.l1)) {
                AtoB.l1 <- Aset.l1 & (n.otu.l1 >= n.rej.stop)
                Aset.l1 <- Aset.l1 & !AtoB.l1
                p.otu.l1[otu.smallp][AtoB.l1] <- n.otu.l1[AtoB.l1]*inv.n.perm.completed
                p.otu.l1[otu.smallp][Aset.l1] <- (n.otu.l1[Aset.l1]+1)*inv.n.perm.completed.1
                q.otu.l1 <- fdr.Sandve(p.otu.l1) 
                if (all(q.otu.l1[otu.smallp][Aset.l1] <= fdr.nominal)) Aset.l1[Aset.l1] <- FALSE
            }
            
            if (any(Aset.l0)) {
                AtoB.l0 <- Aset.l0 & (n.otu.l0 >= n.rej.stop)
                Aset.l0 <- Aset.l0 & !AtoB.l0
                p.otu.l0[otu.smallp][AtoB.l0] <- n.otu.l0[AtoB.l0]*inv.n.perm.completed
                p.otu.l0[otu.smallp][Aset.l0] <- (n.otu.l0[Aset.l0]+1)*inv.n.perm.completed.1
                q.otu.l0 <- fdr.Sandve(p.otu.l0) 
                if (all(q.otu.l0[otu.smallp][Aset.l0] <= fdr.nominal)) Aset.l0[Aset.l0] <- FALSE
            }
            
            if (any(Aset.omni)) { 
                pmin.otu.omni <- pmin(n.otu.l1, n.otu.l0)
                which.pmin.otu[otu.smallp] <- 1*(pmin.otu.omni == n.otu.l1 & pmin.otu.omni != n.otu.l0) + 2*(n.otu.l1 == n.otu.l0)
                
                pnull.otu.l1 <- n.perm.completed + 0.5 - rowRanks(abs(otu.l1.perm[,1:n.perm.completed,drop=FALSE]), ties.method="average")
                pnull.otu.l0 <- n.perm.completed + 0.5 - rowRanks(otu.l0.perm[,1:n.perm.completed,drop=FALSE], ties.method="average")
                pnullmin.otu.omni <- pmin(pnull.otu.l1, pnull.otu.l0)
                
                n.otu.omni <- rowSums( (pnullmin.otu.omni < c(pmin.otu.omni) - tol.eq)) + 0.5 * rowSums( (abs(pnullmin.otu.omni - c(pmin.otu.omni)) < tol.eq))
                    
                AtoB.omni <- Aset.omni & (n.otu.omni >= n.rej.stop)
                Aset.omni <- Aset.omni & !AtoB.omni
                    
                p.otu.omni[otu.smallp][AtoB.omni] <- n.otu.omni[AtoB.omni]*inv.n.perm.completed
                p.otu.omni[otu.smallp][Aset.omni] <- (n.otu.omni[Aset.omni]+1)*inv.n.perm.completed.1
                    
                q.otu.omni <- fdr.Sandve(p.otu.omni) 
                if (all(q.otu.omni[otu.smallp][Aset.omni] <= fdr.nominal)) Aset.omni[Aset.omni] <- FALSE
            }
                
            if (!any(Aset.l1) & !any(Aset.l0) & !any(Aset.omni)) {
                otu.tests.stopped = TRUE 
                cat("otu test stopped at permutation", i.sim, "\n")
                break
            }
        } 
        
    } # i.sim in 1:n.perm.max
    
    p.linear <- matrix(NA, J1, J1)
    q.linear <- matrix(NA, J1, J1)
    p.linear[lower.tri(p.linear)] <- p.otu.l1
    q.linear[lower.tri(q.linear)] <- p.adjust(p.otu.l1, method = "BH")
    rownames(p.linear) <- colnames(otu.tab)
    colnames(p.linear) <- rownames(p.linear)
    rownames(q.linear) <- rownames(p.linear)
    colnames(q.linear) <- colnames(p.linear)
    
    p.nonlinear <- matrix(NA, J1, J1)
    q.nonlinear <- matrix(NA, J1, J1)
    p.nonlinear[lower.tri(p.nonlinear)] <- p.otu.l0
    q.nonlinear[lower.tri(q.nonlinear)] <- p.adjust(p.otu.l0, method = "BH")
    rownames(p.nonlinear) <- colnames(otu.tab)
    colnames(p.nonlinear) <- rownames(p.nonlinear)
    rownames(q.nonlinear) <- rownames(p.nonlinear)
    colnames(q.nonlinear) <- colnames(p.nonlinear)
    
    p.omni <- matrix(NA, J1, J1)
    q.omni <- matrix(NA, J1, J1)
    p.omni[lower.tri(p.omni)] <- p.otu.omni
    q.omni[lower.tri(q.omni)] <- p.adjust(p.otu.omni, method = "BH")
    rownames(p.omni) <- colnames(otu.tab)
    colnames(p.omni) <- rownames(p.omni)
    rownames(q.omni) <- rownames(p.omni)
    colnames(q.omni) <- colnames(p.omni)
    
    which.pmin <- matrix(NA, J1, J1)
    which.pmin[lower.tri(which.pmin)] <- which.pmin.otu
    
    return(list(p.linear=p.linear, q.linear=q.linear, 
                p.nonlinear=p.nonlinear, q.nonlinear=q.nonlinear, 
                p.omni=p.omni, q.omni=q.omni, 
                which.pmin=which.pmin)) 
    
} # TestNet


fdr.Sandve = function(p.otu) {
    
    m = length(p.otu)
    
    p.otu.sort = sort(p.otu)
    n.otu.detected = seq(1, m)
    pi0 = min(1, 2/m*sum(p.otu))
    
    qval.sort = m * pi0 * p.otu.sort / n.otu.detected
    j.min.q = 1
    while (j.min.q < m) {
        min.q = min( qval.sort[j.min.q:m] )
        new.j.min.q = (j.min.q-1) + max( which(qval.sort[j.min.q:m]==min.q) )
        qval.sort[j.min.q:new.j.min.q] = qval.sort[new.j.min.q]
        j.min.q = new.j.min.q+1
    }
    mat = match(p.otu, p.otu.sort)   
    qval.orig = qval.sort[mat]
    results = qval.orig
    return(results)
    
} # fdr.Sandve


