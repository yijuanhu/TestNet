options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args) 

if (length(args)==0) {
    print("No arguments supplied")
    
    n_subj <- 100 # 50, 100, 200, 500, 1000
    J <- 100 # 50, 100, 200
    
    clustered_data <- 0 # 0, 1
    abs_dist <- "Normal" # Normal, Gamma
    cov_model <- "Identity" # Identity, AR1pos, AR1, AR4, Hub, Block, Sparse, Ushape
    linear <- 1 # 1, 0
    
    lib_mean <- 10000  # 10000, 30000
    filter <- "stringent" # lenient, stringent
    
    n_sim <- 100
    i_seed <- 1
    
} else {
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }
}

filename <- paste("C", clustered_data, "_", abs_dist, "_nsubj", n_subj, "_J", J, "_", cov_model, "_ln", linear, "_", lib_mean, "_", filter, "_seed", i_seed, ".txt", sep="")

library(TestNet)
library(dplyr)
library(PRROC) # library(pROC)
library(MASS) # rnegbin
library(dirmult) # simPop


set.seed(123)
mu <- sort(runif(J, 0, 10), decreasing=TRUE)

sim_upper_tri <- function(size, prob, rho) {
    n_vec <- size*(size-1)/2
    n_choose <- floor(n_vec*prob)
    vec <- rep(0, n_vec)
    choose <- sort(sample(1:n_vec, n_choose))
    vec[choose] <- sample(c(rho, -rho), n_choose, replace=TRUE)
    mat <- diag(size)
    mat[upper.tri(mat)] <- vec
    return(mat)
}


metric <- function(R_hat, R0, q.new, thresh = 0.3, stat = "rho", fdr_nominal = 0.1) {
    
    R0_vec <- R0[lower.tri(R0)]
    R0_bin <- (abs(R0_vec) > 0.0)
    R_hat_vec <- R_hat[lower.tri(R_hat)]
    if (stat=="rho") {
        R_hat_bin <- (abs(R_hat_vec) > thresh)
        FDR <- sum(R_hat_bin==TRUE & R0_bin==FALSE, na.rm=TRUE)/max(sum(R_hat_bin==TRUE, na.rm=TRUE), 1)
        sen <- NA
        n <- sum(R_hat_bin)
    } else {
        R_hat_bin <- (abs(R_hat_vec) < thresh)  # changed from <= to <
        q.new_vec <- q.new[lower.tri(q.new)]
        q.new_bin <- (q.new_vec < fdr_nominal)  # changed from <= to <
        FDR <- sum(q.new_bin==TRUE & R0_bin==FALSE, na.rm=TRUE)/max(sum(q.new_bin==TRUE, na.rm=TRUE), 1)
        sen <- ifelse(sum(R0_bin)==0, NA, sum(q.new_bin==TRUE & R0_bin==TRUE, na.rm=TRUE)/max(sum(R0_bin==TRUE, na.rm=TRUE), 1))
        n <- sum(q.new_bin)
    }
    FPR <- sum(R_hat_bin==TRUE & R0_bin==FALSE, na.rm=TRUE)/sum(R0_bin==FALSE, na.rm=TRUE)
    TPR <- ifelse(sum(R0_bin)==0, NA, sum(R_hat_bin==TRUE & R0_bin==TRUE, na.rm=TRUE)/sum(R0_bin==TRUE, na.rm=TRUE))
    AUC <- ifelse(sum(R0_bin)==0, NA, roc.curve(scores.class0 = abs(R_hat_vec)[R0_bin==0], scores.class1 = abs(R_hat_vec)[R0_bin==1])$auc)
    AUC <- ifelse(AUC < 0.5, 1-AUC, AUC)
    PRAUC <- ifelse(sum(R0_bin)==0, NA, pr.curve(scores.class0 = abs(R_hat_vec)[R0_bin==0], scores.class1 = abs(R_hat_vec)[R0_bin==1])$auc.integral)
    return(c(FPR, TPR, FDR, sen, n, AUC, PRAUC))
}


if (cov_model == "Identity") {
    
    seed.var <- 1
    var_cov <- diag(J) 
    
} else if (cov_model == "AR1pos") {

    seed.var <- 4
    set.seed(seed.var)
    var_cov <- diag(J)
    for (j in 1:(J-1)) var_cov[j,j+1] <- 0.5
    var_cov[lower.tri(var_cov)] <- t(var_cov)[lower.tri(var_cov)]

} else if (cov_model == "AR1") {
    
    seed.var <- 4
    set.seed(seed.var)
    var_cov <- diag(J)
    for (j in 1:(J-1)) var_cov[j,j+1] <- sample(c(0.5, -0.5), 1)
    var_cov[lower.tri(var_cov)] <- t(var_cov)[lower.tri(var_cov)] 
    
} else if (cov_model == "AR4") {
    
    seed.var <- 1
    set.seed(seed.var)
    var_cov <- diag(J)
    for (j in 1:(J-1)) var_cov[j,j+1] <- 0.4
    for (j in 1:(J-2)) var_cov[j,j+2] <- 0.2
    for (j in 1:(J-3)) var_cov[j,j+3] <- 0.2
    for (j in 1:(J-4)) var_cov[j,j+4] <- 0.1
    sign <- sample(c(1, -1), size=J, replace=TRUE)
    var_cov <- var_cov*sign
    var_cov[lower.tri(var_cov)] <- t(var_cov)[lower.tri(var_cov)]
    diag(var_cov) <- 1.5

} else if (cov_model == "Hub") {
    
    seed.var <- 1
    set.seed(seed.var)
    n_hub <- 3
    n_nonhub <- J - n_hub
    hub <- sort(sample(1:20, n_hub))
    nonhub <- setdiff(1:J, hub)
    var_cov <- diag(J)
    var_cov <- sim_upper_tri(size=J, prob=0.7, rho=0.3)
    var_cov[nonhub, nonhub] <- sim_upper_tri(size=n_nonhub, prob=0.2, rho=0.3)
    var_cov <- var_cov + t(var_cov)
    if (J==50) diag(var_cov) <- 2.5
    if (J==100) diag(var_cov) <- 3.1
    if (J==200) diag(var_cov) <- 4.4
    
} else if (cov_model == "Block") {
    
    seed.var <- 4
    set.seed(seed.var)
    n.block <- 10
    J.perblock <- J/n.block
    remain <- 1:J
    var_cov <- diag(J)
    var_cov <- sim_upper_tri(size=J, prob=0.2, rho=0.3)
    for (k in 1:n.block) {
        block <- sort(sample(remain, J.perblock))
        var_cov[block, block] <- sim_upper_tri(size=J.perblock, prob=0.5, rho=0.3)
        remain <- setdiff(remain, block)
    }
    var_cov <- var_cov + t(var_cov)
    if (J==100) diag(var_cov) <- 2.8
    if (J==200) diag(var_cov) <- 4.2
    
} else if (cov_model == "Sparse") {
    
    seed.var <- 3
    set.seed(seed.var)
    p1 <- floor(sqrt(J))   
    n.lower.tri <- (p1-1)*p1/2
    choose <- rbinom(n.lower.tri, 1, 0.3) 
    unif <- runif(n.lower.tri)
    unif1 <- ifelse(unif<0.5, unif-1, unif)
    B <- matrix(0, nrow=p1, ncol=p1)
    B[lower.tri(B)] <- choose*unif1
    B <- B + t(B)
    epsilon <- max(-min(eigen(B)$values), 0) + 0.01
    A1 <- B + epsilon*diag(p1)
    A2 <- 4*diag(J-p1)
    var_cov <- diag(J)
    var_cov[1:p1,1:p1] <- A1
    var_cov[(p1+1):J, (p1+1):J] <- A2
    
} else if (cov_model == "Ushape") {
    
    seed.var <- 1
    var_cov <- diag(J)
    for (j in seq(2, J, 2)) {
        var_cov[j-1, j] <- 1
        var_cov[j, j-1] <- 1
    }
}

if (cov_model != "Ushape") {
    eigen_res <- eigen(var_cov)
    w <- which(colSums(eigen_res$vector) < 0)
    eigen_res$vector[,w] <- -eigen_res$vector[,w]
    F_mat <- eigen_res$vector %*% diag(sqrt(eigen_res$values))
}

#-------------
# simulation
#-------------

for (sim in 1:n_sim) {
    
    set.seed(i_seed*1000+sim)
    
    # Z (absolute abundance) => Y (true relative abundance)
    
    const <- 1/sqrt(20) 
    
    if (!clustered_data) {
        
        Z <- matrix(NA, nrow = n_subj, ncol = J)
        
        if (cov_model != "Ushape") {
            for (i in 1:n_subj) Z[i, ] <- mu + F_mat %*% ifelse(rep(abs_dist, J)=="Normal", rnorm(J, mean=0, sd=1), rgamma(J, shape=10, scale=1) * const)
        } else if (cov_model=="Ushape") { # nonlinear
            for (i in 1:n_subj) Z[i, ] <- mu + ifelse(rep(abs_dist, J)=="Normal", rnorm(J, mean=0, sd=1), rgamma(J, shape=10, scale=1) * const)
            for (j in seq(2, J, 2)) {
                raw <- poly(x = Z[, j-1], degree = 2, raw = FALSE)[, 2] + rnorm(n_subj, mean=0, sd=0.05)
                Z[, j] <- raw / sd(raw) + mu[j]
            }
        } 
        
        n_sam <- n_subj
        
    } else if (clustered_data) {
        
        cluster.size <- 4
        cluster.id <- rep(1:n_subj, each=cluster.size)
        
        Z <- matrix(NA, nrow = n_subj*cluster.size, ncol = J)
        
        if (cov_model != "Ushape") {
            for (i in 1:n_subj) {
                mu_i <- mu + rnorm(J, mean=0, sd=1) # sd is the hetero parameter
                for (s in (cluster.size-1):0) {
                    Z[cluster.size*i-s, ] <- mu_i +  F_mat %*% ifelse(rep(abs_dist, J)=="Normal", rnorm(J, mean=0, sd=1), rgamma(J, shape=10, scale=1)*const)
                }
                s <- (cluster.size*(i-1)+1) : (cluster.size*i)
                o <- order(Z[s,1])
                Z[s,] <- Z[s[o],]
            }
        } else if (cov_model=="Ushape") { # nonlinear
            
            for (i in 1:n_subj) {
                mu_i <- mu + rnorm(J, mean=0, sd=1) # sd is the hetero parameter
                for (s in (cluster.size-1):0) {
                    Z[cluster.size*i-s, ] <- mu_i +  ifelse(rep(abs_dist, J)=="Normal", rnorm(J, mean=0, sd=1), rgamma(J, shape=10, scale=1)*const)
                }
                s <- (cluster.size*(i-1)+1) : (cluster.size*i)
                o <- order(Z[s,1])
                Z[s,] <- Z[s[o],]
            }
            for (j in seq(2, J, 2)) {
                raw <- poly(x = Z[, j-1], degree = 2, raw = FALSE)[, 2] + rnorm(n_subj*cluster.size, mean=0, sd=0.05)
                Z[, j] <- raw / sd(raw) + mu_i[j]
            }
        } 
        
        n_sam <- cluster.size * n_subj
    }
    
    Z <- exp(Z)
    colnames(Z) <- paste("taxon", 1:J, sep="")
    Y <- Z / rowSums(Z)
    
    
    # add experimental bias to Z
    
    biased.data <- FALSE
    if (biased.data) {
        bias.factor <- exp(rnorm(J, 0, 0.5))
        Z.biased <- t(t(Z) * bias.factor)
        Y.biased <- Z.biased / rowSums(Z.biased)
    } else {
        Y.biased <- Z / rowSums(Z)
    }
    
    
    # X (observed count)
    
    lib.size <- round(rnorm(n=n_sam, mean=lib_mean, sd=lib_mean/3))
    lib.size <- ifelse(lib.size < 2000, 2000, lib.size)
    X <- matrix(NA, nrow = n_sam, ncol = J)
    
    theta <- 1
    for (j in 1:J) X[,j] <- rnegbin(n=n_sam, mu = lib.size*Y.biased[,j], theta = theta) # var = u+u^2/theta 

    colnames(X) <- paste("taxon", 1:J, sep="")
    sam_ID <- paste("sam", 1:n_sam, sep="")
    rownames(X) <- sam_ID
    metadata <- data.frame(sam_ID = sam_ID)
    rownames(metadata) <- sam_ID
    
    
    # P (observed relative abundance)
    # filtering and preparing statistics
    
    if (filter == "stringent") {
        P <- (X + 0.5)/rowSums(X + 0.5)
        P_mean <- colMeans(P)
        w <- which(P_mean >= 0.001) 
    } else {
        w <- which(colSums(X>0) >= n_sam*0.05) 
    }
    (J1 <- length(w))

    var_cov_filter <- var_cov[w, w]
    X_filter <- X[,w]

    #-------------------------------------------------
    # TestNet
    #-------------------------------------------------
    
    if (!clustered_data) {
        res.testnet <- TestNet(otu.tab = X_filter)
    } else if (clustered_data) {
        res.testnet <- TestNet(otu.tab = X_filter, clustered_data = TRUE, cluster.id = cluster.id)
    }


    ##############
    # summarize
    ##############

    res.sum <- metric(res.testnet$p.omni, var_cov_filter, res.testnet$q.omni, stat = "pvalue", thresh = 0.005)
    write.table(t(signif(res.sum, 3)), filename, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    
} # simulation iteration

