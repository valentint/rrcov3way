##
##  X is an n x m x p array
##
##  - the function ALS() expects a n x m*p matrix
##
t3_als <- function(X, P=2, Q=2, R=2, conv=1e-10, crit=0.975)
{
    di <- dim(X)
    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X)

    Y <- unfold(X)     # default mode is "A", returns n x m*p matrix
    tt <- ThreeWay::T3funcrep(Y, n, m, p, P, Q, R, start=0, conv=conv)

    ## Compute orthogonal distances; flag observations as outliers
    Xfit <- tt$A %*% tt$H %*% t(kronecker(tt$C, tt$B))
    out.Xhat <- array(Xfit, c(n, m, p))
    odsq <- apply((Y-Xfit)^2, 1, sum)
    RD <- sqrt(odsq)
    critRD <- .cutoff.rd(RD, crit=crit, robust=FALSE)   # (mean(RD^(2/3)) + sd(RD^(2/3)) * qnorm(crit))^(3/2)
    flag <- RD <= critRD

    ## dimnames back
    dimnames(tt$A) <- list(dn[[1]], paste("F", 1:P, sep=""))
    dimnames(tt$B) <- list(dn[[2]], paste("F", 1:Q, sep=""))
    dimnames(tt$C) <- list(dn[[3]], paste("F", 1:R, sep=""))
    dimnames(tt$H) <- list(paste("F", 1:P, sep=""), paste("F", 1:(Q*R), sep=""))
    names(RD) <- dn[[1]]
    names(flag) <- dn[[1]]

    ret <- list(fit=tt$f, fp=tt$fp,
                A=tt$A, B=tt$B, C=tt$C, GA=tt$H,
                La=tt$La, Lb=tt$Lb, Lc=tt$Lc,
                iter=tt$iter, flag=flag, rd=RD, cutoff.rd=critRD)
}
