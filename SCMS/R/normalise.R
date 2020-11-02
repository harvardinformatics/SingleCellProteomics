normalizeSingleRuns <- function(mss, channels, type='quantiles') {
    require(MSnbase)
    
    xlst <- sapply(channels, function(x) {
        xms <- mss[, grep(x, sampleNames(mss))]
        xms <- normalise(xms, type)
        return(xms)
    })
    
    nmss <- Reduce(function(x, y) combine(x, y), xlst)
    return(nmss)
}

doNormAndBatch <- function(mss) {
    require(sva)
    
    mss.vsn <- normalise(mss, 'vsn')
    
    pd <- phenoData(mss.vsn)$TreatmentGroup
    names(pd) <- sampleNames(mss.vsn)
    pd.df <- as.data.frame(pd)
    batch <- as.integer(sub('F', '', sub('_.*', '', sampleNames(mss.vsn))))
    pd.df$batch  <- batch

    res.cb <- ComBat(exprs(mss.vsn), pd.df$batch)
    res.cb <- rmCompPCA(res.cb, 1)
    
    return(list(pd, res.cb))
}

rmCompPCA <- function(X, rm=1) {
    rm <- as.integer(rm)
    
    mu <- colMeans(X)

    Xpca <- prcomp(X)
    rot <- Xpca$rotation
    retx <- Xpca$x

    Xhat = retx[,seq(ncol(retx))[-rm]] %*% t(rot[,seq(ncol(rot))[-rm]])
    Xhat = scale(Xhat, center = -mu, scale = FALSE)
    
    return(Xhat)
}

removeVar <- function(mses, fnum) {
    fnum <- as.integer(fnum)

    dsgn <- model.matrix(~ 0+factor(pData(mses)$TreatmentGroup))
    colnames(dsgn) <- c('A', 'B')

    fit <- lmFit(mses, dsgn)
    residuals.m <- residuals.MArrayLM(fit, exprs(mses))
    e <- exprs(mses)
    te <- t(e)
    E <- t(residuals.m)
    svdWa <- svd(E)
    W <- svdWa$u[, fnum]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% te
    cY <- te - W %*% alpha
    cY <- t(cY)

    return(cY)
}
