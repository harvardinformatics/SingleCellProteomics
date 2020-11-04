normalizeData <- function(mss, byrun=TRUE, vsn=TRUE) {
    require(MSnbase)

    if (byrun) {
        channels <- do.call(rbind, strsplit(sampleNames(mss), split='_'))[, 1]
        xlst <- sapply(channels, function(x) {
            xms <- mss[, grep(paste0(x, '_'), sampleNames(mss))]
            xms <- normalise(xms, 'quantiles',)
            return(xms)
        })
        mss <- Reduce(function(x, y) combine(x, y), xlst)
    }

    if (vsn) {
        mss <- normalise(mss, 'vsn')
    }
    
    return(mss)
}

correctBatch <- function(mss) {
    require(sva)
    
    pheno <- pData(mss)
    batch <- pheno$batch
    edata <- exprs(mss)
    combat.dat <- ComBat(edata, batch)
    exprs(mss) <- combat.dat
    
    return(mss)
}

rmCompPCA <- function(mss, rmcomp=1) {
    rmcomp <- as.integer(rmcomp)
    X <- exprs(mss)
    
    Xpca <- prcomp(X)
    loadings <- Xpca$rotation
    scores <- Xpca$x

    mu <- colMeans(X)
    Xhat = scores[,seq(ncol(scores))[-rmcomp]] %*% t(loadings[,seq(ncol(loadings))[-rmcomp]])
    Xhat = scale(Xhat, center = -mu, scale = FALSE)
    exprs(mss) <- Xhat
    
    return(mss)
}

removeVar <- function(mss, numf) {
    numf <- as.integer(numf)

    dsgn <- model.matrix(~ 0+factor(pData(mss)$TreatmentGroup))
    fit <- lmFit(mss, dsgn)
    residuals.m <- residuals.MArrayLM(fit, exprs(mss))
    e <- exprs(mss)
    te <- t(e)
    E <- t(residuals.m)
    svdWa <- svd(E)
    W <- svdWa$u[, numf]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% te
    cY <- te - W %*% alpha
    cY <- t(cY)
    exprs(mss) <- cY
    
    return(mss)
}
