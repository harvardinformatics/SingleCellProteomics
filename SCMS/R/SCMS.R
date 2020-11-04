prepareRunList <- function(dat) {
    xdf <- filterIsolationInterference(dat)
    xdf <- filterChannels(xdf, 126, 127)
    xlst <- partitionMSruns(xdf)

    return(xlst)
}

doExtrapolateMissing <- function(df, bgdL, bgdR) {
    dfL <- df[-grep('ProteinAccessions', colnames(df))][, 1:4]
    dfR <- df[-grep('ProteinAccessions', colnames(df))][, 5:8]
    resL <- extrapolateMissing(dfL, bgdL)
    resR <- extrapolateMissing(dfR, bgdR)
    res <- cbind(resL, resR)
    cbind(df['ProteinAccessions'], res)
}

processSCMS <- function(dat) { # processSCMS(xlst) list of runs
    runlst <- prepareRunList(dat)
    reslst <- lapply(runlst, function(run) {
        backgroundL <- min(run[, 2:5], na.rm=TRUE)
        backgroundR <- min(run[, 6:9], na.rm=TRUE)
        runPSMlst <- partitionProteins(run)
        resPSMlst <- lapply(runPSMlst, function(PSMdf) doExtrapolateMissing(PSMdf, backgroundL, backgroundR))
        resPSM <- Reduce(rbind, resPSMlst)
        resProt <- aggregate(resPSM[-grep('ProteinAccessions', colnames(resPSM))], resPSM[grep('ProteinAccessions', colnames(resPSM))], median)
    })
    res <- Reduce(function(x, y) merge(x, y, by='ProteinAccessions', all=FALSE), reslst)
    rownames(res) <- res$ProteinAccessions
    res <- res[-grep('ProteinAccessions', colnames(res))]
    res
}

makeMSnSet <- function(df) {
    require(MSnbase)
    
    features.df <- data.frame(ID=rownames(df), Acc=rownames(df))
    rownames(features.df) <- features.df$ID

    group <- sapply(colnames(df), function(nm) ifelse(grepl(paste(c('128', '129'), collapse='|'), nm), '0', '1'))
    batch <- do.call(rbind, strsplit(colnames(df), split='_'))[, 1]
    phenoData <- data.frame(TreatmentGroup=group, batch=batch)

    MSnSet(as.matrix(df), features.df, phenoData)
}
