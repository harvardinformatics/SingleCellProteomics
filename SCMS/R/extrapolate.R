extrapolateMissing <- function(df, bgd) {
    # df is ctrl or treatment partial data.frame
    if (!anyNA(df)) {
        return(df)
    }

    while (TRUE) {
        if (nrow(df) == 1) {
            df <- fillRowColumn(df, bgd)
            break
        }

        startNA <- sum(is.na(df))
        df <- fillRowColumn(df, bgd)
        df <- fillRowColumn(t(df), bgd)
        df <- t(df)
        endNA <- sum(is.na(df))
        
        if (startNA == endNA) {
            break
        }
        
    }
    return(df)
}
    
fillRowColumn <- function(df, bgd) {
    res <- apply(df, 1, function(v) {
        if (sum(!is.na(v)) > 1) {
            if (sum(is.na(v)) == 1) {
                v[which(is.na(v))] <- mean(v, na.rm=TRUE)
                v
            } else {
                sampleValuesAndReplace(v)
            }
        } else {
            sampleBaselineAndReplace(v, bgd)
        }
    })
    return(t(res))
}

sampleValuesAndReplace <- function(x) {
    numbermissing <- sum(is.na(x))
    mn <- mean(x[!is.na(x)])
    stdev <- mn * 0.1
    x[which(is.na(x))] <- rnorm(numbermissing, mean=mn, sd=stdev)
    return(x)
}

sampleBaselineAndReplace <- function(x, bgd) {
    numbermissing <- sum(is.na(x))
    if (all(is.na(x))) {
        x[which(is.na(x))] <- runif(numbermissing, min=bgd*0.75, max=bgd)
        x
    } else {
        anchor <- min(x, na.rm=TRUE)
        x[which(is.na(x))] <- runif(numbermissing, min=bgd*0.9, max=anchor)
        x
    }
}
