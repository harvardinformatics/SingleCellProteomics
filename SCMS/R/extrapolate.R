extrapolateMissing <- function(df) {
    # df is ctrl or treatment partial data.frame
    if (!anyNA(df)) {
        return(df)
    }

    while (TRUE) {
        if (nrow(df) == 1) {
            df <- fillRowColumn(df)
            break
        }

        startNA <- sum(is.na(df))
        df <- fillRowColumn(df)
        df <- fillRowColumn(t(df))
        df <- t(df)
        endNA <- sum(is.na(df))
        
        if (startNA == endNA) {
            break
        }
        
    }
    return(df)
}

fillRowColumn <- function(df) {
    res <- apply(df, 1, function(v) {
        if (sum(is.na(v)) == 1) {
            v[which(is.na(v))] <- mean(v, na.rm=TRUE)
            v
        } else if (sum(is.na(v)) <= 2) {
            sampleAndReplace(v)
        } else {
            v
        }
    })
    return(t(res))
}

sampleAndReplace <- function(x) {
    numbermissing <- sum(is.na(x))
    mn <- mean(x[!is.na(x)])
    stdev <- mn * 0.1
    x[which(is.na(x))] <- rnorm(numbermissing, mean=mn, sd=stdev)
    return(x)
}
