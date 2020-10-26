readData <- function(file.xlsx) {
    require(readxl)
    require(stringr)

    nms <- names(read_excel(file.xlsx, n_max=0))
    datix <- which(grepl('Annotated Sequence|^Protein Accessions|File ID|Isolation Interference|Abundance', nms))
    col_types_data <- ifelse(seq(length(nms)) %in% datix, 'guess', 'skip')
    
    dat.df <- as.data.frame(read_excel(file.xlsx, col_types=col_types_data))
    colnames(dat.df) <- str_replace_all(colnames(dat.df), ' ', '')
    rownames(dat.df) <- make.unique(dat.df$AnnotatedSequence, sep=';')
    
    dat.df <- dat.df[-grep('AnnotatedSequence', colnames(dat.df))]
    dat.df$ProteinAccessions <- as.character(sapply(dat.df$'ProteinAccessions', function(s) unlist(strsplit(s, split='; '))[1]))

    return(dat.df)
}

readAnnotation <- function(file.xlsx) {
    require(readxl)
    require(stringr)

    nms <- names(read_excel(file.xlsx, n_max=0))
    annotix <- which(grepl('Annotated Sequence|^Protein Accessions', nms))
    col_types_annot <- ifelse(seq(length(nms)) %in% annotix, 'guess', 'skip')
    
    annot.df <- as.data.frame(read_excel(file.xlsx, col_types=col_types_annot))
    colnames(annot.df) <- str_replace_all(colnames(annot.df), ' ', '')
    annot.df$ProteinAccessions <- as.character(sapply(annot.df$'ProteinAccessions', function(s) unlist(strsplit(s, split='; '))[1]))
    annot.df$uniqAnnotatedSequence <- make.unique(annot.df$AnnotatedSequence, sep=';')

    return(annot.df)
}

filterIsolationInterference <- function(df, percent=70) {
    percent <- as.numeric(percent)
    ix <- which(grepl('IsolationInterference', colnames(df)))
    isolationinterference <- df[, ix]
    df <- df[isolationinterference < percent & !is.na(isolationinterference), ]
    df[-ix]
}

filterChannels <- function(df, carrierORempty, ...) {
    xlst <- list(...)
    carrierORempty <- c(carrierORempty, unlist(xlst))
    cix <- grep(paste(carrierORempty, collapse='|'), colnames(df))
    df <- df[, -cix]
}

filterRows <- function(df) {
    emptyrows <- apply(df, 1, function(x) all(is.na(x)))
    df[!emptyrows, ]
}

partitionMSruns <- function(df) {
    runix <- which(grepl('FileID', colnames(df)))
    partition.lst <- split(df, factor(df[, runix, drop=TRUE]))

    partition.lst <- lapply(partition.lst, function(xdf) {
        prefix <- unique(xdf[, runix, drop=TRUE])
        abundix <- which(grepl('Abundance:', colnames(xdf)))
        colnames(xdf)[abundix] <- gsub('Abundance:', paste0(prefix, '_'), colnames(xdf)[abundix])
        xdf <- xdf[-runix]
    })
}

partitionProteins <- function(df) {
    split(df, factor(df$ProteinAccessions))
}

extrapolateMissing <- function(df) {
    # df is ctrl or treatment partial data.frame
    if (!anyNA(df)) {
        return(df)
    } else if (nrow(df)  == 1) {
        fillRowColumn(df)
    }
}

fillRowColumn <- function(df) {
    res <- apply(df, 1, function(v) {
        if (sum(!is.na(v)) >= 2) {
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








