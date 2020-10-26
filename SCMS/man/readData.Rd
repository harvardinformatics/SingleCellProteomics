\name{readData}
\alias{readData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
readData(file.xlsx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.xlsx}{
%%     ~~Describe \code{file.xlsx} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (file.xlsx) 
{
    require(readxl)
    require(stringr)
    nms <- names(read_excel(file.xlsx, n_max = 0))
    datix <- which(grepl("Annotated Sequence|^Protein Accessions|File ID|Isolation Interference|Abundance", 
        nms))
    col_types_data <- ifelse(seq(length(nms)) \%in\% datix, "guess", 
        "skip")
    dat.df <- as.data.frame(read_excel(file.xlsx, col_types = col_types_data))
    colnames(dat.df) <- str_replace_all(colnames(dat.df), " ", 
        "")
    rownames(dat.df) <- make.unique(dat.df$AnnotatedSequence, 
        sep = ";")
    dat.df <- dat.df[-grep("AnnotatedSequence", colnames(dat.df))]
    dat.df$ProteinAccessions <- as.character(sapply(dat.df$ProteinAccessions, 
        function(s) unlist(strsplit(s, split = "; "))[1]))
    return(dat.df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
