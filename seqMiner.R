#!/usr/bin/env Rscript 

collapseListResult <- function(res) {
    matIdx <- which(unlist(lapply(res, is.matrix)))

    sampleId <- res[[length(res)]]
    res[[length(res)]] <- NULL
    mat <- res[matIdx]
    res[matIdx] <- NULL
    part1 <- do.call(cbind, res)

    if (length(matIdx) == 0) {
        return(part1)
    } else {
        matNames <- names(mat)
        n <- as.matrix(expand.grid(sampleId, matNames, stringsAsFactors=FALSE))
        mat.names <- apply(n, 1, function(x) {str_c(x, "", collapse = ":")})      
        mat.names
        part2 <- do.call(cbind, lapply(mat, function(a) {a <- t(a)}))
        colnames(part2) <- mat.names
        part2
        cbind(part1, part2)
    }
}

printMatrix <- function(res) {
    if (is.null(res)) {
        return(NULL);
    }
    cat (colnames(res), sep = "\t")
    cat ("\n")
    d <- apply(res, 1, function(x) { cat(x, sep = "\t"); cat("\n") })
}

#res <- readVCFToMatrixByRange(opt$args, opt$r, "")
printListOfMatrix <- function(res) {
    lres <- length(res)
    nres <- names(res)
    for (i in 1:lres) {
        cat(gettextf("----- %s -----\n", nres[i]))
        # save(res, file = "a.Rdata")
        m <- (cbind(rownames(res[[i]]), res[[i]]))
        colnames(m) <- c("Position", colnames(res[[i]]))
        printMatrix(m)
    }
}

#' Parse command line options
#'
#' @param args a character vector, usually put commandArgs(TRUE) here
#' @param shortArgs is for short options, e.g "ab:"
#' @param longArgs is for long options, e.g. c("opt1", "opt2=")
#' @return a list, where list names are parsed options names and list values are parsed results. Extra arguemnt are put in 'args' and 'nargs'
#' @examples
#' ret <- getopt(c("a", "-a", "-b",  "bval", "c", "--long1", "--long2", "long2Val"), "ab:c", c("long1", "long2=", "long3"))
getopt <- function(args, shortArgs, longArgs = NULL) {
    options(stringsAsFactors=FALSE)
    ## preprocess option strings
    opt <- list()
    ## parse short args
    if (!is.null(shortArgs)) {
        n <- nchar(shortArgs)
        i <- 1
        while (i <= n) {
            needArgs <- FALSE
            optChar <- substr(shortArgs, i, i)
            if (i < n && substr(shortArgs, i+1, i+1) == ":") {
                needArgs <- TRUE
                i <- i+1
            }
            optFull <- paste("-", optChar, sep = "")
            opt[[length(opt)+1]] <- list(optChar = optChar, optFull = optFull, needArgs = needArgs)
            i <- i + 1
        }
    }
    ## parse long args
    if (!is.null(longArgs)) {
        n <- length(longArgs)
        i <- 1
        while (i <= n) {
            optChar <- longArgs[i]
            last <- nchar(optChar)
            needArgs <- FALSE
            if ( substr(optChar, last, last) == "=") {
                optChar <- substr(optChar, 1, last - 1)
                needArgs <- TRUE
            }
            optFull <- paste("--", optChar, sep = "")
            opt[[length(opt) + 1]] <- list(optChar = optChar, optFull = optFull, needArgs = needArgs)
            i <- i + 1
        }
    }
    f = function(x) function(i) unlist(lapply(x, `[[`, i), use.names=FALSE)
    opt <- as.data.frame(Map(f(opt), names(opt[[1]])))
    #opt <- do.call(rbind, opt)
    #opt <- as.data.frame(opt)
    #print(opt)
    # sanity check
    if (any(grepl("^-", opt$optChar) | grepl("=$", opt$optChar))) {
        stop("Some options have incorrect '-' prefix or '=' suffix")
    }
    # print(opt)
    # go over each args
    res <- list(args = vector(mode = "character"), nargs = 0)
    i <- 1
    n <- length(args)
    while (i <= n) {
        loc <- which( args[i] == opt$optFull)
        if (length(loc) == 0) {
            res$nargs <- res$nargs + 1
            res$args <- c(res$args, args[i])
        } else if (length(loc) >= 1) {
            if (length(loc) > 1) {
                warning(gettextf("Option [ %s ] used more than once", opt$optFull))
                loc <- loc[length(loc)]
            }
            if (opt$needArgs[loc]) {
                if (i+1 > n) {
                    stop(gettextf("Don't have argument after option [ %s ]", args[i]))
                }
                res[[opt$optChar[loc]]] <- args[i+1]
                i <- i + 1
            } else {
                res[[opt$optChar[loc]]] <- "TRUE"
            }
        }
        i <- i + 1
    }
    class(res) <- "getopt"
    return(res)
}
has <- function(x, key) {UseMethod("has", x) }
has.getopt <- function(x, key) {
    return (!is.null(x[[key, exact = TRUE]]))
}
as.logical.getopt<- function(x, key) {
    return (as.logical(x[[key, exact = TRUE]]))
}
as.numeric.getopt <- function(x, key) {
    return (as.numeric(x[[key, exact = TRUE]]))
}
as.integer.getopt <- function(x, key) {
    return (as.integer(x[[key, exact = TRUE]]))
}
get <- function(x, key) {UseMethod("get", x) }
get.getopt <- function(x, key) {
    return (x[[key, exact = TRUE]])
}

args <- commandArgs(TRUE)
opt <- getopt(args, "r:he:a:n:", "geneFile=")
if (has(opt, "h") && as.logical(opt, "h")) {
    cat("seqMiner [[-r region] | [-g geneFile -n geneNames]] [-e vcfFields] [-a anno] inputFile\n")
    cat(" region: e.g. 1:100-200\n")
    cat(" vcfFields: VCF fields to extract from first 8 column, tags in the INFO field or individual fields. e.g. CHROM,POS:AF:GT,GD\n")
    cat(" annotation: e.g. Nonsynonymous, synonymous. You can use '|' to match more than one type of annotations\n")
    cat("\n")
    q("no")
}

suppressPackageStartupMessages(loadOK <- require(stringr, quietly = TRUE))
if (!loadOK) {
    cat("Please install seqminer using:\ninstall.packages(\"seqminer\")\n")
    q("no")
}
suppressPackageStartupMessages(loadOK <- require(seqminer, quietly = TRUE))
if (!loadOK) {
    cat("Please install seqminer using:\ninstall.packages(\"seqminer\")\n")
    q("no")
}

if (opt$nargs != 1) {
    stop(gettextf("Need to specify one VCF/BCF file, but %d given", opt$nargs))
}
if (has(opt, "e") > 0) { # use low level extract
    vcfField <- str_split(get(opt, "e"), ':')[[1]]
    vcfField <- lapply(vcfField, function(x) {str_split(x, ',')[[1]]})
    while (length(vcfField) < 3) {
        vcfField[[length(vcfField) + 1]] <- ""
    }
} else {
    vcfField <- NULL
}
if (opt$nargs > 1) {
    warning("more than one file name provided, only first will be processed")
}
if (has(opt, "a")) {
    anno = get(opt, "a")
} else {
    anno = ""
}

fileName <- opt$args[1]
if (has(opt, "r")) {
    if (!is.null(vcfField)) {
        res <- readVCFToListByRange(fileName, opt$r, anno, vcfField[[1]], vcfField[[2]], vcfField[[3]])
        printMatrix(collapseListResult(res))
        #save(res, file = "a.Rdata")
    } else {
        res <- readVCFToMatrixByRange(fileName, opt$r, anno)
        printListOfMatrix(res)
    }
    q("no")
} 
if (has(opt, "geneFile")) {
    if (!is.null(vcfField)) {
        res <- readVCFToListByGene(fileName, opt$geneFile, opt$n, anno, vcfField[[1]], vcfField[[2]], vcfField[[3]])
        printMatrix(collapseListResult(res))
    } else {
        res <- readVCFToMatrixByGene(fileName, opt$geneFile, opt$n, anno)
        # print(c(fileName, opt$geneFile, opt$n, anno))
        printListOfMatrix(res)
    }
    q("no")
}

