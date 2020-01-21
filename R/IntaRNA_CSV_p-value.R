#!/usr/bin/env Rscript

####################################################################
# Computes p-values and false discovery rates (fdr following Benjamin+Hochberg)
# by fitting a GEV on the energy values computed by IntaRNA. 
# Note, such p-value estimates are only useful for genome-wide predictions.
#
# arguments: <IntaRNA-output-CSV> [<out-CSV> = <intarna-csv-output>] [<col-name-E> = E]
#
# 1 <IntaRNA-output-CSV> = ";"-separated CSV output of IntaRNA
# 2 <out-CSV> = file name to write the extended CSV output to (2 new columns)
# 3 <col-name-E> = the column name that holds the energy values to be fitted
#
# example call:
#
# Rscript --vanilla IntaRNA_CSV_p-value.R predictions.csv
#
# This script is part of the IntaRNA source code package. See
# respective licence and documentation for further information.
#
# https://github.com/BackofenLab/IntaRNA
#
####################################################################


####################################################################
# get command line arguments
####################################################################

args = commandArgs(trailingOnly=TRUE)
# check and parse
if (length(args)<1) { stop("call with <intarna-csv-output> [<out-file-with-p-values> = <intarna-csv-output>] [<col-name-E> = E]", call.=FALSE) }

# get input file = IntaRNA csv output
intarnaOutputFile = args[1];
if (!file.exists(intarnaOutputFile )) { stop("intarna-csv-output file '", intarnaOutputFile, "' does not exist!", call.=FALSE) }

# get output file
outFile = intarnaOutputFile;
if (length(args)>=2) {
	outFile = args[2]
}

# set column to get energies from
colNameE = "E"
# get column name from argument if present
if (length(args)>=3) {
	colNameE = args[3]
}

# column delimiter used in CSV input / output
csvColSep = ";"

# number of digits of p-values
pValPrec = 7

####################################################################
# fits a generalized extreme value distribution to the given energy data
# adopted from 'gev' function of 'evir' library
# @param energy the IntaRNA energy values to fit (a vector)
# @return the fitting parameters xi, mu, and sigma
gevFitting <- function (energy)
####################################################################
{
    n.all <- NA
    energy <- as.numeric(-energy)
    n <- length(energy)
    sigma0 <- sqrt(6 * var(energy))/pi
    mu0 <- mean(energy) - 0.57722 * sigma0
    xi0 <- 0.1
    theta <- c(xi0, sigma0, mu0)
    negloglik <- function(theta, tmp) {
        y <- 1 + (theta[1] * (tmp - theta[3]))/theta[2]
        if ((theta[2] < 0) || (min(y) < 0))
            out <- 1e+06
        else {
            term1 <- length(tmp) * logb(theta[2])
            term2 <- sum((1 + 1/theta[1]) * logb(y))
            term3 <- sum(y^(-1/theta[1]))
            out <- term1 + term2 + term3
        }
        out
    }
	# compute fit
    fit <- optim(theta, negloglik, hessian = TRUE, tmp = energy)
    if (fit$convergence)
        warning("gev fit optimization may not have succeeded")
	
	return( list( xi=fit$par[1], sigma=fit$par[2], mu=fit$par[3] ) )
}


####################################################################
# computes p-values for the given energy values and GEV distribution
# adopted from 'pgev' function of 'evir' library
# @param energy IntaRNA energy values
# @param gev GEV parameters
# @return p-values for each energy value
gevPvalue <- function (energy, gev=list( xi = 1, mu = 0, sigma = 1) )
####################################################################
{
	return ( 1 - exp( - (1  + (gev$xi * ((-energy) - gev$mu))/gev$sigma)^(-1 /gev$xi)))
}



####################################################################
# parse IntaRNA CSV
####################################################################

d = read.csv( intarnaOutputFile, sep=csvColSep )

# check if energy column present
if (!is.element(colNameE, colnames(d))) { 
	stop("'",colNameE,"' is not among the column names of '",intarnaOutputFile,"'", call.=FALSE); 
}
# check if unique
if (sum(colnames(d) == colNameE)>1) {
	stop("column name '",colNameE,"' occurs more than once in '",intarnaOutputFile,"'", call.=FALSE); 
}


####################################################################
# fit p-values
####################################################################

# get energies to fit
E = d[,colnames(d) == colNameE]


# fit negated energies 
gevfit <- gevFitting(E) # fitten

# get rounded pValue
pVal <- round( gevPvalue( E, gevfit ), digits=pValPrec )
# get rounded fdr
fdr <- round( p.adjust(pVal, method="BH"), digits=pValPrec )
		
####################################################################
# write output
####################################################################

o = cbind( d, pVal, fdr )
colnames(o)[ncol(o)-1] = "p-value"
colnames(o)[ncol(o)] = "fdr"

write.table( o, outFile, sep=csvColSep, row.names=FALSE, col.names = TRUE, quote=FALSE )


#################################################################EOF
