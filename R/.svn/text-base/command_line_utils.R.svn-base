#' given a command line object (a list with named items for each commandline arg)
#' from optparse (way worse than argparse), parse out boolean values 
#' @param an optparse object, with variables as named strings
#' @return a post processed version of the object, with boolean represented as such,
#' 		integers converted, etc
#' @keywords command line parse
post.process.commandline <- function(optparse.args) {
  ret <- list()
  for (nm in names(optparse.args)) {
    # if this returns true then we know the variable is set
    if (toupper(optparse.args[[nm]]) == "TRUE") {ret[[nm]] = TRUE}
    else if (toupper(optparse.args[[nm]]) == "FALSE") {ret[[nm]] = FALSE}
    else {ret[[nm]] = optparse.args[[nm]]}  
  }
  return(ret)
}
