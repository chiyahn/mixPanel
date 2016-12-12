#' @description Returns a list with original column and columns that are lagged up to s.
#' @title GetLaggedAndSample
#' @name GetLaggedAndSample
#' @param y n by 1 column
#' @param s The maximum order of lag
#' @return A list that consists of the followings:
#' \item{y.lagged}{(n-s) by 1 column where the first s data are cut}
#' \item{y.sample}{(n-s) by s matrix whose ith column represents ith lagged column.}
GetLaggedAndSample <- function(y, s)
{
  y <- as.numeric(y)
  y.lagged <- sapply(seq(0,s), GetLaggedColumn, y, s) # (n-s) by s matrix
  y.sample <- as.matrix(y.lagged[,1])
  if (s > 0)
    y.lagged <- as.matrix(y.lagged[,-1])
  else
    y.lagged <- matrix(0,ncol=1,nrow=length(y))
  return (list (y.lagged = y.lagged, y.sample = y.sample))
}

#' @description Get a column of 0 \leq j \leq s lagged variable
#' @title GetLaggedColumn
#' @name GetLaggedColumn
#' @param j Lag order of interest
#' @param column Column to be lagged
#' @param s Maximum lag
#' @return A column that is lagged by j.
GetLaggedColumn <- function (j, col, s) {
  if (j != s)
    col <- col[-(0:(s-j))] # destroy first s-j elements
  return (col[1:(length(col)-j)])
}
