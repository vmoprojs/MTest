####################################################
### File name: pairwiseKStets.r
####################################################
pairwiseKStest <- function(X,alternative="greater")
{
  #Returns the p value of the pairwise KS test of X columns
  n <- ncol(X)
  sol <- matrix(NA, ncol = n, nrow = n)
  for(i in 1:(n))
  {
    for(j in (1):n)
    {
      a <- suppressWarnings(ks.test(X[,i],X[,j],alternative = alternative))
      sol[i,j] <- a$p.value
    }
  }
  if(alternative=="less") {mes = ("alternative hypothesis: the CDF of x lies below that of y. Rows are `x` and Columns are `y`")}
  if(alternative=="greater") {mes = ("alternative hypothesis: the CDF of x lies above that of y. Rows are `x` and Columns are `y`")}
  if(alternative=="two.sided") {mes = ("alternative hypothesis: two-sided")}
  colnames(sol) <- colnames(X)
  rownames(sol) <- colnames(X)
  return(list(KSpwMatrix = sol,mes=mes))
}