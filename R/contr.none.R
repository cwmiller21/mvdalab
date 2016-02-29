contr.none <- function (n, contrasts) 
{
  if (length(n) == 1) 
    contr.treatment(n, contrasts = n <= 2)
  else contr.treatment(n, contrasts = length(unique(n)) <= 2)
}