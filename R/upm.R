UPM <- function(w, a = 0, b = 1, a1 = 1, b1 = 1, a2 = 1, b2 = 1)
{
  if(w < 0 || w >1)
  {
    stop("w iw weight taken on first prior (informative), which can lie between 0 and 1")
  }
  if(a < 0)
  {
    a = 0
    warning("Domain of Beta distribution is (0,1)")
  }
  if(b > 1)
  {
    b = 1
    warning("Domain of Beta distribution is (0,1)")
  }
  if(a >= b) stop("a must be less than b")

  val =  w * ((pbeta(b, shape1 = a1, shape2 = b1) - pbeta(a, shape1 = a1, shape2 = b1)) / (b - a)) +
    (1 - w) * ((pbeta(b, shape1 = a2, shape2 = b2) - pbeta(a, shape1 = a2, shape2 = b2)) / (b - a))
  return(val)
}
