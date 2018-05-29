.pnfwunorm=function(q, con=5){
  return(log(1 + q*con)-(con*q)/(1 + con*q))
}

dnfw = function(x, con=5, log=FALSE){
  if(log){
    d = log(con^2*x/((con*x+1)^2*(1/(con+1)+log(con+1)-1)))
  }else{
    d = con^2*x/((con*x+1)^2*(1/(con+1)+log(con+1)-1))
  }
  d[x>1]=0
  d[x<=0]=0
  return(d)
}

pnfw = function(q, con=5, log.p=FALSE){
  p = .pnfwunorm(q=q, con=con)/.pnfwunorm(q=1, con=con)
  p[q>1]=1
  p[q<=0]=0
  if(log.p){
    return(log(p))
  }else{
    return(p)
  }
}

qnfw = function(p, con=5, log.p=FALSE){
  if(log.p){
    p=exp(p)
  }
  p[p>1]=1
  p[p<=0]=0
  p=p*.pnfwunorm(q=1, con=con)
  if(requireNamespace("lamW", quietly = TRUE)){
    return((-(1/lamW::lambertW0(-exp(-p-1)))-1)/con)
  }else if(requireNamespace("gsl", quietly = TRUE)){
    return((-(1/gsl::lambert_W0(-exp(-p-1)))-1)/con)
  }else{
    message('Using slow internal Lambert W0, ideally you shouls install the lamW or the gsl package.')
    return((-(1/ .lambertWp(-exp(-p-1)))-1)/con)
  }
}

rnfw = function(n, con=5){
  if(length(n)>1){
    n=length(n)
  }
  return(qnfw(p=runif(n), con=con))
}

.lambertWp=function(x){
  #taken from the pracma package under the GPL 3+ license
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric (real) vector.")
    if (length(x) == 1) {
        if (x < -1/exp(1))
            return(NaN)
        if (x == -1/exp(1))
            return(-1)
        if (x <= 1) {
            eta <- 2 + 2 * exp(1) * x
            f2 <- 3 * sqrt(2) + 6 - (((2237 + 1457 * sqrt(2)) *
                exp(1) - 4108 * sqrt(2) - 5764) * sqrt(eta))/((215 +
                199 * sqrt(2)) * exp(1) - 430 * sqrt(2) - 796)
            f1 <- (1 - 1/sqrt(2)) * (f2 + sqrt(2))
            w0 <- -1 + sqrt(eta)/(1 + f1 * sqrt(eta)/(f2 + sqrt(eta)))
        }
        else {
            w0 = log(6 * x/(5 * log(12/5 * (x/log(1 + 12 * x/5)))))
        }
        w1 <- w0 - (w0 * exp(w0) - x)/((w0 + 1) * exp(w0) - (w0 +
            2) * (w0 * exp(w0) - x)/(2 * w0 + 2))
        while (abs(w1 - w0) > 1e-15) {
            w0 <- w1
            w1 <- w0 - (w0 * exp(w0) - x)/((w0 + 1) * exp(w0) -
                (w0 + 2) * (w0 * exp(w0) - x)/(2 * w0 + 2))
        }
        return(w1)
    }
    else {
      sapply(x, .lambertWp)
    }
}

