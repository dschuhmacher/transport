wasserstein1d <- function(a, b, p=1, wa=NULL, wb=NULL) {
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  if (m == n && is.null(wa) && is.null(wb)) {
    return(mean(abs(sort(b)-sort(a))^p)^(1/p))
  }
  stopifnot(is.null(wa) || length(wa) == m)
  stopifnot(is.null(wb) || length(wb) == n)
  if (is.null(wa)) {
    wa <- rep(1,m)
  } else { # remove points with zero weight
    wha <- wa > 0
    wa <- wa[wha]
    a <- a[wha]
    m <- length(a)
  }
  if (is.null(wb)) {
    wb <- rep(1,n)
  } else { # remove points with zero weight
    whb <- wb > 0
    wb <- wb[whb]
    b <- b[whb]
    n <- length(b)
  }

  orda <- order(a)
  ordb <- order(b)
  a <- a[orda]
  b <- b[ordb]
  wa <- wa[orda]
  wb <- wb[ordb]
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))  
  cub <- c(cumsum(ub))  
  arep <- hist(cub, breaks = c(-Inf, cua, Inf), plot = FALSE)$counts + 1
  brep <- hist(cua, breaks = c(-Inf, cub, Inf), plot = FALSE)$counts + 1
  # We sum over rectangles with cuts on the vertical axis each time one of the two ecdfs makes a jump.
  # arep and brep tell us how many times each of the a and b data have to be repeated in order to get
  # the points on the horizontal axis.
  # Note that sum(arep)+sum(brep) = m+n-1 (we do not count the height-zero final rectangle where both ecdfs jump to 1)

  aa <- rep(a, times=arep)
  bb <- rep(b, times=brep)

  uu <- sort(c(cua,cub))
  uu0 <- c(0,uu)
  uu1 <- c(uu,1)
  areap <- sum((uu1-uu0)*abs(bb-aa)^p)^(1/p)
  return(areap)
}


# plots without weights
# 
# plot(ecdf(a),verticals=TRUE, xlim=c(-2,2))
# lines(ecdf(b),col=2,verticals=TRUE)


# plots with weights
#
# ua <- (wa/sum(wa))[-length(wa)]
# ub <- (wb/sum(wb))[-length(wb)]
# cua <- c(cumsum(ua))  
# cub <- c(cumsum(ub))  
# plot(sort(c(-2,a,2)),c(0,cua,1,1),type="s",xlim=c(-2,2))
# lines(sort(c(-2,b,2)),c(0,cub,1,1),type="s",col=2)
