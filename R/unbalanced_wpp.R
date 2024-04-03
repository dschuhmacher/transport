#' @rdname unbalanced
#' @export
# default C minimal so that no cutoffs take place under the assumption that points in R^2
unbalanced.wpp <- function(a, b, p = 1, C = NULL, method = c("networkflow", "revsimplex"),
                           output = c("dist", "all", "rawres"), threads=1, ...) {
  stopifnot(is(a, "wpp") && is(b, "wpp"))
  stopifnot(a$dimension == b$dimension)
  if (a$dimension < 2) stop("dimension must be >=2")
  
  method <- match.arg(method)
  output <- match.arg(output)  
  
  if (is.null(C)) {
    coords <- rbind(a$coordinates, b$coordinates)
    span <- apply(coords, 2, \(x) {diff(range(x))})
    C <- (sqrt(sum(span^2)) / 2)^(1/p)
    # so that tunneling of points (deleting units of mass and adding them again somewhere else) is
    # never worth it (except in extreme cases due to rounding issues).
  }
  
  m <- a$N
  n <- b$N
  amass <- a$mass
  atotmass <- a$totmass
  bmass <- b$mass
  btotmass <- b$totmass
 
  if (m == 0 && n == 0) return(0)  # this line is not needed (if a is consistent) but makes things clearer
 
  if (m == 0) return(btotmass)  # atotmass = 0
    
  if (n == 0) return(atotmass)  # btotmass = 0
    
  costm <- gen_cost(a$coordinates, b$coordinates, threads=threads)^(p/2) 
  # could take pmin with 2*C^p, but not needed
  costm <- rbind(costm, C^p)  # we do *not* divide by two! (MSM not HKM)
  costm <- cbind(costm, C^p)
  costm[m+1, n+1] <- 0
  amassplus <- c(amass, btotmass)  
  bmassplus <- c(bmass, atotmass)

  if (method == "networkflow") {
    rawres <- networkflow(matrix(amassplus), matrix(bmassplus), costm, threads=threads) # m x 1 and n x 1 matrices for the masses 
    # sanity check (not sure if already performed inside networkflow):
    primalcost <- sum(costm * rawres$plan)   
    dualcost <- sum(rawres$potential * c(amassplus, bmassplus))
    if (!isTRUE(all.equal(primalcost, dualcost)) || !isTRUE(all.equal(primalcost, rawres$dist))) {  # dist is cost (dist^p)
      warning("Primal-dual gap is ", rawres$dist - dualcost, "\n", "Primal cost: ", primalcost, 
              "; dual cost: ", dualcost, "; rawres$cost: ", rawres$dist)
    }
  } else {
    rawres <- unbalanced_revsimplex_core(amassplus, bmassplus, costm, p, C) 
    # unbalanced_revsimplex_core --> unbalanced.R
    # this rawres does not have a component potential (to be fixed)
  }
  
  if (output == "dist") {
    return(rawres$dist^(1/p)) # rawres$dist is the p-th power of the unbalanced Wasserstein dist
  }
  
  # emulates the output of transport with networkflow and fullreturn=TRUE (trashcan states added)
  if (output == "rawres") { 
    rawres$frame <- rawres$frame[rawres$frame[,3]>0,,drop=FALSE]
    df <- data.frame(from=rawres$frame[,1], to=rawres$frame[,2], mass=rawres$frame[,3])
    out <- list(default=df, primal=rawres$plan, dual=rawres$potential, cost=rawres$dist)
    return(out)
  }
  
  # output = "all"
  rawres$aextra <- rawres$plan[1:m,n+1,drop=FALSE]
  rawres$bextra <- rawres$plan[m+1,1:n,drop=FALSE]
  select <- (rawres$frame[,3] > 0) & # we keep transports over dist 0 in current version (rawres$frame[,1] == rawres$frame[,2]) 
            (rawres$frame[,1] <= m) & (rawres$frame[,2] <= n)  # removes the trashcan states from frame
  rawres$frame <- rawres$frame[select, , drop=FALSE]
  
  result <- list(dist=rawres$dist^(1/p), plan=rawres$frame)
  result$aextra <- as.numeric(rawres$aextra)
  result$bextra <- as.numeric(rawres$bextra)
 
  # there is no component inplace here (for now at least). If a location of a matches a location of b
  # we simply have transport of distance zero, but we do not detect/mark this in any way.
  # The reason why atrans,btrans is computed after aextra,bextra (and then switched) is only to keep
  # the same order as in unbalanded.pgrid
  
  # fill result$atrans: 
  atemp <- rowsum(result$plan[,3], result$plan[,1])
  where <- as.numeric(attr(atemp, "dimnames")[[1]])
  result$atrans <- rep(0, m)
  result$atrans[where] <- atemp
  # fill result$btrans: 
  btemp <- rowsum(result$plan[,3], result$plan[,2])
  where <- as.numeric(attr(btemp, "dimnames")[[1]])
  result$btrans <- rep(0, n)
  result$btrans[where] <- btemp
  
  tol <- ifelse(method == "networkflow", sqrt(.Machine$double.eps), 1e-7)  
  # for networkflow standard tolerance, for revsimplex somewhat smaller due to
  # smaller precision (too many spurious warnings otherwise)
  if (!isTRUE(all.equal(result$atrans + result$aextra, a$mass, tolerance=tol, check.attributes = FALSE))) {
    warning("atrans and aextra do not sum up to a$mass. ", 
            all.equal(result$atrans + result$aextra, a$mass))
    # gives mean relative difference |left-right|/left, no sign!!
  }
  if (!isTRUE(all.equal(result$btrans + result$bextra, b$mass, tolerance=tol, check.attributes = FALSE))) {
    warning("btrans and bextra do not sum up to b$mass. ", 
            all.equal(result$btrans + result$bextra, b$mass))
    # gives mean relative difference |left-right|/left, no sign!!
  }
  
  result$plan <- data.frame(from=result$plan[,1], to=result$plan[,2], mass=result$plan[,3])
  result <- result[c("dist", "plan", "atrans", "btrans", "aextra", "bextra")]
  attr(result, "a") <- a
  attr(result, "b") <- b
  attr(result, "p") <- p
  attr(result, "C") <- C
  class(result) <- "ut_wpp"
  return(result)
}



#' Plot Unbalanced Transport Information
#'
#' Graphic representation of components of the list returned by \code{\link{unbalanced}}.
#'
#' @param x the list returned by \code{\link{unbalanced}} with option \code{output="all"}.
#' @param what character. The aspect of the unbalanced transport information to display.
#' @param axes logical. Whether to plot axes (ignored for \code{what="plan"}).
#' @param xlim,ylim numeric vectors of length 2. The x- and y-limits of the plot.
#' @param ... further graphics parameters passed to \code{\link{plot.pgrid}} for 
#' \code{what="plan"} and passed to \code{\link{matimage}} in all other cases.
#'
#' @return Nothing. Used for the side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(33)
#' m <- 50
#' n <- 20
#' massa <- rexp(m)
#' massb <- rexp(n)
#' a <- wpp( matrix(runif(2*m), m, 2), massa)
#' b <- wpp( matrix(runif(2*n), n, 2), massb)
#' res <- unbalanced(a,b,1,0.3,output="all")
#' plot(res, what="plan")
#' plot(res, what="trans")
#' plot(res, what="extra")}
# possible improvements: also draw the information that is not in the focus of the plot (as specified by what); 
# also the cex is independent for every choice of what (so the areas of the dots don't add up between different "what"s) 
# in the long run we might also directly depicts all the info in one plot (without focus?) and
# and we should definitely switch to ggplot2
plot.ut_wpp <- function(x, what=c("plan", "extra", "trans"), axes=FALSE, xlim=c(0,1), ylim=c(0,1), ...) {
  stopifnot(is(x, "ut_wpp"))
  what <- match.arg(what)
  
  a <- attr(x, "a")
  b <- attr(x, "b")
  
  if (what == "plan") {
    plot.wpp(a, b, x$plan, xlim=xlim, ylim=ylim, ...)
  }
  
  if (what == "extra") {
    # since wpp objects must not have total mass = 0 (this should probably be fixed in the long run)
    # we must proceed in a more complicated way
    if (sum(x$aextra) > 0) {
      aextra <- wpp(a$coordinates, x$aextra)  # wpp removes zero-mass points anyway (not visible as "toosmall")
    } else {
      aextra <- wpp(matrix(c(xlim[1],ylim[1]),1,2)-1,1)  # wpp not visible in plot
    }
    if (sum(x$bextra) > 0) {
      bextra <- wpp(b$coordinates, x$bextra)  # wpp removes zero-mass points anyway (not visible as "toosmall")
    } else {
      bextra <- wpp(matrix(c(xlim[1],ylim[1]),1,2)-1,1)  # wpp not visible in plot
    }
    
    plot.wpp(aextra, bextra, xlim=xlim, ylim=ylim, ...)
  }
  
  if (what == "trans") {
    transsum <- sum(x$atrans)
    if (!isTRUE(all.equal(transsum, sum(x$btrans), tol=1e-5))) {  # for plotting we can be tolerant 
      warning("sums over atrans and btrans are not the same. ", 
              all.equal(x$atrans + x$aextra, a$mass))
      # gives mean relative difference |left-right|/left, no sign!!
    }
    # since wpp objects must not have total mass = 0 (this should probably be fixed in the long run)
    if (transsum == 0) {
      plot.wpp(wpp(matrix(0,1,2),1), wpp(matrix(0,1,2),1), xlim=xlim, ylim=ylim, cex=0, ...)
      warning("no mass transported. Plot is empty.")
    } else {
      wha <- which(x$atrans > 0)
      whb <- which(x$btrans > 0)
      atrans <- wpp(a$coordinates[wha,], x$atrans[wha])  # because wpp removes zero-mass points anyway
      btrans <- wpp(b$coordinates[whb,], x$btrans[whb])  # because wpp removes zero-mass points anyway
      plan <- x$plan
      plan$from  <- match(plan$from, wha)
      plan$to <- match(plan$to, whb)
      plot.wpp(atrans, btrans, plan, xlim=xlim, ylim=ylim, ...) 
    }
  }
  
  invisible()
}

# TO DO: print method for ut_wpp

