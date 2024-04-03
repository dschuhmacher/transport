#' Unbalanced Optimal Transport Between Two Objects
#' 
#' Compute optimal transport between unnormalized images / mass distributions on grids 
#' (\code{pgrid} objects) or between mass distributions on general point patterns
#' (\code{wpp} objects) under the option that mass can be dispose of. Transport cost
#' per unit is the Euclidean distance of the transport to the \code{p}-th power.
#' Disposal cost per unit is \code{C^p}.\cr 
#'
#' @param a,b objects of class \code{\link{pgrid}} or \code{\link{wpp}} that are compatible.
#' @param p a power \eqn{\geq 1} applied to the transport and disposal costs. The order
#' of the resulting unbalanced Wasserstein metric.
#' @param C The base disposal cost (without the power \code{p}) 
#' @param method one of \code{"networkflow"} and \code{"revsimplex"}, specifing the algorithm used. See details.
#' @param output character. One of "dist", "all" and "rawres". Determines what the function
#' returns: only the unbalanced Wasserstein distance; all available information about the 
#' transport plan and the extra mass; or the raw result obtained by the networkflow algorithm.
#' The latter is the same format as in the \code{transport} function with option \code{fullreturn=TRUE}.
#' The choice \code{output = "rawres"} is mainly intended for internal use.
#' @param threads an integer specifying the number of threads for parallel computing in connection with
#' the networkflow method.
#' @param ... other arguments.
#' 
#' @details Given two non-negative mass distributions \eqn{a=(a_x)_{x \in S}}, \eqn{b=(a_y)_{y \in S}}
#' on a set \eqn{S} (a pixel grid / image if \code{a}, \code{b} are of class \code{pgrid} or a more
#' general weighted point pattern if \code{a}, \code{b} are of class \code{wpp}), this function minimizes the functional
#' \deqn{\sum_{x,y \in S} \pi_{x,y} d(x,y)^p + C^p \bigl( \sum_{x \in S} (a_x - \pi^{(1)}_x) + \sum_{y \in S} (b_y - \pi^{(2)}_y) \bigr)}
#' over all \eqn{(\pi_{x,y})_{x,y \in S}} satisfying
#' \deqn{0 \leq \pi^{(1)}_x := \sum_{y \in S} \pi_{x,y} \leq a_x \ \textrm{and} \ 0 \leq \pi^{(2)}_y := \sum_{x \in S} \pi_{x,y} \leq b_y.}
#' 
#' Thus \eqn{\pi_{x,y}} denotes the amount of mass transported from \eqn{x} to \eqn{y}, whereas \eqn{\pi^{(1)}_x}
#' and \eqn{\pi^{(2)}_y} are the total mass transported away from \eqn{x} and total mass transported to \eqn{y}, respectively.
#' Accordingly \eqn{\sum_{x \in S} (a_x - \pi^{(1)}_x)} and \eqn{\sum_{y \in S} (b_y - \pi^{(2)}_y)} are the total
#' amounts of mass of \eqn{a} and \eqn{b}, respectively, that need to be disposed of.
#' 
#' The minimal value of the functional above taken to the \eqn{1/p} is what we refer to as unbalanced
#' \eqn{(p,C)}-Wasserstein metric. This metric is used, in various variants, in an number of research papers.
#' See Heinemann et al. (2022) and the references therein and Müller et al. (2022), Remark 3. We follow the
#' convention of the latter paper regarding the parametrization and the use of the term \emph{unbalanced Wasserstein metric}.
#' 
#' The practical difference between the two methods "networkflow" and "revsimplex" can 
#' roughly described as follows. The former is typically faster for large examples (for \code{pgrid} objects 64x64
#' and beyond), especially if several threads are used. The latter is typically faster
#' for smaller examples (which may be relevant if pairwise transports between many objects
#' are computed) and it guarantees a sparse(r) solution, i.e. at most \eqn{m+n+1} individual
#' transports, where \eqn{m} and \eqn{n} are the numbers of non-zero masses in \code{a} and \code{b}, respectively).
#' Note however that due to the implementation the revsimplex algorithm is a little less
#' precise (roughly within 1e-7 tolerance). For more details on the algorithms see \code{\link{transport}}.
#'
#' @return If \code{output = "dist"} a single numeric, the unbalanced \eqn{(p,C)}-Wasserstein distance.
#' Otherwise a list. If \code{output = "all"} the list is of class \code{ut_pgrid} or \code{ut_wpp} according
#' to the class of the objects \code{a} and \code{b}. It has \code{a}, \code{b}, \code{p}, \code{C} as attributes and 
#' the following components:
#' \item{dist}{same as for \code{output = "dist"}.}
#' \item{plan}{an optimal transport plan. This is a data frame with columns \code{from}, \code{to} and \code{mass}
#'  that specifies from which element of \code{a} to which element of \code{b} what amount of mass is sent.
#'  \code{from} and \code{to} are specified as vector indices in terms of the usual column major enumeration
#'  of the matrices a$mass and b$mass. The plan can be plotted via \code{plot.pgrid(a, b, plan)}.}
#' \item{atrans, btrans}{matrices (pgrid) or vectors (wpp) specifying the masses transported from each point and to each point,
#'  respectively. Corresponds to \eqn{(\pi^{(1)}_x)_{x \in S}} and \eqn{(\pi^{(2)}_y)_{y \in S}} above.}
#' \item{aextra, bextra}{matrices (pgrid) or vectors (wpp) specifying the amount of mass at each point of \code{a} and \code{b},
#' respectively, that cannot be transported and needs to be disposed of. Corresponds to
#' \eqn{(a_x - \pi^{(1)}_x)_{x \in S}} and \eqn{(b_y - \pi^{(2)}_y)_{y \in S}}.}
#' \item{inplace}{(pgrid only) a matrix specifying the amount of mass at each point that can stay in place. Corresponds
#' to \eqn{(\pi_{x,x})_{x \in S}}.}
#'  
#' Note that \code{atrans + aextra + inplace} (pgrid) or \code{atrans + aextra} (wpp)must be equal
#' to \code{a$mass} and likewise for b.
#' A warning occurs if this is not the case (which may indeed happen from time to time for method
#' revsimplex, but the error reported should be very small).
#' 
#' @references Florian Heinemann, Marcel Klatt and Axel Munk (2022).\cr
#'             Kantorovich-Rubinstein distance and barycenter for finitely supported measures: Foundations and Algorithms.\cr
#'             Arxiv preprint.\cr
#'             \doi{10.48550/arXiv.2112.03581}\cr
#'             \cr
#'             Raoul Müller, Dominic Schuhmacher and Jorge Mateu (2020).\cr
#'             Metrics and barycenters for point pattern data
#'             Statistics and Computing 30, 953-972.\cr
#'             \doi{10.1007/s11222-020-09932-y}
#'             
#' @seealso \code{\link{plot.ut_pgrid}} and \code{\link{plot.ut_wpp}}, which can plot the various components of the list obtained for \code{output="all"}.
#'
#' @export
#'
#' @examples
#' a <- pgrid(matrix(1:12, 3, 4))
#' b <- pgrid(matrix(c(9:4, 12:7), 3, 4))
#' res1 <- unbalanced(a, b, 1, 0.5, output="all")
#' res2 <- unbalanced(a, b, 1, 0.3, output="all")
#' plot(a, b, res1$plan, angle=20, rot=TRUE)
#' plot(a, b, res2$plan, angle=20, rot=TRUE)
#' par(mfrow=c(1,2))
#' matimage(res2$aextra, x = a$generator[[1]], y = a$generator[[2]])
#' matimage(res2$bextra, x = b$generator[[1]], y = b$generator[[2]])
#' 
#' set.seed(31)
#' a <- wpp(matrix(runif(8),4,2), 3:6)
#' b <- wpp(matrix(runif(10),5,2), 1:5)
#' res1 <- unbalanced(a, b, 1, 0.5, output="all")
#' res2 <- unbalanced(a, b, 1, 0.3, output="all")
#' plot(a, b, res1$plan)
#' plot(a, b, res2$plan)
#'
#' @export
unbalanced <- function(a, b, ...) {
  stopifnot(class(a) == class(b))	
  UseMethod("unbalanced")
}


#' @rdname unbalanced
#' @export
unbalanced.pgrid <- function(a, b, p = 1, C = NULL, method = c("networkflow", "revsimplex"),
                       output = c("dist", "all", "rawres"), threads=1, ...) {
  stopifnot(is(a, "pgrid") && is(b, "pgrid"))
  stopifnot(compatible(a,b))
  if (a$dimension < 2) stop("pixel grids of dimension >= 2 required")
  if (!(a$structure %in% c("square", "rectangular")))
    stop("transport.pgrid is currently only implemented for rectangular pixel grids")
  
  method <- match.arg(method)
  output <- match.arg(output)
  
  if (is.null(C)) {
    d <- length(a$boundary)
    span <- a$boundary[seq(2, d, 2)] - a$boundary[seq(1, d, 2)]  # for some reason boundary is sequential
    C <- (sqrt(sum(span^2)) / 2)^(1/p)
    # so that tunneling of points (deleting units of mass and adding them again somewhere else) is
    # never worth it. Notice that since pixels are inside the boundary C is a bit larger than 
    # theoretically necessary (for numerical reasons)
  }
  
  # if p==1 remove pointwise minimum to legally get more zero mass points
  if (p == 1) {
    minab <- pmin(a$mass,b$mass)
    ared <- a$mass - minab  
    bred <- b$mass - minab
    atotmass <- sum(ared)
    btotmass <- sum(bred)
  } else {
    ared <- a$mass
    bred <- b$mass
    atotmass <- a$totmass
    btotmass <- b$totmass
  }
  wha <- ared > 0
  whb <- bred > 0
  apos <- ared[wha]
  bpos <- bred[whb]
  m <- sum(wha)
  n <- sum(whb)
  # The following catches the case that after the reduction procedure nothing is left, i.e. the two measures were equal
  if (m==0 && n==0) {
    if (output == "dist") {
      return(0)
    }
    if (output == "rawres") {
      NA
    }
    n1 <- a$n[1]
    n2 <- a$n[2]
    if (p != 1) {
      minab <- matrix(0, n1, n2)
    }
    return(outputallzero(ared, bred, n1, n2, inplace=minab, p=p, C=C)) # I think for p!=0 this is only called if a or b = 0-measure
  }

  gg <- as.matrix(expand.grid(a$generator)) 
    # think of x/y-coordinates given in (colmajor) *matrix enumeration* order
    # since the grid is always regular this wrong order does not matter, even with the selection we do in the next command 
  
  if (p == 1) {
    costm <- gen_cost(gg[wha, , drop=FALSE], gg[whb, , drop=FALSE], threads=threads)^(1/2)
    ltunnel <- costm >= 2*C
    costm[ltunnel] <- 2*C
    acan <- bcan <- FALSE # are there trashcan states at the end of a and b (and the rows and cols of costm)
    aplus <- c(apos)  # matrix to vector
    bplus <- c(bpos)
    wha <- c(wha)
    whb <- c(whb)
    if (!isTRUE(all.equal(atotmass, btotmass))) {
      if (atotmass > btotmass) {
        costm <- cbind(costm, C^p)  # we do *not* divide by two! (MSM not HKM)
        bplus <- c(bplus, atotmass-btotmass)  # matrix to vector
        bcan <- TRUE
        whb <- c(whb, TRUE)  
        ltunnel <- cbind(ltunnel, TRUE)
      } else {
        costm <- rbind(costm, C^p)  # we do *not* divide by two! (MSM not HKM)
        aplus <- c(aplus, btotmass-atotmass)
        acan <- TRUE
        wha <- c(wha, TRUE)
        ltunnel <- rbind(ltunnel, TRUE)
      }
    }
  } else {
    costm <- gen_cost(gg[wha, , drop=FALSE], gg[whb, , drop=FALSE], threads=threads)^(p/2)  
      # could take pmin with 2*C^p, but not needed, let's see what is easier (since most probably the pmin is cheap)
    costm <- rbind(costm, C^p)  # we do *not* divide by two! (MSM not HKM)
    costm <- cbind(costm, C^p)
    costm[m+1, n+1] <- 0
    aplus <- c(apos, btotmass)  # everything gets linearized in colmajor order
    bplus <- c(bpos, atotmass)
    if (!isTRUE(all.equal(sum(aplus), sum(bplus)))) {
      warning("Substantially different mass vectors after extending. ", sum(aplus), " != ", sum(bplus))
    }
    acan <- bcan <- TRUE
    wha <- c(wha, TRUE)  # matrix to vector
    whb <- c(whb, TRUE)  # matrix to vector
  }

  if (method == "networkflow") {
    rawres <- networkflow(matrix(aplus), matrix(bplus), costm, threads=threads) # m x 1 and n x 1 matrices for the masses 
    # sanity check (not sure if already performed inside networkflow)
    primalcost <- sum(costm * rawres$plan)   
    dualcost <- sum(rawres$potential * c(aplus, bplus))
    if (!isTRUE(all.equal(primalcost, dualcost)) || !isTRUE(all.equal(primalcost, rawres$dist))) {  # dist is cost (dist^p)
      warning("Primal-dual gap is ", rawres$dist - dualcost, "\n", "Primal cost: ", primalcost, 
              "; dual cost: ", dualcost, "; rawres$cost: ", rawres$dist)
    }
  } else {
    rawres <- unbalanced_revsimplex_core(aplus, bplus, costm, p, C)   
    # this rawres does not have a component potential (to be fixed)
  }
  
  if (output == "dist") {
    return(rawres$dist^(1/p)) # rawres$dist is the p-th power of the unbalanced Wasserstein dist
  }
  
  # emulates the output of transport with networkflow and fullreturn=TRUE (trashcan states added)
  if (output == "rawres") { 
    rawres$frame <- rawres$frame[rawres$frame[,3]>0,,drop=FALSE]
    if (a$N > m || b$N > n) {
      rawres <- zero_transform(rawres, wha, whb, wha, whb) 
         # this is a bit of a hack, so we can use the same zero_transform function (wha, whb just
         # happen to have the right length and nothing else is needed)
    }
    df <- data.frame(from=rawres$frame[,1], to=rawres$frame[,2], mass=rawres$frame[,3])
    out <- list(default=df, primal=rawres$plan, dual=rawres$potential, cost=rawres$dist)
    return(out)
  }
  
  # output = "all"
  if (p == 1) {
    temp <- rawres$plan * ltunnel
    rawres$aextra <- rowSums(temp[1:m,,drop=FALSE])  # there may or may not be a (m+1)-st row
    rawres$bextra <- colSums(temp[,1:n,drop=FALSE])  # there may or may not be a (n+1)-st column
    rawres$plan[ltunnel] <- 0
    ind <- which(rawres$plan > 0, arr.ind=TRUE) 
    rawres$frame <- cbind(ind, rawres$plan[ind])   # trashcan states (if there were any) will not appear here
                                                   # but we leave them in plan and in potential (currently forever)
    colnames(rawres$frame) <- NULL
  } else { # note: for p >= 2 we may or may not toss if transport is at distance exactly 2*C^p (for p=1 we always toss)
    rawres$aextra <- rawres$plan[1:m,n+1,drop=FALSE]
    rawres$bextra <- rawres$plan[m+1,1:n,drop=FALSE]
    select <- (rawres$frame[,3] > 0) & # we keep transports over dist 0 in current version (rawres$frame[,1] == rawres$frame[,2]) 
              (rawres$frame[,1] <= m) & (rawres$frame[,2] <= n)  # removes the trashcan states from frame
    rawres$frame <- rawres$frame[select, , drop=FALSE]
  }
  
  if (a$N > m || b$N > n){ # if any omission of zero mass points took place
    # (it seems that changing whx to indx <- which(whx) from the beginning would be
    # clearer and uses less memory for sparse images) 
    result <- zero_transform_unbalanced(rawres, wha, whb, a$n[1], a$n[2], p)  
  } else {
    result <- list(dist=rawres$dist^(1/p), plan=rawres$frame)
    result$aextra <- matrix(rawres$aextra, a$n[1], a$n[2])
    result$bextra <- matrix(rawres$bextra, a$n[1], a$n[2])
  }
  
  # fill result$inplace 
  if (p == 1) {
    result$inplace <- minab
  } else {
    result$inplace <- matrix(0,a$n[1],a$n[2])
    ind <- result$plan[,1] == result$plan[,2]
    result$inplace[result$plan[ind,1]] <- result$plan[ind,3]
    result$plan <- result$plan[result$plan[,1] != result$plan[,2], , drop=FALSE]  
    # for p != 1 we lose the transports over dist 0 only here, that's why atrans and btrans have to come after
  }
  
  # fill result$atrans 
  atemp <- rowsum(result$plan[,3], result$plan[,1])
  where <- as.numeric(attr(atemp, "dimnames")[[1]])
  result$atrans <- matrix(0, a$n[1], a$n[2])
  result$atrans[where] <- atemp
  # fill result$btrans 
  btemp <- rowsum(result$plan[,3], result$plan[,2])
  where <- as.numeric(attr(btemp, "dimnames")[[1]])
  result$btrans <- matrix(0, b$n[1], b$n[2])
  result$btrans[where] <- btemp
  
  tol <- ifelse(method == "networkflow", sqrt(.Machine$double.eps), 1e-7)  
    # for networkflow standard tolerance, for revsimplex somewhat smaller due to
    # smaller precision (too many spurious warnings otherwise)
  if (!isTRUE(all.equal(result$atrans + result$aextra + result$inplace, a$mass, tolerance=tol, check.attributes = FALSE))) {
    warning("atrans, aextra and inplace do not sum up to a$mass. ", 
            all.equal(result$atrans + result$aextra + result$inplace, a$mass))
            # gives mean relative difference |left-right|/left, no sign!!
  }
  if (!isTRUE(all.equal(result$btrans + result$bextra + result$inplace, b$mass, tolerance=tol, check.attributes = FALSE))) {
    warning("btrans, bextra and inplace do not sum up to b$mass. ", 
            all.equal(result$btrans + result$bextra + result$inplace, b$mass))
            # gives mean relative difference |left-right|/left, no sign!!
  }
    
  result$plan <- data.frame(from=result$plan[,1], to=result$plan[,2], mass=result$plan[,3])
  result <- result[c("dist", "plan", "atrans", "btrans", "aextra", "bextra", "inplace")]
  attr(result, "a") <- a
  attr(result, "b") <- b
  attr(result, "p") <- p
  attr(result, "C") <- C
  class(result) <- "ut_pgrid"
  return(result)
}



# revsimplex (also used for unbalanced.wpp)
unbalanced_revsimplex_core <- function(aplus, bplus, costm, p, C) {
  m <- length(aplus)
  n <- length(bplus)
  asum <- sum(aplus)
  bsum <- sum(bplus)
  stopifnot(isTRUE(all.equal(asum, bsum)))   # this should have been arranged by putting mass at trashcan state(s)
  
  # turn aplus, bplus into (pseudo)integer vectors if they aren't already
  is.naturalzero <-
  function(x, tol = .Machine$double.eps^0.5)  all((abs(x - round(x)) < tol) & x > -0.5)
  fudgeN <- fudgesum <- 1 
  if (!is.naturalzero(aplus) || !is.naturalzero(bplus)) {
    fudgeN <- 1e9
    fudgesum <- asum
    aplus <- round(aplus/asum * fudgeN)
    bplus <- round(bplus/bsum * fudgeN)
    aplus <- fudge(aplus,fudgeN)
    bplus <- fudge(bplus,fudgeN)
  }
  
  # initialization for computing starting solution in C code (by modrowmin method)
  initassig <- rep(0L,m*n)
  initbasis <- rep(0L,m*n)
  startgiven <- 0
  temp <- .C("revsimplex", as.integer(m), as.integer(n), as.integer(aplus), as.integer(bplus),
                  as.double(costm), assignment = as.integer(initassig), basis = as.integer(initbasis),
                  startgiven = as.integer(0), DUP=TRUE, PACKAGE="transport")

  temp$assignment <- temp$assignment * fudgesum / fudgeN
  nbasis <- sum(temp$basis)
  rawres <- list(frame = data.frame(from = rep(0,nbasis), to = rep(0,nbasis), mass = rep(0,nbasis)),
              plan = matrix(temp$assignment, m, n),
              potential = NA, dist = sum(temp$assignment * costm))  
    # dist is actually cost (no ^(1/p)) as always in rawres
    # we do not have easy access to the dual solution, but we simply do not
    # return it for rawres with revsimplex, and for computing the other output options
    # in unbalanced, we do not need it 
  
  ind <- which(matrix(as.logical(temp$basis),m,n), arr.ind=TRUE) 
  rawres$frame$from <- ind[,1]
  rawres$frame$to <- ind[,2]
  rawres$frame$mass <- temp$assignment[(ind[,2]-1)*m + ind[,1]]

  return(rawres)
}



#' Plot Unbalanced Transport Information
#'
#' Graphic representation of components of the list returned by \code{\link{unbalanced}}.
#'
#' @param x the list returned by \code{\link{unbalanced}} with option \code{output="all"}.
#' @param what character. The aspect of the unbalanced transport information to display.
#' @param axes logical. Whether to plot axes (ignored for \code{what="plan"}).
#' @param ... further graphics parameters passed to \code{\link{plot.pgrid}} for 
#' \code{what="plan"} and passed to \code{\link{matimage}} in all other cases.
#'
#' @return Nothing. Used for the side effect.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- unbalanced( random32a, random32b, p=1, C=0.2, output="all" )
#' plot( res, what="plan", lwd=1.5, angle=20 )
#' plot( res, what="trans" )
#' plot( res, what="extra" )
#' plot( res, what="inplace" )}
# in the long run we might do something fancier here (e.g. that depicts all the info in one plot,
# and we should definitely switch to ggplot2)
plot.ut_pgrid <- function(x, what=c("plan", "extra", "trans", "inplace"), axes=FALSE, ...) {
  stopifnot(is(x, "ut_pgrid"))
  what <- match.arg(what)
  
  a <- attr(x, "a")
  b <- attr(x, "b")
  
  if (what == "plan") {
    plot.pgrid(a, b, x$plan, rot=TRUE, ...)
  }
  
  if (what == "extra") {
    extramax <- max(max(x$aextra), max(x$bextra))
    matimage(x$aextra, x=a$generator[[1]], y=a$generator[[2]], zlim=c(1e-7*extramax, extramax), 
             col = hcl.colors(128, "Reds", alpha=1, rev=TRUE), xlab="", ylab="", axes=axes, ...)
    title(expression("Extra mass (" * phantom("a red") * ", "*phantom("b blue") * ")"), col.main = "black")
    title(expression(phantom("Extra mass (") * "a red" * phantom(", b blue)")), col.main = "#CD1B36FF")
    title(expression(phantom("Extra mass (a red, ") * "b blue", phantom(")")), col.main = "#3674B9FF")
    matimage(x$bextra, x=a$generator[[1]], y=a$generator[[2]], zlim=c(1e-7*extramax, extramax),
             col = hcl.colors(128, "Blues", alpha=1, rev=TRUE), add=TRUE)
  }
  
  if (what == "trans") {
    # here we only plot the difference (a-b pos -> red, b-a pos -> blue,
    # for p=1, it could be (or is!?) arranged that at most one of atrans, btrans is non-null at any given position
    if (any(x$atrans > 0 & x$btrans > 0)) {   
      temp <- x$atrans - x$btrans
      atrans <- pmax(temp, 0)
      btrans <- pmax(-temp, 0)
      net <- TRUE
    } else {
      atrans <- x$atrans
      btrans <- x$btrans
      net <- FALSE
    }
    transmax <- max(max(atrans), max(btrans))
    matimage(atrans, x=a$generator[[1]], y=a$generator[[2]], zlim=c(1e-7*transmax, transmax),
             col = hcl.colors(128, "Reds", alpha=1, rev=TRUE), xlab="", ylab="", axes=axes, ...)
    if (net) {
      title(expression("Net transported mass (" * phantom("a red") * ", "*phantom("b blue") * ")"), col.main = "black")
      title(expression(phantom("Net transported mass (") * "a red" * phantom(", b blue)")), col.main = "#CD1B36FF")
      title(expression(phantom("Net transported mass (a red, ") * "b blue", phantom(")")), col.main = "#3674B9FF")
    } else {
      title(expression("Transported mass (" * phantom("a red") * ", "*phantom("b blue") * ")"), col.main = "black")
      title(expression(phantom("Transported mass (") * "a red" * phantom(", b blue)")), col.main = "#CD1B36FF")
      title(expression(phantom("Transported mass (a red, ") * "b blue", phantom(")")), col.main = "#3674B9FF")
    }
    matimage(btrans, x=a$generator[[1]], y=a$generator[[2]], zlim=c(1e-7*transmax, transmax),
             col = hcl.colors(128, "Blues", alpha=1, rev=TRUE), add=TRUE)
  }
  
  if (what == "inplace") {
    inplacemax <- max(x$inplace)
    matimage(x$inplace, x=a$generator[[1]], y=a$generator[[2]], zlim=c(1e-8, inplacemax),
             col = hcl.colors(128, "Purples", alpha=1, rev=TRUE), xlab="", ylab="", axes=axes, ...)
    title(expression("Mass staying in place"))
  }
  invisible()
}

# TO DO: print method for ut_pgrid
