test_that("wasserstein for wpp objects with i.i.d. weights", { 

  balanced_random <- function(nrep, method, p) {
    m <- 30
    n <- 40
    res <- rep(0,nrep)
    for (i in 1:nrep) {
      set.seed(190611+p*100+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      massx <- rexp(m)
      massx <- massx/sum(massx)
      set.seed(190611+p*400+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      massy <- rexp(n)
      massy <- massy/sum(massy)
      set.seed(190611+p*700+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      x <- wpp(matrix(runif(2*m),m,2),massx)
      set.seed(190611+p*1000+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      y <- wpp(matrix(runif(2*n),n,2),massy)
      res[i] <- wasserstein(x,y,method=method,p=p)
    }
    return(res)
  }

  # 20 each, p1 checked for revsimplex, shortsimplex, primaldual (exactly equal!)
  precomp = data.frame(
    p1 = c( 0.2204111080, 0.1844101315, 0.2971009077, 0.1886649427, 0.2790101370, 0.2546042195, 0.1611890814, 0.2074879170,
            0.2688192942, 0.2676081249, 0.2734300136, 0.2042014557, 0.1737081418, 0.4119017624, 0.1586022259, 0.2059735692,
            0.2162859769, 0.1514188963, 0.2699822686, 0.1859990710),
    p2 = c(0.4047509243, 0.2178269822, 0.1914773931, 0.2721772849, 0.2302878681, 0.3246175470, 0.3040419122, 0.2468305838,
           0.1838419748, 0.2363982393, 0.2071055851, 0.2193084666, 0.1995616592, 0.2743842796, 0.1683908060, 0.2665572555,
           0.4027684820, 0.2988371087, 0.2199253236, 0.2574573042)
  )

  expect_equal(balanced_random(nrep=20, method="revsimplex", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="shortsimplex", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="primaldual", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="networkflow", p=1), precomp$p1)
  expect_equal(balanced_random(nrep=20, method="revsimplex", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="shortsimplex", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="primaldual", p=2), precomp$p2)
  expect_equal(balanced_random(nrep=20, method="networkflow", p=2), precomp$p2)
})



test_that("wasserstein for wpp objects with non-i.d. weights", {  # the objects are unbalanced not the wasserstein dist.
  
  skip_if(getRversion() < "3.6.0", "discrete uniform generation was different for R versions prior to 3.6.0")
  
  unbalanced_random <- function(nrep, method, p) {
    m <- 30
    n <- 40
    res <- rep(0,nrep)
    for (i in 1:nrep) {
      set.seed(190611+p*100+i,kind="Mersenne-Twister",normal.kind = "Inversion",sample.kind="Rejection")
      temp <- sample(30,m,replace=TRUE)
      massx <- runif(m)*10^temp
      massx <- massx/sum(massx)
      set.seed(190611+p*400+i,kind="Mersenne-Twister",normal.kind = "Inversion",sample.kind="Rejection")
      temp <- sample(30,n,replace=TRUE)
      massy <- runif(n)*10^temp
      massy <- massy/sum(massy)
      set.seed(190611+p*700+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      x <- wpp(matrix(runif(2*m),m,2),massx)
      set.seed(190611+p*1000+i,kind="Mersenne-Twister",normal.kind = "Inversion")
      y <- wpp(matrix(runif(2*n),n,2),massy)
      res[i] <- wasserstein(x,y,method=method,p=p)
    }
    return(res)
  }

  # 20 each
  precomp = data.frame(
      p1 = c(0.408712501991695, 0.265729468180307, 0.516829831105873, 0.882572832422254, 
             0.362460109564705, 0.44352849614055, 0.524396418818738, 0.567140446957284, 
             0.199605145672852, 0.500836289468377, 0.608465149250878, 0.526292289496797, 
             0.295525516429156, 0.75368039100222, 0.539846547719169, 0.578900881963893, 
             0.512086614860834, 0.460683637826441, 0.364161062855492, 0.641518533510968),
      p2 = c(0.469550449073627, 0.490510114076552, 0.846518466977827, 0.504892622550457, 
             0.347462537086477, 0.335342245751588, 0.504182118170615, 0.449929012598699, 
             0.311890213777329, 0.225573702901285, 0.260565764246565, 0.372631582632826, 
             0.340233336357178, 0.785096673418506, 0.483406022771683, 0.855022095898563, 
             0.522595036560708, 0.600555445021545, 0.666915983059887, 0.276794690989451)
  )

  # shortlist creates warnings that are no problem, primaldual deprecated
  # --> we don't test these currently (but tests would be ok otherwise)
  expect_equal(unbalanced_random(nrep=20, method="revsimplex", p=1), precomp$p1)
  #expect_equal(unbalanced_random(nrep=20, method="shortsimplex", p=1), precomp$p1)
  #expect_equal(unbalanced_random(nrep=20, method="primaldual", p=1), precomp$p1)
  expect_equal(unbalanced_random(nrep=20, method="revsimplex", p=2), precomp$p2)
  #expect_equal(unbalanced_random(nrep=20, method="shortsimplex", p=2), precomp$p2)
  #expect_equal(unbalanced_random(nrep=20, method="primaldual", p=2), precomp$p2)
})


test_that("non-square pgrids (one toy example)", { 
  a <- pgrid(matrix(1:6, 2 ,3))
  b <- pgrid(matrix(6:1, 2 ,3))
  expect_equal( wasserstein(a, b, p=1), 0.397814379345 )
  expect_equal( wasserstein(a, b, p=1, method="revsimplex"), 0.397814379345 )
})


test_that("unbalanced transport (one toy example)", { 
  a <- pgrid(matrix(1:6, 2 ,3))
  b <- pgrid(matrix(6:1, 2 ,3))
  expect_equal( unbalanced(a, b, p=1, output="dist")/a$totmass, 0.397814379345 )
  expect_equal( unbalanced(a, b, p=1, output="dist"), 8.354101966249 )
  expect_equal( unbalanced(a, b, p=1, C=1.118034/2, output="dist"), 8.354101966249 )
  expect_equal( unbalanced(a, b, p=1, C=1.118033/2, output="dist"), 8.354099 )
  expect_equal( unbalanced(a, b, p=1, C=0.501/2, output="dist"), 4.507 )
  expect_equal( unbalanced(a, b, p=1, C=0.500/2, output="dist"), 18*0.500/2 )
  expect_equal( unbalanced(a, b, p=1, C=0.499/2, output="all")$plan, 
                data.frame(from=numeric(0), to=numeric(0), mass=numeric(0)) ) 
  
  expect_equal( unbalanced(a, b, p=2, output="dist"), 2.29128784748 )
  expect_equal( unbalanced(a, b, p=2, C=1.118034/sqrt(2), output="dist"), 2.29128784748 )
  expect_equal( unbalanced(a, b, p=2, C=1.118033/sqrt(2), output="dist"), 2.29128736502 )  
  expect_equal( unbalanced(a, b, p=2, C=0.501/sqrt(2), output="dist"), 1.502333851 )
  expect_equal( unbalanced(a, b, p=2, C=0.500/sqrt(2), output="dist")^2, 18*(0.500/sqrt(2))^2 ) 
  expect_equal( unbalanced(a, b, p=2, C=0.499/sqrt(2), output="all")$plan,
                data.frame(from=numeric(0), to=numeric(0), mass=numeric(0)) )
  
})

# a somewhat more serious example
#set.seed(221010)
#amat <- pmax(cbind(rbind(matrix(50, 3, 3), matrix(100, 3, 3)), rbind(matrix(80, 3, 3), matrix(30, 3, 3))) + rnorm(36, 0, 30), 0)
#bmat <- pmax(cbind(rbind(matrix(80, 3, 3), matrix(100, 3, 3)), rbind(matrix(30, 3, 3), matrix(100, 3, 3))) + rnorm(36, 0, 30), 0)
amat <- matrix(c(0.00000, 11.15837, 47.11794, 110.68439, 74.07142, 96.41374, 69.52885, 34.19806, 113.82058, 126.79252, 96.75680, 129.98471,
                 0.00000, 62.19007, 22.35849, 87.95412, 72.46086 ,111.66611, 56.46518, 86.01870, 87.39139, 13.56717,  0.00000, 70.25962,
                43.88404, 87.66647, 77.25338, 59.03549, 63.91189, 89.74809, 67.64778, 77.11587, 52.47348, 63.05337, 43.71204, 21.99423), 6, 6)
bmat <- matrix(c(148.10652, 77.14895, 79.83697, 121.31452, 120.33347, 82.87338, 55.72391, 56.09614, 117.18499, 97.70056, 87.84207, 107.50296,
                  72.55216, 93.82185, 89.43807, 132.24490, 62.11505, 33.63378, 50.69261, 30.38996,  0.00000, 106.08027, 62.37144, 95.76406,
                  62.44187, 74.58977, 42.91589, 152.21315, 72.02600, 89.41609, 44.87053, 64.64356, 10.08390, 108.00301, 88.42050, 131.40796), 6, 6)
a <- pgrid(amat)
b <- pgrid(bmat)
a1 <- pgrid(amat/sum(amat))
b1 <- pgrid(bmat/sum(bmat))

test_that("transport (more serious example)", { 
  expect_equal( wasserstein(a1, b1, p=1), 0.0936004132159 )
  expect_equal( wasserstein(a1, b1, p=1, method="revsimplex"), 0.093600413215893 )
  oldopt <- options("transport-CPLEX_no_warn" = TRUE)
  expect_equal( suppressWarnings(wasserstein(a1, b1, p=2, method="shielding")), 0.134995568876 )
  options(oldopt)
  expect_equal( wasserstein(a1, b1, p=2, method="networkflow"), 0.134995568876 )
  expect_equal( wasserstein(a1, b1, p=2, method="revsimplex"), 0.134995568876 )
})

test_that("unbalanced pgrid dist, same total mass (more serious example)", {
  expect_equal( unbalanced(a1, b1, p=1), 0.0936004132159 )
  expect_equal( unbalanced(a1, b1, p=2), 0.134995568876 )
  expect_equal( unbalanced(a1, b1, p=1, C=0.3), 0.0870010219265 )
  expect_equal( unbalanced(a1, b1, p=1, C=0.1), 0.0536390755165 )
  expect_equal( unbalanced(a1, b1, p=2, C=0.3), 0.134784546883 )
  expect_equal( unbalanced(a1, b1, p=2, C=0.1), 0.0764973224334 )
})


test_that("unbalanced pgrid dist, same total mass, revsimplex (more serious example)", {
  skip_on_cran()  # due to smaller precision and random initialization revsimplex 
                  # fails just marginally (probably tolerance 1e-7 would already fix it)
  expect_equal( unbalanced(a1, b1, p=1, method="revsimplex"), 0.0936004132159 )
  expect_equal( unbalanced(a1, b1, p=2, method="revsimplex"), 0.134995568876 )
  expect_equal( unbalanced(a1, b1, p=1, C=0.3, method="revsimplex"), 0.0870010219265 )
  expect_equal( unbalanced(a1, b1, p=1, C=0.1, method="revsimplex"), 0.0536390755165 )
  expect_equal( unbalanced(a1, b1, p=2, C=0.3, method="revsimplex"), 0.134784546883 )
  expect_equal( unbalanced(a1, b1, p=2, C=0.1, method="revsimplex"), 0.0764973224334 )
})
# plans are not necessarily the same, so all entries except the dist might fail here
# the only thing one can really verify is the dist
#expect_equal( unbalanced(a1, b1, p=1, C=0.3, method="revsimplex", output="all")[-2],
#              unbalanced(a1, b1, p=1, C=0.3, output="all")[-2])
#expect_equal( unbalanced(a1, b1, p=2, C=0.3, method="revsimplex", output="all")[-2],
#              unbalanced(a1, b1, p=2, C=0.3, output="all")[-2])
#expect_equal( unbalanced(a, b, p=1, C=0.3, method="revsimplex", output="all")[-2],
#              unbalanced(a, b, p=1, C=0.3, output="all")[-2])
#expect_equal( unbalanced(a, b, p=2, C=0.3, method="revsimplex", output="all")[-(2:4)],
#             unbalanced(a, b, p=2, C=0.3, output="all")[-(2:4)])


test_that("unbalanced pgrid all, different total mass (more serious example)", {
  expect_snapshot_value( unbalanced(a, b, p=1, C=0.2, output="all"), style="json2") # other styles don't seem to work
  expect_snapshot_value( unbalanced(a, b, p=2, C=0.2, output="all"), style="json2") 
  # because for total mass of a > total mass of b and p=1 used to be a problem:
  expect_snapshot_value( unbalanced(b, a, p=1, C=0.2, output="all"), style="json2") # other styles don't seem to work
  expect_snapshot_value( unbalanced(b, a, p=2, C=0.2, output="all"), style="json2") 
  # tests 1,3,4 above had recent minimal changes on my system that did not change the output distance
  # at all (accepted new snapshot on 24/03/09)
})

test_that("unbalanced pgrid all, different total mass, revsimplex (more serious example)", { 
  # expected values computed with networkflow
  skip_on_cran()  # due to smaller precision and random initialization revsimplex 
                  # fails just marginally (probably tolerance 1e-7 would already fix it)
  expect_equal( unbalanced(a, b, p=1, C=0.08, method="revsimplex"), 119.5380336 )
  expect_equal( unbalanced(a, b, p=1, C=0.2, method="revsimplex"), 206.2217782157 )
  expect_equal( unbalanced(a, b, p=1, method="revsimplex"), 507.5840082411 )
  expect_equal( unbalanced(a, b, p=2, C=0.11, method="revsimplex"), 4.25207332745, tolerance=1e-7 )
  expect_equal( unbalanced(a, b, p=2, C=0.2, method="revsimplex"), 6.33807382587 )
  expect_equal( unbalanced(a, b, p=2, method="revsimplex"), 20.88324374742 )

  # because for total mass of a > total mass of b and p=1 it's easy to introduce errors
  # due to all the tricks we do with zero mass (there used to be a problem with the tplan not the dist)
  expect_equal( unbalanced(b, a, p=1, C=0.08, method="revsimplex"), 119.5380336 )
  expect_equal( unbalanced(b, a, p=1, C=0.2, method="revsimplex"), 206.2217782157 )
  expect_equal( unbalanced(b, a, p=1, method="revsimplex"), 507.5840082411 )
  expect_equal( unbalanced(b, a, p=2, C=0.11, method="revsimplex"), 4.25207332745, tolerance=1e-7 )
  expect_equal( unbalanced(b, a, p=2, C=0.2, method="revsimplex"), 6.33807382587 )
  expect_equal( unbalanced(b, a, p=2, method="revsimplex"), 20.88324374742 )
})  

# profvis(unbalanced(random32a, random32b, p=2, output="all")$cost) # it's all in networkflow



#
test_that("semidiscrete, p=2", { 
  # introduced after reported segfault
  skip_on_cran()   # although we set a seed, the result is not consistent on CRAN
                   # *for the additional configurations*; apparently not even between
                   # two checks of flavor "OpenBLAS". Deviations are small but 
                   # it is not clear if there is a strict upper limit and what it is
  a = pgrid(matrix(c(0.7,1.5,0.8,1),2,2))
  n = 100
  set.seed(1111)
  b = wpp(matrix(runif(2*n),nrow=n), mass = rep(1/n,n))
  res <- semidiscrete(a,b,p=2)
  # For some reason things got very inexact on my system since 2024, hence the tolerances
  # On rhub all the tests are fine with the default tolerances
  # Also 
  expect_equal( res$wasserstein_dist, 0.112323776, tolerance=1e-5)
     # not only if !capabilities("long.double"), but also some other specially compiled variants
     # (including the setup CRAN uses for M1mac does not satisfy the usual tolerance limits.)
     # usually the "mean rel. difference" is of order of magnitude of (a little more than) 1e-6,
     # which is surprising (but not of concern). It goes up to 8.902738e-06 for if !capabilities("long.double") 
  expect_snapshot_value( unclass(res), style="json2", tolerance=1e-3) # other styles don't seem to work, 
                 # note that snapshot tests by default are set up to not run on cran (that's why we can keep
                 # the tolerance smaller)
})  


# to prevent further cost matrix desasters
test_that("gen_cost -- gen_cost0 -- gen_cost0d, 2d", { 
  set.seed(240509, kind="Mersenne-Twister", normal.kind = "Inversion")
  xx <- matrix(runif(20), 10, 2) 
  yy <- matrix(runif(18), 9, 2) 
  c0 <- as.matrix(dist(rbind(xx,yy)))[1:10, 11:19]^2
  dimnames(c0) <- NULL
  c1 <- transport:::gen_cost(xx, yy, 1)
  c2 <- transport:::gen_cost0(xx, yy)
  c3 <- transport:::gen_cost0d(xx, yy)
  expect_equal(c1, c0)
  expect_equal(c2, c0)
  expect_equal(c3, c0)
})  


test_that("gen_cost -- gen_cost0d, 5d", { 
  set.seed(240510, kind="Mersenne-Twister", normal.kind = "Inversion")
  xx <- matrix(runif(30), 6, 5) 
  yy <- matrix(runif(40), 8, 5) 
  c0 <- as.matrix(dist(rbind(xx,yy)))[1:6, 7:14]^2
  dimnames(c0) <- NULL
  c1 <- transport:::gen_cost(xx, yy, 1)
  c2 <- transport:::gen_cost0d(xx, yy)
  expect_equal(c1, c0)
  expect_equal(c2, c0)
})  

