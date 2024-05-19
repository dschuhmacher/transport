load(file.path("_data/unbalanced.wpp-breaker.rda"))

# the following example resulted in an infinite loop before we reintroduced maxiters
test_that("wasserstein for wpp objects with i.i.d. weights", { 
  suppressWarnings({
    wasbroken <- unbalanced.wpp(a, b, p = 1, C = 0.2, method="networkflow", output="dist", threads=1)
  })   # warning about max number of iterations reached currently on my Mac
       # which examples trigger this warning depends between systems (so no expected warning)
  expect_equal(wasbroken, 0.367133547132)
  neverbroken <- unbalanced.wpp(a, b, p = 1, C = 0.2, method="revsimplex", output="dist", threads=1)
  expect_equal(wasbroken, 0.367133549274)  # not the same, but equal
    # (primaldual avoids numerical problems by using integers -> *can* only be exact for 8-9 digits due to rounding) 
})  



# the previous tests for unbalanced.wpp
# set.seed(99)
# m <- 15; n <- 20
# amass <- rexp(m); bmass <- rexp(n)
# aloc <- runif(2*m)
# bloc <- runif(2*n)
aloc <- matrix(c(0.14469, 0.31873, 0.44783, 0.93441, 0.38415, 0.14880, 0.08677, 0.56249,
                 0.23284, 0.35738, 0.38091, 0.32375, 0.90057, 0.47434, 0.59988, 0.17468,
                 0.78806, 0.94680, 0.17742, 0.57412, 0.03287, 0.17819, 0.53343, 0.73964, 
                 0.94150, 0.60528, 0.92914, 0.23820, 0.25526, 0.56154), 15, 2)
amass <- c(1.023, 3.029, 0.411, 0.260, 0.157, 0.233, 0.998, 1.444, 1.632, 0.189, 0.561, 0.122, 0.343, 0.104, 1.783)
bloc <- matrix(c(0.36033, 0.43466, 0.59690, 0.80577, 0.06677, 0.12559, 0.19023, 0.47986,
                 0.18065, 0.46663, 0.86854, 0.81021, 0.22948, 0.03031, 0.64434, 0.02060,
                 0.63880, 0.55948, 0.96709, 0.24005, 0.32036, 0.33440, 0.79957, 0.69029,
                 0.54307, 0.24627, 0.17143, 0.99342, 0.47171, 0.43714, 0.14893, 0.47446,
                 0.33290, 0.50481, 0.21834, 0.98259, 0.64568, 0.17477, 0.62269, 0.06644), 20, 2)
bmass <- c(0.115, 1.085, 0.806, 0.897, 3.509, 1.358, 0.303, 0.329, 1.787, 0.695,
           0.070, 0.783, 0.464, 0.360, 0.012, 1.410, 2.428, 0.066, 0.172, 0.304)
a <- wpp(aloc, amass)
b <- wpp(bloc, bmass)

# par(mfrow=c(1,2), mai=c(0.4,0.4,0.05,0.05))
# temp <- unbalanced(a,b, p=1, C=0, output="all")
# format(temp$dist, digits=12)
# temp2 <- unbalanced(b,a, output="all")
# plot(temp); plot(temp)
# plot(temp, what="trans"); plot(temp, what="extra")
# plot(temp); plot(temp2)
test_that("unbalanced wpp dist, different total mass", {
  expect_equal( unbalanced(a, b, p=1), 5.61509324299 )
  expect_equal( unbalanced(b, a, p=1), 5.61509324299 )
  expect_equal( unbalanced(a, b, p=2), 1.9431244082 )
  expect_equal( unbalanced(a, b, p=1, C=0.15), 3.04580570115 )
  expect_equal( unbalanced(a, b, p=1, C=0.075), 1.88499325947 )
  expect_equal( unbalanced(a, b, p=2, C=0.15), 0.676202653091 )
  expect_equal( unbalanced(a, b, p=2, C=0.075), 0.381731142931 )
})

test_that("unbalanced wpp dist, different total mass, revsimplex", {
  expect_equal( unbalanced(a, b, p=1, method="revsimplex"), 5.61509324299 )
  expect_equal( unbalanced(b, a, p=1, method="revsimplex"), 5.61509324299 )
  expect_equal( unbalanced(a, b, p=2, method="revsimplex"), 1.9431244082 )
  expect_equal( unbalanced(a, b, p=1, C=0.15, method="revsimplex"), 3.04580570115 )
  expect_equal( unbalanced(a, b, p=1, C=0.075, method="revsimplex"), 1.88499325947 )
  expect_equal( unbalanced(a, b, p=2, C=0.15, method="revsimplex"), 0.676202653091 )
  expect_equal( unbalanced(a, b, p=2, C=0.075, method="revsimplex"), 0.381731142931 )
})

test_that("unbalanced wpp dist, all-comparison, networkflow vs. revsimplex", {
  resnet <- unbalanced(a, b, p=1, output="all")
  resrev <- unbalanced(a, b, p=1, output="all", method="revsimplex")
  resrev$plan <- resrev$plan[order(resrev$plan$from, resrev$plan$to), ]
  row.names(resrev$plan) <- seq_len(dim(resrev$plan)[1])
  # the plans will usually differ
  expect_equal(resnet$plan, resrev$plan, tolerance=1e-6)
})

