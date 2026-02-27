make_final_scope_fixture <- function(seed = 123) {
  set.seed(seed)

  time <- seq(0, 5, by = 0.05)
  nA <- 140
  nB <- 85

  lambdaA <- rgamma(nA, shape = 2, rate = 6) + 0.12
  lambdaB <- rgamma(nB, shape = 2, rate = 6) + 0.08

  S_A <- outer(lambdaA, time, function(l, t) exp(-l * t))
  S_B <- outer(lambdaB, time, function(l, t) exp(-l * t))

  xA <- as_survmat(S_A, time, group = rep("A", nA))
  xB <- as_survmat(S_B, time, group = rep("B", nB))
  x_all <- bind_survmat(xA, xB)

  list(
    time = time,
    nA = nA,
    nB = nB,
    lambdaA = lambdaA,
    lambdaB = lambdaB,
    xA = xA,
    xB = xB,
    x_all = x_all
  )
}
