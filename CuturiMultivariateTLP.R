CuturiMultivariateTLP <- function(mu,
                                  f,
                                  nu,
                                  g,
                                  p = 2,
                                  eps = 0.01,
                                  iter = 200) {

# Domain Size
M <- ncol(f)  

# Discretize domain
timeDomain <- seq.int(M)/M


# Compute function difference
functionDiff <- matrix(0, ncol = ncol(f), nrow = ncol(f))
for (i in seq_len(nrow(f))) {
  functionDiff <- functionDiff + outer(f[i, ], g[i, ], "-") 
}

# Compute cost matrix
C <- as.matrix(abs(dist(timeDomain, method = "euclidean", upper = T, diag = T, p = p))) + (functionDiff ^ p)/M 

# Gibbs energu
K <-  exp(- C/eps)

v <- matrix(1, nrow = M, ncol = 1)
u <- matrix(0, nrow = M, ncol = 1)
crit <- vector(mode = "numeric", length = M)

phi <- matrix(1/M, nrow = M, ncol = M)
entropy <- - sum(phi * log(phi))
crit[1] <- sum(K * phi) - eps * entropy

for (j in 2:iter) {
  
  u <- mu / (K %*% v)
  v <- nu / (t(K) %*% u)
  
  phi <- diag(as.vector(u)) %*% K %*% diag(as.vector(v))
  
  entropy <- - sum(phi * log(phi))
  crit[j] <- sum(K * phi) - eps * entropy
  
  if (is.na(crit[j])) {
    message("Numerical issues, increase epsilon")
    return(NULL)
  }
  
  message(j)
  if (abs((crit[j]/crit[j-1] - 1)) < 10^(-4)) {
    Sinkhorn <-  eps * sum(phi * log(phi/K))
    TransportMap <-  phi/mu
    d <-  sum(u * (K * C) %*% v)
    return(list(Sinkhorn = Sinkhorn, TransportMap = TransportMap, distTLP = d, crit = crit, phi = phi))
  } else if (j == iter) {
    message("Not converged, increase iter")
    return(NULL)
  }

  
}

}

linearEmbedTLP <- function(listTransportMap,
                           data,
                           sigma,
                           h,
                           p = 2) {

  N <- nrow(listTransportMap[[1]]$TransportMap)
  t <- as.vector(seq.int(0, N - 1)/N) # domain
  m <- nrow(data[[1]])

  E <- sapply(listTransportMap, function(x) matrix((x$TransportMap - diag(1, N)), N, N) %*% as.vector(t) * sigma ^ (1/p))
  coord <- lapply(listTransportMap, function(x) t(matrix(rep(x$TransportMap %*% as.vector(t), m), N, m)))
  data <- Map("*", data, coord)
  funE <- sapply(data, function(x) (x  - h) * sigma ^ (1/p))
  linearEmbed <- rbind(E, funE)

  return(linearEmbed)
}    


## simulate from mulativariate AR models
require(mAr)

armodel1 <- function(numSamples){
  y <- array(0, c(3, 50, numSamples))
  
  for (i in seq.int(numSamples)) {
  w = c(0.25, 0.1, 0)
  C = rbind(c(0.2, 0.1, 0.1), c(0.1, 0.2, 0.1), c(0.1, 0.1, 0.2))
  A = rbind(c(0.4, 1.2, 0.35, -0.3, 0.2, 0.3),c(0.3, 0.7, -0.4, -0.5, 0.1, 0.1), c(0.1, 0.5, -0.2, -0.9, 0.3, 0.2))
  x = mAr.sim(w, A, C, N=50)
  
  y[ , ,i] <- t(as.matrix(x))
  }
  return(y)
}

armodel2 <- function(numSamples){
  y <- array(0, c(3, 50, numSamples))
  
  for (i in seq.int(numSamples)) {
  w = c(0.5, -0.5, 1)
  C = rbind(c(0.2, 0.1, 0.1), c(0.1, 0.2, 0.1), c(0.1, 0.1, 0.2))
  A = rbind(c(0.5, 1.2, 1, -0.3, 0.1, 0.1),c(0.5, 0.7, -0.4, -0.5, 0.1, 0.5), c(0.2, 0.5, -0.2, -0.9, 0.3, 0.2))
  x = mAr.sim(w, A, C, N=50)
  
  y[, , i] <- t(as.matrix(x))
  }
  return(y)
}

addNoise <- function(numSamples,
                     myFun) {
      y <- array(0, c(nrow(myFun), ncol(myFun), numSamples))
      
      for (i in seq.int(numSamples)) {
        y[, , i] <- myFun + rnorm(n = length(myFun), 0, 5)
      }
  return(y)
}