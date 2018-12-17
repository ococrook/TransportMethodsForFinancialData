CuturiOT <- function(mu,
                     nu,
                     p = 2,
                     eps = 0.01,
                     iter = 200) {
  
  # Domain Size
  M <- length(mu)  
  
  # Discretize domain
  timeDomain <- seq.int(M)/M

  # Compute cost matrix
  C <- as.matrix(abs(dist(timeDomain, method = "euclidean", upper = T, diag = T, p = p)))

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
      return(list(Sinkhorn = Sinkhorn, TransportMap = TransportMap, distOT = d, crit = crit, phi = phi))
    } else if (j == iter) {
      message("Not converged, increase iter")
      return(NULL)
    }
    
    
  }
  
}


linearEmbedOT <- function(listTransportMap,
                           sigma,
                           p = 2) {
  
  N <- nrow(listTransportMap[[1]]$TransportMap)
  t <- as.vector(seq.int(0, N - 1)/N) # domain

  E <- sapply(listTransportMap, function(x) matrix((x$TransportMap - diag(1, N)), N, N) %*% as.vector(t) * sigma ^ (1/p))
  linearEmbed <- E
  
  return(linearEmbed)
} 
