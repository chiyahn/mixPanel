#' @description Returns a valid list that represents parameters of a random MS-AR model with
#' non-switching beta and switching sigma; ideal for creating a sample for testing.
#' @export
#' @title GenerateMDPTheta
#' @name GenerateMDPTheta
#' @param M The number of components.
#' @param s The number of terms used for AR(s)
#' @param p The dimension of regressors that depend on components.
#' @param q The dimension of regressors that are independent of components.
#' @return A list that represents the parameters of a model with items:
#' \item{alpha}{M by 1 column that represents a mixing probability}
#' \item{rho}{s by 1 column for state-independent coefficients on AR(s)}
#' \item{beta}{p by M matrix that contains switching
#' coefficients for state-dependent exogenous variables}
#' \item{gamma}{q by 1 column that contains non-switching
#' coefficients for state-independent exogenous variables}
#' \item{mu}{M by 1 column that contains state-dependent mu}
#' \item{sigma}{M by 1 column that contains state-dependent sigma}
GenerateMDPTheta <- function(M = 2, s = 0, p = 0, q = 0)
{
  alpha <- runif(M)
  alpha <- alpha / sum(alpha)

  rho <- NULL
  beta <- NULL
  gamma <- NULL
  mu <- runif(M, 0.4, 0.8) * sample(c(1,-1), M, replace = T)
  sigma <- runif(M, 0.3, 1.5)
  if (s > 0)
  {
    rho <- runif(s, 0.3, 0.8) * sample(c(1,-1), s, replace = T)
    rho <- rho / (1.2 * sum(abs(rho)))
    if (M > 1)
      for (i in 1:(M-1))
      {
        rho.col <- runif(s, 0.3, 0.8) * sample(c(1,-1), s, replace = T)
        rho.col <- rho.col / (1.2 * sum(abs(rho.col)))
        rho <- cbind(rho, rho.col)
      }
    rho <- as.matrix(rho)
  }
  if (p > 0)
  {
    beta <- matrix(runif(p * M, -0.4, 0.4), ncol = M)
    beta <- 0.8 * (beta / sum(abs(beta)))
  }
  if (q > 0)
  {
    gamma <- runif(q, -0.4, 0.4)
    gamma <- 0.8 * (gamma / sum(abs(gamma)))
  }
  theta <- list(alpha = alpha,
                rho = rho,
                beta = beta,
                gamma = gamma,
                mu = mu,
                sigma = sigma)
  return (theta)
}

#' Generates a sample of N individuals with length T from a model specification given.
#' @export
#' @title GenerateMDPSample
#' @name GenerateMDPSample
#' @param theta A list that represents the parameters of a model.
#' @param N The number of individuals.
#' @param T The number of sample observations for each individual.
#' @param initial.y.set n_initial by N matrix that represents previous samples;
#' n_initial must be larger than/equal to s, the number of autoregressive terms.
#' By default, initial.y.set is going to be determined by rnorm(s) for
#' each individual.
#' @param x n by q matrix of data for exogenous variables that
#' have switching coefficeints.
#' @param z n by p matrix of data for exogenous variables that
#' have non-switching coefficeints.
#' @return A list with items:
#' \item{y}{(T + length(initial.y.set)) by N matrix that represents a sample
#' appended with previous values used to estimate autoregressive terms}
#' \item{y.sample}{T by N matrix that represents a sample of the model}
#' \item{y.lagged}{T by (s * N) blocked matrix that represents a lagged sample
#' of the model, whose ith block consisting of s columns gives a lagged sample
#' of ith individual, whose kth column in ith block represents a kth lagged column}
#' \item{components}{M by 1 column that contains state-dependent sigma}
#' \item{MDP.model}{An instance in MDP.model that represents the actual model
#' used to create the sample.}
#' @examples
#' GenerateMDPSample()
#' theta <- GenerateMDPTheta(M = 2, s = 3)
#' GenerateMDPSample(theta)
#' GenerateMDPSample(theta, N = 200)
GenerateMDPSample <- function(theta = NULL, N = 40, T = 5,
                           initial.y.set = NULL,
                           x = NULL, z = NULL)
{
  if (is.null(theta))
    theta <- GenerateMDPTheta()
  M <- length(theta$alpha)
  mu <- as.matrix(theta$mu)
  sigma <- as.matrix(theta$sigma)
  rho <- matrix(rep(0,M), ncol = M)
  beta <- matrix(rep(0,M), ncol = M)
  gamma <- as.matrix(0)
  
  if (!is.null(theta$rho))
    rho <- as.matrix(theta$rho)
  
  s <- nrow(rho)
  p <- 1
  q <- 1

  if (is.null(initial.y.set))
    initial.y.set <- matrix(rnorm((N * s)), ncol = N)

  if (nrow(initial.y.set) < s)
    stop ("EXCEPTION: The length of initial.y.set cannot be smaller than s.")
  if (nrow(initial.y.set) > s)
    warning ("The length of initial.y.set is greater than s;
             only the last s observations are going to be used for
             sample generation.")


  if (is.null(x))
    x <- matrix(rep(0,((s + T) * N)), ncol = N)
  else
  {
    x <- as.matrix(x)
    if (T > ncol(x))
      stop("EXCEPTION: the number of observations for each individual in
            x cannot be smaller than the length of samples, T.")
    beta <- as.matrix(theta$beta)
    
    p <- nrow(beta)
    pN <- p * N
    x <- rbind(matrix(rep(NaN, (nrow(initial.y.set)*pN)), ncol = pN),
               as.matrix(x[(nrow(x) - T + 1):(nrow(x)),]))

  }
  if (is.null(z))
    z <- matrix(rep(0,((s + T) * N)), ncol = N)
  else
  {
    z <- as.matrix(z)
    if (T > ncol(z))
      stop("EXCEPTION: the number of observations for each individual in
            z cannot be smaller than the length of samples, T.")
    gamma <- as.matrix(theta$gamma)

    q <- nrow(gamma)
    qN <- q * N
    z <- rbind(matrix(rep(NaN, (nrow(initial.y.set)*qN)), ncol = qN),
               as.matrix(z[(nrow(z) - T + 1):(nrow(z)),]))
  }

  # initialization
  initial.y.set <- as.matrix(initial.y.set)
  initial.y.set <- matrix(initial.y.set[(nrow(initial.y.set) - s + 1):
                                  nrow(initial.y.set),],
                          ncol = N) # only last s
  y <- rbind(initial.y.set, matrix(rep(-Inf, (N * T)), ncol = N))

  # determines components
  trans.cumsum <- cumsum(theta$alpha)
  components <- rep(1, N)
  probs <- runif(N)
  for (j in 2:M)
    for (i in 1:N)
      if (probs[i] > trans.cumsum[j-1] && probs[i] <= trans.cumsum[j])
        components [i] <- j

  initial.index <- nrow(initial.y.set) + 1
  last.index <- nrow(initial.y.set) + T
  for (i in 1:N)
  {
    component <- components[i]
    x.block.first <- 1 + (i - 1) * p
    z.block.first <- 1 + (i - 1) * q
    x.block <- matrix(x[,x.block.first:(x.block.first + p - 1)], ncol = p)
    z.block <- matrix(z[,z.block.first:(z.block.first + q - 1)], ncol = q)
    for (t in initial.index:last.index)
      y[t,i] <- rnorm(1, mu[component,1], sd=sigma[component,1]) +
        t(rev(y[(t-s):(t-1),i])) %*% as.numeric(rho[,component]) +
        x.block[t,] %*% as.matrix(beta[,component]) +
        z.block[t,] %*% as.matrix(gamma)
  }

  model <- list(theta = theta,
                y = y,
                 log.likelihood = 1,
                 aic = Inf, bic = Inf,
                 components = components,
                 label = "MDP.model")
  
  s <- ifelse(is.null(theta$rho), 0, nrow(theta$rho)) # take original rho back

  if (!(s > 0))
    y <- y[2:nrow(y),]
  lagged.and.sample <- apply(y, 2, GetLaggedAndSample, s)
  y.sample <- sapply(lagged.and.sample, "[[", "y.sample")
  y.lagged <- matrix(sapply(lagged.and.sample, "[[", "y.lagged"), nrow = T)
  x <- as.matrix(x[(nrow(x) - T + 1):(nrow(x)),])
  z <- as.matrix(z[(nrow(z) - T + 1):(nrow(z)),])
  
  return (list(y = y,
               x = x,
               z = z,
               y.sample = y.sample,
               y.lagged = y.lagged,
               components = components,
               MDP.model = model))
}
