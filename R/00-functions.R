library(tidyverse)
library(Matrix)

# linearly scale matrix elements to have a specified minimum and maximum value
scale_matrix = function(mat, scale_range = c(-1,1)) {
  stopifnot(length(scale_range) == 2,
            scale_range[1] < scale_range[2])
  minval = scale_range[1]
  maxval = scale_range[2]
  mat = (mat - min(mat)) / diff(range(mat))
  return(mat * (maxval - minval) + minval)
}

# simulate spatially correlated 2d field by some random fourier modes
simulate_spatial_noise = function(nx, ny, 
                                  nmodes_x=nx, nmodes_y=ny, 
                                  scale_range=NULL) {
  stopifnot(nmodes_x > 0, nmodes_x <= nx, 
            nmodes_y > 0, nmodes_y <= ny, 
            nx > 0, ny > 0)
  mat  = matrix(0, nx, ny)
  mat[1:nmodes_x, 1:nmodes_y] =      rnorm(nmodes_x * nmodes_y) + 
                                1i * rnorm(nmodes_x * nmodes_y)

  mat = Re(fft(mat, inverse=TRUE))
  if (!is.null(scale_range)) {
    mat = scale_matrix(mat, scale_range)
  }
  return(mat)
}

# structure matrix of 2d rw on latitude longitude grid
make_rw2d_structure_matrix = function(nlon, nlat) {
  Tcol = bandSparse(n=nlat, k=0:1, sym=TRUE, 
                    diag=list(c(1, rep(2, nlat-2), 1),
                              rep(-1, nlat-1)))
  Trow = bandSparse(n=nlon, k=0:1, sym=TRUE, 
                    diag=list(c(1, rep(2, nlon-2), 1),
                              rep(-1, nlon-1)))
  TT = kronecker(Diagonal(nlon), Tcol) + kronecker(t(Trow), Diagonal(nlat))
  RR = crossprod(TT)
  return(RR)
}


################################################################################
# functions to infer hyperparameters kappa
################################################################################

# log p(y|x,kappa) as function of x
lp_y_I_x = function(x, params) with(params, {
  alpha = x[1:SS + 0 * SS]
  beta  = x[1:SS + 1 * SS]
  tau   = x[1:SS + 2 * SS]
  - TT * sum(tau) +
  - TT / 2 * sum(exp(-2*tau) * (vyy + beta^2*vzz - 2 * beta * vyz + 
                                (alpha - ybar)^2))
})

# log p(x|kappa) as function of x
lp_x_I_kappa = function(x, kappa, params) with(params, {
  alpha = x[1:SS + 0 * SS]
  beta  = x[1:SS + 1 * SS]
  tau   = x[1:SS + 2 * SS]
  ans = SS /2 * sum(log(kappa)) +
        - 1/2 * (kappa['alpha'] * drop(alpha %*% crossprod(RR, alpha))) +
        - 1/2 * (kappa['beta'] * drop(beta %*% crossprod(RR, beta))) +
        - 1/2 * (kappa['tau'] * drop(tau %*% crossprod(RR, tau)))
  return(unname(ans))
})

# Cholesky of posterior precision
chol_Qxy = function(kappa, params) with(params, {
  Q_x = bdiag(list(kappa['alpha'] * RR, 
                   kappa['beta'] * RR, 
                   kappa['tau'] * RR))
  cQ_xy = Cholesky(Q_x + nH_y, LDL=FALSE) 
  return(cQ_xy)
})

# log posterior of kappa given y as function of kappa (up to an additive constant)
lp_kappa_I_y = function(log_kappa, params) {
  kappa = exp(log_kappa)
  cQ_xy = chol_Qxy(kappa, params)
  xstar = drop(solve(cQ_xy, params$nH_y %*% params$x_hat))
  ldethalf = sum(log(diag(as(cQ_xy, 'sparseMatrix'))))
  lp_y_I_x(xstar, params) + lp_x_I_kappa(xstar, kappa, params) - ldethalf
}


# function that generates samples from p(kappa | y)
# params = list(SS=n_lon * n_lat, TT=n_t, 
#               vyy=hindcast_summary$vyy, vzz=hindcast_summary$vzz, 
#               vyz=hindcast_summary$vyz, ybar=hindcast_summary$ybar, 
#               RR=RR, nH_y=nH_y, x_hat=x_hat)
#
# lp_lkappa = function(lkappa) {
#   lp = sum(-lambda * exp(lkappa))
#   return(lp)
# }

makefunction_sample_lkappa_I_y = function(params, lp_lkappa) {

  force(params)
  force(lp_lkappa)

  # function to optimize 
  fn_opt = function(lkappa) {
    lp = lp_kappa_I_y(log_kappa = lkappa, params = params) + 
         lp_lkappa(lkappa, params$lambda)
    return(lp)
  }
  
  opt_lkappa = optim(par=c(alpha=-10, beta=-10, tau=-10), 
                     fn=fn_opt, control=list(fnscale=-1))
  
  # find hessian of the log posterior at the mode for laplace approximation
  # (multiply by minus because fn_opt is the negative ps
  hess_lkappa = optimHess(par=opt_lkappa[['par']], fn = fn_opt)
  
  # sample kappa values from approximated posterior
  mu_lkappa = opt_lkappa$par
  sigma_lkappa = solve(-hess_lkappa)

  f = function() {
    drop(mvtnorm::rmvnorm(1, mean=mu_lkappa, sigma=sigma_lkappa))
  }

  return(f)
}

# function factory to generate a function that samples from the posterior of x
# given log kappa and y
#
# params = list(SS=n_lon * n_lat, TT=n_t, 
#               vyy=hindcast_summary$vyy, vzz=hindcast_summary$vzz, 
#               vyz=hindcast_summary$vyz, ybar=hindcast_summary$ybar, 
#               RR=RR, nH_y=nH_y, x_hat=x_hat)
#
makefunction_sample_x_I_lkappa_y = function(params) {

  force(params)

  # function to generate one posterior sample of x for one value of log kappa
  f = function(lkappa) with(params, {
    kappa = exp(lkappa)
    # define precision matrix
    Q_x = bdiag(list(kappa['alpha'] * RR, 
                     kappa['beta'] * RR, 
                     kappa['tau'] * RR))
    # calculate posterior mean 
    cQ_xy = Cholesky(Q_x + nH_y, LDL=FALSE) # Cholesky of posterior precision
    x_star = drop(solve(cQ_xy, nH_y %*% x_hat))
    # draw a sample from the conditional given kappa_ 
    zz = solve(cQ_xy, rnorm(n_x), system='Lt')
    x_sampl = drop(x_star + solve(cQ_xy, zz, system='Pt'))
    return(x_sampl)  
  })

  return(f)

}

