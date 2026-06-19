functions {
  matrix L_cov_matern32_warp(matrix X,
                             real eta_sq,
                             vector psi,
                             real a1, real b1,
                             real a2, real b2,
                             real sigma_sq,
                             real delta,
                             int N) {
    matrix[N, N] K;
    vector[N] w1;
    vector[N] w2;

    // warp inputs to (0,1) via Beta CDF
    for (i in 1:N) {
      w1[i] = beta_cdf(X[i,1] | a1, b1);
      w2[i] = beta_cdf(X[i,2] | a2, b2);
    }

    for (i in 1:(N - 1)) {
      K[i, i] = eta_sq + sigma_sq + delta;
      for (j in (i + 1):N) {
        real dsq;
        real r;

        dsq = psi[1] * square(w1[i] - w1[j])
            + psi[2] * square(w2[i] - w2[j]);

        r = sqrt(dsq);

        // Matérn 3/2: eta^2 * (1 + sqrt(3) r) exp(-sqrt(3) r)
        K[i, j] = eta_sq * (1 + sqrt(3) * r) * exp(-sqrt(3) * r);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = eta_sq + sigma_sq + delta;

    return cholesky_decompose(K);
  }
}

data {
  int<lower=1> N;
  matrix[N, 2] X;      // normalized inputs in [0,1] (dev_norm, ay_norm)
  vector[N] z;         // standardized losses
}

parameters {
  vector<lower=0>[2] psi;
  real<lower=0> eta_sq;
  real<lower=0> sigma_sq;

  // warping params
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
}

model {
  matrix[N, N] L_K;
  vector[N] mu = rep_vector(0, N);

  // priors (Lally-style)
  psi      ~ gamma(4, 4);
  eta_sq   ~ student_t(4, 0, 1);
  sigma_sq ~ student_t(4, 0, 1);

  a1 ~ lognormal(0, 0.5);
  b1 ~ lognormal(0, 0.5);
  a2 ~ lognormal(0, 0.5);
  b2 ~ lognormal(0, 0.5);

  // GP likelihood
  L_K = L_cov_matern32_warp(X, eta_sq, psi, a1, b1, a2, b2, sigma_sq, 1e-6, N);
  z ~ multi_normal_cholesky(mu, L_K);
}

generated quantities {
  vector[N] z_rep;
  {
    matrix[N, N] L_K;
    vector[N] mu = rep_vector(0, N);
    L_K = L_cov_matern32_warp(X, eta_sq, psi, a1, b1, a2, b2, sigma_sq, 1e-6, N);
    z_rep = multi_normal_cholesky_rng(mu, L_K);
  }
}