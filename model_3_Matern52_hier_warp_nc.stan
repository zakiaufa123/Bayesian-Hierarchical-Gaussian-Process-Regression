functions {
  matrix cov_matern52_warp(matrix X,
                           vector psi, real eta_sq,
                           real a1, real b1,
                           real a2, real b2,
                           real sigma_sq,
                           real delta) {
    int N = rows(X);
    matrix[N, N] K;
    vector[N] w1;
    vector[N] w2;

    for (i in 1:N) {
      w1[i] = beta_cdf(X[i,1] | a1, b1);
      w2[i] = beta_cdf(X[i,2] | a2, b2);
    }

    for (i in 1:(N - 1)) {
      K[i,i] = eta_sq + sigma_sq + delta;
      for (j in (i + 1):N) {
        real dsq;
        real r;
        dsq = psi[1] * square(w1[i] - w1[j])
            + psi[2] * square(w2[i] - w2[j]);
        r = sqrt(dsq);
        K[i,j] = eta_sq
               * (1 + sqrt(5.0) * r + (5.0 / 3.0) * square(r))
               * exp(-sqrt(5.0) * r);
        K[j,i] = K[i,j];
      }
    }
    K[N,N] = eta_sq + sigma_sq + delta;
    return K;
  }
}

data {
  int<lower=1> L;                 // number of LoBs
  int<lower=1> N_total;           // total observations across LoBs
  array[L] int<lower=1> N_l;      // observations per LoB
  array[L] int<lower=1> start;    // start position in stacked arrays
  matrix[N_total, 2] X;           // stacked inputs (dev_norm, ay_norm)
  vector[N_total] z;              // stacked standardized targets
}

parameters {
  // global hierarchical hyperparameters
  real mu_psi1;
  real mu_psi2;
  real mu_eta;
  real mu_sigma;

  real<lower=0> tau_psi1;
  real<lower=0> tau_psi2;
  real<lower=0> tau_eta;
  real<lower=0> tau_sigma;

  // non-centered latent standard normals
  vector[L] z_psi1;
  vector[L] z_psi2;
  vector[L] z_eta;
  vector[L] z_sigma;

  // LoB-specific warping params
  vector<lower=0>[L] a1;
  vector<lower=0>[L] b1;
  vector<lower=0>[L] a2;
  vector<lower=0>[L] b2;
}

transformed parameters {
  // log-scale LoB-specific params
  vector[L] log_psi1;
  vector[L] log_psi2;
  vector[L] log_eta_sq;
  vector[L] log_sigma_sq;

  // natural scale
  vector<lower=0>[L] psi1;
  vector<lower=0>[L] psi2;
  vector<lower=0>[L] eta_sq;
  vector<lower=0>[L] sigma_sq;

  log_psi1     = mu_psi1  + tau_psi1  * z_psi1;
  log_psi2     = mu_psi2  + tau_psi2  * z_psi2;
  log_eta_sq   = mu_eta   + tau_eta   * z_eta;
  log_sigma_sq = mu_sigma + tau_sigma * z_sigma;

  psi1     = exp(log_psi1);
  psi2     = exp(log_psi2);
  eta_sq   = exp(log_eta_sq);
  sigma_sq = exp(log_sigma_sq);
}

model {
  // hyperpriors
  mu_psi1  ~ normal(0, 1);
  mu_psi2  ~ normal(0, 1);
  mu_eta   ~ normal(0, 1);
  mu_sigma ~ normal(0, 1);

  tau_psi1  ~ student_t(4, 0, 1);
  tau_psi2  ~ student_t(4, 0, 1);
  tau_eta   ~ student_t(4, 0, 1);
  tau_sigma ~ student_t(4, 0, 1);

  // non-centered standard normals
  z_psi1  ~ normal(0, 1);
  z_psi2  ~ normal(0, 1);
  z_eta   ~ normal(0, 1);
  z_sigma ~ normal(0, 1);

  // warping priors
  a1 ~ lognormal(0, 0.5);
  b1 ~ lognormal(0, 0.5);
  a2 ~ lognormal(0, 0.5);
  b2 ~ lognormal(0, 0.5);

  // likelihood per LoB
  for (l in 1:L) {
    int N = N_l[l];
    int s = start[l];

    matrix[N,2] Xl = block(X, s, 1, N, 2);
    vector[N] zl   = segment(z, s, N);

    vector[2] psi_l;
    matrix[N,N] K;
    matrix[N,N] L_K;

    psi_l[1] = psi1[l];
    psi_l[2] = psi2[l];

    K = cov_matern52_warp(Xl, psi_l, eta_sq[l],
                          a1[l], b1[l], a2[l], b2[l],
                          sigma_sq[l], 1e-6);
    L_K = cholesky_decompose(K);

    zl ~ multi_normal_cholesky(rep_vector(0, N), L_K);
  }
}

generated quantities {
  vector[N_total] z_rep;

  {
    for (l in 1:L) {
      int N = N_l[l];
      int s = start[l];

      matrix[N,2] Xl = block(X, s, 1, N, 2);
      vector[2] psi_l;
      matrix[N,N] K;
      matrix[N,N] L_K;

      psi_l[1] = psi1[l];
      psi_l[2] = psi2[l];

      K = cov_matern52_warp(Xl, psi_l, eta_sq[l],
                            a1[l], b1[l], a2[l], b2[l],
                            sigma_sq[l], 1e-6);
      L_K = cholesky_decompose(K);

      z_rep[s:(s + N - 1)] =
        multi_normal_cholesky_rng(rep_vector(0, N), L_K);
    }
  }
}