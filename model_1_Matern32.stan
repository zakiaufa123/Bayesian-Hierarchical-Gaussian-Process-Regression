functions {
  matrix L_cov_matern32_nowarp(matrix X,
                               real eta_sq,
                               vector psi,
                               real sigma_sq,
                               real delta,
                               int N) {
    matrix[N, N] K;

    for (i in 1:(N - 1)) {
      K[i, i] = eta_sq + sigma_sq + delta;

      for (j in (i + 1):N) {
        real dsq;
        real r;

        dsq = psi[1] * square(X[i,1] - X[j,1])
            + psi[2] * square(X[i,2] - X[j,2]);

        // r = sqrt(dsq)
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
  vector<lower=0>[2] psi;   // bandwidth params
  real<lower=0> eta_sq;     // signal variance
  real<lower=0> sigma_sq;   // noise variance
}

model {
  matrix[N, N] L_K;
  vector[N] mu = rep_vector(0, N);

  // Priors (match Lally)
  psi      ~ gamma(4, 4);
  eta_sq   ~ student_t(4, 0, 1);  // half-t due to lower=0
  sigma_sq ~ student_t(4, 0, 1);  // half-t due to lower=0

  // GP likelihood
  L_K = L_cov_matern32_nowarp(X, eta_sq, psi, sigma_sq, 1e-6, N);
  z ~ multi_normal_cholesky(mu, L_K);
}

generated quantities {
  vector[N] z_rep;
  {
    matrix[N, N] L_K;
    vector[N] mu = rep_vector(0, N);
    L_K = L_cov_matern32_nowarp(X, eta_sq, psi, sigma_sq, 1e-6, N);
    z_rep = multi_normal_cholesky_rng(mu, L_K);
  }
}