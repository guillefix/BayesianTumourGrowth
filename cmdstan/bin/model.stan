functions {
  real[] ode_diff(real t,real[] V,
      real[] theta,real[] x_r,int[] x_i){
    real dVdt[2];
    dVdt[1] = V[1] * (theta[1] * (1 - V[1] / theta[2])
      - theta[3] * V[2]);
    dVdt[2] = V[2] * (theta[4] * (1 - V[2] / theta[5])
      - theta[6] * V[1]);
    return dVdt;
  }
}

data {
  int<lower=1> T;
  real V1[T];
  real V2[T];
  real ts[T];
  ## real V0[2];
  real t0;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real<lower=0> theta1[4];
  real lambda[2];
  real<lower=0> sigma;
  real<lower=0> V0[2];
}
transformed parameters {
  real theta[6];
  theta[1] = theta1[1];
  theta[2] = theta1[2];
  theta[3] = lambda[1];
  theta[4] = theta1[3];
  theta[5] = theta1[4];
  theta[6] = lambda[2];
}

model {
  real V_hat[T,2];
  sigma ~ cauchy(0,1); ## Prior
  theta ~ normal(1,2); ## Prior for alpha, beta
  V0 ~ normal(1,2); ## Prior for initial #
  ## Solve ODE at current parameter values
  V_hat = integrate_ode_rk45(ode_diff, V0, t0, ts,
    theta, x_r, x_i);
  ## Likelihood
  for (t in 1:T) {
    V1[t] ~ normal(V_hat[t,1],sigma);
    V2[t] ~ normal(V_hat[t,2],sigma);
  }
}
