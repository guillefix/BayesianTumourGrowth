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
  real V1[T,1];
  real V2[T,1];
  real ts[T];
  real V0[2];
  real t0;
}

transformed data {
  real x_r[0];
  int x_i[0];
}

parameters {
  real theta[6];
  real sigma;
}

model {
  real V_hat[T,1];
  sigma ~ cauchy(0,1); ## Prior
  theta ~ normal(0,2); ## Prior for alpha, beta
  V0 ~ normal(5,2); ## Prior for initial #
  ## Solve ODE at current parameter values
  V_hat = integrate_ode_rk45(ode_diff, V0, t0, ts,
    theta, x_r, x_i);
  ## Likelihood
  for (t in 1:T) {
    V1[t] ~ normal(V_hat[t,1],sigma);
    V2[t] ~ normal(V_hat[t,2],sigma);
  }
}
