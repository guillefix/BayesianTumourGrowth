real[] ode_diff(real t,real V[2],
    real theta[6],real[] x_r,int[] x_i){
  real dVdt[2];
  dVdt[1] = V[1] * (theta[1] * (1 - V[1] / theta[2])
    - theta[3] * V[2]);
  dVdt[2] = V[2] * (theta[4] * (1 - V[2] / theta[5])
    - theta[6] * V[1])
  return dVdt;
}

data {
  real V[T,1];
}

parameters {
  real theta[6]
  real sigma
}

model {
  sigma ~ cauchy(0,1); ## Prior
  theta ~ normal(0,2); ## Prior for alpha, beta
  V0 ~ normal(5,2); ## Prior for initial #
  ## Solve ODE at current parameter values
  real V_hat[T,1];
  V_hat = integrate_ode(ode_diff, V0, t0, ts,
    theta, x_r, x_i);
  ## Likelihood
  for (t in 1:T) {
    V[t] ~ normal(V_hat[t,1],sigma);
  }
}
