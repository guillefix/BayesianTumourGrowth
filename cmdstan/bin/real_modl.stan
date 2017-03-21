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
  int N;
  int odds[N,2];
  int<lower=1> Tind;
  int<lower=1> Tsums;
  real V1[N,Tind];
  real V1stds[N,Tind];
  real V2[N,Tind];
  real V2stds[N,Tind];
  real Vsums[N,Tsums];
  real Vsumstds[N,Tsums];
  real ts_ind[Tind];
  real ts_sums[Tsums];
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
  real<lower=0> V0;
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
  for i in 1:N {
    real V_hat[T,2];
    sigma ~ cauchy(0,1); ## Prior
    theta ~ normal(1,2); ## Prior for alpha, beta
    V0 ~ normal(1,2); ## Prior for initial #
    ## Solve ODE at current parameter values
    V_hat = integrate_ode_rk45(ode_diff, c(V0*odds[i,1],V0*odds[i,2]), t0, ts,
      theta, x_r, x_i);
      ## Likelihood
      for (t in Tinds) {
        V1[i,t] ~ normal(V_hat[t,1],sigma);
        V1stds[i,t] ~ normal(0,sigma);
        V2[i,t] ~ normal(V_hat[t,2],sigma);
        V2stds[i,t] ~ normal(0,sigma);
      }
      for (t in Tsums) {
        Vsums[i,t] ~ normal(V_hat[t,1]+V_hat[t,2],sigma/sqrt(n));
        Vsumstds[i,t] ~ normal(0,sigma/sqrt(n));
      }
  }
}
