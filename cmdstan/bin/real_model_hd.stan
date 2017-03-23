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
  int Ntrain;
  int Ntest;
  int n;
  int odds[Ntrain,2];
  int<lower=1> Tinds;
  int<lower=1> Tsums;
  real V1[Ntrain,Tinds];
  real V1stds[Ntrain,Tinds];
  real V2[Ntrain,Tinds];
  real V2stds[Ntrain,Tinds];
  real Vsums[Ntrain,Tsums];
  real Vsumstds[Ntrain,Tsums];

  int hdodds[Ntest,2];
  real hdV1[Ntest,Tinds];
  real hdV1stds[Ntest,Tinds];
  real hdV2[Ntest,Tinds];
  real hdV2stds[Ntest,Tinds];
  real hdVsums[Ntest,Tsums];
  real hdVsumstds[Ntest,Tsums];
  ## real ts_inds[Tinds];
  int ts_inds_inds[Tinds];
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
  real V_hat[Tsums,2];
  real V0s[2];
  int j;
  sigma ~ cauchy(0,1); ## Prior
  theta ~ normal(0.5,0.5); ## Prior for alpha, beta
  V0 ~ normal(0,0.1); ## Prior for initial #
  for (i in 1:Ntrain) {
    V0s[1]=V0*odds[i,1];
    V0s[2]=V0*odds[i,2];
    ## Solve ODE at current parameter values
    V_hat = integrate_ode_rk45(ode_diff, V0s, t0, ts_sums,
      theta, x_r, x_i);
      ## Likelihood
      for (jj in 1:Tinds) {
        j=ts_inds_inds[jj];
        V1[i,jj] ~ normal(V_hat[j,1],sigma/sqrt(n));
        V1stds[i,jj] ~ normal(0,sigma/sqrt(n));
        V2[i,jj] ~ normal(V_hat[j,2],sigma/sqrt(n));
        V2stds[i,jj] ~ normal(0,sigma/sqrt(n));
      }
      for (jj in 1:Tsums) {
        j = jj;
        Vsums[i,jj] ~ normal(V_hat[j,1]+V_hat[j,2],sigma/sqrt(n));
        Vsumstds[i,jj] ~ normal(0,sigma/sqrt(n));
      }
  }
}

generated quantities {
  real V_hat[Tsums,2];
  real V0s[2];
  int j;
  int i;
  real LL;
  for (ii in 1:Ntest) {
    i = ii;
    V0s[1]=V0*hdodds[i,1];
    V0s[2]=V0*hdodds[i,2];
    ## Solve ODE at current parameter values
    V_hat = integrate_ode_rk45(ode_diff, V0s, t0, ts_sums,
      theta, x_r, x_i);
      ## Likelihood
      LL = 0;
      for (jj in 1:Tinds) {
        j=ts_inds_inds[jj];
        LL = LL + normal_lpdf(hdV1[i,jj] | V_hat[j,1],sigma/sqrt(n));
        LL = LL + normal_lpdf(hdV1stds[i,jj] | 0,sigma/sqrt(n));
        LL = LL + normal_lpdf(hdV2[i,jj] | V_hat[j,2],sigma/sqrt(n));
        LL = LL + normal_lpdf(hdV2stds[i,jj] | 0,sigma/sqrt(n));
      }
      for (jj in 1:Tsums) {
        j = jj;
        LL = LL + normal_lpdf(hdVsums[i,jj] | V_hat[j,1]+V_hat[j,2],sigma/sqrt(n));
        LL = LL + normal_lpdf(hdVsumstds[i,jj] | 0,sigma/sqrt(n));
      }
  }
}

## lpdf
## try hold out
## waic / loo-cv "loo" in R need to feed sampled log likelihood on real data
