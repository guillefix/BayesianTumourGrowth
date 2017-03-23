functions {
  real[] ode_diff(real t,real[] V,
      real[] theta,real[] x_r,int[] x_i){
    real dVdt[2];
    dVdt[1] = V[1] * (theta[1] * log(theta[2]/V[1])//^theta[7])
      - theta[3] * V[2]);
    dVdt[2] = V[2] * (theta[4] * log(theta[5]/V[2])//^theta[7])
      - theta[6] * V[1]);
    return dVdt;
  }
}

data {
  int N;
  // int Ntest;
  int n;
  int odds[N,2];
  int<lower=1> Tinds;
  int<lower=1> Tsums;
  real V1[N,Tinds];
  real V1stds[N,Tinds];
  real V2[N,Tinds];
  real V2stds[N,Tinds];
  real Vsums[N,Tsums];
  real Vsumstds[N,Tsums];

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
  // real<lower=1> alpha;
}
transformed parameters {
  real theta[6];
  theta[1] = theta1[1];
  theta[2] = theta1[2];
  theta[3] = lambda[1];
  theta[4] = theta1[3];
  theta[5] = theta1[4];
  theta[6] = lambda[2];
  // theta[7] = alpha;
}

model {
  real V_hat[Tsums,2];
  real V0s[2];
  int j;
  sigma ~ cauchy(0,1); ## Prior
  // for (k in 1:6) {
  //   theta[k] ~ normal(0.5,0.5); ## Prior for thetas
  // }
  // theta[1] ~ lognormal(log(0.5)-0.5,1);
  // theta[2] ~ lognormal(log(0.5)-0.5,1);
  // theta[3] ~ normal(0.5,1);
  // theta[4] ~ lognormal(log(0.5)-0.5,1);
  // theta[5] ~ lognormal(log(0.5)-0.5,1);
  // theta[6] ~ normal(0.5,1);
  // V0 ~ lognormal(log(0.5)-0.5,1); ## Prior for initial #
  theta[1] ~ lognormal(0,0.5);
  theta[2] ~ lognormal(0,0.5);
  theta[3] ~ normal(0.5,2);
  theta[4] ~ lognormal(0,0.5);
  theta[5] ~ lognormal(0,0.5);
  theta[6] ~ normal(0.5,2);
  V0 ~ normal(0,0.25); ## Prior for initial #

  for (i in 1:N) {
    V0s[1]=V0*odds[i,1]+0.001;
    V0s[2]=V0*odds[i,2]+0.001;
    ## Solve ODE at current parameter values
    V_hat = integrate_ode_rk45(ode_diff, V0s, t0, ts_sums,
      theta, x_r, x_i);
      print(theta);
      print(V0s);
      print(V_hat);
      ## Likelihood
      for (jj in 1:Tinds) {
        j=ts_inds_inds[jj];
        print(V_hat[j,1]);
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
  real LL[Tsums];
  for (ii in 1:N) {
    i = ii;
    V0s[1]=V0*odds[i,1]+0.001;
    V0s[2]=V0*odds[i,2]+0.001;
    ## Solve ODE at current parameter values
    V_hat = integrate_ode_rk45(ode_diff, V0s, t0, ts_sums,
      theta, x_r, x_i);
      ## Likelihood
      for (jj in 1:Tsums) {
        j = jj;
        LL[j] = 0;
        LL[j] = LL[j] + normal_lpdf(Vsums[i,jj] | V_hat[j,1]+V_hat[j,2],sigma/sqrt(n));
        LL[j] = LL[j] + normal_lpdf(Vsumstds[i,jj] | 0,sigma/sqrt(n));
        // print(LL[j]);
      }
      for (jj in 1:Tinds) {
        j=ts_inds_inds[jj];
        LL[j] = LL[j] + normal_lpdf(V1[i,jj] | V_hat[j,1],sigma/sqrt(n));
        LL[j] = LL[j] + normal_lpdf(V1stds[i,jj] | 0,sigma/sqrt(n));
        LL[j] = LL[j] + normal_lpdf(V2[i,jj] | V_hat[j,2],sigma/sqrt(n));
        LL[j] = LL[j] + normal_lpdf(V2stds[i,jj] | 0,sigma/sqrt(n));
        // print(LL[j]);

      }
  }
}
