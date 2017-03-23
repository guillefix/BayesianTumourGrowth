library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data <- list(N = 5,
             n = 12,
             Tsums = 6,
             Tinds = 2,
             ts_inds_inds = c(4,6),
             ts_sums = c(3,5,7,10,13,15),
             t0 = 0,
             odds = structure(c(1,0,3,1,1,1,1,3,0,1), .Dim = c(5,2)),
             Vsums = structure(c(0.0636679167,0.0661274167,0.0987616667,0.16010675,0.25829375,0.280966,0.0574364167,0.0799836667,0.1146740833,0.2577545833,0.41375175,0.4119126667,0.0589646667,0.0776658333,0.1201630833,0.27822025,0.5145295,0.557941,0.0544568333,0.08331,0.1326801667,0.35862575,0.5359255,0.523286,0.0503110833,0.0921101667,0.1695063333,0.4100345833,0.5936558333,0.6179516667), .Dim = c(5,6)),
             Vsumstds = structure(c(0.009479858,0.016532594,0.0242550332,0.0714584987,0.0963189337,0.0344588741,0.0073424551,0.0101382147,0.0258780892,0.0847120964,0.0766670382,0.0301724037,0.0056577373,0.0083716803,0.0113536203,0.0565088006,0.0929378678,0.0378267744,0.0077258534,0.0156706135,0.0118896548,0.0923792651,0.1044882296,0.0394104073,0.0080204573,0.0246143359,0.0285837478,0.1238795509,0.1385946846,0.0198155617), .Dim = c(5,6)),
             V1 = structure(c(0.99,0.98,0.56,0.45,0.32,0.18,0.12,0.06,0.00,0.00), .Dim = c(5,2)),
             V1stds = structure(c(0.001,0.01,0.01,0.05,0.02,0.01,0.01,0.01,0.001,0.001), .Dim = c(5,2)),
             V2 = structure(c(0.00,0.00,0.43,0.53,0.67,0.80,0.87,0.92,0.99,0.97), .Dim = c(5,2)),
             V2stds = structure(c(0.001,0.001,0.01,0.05,0.02,0.01,0.01,0.01,0.001,0.01), .Dim = c(5,2)))

data_hd <- list(Ntrain = 4,
              Ntest = 1,
             n = 12,
             Tsums = 6,
             Tinds = 2,
             ts_inds_inds = c(4,6),
             ts_sums = c(3,5,7,10,13,15),
             t0 = 0,
             odds = structure(c(1,0,3,1,1,1,1,3), .Dim = c(4,2)),
             Vsums = structure(c(0.0636679167,0.0661274167,0.0987616667,0.16010675,0.25829375,0.280966,0.0574364167,0.0799836667,0.1146740833,0.2577545833,0.41375175,0.4119126667,0.0589646667,0.0776658333,0.1201630833,0.27822025,0.5145295,0.557941,0.0544568333,0.08331,0.1326801667,0.35862575,0.5359255,0.523286), .Dim = c(4,6)),
             Vsumstds = structure(c(0.009479858,0.016532594,0.0242550332,0.0714584987,0.0963189337,0.0344588741,0.0073424551,0.0101382147,0.0258780892,0.0847120964,0.0766670382,0.0301724037,0.0056577373,0.0083716803,0.0113536203,0.0565088006,0.0929378678,0.0378267744,0.0077258534,0.0156706135,0.0118896548,0.0923792651,0.1044882296,0.0394104073), .Dim = c(4,6)),
             V1 = structure(c(0.99,0.98,0.56,0.45,0.32,0.18,0.12,0.06), .Dim = c(4,2)),
             V1stds = structure(c(0.001,0.01,0.01,0.05,0.02,0.01,0.01,0.01), .Dim = c(4,2)),
             V2 = structure(c(0.00,0.00,0.43,0.53,0.67,0.80,0.87,0.92), .Dim = c(4,2)),
             V2stds = structure(c(0.001,0.001,0.01,0.05,0.02,0.01,0.01,0.01), .Dim = c(4,2)),
             hdodds = structure(c(0,1), .Dim = c(1,2)),
             hdVsums = structure(c(0.0503110833,0.0921101667,0.1695063333,0.4100345833,0.5936558333,0.6179516667), .Dim = c(1,6)),
             hdVsumstds = structure(c(0.0080204573,0.0246143359,0.0285837478,0.1238795509,0.1385946846,0.0198155617), .Dim = c(1,6)),
             hdV1 = structure(c(0.00,0.00), .Dim = c(1,2)),
             hdV1stds = structure(c(0.001,0.001), .Dim = c(1,2)),
             hdV2 = structure(c(0.99,0.97), .Dim = c(1,2)),
             hdV2stds = structure(c(0.001,0.01), .Dim = c(1,2)))


fit <- stan(file = 'real_model_hd.stan', data = data_hd, iter = 10000, chains = 4)
plot(fit)
print(fit)


library("loo")

logliks = extract_log_lik(fit, parameter_name="LL")

loo1=loo(logliks)

sumlogliks = rowSums(logliks);

hist(sumlogliks, breaks=50)
hist(exp(sumlogliks), breaks=50)

fit_gomp <- stan(file = 'real_model_gomp_hd.stan', data = data_hd, iter = 10000, chains = 4)
plot(fit_gomp)
print(fit_gomp)

logliks_gomp = extract_log_lik(fit_gomp, parameter_name="LL")

loo_gomp=loo(logliks_gomp)

sumlogliks_gomp = rowSums(logliks_gomp);

hist(sumlogliks_gomp, breaks=50)
hist(exp(sumlogliks_gomp), breaks=50)

compare(loo1, loo_gomp)


rr = extract(fit, "theta[1]")
rr=unlist(rr, use.names=FALSE)
write.csv(rr, "rr.csv", row.names=FALSE)

kr = extract(fit, "theta[2]")
kr=unlist(kr, use.names=FALSE)
write.csv(rr, "kr.csv", row.names=FALSE)

lr = extract(fit, "theta[3]")
lr=unlist(lr, use.names=FALSE)
write.csv(rr, "lr.csv", row.names=FALSE)

rc = extract(fit, "theta[4]")
rc=unlist(rc, use.names=FALSE)
write.csv(rr, "rc.csv", row.names=FALSE)

kc = extract(fit, "theta[5]")
kc=unlist(kc, use.names=FALSE)
write.csv(rr, "kc.csv", row.names=FALSE)

lc = extract(fit, "theta[6]")
lc=unlist(lc, use.names=FALSE)
write.csv(rr, "lc.csv", row.names=FALSE)

sig = extract(fit, "sigma")
sig=unlist(sig, use.names=FALSE)
write.csv(sig, "sig.csv", row.names=FALSE)

v0 = extract(fit, "V0")
v0=unlist(v0, use.names=FALSE)
write.csv(v0, "v0.csv", row.names=FALSE)

#####

rr = extract(fit_gomp, "theta[1]")
rr=unlist(rr, use.names=FALSE)
write.csv(rr, "rr.csv", row.names=FALSE)

kr = extract(fit_gomp, "theta[2]")
kr=unlist(kr, use.names=FALSE)
write.csv(rr, "kr.csv", row.names=FALSE)

lr = extract(fit_gomp, "theta[3]")
lr=unlist(lr, use.names=FALSE)
write.csv(rr, "lr.csv", row.names=FALSE)

rc = extract(fit_gomp, "theta[4]")
rc=unlist(rc, use.names=FALSE)
write.csv(rr, "rc.csv", row.names=FALSE)

kc = extract(fit_gomp, "theta[5]")
kc=unlist(kc, use.names=FALSE)
write.csv(rr, "kc.csv", row.names=FALSE)

lc = extract(fit_gomp, "theta[6]")
lc=unlist(lc, use.names=FALSE)
write.csv(rr, "lc.csv", row.names=FALSE)

sig = extract(fit_gomp, "sigma")
sig=unlist(sig, use.names=FALSE)
write.csv(sig, "sig.csv", row.names=FALSE)

v0 = extract(fit_gomp, "V0")
v0=unlist(v0, use.names=FALSE)
write.csv(v0, "v0.csv", row.names=FALSE)

