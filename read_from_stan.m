function [rc,kc,lc,rr,kr,lr,v0,sig] = read_from_stan()
  rc = csvread('cmdstan/bin/rc.csv');
  kc = csvread('cmdstan/bin/kc.csv');
  lc = csvread('cmdstan/bin/lc.csv');
  rr = csvread('cmdstan/bin/rr.csv');
  kr = csvread('cmdstan/bin/kr.csv');
  lr = csvread('cmdstan/bin/lr.csv');
  v0 = csvread('cmdstan/bin/v0.csv');
  sig = csvread('cmdstan/bin/sig.csv');
  % rc=rc(1:100);
  % kc=kc(1:100);
  % lc=lc(1:100);
  % rr=rr(1:100);
  % kr=kr(1:100);
  % lr=lr(1:100);
  % v0=v0(1:100);
  % sig=sig(1:100);
end
