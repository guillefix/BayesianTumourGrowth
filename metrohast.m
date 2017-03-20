x0 = ones(1,7);
x=x0;
V0=[0.031;0.015];
tspan=[0,15];
% sampled_times=1:0.5:15;
% sampled_times_individual_data=[10,15];
sampled_times_individual_data=1:0.5:15;
sampled_times_sum_data=1:0.5:15;
standard_deviation_noise=0.05;
% [Vdata1, Vdata2] = generate_fake_data(V0, tspan, sampled_times,standard_deviation_noise);
[V_data_sum,Vdata1, Vdata2] = generate_fake_realistic_data(V0, tspan,...
    sampled_times_individual_data, sampled_times_sum_data, standard_deviation_noise);
% neg log likelihood
% NLL = @(x) findLogLikData(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan);
NLL = @(x) findLogLikRealData(x,V_data_sum,Vdata1,Vdata2,V0,...
    sampled_times_individual_data, sampled_times_sum_data,standard_deviation_noise,tspan);
NLLsig = @(x) findLogLikRealData(x,V_data_sum,Vdata1,Vdata2,V0,...
    sampled_times_individual_data, sampled_times_sum_data,x(7),tspan);

NumSample=10000;

xs = [];
sigma=0.02;

for i=1:NumSample
   if (mod(i,NumSample/100) == 0) disp(i);
   end
   xnew = proposalSample(x,sigma);
   if (xnew(2) < 0 || xnew(5) < 0 || xnew(1) < 0 || xnew(4) < 0) continue;
   end
   r = exp(-NLLsig(xnew))*proposalDist(xnew,x0,1)*proposalDist(x,xnew,sigma)...
       /(exp(-NLLsig(x))*proposalDist(x,x0,1)*proposalDist(xnew,x,sigma));
   if (rand < r)
       x = xnew;
       xs = [xs; x];
   end
end
xs=xs(end/2:end,:);
sizes=size(xs);
sum(xs,1)/(sizes(1))