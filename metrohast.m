x= [0.5,0.5,0.5,0.5,0.5,0.5];
V0=[0.031;0.015];
tspan=[0,15];
sampled_times=1:0.5:15;
standard_deviation_noise=0.01;
[Vdata1, Vdata2] = generate_fake_data(V0, tspan, sampled_times,standard_deviation_noise);
% neg log likelihood
NLL = @(x) findLogLikData(x,Vdata1,Vdata2,V0,sampled_times,standard_deviation_noise,tspan);

NumSample=5000;

x = 0.5*ones(1,6);
xs = [];
sigma=0.1;

for i=1:NumSample
   if (mod(i,100) == 0) disp(i);
   end
   xnew = proposalSample(x,sigma);
   if (x(2) < 0 || x(5) < 0 || x(1) < 0 || x(4) < 0) continue;
   end
   r = exp(-NLL(xnew))*proposalDist(x,xnew,sigma)/(exp(-NLL(x))*proposalDist(xnew,x,sigma));
   if (rand < r)
       x = xnew;
       xs = [xs; x];
   end
   
end
xs=xs(end/2:end,:);
sizes=size(xs);
sum(xs,1)/(sizes(1))