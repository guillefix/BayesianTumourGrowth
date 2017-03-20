x0 = ones(1,8);
x=x0;
V0=[0.031;0.015];
tspan=[0,15];
% sampled_times=1:0.5:15;
% sampled_times_individual_data=[10,15];
sampled_times_individual_data=[10,15];
sampled_times_sum_data=[3 5 7 10 13 15];
standard_deviation_noise=0.05;

odds_of_initial_values = [1 0;
        3 1;
        1 1;
        1 3;
        0 1];

% neg log likelihood
V_sum_data_average_PC3_0Gy = [0.0636679167 0.0661274167 0.0987616667 0.16010675 0.25829375 0.280966;
    0.0574364167 0.0799836667 0.1146740833 0.2577545833 0.41375175 0.4119126667;
    0.0589646667 0.0776658333 0.1201630833 0.27822025 0.5145295 0.557941;
    0.0544568333 0.08331 0.1326801667 0.35862575 0.5359255 0.523286;
    0.0503110833 0.0921101667 0.1695063333 0.4100345833	0.5936558333 0.6179516667];

V_individual_control_data_PC3_0Gy = [0.99 0.98;
    0.56 0.45;
    0.32 0.18;
    0.12 0.06;
    0.00 0.00];

V_individual_resistant_data_PC3_0Gy = [0.00 0.00;
    0.43 0.53;
    0.67 0.80;
    0.87 0.92;
    0.99 0.97];

V_sum_data_sd_PC3_0Gy = [0.009479858 0.016532594 0.0242550332 0.0714584987 0.0963189337 0.0344588741;
    0.0073424551 0.0101382147 0.0258780892 0.0847120964 0.0766670382 0.0301724037;
    0.0056577373 0.0083716803 0.0113536203 0.0565088006 0.0929378678 0.0378267744;
    0.0077258534 0.0156706135 0.0118896548 0.0923792651 0.1044882296 0.0394104073;
    0.0080204573 0.0246143359 0.0285837478 0.1238795509 0.1385946846 0.0198155617];

V_individual_control_sd_PC3_0Gy = [0.001 0.01;
    0.01 0.05;
    0.02 0.01;
    0.01 0.01;
    0.001 0.001];

V_individual_resistant_sd_PC3_0Gy = [0.001 0.001;
    0.01 0.05;
    0.02 0.01;
    0.01 0.01;
    0.001 0.01];

ratiosNumber = 5;

NLLsig = @(x) findLogLikRealerData(x,odds_of_initial_values,V_sum_data_average_PC3_0Gy,V_sum_data_sd_PC3_0Gy,...
    V_individual_control_data_PC3_0Gy,V_individual_control_sd_PC3_0Gy,...
    V_individual_resistant_data_PC3_0Gy,V_individual_resistant_sd_PC3_0Gy,...
    ratiosNumber,x(8),sampled_times_individual_data, ...
    sampled_times_sum_data, x(7),tspan);

NumSample=100000;

xs = [];
sigma=0.02;

for i=1:NumSample
   if (mod(i,NumSample/100) == 0) disp(i);
   end
   xnew = proposalSample(x,sigma);
   if (xnew(2) < 0 || xnew(5) < 0 || xnew(1) < 0 || xnew(4) < 0 || xnew(8) <0) continue;
   end
   r = exp(-NLLsig(xnew))*proposalDist(xnew,0.5*x0,2)*proposalDist(x,xnew,sigma)...
       /(exp(-NLLsig(x))*proposalDist(x,0.5*x0,2)*proposalDist(xnew,x,sigma));
   if (rand < r)
       x = xnew;
       xs = [xs; x];
   end
end
xs=xs(end/2:end,:);
sizes=size(xs);
sum(xs,1)/(sizes(1))