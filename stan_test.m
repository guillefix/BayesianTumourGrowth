
x = ones(1,7);
V0=[0.031;0.015];
tspan=[0,15];
% sampled_times=1:0.5:15;
% sampled_times_individual_data=[10,15];
sampled_times_individual_data=1:0.5:15;
sampled_times_sum_data=1:0.5:15;
standard_deviation_noise=0.05;
% [Vdata1, Vdata2] = generate_fake_data(V0, tspan, sampled_times,standard_deviation_noise);
[V_data_sum,Vdata1, Vdata2] = generate_fake_realistic_data(V0, tspan, sampled_times_individual_data, sampled_times_sum_data, standard_deviation_noise);

fit1 = stan('file','model.stan','data',[Vdata1,Vdata2],'iter',1000,'chains',4);