function [Vdata1, Vdata2] = generate_fake_data(V0, tspan, sampled_times,...
    standard_deviation_noise)

rc=0.236;
Kc=0.473;
lc=0.655;
rr=0.401;
Kr=0.773;
lr=0.480;

dvdt = @(t,V) [V(1)*(rc*(1-V(1)/Kc)-lr*V(2));V(2)*(rr*(1-V(2)/Kr)-lc*V(1))] ;
% V0=[1;2];
% tspan=[0,5];
[t,V]=ode45(dvdt,tspan,V0);

<<<<<<< HEAD:generate_fake_data.m
% plot(t,V(:,1))
% hold on
% plot(t,V(:,2))
=======
plot(t,V(:,1))
hold on
plot(t,V(:,2))
hold off
>>>>>>> 362f48293cf59e70b38ca9fc90c4bd9755bdf859:Old/generate_fake_data.m

% ts=1:0.01:5;

Vdata1=interp1(t,V(:,1),sampled_times);
Vdata1=Vdata1+standard_deviation_noise*randn(1,length(Vdata1));
Vdata2=interp1(t,V(:,2),sampled_times);
Vdata2=Vdata2+standard_deviation_noise*randn(1,length(Vdata2));
end