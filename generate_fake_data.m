function [Vdata1, Vdata2] = generate_fake_data()

rc=1;
Kc=1;
lc=1;
rr=2;
Kr=2;
lr=2;

dvdt = @(t,V) [V(1)*(rc*(1-V(1)/Kc)-lr*V(2));V(2)*(rr*(1-V(2)/Kr)-lc*V(1))] ;
V0=[1;2];
tspan=[0,5];
[t,V]=ode45(dvdt,tspan,V0);
plot(t,V(:,1))
hold on
plot(t,V(:,2))

ts=1:0.01:5;
Vdata1=interp1(t,V(:,1),ts);
Vdata1=Vdata1+0.01*randn(1,length(Vdata1));
Vdata2=interp1(t,V(:,2),ts);
Vdata2=Vdata2+0.01*randn(1,length(Vdata2));
end